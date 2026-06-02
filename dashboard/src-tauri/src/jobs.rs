use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use tauri::{AppHandle, Emitter, Manager};
use tauri_plugin_shell::process::CommandEvent;
use tauri_plugin_shell::ShellExt;

use crate::db;
use crate::{ActiveChild, DbState};

#[derive(Clone, serde::Serialize)]
pub struct ProgressPayload {
    pub job_id: i64,
    pub line: String,
}

#[derive(Clone, serde::Serialize)]
pub struct ProgressBarPayload {
    pub job_id: i64,
    pub stage: String,
    pub current: u32,
    pub total: u32,
}

#[derive(Clone, serde::Serialize)]
pub struct DonePayload {
    pub job_id: i64,
    pub output_dir: String,
}

fn now_secs() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs() as i64
}

/// Classify a single file by extension into its Yleaf CLI flag.
fn flag_for_file(p: &std::path::Path) -> Option<&'static str> {
    let name = p.file_name().unwrap_or_default().to_string_lossy();
    if name.ends_with(".ped") {
        return Some("-plink");
    }
    if name.ends_with(".bed") && p.with_extension("bim").exists() {
        return Some("-plink");
    }
    if name.ends_with(".vcf") || name.ends_with(".vcf.gz") {
        return Some("-vcf");
    }
    if name.ends_with(".cram") {
        return Some("-cram");
    }
    if name.ends_with(".fastq")
        || name.ends_with(".fastq.gz")
        || name.ends_with(".fq")
        || name.ends_with(".fq.gz")
    {
        return Some("-fastq");
    }
    if name.ends_with(".bam") {
        return Some("-bam");
    }
    None
}

/// Detect ALL eligible input types for a file or directory path.
/// File: returns one element. Directory: returns one element per type found
/// (BAM, CRAM, VCF, FASTQ, PLINK). For per-file flags (-vcf, -plink) the
/// resolved_path is the first matching file in the dir; for directory-aware
/// flags (-bam, -cram, -fastq) it is the directory itself.
fn resolve_inputs_all(input: &str) -> Result<Vec<(&'static str, String)>, String> {
    let p = std::path::Path::new(input);

    if !p.is_dir() {
        // Single file. Unknown extension still falls back to -bam for back-compat.
        let flag = flag_for_file(p).unwrap_or("-bam");
        return Ok(vec![(flag, input.to_string())]);
    }

    let mut entries: Vec<std::path::PathBuf> = std::fs::read_dir(p)
        .map_err(|e| e.to_string())?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect();
    entries.sort();

    // Stable order: plink > bam > cram > vcf > fastq (same as old single-pick priority).
    let mut out: Vec<(&'static str, String)> = Vec::new();
    let mut seen: std::collections::BTreeSet<&'static str> = std::collections::BTreeSet::new();
    let push = |flag: &'static str, path: String, seen: &mut std::collections::BTreeSet<&'static str>, out: &mut Vec<(&'static str, String)>| {
        if seen.insert(flag) {
            out.push((flag, path));
        }
    };

    for f in &entries {
        let Some(flag) = flag_for_file(f) else { continue };
        let path = match flag {
            // Per-file flags: pass the matching file path.
            "-plink" | "-vcf" => f.to_string_lossy().into_owned(),
            // Directory-aware flags: pass the parent directory once.
            _ => input.to_string(),
        };
        push(flag, path, &mut seen, &mut out);
    }

    if out.is_empty() {
        return Err(format!("No recognized input files found in: {input}"));
    }
    // Preserve the historical priority order for predictability.
    let rank = |f: &&'static str| match *f {
        "-plink" => 0,
        "-bam" => 1,
        "-cram" => 2,
        "-vcf" => 3,
        "-fastq" => 4,
        _ => 99,
    };
    out.sort_by_key(|(flag, _)| rank(flag));
    Ok(out)
}

/// Spawn the Yleaf sidecar — possibly more than once when the user-selected
/// directory contains several eligible input types (e.g. BAM + VCF). Each
/// type becomes its own job row with its own output subdirectory so the
/// per-type `hg_prediction.hg` / `report.json` outputs don't collide.
///
/// Returns the list of newly-created job IDs (one element in the common
/// single-type case, multiple when the directory had mixed types).
#[tauri::command]
pub async fn run_yleaf(
    app: AppHandle,
    db: tauri::State<'_, DbState>,
    active_child: tauri::State<'_, ActiveChild>,
    bam_path: String,
    output_dir: String,
    reference_genome: String,
    tree: Vec<String>,
    threads: u32,
    reads_threshold: u32,
    quality_thresh: u32,
    base_majority: u32,
    prediction_quality: f64,
    draw_haplogroups: bool,
    ancient_dna: bool,
    private_mutations: bool,
    collapsed_draw_mode: bool,
    mixture_mode: bool,
) -> Result<Vec<i64>, String> {
    let inputs = resolve_inputs_all(&bam_path)?;
    let multi = inputs.len() > 1;

    // Kill any zombie that targeted the base output_dir (previous-session leftover).
    // Do this ONCE up front so it doesn't clobber siblings we're about to spawn.
    {
        let mut map = active_child.0.lock().unwrap();
        let to_kill: Vec<i64> = map
            .iter()
            .filter(|(_, (dir, _))| {
                dir == &output_dir || std::path::Path::new(dir).starts_with(&output_dir)
            })
            .map(|(id, _)| *id)
            .collect();
        for id in to_kill {
            if let Some((_, child)) = map.remove(&id) {
                let _ = child.kill();
            }
        }
    }

    let mut spawned_ids: Vec<i64> = Vec::with_capacity(inputs.len());
    let trees_str = tree.join(",");
    let app_data_dir = app.path().app_data_dir().map_err(|e| e.to_string())?;

    for (input_flag, resolved_path) in inputs {
        // Per-type output subdir when multiple types were detected; otherwise
        // keep the user-chosen output_dir unchanged (back-compat).
        let job_output_dir = if multi {
            let sub = input_flag.trim_start_matches('-');
            format!("{}/{}", output_dir.trim_end_matches('/'), sub)
        } else {
            output_dir.clone()
        };
        if multi {
            std::fs::create_dir_all(&job_output_dir)
                .map_err(|e| format!("create {job_output_dir}: {e}"))?;
        }

        // Build CLI args for this type
        let mut args: Vec<String> = vec![
            input_flag.into(),
            resolved_path.clone(),
            "-o".into(),
            job_output_dir.clone(),
            "-rg".into(),
            reference_genome.clone(),
            "-t".into(),
            threads.to_string(),
            "-r".into(),
            reads_threshold.to_string(),
            "-q".into(),
            quality_thresh.to_string(),
            "-b".into(),
            base_majority.to_string(),
            "-pq".into(),
            prediction_quality.to_string(),
            "-tree".into(),
        ];
        args.extend(tree.clone());
        if draw_haplogroups {
            args.push("-dh".into());
        }
        if draw_haplogroups && collapsed_draw_mode {
            args.push("-hc".into());
        }
        if ancient_dna {
            args.push("-aDNA".into());
        }
        if private_mutations {
            args.push("-p".into());
        }
        if mixture_mode {
            args.push("-mix".into());
        }
        args.push("--report-json".into());
        args.push(format!("{}/report.json", job_output_dir));
        args.push("-force".into());

        // Display name: file stem; for mixed dirs annotate with the type tag
        // so the job list distinguishes "<dir>" entries.
        let base_name = std::path::Path::new(&resolved_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        let sample_name = if multi {
            format!("{} [{}]", base_name, input_flag.trim_start_matches('-'))
        } else {
            base_name
        };

        // Insert job record
        let db_arc = Arc::clone(&db.0);
        let job_id = {
            let conn = db_arc.lock().map_err(|e| e.to_string())?;
            db::insert_job(
                &conn, &sample_name, &resolved_path, &job_output_dir, &trees_str, now_secs(),
                &reference_genome, threads as i64, reads_threshold as i64,
                quality_thresh as i64, base_majority as i64, prediction_quality,
                draw_haplogroups, ancient_dna, private_mutations, collapsed_draw_mode,
            ).map_err(|e| e.to_string())?
        };

        // Spawn sidecar — pass persistent data dir so downloaded references survive restarts
        let sidecar = app
            .shell()
            .sidecar("yleaf")
            .map_err(|e| e.to_string())?
            .env("YLEAF_DATA_DIR", app_data_dir.to_string_lossy().as_ref())
            .env("PYTHONUNBUFFERED", "1")
            .args(&args);

        app.emit("yleaf-progress", ProgressPayload {
            job_id,
            line: "Starting Yleaf (first run may take a few minutes while the reference genome is prepared)...".to_string(),
        }).ok();
        let (mut rx, child) = sidecar.spawn().map_err(|e| e.to_string())?;
        active_child.0.lock().unwrap().insert(job_id, (job_output_dir.clone(), child));

        // Background task: stream events + update DB on completion
        let app_clone = app.clone();
        let active_child_arc = active_child.0.clone();
        let output_dir_clone = job_output_dir.clone();
        let trees_clone = tree.clone();
        let db_arc_for_task = Arc::clone(&db.0);
        tauri::async_runtime::spawn(async move {
            let db_arc = db_arc_for_task;
        while let Some(event) = rx.recv().await {
            match event {
                CommandEvent::Stdout(bytes) | CommandEvent::Stderr(bytes) => {
                    let line = String::from_utf8_lossy(&bytes).to_string();
                    // Intercept structured progress lines: "[PROGRESS] <stage> <n>/<total>"
                    if let Some(rest) = line.trim().strip_prefix("[PROGRESS] ") {
                        let parts: Vec<&str> = rest.splitn(2, ' ').collect();
                        if parts.len() == 2 {
                            let nums: Vec<&str> = parts[1].split('/').collect();
                            if nums.len() == 2 {
                                if let (Ok(current), Ok(total)) =
                                    (nums[0].parse::<u32>(), nums[1].parse::<u32>())
                                {
                                    app_clone.emit("yleaf-progress-bar", ProgressBarPayload {
                                        job_id,
                                        stage: parts[0].to_string(),
                                        current,
                                        total,
                                    }).ok();
                                    continue;
                                }
                            }
                        }
                    }
                    app_clone
                        .emit("yleaf-progress", ProgressPayload { job_id, line })
                        .ok();
                }
                CommandEvent::Error(e) => {
                    active_child_arc.lock().unwrap().remove(&job_id);
                    app_clone.emit("yleaf-error", (job_id, e)).ok();
                    let conn = db_arc.lock().unwrap();
                    db::update_job_error(&conn, job_id, now_secs()).ok();
                    app_clone.emit("job-updated", job_id).ok();
                    return;
                }
                CommandEvent::Terminated(status) => {
                    active_child_arc.lock().unwrap().remove(&job_id);
                    let end = now_secs();
                    if status.code == Some(0) {
                        let tree_refs: Vec<&str> = trees_clone.iter().map(|s| s.as_str()).collect();
                        if let Some((hg, qc, qc1, qc2, qc3, reads, markers)) =
                            crate::results::extract_primary_result(
                                std::path::Path::new(&output_dir_clone),
                                &tree_refs,
                            )
                        {
                            let conn = db_arc.lock().unwrap();
                            db::update_job_done(
                                &conn, job_id, &hg, qc, qc1, qc2, qc3, reads, markers, end,
                            )
                            .ok();
                        }
                        app_clone
                            .emit(
                                "yleaf-done",
                                DonePayload {
                                    job_id,
                                    output_dir: output_dir_clone,
                                },
                            )
                            .ok();
                    } else {
                        let conn = db_arc.lock().unwrap();
                        db::update_job_error(&conn, job_id, end).ok();
                        app_clone
                            .emit(
                                "yleaf-error",
                                (job_id, format!("Yleaf exited with code {:?}", status.code)),
                            )
                            .ok();
                    }
                    app_clone.emit("job-updated", job_id).ok();
                    return;
                }
                _ => {}
            }
        }
    });

        spawned_ids.push(job_id);
    }

    Ok(spawned_ids)
}

#[tauri::command]
pub fn kill_job(
    db: tauri::State<'_, DbState>,
    active_child: tauri::State<'_, ActiveChild>,
    job_id: i64,
) -> Result<(), String> {
    let entry = active_child.0.lock().unwrap().remove(&job_id);
    if let Some((_, child)) = entry {
        child.kill().map_err(|e| e.to_string())?;
    }
    let conn = db.0.lock().unwrap();
    db::update_job_error(&conn, job_id, now_secs()).ok();
    Ok(())
}
