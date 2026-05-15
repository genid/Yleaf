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

/// Detect input type from a file or directory path.
/// Returns (cli_flag, resolved_file_path).
/// Priority for directories: plink > bam > cram > vcf > fastq.
fn resolve_input(input: &str) -> Result<(&'static str, String), String> {
    let p = std::path::Path::new(input);

    if p.is_dir() {
        let mut entries: Vec<std::path::PathBuf> = std::fs::read_dir(p)
            .map_err(|e| e.to_string())?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .collect();
        entries.sort();

        // PLINK text (.ped)
        if let Some(f) = entries.iter().find(|f| {
            f.extension().and_then(|s| s.to_str()) == Some("ped")
        }) {
            return Ok(("-plink", f.to_string_lossy().into_owned()));
        }
        // PLINK binary (.bed + sibling .bim)
        if let Some(f) = entries.iter().find(|f| {
            f.extension().and_then(|s| s.to_str()) == Some("bed")
                && f.with_extension("bim").exists()
        }) {
            return Ok(("-plink", f.to_string_lossy().into_owned()));
        }
        // BAM — pass the directory; Yleaf processes all BAMs inside it
        if entries.iter().any(|f| {
            f.extension().and_then(|s| s.to_str()) == Some("bam")
        }) {
            return Ok(("-bam", input.to_string()));
        }
        // CRAM — same: pass the directory
        if entries.iter().any(|f| {
            f.extension().and_then(|s| s.to_str()) == Some("cram")
        }) {
            return Ok(("-cram", input.to_string()));
        }
        // VCF / VCF.GZ
        if let Some(f) = entries.iter().find(|f| {
            let name = f.file_name().unwrap_or_default().to_string_lossy();
            name.ends_with(".vcf") || name.ends_with(".vcf.gz")
        }) {
            return Ok(("-vcf", f.to_string_lossy().into_owned()));
        }
        // FASTQ — pass the directory; Yleaf processes all FASTQs inside it
        if entries.iter().any(|f| {
            let name = f.file_name().unwrap_or_default().to_string_lossy();
            name.ends_with(".fastq") || name.ends_with(".fastq.gz")
                || name.ends_with(".fq") || name.ends_with(".fq.gz")
        }) {
            return Ok(("-fastq", input.to_string()));
        }

        Err(format!("No recognized input files found in: {input}"))
    } else {
        let name = p.file_name().unwrap_or_default().to_string_lossy();
        if name.ends_with(".ped")
            || (name.ends_with(".bed") && p.with_extension("bim").exists())
        {
            Ok(("-plink", input.to_string()))
        } else if name.ends_with(".vcf") || name.ends_with(".vcf.gz") {
            Ok(("-vcf", input.to_string()))
        } else if name.ends_with(".cram") {
            Ok(("-cram", input.to_string()))
        } else if name.ends_with(".fastq") || name.ends_with(".fastq.gz")
            || name.ends_with(".fq") || name.ends_with(".fq.gz")
        {
            Ok(("-fastq", input.to_string()))
        } else {
            Ok(("-bam", input.to_string()))
        }
    }
}

/// Spawn the Yleaf sidecar with the given arguments.
///
/// Returns immediately with the new job ID. Progress is streamed via
/// "yleaf-progress" events; completion via "yleaf-done"; errors via "yleaf-error".
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
) -> Result<i64, String> {
    let (input_flag, resolved_path) = resolve_input(&bam_path)?;

    // Build CLI args
    let mut args: Vec<String> = vec![
        input_flag.into(),
        resolved_path.clone(),
        "-o".into(),
        output_dir.clone(),
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
    args.push(format!("{}/report.json", output_dir));
    args.push("-force".into());

    // Derive a display name from the input filename
    let sample_name = std::path::Path::new(&resolved_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();
    let trees_str = tree.join(",");

    // Insert job record
    let db_arc = Arc::clone(&db.0);
    let job_id = {
        let conn = db_arc.lock().map_err(|e| e.to_string())?;
        db::insert_job(
            &conn, &sample_name, &resolved_path, &output_dir, &trees_str, now_secs(),
            &reference_genome, threads as i64, reads_threshold as i64,
            quality_thresh as i64, base_majority as i64, prediction_quality,
            draw_haplogroups, ancient_dna, private_mutations, collapsed_draw_mode,
        ).map_err(|e| e.to_string())?
    };

    // Kill any existing process targeting the same output directory (zombie from a previous session)
    {
        let mut map = active_child.0.lock().unwrap();
        let to_kill: Vec<i64> = map
            .iter()
            .filter(|(_, (dir, _))| dir == &output_dir)
            .map(|(id, _)| *id)
            .collect();
        for id in to_kill {
            if let Some((_, child)) = map.remove(&id) {
                let _ = child.kill();
            }
        }
    }

    // Spawn sidecar — pass persistent data dir so downloaded references survive restarts
    let app_data_dir = app.path().app_data_dir().map_err(|e| e.to_string())?;
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
    active_child.0.lock().unwrap().insert(job_id, (output_dir.clone(), child));

    // Background task: stream events + update DB on completion
    let app_clone = app.clone();
    let active_child_arc = active_child.0.clone();
    let output_dir_clone = output_dir.clone();
    let trees_clone = tree.clone();
    tauri::async_runtime::spawn(async move {
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

    Ok(job_id)
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
