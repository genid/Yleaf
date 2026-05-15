use serde::Serialize;
use std::fs;
use std::path::Path;

#[derive(Serialize, Clone)]
pub struct HgPrediction {
    pub tree: String,
    pub sample_name: String,
    pub haplogroup: String,
    pub hg_marker: String,
    pub total_reads: i64,
    pub valid_markers: i64,
    pub qc_score: f64,
    pub qc1: f64,
    pub qc2: f64,
    pub qc3: f64,
}

#[derive(Serialize, Clone)]
pub struct MixtureContributor {
    pub rank: u32,
    pub haplogroup: String,
    pub ratio: f64,
    pub markers: u32,
    pub qc: f64,
}

#[derive(Serialize, Clone)]
pub struct MixtureSample {
    pub tree: String,
    pub sample_name: String,
    pub common_ancestor: String,
    pub contributors: Vec<MixtureContributor>,
}

#[derive(Serialize)]
pub struct MarkerStats {
    pub derived: usize,
    pub ancestral: usize,
    pub not_covered: usize,
}

#[derive(Serialize)]
pub struct JobResults {
    pub predictions: Vec<HgPrediction>,
    pub marker_stats: Option<MarkerStats>,
    pub mixtures: Vec<MixtureSample>,
}

/// Parse a single-tree .hg file (Sample_name, Hg, Hg_marker, Total_reads,
/// Valid_markers, QC-score, QC-1, QC-2, QC-3).
fn parse_hg_file(path: &Path, tree: &str) -> Result<Vec<HgPrediction>, String> {
    let content = fs::read_to_string(path).map_err(|e| format!("{path:?}: {e}"))?;
    let mut out = Vec::new();
    for (i, line) in content.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 9 {
            continue;
        }
        out.push(HgPrediction {
            tree: tree.to_string(),
            sample_name: c[0].to_string(),
            haplogroup: c[1].to_string(),
            hg_marker: c[2].to_string(),
            total_reads: c[3].parse().unwrap_or(0),
            valid_markers: c[4].parse().unwrap_or(0),
            qc_score: c[5].parse().unwrap_or(0.0),
            qc1: c[6].parse().unwrap_or(0.0),
            qc2: c[7].parse().unwrap_or(0.0),
            qc3: c[8].parse().unwrap_or(0.0),
        });
    }
    Ok(out)
}

/// Scan sample .out file and count D / A / not-covered markers.
/// Column 10 (0-indexed) holds the state: D, A, or anything else (NA / no coverage).
fn parse_marker_stats(sample_dir: &Path) -> Option<MarkerStats> {
    let dir_name = sample_dir.file_name()?.to_str()?;
    // Try yfull-style .out first, then fallback to any .out in the dir
    let out_path = sample_dir.join(format!("{dir_name}.out"));
    let content = fs::read_to_string(&out_path).ok()?;

    let (mut d, mut a, mut nc) = (0usize, 0usize, 0usize);
    for (i, line) in content.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let mut cols = line.splitn(12, '\t');
        // advance to column 10
        for _ in 0..10 {
            cols.next();
        }
        match cols.next() {
            Some("D") => d += 1,
            Some("A") => a += 1,
            _ => nc += 1,
        }
    }
    Some(MarkerStats {
        derived: d,
        ancestral: a,
        not_covered: nc,
    })
}

/// Parse a .mix file produced by Yleaf mixture mode.
fn parse_mix_file(path: &Path) -> Option<MixtureSample> {
    let content = fs::read_to_string(path).ok()?;
    let mut tree = String::new();
    let mut common_ancestor = String::new();
    let mut contributors: Vec<MixtureContributor> = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') {
            let body = line.trim_start_matches('#');
            let mut kv = body.splitn(2, '\t');
            let key = kv.next().unwrap_or("").trim();
            let val = kv.next().unwrap_or("").trim();
            match key {
                "tree" => tree = val.to_string(),
                "common_ancestor" => common_ancestor = val.to_string(),
                _ => {}
            }
            continue;
        }
        if line.starts_with("contributor\t") {
            continue; // header row
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() >= 5 {
            contributors.push(MixtureContributor {
                rank: c[0].parse().unwrap_or(0),
                haplogroup: c[1].to_string(),
                ratio: c[2].parse().unwrap_or(0.0),
                markers: c[3].parse().unwrap_or(0),
                qc: c[4].parse().unwrap_or(0.0),
            });
        }
    }

    if contributors.is_empty() {
        return None;
    }

    // Derive sample name from the file stem (strip trailing .{tree}.mix or .mix)
    let stem = path.file_stem()?.to_str()?;
    let sample_name = if !tree.is_empty() && stem.ends_with(&format!(".{tree}")) {
        stem[..stem.len() - tree.len() - 1].to_string()
    } else {
        stem.to_string()
    };

    Some(MixtureSample { tree, sample_name, common_ancestor, contributors })
}

/// Collect .mix files from a single directory level (non-recursive).
fn scan_mix_dir(dir: &Path) -> Vec<MixtureSample> {
    let Ok(rd) = fs::read_dir(dir) else { return Vec::new() };
    rd.filter_map(|e| e.ok())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("mix"))
        .filter_map(|e| parse_mix_file(&e.path()))
        .collect()
}

/// Collect all .mix files under output_dir (top level + one level of subdirectories).
fn collect_mix_files(output_dir: &Path) -> Vec<MixtureSample> {
    let mut results = scan_mix_dir(output_dir);
    if let Ok(rd) = fs::read_dir(output_dir) {
        for entry in rd.filter_map(|e| e.ok()) {
            if entry.path().is_dir() {
                results.extend(scan_mix_dir(&entry.path()));
            }
        }
    }
    results
}

/// Pick the first sample subdirectory in output_dir for marker stats.
fn first_sample_dir(output_dir: &Path) -> Option<std::path::PathBuf> {
    fs::read_dir(output_dir)
        .ok()?
        .filter_map(|e| e.ok())
        .find(|e| e.path().is_dir())
        .map(|e| e.path())
}

/// Parse the best primary haplogroup + QC from an output directory.
/// Used internally when updating the DB after a run.
pub fn extract_primary_result(
    output_dir: &Path,
    trees: &[&str],
) -> Option<(String, f64, f64, f64, f64, i64, i64)> {
    // In single-tree mode, hg_prediction.hg holds the result.
    // In multi-tree mode, use the first tree's per-tree file.
    let hg_path = if trees.len() == 1 {
        output_dir.join("hg_prediction.hg")
    } else {
        output_dir.join(format!("hg_prediction_{}.hg", trees[0]))
    };
    let preds = parse_hg_file(&hg_path, trees[0]).ok()?;
    let p = preds.into_iter().next()?;
    Some((p.haplogroup, p.qc_score, p.qc1, p.qc2, p.qc3, p.total_reads, p.valid_markers))
}

/// Tauri command: return full results for a finished job.
#[tauri::command]
pub fn get_job_results(output_dir: String, trees: String) -> Result<JobResults, String> {
    let out_dir = Path::new(&output_dir);
    let tree_list: Vec<&str> = trees.split(',').filter(|s| !s.is_empty()).collect();

    let predictions = if tree_list.len() <= 1 {
        // Single tree — one hg_prediction.hg file
        let tree = tree_list.first().copied().unwrap_or("yfull");
        parse_hg_file(&out_dir.join("hg_prediction.hg"), tree)?
    } else {
        // Multi-tree — one per-tree file per tree
        let mut all = Vec::new();
        for tree in &tree_list {
            let path = out_dir.join(format!("hg_prediction_{tree}.hg"));
            if path.exists() {
                all.extend(parse_hg_file(&path, tree)?);
            }
        }
        all
    };

    let mut predictions = predictions;
    predictions.sort_by(|a, b| a.sample_name.cmp(&b.sample_name));

    let marker_stats = first_sample_dir(out_dir).and_then(|d| parse_marker_stats(&d));
    let mixtures = collect_mix_files(out_dir);

    Ok(JobResults {
        predictions,
        marker_stats,
        mixtures,
    })
}

/// Tauri command: return the ROOT→leaf path for each tree.
/// Looks in {output_dir}/{sample_name}/haplogroup_path_{tree}.json (per-sample, new format).
/// Falls back to {output_dir}/haplogroup_path_{tree}.json for older outputs.
#[tauri::command]
pub fn get_haplogroup_path(
    output_dir: String,
    sample_name: String,
    trees: Vec<String>,
) -> Result<std::collections::HashMap<String, serde_json::Value>, String> {
    let base = Path::new(&output_dir);
    let mut result = std::collections::HashMap::new();
    for tree in &trees {
        let per_sample = base.join(&sample_name).join(format!("haplogroup_path_{}.json", tree));
        let fallback   = base.join(format!("haplogroup_path_{}.json", tree));
        let path_file  = if per_sample.exists() { per_sample } else { fallback };
        if path_file.exists() {
            let content = fs::read_to_string(&path_file).map_err(|e| e.to_string())?;
            let value: serde_json::Value = serde_json::from_str(&content).map_err(|e| e.to_string())?;
            result.insert(tree.clone(), value);
        }
    }
    Ok(result)
}

/// Tauri command: rename a job in the DB.
#[tauri::command]
pub fn rename_job(
    state: tauri::State<'_, crate::DbState>,
    job_id: i64,
    name: String,
) -> Result<(), String> {
    let conn = state.0.lock().map_err(|e| e.to_string())?;
    crate::db::rename_job(&conn, job_id, &name).map_err(|e| e.to_string())
}

/// Tauri command: remove a job record from the DB.
#[tauri::command]
pub fn delete_job(
    state: tauri::State<'_, crate::DbState>,
    job_id: i64,
) -> Result<(), String> {
    let conn = state.0.lock().map_err(|e| e.to_string())?;
    crate::db::delete_job(&conn, job_id).map_err(|e| e.to_string())
}

/// Tauri command: list all jobs from the DB.
#[tauri::command]
pub fn list_jobs(
    state: tauri::State<'_, crate::DbState>,
) -> Result<Vec<crate::db::Job>, String> {
    let conn = state.0.lock().map_err(|e| e.to_string())?;
    crate::db::list_jobs(&conn).map_err(|e| e.to_string())
}
