//! Stage multiple user-selected input files into a single directory so the
//! existing Yleaf directory-mode pipeline can pick them up. Files are hard-
//! linked when possible (instant, zero disk cost), otherwise symlinked, and
//! only copied as a last resort. Index files (.bai/.crai) next to BAM/CRAM
//! inputs are also linked over so samtools doesn't trip on missing indexes.

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

fn link_one(src: &Path, dst: &Path) -> Result<(), String> {
    if fs::hard_link(src, dst).is_ok() {
        return Ok(());
    }
    #[cfg(unix)]
    {
        if std::os::unix::fs::symlink(src, dst).is_ok() {
            return Ok(());
        }
    }
    #[cfg(windows)]
    {
        if std::os::windows::fs::symlink_file(src, dst).is_ok() {
            return Ok(());
        }
    }
    fs::copy(src, dst).map(|_| ()).map_err(|e| format!("copy {:?} -> {:?}: {}", src, dst, e))
}

fn link_index_if_present(src: &Path, dst_dir: &Path) {
    for ext in ["bai", "crai", "csi"] {
        let candidates = [
            src.with_extension(format!("{}.{}", src.extension().and_then(|s| s.to_str()).unwrap_or(""), ext)),
            src.with_extension(ext),
        ];
        for cand in candidates {
            if cand.exists() {
                let name = cand.file_name().unwrap();
                let _ = link_one(&cand, &dst_dir.join(name));
            }
        }
    }
}

#[tauri::command]
pub fn stage_batch_files(paths: Vec<String>) -> Result<String, String> {
    if paths.is_empty() {
        return Err("no files supplied".into());
    }
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let dir = std::env::temp_dir().join(format!("yleaf_batch_{}", ts));
    fs::create_dir_all(&dir).map_err(|e| format!("create {:?}: {}", dir, e))?;

    for p in &paths {
        let src = PathBuf::from(p);
        let name = src
            .file_name()
            .ok_or_else(|| format!("no filename: {}", p))?;
        link_one(&src, &dir.join(name))?;
        link_index_if_present(&src, &dir);
    }

    dir.to_str()
        .map(|s| s.to_string())
        .ok_or_else(|| "non-UTF8 path".into())
}

#[tauri::command]
pub fn cleanup_batch_dir(dir: String) -> Result<(), String> {
    let p = PathBuf::from(&dir);
    // Sanity guard: only delete dirs we created (prefix in temp).
    if !p.starts_with(std::env::temp_dir())
        || !p.file_name().and_then(|n| n.to_str()).map_or(false, |n| n.starts_with("yleaf_batch_"))
    {
        return Err(format!("refuse to delete {:?} — not a yleaf batch dir", p));
    }
    if p.exists() {
        fs::remove_dir_all(&p).map_err(|e| format!("rm {:?}: {}", p, e))?;
    }
    Ok(())
}
