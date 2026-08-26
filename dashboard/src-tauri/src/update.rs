//! Startup update check for the dashboard.
//!
//! The app has no self-updater, and Yleaf is distributed from GitHub rather
//! than an app store, so this only reports that a newer release exists and
//! links to the releases page — installing stays the user's deliberate step.
//!
//! Best-effort by design: short timeout, every failure returns `None`. An
//! offline machine must behave exactly as it did before.

use serde::Serialize;

const RELEASES_API_URL: &str = "https://api.github.com/repos/genid/Yleaf/releases/latest";
const RELEASES_PAGE_URL: &str = "https://github.com/genid/Yleaf/releases/latest";
const REQUEST_TIMEOUT_SECS: u64 = 3;

#[derive(Serialize, Clone)]
pub struct UpdateInfo {
    pub current: String,
    pub latest: String,
    pub url: String,
}

#[derive(serde::Deserialize)]
struct Release {
    tag_name: String,
}

/// Comparable key for a dotted version. Trailing non-numeric parts (rc, beta)
/// are ignored, so 4.2.0-rc1 never compares newer than 4.2.0.
fn version_key(version: &str) -> Vec<u32> {
    let mut parts = Vec::new();
    for part in version.split('.') {
        let digits: String = part.chars().take_while(|c| c.is_ascii_digit()).collect();
        match digits.parse::<u32>() {
            Ok(number) => parts.push(number),
            Err(_) => break,
        }
    }
    parts
}

fn is_newer(latest: &str, current: &str) -> bool {
    let (latest_key, current_key) = (version_key(latest), version_key(current));
    if latest_key.is_empty() || current_key.is_empty() {
        return false;
    }
    latest_key > current_key
}

/// Newer release than the running app, or None when up to date or unreachable.
#[tauri::command]
pub async fn check_for_update(app: tauri::AppHandle) -> Option<UpdateInfo> {
    let current = app.package_info().version.to_string();

    let client = reqwest::Client::builder()
        .use_rustls_tls()
        .timeout(std::time::Duration::from_secs(REQUEST_TIMEOUT_SECS))
        .build()
        .ok()?;

    // GitHub rejects requests without a User-Agent.
    let response = client
        .get(RELEASES_API_URL)
        .header("User-Agent", "Yleaf-Dashboard")
        .header("Accept", "application/vnd.github+json")
        .send()
        .await
        .map_err(|e| eprintln!("update check skipped: {e}"))
        .ok()?;

    // reqwest is built without its `json` feature here, so parse the body
    // ourselves — the same pattern uysd.rs uses.
    let body = response
        .text()
        .await
        .map_err(|e| eprintln!("update check skipped: {e}"))
        .ok()?;
    let release: Release = serde_json::from_str(&body)
        .map_err(|e| eprintln!("update check skipped: {e}"))
        .ok()?;

    let latest = release.tag_name.trim_start_matches(['v', 'V']).to_string();
    if !is_newer(&latest, &current) {
        return None;
    }
    Some(UpdateInfo { current, latest, url: RELEASES_PAGE_URL.to_string() })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newer_versions_are_detected() {
        assert!(is_newer("4.2.0", "4.1.4"));
        assert!(is_newer("5.0", "4.1.4"));
        // numeric compare, not lexicographic
        assert!(is_newer("4.10.0", "4.9.9"));
    }

    #[test]
    fn same_or_older_versions_are_not_offered() {
        assert!(!is_newer("4.1.4", "4.1.4"));
        assert!(!is_newer("4.1.3", "4.1.4"));
        assert!(!is_newer("4.9.9", "4.10.0"));
    }

    #[test]
    fn unparseable_versions_never_offer_an_update() {
        assert!(!is_newer("", "4.1.4"));
        assert!(!is_newer("4.1.4", ""));
        assert!(!is_newer("nightly", "4.1.4"));
        // a pre-release of the version we already run is not an upgrade
        assert!(!is_newer("4.2.0-rc1", "4.2.0"));
    }
}
