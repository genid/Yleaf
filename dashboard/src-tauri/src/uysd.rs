use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use reqwest::cookie::CookieStore;
use tauri::{command, AppHandle, Manager, State, WebviewUrl, WebviewWindowBuilder};

use crate::{db::SampleMeta, DbState};

#[command]
pub fn open_uysd_window(app: AppHandle, url: String) {
    let parsed = match url.parse() {
        Ok(u) => u,
        Err(_) => return,
    };
    // Reuse an existing UYSD window if open; otherwise create one.
    if let Some(win) = app.get_webview_window("uysd") {
        let _ = win.navigate(parsed);
        let _ = win.set_focus();
    } else {
        let _ = WebviewWindowBuilder::new(&app, "uysd", WebviewUrl::External(parsed))
            .title("UYSD — Y-SNP Distribution")
            .inner_size(1280.0, 900.0)
            .build();
    }
}

#[command]
pub fn get_uysd_embedded(db: State<DbState>) -> bool {
    let conn = db.0.lock().unwrap();
    conn.query_row(
        "SELECT value FROM uysd_settings WHERE key='embedded_map'",
        [],
        |row| row.get::<_, i64>(0),
    )
    .unwrap_or(0)
        != 0
}

#[command]
pub fn set_uysd_embedded(db: State<DbState>, enabled: bool) {
    let conn = db.0.lock().unwrap();
    let _ = conn.execute(
        "INSERT INTO uysd_settings (key, value) VALUES ('embedded_map', ?1)
         ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        [enabled as i64],
    );
}

#[command]
pub fn get_sample_meta(db: State<DbState>, job_id: i64) -> Vec<SampleMeta> {
    let conn = db.0.lock().unwrap();
    crate::db::get_sample_meta(&conn, job_id).unwrap_or_default()
}

#[command]
pub fn upsert_sample_meta(db: State<DbState>, meta: SampleMeta) {
    let conn = db.0.lock().unwrap();
    let _ = crate::db::upsert_sample_meta(&conn, &meta);
}

// ── HTTP: login / submit / poll ──────────────────────────────────────────────

const UYSD_BASE: &str = "https://ysnp.erasmusmc.nl";

/// Log in to UYSD. Returns a session token string (sessionid=…;csrftoken=…)
/// that the frontend passes back to submit/poll, or an error message.
#[command]
pub async fn uysd_login(username: String, password: String) -> Result<String, String> {
    // Explicit cookie jar so we can read cookies set during the redirect chain
    // (sessionid is set on the 302 after login; reqwest::Response::cookies() on
    // the final landing page returns nothing useful).
    let jar = std::sync::Arc::new(reqwest::cookie::Jar::default());
    let client = reqwest::Client::builder()
        .cookie_provider(jar.clone())
        .use_rustls_tls()
        .build()
        .map_err(|e| e.to_string())?;

    let base = reqwest::Url::parse(UYSD_BASE).map_err(|e| e.to_string())?;
    let login_url = format!("{}/login/", UYSD_BASE);

    // GET /login/ — sets csrftoken cookie in the jar; we also need the form's
    // csrfmiddlewaretoken value (Django checks both).
    let get_resp = client.get(&login_url).send().await.map_err(|e| e.to_string())?;
    let body = get_resp.text().await.map_err(|e| e.to_string())?;
    let csrf_token = extract_csrf_token(&body)
        .ok_or_else(|| "Could not find csrfmiddlewaretoken on /login/ page".to_string())?;

    // POST /login/
    let _post_resp = client
        .post(&login_url)
        .header("Referer", &login_url)
        .form(&[
            ("username", username.as_str()),
            ("password", password.as_str()),
            ("csrfmiddlewaretoken", csrf_token.as_str()),
        ])
        .send()
        .await
        .map_err(|e| e.to_string())?;

    // Authoritative signal: sessionid present in the jar for this domain.
    let session = jar
        .cookies(&base)
        .map(|h| h.to_str().unwrap_or("").to_string())
        .unwrap_or_default();
    if session.contains("sessionid") {
        return Ok(session);
    }
    Err("Invalid credentials or login failed".into())
}

/// Submit to UYSD. Returns the submission result key URL.
#[command]
pub async fn uysd_submit(
    session_token: String,
    country_csv_path: String,
    yleaf_zip_path: String,
) -> Result<String, String> {
    let client = reqwest::Client::builder()
        .cookie_store(true)
        .use_rustls_tls()
        .build()
        .map_err(|e| e.to_string())?;

    // Parse session_token back into cookies, extract csrftoken
    let mut csrf = String::new();
    for part in session_token.split(';') {
        let part = part.trim();
        if let Some(v) = part.strip_prefix("csrftoken=") {
            csrf = v.to_string();
        }
    }

    let csv_bytes = std::fs::read(&country_csv_path).map_err(|e| e.to_string())?;
    let zip_bytes = std::fs::read(&yleaf_zip_path).map_err(|e| e.to_string())?;

    let submit_url = format!("{}/submission/", UYSD_BASE);
    let form = reqwest::multipart::Form::new()
        .part("country_file", reqwest::multipart::Part::bytes(csv_bytes).file_name("country_file.csv").mime_str("text/csv").map_err(|e| e.to_string())?)
        .part("yleaf_zip", reqwest::multipart::Part::bytes(zip_bytes).file_name("yleaf.zip").mime_str("application/zip").map_err(|e| e.to_string())?)
        .text("terms_conditions", "on")
        .text("yleaf_submit", "1")
        .text("csrfmiddlewaretoken", csrf.clone());

    let resp = client
        .post(&submit_url)
        .header("Referer", &submit_url)
        .header("Cookie", &session_token)
        .header("X-CSRFToken", &csrf)
        .multipart(form)
        .send()
        .await
        .map_err(|e| e.to_string())?;

    let final_url = resp.url().to_string();
    if final_url.contains("submission_result") {
        return Ok(final_url);
    }
    Err(format!("Submission failed — landed at: {final_url}"))
}

/// Poll a submission_result URL using the authenticated session.
///
/// Returns:
///   Ok("pending")          — submission is still being processed.
///   Ok("ok")               — submission succeeded and is in the database.
///   Err("<message>")       — submission failed (auth, server, validation) — message is the reason.
///
/// Note: the result page requires authentication.  Without the session cookie
/// UYSD redirects to its login page, which historically caused the poller to
/// misclassify every submission as success.  Always pass the same session
/// token that uysd_login returned.
#[command]
pub async fn uysd_poll_result(
    result_url: String,
    session_token: String,
) -> Result<String, String> {
    let client = reqwest::Client::builder()
        .use_rustls_tls()
        .build()
        .map_err(|e| e.to_string())?;
    let resp = client
        .get(&result_url)
        .header("Cookie", &session_token)
        .send()
        .await
        .map_err(|e| e.to_string())?;
    let body = resp.text().await.map_err(|e| e.to_string())?;

    // Authentication was rejected — UYSD redirected us to the login page.
    if body.contains("<h1>Login</h1>") || body.contains("\"id_username\"") {
        return Err("Session expired or not authenticated. Please log in again.".into());
    }

    // Explicit failure messages observed in UYSD's submission_result page.
    if body.contains("Server failure")
        || body.contains("Server encountered an error")
        || body.contains("could not be processed")
    {
        // Try to extract the human message inside the result-card area.
        let msg = extract_failure_message(&body)
            .unwrap_or_else(|| "Server failure — UYSD reports the submission could not be processed.".to_string());
        return Err(msg);
    }

    // Still processing.
    if body.contains("pending") || body.contains("Processing") {
        return Ok("pending".into());
    }

    // Any other shape we treat as success.  TODO: once UYSD exposes a
    // machine-readable `data-status` / status JSON, switch to that.
    Ok("ok".into())
}

fn extract_failure_message(html: &str) -> Option<String> {
    // Look for the first heading + paragraph in the content area.
    // Typical shape: <h2>Server failure</h2><p>… message …</p>
    let re_blocks = regex_lite_extract(html, "<h[1-6][^>]*>", "</h");
    let title = re_blocks.first().cloned();
    let para = regex_lite_extract(html, "<p[^>]*>", "</p>")
        .into_iter()
        .find(|p| {
            let s = p.to_lowercase();
            s.contains("server") || s.contains("error") || s.contains("could not") || s.contains("inconvenience")
        });
    match (title, para) {
        (Some(t), Some(p)) => Some(format!("{t} — {p}")),
        (Some(t), None) => Some(t),
        (None, Some(p)) => Some(p),
        _ => None,
    }
}

/// Very small "extract text between known delimiters" helper — avoids pulling
/// in a full regex crate for one-off content scraping.  Returns the texts
/// between consecutive matches of `open` and `close`, with HTML tags stripped.
fn regex_lite_extract(html: &str, open: &str, close: &str) -> Vec<String> {
    // open is something like "<h2>" with a wildcard for attributes — we just
    // anchor on the first 3 chars ("<h2", "<p")
    let anchor = &open[..open.len().min(3)];
    let mut out = Vec::new();
    let mut cur = html;
    while let Some(start) = cur.find(anchor) {
        // skip until '>'
        let after_tag = match cur[start..].find('>') {
            Some(p) => start + p + 1,
            None => break,
        };
        let end_close = match cur[after_tag..].find(close) {
            Some(p) => after_tag + p,
            None => break,
        };
        let raw = &cur[after_tag..end_close];
        // strip nested tags
        let mut text = String::new();
        let mut in_tag = false;
        for c in raw.chars() {
            match c {
                '<' => in_tag = true,
                '>' => in_tag = false,
                _ if !in_tag => text.push(c),
                _ => {}
            }
        }
        let text = text.trim().to_string();
        if !text.is_empty() {
            out.push(text);
        }
        cur = &cur[end_close + close.len()..];
    }
    out
}

fn base_haplogroup(hg: &str) -> &str {
    // Strip "*(…)" wildcard suffix; "E-Z15929*(xE-Y25504)" -> "E-Z15929"
    hg.split('*').next().unwrap_or(hg)
}

fn page_has_frequencies(html: &str) -> bool {
    // The frequency-data array is inlined as:
    //   var frequencies = JSON.parse("[{...");   (has data)
    //   var frequencies = JSON.parse("[]");       (empty -- valid hg, no samples)
    //   (line absent)                             (error page or unknown hg)
    //
    // Other JSON.parse calls in the same page (slug_map, small_tree, …) must
    // NOT be matched -- match the literal variable assignment only.
    html.contains(r#"var frequencies = JSON.parse("[{"#)
}

fn percent_encode_segment(s: &str) -> String {
    // RFC 3986 unreserved set; everything else gets percent-encoded.
    // Url::path_segments_mut().push() leaves sub-delims (`*`, `(`, `)`, `,`)
    // unencoded, which WebKit then refuses to load in an iframe src.
    let mut out = String::with_capacity(s.len() * 3);
    for byte in s.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => {
                out.push(byte as char);
            }
            _ => out.push_str(&format!("%{:02X}", byte)),
        }
    }
    out
}

fn full_url_for(haplogroup: &str) -> String {
    format!("{UYSD_BASE}/haplogroup/{}", percent_encode_segment(haplogroup))
}

fn embed_url_for(haplogroup: &str) -> String {
    // ?embed=1 triggers UYSD's minimal layout (chrome hidden, map+tree-nav only)
    format!("{}?embed=1", full_url_for(haplogroup))
}

#[derive(serde::Serialize)]
pub struct UysdMapUrls {
    /// Bare URL — opened externally when the user clicks the preview.
    pub full_url: String,
    /// Minimal-layout URL used as the iframe preview src.
    pub embed_url: String,
}

#[derive(serde::Serialize, serde::Deserialize)]
pub struct KnownLocations {
    pub countries: Vec<String>,
    pub regions: Vec<String>,
}

/// Fetch the list of country/region names UYSD accepts for submission.
/// Caching is handled on the frontend (one fetch per session is enough; the
/// list rarely changes and UYSD ships a `Cache-Control: max-age=86400` header).
#[command]
pub async fn uysd_get_known_locations() -> Result<KnownLocations, String> {
    let client = reqwest::Client::builder()
        .use_rustls_tls()
        .build()
        .map_err(|e| e.to_string())?;
    let resp = client
        .get(format!("{UYSD_BASE}/api/known_locations/"))
        .send()
        .await
        .map_err(|e| e.to_string())?;
    if !resp.status().is_success() {
        return Err(format!("UYSD returned HTTP {} for /api/known_locations/", resp.status()));
    }
    let body = resp.text().await.map_err(|e| e.to_string())?;
    serde_json::from_str::<KnownLocations>(&body).map_err(|e| e.to_string())
}

/// Resolve which UYSD haplogroup-map URL to load.  Tries the full wildcard form
/// first; if UYSD reports no frequency data, falls back to the base haplogroup.
#[command]
pub async fn uysd_resolve_map_url(haplogroup: String) -> UysdMapUrls {
    let base_hg = base_haplogroup(&haplogroup);
    let make = |hg: &str| UysdMapUrls {
        full_url: full_url_for(hg),
        embed_url: embed_url_for(hg),
    };

    if base_hg == haplogroup {
        return make(&haplogroup); // nothing to disambiguate
    }

    let client = match reqwest::Client::builder().use_rustls_tls().build() {
        Ok(c) => c,
        Err(_) => return make(base_hg),
    };
    // Probe the full wildcard URL (without ?embed=1; the frequency JSON is the
    // same either way and the bare URL keeps things simple).
    let resp = match client.get(full_url_for(&haplogroup)).send().await {
        Ok(r) => r,
        Err(_) => return make(base_hg),
    };
    let body = match resp.text().await {
        Ok(b) => b,
        Err(_) => return make(base_hg),
    };
    if page_has_frequencies(&body) {
        make(&haplogroup)
    } else {
        make(base_hg)
    }
}

fn extract_csrf_token(html: &str) -> Option<String> {
    // <input type="hidden" name="csrfmiddlewaretoken" value="...">
    let needle = "name=\"csrfmiddlewaretoken\" value=\"";
    let start = html.find(needle)? + needle.len();
    let end = html[start..].find('"')? + start;
    Some(html[start..end].to_string())
}

// ── CSV / ZIP builders ───────────────────────────────────────────────────────

#[derive(serde::Deserialize)]
pub struct CsvRow {
    pub sample_name: String,
    pub country: String,
    pub region: String,
    pub comment: String,
    pub publication: String,
}

fn csv_escape(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "\"\""))
}

/// Write rows to a temp CSV file and return its path.
#[command]
pub fn build_country_csv(rows: Vec<CsvRow>) -> Result<String, String> {
    let path = std::env::temp_dir().join("yleaf_country_file.csv");
    let mut f = File::create(&path).map_err(|e| e.to_string())?;
    for row in &rows {
        let line = format!(
            "{},{},{},{},{}\r\n",
            csv_escape(&row.sample_name),
            csv_escape(&row.country),
            csv_escape(&row.region),
            csv_escape(&row.comment),
            csv_escape(&row.publication),
        );
        f.write_all(line.as_bytes()).map_err(|e| e.to_string())?;
    }
    Ok(path.to_string_lossy().into_owned())
}

/// Zip every <sample>/<sample>.out and hg_prediction*.hg under output_dir.
/// Returns the path to the created archive.
/// Build the submission zip in the legacy single-database shape that UYSD's
/// pipeline expects:  one `hg_prediction.hg` and one `<sample>/<sample>.out`.
///
/// Yleaf v4 multi-database runs produce per-tree files (`hg_prediction_<tree>.hg`,
/// `<sample>.<tree>.out`) plus a combined `hg_prediction_combined.hg` with 13
/// columns — none of which UYSD knows how to parse.  UYSD's haplogroup map
/// is YFull-only anyway, so we pick yfull and rename to the legacy shape.
#[command]
pub fn build_yleaf_zip(output_dir: String) -> Result<String, String> {
    let base = PathBuf::from(&output_dir);
    let zip_path = std::env::temp_dir().join("yleaf_submission.zip");
    let zip_file = File::create(&zip_path).map_err(|e| e.to_string())?;
    let mut zip = zip::ZipWriter::new(zip_file);
    let options = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);

    // Pick exactly one .hg file at the root: prefer hg_prediction_yfull.hg
    // (v4 multi-tree), fall back to hg_prediction.hg (v3 / single-tree yfull),
    // and rename either to `hg_prediction.hg` in the zip.
    let mut chose_hg: Option<PathBuf> = None;
    let mut fallback_hg: Option<PathBuf> = None;
    for entry in std::fs::read_dir(&base).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name();
        let n = name.to_string_lossy();
        if n == "hg_prediction_yfull.hg" {
            chose_hg = Some(entry.path());
            break;
        }
        if n == "hg_prediction.hg" {
            fallback_hg = Some(entry.path());
        }
    }
    let hg_src = chose_hg.or(fallback_hg)
        .ok_or_else(|| "No yfull hg_prediction file found in output dir".to_string())?;
    zip.start_file("hg_prediction.hg", options).map_err(|e| e.to_string())?;
    let data = std::fs::read(&hg_src).map_err(|e| e.to_string())?;
    zip.write_all(&data).map_err(|e| e.to_string())?;

    // For each sample subdir, pick one .out: prefer <sample>.yfull.out,
    // fall back to <sample>.out, and store it as `<sample>/<sample>.out`.
    for entry in std::fs::read_dir(&base).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        if !entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
            continue;
        }
        let sample_dir = entry.path();
        let sample_name = entry.file_name();
        let sample_str = sample_name.to_string_lossy().into_owned();

        let yfull_out = sample_dir.join(format!("{sample_str}.yfull.out"));
        let plain_out = sample_dir.join(format!("{sample_str}.out"));
        let src = if yfull_out.exists() {
            yfull_out
        } else if plain_out.exists() {
            plain_out
        } else {
            continue; // sample without a usable yfull/plain .out — skip
        };

        let zip_entry = format!("{sample_str}/{sample_str}.out");
        zip.start_file(&zip_entry, options).map_err(|e| e.to_string())?;
        let data = std::fs::read(&src).map_err(|e| e.to_string())?;
        zip.write_all(&data).map_err(|e| e.to_string())?;
    }

    zip.finish().map_err(|e| e.to_string())?;
    Ok(zip_path.to_string_lossy().into_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_escape_quotes() {
        assert_eq!(csv_escape("he said \"hi\""), "\"he said \\\"hi\\\"\"".replace("\\\"", "\"\""));
        assert_eq!(csv_escape(r#"a "b" c"#), r#""a ""b"" c""#);
    }

    #[test]
    fn build_country_csv_roundtrip() {
        let rows = vec![
            CsvRow {
                sample_name: "S1".into(),
                country: "Netherlands".into(),
                region: "South Holland".into(),
                comment: "test".into(),
                publication: "https://doi.org/10.1000/xyz".into(),
            },
            CsvRow {
                sample_name: "S2".into(),
                country: "Germany".into(),
                region: "".into(),
                comment: "".into(),
                publication: "".into(),
            },
        ];
        let path = build_country_csv(rows).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("\"S1\",\"Netherlands\",\"South Holland\""));
        assert!(content.contains("\"S2\",\"Germany\""));
        assert!(content.contains("\r\n"));
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn build_yleaf_zip_collects_files() {
        let dir = tempfile::tempdir().unwrap();
        let base = dir.path();
        // hg_prediction file
        std::fs::write(base.join("hg_prediction_yfull.hg"), "hg_data").unwrap();
        // sample subdir with .out file
        let sample_dir = base.join("Sample_1");
        std::fs::create_dir(&sample_dir).unwrap();
        std::fs::write(sample_dir.join("Sample_1.out"), "out_data").unwrap();

        let zip_path = build_yleaf_zip(base.to_string_lossy().into_owned()).unwrap();
        let zip_file = File::open(&zip_path).unwrap();
        let mut archive = zip::ZipArchive::new(zip_file).unwrap();
        let names: Vec<String> = (0..archive.len())
            .map(|i| archive.by_index(i).unwrap().name().to_string())
            .collect();
        assert!(names.iter().any(|n| n == "hg_prediction_yfull.hg"), "missing hg file: {names:?}");
        assert!(names.iter().any(|n| n == "Sample_1/Sample_1.out"), "missing .out file: {names:?}");
        std::fs::remove_file(zip_path).ok();
    }
}
