use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

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
    let client = reqwest::Client::builder()
        .cookie_store(true)
        .use_rustls_tls()
        .build()
        .map_err(|e| e.to_string())?;

    // GET /login/ — grab csrftoken from Set-Cookie and csrfmiddlewaretoken from form
    let login_url = format!("{}/login/", UYSD_BASE);
    let get_resp = client.get(&login_url).send().await.map_err(|e| e.to_string())?;
    let csrf_cookie = get_resp
        .cookies()
        .find(|c| c.name() == "csrftoken")
        .map(|c| c.value().to_string())
        .unwrap_or_default();
    let body = get_resp.text().await.map_err(|e| e.to_string())?;
    let csrf_token = extract_csrf_token(&body).unwrap_or_else(|| csrf_cookie.clone());

    // POST /login/
    let post_resp = client
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

    // Success = redirect to /account/
    if post_resp.url().path().contains("account") || post_resp.status().is_success() {
        // Collect session cookies into a compact token string for the frontend
        let session = post_resp
            .cookies()
            .map(|c| format!("{}={}", c.name(), c.value()))
            .collect::<Vec<_>>()
            .join(";");
        if session.contains("sessionid") {
            return Ok(session);
        }
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

/// Poll a submission_result URL. Returns "pending", "ok:<html>", or "error:<msg>".
#[command]
pub async fn uysd_poll_result(result_url: String) -> Result<String, String> {
    let client = reqwest::Client::builder()
        .use_rustls_tls()
        .build()
        .map_err(|e| e.to_string())?;
    let resp = client.get(&result_url).send().await.map_err(|e| e.to_string())?;
    let body = resp.text().await.map_err(|e| e.to_string())?;
    if body.contains("pending") || body.contains("Processing") {
        return Ok("pending".into());
    }
    Ok(format!("ok:{body}"))
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
#[command]
pub fn build_yleaf_zip(output_dir: String) -> Result<String, String> {
    let base = PathBuf::from(&output_dir);
    let zip_path = std::env::temp_dir().join("yleaf_submission.zip");
    let zip_file = File::create(&zip_path).map_err(|e| e.to_string())?;
    let mut zip = zip::ZipWriter::new(zip_file);
    let options = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);

    // hg_prediction*.hg files at root of output_dir
    for entry in std::fs::read_dir(&base).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if name_str.starts_with("hg_prediction") && name_str.ends_with(".hg") {
            zip.start_file(name_str.as_ref(), options).map_err(|e| e.to_string())?;
            let data = std::fs::read(entry.path()).map_err(|e| e.to_string())?;
            zip.write_all(&data).map_err(|e| e.to_string())?;
        }
    }

    // <sample>/<sample>.out files in subdirs
    for entry in std::fs::read_dir(&base).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        if !entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
            continue;
        }
        let sample_dir = entry.path();
        let sample_name = entry.file_name();
        let sample_str = sample_name.to_string_lossy();
        for inner in std::fs::read_dir(&sample_dir).map_err(|e| e.to_string())? {
            let inner = inner.map_err(|e| e.to_string())?;
            let iname = inner.file_name();
            let istr = iname.to_string_lossy();
            if istr.ends_with(".out") {
                let zip_entry = format!("{}/{}", sample_str, istr);
                zip.start_file(&zip_entry, options).map_err(|e| e.to_string())?;
                let data = std::fs::read(inner.path()).map_err(|e| e.to_string())?;
                zip.write_all(&data).map_err(|e| e.to_string())?;
            }
        }
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
