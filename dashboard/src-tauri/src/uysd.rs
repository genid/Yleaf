use tauri::{command, AppHandle, Manager, State, WebviewUrl, WebviewWindowBuilder};

use crate::DbState;

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
