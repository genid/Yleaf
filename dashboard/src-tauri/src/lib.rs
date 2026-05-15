mod db;
mod jobs;
mod results;

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use tauri::Manager;

pub struct DbState(pub Arc<Mutex<rusqlite::Connection>>);

pub struct ActiveChild(pub Arc<Mutex<HashMap<i64, (String, tauri_plugin_shell::process::CommandChild)>>>);

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            let icon = tauri::image::Image::from_bytes(include_bytes!("../icons/128x128.png"))
                .expect("failed to load icon");
            app.get_webview_window("main")
                .unwrap()
                .set_icon(icon)
                .expect("failed to set window icon");
            let db_path = app
                .path()
                .app_data_dir()
                .expect("no app data dir")
                .join("jobs.db");
            std::fs::create_dir_all(db_path.parent().unwrap())?;
            let conn = rusqlite::Connection::open(&db_path)
                .expect("failed to open SQLite DB");
            db::init_schema(&conn).expect("failed to init DB schema");
            app.manage(DbState(Arc::new(Mutex::new(conn))));
            app.manage(ActiveChild(Arc::new(Mutex::new(HashMap::new()))));
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            jobs::run_yleaf,
            jobs::kill_job,
            results::list_jobs,
            results::get_job_results,
            results::get_haplogroup_path,
            results::delete_job,
            results::rename_job,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
