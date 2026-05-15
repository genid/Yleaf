use rusqlite::{params, Connection, Result as SqlResult};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Job {
    pub id: i64,
    pub sample_name: String,
    pub sample_path: String,
    pub output_dir: String,
    pub trees: String,
    pub status: String,
    pub start_ts: i64,
    pub end_ts: Option<i64>,
    pub haplogroup: Option<String>,
    pub qc_score: Option<f64>,
    pub qc1: Option<f64>,
    pub qc2: Option<f64>,
    pub qc3: Option<f64>,
    pub total_reads: Option<i64>,
    pub valid_markers: Option<i64>,
    // run settings (for reload)
    pub reference_genome: String,
    pub threads: i64,
    pub reads_threshold: i64,
    pub quality_thresh: i64,
    pub base_majority: i64,
    pub prediction_quality: f64,
    pub draw_haplogroups: bool,
    pub ancient_dna: bool,
    pub private_mutations: bool,
    pub collapsed_draw_mode: bool,
}

pub fn init_schema(conn: &Connection) -> SqlResult<()> {
    conn.execute_batch(
        "CREATE TABLE IF NOT EXISTS jobs (
            id                 INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_name        TEXT    NOT NULL,
            sample_path        TEXT    NOT NULL,
            output_dir         TEXT    NOT NULL,
            trees              TEXT    NOT NULL DEFAULT '',
            status             TEXT    NOT NULL DEFAULT 'running',
            start_ts           INTEGER NOT NULL,
            end_ts             INTEGER,
            haplogroup         TEXT,
            qc_score           REAL,
            qc1                REAL,
            qc2                REAL,
            qc3                REAL,
            total_reads        INTEGER,
            valid_markers      INTEGER,
            reference_genome   TEXT    NOT NULL DEFAULT 'hg38',
            threads            INTEGER NOT NULL DEFAULT 4,
            reads_threshold    INTEGER NOT NULL DEFAULT 10,
            quality_thresh     INTEGER NOT NULL DEFAULT 20,
            base_majority      INTEGER NOT NULL DEFAULT 90,
            prediction_quality REAL    NOT NULL DEFAULT 0.95,
            draw_haplogroups   INTEGER NOT NULL DEFAULT 0,
            ancient_dna        INTEGER NOT NULL DEFAULT 0,
            private_mutations  INTEGER NOT NULL DEFAULT 0,
            collapsed_draw_mode INTEGER NOT NULL DEFAULT 0
        );",
    )?;
    // Migrate existing DBs — ignore errors if columns already exist
    for col in &[
        "ALTER TABLE jobs ADD COLUMN reference_genome   TEXT    NOT NULL DEFAULT 'hg38'",
        "ALTER TABLE jobs ADD COLUMN threads            INTEGER NOT NULL DEFAULT 4",
        "ALTER TABLE jobs ADD COLUMN reads_threshold    INTEGER NOT NULL DEFAULT 10",
        "ALTER TABLE jobs ADD COLUMN quality_thresh     INTEGER NOT NULL DEFAULT 20",
        "ALTER TABLE jobs ADD COLUMN base_majority      INTEGER NOT NULL DEFAULT 90",
        "ALTER TABLE jobs ADD COLUMN prediction_quality REAL    NOT NULL DEFAULT 0.95",
        "ALTER TABLE jobs ADD COLUMN draw_haplogroups    INTEGER NOT NULL DEFAULT 0",
        "ALTER TABLE jobs ADD COLUMN ancient_dna         INTEGER NOT NULL DEFAULT 0",
        "ALTER TABLE jobs ADD COLUMN private_mutations   INTEGER NOT NULL DEFAULT 0",
        "ALTER TABLE jobs ADD COLUMN collapsed_draw_mode INTEGER NOT NULL DEFAULT 0",
    ] {
        let _ = conn.execute_batch(col);
    }
    Ok(())
}

pub fn insert_job(
    conn: &Connection,
    sample_name: &str,
    sample_path: &str,
    output_dir: &str,
    trees: &str,
    start_ts: i64,
    reference_genome: &str,
    threads: i64,
    reads_threshold: i64,
    quality_thresh: i64,
    base_majority: i64,
    prediction_quality: f64,
    draw_haplogroups: bool,
    ancient_dna: bool,
    private_mutations: bool,
    collapsed_draw_mode: bool,
) -> SqlResult<i64> {
    conn.execute(
        "INSERT INTO jobs (sample_name, sample_path, output_dir, trees, status, start_ts,
                           reference_genome, threads, reads_threshold, quality_thresh,
                           base_majority, prediction_quality, draw_haplogroups, ancient_dna,
                           private_mutations, collapsed_draw_mode)
         VALUES (?1,?2,?3,?4,'running',?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15)",
        params![
            sample_name, sample_path, output_dir, trees, start_ts,
            reference_genome, threads, reads_threshold, quality_thresh,
            base_majority, prediction_quality, draw_haplogroups as i64, ancient_dna as i64,
            private_mutations as i64, collapsed_draw_mode as i64
        ],
    )?;
    Ok(conn.last_insert_rowid())
}

pub fn update_job_done(
    conn: &Connection,
    id: i64,
    haplogroup: &str,
    qc_score: f64,
    qc1: f64,
    qc2: f64,
    qc3: f64,
    total_reads: i64,
    valid_markers: i64,
    end_ts: i64,
) -> SqlResult<()> {
    conn.execute(
        "UPDATE jobs SET status='done', haplogroup=?1, qc_score=?2, qc1=?3, qc2=?4,
         qc3=?5, total_reads=?6, valid_markers=?7, end_ts=?8 WHERE id=?9",
        params![haplogroup, qc_score, qc1, qc2, qc3, total_reads, valid_markers, end_ts, id],
    )?;
    Ok(())
}

pub fn update_job_error(conn: &Connection, id: i64, end_ts: i64) -> SqlResult<()> {
    conn.execute(
        "UPDATE jobs SET status='error', end_ts=?1 WHERE id=?2",
        params![end_ts, id],
    )?;
    Ok(())
}

pub fn delete_job(conn: &Connection, id: i64) -> SqlResult<()> {
    conn.execute("DELETE FROM jobs WHERE id=?1", params![id])?;
    Ok(())
}

pub fn rename_job(conn: &Connection, id: i64, name: &str) -> SqlResult<()> {
    conn.execute("UPDATE jobs SET sample_name=?1 WHERE id=?2", params![name, id])?;
    Ok(())
}

pub fn list_jobs(conn: &Connection) -> SqlResult<Vec<Job>> {
    let mut stmt = conn.prepare(
        "SELECT id, sample_name, sample_path, output_dir, trees, status,
                start_ts, end_ts, haplogroup, qc_score, qc1, qc2, qc3,
                total_reads, valid_markers,
                reference_genome, threads, reads_threshold, quality_thresh,
                base_majority, prediction_quality, draw_haplogroups, ancient_dna,
                private_mutations, collapsed_draw_mode
         FROM jobs ORDER BY start_ts DESC",
    )?;
    let rows = stmt.query_map([], |row| {
        Ok(Job {
            id: row.get(0)?,
            sample_name: row.get(1)?,
            sample_path: row.get(2)?,
            output_dir: row.get(3)?,
            trees: row.get(4)?,
            status: row.get(5)?,
            start_ts: row.get(6)?,
            end_ts: row.get(7)?,
            haplogroup: row.get(8)?,
            qc_score: row.get(9)?,
            qc1: row.get(10)?,
            qc2: row.get(11)?,
            qc3: row.get(12)?,
            total_reads: row.get(13)?,
            valid_markers: row.get(14)?,
            reference_genome: row.get::<_, Option<String>>(15)?.unwrap_or_else(|| "hg38".into()),
            threads: row.get::<_, Option<i64>>(16)?.unwrap_or(4),
            reads_threshold: row.get::<_, Option<i64>>(17)?.unwrap_or(10),
            quality_thresh: row.get::<_, Option<i64>>(18)?.unwrap_or(20),
            base_majority: row.get::<_, Option<i64>>(19)?.unwrap_or(90),
            prediction_quality: row.get::<_, Option<f64>>(20)?.unwrap_or(0.95),
            draw_haplogroups: row.get::<_, Option<i64>>(21)?.unwrap_or(0) != 0,
            ancient_dna: row.get::<_, Option<i64>>(22)?.unwrap_or(0) != 0,
            private_mutations: row.get::<_, Option<i64>>(23)?.unwrap_or(0) != 0,
            collapsed_draw_mode: row.get::<_, Option<i64>>(24)?.unwrap_or(0) != 0,
        })
    })?;
    rows.collect()
}
