import sqlite3
import pandas as pd
from datetime import datetime


class ResultsDB:
    def __init__(self, db_path="data/processed/results.db"):
        self.db_path = db_path
        self._init_db()

    def _init_db(self):
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS experiments (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    run_date TEXT NOT NULL,
                    sample_size INTEGER NOT NULL
                )
            """)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS results (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    experiment_id INTEGER NOT NULL,
                    protein TEXT NOT NULL,
                    raw_count INTEGER NOT NULL,
                    harmonized_count INTEGER NOT NULL,
                    ic_raw REAL NOT NULL,
                    ic_harmonized REAL NOT NULL,
                    improvement REAL NOT NULL,
                    FOREIGN KEY (experiment_id) REFERENCES experiments(id)
                )
            """)

    def save_results(self, df: pd.DataFrame, sample_size: int) -> int:
        # Support both old column name (Improvement) and new (Improvement_Pct)
        imp_col = 'Improvement_Pct' if 'Improvement_Pct' in df.columns else 'Improvement'

        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "INSERT INTO experiments (run_date, sample_size) VALUES (?, ?)",
                (datetime.now().isoformat(), sample_size),
            )
            experiment_id = cursor.lastrowid

            for _, row in df.iterrows():
                conn.execute(
                    """INSERT INTO results
                       (experiment_id, protein, raw_count, harmonized_count,
                        ic_raw, ic_harmonized, improvement)
                       VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    (
                        experiment_id,
                        row["Protein"],
                        int(row["Raw_Count"]),
                        int(row["Harmonized_Count"]),
                        float(row["IC_Raw"]),
                        float(row["IC_Harmonized"]),
                        float(row[imp_col]),
                    ),
                )
        return experiment_id

    def get_experiment(self, experiment_id: int) -> pd.DataFrame:
        with sqlite3.connect(self.db_path) as conn:
            return pd.read_sql_query(
                "SELECT * FROM results WHERE experiment_id = ?",
                conn,
                params=(experiment_id,),
            )

    def list_experiments(self) -> pd.DataFrame:
        with sqlite3.connect(self.db_path) as conn:
            return pd.read_sql_query("SELECT * FROM experiments", conn)
