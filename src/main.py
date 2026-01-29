import pandas as pd
from src.data_loader import BioDataLoader
from src.harmonizer import GOHarmonizer
from src.database import ResultsDB

def main():
    # 1. Завантаження
    print("--- Етап 1: Завантаження даних ---")
    loader = BioDataLoader(
        obo_path="data/raw/go-basic.obo",
        gaf_path="data/raw/goa_human_plus.gaf"
    )
    dag = loader.load_ontology()
    loader.load_annotations()
    df = loader.get_df()
    
    print(f"Всього анотацій: {len(df)}")

    # 2. Ініціалізація гармонізатора
    print("\n--- Етап 2: Ініціалізація алгоритму ---")
    harmonizer = GOHarmonizer(dag, df)

    # 3. Експеримент
    print("\n--- Етап 3: Запуск гармонізації (True Path Rule) ---")
    results_table = harmonizer.run_experiment(sample_size=100)
    
    print("\nРезультати експерименту:")
    print(results_table.to_string())

    # Зберегти у CSV для диплому
    results_table.to_csv("data/processed/harmonization_results.csv", index=False)
    print("\nРезультати збережено у data/processed/harmonization_results.csv")

    # Зберегти у SQLite
    db = ResultsDB()
    exp_id = db.save_results(results_table, sample_size=100)
    print(f"Результати збережено у SQLite (experiment_id={exp_id})")

if __name__ == "__main__":
    main()