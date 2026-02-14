"""
Pre-compute all protein data so the Protein Explorer works on Streamlit Cloud
without the raw GO/GAF files.  Run once locally:

    python precompute_explorer.py

Generates:
    data/processed/all_protein_metrics.csv    (~2 MB)
    data/processed/protein_terms.csv.gz       (~5-10 MB compressed)
    data/processed/go_terms_info.csv          (~500 KB)
"""

import pandas as pd
from src.data_loader import BioDataLoader
from src.harmonizer import GOHarmonizer


def main():
    # 1. Load data
    print("Loading data...")
    loader = BioDataLoader(
        obo_path="data/raw/go-basic.obo",
        gaf_path="data/raw/goa_human_plus.gaf",
    )
    dag = loader.load_ontology()
    loader.load_annotations()
    df = loader.get_df()

    harmonizer = GOHarmonizer(dag, df)

    # 2. GO terms info table (go_id → name, namespace, IC)
    print("Building GO terms info table...")
    term_rows = []
    for go_id in harmonizer.term_counts:
        name = dag[go_id].name if go_id in dag else ""
        ns = harmonizer.get_term_namespace(go_id)
        ic = harmonizer.get_information_content(go_id)
        term_rows.append({"go_id": go_id, "name": name, "namespace": ns, "ic": round(ic, 4)})
    pd.DataFrame(term_rows).to_csv("data/processed/go_terms_info.csv", index=False)
    print(f"  Saved go_terms_info.csv ({len(term_rows)} terms)")

    # 3. Harmonize ALL proteins
    unique_proteins = df["DB_Object_ID"].unique()
    print(f"Harmonizing all {len(unique_proteins)} proteins...")

    metrics_rows = []
    terms_rows = []

    for i, prot in enumerate(unique_proteins):
        if (i + 1) % 5000 == 0:
            print(f"  {i + 1}/{len(unique_proteins)}...")

        raw, harmonized = harmonizer.harmonize_protein(prot)
        ic_raw = harmonizer.evaluate_informativeness(raw)
        ic_harm = harmonizer.evaluate_informativeness(harmonized)
        ns_ic_raw, ns_count_raw = harmonizer.evaluate_by_namespace(raw)
        ns_ic_harm, ns_count_harm = harmonizer.evaluate_by_namespace(harmonized)
        improvement = ((ic_harm - ic_raw) / ic_raw * 100) if ic_raw > 0 else 0.0

        metrics_rows.append({
            "Protein": prot,
            "Raw_Count": len(raw),
            "Harmonized_Count": len(harmonized),
            "IC_Raw": round(ic_raw, 4),
            "IC_Harmonized": round(ic_harm, 4),
            "Improvement_Pct": round(improvement, 2),
            "BP_Raw": ns_count_raw["BP"], "BP_Harmonized": ns_count_harm["BP"],
            "MF_Raw": ns_count_raw["MF"], "MF_Harmonized": ns_count_harm["MF"],
            "CC_Raw": ns_count_raw["CC"], "CC_Harmonized": ns_count_harm["CC"],
            "IC_BP_Raw": round(ns_ic_raw["BP"], 4),
            "IC_BP_Harmonized": round(ns_ic_harm["BP"], 4),
            "IC_MF_Raw": round(ns_ic_raw["MF"], 4),
            "IC_MF_Harmonized": round(ns_ic_harm["MF"], 4),
            "IC_CC_Raw": round(ns_ic_raw["CC"], 4),
            "IC_CC_Harmonized": round(ns_ic_harm["CC"], 4),
        })

        # Store term-level data
        for t in raw:
            terms_rows.append({"protein": prot, "go_id": t, "is_original": True})
        for t in harmonized - raw:
            terms_rows.append({"protein": prot, "go_id": t, "is_original": False})

    # 4. Save
    pd.DataFrame(metrics_rows).to_csv(
        "data/processed/all_protein_metrics.csv", index=False)
    print(f"  Saved all_protein_metrics.csv ({len(metrics_rows)} proteins)")

    pd.DataFrame(terms_rows).to_csv(
        "data/processed/protein_terms.csv.gz", index=False, compression="gzip")
    print(f"  Saved protein_terms.csv.gz ({len(terms_rows)} term entries)")

    print("Done!")


if __name__ == "__main__":
    main()
