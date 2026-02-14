import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from src.data_loader import BioDataLoader
from src.harmonizer import GOHarmonizer
from src.metrics import SemanticSimilarity
from src.visualizations import ResultsVisualizer
from src.database import ResultsDB


# Well-known drug targets (UniProt IDs) for the case study.
# These are all cancer-pathway proteins targeted by approved drugs.
DRUG_TARGETS = {
    'P00533': 'EGFR',       # erlotinib, gefitinib
    'P04626': 'ERBB2',      # trastuzumab, lapatinib
    'P31749': 'AKT1',       # MK-2206
    'P42336': 'PIK3CA',     # alpelisib
    'P15056': 'BRAF',       # vemurafenib, dabrafenib
    'Q16539': 'MAPK14',     # losmapimod
    'P42345': 'MTOR',       # rapamycin, everolimus
    'P07900': 'HSP90AA1',   # geldanamycin
}


def main():
    np.random.seed(42)
    os.makedirs("data/processed/figures", exist_ok=True)

    # ==================================================================
    # STAGE 1: Load data
    # ==================================================================
    print("=" * 60)
    print("STAGE 1: Loading Data")
    print("=" * 60)

    loader = BioDataLoader(
        obo_path="data/raw/go-basic.obo",
        gaf_path="data/raw/goa_human_plus.gaf",
    )
    dag = loader.load_ontology()
    loader.load_annotations()
    df = loader.get_df()

    print(f"Total annotations: {len(df)}")
    print(f"Unique proteins:   {df['DB_Object_ID'].nunique()}")
    print(f"Unique GO terms:   {df['GO_ID'].nunique()}")
    print(f"Evidence codes:    {sorted(df['Evidence'].unique())}")

    # ==================================================================
    # STAGE 2: Initialise harmoniser (propagated IC)
    # ==================================================================
    print("\n" + "=" * 60)
    print("STAGE 2: Initialising Harmoniser (Propagated IC)")
    print("=" * 60)

    harmonizer = GOHarmonizer(dag, df)

    # ==================================================================
    # STAGE 3: Main experiment — 1000 proteins
    # ==================================================================
    print("\n" + "=" * 60)
    print("STAGE 3: Running Harmonisation Experiment (n=1000)")
    print("=" * 60)

    results = harmonizer.run_experiment(sample_size=1000)

    results.to_csv("data/processed/harmonization_results.csv", index=False)
    print(f"\nResults saved to CSV ({len(results)} proteins)")

    db = ResultsDB()
    exp_id = db.save_results(results, sample_size=len(results))
    print(f"Results saved to SQLite (experiment_id={exp_id})")

    print("\n--- Summary Statistics ---")
    print(f"  Mean IC (Raw):          {results['IC_Raw'].mean():.2f}")
    print(f"  Mean IC (Harmonised):   {results['IC_Harmonized'].mean():.2f}")
    print(f"  Mean Improvement:       {results['Improvement_Pct'].mean():.1f}%")
    print(f"  Median Improvement:     {results['Improvement_Pct'].median():.1f}%")
    print(f"  Mean Raw Count:         {results['Raw_Count'].mean():.1f}")
    print(f"  Mean Harmonised Count:  {results['Harmonized_Count'].mean():.1f}")

    # ==================================================================
    # STAGE 4: Statistical significance tests
    # ==================================================================
    print("\n" + "=" * 60)
    print("STAGE 4: Statistical Significance Tests")
    print("=" * 60)

    valid = results[results['IC_Raw'] > 0]
    stat, p_value = wilcoxon(valid['IC_Raw'], valid['IC_Harmonized'])
    print(f"\nWilcoxon signed-rank test (overall):")
    print(f"  Test statistic: {stat:.4f}")
    print(f"  P-value:        {p_value:.2e}")
    print(f"  Significant (p < 0.001): {'YES' if p_value < 0.001 else 'NO'}")

    for ns in ['BP', 'MF', 'CC']:
        ns_valid = results[results[f'IC_{ns}_Raw'] > 0]
        if len(ns_valid) > 20:
            stat_ns, p_ns = wilcoxon(ns_valid[f'IC_{ns}_Raw'],
                                     ns_valid[f'IC_{ns}_Harmonized'])
            print(f"\n  {ns} sub-ontology:")
            print(f"    Wilcoxon statistic: {stat_ns:.4f}")
            print(f"    P-value: {p_ns:.2e}")

    stats_summary = pd.DataFrame({
        'metric': [
            'mean_ic_raw', 'mean_ic_harmonized',
            'mean_improvement_pct', 'median_improvement_pct',
            'wilcoxon_statistic', 'wilcoxon_p_value',
            'sample_size', 'significant_p001',
        ],
        'value': [
            results['IC_Raw'].mean(), results['IC_Harmonized'].mean(),
            results['Improvement_Pct'].mean(), results['Improvement_Pct'].median(),
            stat, p_value, len(results), p_value < 0.001,
        ],
    })
    stats_summary.to_csv("data/processed/statistics_summary.csv", index=False)

    # ==================================================================
    # STAGE 5: Evidence-code breakdown
    # ==================================================================
    print("\n" + "=" * 60)
    print("STAGE 5: Evidence Code Breakdown")
    print("=" * 60)

    sample_for_evidence = results['Protein'].values[:200]
    evidence_results = harmonizer.get_evidence_breakdown(sample_for_evidence)
    evidence_results.to_csv("data/processed/evidence_breakdown.csv", index=False)

    print("\nEvidence category summary (mean values):")
    if not evidence_results.empty:
        evidence_summary = evidence_results.groupby('Evidence_Category').agg({
            'IC_Raw': 'mean',
            'IC_Harmonized': 'mean',
            'Improvement_Pct': 'mean',
            'Raw_Count': 'mean',
        }).round(2)
        print(evidence_summary.to_string())

    # ==================================================================
    # STAGE 6: Visualisations
    # ==================================================================
    print("\n" + "=" * 60)
    print("STAGE 6: Creating Visualisations")
    print("=" * 60)

    viz = ResultsVisualizer()
    viz.plot_ic_distribution(results)
    viz.plot_namespace_breakdown(results)
    viz.plot_evidence_breakdown(evidence_results)
    viz.plot_scatter_raw_vs_improvement(results)

    # DAG subgraph — pick a protein with a moderate number of annotations
    example_protein = results.sort_values('Raw_Count').iloc[len(results) // 2]['Protein']
    raw_terms, harm_terms = harmonizer.harmonize_protein(example_protein)
    viz.plot_dag_subgraph(dag, raw_terms, harm_terms, protein_name=example_protein)

    # ==================================================================
    # STAGE 7: Drug-repurposing case study
    # ==================================================================
    print("\n" + "=" * 60)
    print("STAGE 7: Drug Repurposing Case Study")
    print("=" * 60)

    available_proteins = set(df['DB_Object_ID'].unique())
    found_targets = {uid: name for uid, name in DRUG_TARGETS.items()
                     if uid in available_proteins}

    print(f"Found {len(found_targets)}/{len(DRUG_TARGETS)} drug targets in dataset:")
    for uid, name in found_targets.items():
        n_terms = len(harmonizer.get_protein_terms(uid))
        print(f"  {uid} ({name}): {n_terms} direct annotations")

    # Need at least 3 targets for a meaningful case study
    if len(found_targets) >= 3:
        target_ids = list(found_targets.keys())[:6]
        target_names = [found_targets[t] for t in target_ids]
    else:
        # Fallback: pick proteins with the most annotations
        print("\nNot enough known targets found. Using top-annotated proteins instead.")
        protein_counts = df.groupby('DB_Object_ID').size().sort_values(ascending=False)
        target_ids = protein_counts.head(5).index.tolist()
        target_names = target_ids

    print(f"\nComputing Resnik BMA similarity for {len(target_ids)} proteins...")
    sim = SemanticSimilarity(harmonizer)

    print("  Raw similarity matrix...")
    sim_raw = sim.pairwise_similarity_matrix(target_ids, use_harmonized=False)

    print("  Harmonised similarity matrix...")
    sim._mica_cache.clear()
    sim_harm = sim.pairwise_similarity_matrix(target_ids, use_harmonized=True)

    sim_raw.to_csv("data/processed/similarity_raw.csv")
    sim_harm.to_csv("data/processed/similarity_harmonized.csv")

    print("\nRaw similarity matrix:")
    print(sim_raw.round(3).to_string())
    print("\nHarmonised similarity matrix:")
    print(sim_harm.round(3).to_string())

    viz.plot_similarity_heatmap(sim_raw, sim_harm, target_names)

    # Pairwise comparison summary
    n = len(target_ids)
    raw_vals = [sim_raw.values[i, j] for i in range(n) for j in range(i + 1, n)]
    harm_vals = [sim_harm.values[i, j] for i in range(n) for j in range(i + 1, n)]

    print(f"\nMean pairwise similarity (raw):        {np.mean(raw_vals):.4f}")
    print(f"Mean pairwise similarity (harmonised): {np.mean(harm_vals):.4f}")

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    case_study = pd.DataFrame({
        'Protein_1': [f"{target_names[i]} ({target_ids[i]})" for i, _ in pairs],
        'Protein_2': [f"{target_names[j]} ({target_ids[j]})" for _, j in pairs],
        'Similarity_Raw': raw_vals,
        'Similarity_Harmonised': harm_vals,
        'Delta': [h - r for r, h in zip(raw_vals, harm_vals)],
    })
    case_study.to_csv("data/processed/drug_repurposing_case_study.csv", index=False)
    print("Case study results saved.")

    # ==================================================================
    # DONE
    # ==================================================================
    print("\n" + "=" * 60)
    print("ALL STAGES COMPLETE")
    print("=" * 60)
    print("\nOutput files:")
    for f in [
        "data/processed/harmonization_results.csv",
        "data/processed/statistics_summary.csv",
        "data/processed/evidence_breakdown.csv",
        "data/processed/similarity_raw.csv",
        "data/processed/similarity_harmonized.csv",
        "data/processed/drug_repurposing_case_study.csv",
        "data/processed/figures/ic_distribution.png",
        "data/processed/figures/namespace_breakdown.png",
        "data/processed/figures/evidence_breakdown.png",
        "data/processed/figures/scatter_raw_vs_improvement.png",
        "data/processed/figures/dag_subgraph.png",
        "data/processed/figures/similarity_heatmap.png",
    ]:
        print(f"  {f}")


if __name__ == "__main__":
    main()
