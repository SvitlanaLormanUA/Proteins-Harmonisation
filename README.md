# Proteins-Harmonisation

**Deployed:** https://proteins-harmonisation.streamlit.app/

Code for a thesis evaluating methods to harmonize predicted protein annotations within the Gene Ontology (GO) DAG hierarchy. By enforcing hierarchical consistency, this model improves the informativeness and accuracy of computational drug repurposing pipelines. Includes data preprocessing, GO-term harmonization, and drug-target interaction analysis.

```text
Proteins-Harmonisation/
├── .gitignore
├── README.md
├── app.py                      # Streamlit dashboard
├── main.py                     # CLI entry point
├── precompute_explorer.py      # Pre-computes protein explorer data
├── requirements.txt            # Production dependencies
├── requirements-dev.txt        # Dev dependencies
├── src/                        # Core library
│   ├── __init__.py
│   ├── config.py
│   ├── data_loader.py
│   ├── database.py
│   ├── harmonizer.py
│   ├── main.py
│   ├── metrics.py
│   └── visualizations.py
├── data/
│   ├── raw/                    # (gitignored except .gitkeep)
│   │   ├── go-basic.obo
│   │   └── goa_human_plus.gaf
│   └── processed/
│       ├── harmonization_results.csv
│       ├── statistics_summary.csv
│       ├── evidence_breakdown.csv
│       ├── similarity_raw.csv
│       ├── similarity_harmonized.csv
│       ├── drug_repurposing_case_study.csv
│       ├── all_protein_metrics.csv   # Pre-computed explorer
│       ├── protein_terms.csv.gz      # Pre-computed explorer
│       ├── go_terms_info.csv         # Pre-computed explorer
│       ├── results.db                # (gitignored)
│       └── figures/
│           ├── dag_subgraph.png
│           ├── evidence_breakdown.png
│           ├── ic_distribution.png
│           ├── namespace_breakdown.png
│           ├── scatter_raw_vs_improvement.png
│           └── similarity_heatmap.png
├── notebooks/                  # (empty)
└── venv/                       # (gitignored)
