# Proteins-Harmonisation

Code for a thesis evaluating methods to harmonize predicted protein annotations within the Gene Ontology DAG hierarchy. By enforcing hierarchical consistency, this model improves the informativeness and accuracy of computational drug repurposing pipelines. Includes data preprocessing, GO-term harmonization, and drug-target interaction analysis.

Project structure:
bio-thesis/
├── venv/  
├── data/
│ ├── raw/  
│ └── processed/  
├── notebooks/ # Jupyter Notebooks
│ └── 01_explore_data.ipynb
├── src/  
│ ├── **init**.py
│ ├── config.py  
│ ├── data_loader.py  
│ ├── harmonizer.py  
│ └── metrics.py # Semantic Similarity -- results evaluation
├── main.py  
├── requirements.txt  
└── README.md
