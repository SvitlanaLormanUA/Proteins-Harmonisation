import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(
    page_title="GO Harmonization",
    page_icon="\U0001f9ec",
    layout="wide",
)

# ── data loading (cached) ─────────────────────────────────────────────

@st.cache_data
def load_results():
    results = pd.read_csv("data/processed/harmonization_results.csv")
    stats = pd.read_csv("data/processed/statistics_summary.csv")
    evidence = pd.read_csv("data/processed/evidence_breakdown.csv")
    try:
        sim_raw = pd.read_csv("data/processed/similarity_raw.csv", index_col=0)
        sim_harm = pd.read_csv("data/processed/similarity_harmonized.csv", index_col=0)
        case_study = pd.read_csv("data/processed/drug_repurposing_case_study.csv")
    except FileNotFoundError:
        sim_raw = sim_harm = case_study = None
    return results, stats, evidence, sim_raw, sim_harm, case_study


@st.cache_resource
def load_harmonizer():
    from src.data_loader import BioDataLoader
    from src.harmonizer import GOHarmonizer

    loader = BioDataLoader(
        obo_path="data/raw/go-basic.obo",
        gaf_path="data/raw/goa_human_plus.gaf",
    )
    dag = loader.load_ontology()
    loader.load_annotations()
    df = loader.get_df()
    harmonizer = GOHarmonizer(dag, df)
    return harmonizer, dag, df


results, stats, evidence, sim_raw, sim_harm, case_study = load_results()

# ── sidebar ────────────────────────────────────────────────────────────

PAGES = [
    "Overview",
    "IC Analysis",
    "Sub-ontology Breakdown",
    "Evidence Codes",
    "Drug Repurposing",
    "Protein Explorer",
]
page = st.sidebar.radio("Navigation", PAGES)

UNIPROT_TO_GENE = {
    "P00533": "EGFR",  "P04626": "ERBB2",  "P31749": "AKT1",
    "P42336": "PIK3CA", "P15056": "BRAF",   "Q16539": "MAPK14",
    "P42345": "MTOR",   "P07900": "HSP90AA1",
}
DRUG_INFO = {
    "EGFR": "erlotinib, gefitinib",
    "ERBB2": "trastuzumab, lapatinib",
    "AKT1": "MK-2206",
    "PIK3CA": "alpelisib",
    "BRAF": "vemurafenib, dabrafenib",
    "MAPK14": "losmapimod",
    "MTOR": "rapamycin, everolimus",
    "HSP90AA1": "geldanamycin",
}


# ══════════════════════════════════════════════════════════════════════
#  PAGE 1 — Overview
# ══════════════════════════════════════════════════════════════════════

if page == "Overview":
    st.title("GO Annotation Harmonization")
    st.markdown(
        "**Methods for Harmonizing Predicted Functional Protein Annotations "
        "within the Gene Ontology DAG Hierarchy to Improve the Informativeness "
        "of Computational Drug Repurposing Approaches**"
    )
    st.markdown("---")

    # key metrics
    stat_val = lambda name: float(stats.loc[stats["metric"] == name, "value"].values[0])

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Proteins analysed", f"{len(results):,}")
    c2.metric("Mean IC improvement", f"{stat_val('mean_improvement_pct'):.1f}%")
    c3.metric("Median IC improvement", f"{stat_val('median_improvement_pct'):.1f}%")
    c4.metric("Wilcoxon p-value", f"{stat_val('wilcoxon_p_value'):.2e}")

    st.markdown("---")

    col_l, col_r = st.columns(2)
    with col_l:
        st.subheader("Summary")
        summary = pd.DataFrame({
            "Metric": [
                "Mean raw annotations / protein",
                "Mean harmonised annotations / protein",
                "Mean IC (raw)",
                "Mean IC (harmonised)",
            ],
            "Value": [
                f"{results['Raw_Count'].mean():.1f}",
                f"{results['Harmonized_Count'].mean():.1f}",
                f"{results['IC_Raw'].mean():.2f}",
                f"{results['IC_Harmonized'].mean():.2f}",
            ],
        })
        st.dataframe(summary, use_container_width=True, hide_index=True)

    with col_r:
        st.subheader("How it works")
        st.markdown(
            "1. Load Gene Ontology DAG and protein annotations (GOA Human)\n"
            "2. For each protein, propagate annotations up the DAG (**True Path Rule**)\n"
            "3. Compute **Information Content** with propagated term frequencies\n"
            "4. Measure improvement in annotation informativeness\n"
            "5. Evaluate impact on **semantic similarity** between drug targets"
        )

    st.subheader("Full results table")
    st.dataframe(results, use_container_width=True, height=400)


# ══════════════════════════════════════════════════════════════════════
#  PAGE 2 — IC Analysis
# ══════════════════════════════════════════════════════════════════════

elif page == "IC Analysis":
    st.title("Information Content Analysis")
    st.markdown(
        "IC measures how *specific* a protein's annotation profile is. "
        "Harmonization adds ancestor terms, increasing total IC because "
        "propagated term frequencies keep IC accurate (root terms get IC ~ 0)."
    )

    col1, col2 = st.columns(2)

    with col1:
        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=results["IC_Raw"], name="Before", marker_color="#e74c3c", opacity=0.7, nbinsx=40,
        ))
        fig.add_trace(go.Histogram(
            x=results["IC_Harmonized"], name="After", marker_color="#2ecc71", opacity=0.7, nbinsx=40,
        ))
        fig.update_layout(
            barmode="overlay", title="IC Distribution: Before vs After",
            xaxis_title="Total Information Content", yaxis_title="Proteins",
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        fig = px.histogram(results, x="Improvement_Pct", nbins=40,
                           title="IC Improvement Distribution",
                           labels={"Improvement_Pct": "IC Improvement (%)"})
        fig.add_vline(
            x=results["Improvement_Pct"].median(), line_dash="dash", line_color="red",
            annotation_text=f"Median: {results['Improvement_Pct'].median():.1f}%",
        )
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Raw annotation count vs IC improvement")
    fig = px.scatter(
        results, x="Raw_Count", y="Improvement_Pct", color="IC_Raw",
        labels={"Raw_Count": "Raw annotations", "Improvement_Pct": "IC Improvement (%)", "IC_Raw": "Raw IC"},
        color_continuous_scale="viridis", opacity=0.6,
    )
    fig.update_layout(title="Proteins with fewer raw annotations benefit more from harmonization")
    st.plotly_chart(fig, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════
#  PAGE 3 — Sub-ontology Breakdown
# ══════════════════════════════════════════════════════════════════════

elif page == "Sub-ontology Breakdown":
    st.title("Sub-ontology Breakdown (BP / MF / CC)")
    st.markdown(
        "Gene Ontology has three sub-ontologies: **Biological Process** (deepest hierarchy), "
        "**Molecular Function**, and **Cellular Component**. BP benefits the most from "
        "harmonization because its DAG is the deepest."
    )

    # compute per-namespace improvement
    ns_rows = []
    for _, row in results.iterrows():
        for ns in ["BP", "MF", "CC"]:
            raw_ic = row[f"IC_{ns}_Raw"]
            harm_ic = row[f"IC_{ns}_Harmonized"]
            if raw_ic > 0:
                imp = (harm_ic - raw_ic) / raw_ic * 100
                ns_rows.append({"Namespace": ns, "Improvement_Pct": imp,
                                "Raw_Count": row[f"{ns}_Raw"],
                                "Harmonized_Count": row[f"{ns}_Harmonized"]})
    ns_df = pd.DataFrame(ns_rows)

    col1, col2 = st.columns(2)

    with col1:
        fig = px.box(ns_df, x="Namespace", y="Improvement_Pct",
                     color="Namespace",
                     color_discrete_map={"BP": "#3498db", "MF": "#e74c3c", "CC": "#2ecc71"},
                     title="IC Improvement by Sub-ontology")
        fig.update_layout(yaxis_title="IC Improvement (%)", showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        bar_data = pd.DataFrame({
            "Namespace": ["BP", "MF", "CC"] * 2,
            "Stage": ["Before"] * 3 + ["After"] * 3,
            "Mean terms": [
                results["BP_Raw"].mean(), results["MF_Raw"].mean(), results["CC_Raw"].mean(),
                results["BP_Harmonized"].mean(), results["MF_Harmonized"].mean(), results["CC_Harmonized"].mean(),
            ],
        })
        fig = px.bar(bar_data, x="Namespace", y="Mean terms", color="Stage", barmode="group",
                     color_discrete_map={"Before": "#e74c3c", "After": "#2ecc71"},
                     title="Average Term Count by Sub-ontology")
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("Per-namespace statistics")
    ns_stats = ns_df.groupby("Namespace")["Improvement_Pct"].describe().round(1)
    st.dataframe(ns_stats, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════
#  PAGE 4 — Evidence Codes
# ══════════════════════════════════════════════════════════════════════

elif page == "Evidence Codes":
    st.title("Evidence Code Analysis")
    st.markdown(
        "Annotations have different evidence types. **IEA** (electronic/computational) "
        "annotations are the shallowest and benefit the most from harmonization. "
        "**Experimental** annotations already have more terms but still gain significantly."
    )

    if evidence.empty:
        st.warning("No evidence breakdown data found.")
    else:
        grouped = evidence.groupby("Evidence_Category").agg({
            "IC_Raw": "mean", "IC_Harmonized": "mean",
            "Improvement_Pct": "mean", "Raw_Count": "mean",
        }).reset_index().round(2)

        col1, col2 = st.columns(2)

        with col1:
            fig = px.bar(grouped, x="Evidence_Category", y="Improvement_Pct",
                         title="Mean IC Improvement by Evidence Category",
                         labels={"Evidence_Category": "Evidence", "Improvement_Pct": "Improvement (%)"},
                         color="Improvement_Pct", color_continuous_scale="blues")
            fig.update_layout(xaxis_tickangle=-25)
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            fig = go.Figure()
            fig.add_trace(go.Bar(name="Before", x=grouped["Evidence_Category"],
                                 y=grouped["IC_Raw"], marker_color="#e74c3c", opacity=0.8))
            fig.add_trace(go.Bar(name="After", x=grouped["Evidence_Category"],
                                 y=grouped["IC_Harmonized"], marker_color="#2ecc71", opacity=0.8))
            fig.update_layout(barmode="group", title="Mean IC by Evidence Category",
                              xaxis_tickangle=-25, yaxis_title="Mean IC")
            st.plotly_chart(fig, use_container_width=True)

        st.subheader("Detailed breakdown")
        st.dataframe(grouped, use_container_width=True, hide_index=True)


# ══════════════════════════════════════════════════════════════════════
#  PAGE 5 — Drug Repurposing Case Study
# ══════════════════════════════════════════════════════════════════════

elif page == "Drug Repurposing":
    st.title("Drug Repurposing Case Study")
    st.markdown(
        "We selected well-known cancer-pathway drug targets and computed "
        "**Resnik Best-Match-Average (BMA)** semantic similarity before and after "
        "harmonization. Higher self-similarity (diagonal) after harmonization shows "
        "that annotation profiles become richer and more informative."
    )

    if sim_raw is not None:
        # Map UniProt IDs to gene names for display
        uniprot_ids = list(sim_raw.columns)
        gene_names = [UNIPROT_TO_GENE.get(uid, uid) for uid in uniprot_ids]

        info_rows = []
        for uid, gene in zip(uniprot_ids, gene_names):
            drugs = DRUG_INFO.get(gene, "")
            info_rows.append({"UniProt ID": uid, "Gene": gene, "Known drugs": drugs})
        st.dataframe(pd.DataFrame(info_rows), use_container_width=True, hide_index=True)

        col1, col2 = st.columns(2)
        vmax = max(sim_raw.values.max(), sim_harm.values.max())

        with col1:
            fig = px.imshow(sim_raw.values, x=gene_names, y=gene_names,
                            color_continuous_scale="YlOrRd", zmin=0, zmax=vmax,
                            text_auto=".2f", title="Before Harmonization")
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            fig = px.imshow(sim_harm.values, x=gene_names, y=gene_names,
                            color_continuous_scale="YlOrRd", zmin=0, zmax=vmax,
                            text_auto=".2f", title="After Harmonization")
            st.plotly_chart(fig, use_container_width=True)

        if case_study is not None:
            st.subheader("Pairwise similarity changes")
            st.dataframe(case_study, use_container_width=True, hide_index=True)

            n = len(gene_names)
            raw_vals = [sim_raw.values[i, j] for i in range(n) for j in range(i + 1, n)]
            harm_vals = [sim_harm.values[i, j] for i in range(n) for j in range(i + 1, n)]
            c1, c2 = st.columns(2)
            c1.metric("Mean pairwise similarity (raw)", f"{np.mean(raw_vals):.4f}")
            c2.metric("Mean pairwise similarity (harmonised)", f"{np.mean(harm_vals):.4f}")
    else:
        st.info("Run main.py first to generate similarity matrices.")


# ══════════════════════════════════════════════════════════════════════
#  PAGE 6 — Protein Explorer
# ══════════════════════════════════════════════════════════════════════

elif page == "Protein Explorer":
    st.title("Protein Explorer")
    st.markdown(
        "Search for any protein in the GOA Human dataset to see its "
        "annotation profile before and after harmonization."
    )

    if "harmonizer_loaded" not in st.session_state:
        st.session_state.harmonizer_loaded = False

    if not st.session_state.harmonizer_loaded:
        if st.button("Load full model (takes ~30 s on first run)"):
            with st.spinner("Loading Gene Ontology and annotations..."):
                harmonizer, dag, ann_df = load_harmonizer()
                st.session_state.harmonizer_loaded = True
                st.session_state.harmonizer = harmonizer
                st.session_state.dag = dag
                st.session_state.ann_df = ann_df
            st.rerun()
        st.stop()

    harmonizer = st.session_state.harmonizer
    ann_df = st.session_state.ann_df

    all_proteins = sorted(ann_df["DB_Object_ID"].unique())

    protein_id = st.selectbox(
        "Select or type a protein ID",
        options=all_proteins,
        index=None,
        placeholder="Start typing a protein ID...",
    )

    if protein_id:
        raw_terms, harm_terms = harmonizer.harmonize_protein(protein_id)
        ic_raw = harmonizer.evaluate_informativeness(raw_terms)
        ic_harm = harmonizer.evaluate_informativeness(harm_terms)
        ns_ic_raw, ns_count_raw = harmonizer.evaluate_by_namespace(raw_terms)
        ns_ic_harm, ns_count_harm = harmonizer.evaluate_by_namespace(harm_terms)

        improvement = ((ic_harm - ic_raw) / ic_raw * 100) if ic_raw > 0 else 0

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Raw terms", len(raw_terms))
        c2.metric("Harmonised terms", len(harm_terms))
        c3.metric("IC (raw)", f"{ic_raw:.2f}")
        c4.metric("IC improvement", f"{improvement:.1f}%")

        # namespace breakdown chart
        ns_data = pd.DataFrame({
            "Namespace": ["BP", "MF", "CC"] * 2,
            "Stage": ["Before"] * 3 + ["After"] * 3,
            "Terms": [ns_count_raw["BP"], ns_count_raw["MF"], ns_count_raw["CC"],
                      ns_count_harm["BP"], ns_count_harm["MF"], ns_count_harm["CC"]],
            "IC": [ns_ic_raw["BP"], ns_ic_raw["MF"], ns_ic_raw["CC"],
                   ns_ic_harm["BP"], ns_ic_harm["MF"], ns_ic_harm["CC"]],
        })

        col_l, col_r = st.columns(2)
        with col_l:
            fig = px.bar(ns_data, x="Namespace", y="Terms", color="Stage", barmode="group",
                         color_discrete_map={"Before": "#e74c3c", "After": "#2ecc71"},
                         title="Term count by sub-ontology")
            st.plotly_chart(fig, use_container_width=True)
        with col_r:
            fig = px.bar(ns_data, x="Namespace", y="IC", color="Stage", barmode="group",
                         color_discrete_map={"Before": "#e74c3c", "After": "#2ecc71"},
                         title="IC by sub-ontology")
            st.plotly_chart(fig, use_container_width=True)

        # term tables
        propagated_only = harm_terms - raw_terms
        dag = st.session_state.dag

        def term_table(terms):
            rows = []
            for t in sorted(terms):
                name = dag[t].name if t in dag else ""
                ns = harmonizer.get_term_namespace(t)
                ic = harmonizer.get_information_content(t)
                rows.append({"GO ID": t, "Name": name, "Namespace": ns, "IC": round(ic, 3)})
            return pd.DataFrame(rows)

        tab1, tab2 = st.tabs(["Original annotations", "Propagated (new) annotations"])
        with tab1:
            st.dataframe(term_table(raw_terms), use_container_width=True, hide_index=True)
        with tab2:
            st.dataframe(term_table(propagated_only), use_container_width=True, hide_index=True)
