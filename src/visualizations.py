import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx


# Clean style for thesis figures
plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
})


class ResultsVisualizer:
    def __init__(self, output_dir='data/processed/figures'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def _save(self, fig, filename):
        path = os.path.join(self.output_dir, filename)
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {path}")

    # ------------------------------------------------------------------
    # 1. IC distribution: before vs after harmonization
    # ------------------------------------------------------------------
    def plot_ic_distribution(self, results_df):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        axes[0].hist(results_df['IC_Raw'], bins=40, alpha=0.7,
                     label='Before harmonization', color='#e74c3c')
        axes[0].hist(results_df['IC_Harmonized'], bins=40, alpha=0.7,
                     label='After harmonization', color='#2ecc71')
        axes[0].set_xlabel('Total Information Content')
        axes[0].set_ylabel('Number of proteins')
        axes[0].set_title('IC Distribution: Before vs After Harmonization')
        axes[0].legend()

        axes[1].hist(results_df['Improvement_Pct'], bins=40, color='#3498db', alpha=0.8)
        axes[1].set_xlabel('IC Improvement (%)')
        axes[1].set_ylabel('Number of proteins')
        axes[1].set_title('Distribution of IC Improvement')
        median_val = results_df['Improvement_Pct'].median()
        axes[1].axvline(median_val, color='red', linestyle='--',
                        label=f'Median: {median_val:.1f}%')
        axes[1].legend()

        fig.tight_layout()
        self._save(fig, 'ic_distribution.png')

    # ------------------------------------------------------------------
    # 2. Sub-ontology (namespace) breakdown
    # ------------------------------------------------------------------
    def plot_namespace_breakdown(self, results_df):
        bp_imp, mf_imp, cc_imp = [], [], []
        for _, row in results_df.iterrows():
            if row['IC_BP_Raw'] > 0:
                bp_imp.append((row['IC_BP_Harmonized'] - row['IC_BP_Raw'])
                              / row['IC_BP_Raw'] * 100)
            if row['IC_MF_Raw'] > 0:
                mf_imp.append((row['IC_MF_Harmonized'] - row['IC_MF_Raw'])
                              / row['IC_MF_Raw'] * 100)
            if row['IC_CC_Raw'] > 0:
                cc_imp.append((row['IC_CC_Harmonized'] - row['IC_CC_Raw'])
                              / row['IC_CC_Raw'] * 100)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Box plot: improvement by namespace
        data = [bp_imp, mf_imp, cc_imp]
        labels = ['Biological\nProcess', 'Molecular\nFunction', 'Cellular\nComponent']
        colors = ['#3498db', '#e74c3c', '#2ecc71']

        bp = axes[0].boxplot(data, labels=labels, patch_artist=True, showfliers=False)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        axes[0].set_ylabel('IC Improvement (%)')
        axes[0].set_title('IC Improvement by GO Sub-ontology')

        # Bar chart: average term count before/after
        ns_labels = ['BP', 'MF', 'CC']
        raw_means = [results_df['BP_Raw'].mean(),
                     results_df['MF_Raw'].mean(),
                     results_df['CC_Raw'].mean()]
        harm_means = [results_df['BP_Harmonized'].mean(),
                      results_df['MF_Harmonized'].mean(),
                      results_df['CC_Harmonized'].mean()]

        x = np.arange(len(ns_labels))
        width = 0.35
        axes[1].bar(x - width / 2, raw_means, width, label='Before', color='#e74c3c', alpha=0.7)
        axes[1].bar(x + width / 2, harm_means, width, label='After', color='#2ecc71', alpha=0.7)
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(ns_labels)
        axes[1].set_ylabel('Average term count')
        axes[1].set_title('Average Term Count by Sub-ontology')
        axes[1].legend()

        fig.tight_layout()
        self._save(fig, 'namespace_breakdown.png')

        return {'BP': bp_imp, 'MF': mf_imp, 'CC': cc_imp}

    # ------------------------------------------------------------------
    # 3. Evidence code breakdown
    # ------------------------------------------------------------------
    def plot_evidence_breakdown(self, evidence_df):
        if evidence_df.empty:
            print("  No evidence data to plot.")
            return

        grouped = evidence_df.groupby('Evidence_Category').agg({
            'IC_Raw': 'mean',
            'IC_Harmonized': 'mean',
            'Improvement_Pct': 'mean',
            'Raw_Count': 'mean',
        }).reset_index()

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        categories = grouped['Evidence_Category']
        x = np.arange(len(categories))

        axes[0].bar(x, grouped['Improvement_Pct'], color='#3498db', alpha=0.8)
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(categories, rotation=25, ha='right')
        axes[0].set_ylabel('Mean IC Improvement (%)')
        axes[0].set_title('IC Improvement by Evidence Category')

        width = 0.35
        axes[1].bar(x - width / 2, grouped['IC_Raw'], width,
                    label='Before', color='#e74c3c', alpha=0.7)
        axes[1].bar(x + width / 2, grouped['IC_Harmonized'], width,
                    label='After', color='#2ecc71', alpha=0.7)
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(categories, rotation=25, ha='right')
        axes[1].set_ylabel('Mean IC')
        axes[1].set_title('Mean IC by Evidence Category')
        axes[1].legend()

        fig.tight_layout()
        self._save(fig, 'evidence_breakdown.png')

    # ------------------------------------------------------------------
    # 4. Scatter: raw annotation count vs IC improvement
    # ------------------------------------------------------------------
    def plot_scatter_raw_vs_improvement(self, results_df):
        fig, ax = plt.subplots(figsize=(10, 7))

        scatter = ax.scatter(
            results_df['Raw_Count'],
            results_df['Improvement_Pct'],
            c=results_df['IC_Raw'],
            cmap='viridis',
            alpha=0.6,
            s=30,
        )
        plt.colorbar(scatter, label='Raw IC')
        ax.set_xlabel('Number of Raw Annotations')
        ax.set_ylabel('IC Improvement (%)')
        ax.set_title('Raw Annotation Count vs IC Improvement')

        fig.tight_layout()
        self._save(fig, 'scatter_raw_vs_improvement.png')

    # ------------------------------------------------------------------
    # 5. DAG subgraph for one protein
    # ------------------------------------------------------------------
    def plot_dag_subgraph(self, godag, raw_terms, harmonized_terms, protein_name='Protein'):
        G = nx.DiGraph()

        for go_id in harmonized_terms:
            if go_id not in godag:
                continue
            term = godag[go_id]
            label = f"{go_id}\n{term.name[:30]}" if hasattr(term, 'name') else go_id
            G.add_node(go_id, label=label)
            for parent in term.parents:
                if parent.id in harmonized_terms:
                    G.add_edge(go_id, parent.id)

        if not G.nodes:
            print("  No valid nodes for DAG visualization.")
            return

        # Limit graph size for readability
        if len(G.nodes) > 35:
            keep = set(raw_terms)
            for t in raw_terms:
                if t in godag:
                    for p in godag[t].parents:
                        keep.add(p.id)
                        for gp in p.parents:
                            keep.add(gp.id)
            G = G.subgraph([n for n in G.nodes if n in keep]).copy()

        fig, ax = plt.subplots(figsize=(16, 11))

        # Use graphviz layout if available, otherwise spring layout
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except Exception:
            pos = nx.spring_layout(G, k=2.5, iterations=80, seed=42)

        colors = ['#e74c3c' if node in raw_terms else '#2ecc71' for node in G.nodes]

        nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=900, alpha=0.8, ax=ax)
        nx.draw_networkx_edges(G, pos, edge_color='#95a5a6', arrows=True,
                               arrowsize=15, width=1.5, ax=ax)

        labels = nx.get_node_attributes(G, 'label')
        nx.draw_networkx_labels(G, pos, labels, font_size=6, ax=ax)

        original_patch = mpatches.Patch(color='#e74c3c', label='Original annotations')
        propagated_patch = mpatches.Patch(color='#2ecc71', label='Propagated (harmonized)')
        ax.legend(handles=[original_patch, propagated_patch], loc='upper left', fontsize=10)

        ax.set_title(
            f'GO DAG Subgraph for {protein_name}\n'
            f'(Original: {len(raw_terms)} terms -> Harmonized: {len(harmonized_terms)} terms)',
            fontsize=13,
        )
        ax.axis('off')

        fig.tight_layout()
        self._save(fig, 'dag_subgraph.png')

    # ------------------------------------------------------------------
    # 6. Similarity heatmaps (before vs after harmonization)
    # ------------------------------------------------------------------
    def plot_similarity_heatmap(self, sim_matrix_raw, sim_matrix_harm, protein_labels):
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        vmax = max(sim_matrix_raw.values.max(), sim_matrix_harm.values.max())
        if vmax == 0:
            vmax = 1.0

        for ax, matrix, title in [
            (axes[0], sim_matrix_raw, 'Before Harmonization'),
            (axes[1], sim_matrix_harm, 'After Harmonization'),
        ]:
            im = ax.imshow(matrix.values, cmap='YlOrRd', vmin=0, vmax=vmax)
            ax.set_xticks(range(len(protein_labels)))
            ax.set_yticks(range(len(protein_labels)))
            ax.set_xticklabels(protein_labels, rotation=45, ha='right', fontsize=8)
            ax.set_yticklabels(protein_labels, fontsize=8)
            ax.set_title(f'Protein Similarity\n({title})')
            plt.colorbar(im, ax=ax, shrink=0.8)

            for i in range(len(protein_labels)):
                for j in range(len(protein_labels)):
                    ax.text(j, i, f'{matrix.values[i, j]:.2f}',
                            ha='center', va='center', fontsize=7)

        fig.suptitle('Drug Repurposing Case Study: Protein Semantic Similarity',
                     fontsize=13, y=1.02)
        fig.tight_layout()
        self._save(fig, 'similarity_heatmap.png')
