import numpy as np
import pandas as pd


class SemanticSimilarity:
    """
    Resnik semantic similarity between GO terms and
    Best Match Average (BMA) for protein-level comparison.

    Resnik similarity: sim(t1, t2) = IC of the Most Informative Common Ancestor (MICA).
    BMA: averages the best pairwise Resnik matches between two proteins' term sets,
         computed separately per namespace (BP, MF, CC).
    """

    def __init__(self, harmonizer):
        self.harmonizer = harmonizer
        self.godag = harmonizer.godag
        self._mica_cache = {}

    def get_ancestors(self, go_id):
        """Get all ancestors of a GO term (including self)."""
        ancestors = {go_id}
        if go_id in self.godag:
            ancestors.update(self.godag[go_id].get_all_parents())
        return ancestors

    def get_common_ancestors(self, term1, term2):
        """Get set of common ancestors of two GO terms."""
        return self.get_ancestors(term1) & self.get_ancestors(term2)

    def resnik_similarity(self, term1, term2):
        """
        Resnik similarity = IC of the Most Informative Common Ancestor (MICA).
        sim_Resnik(t1, t2) = max{ IC(t) : t in ancestors(t1) ∩ ancestors(t2) }
        """
        key = (min(term1, term2), max(term1, term2))
        if key in self._mica_cache:
            return self._mica_cache[key]

        common = self.get_common_ancestors(term1, term2)
        if not common:
            self._mica_cache[key] = 0.0
            return 0.0

        max_ic = max(self.harmonizer.get_information_content(t) for t in common)
        self._mica_cache[key] = max_ic
        return max_ic

    def _bma_for_namespace(self, terms1, terms2):
        """
        BMA between two term lists (assumed to be from the same namespace).
        For each term in set A, find best Resnik match in set B, and vice versa.
        BMA = (mean(best_A→B) + mean(best_B→A)) / 2
        """
        if not terms1 or not terms2:
            return 0.0

        best_1to2 = []
        for t1 in terms1:
            best = max(self.resnik_similarity(t1, t2) for t2 in terms2)
            best_1to2.append(best)

        best_2to1 = []
        for t2 in terms2:
            best = max(self.resnik_similarity(t2, t1) for t1 in terms1)
            best_2to1.append(best)

        return (np.mean(best_1to2) + np.mean(best_2to1)) / 2

    def protein_similarity_bma(self, terms1, terms2):
        """
        BMA similarity between two proteins given their GO term lists.
        Computed separately per namespace and averaged.
        """
        # Filter to valid terms
        terms1 = [t for t in terms1 if t in self.godag]
        terms2 = [t for t in terms2 if t in self.godag]

        if not terms1 or not terms2:
            return 0.0

        ns_sims = []
        for ns in ['biological_process', 'molecular_function', 'cellular_component']:
            t1_ns = [t for t in terms1 if self.godag[t].namespace == ns]
            t2_ns = [t for t in terms2 if self.godag[t].namespace == ns]
            if t1_ns and t2_ns:
                ns_sims.append(self._bma_for_namespace(t1_ns, t2_ns))

        return np.mean(ns_sims) if ns_sims else 0.0

    def pairwise_similarity_matrix(self, protein_ids, use_harmonized=True):
        """
        Compute pairwise BMA similarity matrix for a list of proteins.
        Pre-computes term sets to avoid redundant harmonization calls.
        """
        n = len(protein_ids)

        # Pre-compute terms for all proteins
        protein_terms = {}
        for pid in protein_ids:
            if use_harmonized:
                _, terms = self.harmonizer.harmonize_protein(pid)
                protein_terms[pid] = list(terms)
            else:
                protein_terms[pid] = self.harmonizer.get_protein_terms(pid)

        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                sim = self.protein_similarity_bma(
                    protein_terms[protein_ids[i]],
                    protein_terms[protein_ids[j]],
                )
                matrix[i, j] = sim
                matrix[j, i] = sim
            print(f"  Similarity row {i + 1}/{n} done.")

        return pd.DataFrame(matrix, index=protein_ids, columns=protein_ids)
