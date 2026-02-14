import math
import numpy as np
import pandas as pd
from collections import Counter, defaultdict


class GOHarmonizer:
    """
    Harmonizes GO annotations using the True Path Rule
    and computes Information Content with propagated term frequencies.

    Key improvement over naive IC: term frequencies are propagated up the DAG,
    so freq(t) = number of proteins annotated to t OR any of its descendants.
    This gives root terms IC ≈ 0 and specific leaf terms high IC.
    """

    NAMESPACES = {
        'biological_process': 'BP',
        'molecular_function': 'MF',
        'cellular_component': 'CC',
    }

    EVIDENCE_CATEGORIES = {
        'Experimental': {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'},
        'High-throughput': {'HTP', 'HDA', 'HMP', 'HGI', 'HEP'},
        'Electronic (IEA)': {'IEA'},
        'Computational': {'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA'},
        'Author/Curator': {'TAS', 'NAS', 'IC', 'ND'},
    }

    def __init__(self, godag, associations_df):
        """
        godag: GODag object from goatools
        associations_df: DataFrame with columns ['DB_Object_ID', 'GO_ID', 'Evidence']
        """
        self.godag = godag
        self.df = associations_df
        self.total_proteins = self.df['DB_Object_ID'].nunique()
        print(f"Unique proteins: {self.total_proteins}")
        print(f"Total annotations: {len(self.df)}")
        self.term_counts = self._calculate_propagated_counts()

    def _calculate_propagated_counts(self):
        """
        Propagated term frequency computation.
        For each GO term t: freq(t) = |{proteins annotated to t or any descendant of t}|.
        This is done by collecting protein sets per term, then propagating up to all ancestors.
        """
        print("Calculating propagated term frequencies (this may take a minute)...")

        # Step 1: for each GO term, gather the set of proteins directly annotated to it
        term_proteins = self.df.groupby('GO_ID')['DB_Object_ID'].apply(set).to_dict()

        # Step 2: propagate protein sets up to all ancestors in the DAG
        propagated = defaultdict(set)
        for go_id, proteins in term_proteins.items():
            propagated[go_id].update(proteins)
            if go_id in self.godag:
                for parent_id in self.godag[go_id].get_all_parents():
                    propagated[parent_id].update(proteins)

        counts = Counter({go_id: len(prots) for go_id, prots in propagated.items()})
        print(f"Propagated counts computed for {len(counts)} terms.")
        return counts

    def get_information_content(self, go_id):
        """
        IC(t) = -ln(P(t)), where P(t) = freq(t) / N.
        N = total number of annotated proteins.
        Root terms have P ≈ 1 → IC ≈ 0.
        Specific leaf terms have P << 1 → high IC.
        """
        if go_id not in self.term_counts or self.term_counts[go_id] == 0:
            return 0.0
        prob = self.term_counts[go_id] / self.total_proteins
        if prob >= 1.0:
            return 0.0
        return -math.log(prob)

    def get_term_namespace(self, go_id):
        """Returns the sub-ontology abbreviation: 'BP', 'MF', or 'CC'."""
        if go_id in self.godag:
            ns = self.godag[go_id].namespace
            return self.NAMESPACES.get(ns, 'unknown')
        return 'unknown'

    def get_protein_terms(self, protein_id):
        """Returns list of GO term IDs directly annotated to a protein."""
        return self.df[self.df['DB_Object_ID'] == protein_id]['GO_ID'].tolist()

    def get_protein_evidence(self, protein_id):
        """Returns DataFrame of (GO_ID, Evidence) for a protein."""
        return self.df[self.df['DB_Object_ID'] == protein_id][['GO_ID', 'Evidence']]

    def harmonize_protein(self, protein_id):
        """
        True Path Rule harmonization.
        For each directly annotated term, add all ancestor terms in the DAG.
        Returns (original_terms, harmonized_terms) as sets.
        """
        original_terms = set(self.get_protein_terms(protein_id))
        harmonized_terms = set(original_terms)

        for go_id in original_terms:
            if go_id in self.godag:
                harmonized_terms.update(self.godag[go_id].get_all_parents())

        return original_terms, harmonized_terms

    def evaluate_informativeness(self, terms):
        """Sum of IC across all terms in the set."""
        return sum(self.get_information_content(t) for t in terms)

    def evaluate_by_namespace(self, terms):
        """
        Returns (ns_ic, ns_count) dicts broken down by sub-ontology.
        ns_ic: {'BP': total_ic, 'MF': ..., 'CC': ...}
        ns_count: {'BP': num_terms, ...}
        """
        ns_ic = {'BP': 0.0, 'MF': 0.0, 'CC': 0.0}
        ns_count = {'BP': 0, 'MF': 0, 'CC': 0}
        for t in terms:
            ns = self.get_term_namespace(t)
            if ns in ns_ic:
                ns_ic[ns] += self.get_information_content(t)
                ns_count[ns] += 1
        return ns_ic, ns_count

    def run_experiment(self, sample_size=1000):
        """
        Run harmonization on a sample of proteins.
        Returns DataFrame with per-protein results including namespace breakdown.
        """
        unique_proteins = self.df['DB_Object_ID'].unique()
        actual_size = min(sample_size, len(unique_proteins))
        sample = np.random.choice(unique_proteins, actual_size, replace=False)

        print(f"Running harmonization on {actual_size} proteins...")
        results = []

        for i, prot in enumerate(sample):
            if (i + 1) % 200 == 0:
                print(f"  Processed {i + 1}/{actual_size}...")

            raw, harmonized = self.harmonize_protein(prot)
            ic_raw = self.evaluate_informativeness(raw)
            ic_harm = self.evaluate_informativeness(harmonized)

            ns_ic_raw, ns_count_raw = self.evaluate_by_namespace(raw)
            ns_ic_harm, ns_count_harm = self.evaluate_by_namespace(harmonized)

            improvement = ((ic_harm - ic_raw) / ic_raw * 100) if ic_raw > 0 else 0.0

            results.append({
                'Protein': prot,
                'Raw_Count': len(raw),
                'Harmonized_Count': len(harmonized),
                'IC_Raw': round(ic_raw, 4),
                'IC_Harmonized': round(ic_harm, 4),
                'Improvement_Pct': round(improvement, 2),
                'BP_Raw': ns_count_raw['BP'],
                'BP_Harmonized': ns_count_harm['BP'],
                'MF_Raw': ns_count_raw['MF'],
                'MF_Harmonized': ns_count_harm['MF'],
                'CC_Raw': ns_count_raw['CC'],
                'CC_Harmonized': ns_count_harm['CC'],
                'IC_BP_Raw': round(ns_ic_raw['BP'], 4),
                'IC_BP_Harmonized': round(ns_ic_harm['BP'], 4),
                'IC_MF_Raw': round(ns_ic_raw['MF'], 4),
                'IC_MF_Harmonized': round(ns_ic_harm['MF'], 4),
                'IC_CC_Raw': round(ns_ic_raw['CC'], 4),
                'IC_CC_Harmonized': round(ns_ic_harm['CC'], 4),
            })

        print("Experiment complete.")
        return pd.DataFrame(results)

    def get_evidence_breakdown(self, protein_ids):
        """
        Analyze harmonization impact grouped by evidence code category.
        For each protein and each evidence category, harmonize only the terms
        from that category and measure IC improvement.
        """
        results = []
        for prot in protein_ids:
            evidence_df = self.get_protein_evidence(prot)

            for category, codes in self.EVIDENCE_CATEGORIES.items():
                cat_terms = set(evidence_df[evidence_df['Evidence'].isin(codes)]['GO_ID'])
                if not cat_terms:
                    continue

                # Harmonize just these terms
                harmonized = set(cat_terms)
                for go_id in cat_terms:
                    if go_id in self.godag:
                        harmonized.update(self.godag[go_id].get_all_parents())

                ic_raw = self.evaluate_informativeness(cat_terms)
                ic_harm = self.evaluate_informativeness(harmonized)
                improvement = ((ic_harm - ic_raw) / ic_raw * 100) if ic_raw > 0 else 0.0

                results.append({
                    'Protein': prot,
                    'Evidence_Category': category,
                    'Raw_Count': len(cat_terms),
                    'Harmonized_Count': len(harmonized),
                    'IC_Raw': round(ic_raw, 4),
                    'IC_Harmonized': round(ic_harm, 4),
                    'Improvement_Pct': round(improvement, 2),
                })

        return pd.DataFrame(results)
