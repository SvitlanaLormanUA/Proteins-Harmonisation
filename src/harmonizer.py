import networkx as nx
import numpy as np
import pandas as pd
import math
from collections import Counter

class GOHarmonizer:
    def __init__(self, godag, associations_df):
        """
        godag: об'єкт GODag з goatools
        associations_df: DataFrame з колонками ['DB_Object_ID', 'GO_ID']
        """
        self.godag = godag
        self.df = associations_df
        # Підрахунок частоти термінів для Information Content (IC)
        self.term_counts = self._calculate_term_counts()
        self.total_annotations = len(self.df)
        
    def _calculate_term_counts(self):
        """Рахує, скільки разів кожен GO термін зустрічається в базі"""
        print("Calculating term frequencies...")
        counts = Counter(self.df['GO_ID'])
        
        # Важливо: Треба пропагувати каунти вгору. 
        # Якщо термін зустрічається 1 раз, його батько теж +1 раз.
        # Для спрощення поки беремо прямі входження, але для точності:
        # (Це місце можна покращити для розділу "Удосконалення")
        return counts

    def get_information_content(self, go_id):
        """
        Розрахунок Information Content (IC) за формулою Шеннона.
        IC(t) = -log(P(t)), де P(t) - ймовірність появи терміну.
        Чим рідкісніший (специфічніший) термін, тим вище IC.
        """
        if go_id not in self.term_counts:
            return 0.0
        
        prob = self.term_counts[go_id] / self.total_annotations
        return -math.log(prob)

    def get_protein_terms(self, protein_id):
        """Повертає список GO_ID для конкретного білка"""
        return self.df[self.df['DB_Object_ID'] == protein_id]['GO_ID'].tolist()

    def harmonize_protein(self, protein_id):
        """
        Основний алгоритм гармонізації (True Path Rule).
        Вхід: "Сирі" анотації білка.
        Вихід: Розширений набір анотацій (включаючи всіх предків).
        """
        original_terms = set(self.get_protein_terms(protein_id))
        harmonized_terms = set(original_terms)
        
        for go_id in original_terms:
            if go_id not in self.godag:
                continue
                
            term = self.godag[go_id]
            # Отримуємо всіх батьків рекурсивно
            parents = term.get_all_parents()
            harmonized_terms.update(parents)
            
        return list(original_terms), list(harmonized_terms)

    def evaluate_informativeness(self, terms_list):
        """
        Оцінка "сили" анотації.
        Сумуємо IC всіх термінів.
        """
        total_ic = sum(self.get_information_content(term) for term in terms_list)
        return total_ic

    def run_experiment(self, sample_size=5):
        """Демонстрація для курсової: беремо 5 білків і порівнюємо До/Після"""
        unique_proteins = self.df['DB_Object_ID'].unique()
        # Беремо випадкові білки
        sample_proteins = np.random.choice(unique_proteins, sample_size, replace=False)
        
        results = []
        for prot in sample_proteins:
            raw, harmonized = self.harmonize_protein(prot)
            
            ic_raw = self.evaluate_informativeness(raw)
            ic_harm = self.evaluate_informativeness(harmonized)
            
            results.append({
                "Protein": prot,
                "Raw_Count": len(raw),
                "Harmonized_Count": len(harmonized),
                "IC_Raw": round(ic_raw, 2),
                "IC_Harmonized": round(ic_harm, 2),
                "Improvement": round(((ic_harm - ic_raw) / ic_raw * 100), 1) if ic_raw > 0 else 0
            })
            
        return pd.DataFrame(results)