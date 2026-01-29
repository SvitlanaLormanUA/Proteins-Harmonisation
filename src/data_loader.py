import os
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader
import pandas as pd

class BioDataLoader:
    def __init__(self, obo_path, gaf_path):
        self.obo_path = obo_path
        self.gaf_path = gaf_path
        self.godag = None
        self.associations = None

    def load_ontology(self):
        """Завантажує структуру Gene Ontology (DAG)"""
        if not os.path.exists(self.obo_path):
            raise FileNotFoundError(f"Файл онтології не знайдено: {self.obo_path}")
        
        print(f"Loading Ontology from {self.obo_path}...")
        self.godag = GODag(self.obo_path)
        print(f"Ontology loaded. Terms count: {len(self.godag)}")
        return self.godag

    def load_annotations(self):
        """Завантажує анотації (GAF файл) і перетворює їх у зручний формат"""
        if not os.path.exists(self.gaf_path):
            raise FileNotFoundError(f"Файл анотацій не знайдено: {self.gaf_path}")

        print(f"Loading Annotations from {self.gaf_path}...")
        ogaf = GafReader(self.gaf_path)
        self.associations = ogaf.associations
        
        print(f"Annotations loaded. Associations found: {len(self.associations)}")
        return self.associations

    def get_df(self):
        """Повертає DataFrame для зручного перегляду даних"""
        data = []
        for assoc in self.associations:
            data.append({
                'DB_Object_ID': assoc.DB_ID, # ID білка (UniProt)
                'GO_ID': assoc.GO_ID,               # ID терміну GO
                'Evidence': assoc.Evidence_Code     # Код доказовості
            })
        return pd.DataFrame(data)

if __name__ == "__main__":
    loader = BioDataLoader(
        obo_path="data/raw/go-basic.obo",
        gaf_path="data/raw/goa_human_plus.gaf"
    )
    loader.load_ontology()