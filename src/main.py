from .data_loader import BioDataLoader

def main():
    # Шляхи відносно кореня проєкту (запускати через python main.py з кореня)
    obo_file = "data/raw/go-basic.obo"
    gaf_file = "data/raw/goa_human_plus.gaf"

    loader = BioDataLoader(obo_file, gaf_file)
    
    dag = loader.load_ontology()
    assocs = loader.load_annotations()
    
    df = loader.get_df()
    print("\nПерші 5 записів даних:")
    print(df.head())

if __name__ == "__main__":
    main()