#-------------------------------------
# file: main.py

import iedb_extractor as iedb
import ncbi_extractor as ncbi

import pandas as pd
import numpy as np

def main():
    print("Extacting IEDB")
    iedb_df = iedb.iedb_extract()
    print("")
    print("Finished Extacting IEDB")
    print("")
    print("Extacting NCBI")
    ncbi_df = ncbi.ncbi_extract(iedb_df)

    df = iedb_df.merge(ncbi_df, how='left', on='protein_id')
    df.to_csv("output/left_join.csv", index = False)

def check_sequencing():
    

if __name__ == "__main__":
    main()
