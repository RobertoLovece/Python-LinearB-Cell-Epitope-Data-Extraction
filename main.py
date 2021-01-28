#-------------------------------------
# file: main.py

import iedb_extractor as iedb
import ncbi_extractor as ncbi

import pandas as pd
import numpy as np

def main():
    print("Extacting IEDB")
    df = iedb.iedb_extract()
    print("")
    print("Finished Extacting IEDB")
    print("")
    print("Extacting NCBI")
    ncbi.ncbi_extract(df)

if __name__ == "__main__":
    main()
