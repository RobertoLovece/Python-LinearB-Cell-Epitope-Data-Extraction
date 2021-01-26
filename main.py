#-------------------------------------
# file: main.py

import iedb_extractor as iedb
import ncbi_extractor as ncbi

import pandas as pd
import numpy as np

def main():
    df = iedb.iedb_extract()
    ncbi.ncbi_extract(df)

if __name__ == "__main__":
    main()
