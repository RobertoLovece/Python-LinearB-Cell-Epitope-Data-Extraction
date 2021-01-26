import iedb_extractor as iedb
import pandas as pd
import numpy as np

def main():
    df = iedb.iedb_extract()
    print(df)

if __name__ == "__main__":
    main()
