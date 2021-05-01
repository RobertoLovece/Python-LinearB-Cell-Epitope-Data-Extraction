#-------------------------------------
# file: main.py

import pandas as pd
import numpy as np

import iedb_extractor as iedb
import ncbi_extractor as ncbi
import post_processing as post

def main():
    print("Extacting IEDB")
    iedb_df = iedb.iedb_extract('input/example_xml/*.xml')
    print("Finished Extacting IEDB")

    print("Extacting NCBI")
    ncbi_df = ncbi.ncbi_extract(iedb_df)
    print("Finished Extacting NCBI")

    # left join the iedb_df to the ncbi_df
    df = iedb_df.merge(ncbi_df, how='left', on='protein_id')
    df.to_csv("output/left_join.csv", index = False)
    print("Checking Sequencing...")
    df = post.check_sequencing(df)
    print("Sequencing Finished")

    print("Starting Windowing Process...")
    post.windowing(df)
    print("Finished Windowing Process")

if __name__ == "__main__":
    main()
