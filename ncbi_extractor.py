#-------------------------------------
# file: ncbi_extractor.py

import pandas as pd
import numpy as np
import eutils
import time

def ncbi_extract(df):
    ec = eutils.Client()
    unique = df.protein_id.unique()
    #test = "P06821.1, P14013.1, CAA32579.1, P14916.2, 30173397, 229595253"
    #ec.efetch(db="protein", id=test)

    queries = protein_ids_to_strings(unique)
    results = []
    db_return = None

    for query in queries:
        db_return = ec.efetch(db="protein", id=query, retmode="xml", rettype="native")
        results.append(db_return)
        db_return = None
        time.sleep(3)
    

def protein_ids_to_strings(unique):
    # splits a list into multiple different list of size equal to the final int
    # suggested limit is 200
    split = [unique[x:x+200] for x in range(0, len(unique), 200)]
    #print(split)

    queries = []
    query = ""

    for list in split:
        query = ','.join(list)
        queries.append(query)
        query = ""

    return queries
