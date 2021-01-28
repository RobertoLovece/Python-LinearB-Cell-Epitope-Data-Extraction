#-------------------------------------
# file: ncbi_extractor.py

import pandas as pd
import numpy as np
import eutils
import time

from Bio import Entrez
from Bio import SeqIO

def ncbi_extract(df):

    Entrez.email = "rlovece@hotmail.co.uk"
    Entrez.api_key = "ee238abcf0e852f21ef1b98fd0d440477a08"

    unique = df.protein_id.unique()
    print("unique protein_id amount = " + str(len(unique)))

    unique_df = pd.DataFrame.from_records(np.vstack(unique), columns=["protein_id"])
    unique_df.to_csv("output/unique_protein_ids.csv")

    processed_array = []
    results = []
    db_return = None

    queries = protein_ids_to_strings(unique)

    for query in queries:
        db_return = Entrez.efetch(db="protein", id="ABD64214.1", retmode="xml", rettype="fasta")
        #db_return = Entrez.efetch(db="protein", id=query, retmode="xml", rettype="fasta")
        results.append(db_return)
        db_return = None
        time.sleep(3)

    count = 0

    for result in results:
        record = Entrez.read(result)
        result.close()

        for d in record:
            extras = np.array([unique[count], "NCBI protein"])
            final_array = np.append(process_dictionary(d), extras)

            processed_array.append(final_array)
            print("processed protein_id " + str(unique[count]) + " (" + str(count) + ")")
            count+=1

    output_array = np.vstack(processed_array)
    df = pd.DataFrame.from_records(output_array,columns=["seqtype","accver","taxid","orgname","defline","length","sequence","sid","UID","DB"])
    df.to_csv("output/ncbi_extractor_output.csv",index=False)

def protein_ids_to_strings(unique):
    # splits a list into multiple different list of size equal to the final int
    # suggested limit is 200
    split = [unique[x:x+200] for x in range(0, len(unique), 200)]

    queries = []
    query = ""

    for list in split:
        query = ','.join(list)
        queries.append(query)
        query = ""

    return queries

def process_dictionary(d):

    if (d.get('TSeq_seqtype') != None and len(str(d.get('TSeq_seqtype'))) != 0):
        seqtype = d.get('TSeq_seqtype')
    else:
        seqtype = "NA"

    if (d.get('TSeq_accver') != None and len(str(d.get('TSeq_accver'))) != 0):
        accver = d.get('TSeq_accver')
    else:
        accver = "NA"

    if (d.get('TSeq_taxid') != None and len(str(d.get('TSeq_taxid'))) != 0):
        taxid = d.get('TSeq_taxid')
    else:
        taxid = "NA"

    if (d.get('TSeq_orgname') != None and len(str(d.get('TSeq_orgname'))) != 0):
        orgname = d.get('TSeq_orgname')
    else:
        orgname = "NA"

    if (d.get('TSeq_defline') != None and len(str(d.get('TSeq_defline'))) != 0):
        defline = d.get('TSeq_defline')
    else:
        defline = "NA"

    if (d.get('TSeq_length') != None and len(str(d.get('TSeq_length'))) != 0):
        length = d.get('TSeq_length')
    else:
        length = "NA"

    if (d.get('TSeq_sequence') != None and len(str(d.get('TSeq_sequence'))) != 0):
        sequence = d.get('TSeq_sequence')
    else:
        sequence = "NA"

    if (d.get('TSeq_sid') != None and len(str(d.get('TSeq_sid'))) != 0):
        sid = d.get('TSeq_sid')
    else:
        sid = "NA"

    lst = np.array([seqtype, accver, taxid, orgname, defline, length, sequence, sid])
    #df = pd.DataFrame([lst], columns=["seqtype","accver","taxid","orgname","defline","length","sequence","sid"])

    return lst
