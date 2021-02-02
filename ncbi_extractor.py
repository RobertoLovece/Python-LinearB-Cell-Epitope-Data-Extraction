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

    processed_arrays = []
    results = []
    db_return_dict = {}

    count = 0
    cooldown = 0

    for protein_id in unique:
        try:
            # after 200 requests take a break
            if (cooldown == 200):
                time.sleep(5)
                cooldown = 0
            db_return = Entrez.efetch(db="protein", id=protein_id, retmode="xml", rettype="fasta")
            db_return_dict[protein_id] = db_return
            print("protein_id '" + protein_id + "' (" + str(count) + ") returned " + str(db_return))
            count += 1
            cooldown += 1
        except:
            print("protein_id '" + protein_id + "' (" + str(count) + ") Bad Request - Protein_Id was not found")
            count += 1
            cooldown += 1

    count = 0

    for key in db_return_dict:
        record = Entrez.read(db_return_dict[key])
        db_return_dict[key].close()

        final_array = np.append(process_dictionary(record), [key, "NCBI protein"])
        processed_arrays.append(final_array)
        print("processed protein_id " + str(key) + " (" + str(count) + ")")
        count+=1

    output_array = np.vstack(processed_arrays)
    df = pd.DataFrame.from_records(output_array,columns=["seqtype","accver","taxid","orgname","defline","length","sequence","sid","protein_id","DB"])
    df.to_csv("output/ncbi_extractor_output.csv",index=False)

    return df

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

    element = d[0]

    if (element.get('TSeq_seqtype') != None and len(str(element.get('TSeq_seqtype'))) != 0):
        seqtype = element.get('TSeq_seqtype')
    else:
        seqtype = "NA"

    if (element.get('TSeq_accver') != None and len(str(element.get('TSeq_accver'))) != 0):
        accver = element.get('TSeq_accver')
    else:
        accver = "NA"

    if (element.get('TSeq_taxid') != None and len(str(element.get('TSeq_taxid'))) != 0):
        taxid = element.get('TSeq_taxid')
    else:
        taxid = "NA"

    if (element.get('TSeq_orgname') != None and len(str(element.get('TSeq_orgname'))) != 0):
        orgname = element.get('TSeq_orgname')
    else:
        orgname = "NA"

    if (element.get('TSeq_defline') != None and len(str(element.get('TSeq_defline'))) != 0):
        defline = element.get('TSeq_defline')
    else:
        defline = "NA"

    if (element.get('TSeq_length') != None and len(str(element.get('TSeq_length'))) != 0):
        length = element.get('TSeq_length')
    else:
        length = "NA"

    if (element.get('TSeq_sequence') != None and len(str(element.get('TSeq_sequence'))) != 0):
        sequence = element.get('TSeq_sequence')
    else:
        sequence = "NA"

    if (element.get('TSeq_sid') != None and len(str(element.get('TSeq_sid'))) != 0):
        sid = element.get('TSeq_sid')
    else:
        sid = "NA"

    lst = np.array([seqtype, accver, taxid, orgname, defline, length, sequence, sid])
    #df = pd.DataFrame([lst], columns=["seqtype","accver","taxid","orgname","defline","length","sequence","sid"])

    return lst
