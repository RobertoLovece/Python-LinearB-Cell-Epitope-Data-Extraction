#-------------------------------------
# file: ncbi_extractor.py

import pandas as pd
import numpy as np

from Bio import Entrez
from Bio import SeqIO

request_size = 200

def ncbi_extract(df):

    # initialize Entrez with details
    Entrez.email = "rlovece@hotmail.co.uk"
    Entrez.api_key = "ee238abcf0e852f21ef1b98fd0d440477a08"

    # get unique protien ids from iedb extraction
    unique = df.protein_id.unique()
    print("unique protein_id amount = " + str(len(unique)))

    # output unique protein ids to file
    unique_df = pd.DataFrame.from_records(np.vstack(unique), columns=["protein_id"])
    unique_df.to_csv("output/unique_protein_ids.csv")

    processed_array = []

    need_repeating = []
    protein_dict = {}

    ids = protein_ids_to_strings(unique)
    unique_split = [unique[x:x+request_size] for x in range(0, len(unique), request_size)]

    count = 0
    match_array = []
    accver_array = []

    while count < len(ids):
        id_request = ids[count]
        unique_item = unique_split[count]

        print("Processing queries " + str(count * request_size) + " - " + str(count * request_size + len(unique_item)))

        try:
            handle = Entrez.efetch(db="protein", id=id_request, retmode="xml", rettype="fasta")
            records = Entrez.read(handle)

            handle.close()

            if (len(records) == len(unique_item)):
                match_array = match_array + unique_item.tolist()
                index = 0
                for record in records:
                    final_dict = process_record(record)

                    final_dict["protein_id"] = unique_item[index]
                    final_dict["DB"] = "NCBI protein"

                    processed_array.append(final_dict)

                    index += 1

            else:
                for record in records:
                    accver = record.get('TSeq_accver')

                    # if the accver matches the protein_id add the record to a
                    # dictonary with the key as the protein_id
                    #print(accver + " = " + str(accver in unique_item))

                    if accver in unique_item:
                        final_dict = process_record(record)

                        accver_array.append(accver)
                        final_dict["protein_id"] = accver
                        final_dict["DB"] = "NCBI protein"

                        processed_array.append(final_dict)

            #need_repeating.extend(list(set(unique_item) - set(accver_array)))
        except:
            print("Bad Request - Query id's not found")
        count += 1

    need_repeating = (set(unique) - set(accver_array) - set(match_array))

    print("Repeating Failed " + str(len(need_repeating)) + " Queries")

    count = 1
    bad_requests = []

    for repeat in need_repeating:
        try:
            handle = Entrez.efetch(db="protein", id=repeat, retmode="xml", rettype="fasta")
            record = Entrez.read(handle)
            handle.close()

            final_dict = process_record(record[0])

            final_dict["protein_id"] = repeat
            final_dict["DB"] = "NCBI protein"

            processed_array.append(final_dict)

            print("(" + str(count) + "/" + str(len(need_repeating)) + ") Protein_Id '" + repeat + "' Returned " + str(handle))
        except:
            print("(" + str(count) + "/" + str(len(need_repeating)) + ") Protein_Id '" + repeat + "' Bad Request - Protein_Id was not found")
            bad_requests.append(repeat)
        count += 1

    bad = pd.DataFrame(bad_requests, columns = ["protein_id"])
    bad.to_csv("output/bad_requests.csv", index = False)

    df = pd.DataFrame.from_records(processed_array,columns = ["accver","taxid","orgname","defline","length","sequence","sid","protein_id","DB"])
    df.to_csv("output/ncbi_extractor_output.csv", index = False)

    return df

def protein_ids_to_strings(unique):
    # splits a list into multiple different list of size equal to the
    # request_size global variable
    split = [unique[x:x+request_size] for x in range(0, len(unique), request_size)]

    queries = []
    query = ""

    for list in split:
        query = ','.join(list)
        queries.append(query)
        query = ""

    return queries

def process_record(element):

    dict = {}

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

    dict["accver"] = accver
    dict["taxid"] = taxid
    dict["orgname"] = orgname
    dict["defline"] = defline
    dict["length"] = length
    dict["sequence"] = sequence
    dict["sid"] = sid

    #lst = np.array([accver, taxid, orgname, defline, length, sequence, sid])

    return dict
