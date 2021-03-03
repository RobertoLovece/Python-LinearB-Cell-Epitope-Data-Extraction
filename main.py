#-------------------------------------
# file: main.py

import iedb_extractor as iedb
import ncbi_extractor as ncbi

import pandas as pd
import numpy as np

def main():
    print("Extacting IEDB")
    iedb_df = iedb.iedb_extract()
    print()
    print("Finished Extacting IEDB")
    print()
    print("Extacting NCBI")
    ncbi_df = ncbi.ncbi_extract(iedb_df)

    df = iedb_df.merge(ncbi_df, how='left', on='protein_id')
    df.to_csv("output/left_join.csv", index = False)
    #df = pd.read_csv("output/left_join.csv")
    check_sequencing(df)
    windowing(df)

def check_sequencing(df):
    passed = 0
    failed = 0

    to_remove = []
    index = 0

    for row in zip(df['protein_id'], df['epit_seq'], df['start_pos'], df['end_pos'], df['sequence']):
        if (check(row) == True):
            passed += 1
        else:
            failed += 1
            to_remove.append(index)

        index += 1

    print("Passed Epitope To Sequence Matching (" + str(passed) + ")")
    print("Failed Epitope To Sequence Matching (" + str(failed) + ")")

    if (len(to_remove) != 0):
        print("Removing Rows Of Data At Indexes " + str(to_remove))
        df.drop(df.index[to_remove], inplace=True)

    df.to_csv("output/sequence_check.csv", index = False)

def check(row):
    epit_seq = str(row[1])
    start_pos = int(row[2])
    end_pos = int(row[3])
    sequence = np.asarray(list(str(row[4])))

    sample_list = sequence[start_pos-1:end_pos]
    sample_string = ''.join(sample_list)

    if (sample_string == epit_seq):
        return True
    else:
        return False

def windowing(df):
    window_left = 7
    window_right = 7



    for row in zip(df['protein_id'], df['epitope_id'], df['epit_seq'], df['start_pos'], df['end_pos'], df['sequence']):

        epit_seq = str(row[2])
        start_pos = int(row[3])
        end_pos = int(row[4])
        sequence = np.asarray(list(str(row[5])))

        check = ""
        print(epit_seq)

        for x in range(start_pos-1,end_pos):
            aa_pos = x + 1
            test = sequence[x-window_left:x+window_right+1]
            print(sequence[x])
            print(str(x) + str(test))

if __name__ == "__main__":
    main()
