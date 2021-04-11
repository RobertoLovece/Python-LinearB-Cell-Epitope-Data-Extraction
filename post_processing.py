#-------------------------------------
# file: post_processing.py

import pandas as pd
import numpy as np

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
    return df

def check(row):
    try:
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
    except:
        return False

def windowing(df):
    window_left = 7
    window_right = 7
    windowed_count = 0

    row_list = []

    no_dominant_class = []

    for row in zip(df['protein_id'], df['epitope_id'], df['epit_seq'], df['start_pos'], df['end_pos'], df['majority_class'], df['sequence']):

        protein_id = str(row[0])
        epitope_id = str(row[1])
        epit_seq = str(row[2])
        start_pos = int(row[3])
        end_pos = int(row[4])
        majority_class = str(row[5])
        sequence = np.asarray(list(str(row[6])))

        #print(epit_seq + " = " + str(len(epit_seq)))

        if (majority_class != "No dominant class"):
            windowed_count += 1
            for pos in range(start_pos-1,end_pos):
                dict = {}
                dict["protein_id"] = protein_id
                dict["epitope_id"] = epitope_id

                AA_pos = pos + 1
                offset_left = 0
                offset_right = 0

                if (AA_pos - window_left <= 0):
                    offset_left = -(AA_pos - window_left) + 1

                if (AA_pos + window_right > len(sequence)):
                    offset_right = AA_pos + window_right - len(sequence)

                AA_window_list = list(sequence[offset_left+pos-window_left:pos+window_right+1])


                if (offset_left != 0):
                    padded_value = sequence[0]

                    padding_left = [padded_value for x in range(0, offset_left)]
                    AA_window_list = padding_left + AA_window_list

                if (offset_right != 0):
                    padded_value = sequence[len(sequence)-1]

                    padding_right = [padded_value for x in range(len(sequence), AA_pos + window_right)]
                    AA_window_list = AA_window_list + padding_right

                AA_window_string = ''.join(AA_window_list)

                dict["AA_position"] = AA_pos
                dict["AA_window"] = AA_window_string

                dict["class"] = majority_class

                row_list.append(dict)
        else:
            dict = {}
            dict["protein_id"] = protein_id
            dict["epitope_id"] = epitope_id

            no_dominant_class.append(dict)

    print("Windowed Epitopes (" + str(windowed_count) + ")")
    print("Epitopes Left Out Of Windowing Due To No Dominant Class (" + str(len(no_dominant_class)) + ")")
    df = pd.DataFrame(no_dominant_class)
    df.to_csv("output/no_dominant_class.csv", index = False)

    df = pd.DataFrame(row_list, columns=['protein_id', 'epitope_id', 'AA_position', 'AA_window', 'class'])
    df.to_csv("output/windowed.csv", index = False)
