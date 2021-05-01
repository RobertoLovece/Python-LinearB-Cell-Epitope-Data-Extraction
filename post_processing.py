#-------------------------------------
# file: post_processing.py

import pandas as pd
import numpy as np

# check the sequence of epitopes from the left join
def check_sequencing(df):
    passed = 0
    failed = 0

    to_remove = []
    index = 0

    # iterate through the dataframe and get these parameters
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
    # check if epitope sequence is inside the protein sequence where it should
    # be using start_pos and end_pos
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
    # windowing sizes
    window_left = 7
    window_right = 7
    # how many epitopes have been windowed
    windowed_count = 0

    # rows of values
    row_list = []

    # epitopes with no dominant class
    no_dominant_class = []

    # iterate through dataframe
    for row in zip(df['protein_id'], df['epitope_id'], df['epit_seq'], df['start_pos'], df['end_pos'], df['majority_class'], df['sequence']):

        protein_id = str(row[0])
        epitope_id = str(row[1])
        epit_seq = str(row[2])
        start_pos = int(row[3])
        end_pos = int(row[4])
        majority_class = str(row[5])
        sequence = np.asarray(list(str(row[6])))


        if (majority_class != "No dominant class"):
            windowed_count += 1
            for pos in range(start_pos-1,end_pos):
                dict = {}
                dict["protein_id"] = protein_id
                dict["epitope_id"] = epitope_id

                AA_pos = pos + 1
                offset_left = 0
                offset_right = 0

                # padding is required at the start of the sequence change the
                # offset left as not to exceed the sequence bounds
                if (AA_pos - window_left <= 0):
                    offset_left = -(AA_pos - window_left) + 1

                # if padding is required at the end of the sequence change the
                # offset right as not to exceed the sequence bounds
                if (AA_pos + window_right > len(sequence)):
                    offset_right = AA_pos + window_right - len(sequence)

                # get the area of the sequence without exceed the bounds
                AA_window_list = list(sequence[offset_left+pos-window_left:pos+window_right+1])


                # if padding is required left
                if (offset_left != 0):
                    # get the value the padding will be
                    padded_value = sequence[0]

                    # get the amount of values need to be padded
                    padding_left = [padded_value for x in range(0, offset_left)]
                    # add it to the AA_window
                    AA_window_list = padding_left + AA_window_list

                # if padding is required right
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
            # if the epitope has no dominant class
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

    return df
