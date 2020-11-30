#-------------------------------------
# file: iedb_extractor.py

import pandas as pd
import numpy as np
import xml.etree.ElementTree as et
import glob
import os

path = "{http://www.iedb.org/schema/CurationSchema}"
xml_path = 'example_XML/*.xml'
file_path = 'output.csv'

def main():
    # checks if output file already exists if so deletes it otherwise it will append to an existing file messing up the output
    if (os.path.exists(file_path)):
        os.remove(file_path)

    # makes sure the first df added to the csv file outputs it's header
    file_empty = True

    # iterates through all xml files in xml_path directory
    for xml in glob.iglob(xml_path):
        try:
            root = et.parse(xml).getroot()
            print(xml)
            if (check_linear_bcell_epitopes(root) == True):
                final_array = None
                # only getting one epitope from each file

                # get the article information first
                article_array = get_article_information(root)
                epitope_array = process_epitope(root)
                assay_array = process_assays(root)

                # combines the epitope and assay data into one array
                final_array = np.append(article_array, epitope_array)
                final_array = np.append(final_array, assay_array)

                # 15 parameters need to get passed to the DataFrame
                df = pd.DataFrame([final_array],columns=["pubmed_id","year","epit_name","epitope_id","evid_code","epit_struc_def","sourceOrg_id","protein_id","epit_seq","start_pos","end_pos","host_id","bcell_id","class","assay_type"])
                # print(df)

                # if the file is empty adds a header otherwise appends it to csv file
                if (file_empty == True):
                    df.to_csv('output.csv', mode='a',index=False)
                    file_empty = False
                else:
                    df.to_csv('output.csv', mode='a',index=False, header=False)

        except et.ParseError:
            print("Parse error in "+ xml +",invalid xml format")

# finds out if a Linear B-Cell epitopes is present to process file
def check_linear_bcell_epitopes(root):

    bcell_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path + "Assays/"+ path +"BCell"
    linear_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule"

    bcell_array = []
    linear_array = []

    # checks if xml contains bcell tag and linear
    for bcell in root.findall(bcell_path):
        bcell_array.append(bcell)

    for linear in root.findall(linear_path):
        linear_array.append(linear)

    # if it has the tags return true otherwise false
    if (len(linear_array) != 0 and len(bcell_array) != 0):
        return True
    else:
        return False

# gets the article part of the xml file
def get_article_information(root):
        # part of the article section
        pubmed_id_path = path +"Reference/"+ path +"Article/"+ path +"PubmedId"
        year_path = path +"Reference/"+ path +"Article/"+ path +"ArticleYear"

        pubmed_id = process_epitope_data(root, pubmed_id_path)
        year = process_epitope_data(root, year_path)

        article_array = np.array([pubmed_id, year])

        return article_array

def process_epitope(root):

    # the paths to all the parameters that need processing
    # part of the the epitope section
    epit_name_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeName"
    epitope_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeId"
    evid_code_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeEvidenceCode"
    epit_struc_def_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructureDefines"
    sourceorg_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"SourceOrganismId"
    protein_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"SourceMolecule/"+ path +"GenBankId"
    epit_seq_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"LinearSequence"
    start_pos_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"StartingPosition"
    end_pos_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"EndingPosition"

    # processes the parameter that are required
    epit_name = process_epitope_data(root, epit_name_path)
    epitope_id = process_epitope_data(root, epitope_id_path)
    evid_code = process_epitope_data(root, evid_code_path)
    epit_struc_def = process_epitope_data(root, epit_struc_def_path)
    sourceorg_id = process_epitope_data(root, sourceorg_id_path)
    protein_id = process_epitope_data(root, protein_id_path)
    epit_seq = process_epitope_data(root, epit_seq_path)

    # If the starting or end positions tags are empty we check another location
    if(root.find(start_pos_path) != None):
        start_pos = root.find(start_pos_path).text
    else:
        start_pos_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"ReferenceStartingPosition"
        if(root.find(start_pos_path) != None):
            start_pos = root.find(start_pos_path).text
        else:
            start_pos = "NA"

    if(root.find(end_pos_path) != None):
        end_pos = root.find(end_pos_path).text
    else:
        end_pos_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"ReferenceEndingPosition"
        if(root.find(end_pos_path) != None):
            end_pos = root.find(end_pos_path).text
        else:
            end_pos = "NA"

    epitope_array = np.array([epit_name, epitope_id, evid_code, epit_struc_def, sourceorg_id, protein_id, epit_seq, start_pos, end_pos])

    return epitope_array

# get the data from each parameter
def process_epitope_data(root, path):
    if(root.find(path) != None):
        processed = root.find(path).text
    else:
        processed = "NA"
    return processed

# process the data from each assay section of the epitope
def process_assays(root):

    # paths for the assay parameters that need processing
    host_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"Immunization/"+ path +"HostOrganism/"+ path +"OrganismId"
    bcell_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"BCellId"
    class_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"QualitativeMeasurement"
    assay_type_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"AssayTypeId"

    # finds all data from each assay parameter as each assay can have multiple tags
    # then appends it to a string so it's outputted correctly
    loop_list = []

    for host_id in root.findall(host_id_path):
        loop_list.append('"'+ host_id.text +'"')
        host_ids = ','.join(loop_list)

    loop_list.clear()

    for bcell_id in root.findall(bcell_id_path):
        loop_list.append('"'+ bcell_id.text +'"')
        bcell_ids = ','.join(loop_list)

    loop_list.clear()

    for assay_class in root.findall(class_path):
        loop_list.append('"'+ assay_class.text +'"')
        classes = ','.join(loop_list)

    loop_list.clear()

    for assay_type in root.findall(assay_type_path):
        loop_list.append('"'+ assay_type.text +'"')
        assay_types = ','.join(loop_list)

    assay_array = np.array([host_ids, bcell_ids, classes, assay_types])

    return assay_array

if __name__ == "__main__":
    main()
