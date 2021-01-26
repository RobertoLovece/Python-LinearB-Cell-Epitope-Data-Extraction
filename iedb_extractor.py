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

def iedb_extract():

    # array of all processed epitopes
    processed_array = []

    # iterates through all xml files in xml_path directory
    for xml in glob.iglob(xml_path):
        try:
            root = et.parse(xml).getroot()
            print(xml)
            if (check_linear_bcell_epitopes(root) == True):
                # path to each epitope
                epitope_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope"

                # final output of each bit of data
                final_array = None

                # get the article information first
                article_array = get_article_information(root)

                # go through all the epitopes in the xml file
                for epitope in root.findall(epitope_path):
                    epitope_array = process_epitope(epitope)

                    if(epitope.find(path + "EpitopeName") != None):
                        print(epitope.find(path + "EpitopeName").text)

                    if (check_assay(epitope) == True):
                        assay_array = process_assays(epitope)

                        # combines the epitope and assay data into one array
                        final_array = np.append(article_array, epitope_array)
                        final_array = np.append(final_array, assay_array)
                    else:
                        empty_array = np.array(["NA","NA","NA","NA"])
                        final_array = np.append(article_array, epitope_array)
                        final_array = np.append(final_array, empty_array)

                    # 15 parameters need to get passed to the DataFrame
                    df = pd.DataFrame([final_array],columns=["pubmed_id","year","epit_name","epitope_id","evid_code","epit_struc_def","sourceOrg_id","protein_id","epit_seq","start_pos","end_pos","host_id","bcell_id","assay_class","assay_type"])

                    if (add_to_dataframe(df) == True):
                        processed_array.append(final_array)

        except et.ParseError:
            print("Parse error in "+ xml +",invalid xml format")

    if (len(processed_array) != 0):
        output_array = np.vstack(processed_array)
        df = pd.DataFrame.from_records(output_array,columns=["pubmed_id","year","epit_name","epitope_id","evid_code","epit_struc_def","sourceOrg_id","protein_id","epit_seq","start_pos","end_pos","host_id","bcell_id","assay_class","assay_type"])
        df.to_csv(file_path,index=False)
        return df

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

def process_epitope(epitope):
    # the paths to all the parameters that need processing
    # part of the the epitope section
    epit_name_path = path +"EpitopeName"
    epitope_id_path = path +"EpitopeId"
    evid_code_path = path +"EpitopeEvidenceCode"
    epit_struc_def_path = path +"EpitopeStructureDefines"
    sourceorg_id_path = path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"SourceOrganismId"
    protein_id_path = path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"SourceMolecule/"+ path +"GenBankId"
    epit_seq_path = path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"LinearSequence"
    start_pos_path = path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"StartingPosition"
    end_pos_path = path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"EndingPosition"

    # processes the parameter that are required
    epit_name = process_epitope_data(epitope, epit_name_path)
    epitope_id = process_epitope_data(epitope, epitope_id_path)
    evid_code = process_epitope_data(epitope, evid_code_path)
    epit_struc_def = process_epitope_data(epitope, epit_struc_def_path)
    sourceorg_id = process_epitope_data(epitope, sourceorg_id_path)
    protein_id = process_epitope_data(epitope, protein_id_path)
    epit_seq = process_epitope_data(epitope, epit_seq_path)

    # If the starting or end positions tags are empty we check another location
    if(epitope.find(start_pos_path) != None):
        start_pos = epitope.find(start_pos_path).text
    else:
        start_pos_path = path +"ReferenceStartingPosition"
        if(epitope.find(start_pos_path) != None):
            start_pos = epitope.find(start_pos_path).text
        else:
            start_pos = "NA"

    if(epitope.find(end_pos_path) != None):
        end_pos = epitope.find(end_pos_path).text
    else:
        end_pos_path = path +"ReferenceEndingPosition"
        if(epitope.find(end_pos_path) != None):
            end_pos = epitope.find(end_pos_path).text
        else:
            end_pos = "NA"

    epitope_array = np.array([epit_name, epitope_id, evid_code, epit_struc_def, sourceorg_id, protein_id, epit_seq, start_pos, end_pos])

    return epitope_array

# check if there is an array present in the epitope
def check_assay(epitope):
    # paths to assay stuff
    assay_path = path + "Assays"
    host_id_path = path +"Assays/"+ path +"BCell/"+ path +"Immunization/"+ path +"HostOrganism/"+ path +"OrganismId"
    bcell_id_path = path +"Assays/"+ path +"BCell/"+ path +"BCellId"
    class_path = path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"QualitativeMeasurement"
    assay_type_path = path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"AssayTypeId"

    # checks if the parametes I need are in the assay
    if(epitope.find(assay_path) == None):
        return False
    else:
        if (len(epitope.findall(host_id_path)) != 0):
            return True
        elif (len(epitope.findall(bcell_id_path)) != 0):
            return True
        elif (len(epitope.findall(class_path)) != 0):
            return True
        elif (len(epitope.findall(assay_type_path)) != 0):
            return True
        else:
            return False

def add_to_dataframe(df):
    if(df.get('protein_id')[0] == "NA"):
        return False
    elif(df.get('epit_seq')[0] == "NA"):
        return False
    elif(df.get('assay_type')[0] == "NA"):
        return False
    else:
        return True

# get the data from each parameter
def process_epitope_data(root, path):
    if(root.find(path) != None):
        processed = root.find(path).text
    else:
        processed = "NA"
    return processed

# process the data from each assay section of the epitope
def process_assays(epitope):
    # paths for the assay parameters that need processing
    host_id_path = path +"Assays/"+ path +"BCell/"+ path +"Immunization/"+ path +"HostOrganism/"+ path +"OrganismId"
    bcell_id_path = path +"Assays/"+ path +"BCell/"+ path +"BCellId"
    class_path = path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"QualitativeMeasurement"
    assay_type_path = path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"AssayTypeId"

    # finds all data from each assay parameter as each assay can have multiple tags
    # then appends it to a string so it's outputted correctly
    loop_list = []


    for host_id in epitope.findall(host_id_path):
        if (len(epitope.findall(host_id_path)) == 1):
            loop_list.append('"'+ host_id.text + '"')
        elif (len(loop_list) == 0):
            loop_list.append('"'+ host_id.text)
        elif (len(loop_list) == len(epitope.findall(host_id_path))-1):
            loop_list.append(host_id.text + '"')
        else:
            loop_list.append(host_id.text)
        host_ids = ','.join(loop_list)

    loop_list.clear()

    for bcell_id in epitope.findall(bcell_id_path):
        if (len(epitope.findall(bcell_id_path)) == 1):
            loop_list.append('"'+ bcell_id.text + '"')
        elif (len(loop_list) == 0):
            loop_list.append('"'+ bcell_id.text)
        elif (len(loop_list) == len(epitope.findall(bcell_id_path))-1):
            loop_list.append(bcell_id.text + '"')
        else:
            loop_list.append(bcell_id.text)
        bcell_ids = ','.join(loop_list)

    loop_list.clear()

    for assay_class in epitope.findall(class_path):

        if(assay_class.text == "Positive"):
            class_value = "1"
        else:
            class_value = "-1"

        if (len(epitope.findall(class_path)) == 1):
            loop_list.append('"'+ class_value +'"')
        elif (len(loop_list) == 0):
            loop_list.append('"'+ class_value)
        elif (len(loop_list) == len(epitope.findall(class_path))-1):
            loop_list.append(class_value + '"')
        else:
            loop_list.append(class_value)
        classes = ','.join(loop_list)

    loop_list.clear()

    for assay_type in epitope.findall(assay_type_path):
        if (len(epitope.findall(assay_type_path)) == 1):
            loop_list.append('"'+ assay_type.text +'"')
        elif (len(loop_list) == 0):
            loop_list.append('"' + assay_type.text)
        elif (len(loop_list) == len(epitope.findall(assay_type_path))-1):
            loop_list.append(assay_type.text + '"')
        else:
            loop_list.append(assay_type.text)
        assay_types = ','.join(loop_list)

    assay_array = np.array([host_ids, bcell_ids, classes, assay_types])

    return assay_array
