#-------------------------------------
# file: iedb_extractor.py

import pandas as pd
import numpy as np
import xml.etree.ElementTree as et
import glob

path = "{http://www.iedb.org/schema/CurationSchema}"

def iedb_extract(xml_path):

    # array of all processed epitopes
    processed_array = []

    # iterates through all xml files in xml_path directory
    for xml in glob.iglob(xml_path):
        try:
            root = et.parse(xml).getroot()
            print("Reading " + str(xml) + " XML File")
            if (check_linear_bcell_epitopes(root) == True):
                # path to each epitope
                epitope_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope"

                # get the article information first
                article_dict = get_article_information(root)

                # go through all the epitopes in the xml file
                for epitope in root.findall(epitope_path):
                    epitope_dict = process_epitope(epitope)

                    # prints the epitope being processed if it has a name
                    if(epitope.find(path + "EpitopeName") != None):
                        print("Processing Epitope " + epitope.find(path + "EpitopeName").text + " From XMLA")

                    if (check_assay(epitope) == True):
                        assay_dict = process_assays(epitope)

                        # combines all dictonaries of data into one
                        final_dict = {**article_dict, **epitope_dict, **assay_dict}

                    # when there's no assay
                    else:
                        empty_dict = {'host_id': "NA",'bcell_id': "NA",'assay_class': "NA",'majority_class': "NA",'assay_type': "NA"}
                        final_dict = {**article_dict, **epitope_dict, **empty_dict}

                    if (check_add_to_dataframe(final_dict) == True):
                        processed_array.append(final_dict)

        except et.ParseError:
            print("Parse error in "+ xml +",invalid xml format")

    # turn the final data extracted into a dataframe
    if (len(processed_array) != 0):
        df = pd.DataFrame(processed_array)
        df.to_csv('output/iedb_extractor_output.csv',index=False)
        return df

# finds out if a Linear B-Cell epitopes is present to process file
def check_linear_bcell_epitopes(root):

    bcell_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path + "Assays/"+ path +"BCell"
    linear_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/" + path + "LinearSequence"

    # if it has the tags return true otherwise false
    if (len(root.findall(bcell_path)) != 0 and len(root.findall(linear_path)) != 0):
        return True
    else:
        return False

# gets the article part of the xml file
def get_article_information(root):
    # dictonary to return values with key as name of field they represent
    dict = {}

    # part of the article section
    pubmed_id_path = path +"Reference/"+ path +"Article/"+ path +"PubmedId"
    year_path = path +"Reference/"+ path +"Article/"+ path +"ArticleYear"

    # adds the values for pubmed_id and year with the respective names as
    # keys
    dict['pubmed_id'] = process_epitope_data(root, pubmed_id_path)
    dict['year'] = process_epitope_data(root, year_path)

    return dict

def process_epitope(epitope):

    # dictonary to return values with key as name of field they represent
    dict = {}

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

    # processes the parameter that are required and add them to dict
    dict['epit_name'] = process_epitope_data(epitope, epit_name_path)
    dict['epitope_id'] = process_epitope_data(epitope, epitope_id_path)
    dict['evid_code'] = process_epitope_data(epitope, evid_code_path)
    dict['epit_struc_def'] = process_epitope_data(epitope, epit_struc_def_path)
    dict['sourceOrg_id'] = process_epitope_data(epitope, sourceorg_id_path)
    dict['protein_id'] = process_epitope_data(epitope, protein_id_path)
    dict['epit_seq'] = process_epitope_data(epitope, epit_seq_path)

    # If the starting or end positions tags are empty we check another location
    if(epitope.find(start_pos_path) != None):
        dict['start_pos'] = epitope.find(start_pos_path).text
    else:
        start_pos_path = path +"ReferenceStartingPosition"

        if(epitope.find(start_pos_path) != None):
            dict['start_pos'] = epitope.find(start_pos_path).text
        else:
            dict['start_pos'] = "NA"

    if(epitope.find(end_pos_path) != None):
        dict['end_pos'] = epitope.find(end_pos_path).text
    else:
        end_pos_path = path +"ReferenceEndingPosition"
        if(epitope.find(end_pos_path) != None):
            dict['end_pos'] = epitope.find(end_pos_path).text
        else:
            dict['end_pos'] = "NA"

    return dict

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

def check_add_to_dataframe(dict):
    if(dict.get('protein_id') == "NA"):
        return False
    elif(dict.get('epit_seq') == "NA"):
        return False
    elif(dict.get('assay_type') == "NA"):
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

    # dictonary to return values with key as name of field they represent
    dict = {}

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
        dict['host_id'] = ','.join(loop_list)

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
        dict['bcell_id'] = ','.join(loop_list)

    loop_list.clear()

    class_array = []

    for assay_class in epitope.findall(class_path):

        if(assay_class.text == "Positive"):
            class_value = "1"
        else:
            class_value = "-1"

        class_array.append(class_value)

        if (len(epitope.findall(class_path)) == 1):
            loop_list.append('"'+ class_value +'"')
        elif (len(loop_list) == 0):
            loop_list.append('"'+ class_value)
        elif (len(loop_list) == len(epitope.findall(class_path))-1):
            loop_list.append(class_value + '"')
        else:
            loop_list.append(class_value)
        dict['assay_class'] = ','.join(loop_list)

    dict['majority_class'] = get_majority_class(class_array)

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
        dict['assay_type'] = ','.join(loop_list)

    return dict

def get_majority_class(array):
    positive_count = array.count("1")
    negative_count = array.count("-1")

    if (positive_count > negative_count):
        majority_class = "1"
    elif(negative_count > positive_count):
        majority_class = "-1"
    else:
        majority_class = "No dominant class"

    return majority_class
