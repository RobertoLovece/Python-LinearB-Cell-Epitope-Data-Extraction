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

    xml_files = glob.glob(xml_path)
    #checks if output file already exists if so deletes it otherwise it will append to an existing file messing up the output
    if (os.path.exists(file_path)):
        os.remove(file_path)

    file_empty = True

    for xml in xml_files:
        try:
            root = et.parse(xml).getroot()
            print(xml)
            if (check_linear_bcell_epitopes(root) == True):
                epitope_array = process_epitopes(root)
                assay_array = process_assays(root)

                epitope_array = np.append(epitope_array,assay_array)

                # 15 parameters need to get passed
                df = pd.DataFrame([epitope_array],columns=["pubmed_id","year","epit_name","epitope_id","evid_code","epit_struc_def","sourceOrg_id","protein_id","epit_seq","start_pos","end_pos","host_id","bcell_id","class","assay_type"])
                #print(df)
                if (file_empty == True):
                    df.to_csv('output.csv', mode='a',index=False)
                    file_empty = False
                else:
                    df.to_csv('output.csv', mode='a',index=False, header=False)

        except et.ParseError:
            print("Parse error in "+ xml +",invalid xml format")

# find out if it's a Linear B-Cell epitopes
def check_linear_bcell_epitopes(root):

    bcell_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path + "Assays/"+ path +"BCell"
    linear_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule"

    bcell_array = []
    linear_array = []

    for bcell in root.findall(bcell_path):
        bcell_array.append(bcell)

    for linear in root.findall(linear_path):
        linear_array.append(linear)

    if (len(linear_array) != 0 and len(bcell_array) != 0):
        return True
    else:
        return False

def process_epitopes(root):

    # the paths to all the parameters that need processing
    pubmed_id_path = path +"Reference/"+ path +"Article/"+ path +"PubmedId"
    year_path = path +"Reference/"+ path +"Article/"+ path +"ArticleYear"
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
    pubmed_id = process_epitope_data(root, pubmed_id_path)
    year = process_epitope_data(root, year_path)
    epit_name = process_epitope_data(root, epit_name_path)
    epitope_id = process_epitope_data(root, epitope_id_path)
    evid_code = process_epitope_data(root, evid_code_path)
    epit_struc_def = process_epitope_data(root, epit_struc_def_path)
    sourceorg_id = process_epitope_data(root, sourceorg_id_path)
    protein_id = process_epitope_data(root, protein_id_path)
    epit_seq = process_epitope_data(root, epit_seq_path)

    # If the starting or end positions are empty we check another location
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

    epitope_array = np.array([pubmed_id, year, epit_name, epitope_id, evid_code, epit_struc_def, sourceorg_id, protein_id, epit_seq, start_pos, end_pos])

    return epitope_array

def process_epitope_data(root, path):
    if(root.find(path) != None):
        processed = root.find(path).text
    else:
        processed = "NA"
    return processed

def process_assays(root):

    # paths for the assay parameters that need processing
    host_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"Immunization/"+ path +"HostOrganism/"+ path +"OrganismId"
    bcell_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"BCellId"
    class_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"QualitativeMeasurement"
    assay_type_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"Assays/"+ path +"BCell/"+ path +"AssayInformation/"+ path +"AssayTypeId"

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
