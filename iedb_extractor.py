#-------------------------------------
# file: iedb_extractor.py

import pandas as pd
import numpy as np
import xml.etree.ElementTree as et

root = et.parse('285.xml').getroot()

path = "{http://www.iedb.org/schema/CurationSchema}"

def main():
    if (check_linear_bcell_epitopes() == True):
        readAll()
        epitope_array = process_epitopes()
        assay_array = process_assays()

        processed_epitope_array = np.append(epitope_array,assay_array)

        df = pd.DataFrame([processed_epitope_array],columns=["pubmed_id","year","epit_name","epitope_id","evid_code","epit_struc_def","sourceOrg_id","protein_id","epit_seq","start_pos","end_pos","host_id","bcell_id","class","assay_type"])
        print(df)

# find out if it's a Linear B-Cell epitopes
def check_linear_bcell_epitopes():

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

def process_epitopes():

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
    if(root.find(pubmed_id_path) != None):
        pubmed_id = root.find(pubmed_id_path).text
    else:
        pubmed_id = "NA"

    if(root.find(year_path) != None):
        year = root.find(year_path).text
    else:
        year = "NA"

    if(root.find(epit_name_path) != None):
        epit_name = root.find(epit_name_path).text
    else:
        epit_name = "NA"

    if(root.find(epitope_id_path) != None):
        epitope_id = root.find(epitope_id_path).text
    else:
        epitope_id = "NA"

    if(root.find(evid_code_path) != None):
        evid_code = root.find(evid_code_path).text
    else:
        evid_code = "NA"

    if(root.find(epit_struc_def_path) != None):
        epit_struc_def = root.find(epit_struc_def_path).text
    else:
        epit_struc_def = "NA"

    if(root.find(sourceorg_id_path) != None):
        sourceorg_id = root.find(sourceorg_id_path).text
    else:
        sourceorg_id = "NA"

    if(root.find(protein_id_path) != None):
        protein_id = root.find(protein_id_path).text
    else:
        protein_id = "NA"

    if(root.find(epit_seq_path) != None):
        epit_seq = root.find(epit_seq_path).text
    else:
        epit_seq = "NA"

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

def process_assays():

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
