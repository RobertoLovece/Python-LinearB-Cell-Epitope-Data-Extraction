#-------------------------------------
# file: iedb_extractor.py

import pandas as pd
import xml.etree.ElementTree as et

root = et.parse('285.xml').getroot()

path = "{http://www.iedb.org/schema/CurationSchema}"

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

    # the paths to all the parameter that need processing
    pubmed_id_path = path +"Reference/"+ path + "Article/" + path + "PubmedId"
    year_path = path +"Reference/"+ path + "Article/" + path + "ArticleYear"
    epit_name_path = path +"Reference/"+ path + "Epitopes/" + path + "Epitope/" + "EpitopeName"
    epitope_id_path = path +"Reference/"+ path + "Epitopes/" + path + "Epitope/" + "EpitopeId"
    evid_code_path = path +"Reference/"+ path + "Epitopes/" + path + "Epitope/" + "EpitopeEvidenceCode"
    epit_struc_def_path = path +"Reference/"+ path + "Epitopes/" + path + "Epitope/" + "EpitopeStructureDefines"
    sourceorg_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"SourceOrganismId"
    protein_id_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"SourceMolecule/"+ path +"GenBankId"
    epit_seq_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"LinearSequence"
    start_pos_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"StartingPosition"
    end_pos_path = path +"Reference/"+ path +"Epitopes/"+ path +"Epitope/"+ path +"EpitopeStructure/"+ path +"FragmentOfANaturalSequenceMolecule/"+ path +"EndingPosition"




print(check_linear_bcell_epitopes())
