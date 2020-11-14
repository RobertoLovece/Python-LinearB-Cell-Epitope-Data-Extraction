#-------------------------------------
# file: iedb_extractor.py

import pandas as pd
import xml.etree.ElementTree as et

root = et.parse('285.xml').getroot()

# find out if it's a Linear B-Cell epitopes
def check_linear_bcell_epitopes():

    path_to_bcell = "{http://www.iedb.org/schema/CurationSchema}Reference/{http://www.iedb.org/schema/CurationSchema}Epitopes/{http://www.iedb.org/schema/CurationSchema}Epitope/{http://www.iedb.org/schema/CurationSchema}Assays/{http://www.iedb.org/schema/CurationSchema}BCell"
    path_to_linear = "{http://www.iedb.org/schema/CurationSchema}Reference/{http://www.iedb.org/schema/CurationSchema}Epitopes/{http://www.iedb.org/schema/CurationSchema}Epitope/{http://www.iedb.org/schema/CurationSchema}EpitopeStructure/{http://www.iedb.org/schema/CurationSchema}FragmentOfANaturalSequenceMolecule"

    bcell_array = []
    linear_array = []

    for bcell in root.findall(path_to_bcell):
        bcell_array.append(bcell)

    for linear in root.findall(path_to_linear):
        linear_array.append(linear)

    if (len(linear_array) != 0 and len(bcell_array) != 0):
        return True
    else:
        return False

print(check_linear_bcell_epitopes())
