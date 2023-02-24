# This is a sample Python script.

##
import pandas as pd
import json
import xml.etree.ElementTree as ET
import dbfunctions
from importlib import reload
#-----
tree = ET.parse('full database.xml')
root = tree.getroot()

##
table = pd.read_excel("data/tbl_Medikation.xlsx")
table = table.loc[:,['case_line_id','case_id','Z_Wirkstoff_Gesamt', 'Z_ATC_Gesamt']]
table = table.set_index('case_line_id')

# filtered the ones out without ATC, checked and proved
err_table = table.loc[pd.isna(table['Z_ATC_Gesamt'])]
#table.dropna(subset=['Z_ATC_Gesamt'], inplace=True)

##
output_list = list()
for index, row in table.iterrows():

    output = list()
    # ATC Feld darf nicht leer sein
    if pd.isna(row['Z_ATC_Gesamt']) == False:
        # Node des jeweiligen Medikaments finden und abspeichern
        __, node_drug, single_error = dbfunctions.find_drug_root(row['Z_ATC_Gesamt'], root)
    else: single_error == 1

    if single_error == 0:
        for elem in node_drug:
            output_temp = dbfunctions.find_enzymes(elem, root)
            output = output + output_temp
        output = list(set(output))
    else:
        output.append('err1')
    output_list.append(output)

table['AllInteraction'] = output_list


## Speichern!
filename = r'base_output/drc_p1_pv2.json'
table.to_json(filename,orient='index', force_ascii='False')
filename = r'base_output/drc_p1_pv2.xlsx'
table.to_excel(filename)


