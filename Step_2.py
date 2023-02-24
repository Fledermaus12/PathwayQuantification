##
import pandas as pd
import xml.etree.ElementTree as ET
import dbfunctions
from importlib import reload
#-----
tree = ET.parse('full database.xml')
root = tree.getroot()

data = pd.read_json(r'base_output/drc_p1_pv2.json', orient='index', convert_axes=False, dtype=False)

##
s = data.AllInteraction.str.contains('err1', regex=False)
err_data = data.loc[s].copy()

mask = err_data.loc[:,'Z_Wirkstoff_Gesamt'].str.rstrip("/")
    # alle / hinten entfernen, weil es sonst leere Medikamente
    # produziert
mask = mask.str.split(pat="/")
err_data.loc[:, 'Z_Wirkstoff_Gesamt'] = mask

err_data_without_substance = err_data.loc[pd.isna(err_data['Z_Wirkstoff_Gesamt'])]
# --> tested: every row has at least substance

##
output_list = list()
for index, row in err_data.iterrows():
    output = list()
    for substance in row['Z_Wirkstoff_Gesamt']:
        node = dbfunctions.runit(substance, root)
        if node == None:
            substance_mod = substance + 'e'
            print('------', substance_mod)
            node = dbfunctions.runit(substance_mod, root)

        if node == None:
            output.append('err1')
        else:
            out = dbfunctions.find_enzymes(node, root)
            output = output + out
    output = list(set(output)) # no doubles!

    # if there is err in combination preperate, it is not usable
    if 'err1' in output:
        output = ['err1']
    output_list.append(output)


##
s = data.AllInteraction.str.contains('err1', regex=False)
err_data = data.loc[s].copy()
err_data['AllInteraction'] = output_list

##
data.update(err_data)
s = data.AllInteraction.str.contains('err1', regex=False)
err_data_leftover = data.loc[s].copy()
##
filename = r'base_output/drc_p2_pv2.json'
data.to_json(filename,orient='index', force_ascii='False')
filename = r'base_output/drc_p2_pv2.xlsx'
data.to_excel(filename)

