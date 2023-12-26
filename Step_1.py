# This is a sample Python script.

##
import pandas as pd
import json
import xml.etree.ElementTree as ET
import dbfunctions
from importlib import reload
#-----
tree = ET.parse('../full database.xml')
root = tree.getroot()

##
table = pd.read_json('v11.11/manual_table_translated1.json')
table = table[table['substances'] != 'error1']
substances = table.substances
substances


# search in drugbank
output_list = list()
for entry in list(substances):
    output = list()
    print('------', entry)
    entry = entry.capitalize()
    node = dbfunctions.runit(entry, root)

    if node == None:
        output.append('error2')
    else:
        out = dbfunctions.find_enzymes(node, root)
        output = output + out
    output = list(set(output)) # no doubles!

    # if there is err in combination preperate, it is not usable
    if 'error2' in output:
        output = ['error2']
    output_list.append(output)

check_error2 = pd.Series(output_list).value_counts()
check_error2.to_excel('check_error2.xlsx')
results = pd.DataFrame({
    'substance': substances,
    'interactions': output_list
})
results



## Speichern!
filename = r'db_results.json'
results.to_json(filename,orient='index', force_ascii='False')
filename = r'db_results.xlsx'
results.to_excel(filename)





# Make matrix 
data = pd.read_json(r'db_results.json', orient='index', convert_axes=False, dtype=False)
data

def search(keyword, liste):
    '''
    :param keyword: the enzyme name, that has to be sorted out
    :param liste: the list of the identified interactions
    :return: a list with 1 & 0, which entry has shown a interaction with given enzyme (keyword)
    '''
    lkeyword = list()
    for elem in liste:
        if keyword in elem: lkeyword.append(1)
        else: lkeyword.append(0)
    return lkeyword
enzymes = ['2d6sub', '2c19sub', '2c9sub', '3a4sub', '1a2sub']
for i in enzymes:
    search_result = search(i, list(data['interactions']))
    data[i] = search_result
data

def classify(row):
    if 'error2' in row['interactions']:
        return 'not_classified'
    else:
        return 'classified'

# Add the new column "classification" based on the "interactions" column
data['db_classified'] = data.apply(classify, axis=1)

##
data.to_excel(r'v11.11/db_results_matrix.xlsx')
data.to_json(r'v11.11/db_results_matrix.json',orient='index', force_ascii='False')
data
