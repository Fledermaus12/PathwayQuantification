import pandas as pd 
import numpy as np
import math



############### IMPORTS #################




# Drugbank

drc = pd.read_json(r'v11.11/db_results_matrix.json', orient='index', convert_axes=False, dtype=False)
drc.drop(columns=['3a4sub', '1a2sub'], inplace=True)
drc = drc[~drc['substance'].str.contains('error1')]
drc = drc.drop_duplicates(subset='substance')


# Manual Table #

man = pd.read_json(r'v11.11/manual_table_translated1.json')
man = man[~man['substances'].str.contains('error1')]
# making sure only to keep the classified ones in case duplicates appear
classified_rows = man[man['man_classified'] == 'classified']
classified_unique = classified_rows.drop_duplicates(subset='substances')
unclassified_rows = man[man['man_classified'] == 'not_classified']
unclassified_unique = unclassified_rows.drop_duplicates(subset='substances')
man = pd.concat([classified_unique, unclassified_unique]).drop_duplicates(subset='substances', keep='first')

# set bisoprolol to 0 
s = 'Bisoprolol'
man.set_index('Z_Wirkstoff_Gesamt', inplace=True)
man.loc[s, 'CYP2D6'] = 0
man.reset_index(inplace=True)
man = man.copy()
del s

man = man[~man['substances'].str.contains('error1')]
duplicates = man[man['substances'].duplicated(keep=False)]
duplicates
man = man.drop_duplicates(subset='substances')

# Flockhardt 

flock = pd.read_excel('FlockhardtTable_converted.xlsx', index_col=0)
#flock = flock[['2C9', '2C19', '2D6']]


# FDA

fda = pd.read_excel('fda_table_clean_converted.xlsx', index_col=0)
fda = fda[['2d6sub', '2c19sub', '2c9sub']]
fda['fda_classified'] = 'classified'
fda




# Load the original Table
man_og_start = pd.read_excel("../data/CYP_Substrate_haendisch.xlsx")
man_og = man_og_start[['Z_Wirkstoff_Gesamt','CYP2D6.1', 'CYP2C19.1', 'CYP2C9.1']]
man_og = man_og.rename(columns={
    'CYP2D6.1':'CYP2D6', 
    'CYP2C19.1': 'CYP2C19', 
    'CYP2C9.1': 'CYP2C9'
})
man_og


# Count basic occurences
def count_occurence_basic(med):
    if man_og['Z_Wirkstoff_Gesamt'].str.contains(med, regex=False).any() == True:
        count = len(man_og[man_og['Z_Wirkstoff_Gesamt'].str.contains(med, regex=False)])
    else: 
        print(med, 'fail to find in table')
        count = 0
    return count

man['counts'] = man['Z_Wirkstoff_Gesamt'].apply(count_occurence_basic)
man.counts.sum()


'''


MAN & DRUGBANK FUSION


'''


# Fusion MANUELL and DRUGBANK
drc_reduced = drc.drop(columns='interactions')
drc_reduced.set_index('substance', inplace=True)
#man_reduced = man.drop(columns='Z_Wirkstoff_Gesamt')
man_reduced = man.copy()
man_reduced.set_index('substances', inplace=True)
df = drc_reduced.join(man_reduced, how='outer')
df = df.rename(columns={
    'index': 'substance',
    '2d6sub':'2d6sub_db',
    '2c19sub':'2c19sub_db',
    '2c9sub': '2c9sub_db',
    'CYP2D6': '2d6sub_man',
    'CYP2C19': '2c19sub_man',
    'CYP2C9': '2c9sub_man'})
df

# count, how many entries each substance causes in the dataset
check1 = df["db_classified"] == "not_db_classified"
check2 = df["man_classified"] == "not_classified"
df[check1].sort_values('counts', ascending=False)
print(df[check1].sort_values('counts', ascending=False)[:20])



'''


MAN & FDA FUSION


'''






man_reduced = man.copy()
man_reduced = man_reduced.rename(columns={
    'substances': 'substance',
    'CYP2D6': '2d6sub',
    'CYP2C19': '2c19sub',
    'CYP2C9': '2c9sub'})
man_reduced = man_reduced.set_index('substance')

df = man_reduced.join(fda, how='left', lsuffix='_man', rsuffix='_fda')
df = df.reset_index()
df = df.rename(columns={'index': 'substance'})
df = df.drop_duplicates(subset='substance')
df = df.set_index('substance')
df.fda_classified = df.fda_classified.fillna('not_classified')
df.man_classified = df.man_classified.fillna('not_classified')
df[['2d6sub_fda', '2c19sub_fda', '2c9sub_fda']] = df[['2d6sub_fda', '2c19sub_fda', '2c9sub_fda']].fillna(0)
df[['2d6sub_man', '2c19sub_man', '2c9sub_man']] = df[['2d6sub_man', '2c19sub_man', '2c9sub_man']].fillna(0)

check1 = df.man_classified == 'classified'
check2 = df.man_classified == 'not_classified'
check3 = df.fda_classified == 'classified'
check4 = df.fda_classified == 'not_classified'
data = df[check4].sort_values('counts', ascending=False)[:20]
print(data)




'''


MAN & FLOCKHART FUSION


'''



man_reduced = man.copy()
man_reduced = man_reduced.rename(columns={
    'substances': 'substance',
    'CYP2D6': '2d6sub',
    'CYP2C19': '2c19sub',
    'CYP2C9': '2c9sub'})
man_reduced = man_reduced.set_index('substance')

flock_copy = flock.set_index('drugs')
df = man_reduced.join(flock_copy, how='left', lsuffix='_man', rsuffix='_flock')
df = df.reset_index()
df = df.rename(columns={'index': 'substance'})
df = df.drop_duplicates(subset='substance')
df = df.set_index('substance')
df.flock_classified = df.flock_classified.fillna('not_classified')
df.man_classified = df.man_classified.fillna('not_classified')
df[['2d6sub_flock', '2c19sub_flock', '2c9sub_flock']] = df[['2d6sub_flock', '2c19sub_flock', '2c9sub_flock']].fillna(0)
df[['2d6sub_man', '2c19sub_man', '2c9sub_man']] = df[['2d6sub_man', '2c19sub_man', '2c9sub_man']].fillna(0)

check1 = df.man_classified == 'classified'
check2 = df.man_classified == 'not_classified'
check3 = df.flock_classified == 'classified'
check4 = df.flock_classified == 'not_classified'
data = df[check4].sort_values('counts', ascending=False)[:20]



'''


MAN: OLD TABLE AND NEW TABLE COMPARISON


'''

man1 = man.set_index('substances')
man1 = man1.loc[man1['man_classified']=='classified']
man1
old_man = pd.read_excel('v10.24/classified_old_man.xlsx', index_col=0)
old_man1 = old_man.iloc[:, [4,5,6,8]]
df = old_man1.join(man1, how='outer')
df.to_excel('temp1.xlsx')



'''


TOTAL AMOUNT OF THE FLOCKHART ANALYSIS


'''

substrate = pd.read_excel("FlockhardtTableExtended.xlsx", sheet_name="Substrat")
inhibitor = pd.read_excel("FlockhardtTableExtended.xlsx", sheet_name="Inhibitor")
inducer = pd.read_excel("FlockhardtTableExtended.xlsx", sheet_name="Inducer")
selection = ['2D6', '2C19', '2C9']

substrate = substrate[selection]
inhibitor = inhibitor[selection]
inducer = inducer[selection]

substrate = substrate.values.flatten().tolist()
inhibitor = inhibitor.values.flatten().tolist()
inducer = inducer.values.flatten().tolist()

li = substrate + inhibitor + inducer

li = list(set(li))
df = pd.Series(li)
df = df.dropna()
li = df.tolist()
li

man[man['substances'].isin(li)]


'''


VENN CHART


'''
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from venny4py.venny4py import *

# PREPARE FLOCKHARDT
# mark properly if a substance is classified or not
substrate = pd.read_excel("FlockhardtTableExtended.xlsx", sheet_name="Substrat")
inhibitor = pd.read_excel("FlockhardtTableExtended.xlsx", sheet_name="Inhibitor")
inducer = pd.read_excel("FlockhardtTableExtended.xlsx", sheet_name="Inducer")
substrate = substrate.values.flatten().tolist()
inhibitor = inhibitor.values.flatten().tolist()
inducer = inducer.values.flatten().tolist()

li = substrate + inhibitor + inducer

li = list(set(li))
df = pd.Series(li)
df = df.dropna()
li = df.tolist()
li
flock_copy = flock.set_index('drugs')

drc_reduced = drc.drop(columns='interactions')
drc_reduced.set_index('substance', inplace=True)
man_reduced = man.drop(columns='Z_Wirkstoff_Gesamt')
man_reduced.set_index('substances', inplace=True)

df = drc_reduced.join(man_reduced, how='outer')
df = df.rename(columns={
    'index': 'substance',
    '2d6sub':'2d6sub_db',
    '2c19sub':'2c19sub_db',
    '2c9sub': '2c9sub_db',
    'CYP2D6': '2d6sub_man',
    'CYP2C19': '2c19sub_man',
    'CYP2C9': '2c9sub_man'})
df = df.join(fda, how='left', lsuffix='_man', rsuffix='_fda')



df = df.join(flock_copy, how='left', lsuffix='_fda', rsuffix='_flock')
df.loc[df.index.isin(li), "flock_classified"] = 'classified'
df.flock_classified = df.flock_classified.fillna('not_classfied')



check = df.db_classified == "classified"
set1 = set(df[check].index)
check = df.man_classified == "classified"
set2 = set(df[check].index)
check = df.fda_classified == "classified"
set3 = set(df[check].index)
check = df.flock_classified == "classified"
set4 = set(df[check].index)

# Create a Venn diagram
#venn_diagram = venn2([set1, set2], set_labels=('Drugbank', 'Manual'))
sets = {
    'Drugbank': set1,
    'Manual': set2, 
    "FDA": set3, 
    "Flockhart":set4
}
    
venny4py(sets=sets)

# Display the plot
plt.savefig('venn_plot.png', dpi=300)
plt.show()
