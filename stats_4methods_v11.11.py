# Imports

import pandas as pd
import numpy as np
from sklearn.metrics import cohen_kappa_score
import os

cwd = os.getcwd()
print(cwd)

# translate
translation = pd.read_excel("../data/CID_DB_directory.xlsx")
translation = translation[['compound_eng','compound_ger']]
translation = translation.set_index('compound_eng')
translation = translation.dropna()
translation = translation[~translation.apply(lambda row: row.str.contains(r"\\")).any(axis=1)]
translation_uni = translation.reset_index()
translation_uni.drop_duplicates(subset='compound_eng', inplace=True)
translation_uni.set_index('compound_eng', inplace=True)

og_table = pd.read_excel("../data/tbl_Medikation.xlsx")
man_og_start = pd.read_excel("../data/CYP_Substrate_haendisch.xlsx")
man_og = man_og_start[['Z_Wirkstoff_Gesamt','CYP2D6.1', 'CYP2C19.1', 'CYP2C9.1']]
man_og = man_og.rename(columns={
    'CYP2D6.1':'CYP2D6', 
    'CYP2C19.1': 'CYP2C19', 
    'CYP2C9.1': 'CYP2C9'
})
man_og = man_og.dropna()

### Function define
def create_proportion_table(name1, name2, se_absolute, se_percent): 
    '''
        se_absolute and se_percent follow the template: 
        2D6_table1, 2C19_table1, 2C9_table1, 2D6_table2, 2C19_table2, 2C9_table2

        names are the names of the table (table1, table2)
    '''
    se = {
        f'{name1} (count)': list(se_absolute[0:3]),
        f'{name2} (count)': list(se_absolute[3:]),
        f'{name1} (percent)': list(se_percent[0:3]),
        f'{name2} (percent)': list(se_percent[3:])
    }
    se = pd.DataFrame(se)
    se['Enzyme'] = ['2D6', '2C19', '2C9']
    se.set_index('Enzyme', inplace=True)
    print(se)
    
def find_diff(df_diff, col1, col2):
    opt = df_diff[['substance',col1, col2]]
    optdiff = opt.loc[opt.loc[:,col1] != opt.loc[:,col2]]
    optgb = optdiff.groupby(by=["substance"]).sum()
    return optgb




def count_occurence2(engmed):
    #print(engmed)
    if engmed in translation_uni.index:
        demed = translation_uni.loc[engmed, 'compound_ger']
    else: 
        print(engmed,'fail to translate')
        demed = 'error1'
    if og_table['Z_Wirkstoff_Gesamt'].str.contains(demed).any() == True:
        count = len(og_table[og_table['Z_Wirkstoff_Gesamt'].str.contains(demed)])
    else: 
        print(engmed, demed, 'fail to find in table')
        count = 0
    return count

def check_listed(engmed, not_in_list):
    if engmed in not_in_list:
        exist = 'not_classified'
    else: 
        exist = 'classified'
    return exist

def translate_and_check_listed(engmed):
    if engmed in translation_uni.index:
        demed = translation_uni.loc[engmed, 'compound_ger']
    else: 
        print(engmed,'fail to translate')
        demed = 'error1'
    if demed == 'error1':
        exist = 'error1'
    elif man_og['Z_Wirkstoff_Gesamt'].str.contains(demed).any() == True:
        exist = 'classified'
    else: 
        exist = 'not_classified'
    return exist



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









# Drugbank vs Manual

drc_reduced = drc.drop(columns='interactions')
drc_reduced.set_index('substance', inplace=True)
man_reduced = man.drop(columns='Z_Wirkstoff_Gesamt')
man_reduced.set_index('substances', inplace=True)
automatic_list = list(drc_reduced.index)
manual_list = list(man_reduced.index)

set1 = set(automatic_list)
set2 = set(manual_list)
not_in_list1 = list(set2 - set1)
not_in_list2 = list(set1 - set2)

print("Elements in list2 but not in list1:", len(not_in_list1), not_in_list1)
print("Elements in list1 but not in list2:", len(not_in_list2), not_in_list2)

# Beispiel von später
df = drc_reduced.join(man_reduced, how='outer')
df = df.rename(columns={
    'index': 'substance',
    '2d6sub':'2d6sub_db',
    '2c19sub':'2c19sub_db',
    '2c9sub': '2c9sub_db',
    'CYP2D6': '2d6sub_man',
    'CYP2C19': '2c19sub_man',
    'CYP2C9': '2c9sub_man'})

df_export = df.copy()
counts = [count_occurence2(e) for i, e in enumerate(list(df_export.index))]
df_export['Count'] = counts
df_export
df_export.to_excel('v11.11/man_db_compared.xlsx')

print(df.man_classified.value_counts())
print(df.db_classified.value_counts())
check1 = df['db_classified']=='classified'
check2 = df['man_classified']=='classified'
df_portion_classified = df[check1 & check2].copy()
df_portion_classified.drop(['db_classified', 'man_classified'], axis=1, inplace=True)
print("Length of the classified only table: ", len(df_portion_classified))

print('Schnittmenge:')
res_absolute = df_portion_classified.iloc[:, :].sum()
res_perc = round(res_absolute / len(df_portion_classified), 3)
create_proportion_table('Drugbank', 'Manual', res_absolute, res_perc)

print()
print('Gesamtes Set:')
df_matrix = df.drop(['db_classified', 'man_classified'], axis=1)
res1 = df_matrix.iloc[:, :].sum()
res2 = round(df_matrix.iloc[:, :].sum() / len(df_matrix), 3)
create_proportion_table('Drugbank', 'Manual', res1, res2)

# Kappa Score
print("\nCohen's Kappa Score - by substance : 2d6, 2c19, 2c9")
print(round(cohen_kappa_score(df_portion_classified['2d6sub_db'], df_portion_classified['2d6sub_man']), 3))
print(round(cohen_kappa_score(df_portion_classified['2c19sub_db'], df_portion_classified['2c19sub_man']), 3))
print(round(cohen_kappa_score(df_portion_classified['2c9sub_db'], df_portion_classified['2c9sub_man']), 3))


# Show differences - Schnittmenge und ohne Schnittmenge in Einem
df_diff = df_matrix.reset_index()
df_diff.rename(columns={'index':'substance'}, inplace=True)
enzyme = '2d6'
table = find_diff(df_diff, f'{enzyme}sub_db', f'{enzyme}sub_man')
counts = [count_occurence2(e) for i, e in enumerate(list(table.index))]
table['Count'] = counts

classified = df[['db_classified','man_classified']]
table = table.join(classified, how='left')

#table['manually_classified'] = [translate_and_check_listed(e) for i, e in enumerate(list(table.index))]
print('\nLength of differences:', len(table))
table.to_excel(f'v11.11/diff_mandb_{enzyme}.xlsx')
table




###############################################




# FDA
man_reduced = man.drop(columns='Z_Wirkstoff_Gesamt')
man_reduced = man_reduced.rename(columns={
    'substances': 'substance',
    'CYP2D6': '2d6sub',
    'CYP2C19': '2c19sub',
    'CYP2C9': '2c9sub'})
man_reduced = man_reduced.set_index('substance')

man_list = list(man_reduced.index)
fda_list = list(fda.index)

set1 = set(man_list)
set2 = set(fda_list)
not_in_list1 = list(set2 - set1)
not_in_list2 = list(set1 - set2)

print("Elements in list2 but not in list1:", len(not_in_list1), not_in_list1)
print("Elements in list1 but not in list2:", len(not_in_list2), not_in_list2)


df = man_reduced.join(fda, how='left', lsuffix='_man', rsuffix='_fda')
df = df.reset_index()
df = df.rename(columns={'index': 'substance'})
df = df.drop_duplicates(subset='substance')
df = df.set_index('substance')
df.fda_classified = df.fda_classified.fillna('not_classified')
df.man_classified = df.man_classified.fillna('not_classified')
df[['2d6sub_fda', '2c19sub_fda', '2c9sub_fda']] = df[['2d6sub_fda', '2c19sub_fda', '2c9sub_fda']].fillna(0)
df[['2d6sub_man', '2c19sub_man', '2c9sub_man']] = df[['2d6sub_man', '2c19sub_man', '2c9sub_man']].fillna(0)
print("The fusion result in this many substances: ", len(df))

check1 = df.man_classified == 'classified'
check2 = df.man_classified == 'not_classified'
check3 = df.fda_classified == 'classified'
check4 = df.fda_classified == 'not_classified'

df[check1]
df_export = df.copy()
counts = [count_occurence2(e) for i, e in enumerate(list(df_export.index))]
df_export['Count'] = counts
df_export
df_export.to_excel('v11.11/man_fda_compared.xlsx')

# Schnittmenge and count them
df_class = df.loc[check1 & check3,:].copy()
df_class

print('This many substances are in both tables:', len(df_class))

# Kappa Score
print("\nCohen's Kappa Score - by substance: 2d6, 2c19, 2c9")
print(round(cohen_kappa_score(df_class['2d6sub_man'], df_class['2d6sub_fda']), 3))
print(round(cohen_kappa_score(df_class['2c19sub_man'], df_class['2c19sub_fda']), 3))
print(round(cohen_kappa_score(df_class['2c9sub_man'], df_class['2c9sub_fda']), 3))

# Show proportion table

df_class.drop(['fda_classified', 'man_classified'], axis=1, inplace=True)
res1 = df_class.sum()
res2 = round(df_class.sum() / len(df_class), 3)
create_proportion_table('Manuell', 'FDA', res1, res2)
"""
# Show differences and count them 
df_diff = df_class.reset_index()
enzyme = '2c19'
table = find_diff(df_diff, f'{enzyme}sub_man', f'{enzyme}sub_fda')
counts = [count_occurence2(e) for i, e in enumerate(list(table.index))]
table['Count'] = counts
print('\nLength of differences:', len(table))
table.to_excel(f'v11.11/diff_manfda_{enzyme}.xlsx')
print(table)
"""

# whole set

df = man_reduced.join(fda, how='left', lsuffix='_man', rsuffix='_fda')
df = df.reset_index()
df = df.rename(columns={'index': 'substance'})
df = df.drop_duplicates(subset='substance')
df = df.set_index('substance')
df.fda_classified = df.fda_classified.fillna('not_classfied')
df.man_classified = df.man_classified.fillna('not_classified')
df[['2d6sub_fda', '2c19sub_fda', '2c9sub_fda']] = df[['2d6sub_fda', '2c19sub_fda', '2c9sub_fda']].fillna(0)
df[['2d6sub_man', '2c19sub_man', '2c9sub_man']] = df[['2d6sub_man', '2c19sub_man', '2c9sub_man']].fillna(0)

print("The fusion result in this many substances: ", len(df))
#print(df)
#df.to_excel('man_fda_onlyadred.xlsx')
classified = df[['fda_classified','man_classified']]
df.drop(['fda_classified', 'man_classified'], axis=1, inplace=True)


# print proportion table
res1 = df.sum()
res2 = round(df.sum() / len(df), 3)
create_proportion_table('Manual', 'FDA', res1, res2)

# print differences
df_diff = df.reset_index()
df_diff = df_diff.rename(columns={'index': 'substance'})

#print(df_diff)
#df_diff = df_diff.fillna(0)
# select enzyme pair
enzyme = '2c19'
table = find_diff(df_diff, f'{enzyme}sub_man', f'{enzyme}sub_fda')
counts = [count_occurence2(e) for i, e in enumerate(list(table.index))]
table['Count'] = counts


table = table.join(classified, how='left')
print('\nLength of differences:', len(table))
table.to_excel(f'diff_manfda_clean_{enzyme}.xlsx')
table


###########################################




# Flockhardt

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





man_reduced = man.drop(columns='Z_Wirkstoff_Gesamt')
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
df
df.loc[df.index.isin(li), "flock_classified"] = 'classified'
df.flock_classified = df.flock_classified.fillna('not_classfied')
df.man_classified = df.man_classified.fillna('not_classified')

df[['2d6sub_flock', '2c19sub_flock', '2c9sub_flock']] = df[['2d6sub_flock', '2c19sub_flock', '2c9sub_flock']].fillna(0)
df[['2d6sub_man', '2c19sub_man', '2c9sub_man']] = df[['2d6sub_man', '2c19sub_man', '2c9sub_man']].fillna(0)
print("The fusion result in this many substances: ", len(df))

check1 = df.man_classified == 'classified'
check2 = df.man_classified == 'not_classified'
check3 = df.flock_classified == 'classified'
check4 = df.flock_classified == 'not_classified'

df_export = df.copy()
counts = [count_occurence2(e) for i, e in enumerate(list(df_export.index))]
df_export['Count'] = counts
df_export
df_export.to_excel('v11.11/man_flock_compared.xlsx')
df[check3]

df_class = df.loc[check1 & check3,:].copy()
print('This many substances are in both tables:', len(df_class))


# Kappa Score
print("\nCohen's Kappa Score - by substances: 2d6, 2c19, 2c9")
print(round(cohen_kappa_score(df_class['2d6sub_man'], df_class['2d6sub_flock']), 3))
print(round(cohen_kappa_score(df_class['2c19sub_man'], df_class['2c19sub_flock']), 3))
print(round(cohen_kappa_score(df_class['2c9sub_man'], df_class['2c9sub_flock']), 3))

# print proportion table
df_matrix = df_class.drop(['flock_classified', 'man_classified'], axis=1)
res1 = df_matrix.iloc[:, :].sum()
res2 = round(df_matrix.iloc[:, :].sum() / len(df_matrix), 3)
create_proportion_table('Manual', 'Flock', res1, res2)



# find differences
df_diff = df_class.reset_index()
enzyme = '2c9'
table = find_diff(df_diff, f'{enzyme}sub_man', f'{enzyme}sub_flock')
counts = [count_occurence2(e) for i, e in enumerate(list(table.index))]
table['Count'] = counts
print(f'\nLength of differences for {enzyme}:', len(table))
#table.to_excel(f'diff_manflock_{enzyme}.xlsx')
print(table)


# Gesamtes Set 

df = man_reduced.join(flock_copy, how='left', lsuffix='_man', rsuffix='_flock')
df = df.reset_index()
df = df.rename(columns={'index': 'substance'})
df = df.drop_duplicates(subset='substance')
df = df.set_index('substance')
df.loc[df.index.isin(li), "flock_classified"] = 'classified'
df.flock_classified = df.flock_classified.fillna('not_classfied')
df.man_classified = df.man_classified.fillna('not_classified')
classified = df[['flock_classified','man_classified']]

df[['2d6sub_flock', '2c19sub_flock', '2c9sub_flock']] = df[['2d6sub_flock', '2c19sub_flock', '2c9sub_flock']].fillna(0)
df[['2d6sub_man', '2c19sub_man', '2c9sub_man']] = df[['2d6sub_man', '2c19sub_man', '2c9sub_man']].fillna(0)
print("The fusion result in this many substances: ", len(df))

# print proportion table
df_matrix = df.drop(['flock_classified', 'man_classified'], axis=1)
res1 = df_matrix.iloc[:, :].sum()
res2 = round(df_matrix.iloc[:, :].sum() / len(df_matrix), 3)
create_proportion_table('Manual', 'Flock', res1, res2)

# find differences
df = df.reset_index()
enzyme = '2c9'
table = find_diff(df, f'{enzyme}sub_man', f'{enzyme}sub_flock')
counts = [count_occurence2(e) for i, e in enumerate(list(table.index))]
table['Count'] = counts


table = table.join(classified, how='left')
print('\nLength of differences:', len(table))
table.to_excel(f'v11.11/diff_manflock_{enzyme}.xlsx')
table