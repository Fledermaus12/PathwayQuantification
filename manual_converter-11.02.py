import pandas as pd 
import numpy as np

manual = pd.read_excel('CYP_Substrate_haendisch.xlsx')
manual = manual[['Z_Wirkstoff_Gesamt','CYP2D6.1', 'CYP2C19.1', 'CYP2C9.1']]
manual = manual.rename(columns={
    'CYP2D6.1':'CYP2D6', 
    'CYP2C19.1': 'CYP2C19', 
    'CYP2C9.1': 'CYP2C9'
})
manual['man_classified'] = manual.apply(lambda row: 'classified' if row.isnull().sum() == 0 else 'not_classified', axis=1)
manual = manual.fillna(0)

# COUNT TOTAL ENTRIES
ec = manual[~manual['Z_Wirkstoff_Gesamt'].str.contains(r"/", regex=False)].copy()
eccs = manual[manual['Z_Wirkstoff_Gesamt'].str.contains(r"/", regex=False)]
df = eccs.assign(Z_Wirkstoff_Gesamt=eccs['Z_Wirkstoff_Gesamt'].str.split('/'))
df = df.explode('Z_Wirkstoff_Gesamt')
df = pd.concat([ec, df])
total_entries = df.copy()
# 3922 + 20257 = 24179

# Hier finden wir nur noch Substanzen anstatt Einträge vor
manual_uni = manual.drop_duplicates(keep='first', inplace=False, ignore_index=False)



'''
Splitting der Kombinationspräparate ---
1) Splitting nach "/" in eine Liste
2) Für jede Substanz eine Zeile
3) sortiere die Substanzen aus, die als Einzelsubstanz bereits KLASSIFIZIERT sind.
Arbeite mit den Substanzen weiter, die als Einzelsubstanz noch nicht klassifiziert waren
4) zeige nur Substanzen an, die als Kombi klassifiziert wurden. 

5) sortiere die Substanzen nach oben, die als Kombi Interaktion verursachen ->
Beachte die Substanzen, die hier als 1 markiert werden. 
Frage: welche von diesen als 1 markierten Substanzen ist wirklich eine 1? 
denn die 1 gehört zum Kombinationspräparat aber nicht zwangsläufig zu den einzelnen Bestandteilen davon

6) die Substanzen, die nicht als Einzelsubstanz und nicht als Kombi klassifiziert sind
werden der Haupttabelle angefügt

7) Substanzen, die nicht als Einzelsubstanz aber als Kombi klassifiziert 
und ZEROS sind, 
werden der Haupttabelle angefügt

8) Substanzen, die nicht als Einzelsubstanz aber als Kombi klassifiziert
und ONES sind, 
werden korrigiert,
werden der Haupttabelle angefügt

9) Duplikate entfernen
'''
# manual_uni_clean enthält nach dieser Zeile nur noch Einzelsubstanzen und keine Kombinationspräparate
manual_uni_clean = manual_uni[~manual_uni['Z_Wirkstoff_Gesamt'].str.contains(r"/", regex=False)].copy()
duplicates = manual_uni_clean[manual_uni_clean['Z_Wirkstoff_Gesamt'].duplicated(keep=False)]
duplicates # must be empty


check_slash = manual_uni[manual_uni['Z_Wirkstoff_Gesamt'].str.contains(r"/", regex=False)]
check_slash_one = check_slash[check_slash[['CYP2D6', 'CYP2C19', 'CYP2C9']].sum(axis=1) != 0] # sort out the interesting ones
check_slash_zero = check_slash[check_slash[['CYP2D6', 'CYP2C19', 'CYP2C9']].sum(axis=1) == 0] # sort out the interesting ones
check_slash_classified = check_slash[check_slash['man_classified'] == 'classified']
 
# 1)
df = check_slash.assign(Z_Wirkstoff_Gesamt=check_slash['Z_Wirkstoff_Gesamt'].str.split('/'))
df.sort_values(by=['CYP2D6', 'CYP2C19', 'CYP2C9'], ascending=False)[:15]
# 2)
df = df.explode('Z_Wirkstoff_Gesamt')
# in klassifiert und unklassifiziert unterteilen und dann die Duplikate in jedem Set löschen,
# als nächstes die zwei sets untereinander zusammenführen 
classified_rows = df[df['man_classified'] == 'classified']
classified_unique = classified_rows.drop_duplicates(subset='Z_Wirkstoff_Gesamt')
unclassified_rows = df[df['man_classified'] == 'not_classified']
unclassified_unique = unclassified_rows.drop_duplicates(subset='Z_Wirkstoff_Gesamt')
df = pd.concat([classified_unique, unclassified_unique]).drop_duplicates(subset='Z_Wirkstoff_Gesamt', keep='first')


# 3)
classified = list(manual_uni_clean[manual_uni_clean['man_classified'] == 'classified']['Z_Wirkstoff_Gesamt'])
def check_if_in_muc(substanz):
    if substanz in classified: return True
    else: return False
df['singleclass'] = df['Z_Wirkstoff_Gesamt'].apply(check_if_in_muc)

# 4)
df_new = df[df.singleclass == False]

# 5)
df_new_class = df_new[df_new['man_classified'] == 'classified']
df_new_class.sort_values(by=['CYP2D6', 'CYP2C19', 'CYP2C9'], ascending=False)

columns_to_check = ['CYP2D6', 'CYP2C19', 'CYP2C9']
ones = df_new_class[df_new_class[columns_to_check].eq(1.0).any(axis=1)].copy()
zeros = df_new_class[~df_new_class[columns_to_check].eq(1.0).any(axis=1)].copy()

# 6) completely unclassified
df_new_unclass = df_new[df_new['man_classified'] != 'classified'].copy()
df_new_unclass.drop_duplicates()

manual_uni_clean = pd.concat([manual_uni_clean, df_new_unclass])

manual_uni_clean[manual_uni_clean['Z_Wirkstoff_Gesamt'].duplicated(keep=False)] 
manual_uni_clean.drop_duplicates(subset='Z_Wirkstoff_Gesamt', inplace=True, keep="last")


# 7) als Einzelsubstanz -> unklassifiziert, Kombi -> klassifiziert, ZEROS
manual_uni_clean = pd.concat([manual_uni_clean, zeros])

manual_uni_clean[manual_uni_clean['Z_Wirkstoff_Gesamt'].duplicated(keep=False)] 
manual_uni_clean.drop_duplicates(subset='Z_Wirkstoff_Gesamt', inplace=True, keep="last")

# 8) als Einzelsubstanz -> unklassifiziert, Kombi -> klassifiziert, ONES
ones
ones.set_index('Z_Wirkstoff_Gesamt', inplace=True)
ones.loc['Naloxon', 'CYP2C19'] = 0.0
ones.loc['Ezetimib', 'CYP2D6'] = 0.0
ones.loc['Olodaterol', 'CYP2D6'] = 0.0
ones.loc['Dutasterid', 'CYP2D6'] = 0.0
ones.loc['Pseudoephedrin', 'CYP2C9'] = 0.0
ones.loc['Clarithromycin', 'CYP2C19'] = 0.0
ones.loc['Amoxicillin', 'CYP2C19'] = 0.0
ones = ones.reset_index()
ones
manual_uni_clean = pd.concat([manual_uni_clean, ones])

manual_uni_clean[manual_uni_clean['Z_Wirkstoff_Gesamt'].duplicated(keep=False)]
manual_uni_clean.drop_duplicates(subset='Z_Wirkstoff_Gesamt', inplace=True, keep="last")


# 9)
manual_uni_clean[manual_uni_clean['Z_Wirkstoff_Gesamt'].duplicated(keep=False)] # empty

manual_uni_clean.drop(columns="singleclass", inplace=True)
manual_uni_clean.reset_index(drop=True, inplace=True)



'''

TRANSLATION

'''
translation = pd.read_excel("../data/CID_DB_directory.xlsx")
translation = translation[['compound_eng','compound_ger']]
translation = translation.set_index('compound_ger')
translation = translation.dropna()
translation = translation[~translation.apply(lambda row: row.str.contains(r"\\")).any(axis=1)]


def checktranslation(med):
    if med in translation.index:
        engl = translation.loc[med, 'compound_eng']
        if type(engl) != str:
            print(med, engl)
            engl = engl[0]
            
        return engl
    else:
        return 'error1'
    
substances_engl = []
for entry in list(manual_uni_clean.Z_Wirkstoff_Gesamt):
    answer = checktranslation(entry)
    substances_engl.append(answer)
print(substances_engl)
print(len(substances_engl))
manual_uni_clean['substances'] = substances_engl
df = manual_uni_clean.copy()

df[df['substances'] == 'error1']
df['man_classified'].value_counts()



errors = df[df['substances'] == 'error1'].copy()
errors = errors['Z_Wirkstoff_Gesamt'].to_list()
total_entries[~total_entries['Z_Wirkstoff_Gesamt'].isin(errors)]

#

df.to_excel('v11.11/manual_table_translated1.xlsx')
df.to_json('v11.11/manual_table_translated1.json')
