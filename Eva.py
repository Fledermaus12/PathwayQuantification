##
import pandas as pd

data = pd.read_json(r'base_output/drc_p2_pv2.json', orient='index', convert_axes=False, dtype=False)


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
    search_result = search(i, list(data['AllInteraction']))
    data[i] = search_result

##
data.to_excel(r'evaluater_output/drc_p3eva.xlsx')
data.to_json(r'evaluater_output/drc_p3eva.json',orient='index', force_ascii='False')