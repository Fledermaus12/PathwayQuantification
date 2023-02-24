import pandas as pd
import xml.etree.ElementTree as ET

def runit(var, root): # submodule of findit
    '''r = root of element-tree
    var = the drug that needs to be found'''
    pfadstart = r"{http://www.drugbank.ca}drug[{http://www.drugbank.ca}name='"
    pfadend = r"'" + r']'
    pfad = pfadstart + var + pfadend
    drugroot = root.find(pfad)
    return drugroot


def find_drug_root(atc, root):
    '''
    :param atc: ATC code that you want to look up
    :param root: root
    :return: a list of nodes
    '''
    varbegin = r"*[@code='"
    varend = r"']"
    var = varbegin + atc + varend
    # create path searching term
    print('searching for ', atc)
    drugroot, atcnode = list(), list()
    for sn in root:
        ssn = sn.find("{http://www.drugbank.ca}atc-codes")
        sssn = ssn.find(var)
        if sssn != None:
            atcnode.append(sssn)
            drugroot.append(sn) # save the drug root of the success
    if len(drugroot) != 0:
        err = 0
        for idx_element in range(0,len(drugroot)):
            drugname = drugroot[idx_element].find("{http://www.drugbank.ca}name").text
            atc_check = atcnode[idx_element].get('code')
            print(idx_element, atc_check, drugname)
    else: err = 1
    return atcnode, drugroot, err


def find_enzymes(node, root):
    identifier = list()
    if node == None:
        print("not found")
        return identifier
    # if node is empty do not continue

    drugname = node.find("{http://www.drugbank.ca}name").text
    print(drugname)
    node_enzymes = node.find('{http://www.drugbank.ca}enzymes')
    if node_enzymes != None:
        # only loop over enzymes if any enzymes are listed
        # if no enzymes are listed then just return an empty identifier

        # enzymes are listed, so we can loop over!
        for enzyme in node_enzymes:
            enzyme_name = enzyme.find("{http://www.drugbank.ca}name").text
            print(drugname, enzyme_name)

            if enzyme_name == "Cytochrome P450 2D6":
                node_enz_sn = node_enzymes.find("{http://www.drugbank.ca}enzyme[{http://www.drugbank.ca}name='Cytochrome P450 2D6']")
                node_enz_sn_action = node_enz_sn.find("{http://www.drugbank.ca}actions")
                for child in node_enz_sn_action:
                    print(drugname, '--> Cytochrome P450 2D6 --->', child.text)
                    if child.text == 'inducer': identifier.append('2d6ind')
                    if child.text == 'substrate': identifier.append('2d6sub')
                    if child.text == 'inhibitor': identifier.append('2d6inh')
            if enzyme_name == 'Cytochrome P450 2C19':
                node_enz_sn = node_enzymes.find("{http://www.drugbank.ca}enzyme[{http://www.drugbank.ca}name='Cytochrome P450 2C19']")
                node_enz_sn_action = node_enz_sn.find("{http://www.drugbank.ca}actions")
                for child in node_enz_sn_action:
                    print(drugname, '--> Cytochrome P450 2C19 --->', child.text)
                    if child.text == 'inducer': identifier.append('2c19ind')
                    if child.text == 'substrate': identifier.append('2c19sub')
                    if child.text == 'inhibitor': identifier.append('2c19inh')
            if enzyme_name == 'Cytochrome P450 2C9':
                node_enz_sn = node_enzymes.find("{http://www.drugbank.ca}enzyme[{http://www.drugbank.ca}name='Cytochrome P450 2C9']")
                node_enz_sn_action = node_enz_sn.find("{http://www.drugbank.ca}actions")
                for child in node_enz_sn_action:
                    print(drugname, '--> Cytochrome P450 2C9 --->', child.text)
                    if child.text == 'inducer': identifier.append('2c9ind')
                    if child.text == 'substrate': identifier.append('2c9sub')
                    if child.text == 'inhibitor': identifier.append('2c9inh')
            if enzyme_name == 'Cytochrome P450 3A4':
                node_enz_sn = node_enzymes.find("{http://www.drugbank.ca}enzyme[{http://www.drugbank.ca}name='Cytochrome P450 3A4']")
                node_enz_sn_action = node_enz_sn.find("{http://www.drugbank.ca}actions")
                for child in node_enz_sn_action:
                    print(drugname, '--> Cytochrome P450 3A4 --->', child.text)
                    if child.text == 'inducer': identifier.append('3a4ind')
                    if child.text == 'substrate': identifier.append('3a4sub')
                    if child.text == 'inhibitor': identifier.append('3a4inh')
            if enzyme_name == 'Cytochrome P450 1A2':
                node_enz_sn = node_enzymes.find("{http://www.drugbank.ca}enzyme[{http://www.drugbank.ca}name='Cytochrome P450 1A2']")
                node_enz_sn_action = node_enz_sn.find("{http://www.drugbank.ca}actions")
                for child in node_enz_sn_action:
                    print(drugname, '--> Cytochrome P450 1A2 --->', child.text)
                    if child.text == 'inducer': identifier.append('1a2ind')
                    if child.text == 'substrate': identifier.append('1a2sub')
                    if child.text == 'inhibitor': identifier.append('1a2inh')
    return identifier

def make_enzyme_list(atc_list, root):
    fehler = list()
    enzli = list()
    for elem in atc_list:
        output = list()
        if pd.isna(elem) != True:
            __, node_drug, single_error = find_drug_root(root, elem)
            output = find_enzymes(node_drug, root)
            if single_error == 1:
                fehler.append(elem) # check, if a medication has not been found
                output.append('err')
        enzli.append(output)
    return enzli, fehler