import mod
import glob
import os
import re
import regex
import itertools
import requests
import networkx as nx
import copy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import openbabel
from mod import *

config.rule.ignoreConstraintsDuringInversion = True

# load compound mol format form KEGG Database
def get_compound_mol(compound_id):
    url = f'http://rest.kegg.jp/get/compound:{compound_id}/mol'
    try:
        response = requests.get(url)
        response.raise_for_status()  # This throws an exception if the status code is not 200 (OK)
        return response.text
    except requests.exceptions.RequestException as e:
        #print(f"Fehler bei der Anfrage: {e}")
        return None

def mol_to_smiles(molfile_content):
    if isinstance(molfile_content, bytes):
        molfile_content = molfile_content.decode('utf-8')
    elif not isinstance(molfile_content, str):
        return False

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "smi")
    mol = openbabel.OBMol()
    success = obConversion.ReadString(mol, molfile_content)
    
    if not success:
        return False
    
    smiles = obConversion.WriteString(mol).strip()
    if smiles != '[H+]':
        smiles_without_charges = re.sub(r'\[\w+[\+\-]\d*\]', lambda x: x.group(0).replace("+", "").replace("-", ""), smiles)
    else:
        return smiles
    return smiles_without_charges


# create the compound form KEGG mol data 
def mol_to_graph(mol_str):
    mol = Chem.MolFromMolBlock(mol_str)
    Chem.rdmolops.AssignStereochemistry(mol, force=True, cleanIt=True)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL - Chem.SanitizeFlags.SANITIZE_PROPERTIES)

    G = nx.Graph()
    for atom in mol.GetAtoms():
        label = atom.GetSymbol()
        G.add_node(atom.GetIdx(), label=label)
    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        bond_type = bond.GetBondTypeAsDouble()
        G.add_edge(u, v, label=bond_type)
    return G  
   
def count_protons(text):
    # Define patterns for the left block and right block
    left_block_pattern = r"left\s*\[(.*?)\]\s*\]"
    right_block_pattern = r"right\s*\[(.*?)\]\s*\]"

    # Extract the left block and right block
    left_block = re.search(left_block_pattern, text, re.DOTALL)
    right_block = re.search(right_block_pattern, text, re.DOTALL)

    # Count the number of 'H+' in the blocks
    left_h_plus_count = len(re.findall(r'label\s*\"H\+\"', left_block.group(1))) if left_block else 0
    right_h_plus_count = len(re.findall(r'label\s*\"H\+\"', right_block.group(1))) if right_block else 0

    return left_h_plus_count, right_h_plus_count
    

def getReactionSmiles(dg,correct_edge):
        origSmiles = {}
        # change the class labels for every atom
        for v in dg.vertices:
                s = v.graph.smilesWithIds
                #s = s.replace(":", ":o")
                s = regex.sub(":([0-9]+)\]",":o\\1]",s)
                origSmiles[v.graph] = s
        res = {}
        
        maps = DGVertexMapper(correct_edge,leftLimit=1,rightLimit=1)

        # replace the original local IDs with global IDs
        eductSmiles = [origSmiles[g] for g in maps.left]
        
        for v in maps.left.vertices:
            s = eductSmiles[v.graphIndex]
            s = s.replace(":o%d]" % v.vertex.id, ":%d]" % v.id)
            eductSmiles[v.graphIndex] = s

        strs = set()
        for r in maps:
            #print('No of Maps:',len(r.map))
            m = r.map
            productSmiles = [origSmiles[g] for g in maps.right]
            for ev in maps.left.vertices:
                pv = m[ev]
                if not pv:
                    continue
                s = productSmiles[pv.graphIndex]
                s = s.replace(f":o{pv.vertex.id}]", f":{ev.id}]")            
                productSmiles[pv.graphIndex] = s
            
            count = maps.left.numVertices
            
            for pv in maps.right.vertices:
                ev = m.inverse(pv) 
                if ev: 
                    continue
                s = productSmiles[pv.graphIndex]
                s = s.replace(f":o{pv.vertex.id}]", f":{count}]")
                count += 1
                productSmiles[pv.graphIndex] = s    
                    
            left = ".".join(eductSmiles)
            right = ".".join(productSmiles)
            s = f"{left}>>{right}"
            assert ":o" not in s
            strs.add(s)

        return list(sorted(strs))

def generate_combinations(input_list):
    # Liste, um die Kombinationen zu speichern
    all_combinations = []
    
    # Erzeuge alle Kombinationen für jede mögliche Länge
    for r in range(1, len(input_list)):
        combinations_r = list(itertools.combinations(input_list, r))
        # Wandel jede Kombination von Tupeln in eine Liste um
        all_combinations.extend([list(comb) for comb in combinations_r])
    
    return all_combinations

def create_Educts_and_Products(reaction_data,path,rule,rule_var,universe):

    educ = [] # List of educts known by the database for finding correct hyperedge
    prod = [] # List of products known by the database for finding correct hyperedge
    h_mol = None # marker if proton molecules were at the Kegg reaction
    subspace = [] # List of possibile substrates for creating Network
    
    # load all compounants of the reaction    
    educts = copy.deepcopy(reaction_data[rule][0])
    products = copy.deepcopy(reaction_data[rule][1])
    
    # Define all educts 
    for e in educts.keys():
        if e == 'C00080':
            h_mol = universe[e]
        while educts[e] > 0:            
            if e == False:
                break
            elif isinstance(universe[e], Graph) and universe[e].id:
                educ.append(universe[e])
                educts[e] = educts[e]-1        
            else:
                break
        if e != False and isinstance(universe[e], Graph) and universe[e].id: 
            subspace.append(universe[e])
    
    # Define all products
    for e in products.keys():
        if e == 'C00080':
            h_mol = universe[e]
        while products[e] > 0:
            if e == False:
                break
            elif isinstance(universe[e], Graph) and universe[e].id:
                prod.append(universe[e])
                products[e] = products[e]-1    
            else:
                break  
        if e != False:       
            iso_check = False
            for g1 in subspace:
                if isinstance(g1, Graph) and g1.id:
                    if isinstance(universe[e], Graph) and universe[e].id:
                        if g1.isomorphism(universe[e]) == True:
                            iso_check = True
            if iso_check == False and isinstance(universe[e], Graph) and universe[e].id: 
                subspace.append(universe[e])                
   
    # add protons to Educts and Products
    with open(path+'/'+rule+'/'+rule_var,'r') as f:
        lines = f.readlines()
    string_rule = " ".join(lines)        
    left_protons,right_protons = count_protons(string_rule)
    if left_protons != 0 or right_protons != 0:
        if h_mol is not None:
            left_protons = left_protons-1
        else:
            h_mol = smiles('[H+]', name='C00080')  
            subspace.append(h_mol)             
        while left_protons > 0:
            educ.append(h_mol)
            left_protons = left_protons-1  
        while right_protons > 0:
            prod.append(h_mol)
            right_protons = right_protons-1    
  
    return educ,prod,subspace

def find_AAMs(dg,found,count_aam,count_aam_list,var_succ):
    # search for HyperEdge which creates all the products -> global atom-to-atom maps
    correct_edge = False       
    print('No. of Edges in Network',dg.numEdges)
    for hyperedge in dg.edges:
        try:
            correct_edge = dg.findEdge(educ,prod)     
        except:
            print('findEdge Faild')

    if correct_edge != False:
        if not correct_edge.isNull():
            print('HyperEdge found')        
            try:          
                # parse graphs to SMILES Format  
                res = getReactionSmiles(dg,correct_edge)        
            except:
                with open(f"./output_WildcardCase.txt", mode="a") as outp:
                    outp.write(rule+"\n")  
            print(res)
            if len(res) > 0:
                if found == False:
                    found = True
                    count_aam += 1
                    count_aam_list.append(rule)
                var_succ = True
                for s in res:
                    with open(f"./output_GlobalReactions.txt", mode="a") as outp:
                        outp.write(rule+' '+str(s)+"\n")              
    if found == False: 
        # check for partial AAMs
        for hyperedge in dg.edges:
            print('Search for Partial AAMs')    
            product_check = True
            for t in dg.products:
                if t not in prod:
                    product_check = False
            print('Result Check',product_check)       
            if product_check == True:
                try:
                    # parse graphs to SMILES Format   
                    res = getReactionSmiles(dg,hyperedge)
                except:
                    with open(f"./output_WildcardCase.txt", mode="a") as outp:
                        outp.write(rule+"\n")  
                print(res)
                if len(res) > 0:
                    if found == False:
                        found = True
                        count_aam += 1
                        count_aam_list.append(rule)
                    var_succ = True
                    for s in res:
                        print('found')
                        with open(f"./output_PartialReactions.txt", mode="a") as outp:
                            outp.write(rule+' '+str(s)+"\n")        
    return found,count_aam,count_aam_list,var_succ

def createMoleculeDatabase(reaction_data,universe):   
    for molecule in reaction_data:
        for m in molecule.keys():
            mol_kegg = get_compound_mol(m)
            s = mol_to_smiles(mol_kegg)
            if s == False:
                universe[m] = s  
            else:
                mol = smiles(s, name=m,allowAbstract=True)   
                if mol not in universe:     
                    universe[m] = mol  
    return universe

############################### START ###################################

# load reaction data form KEGG database
# reaction_data = {'reaction_name': [educts = {},products = {}]}
with open('/homes/biertank/nora/Githubs/AAMRuleTool/RCLASS2DPORules_version2.0/03/Filter Reactions/Reaction_should_be_generated.txt','r') as in_f:
#with open('/homes/biertank/nora/Githubs/AAMRuleTool/RCLASS2DPORules_version2.0/03/Filter Reactions/Test.txt','r') as in_f:
    lines = in_f.readlines()

reaction_data = {}
for line in lines:
    # Extract Reaction Name
    if line.startswith('Re'):
        line = line.split(':')
        line = line[1].split('\n')
        rn = line[0]
        continue
    # Extract Compounds an Stöchiomatric Data    
    if line.startswith('C'):
        compound_ids_str = re.search(r"Compound IDs:\[(.*?)\]", line).group(1)
        compound_ids = eval(compound_ids_str)  
        educts = {}
        products = {}
        delimiter = '<=>'
        left_side = True

        for i in range(len(compound_ids)):
            item = compound_ids[i]
            if item == delimiter:
                left_side = False
                continue
            if item.startswith('C'):
                if left_side:
                    count = 1
                    if i > 0 and compound_ids[i-1].isdigit():
                        count = int(compound_ids[i-1])
                    educts[item] = count
                else:
                    count = 1
                    if i > 0 and compound_ids[i-1].isdigit():
                        count = int(compound_ids[i-1])
                    products[item] = count
        reaction_data[rn] = [educts,products]

########################## STATS ##############################################
count_aam = 0 # counts the number of generated aams
count_miss = 0 # counts the number of reactions which could not be generated
count_aam_list = [] 
count_miss_list = [] 
##############################################################################


# load all Reaction rules witch AAMs should be generated
#path = '/homes/biertank/nora/Githubs/AAMRuleTool/RCLASS2DPORules_version2.0/03/03_DPO_bigRun'
path = '/homes/biertank/nora/Githubs/AAMRuleTool/RCLASS2DPORules_version2.0/03/03_DPO_Rules'
#path = '/homes/biertank/nora/Githubs/AAMRuleTool/RCLASS2DPORules_version2.0/03/Test'
rules_names = os.listdir(path)        
universe = {} # list of all known Molecules by KEGG DB  

# create a list of DPO-Rules which are already tried
done_list = []
with open('output_rxns.txt','r') as f:
    lines = f.readlines()    
for i in lines:
    done_list.append(i.strip())


for rule in rules_names:  
    if rule in done_list:
        continue
    # complete universe       
    universe = createMoleculeDatabase(reaction_data[rule],universe)
    print('universe',universe)
    
    # write logbook
    with open(f"./output_rxns.txt", mode="a") as outp:
        outp.write(rule+'\n')        
    
    found = False # marker if a aam could be generated    
    rule_variations = os.listdir(path+'/'+rule) # find all rules for one reaction     
    for rule_var in rule_variations:
        print('TRY ', rule_var)
        var_succ = False # marker if a variation of a rule could generated a aam
        
        # add rule
        b = [ruleGML(path+'/'+rule+'/'+rule_var, add=False,printStereoWarnings=False)]


        # generated educts and products
        educ,prod,subspace = create_Educts_and_Products(reaction_data,path,rule,rule_var,universe)
        print('Educts of the Reaction: ',educ)
        print('Products of the Reaction: ',prod)  

        ## print in file
        #post.summarySection("Educt Graph")
        #for a in educ:   
        #    a.print()      
        #post.summarySection("Product Graph")
        #for a in prod:   
        #    a.print()
        
        # Create Network        
        dg = DG(labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation),graphDatabase=educ)
        dg.build().apply(subspace, ruleGML(path+'/'+rule+'/'+rule_var, add=False,printStereoWarnings=False), onlyProper=False)
            
        # print in file
        #post.summarySection("Network")     
        #dg.print()
        
        #post.summarySection("Database")     
        #for g in dg.graphDatabase:
        #    g.print()
            
        #post.summarySection("Rule")  
        #for r in b:
        #    r.print()                
        
        # find atom-to-atom maps in Network
        found,count_aam,count_aam_list,var_succ = find_AAMs(dg,found,count_aam,count_aam_list,var_succ)

        if found == False:
            print('### HyperEdge not found! Check reverse Modus ###') 
            # try the opposite direction            
            a = [ruleGML(path+'/'+rule+'/'+rule_var,invert=True, add=False,printStereoWarnings=False)]
            
            # generated educts and products
            educ,prod,subspace = create_Educts_and_Products(reaction_data,path,rule,rule_var,universe)
            print('Educts of the Reaction: ',educ)
            print('Products of the Reaction: ',prod)  
            print('Database for Network: ',subspace)          
            
            # Create Network
            dg = DG(labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation),graphDatabase=educ)
            dg.build().apply(subspace, ruleGML(path+'/'+rule+'/'+rule_var,invert=True, add=False,printStereoWarnings=False), onlyProper=False)
                                
            #post.summarySection("Network Reverse")     
            #dg.print()
            correct_edge = False         
            
            # find atom-to-atom maps in Network
            found,count_aam,count_aam_list,var_succ = find_AAMs(dg,found,count_aam,count_aam_list,var_succ)

        if found == False:
            with open(f"./output_FalseRuleVariation.txt", mode="a") as outp:
                outp.write(rule_var+'\n')     
                #post.summarySection("Network False")     
                #dg.print()
    if found == False:
        count_miss += 1
        count_miss_list.append(rule)
        with open(f"./output_FalseReactions.txt", mode="a") as outp:
            outp.write(rule+'\n')            

print('######################## SUMMARY #############################')
print('Total No. of Reactions: ',len(rules_names))
print('No. of generated AAMs: ',count_aam, count_aam_list)
print('No. of RXN which faild: ',count_miss, count_miss_list)
print('#####################################################')
