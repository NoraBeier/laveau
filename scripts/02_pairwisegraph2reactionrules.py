import os
import networkx as nx
import matplotlib.pyplot as plt
import requests
import re
from itertools import product
from rdkit import Chem
from rdkit.Chem import AllChem
from networkx.algorithms import isomorphism as iso
from networkx.algorithms.isomorphism import GraphMatcher
from collections import OrderedDict
from collections import defaultdict
import time

def print_Rule(graph,label):
    pos = nx.spring_layout(graph)
    labels = nx.get_node_attributes(graph, label)
    nx.draw_networkx_labels(graph, pos,labels)
    nx.draw(graph, pos)
    nx.draw_networkx_edge_labels(graph, pos)
    plt.show() 

def node_match1(a, b):
    return a['label'] == b['label']

def node_match2(a, b):
    if b['label'] == 'R':
        return True
    else:
        return a['label'] == b['label']    

def node_match_typ(a, b):
    return a['atomtype'] == b['atomtype']
         
def edge_match(a, b):
    return a['label'] == b['label']

# check for perutaion of matching nodes (most because of hydrogens)       
def delete_sub_mappings(mapping):                
    unique_mappings = []
    seen_keys = set()
    for d in mapping:
        sorted_keys = tuple(sorted(d.keys()))
        if sorted_keys in seen_keys:
            continue
        seen_keys.add(sorted_keys)
        unique_mappings.append(d)                                                                
    return unique_mappings  

# finds all molecular subgraphs which were earlier created an load them in left an right lists
def load_graphs_from_folders(root_path, keyword,comp_pair):
    
    # translate the MOD bond names to the mol bond names 
    def update_bound_attribute(G):
        for u, v, data in G.edges(data=True):
            bound_value = data['bound']
            if bound_value == '-' or bound_value == 'r':
                data['label'] = 1.0
            elif bound_value == '=' or bound_value == 'dr':
                data['label'] = 2.0
            elif bound_value == ':':
                data['label'] = 1.5
            elif bound_value == '#':
                data['label'] = 3
                       
        return G
              
    rclass_left_comp1 = {}
    rclass_right_comp1 = {} 
    rclass_left_comp2 = {}
    rclass_right_comp2 = {} 

    sub_count = 0    
    for dirpath, _, filenames in os.walk(root_path):
        for filename in filenames:
            if keyword in filename:
                # load subgraph    
                filepath = os.path.join(dirpath, filename)
                sub_no = re.search(r'_(\d+)\.gml', filename)
                list_name = sub_no.group(1)
                G = nx.read_gml(filepath)              
                
                # for counting no of subgraphs per side 
                if 'left' in filename:
                    if int(sub_no.group(1)) > sub_count:
                        sub_count = int(sub_no.group(1))
                                       
                # rename all nodes so that every id is unique
                rename = {node : list_name +'_'+str(node) for node in G.nodes()} 
                nx.relabel_nodes(G, rename, copy=False)
                if True:
                    #  modify the edges and nodes so that it fits to the mol data
                    update_bound_attribute(G)

                    for node in G.nodes():
                        if G.has_node(node) and 'atom' in G.nodes[node]:
                            if '+' in G.nodes[node]['atom'] or '-' in G.nodes[node]['atom']:
                                G.nodes[node]['label'] = G.nodes[node]['atom'].replace('+', '').replace('-', '') 
                            else:
                                G.nodes[node]['label'] = G.nodes[node]['atom']
                                       
                    # check compound 1                
                    mapping = subgraph_isomorphism_with_attributes(comp_pair[0],G,node_match2,edge_match)
                    if len(mapping) != 0:
                        mapping = delete_sub_mappings(mapping)
                        if 'left' in filename:       
                            if list_name not in rclass_left_comp1:
                                rclass_left_comp1[list_name] = []
                            rclass_left_comp1[list_name].append((G,mapping))  
                        elif 'right' in filename:         
                            if list_name not in rclass_right_comp1:
                                rclass_right_comp1[list_name] = []
                            rclass_right_comp1[list_name].append((G,mapping))    
                    # check comound 2 
                    mapping = subgraph_isomorphism_with_attributes(comp_pair[1],G,node_match2,edge_match)        
                    if len(mapping) != 0:
                        mapping = delete_sub_mappings(mapping)
                        if 'left' in filename:
                            if list_name not in rclass_left_comp2:
                                rclass_left_comp2[list_name] = []
                            rclass_left_comp2[list_name].append((G,mapping))  
                        elif 'right' in filename:                           
                            if list_name not in rclass_right_comp2:
                                rclass_right_comp2[list_name] = []
                            rclass_right_comp2[list_name].append((G,mapping))  

    # check if all subgraphs where found in the compound
    if len(rclass_left_comp1) != sub_count:
        rclass_left_comp1 = {}
    if len(rclass_left_comp2) != sub_count:
        rclass_left_comp2 = {}
    if len(rclass_right_comp1) != sub_count:
        rclass_right_comp1 = {}
    if len(rclass_right_comp2) != sub_count:
        rclass_right_comp2 = {}

    return rclass_left_comp1, rclass_right_comp1, rclass_left_comp2, rclass_right_comp2

# load compound mol format form KEGG Database
def get_compound_mol(compound_id):
    url = f'http://rest.kegg.jp/get/compound:{compound_id}/mol'
    try:
        response = requests.get(url)
        response.raise_for_status()  # This throws an exception if the status code is not 200 (OK)
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Fehler bei der Anfrage: {e}")
        return None

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

# find all matches of the RDM Patterns in the Molecule  
def subgraph_isomorphism_with_attributes(G, s,node_match,edge_match):  

    # check if Subgraphs has polaritis if yes then remove then, otherwise isomorphism test will not work because mol format hasnt it                      
    GM = iso.GraphMatcher(G, s,node_match=node_match,edge_match=edge_match)
    subgraph_matches = list(GM.subgraph_monomorphisms_iter())       
    
    return subgraph_matches

                    
def connect_patterns(comp,pattern):
    conected_graphs = [] # List of already generated patterns, avoid duplications
    if len(pattern) != 1:
        # try all combinaions of subgraph mappings
        mapping_combis = list(product(pattern[0][1],pattern[1][1]))
        mapping_combis = list(product(*[pattern[i][1] for i in range(len(pattern))]))

        for combi in mapping_combis:       
            # time check for to comlex RCLASSES
            current_time = time.time()
            elapsed_time = current_time - start_time
            time_limit = 3600
            if elapsed_time >= time_limit:
                return conected_graphs       
            # find all node in the compound which should be in the rule
            nodes_inRule =  {}
            rdm_pattern = comp.copy()
            for maps in combi:
                for key in maps.keys():
                    nodes_inRule[key] = []
            for maps in combi:
                for key in maps.keys():
                    nodes_inRule[key].append(maps[key])

            # cut RDM Pattern out of molecule and transfer node attribut from subgraphs to node
            for node in comp.nodes():
                if node in nodes_inRule.keys():
                    for sub_node in nodes_inRule[node]: 
                        for g in pattern:                     
                            if sub_node in g[0].nodes(data=True):
                                label_att = g[0].nodes[sub_node]['atom']
                                if 'map' in rdm_pattern.nodes[node]:
                                    map_att_list = [ g[0].nodes[sub_node]['map'] ]
                                    map_att_list.extend(rdm_pattern.nodes[node]['map'])
                                else:
                                    map_att_list = [g[0].nodes[sub_node]['map']]
                                if 'atomtype' in rdm_pattern.nodes[node]:
                                    atomtype_att = [g[0].nodes[sub_node]['atomtype']]
                                    atomtype_att.extend(rdm_pattern.nodes[node]['atomtype'])
                                else:
                                    atomtype_att = [g[0].nodes[sub_node]['atomtype']]
                                
                                rdm_pattern.nodes[node]['map'] = map_att_list
                                rdm_pattern.nodes[node]['atomtype'] = atomtype_att
                                rdm_pattern.nodes[node]['label'] = label_att
                else:
                    rdm_pattern.remove_node(node)
      
            # if there are R-Atoms with highter degrees then 1, rewrite the atom label with the atom of the compound         
            for node,data in rdm_pattern.nodes(data=True):        
                if data['label'] == 'R':
                    if rdm_pattern.degree(node) > 1:
                        label_att = comp.nodes[node]['label']
                        rdm_pattern.nodes[node]['label'] = label_att
                                                   
            # check if identical graphs was already generated
            iso_check = True         
            for i in conected_graphs:
                if nx.is_isomorphic(i,rdm_pattern,node_match=node_match1,edge_match=edge_match): 
                    if nx.is_isomorphic(i,rdm_pattern,node_match=node_match_typ,edge_match=edge_match): 
                        iso_heck = False
            if iso_check == True:
                conected_graphs.append(rdm_pattern)
                     
    else:
            # cut RDM Pattern out of molecule and transfer node attribut from subgraphs to node
            rdm_pattern = comp.copy()
            for node in comp.nodes():
                if node in pattern[0][1][0].keys():  
                    label_att = pattern[0][0].nodes[pattern[0][1][0][node]]['atom']
                    map_att = pattern[0][0].nodes[pattern[0][1][0][node]]['map']
                    map_list = [map_att]
                    atomtype_att = pattern[0][0].nodes[pattern[0][1][0][node]]['atomtype']
                    rdm_pattern.nodes[node]['map'] = map_list
                    rdm_pattern.nodes[node]['atomtype'] = atomtype_att
                    rdm_pattern.nodes[node]['label'] = label_att
                else:
                    rdm_pattern.remove_node(node)
                   
            conected_graphs.append(rdm_pattern)   
          
    return conected_graphs            
   
# check if both side are represented in the compound pair
def matchingMolecule(comp_pair,variation):

    filtered_list_var1 = {}
    filtered_list_var1['left'] = []
    filtered_list_var1['right'] = []    
    filtered_list_var2 = {}
    filtered_list_var2['left'] = []
    filtered_list_var2['right'] = []
    left_mol1 = []
    left_mol2 = []
    right_mol1 = []
    right_mol2 = []       
  
    # try it for both compounds 
    for com in range(len(comp_pair)):              
            
        # try left and right reaction side
        for side in variation.keys():
            # try it for all created patterns            
            for g in range(len(variation[side])):     
                
                s = variation[side][g].copy()          
                gm = GraphMatcher(comp_pair[com],s,node_match=node_match2,edge_match=edge_match)
                mol_check = gm.subgraph_is_isomorphic() 

                if side == 'left':
                    if com == 0 and mol_check == True:    
                        left_mol1.append(g)
                    if com == 1 and mol_check == True:                                
                        left_mol2.append(g)
                if side == 'right':    
                    if com == 0 and mol_check == True:    
                        right_mol1.append(g)
                    if  com == 1 and mol_check == True:                                
                        right_mol2.append(g)                 
        
    if len(left_mol1) == 0: 
        right_mol2 = []
    if len(left_mol2) == 0: 
        right_mol1 = []
    if len(right_mol1) == 0: 
        left_mol2 = []
    if len(right_mol2) == 0: 
        left_mol1 = []
           
    for l1 in left_mol1:
        if variation['left'][l1] not in filtered_list_var1['left']:
            filtered_list_var1['left'].append(variation['left'][l1])
    for r2 in right_mol2:
        if variation['right'][r2] not in filtered_list_var1['right']:
            filtered_list_var1['right'].append(variation['right'][r2])      
    for r1 in right_mol1:
        if variation['right'][r1] not in filtered_list_var2['left']:
            filtered_list_var2['left'].append(variation['right'][r1])
    for l2 in left_mol2:
        if variation['left'][l2] not in filtered_list_var2['right']:
            filtered_list_var2['right'].append(variation['left'][l2])          
    # Check for R-M-D-Atom consistency
    if len(filtered_list_var1['left']) != 0 and len(filtered_list_var2['left']) == 0: 
    
        variation = consistency_check(filtered_list_var1,comp_pair,True)       
    elif len(filtered_list_var1['left']) == 0 and len(filtered_list_var2['left']) != 0: 
        variation = consistency_check(filtered_list_var2,comp_pair,True)       
    else:   
        filtered_list_var1 = consistency_check(filtered_list_var1,comp_pair,True) 
        filtered_list_var2 = consistency_check(filtered_list_var2,comp_pair,True)  

        # merge both variations together
        variation= filtered_list_var1  
        for l2 in filtered_list_var2['left']:
            variation['left'].append(l2)
        for r2 in filtered_list_var2['right']:
            variation['right'].append(r2)            
    return variation
    
# checks if a rule variation is only a sub graph of a other rule variation     
def subgraph_check(rule):
    for side in rule.keys():
        filtered_list = []
        filtered_list.append(rule[side][0]) 
        for i in range(1,len(rule[side])):
            check = True
            for j in range(len(filtered_list)):
                if nx.is_isomorphic(rule[side][i],filtered_list[j],node_match=node_match2,edge_match=edge_match):
                   check = False 
                   break
            
            if check == True:                     
                filtered_list.append(rule[side][i])
            
            # if the new Graph has more information (means less R-Atoms) than use that
            if check == False: 
                R_new = sum(1 for node, data in rule[side][i].nodes(data=True) if data.get('label') == 'R')
                R_old = sum(1 for node, data in filtered_list[j].nodes(data=True) if data.get('label') == 'R') 
                if R_new < R_old:
                    filtered_list[j] = rule[side][i]              

        rule[side] = filtered_list
    return rule
    
def consistency_check(rule,comp,new_mapping):
    filtered_list = {}
    filtered_list['left'] = []
    filtered_list['right'] = []       

    # create a list of the atoms with there bonding from the neighborhood
    def get_edge_attributes_for_node(G, node):
        edge_attributes = []
        for edge in G.edges(node, data=True):
            u, v, attr = edge
            if u == node:
                edge_attributes.append([G.nodes[v]['label'],attr['label']])
            if v == node:
                edge_attributes.append([G.nodes[u]['label'],attr['label']])   
        return edge_attributes

    # check for same neighbor nodes
    def neighborhood_check(neighborhood1,neighborhood2):
        combined_hist = {}

        itemset = set()
        for n1 in neighborhood1:
            if n1[0] != 'H':
                itemset.add((n1[0], n1[1]))
        for n2 in neighborhood2:
            if n2[0] != 'H':
                itemset.add((n2[0], n2[1]))  
        for item in itemset:
            combined_hist[item] = 0
            
        for n1 in neighborhood1:
             if n1[0] != 'H':
                combined_hist[(n1[0], n1[1])] = combined_hist[(n1[0], n1[1])] + 1
        for n2 in neighborhood2:
            if n2[0] != 'H':
                combined_hist[(n2[0], n2[1])] = combined_hist[(n2[0], n2[1])] - 1

        no_match_n1 = {}
        no_match_n2 = {}
        for n in combined_hist:
            if combined_hist[n] > 0:
                no_match_n1[n] = combined_hist[n]
            elif combined_hist[n] < 0:
                no_match_n2[n] = -1*combined_hist[n]

        has_atom_n1 = any( [ n[0] != "R" for n in no_match_n1 ] )
        has_atom_n2 = any( [ n[0] != "R" for n in no_match_n2 ] )

        if has_atom_n1 and has_atom_n2:
            return False

        bond_hist = {1.0: 0, 2.0: 0, 1.5: 0,3.0: 0}
        for n1 in no_match_n1:
            bond_hist[n1[1]] = bond_hist[n1[1]] + no_match_n1[n1]
        for n2 in no_match_n2:
            bond_hist[n2[1]] = bond_hist[n2[1]] - no_match_n2[n2]     
        for b in bond_hist:
            if bond_hist[b] != 0:
                return False     
        return True
   
    for var in rule['left']:
        neighborhood_check_left = False
        for r_node in var.nodes():
            if 'r' in var.nodes[r_node]['atomtype']:
                neighborhood_pattern = get_edge_attributes_for_node(var,r_node)
                if new_mapping == True:
                    mapping = subgraph_isomorphism_with_attributes(comp[0],var,node_match2,edge_match)    
                    if len(mapping) != 0:
                        mapping = delete_sub_mappings(mapping)
                    for m in mapping[0].keys():
                        if mapping[0][m] == r_node:
                            neighborhood_molecule = get_edge_attributes_for_node(comp[0],m) 
                else:
                    neighborhood_molecule = get_edge_attributes_for_node(comp[0],r_node)

                if  neighborhood_check(neighborhood_pattern,neighborhood_molecule) == True:
                        neighborhood_check_left = True
                else:
                    neighborhood_check_left = False 
                    break               
        if neighborhood_check_left == True:
            filtered_list['left'].append(var)
            
    for var in rule['right']:                    #        
        neighborhood_check_right = False
        for r_node in var.nodes():

            if 'r' in var.nodes[r_node]['atomtype']:
                neighborhood_pattern = get_edge_attributes_for_node(var,r_node)

                if new_mapping == True:
                    mapping = subgraph_isomorphism_with_attributes(comp[1],var,node_match2,edge_match)    
                    if len(mapping) != 0:
                        mapping = delete_sub_mappings(mapping)
                    for m in mapping[0].keys():
                        if mapping[0][m] == r_node:
                            neighborhood_molecule = get_edge_attributes_for_node(comp[1],m)       
                else:           
                    neighborhood_molecule = get_edge_attributes_for_node(comp[1],r_node)    

                if neighborhood_check(neighborhood_pattern,neighborhood_molecule) == True:  
                    neighborhood_check_right = True
                else:
                    neighborhood_check_right = False
                    break               
 
        if neighborhood_check_right == True:
            filtered_list['right'].append(var)
     
    return filtered_list

# checks if the diffrent mapping only comes from diffrent order of keys
def check_mappings(rclass):
    # try for every subgraph
    for subgraph in rclass.keys():  

        # try for every variation of a subgraph
        for var in range(len(rclass[subgraph])-1):
            filtered_maps = []
            filtered_maps.append(rclass[subgraph][var][1][0])
            # check all mappings
            for m in range(len(rclass[subgraph][var][1])-1):    
                if sorted(rclass[subgraph][var][1][m].keys()) != sorted(rclass[subgraph][var][1][m+1].keys()): 
                    filtered_maps.append(rclass[subgraph][var][1][m+1])   
                        
            rclass[subgraph][var] = (rclass[subgraph][var][0],filtered_maps)                
 
    return rclass
    
# checks if of both sides are the same amount of atoms    
def check_atom_number(variation):
 
    def count_nodes_with_condition(graph):
        count = 0
        for node, data in graph.nodes(data=True):
            if 'label' in data and data['label'] not in ['H', 'R'] and ('atomtype' not in data or 'd' not in data['atomtype']):
                count += 1
        return count

    def count_nodes_with_DAtom(graph):
        count = 0
        for node, data in graph.nodes(data=True):
            if 'label' in data and data['label'] not in ['H', 'R']:
                count += 1
        return count
        
    filtered_leftside = []
    filtered_rightside = []
    
    left_atoms = []
    right_atoms = []
    
    for left in variation['left']:
         atom_num = count_nodes_with_condition(left)   
         left_atoms.append(atom_num)
    for right in variation['right']:
         atom_num = count_nodes_with_condition(right)   
         right_atoms.append(atom_num)

    for i in range(len(left_atoms)):
        if left_atoms[i] in right_atoms:
            filtered_leftside.append(variation['left'][i])
    for i in range(len(right_atoms)):
        if right_atoms[i] in left_atoms:
            filtered_rightside.append(variation['right'][i])

    if len(filtered_leftside) == 0 or len(filtered_rightside) == 0:
        left_atoms = []
        right_atoms = []
        for left in variation['left']:
             atom_num = count_nodes_with_DAtom(left)   
             left_atoms.append(atom_num)
        for right in variation['right']:
             atom_num = count_nodes_with_DAtom(right)   
             right_atoms.append(atom_num)

        for i in range(len(left_atoms)):
            if left_atoms[i] in right_atoms:
                filtered_leftside.append(variation['left'][i])
        for i in range(len(right_atoms)):
            if right_atoms[i] in left_atoms:
                filtered_rightside.append(variation['right'][i])
    else:         
        variation['left'] = filtered_leftside
        variation['right'] = filtered_rightside
    
    return variation
    
# counter for stats
counter_work = 0
counter_notwork = 0   
counter_input = 0

# read the linked names of componands for each RCLASS
with open('./Additional_ Files/RCLASS_RPAIR.txt', 'r') as f:
    lines = f.readlines()
comp_list = {}    
for line in lines:
    if line.startswith("RCLASS:"):
        rxn = line.split("RCLASS:")[1].strip()
    elif line.startswith("RPAIRs:"):
        rpair_str = line.split("RPAIRs:")[1].strip()
        rpair_l = rpair_str.split(", ")
        rpair_list = [l.split("_") for l in rpair_l]
        if rxn in comp_list:
            comp_list[rxn].append(rpair_list)
        else:
            comp_list[rxn] = rpair_list

# read list with undefined Atom RCLASSES
with open('./Additional_ Files/List_UndefindAoms.txt', 'r') as f:
    lines = f.readlines()
undefind_list = [s.strip() for s in lines]

root_path = sys.argv[2]
output_path = './02_Reaction Rules/'

# try all RCLASSES from the comound_pairs
for rxn in comp_list:
        # time counter which continue the process if the RCLASS take to much time
        start_time = time.time()
        # pass the RCLASSES which as a undefind Atom in the reaction
        if rxn not in undefind_list: 
            counter_input = counter_input + 1

            # translate all mol compound data from KEGG in graphs
            comp_graphs = [] # Compound Graphs as Tuples 
            compound_ids = comp_list[rxn]
            for mol in compound_ids:
                m1 = get_compound_mol(mol[0])
                g1 = mol_to_graph(m1)
                m2 = get_compound_mol(mol[1])
                g2 = mol_to_graph(m2)
                comp_graphs.append([g1,g2])

            # load the molecular subgraphs 
            rclass_left_comp1, rclass_right_comp1, rclass_left_comp2, rclass_right_comp2 = load_graphs_from_folders(root_path, rxn, comp_graphs[0])

            # if one side could only be found on one compound the other side matches have to be false positiv -> dont check further 
            if len(rclass_left_comp1) != 0 and len(rclass_left_comp2) == 0:
                 rclass_right_comp1 = []    
                 rclass_left_comp1 = check_mappings(rclass_left_comp1 )
            if len(rclass_left_comp1) == 0 and len(rclass_left_comp2) != 0:
                 rclass_right_comp2 = [] 
                 rclass_left_comp2 = check_mappings(rclass_left_comp2)
            if len(rclass_right_comp1) == 0 and len(rclass_right_comp2) != 0:
                 rclass_left_comp2 = []      
                 rclass_right_comp2 =  check_mappings(rclass_right_comp2)
            if len(rclass_right_comp1) != 0 and len(rclass_right_comp2) == 0:
                 rclass_left_comp1 = []      
                 rclass_right_comp1 =  check_mappings(rclass_right_comp1)
            if len(rclass_right_comp1) != 0 and len(rclass_right_comp2) != 0:
                rclass_left_comp1 = check_mappings(rclass_left_comp1)
                rclass_left_comp2 = check_mappings(rclass_left_comp2)
                rclass_right_comp2 =  check_mappings(rclass_right_comp2)
                rclass_right_comp1 =  check_mappings(rclass_right_comp1)
            if (len(rclass_right_comp1) == 0 and len(rclass_right_comp2) == 0) or (len(rclass_left_comp1) == 0 and len(rclass_left_comp2) == 0):
                counter_notwork = counter_notwork + 1
                with open('List_ErrorRCLASS', 'a') as file:
                    file.write(rxn + ' Compound checks faild! \n')
                continue   

            # combine all possibile combinations of Subgraphs
            if len(rclass_left_comp1) != 0:
                subcombis_left_comp1 = list(product(*rclass_left_comp1.values()))
            else:
                subcombis_left_comp1 = []
            if len(rclass_left_comp2) != 0:
                subcombis_left_comp2 = list(product(*rclass_left_comp2.values()))
            else:
                subcombis_left_comp2 = []
            if len(rclass_right_comp1) != 0:
                subcombis_right_comp1 = list(product(*rclass_right_comp1.values())) 
            else: 
                subcombis_right_comp1 = []
            if len(rclass_right_comp2) != 0:
                subcombis_right_comp2 = list(product(*rclass_right_comp2.values())) 
            else:
                subcombis_right_comp2 = []
            
            # check if the amount of one combi part is to big to handle
            continue_check = [len(subcombis_left_comp1),len(subcombis_left_comp2),len(subcombis_right_comp1),len(subcombis_right_comp2)]
           
            # try to puzzle the subgraphs to one RCLASS
            rdm_pattern_left_comp1 =  []
            rdm_pattern_right_comp1 =  []
            rdm_pattern_left_comp2 =  []
            rdm_pattern_right_comp2 =  []       
            
            if len(subcombis_left_comp1) != 0:
                for pattern in subcombis_left_comp1:
                    pattern_list = connect_patterns(comp_graphs[0][0],pattern)
                    for i in pattern_list:
                        rdm_pattern_left_comp1.append(i)
               
            if len(subcombis_left_comp2) != 0:
                for pattern in subcombis_left_comp2:
                    pattern_list = connect_patterns(comp_graphs[0][1],pattern)
                    for i in pattern_list:
                        rdm_pattern_left_comp2.append(i)
            if len(subcombis_right_comp1) != 0:
                for pattern in subcombis_right_comp1:
                    pattern_list = connect_patterns(comp_graphs[0][0],pattern)
                    for i in pattern_list:
                        rdm_pattern_right_comp1.append(i)            
            if len(subcombis_right_comp2) != 0:
                for pattern in subcombis_right_comp2:
                    pattern_list = connect_patterns(comp_graphs[0][1],pattern)
                    for i in pattern_list:
                        rdm_pattern_right_comp2.append(i)

            # merge the partner RDM Pattern together
            variation1 = {}
            variation1['left'] = rdm_pattern_left_comp1
            variation1['right'] = rdm_pattern_right_comp2
            variation2 = {}
            variation2['left'] = rdm_pattern_right_comp1
            variation2['right'] = rdm_pattern_left_comp2

            # Remove the polarity of the atoms
            for graph in variation1['right']: 
                for node, data in graph.nodes(data=True):
                    if 'label' in data:
                        data['label'] = data['label'].replace('+', '').replace('-', '')
            for graph in variation1['left']: 
                for node, data in graph.nodes(data=True):
                    if 'label' in data:
                        data['label'] = data['label'].replace('+', '').replace('-', '')  
            for graph in variation2['right']: 
                for node, data in graph.nodes(data=True):
                    if 'label' in data:
                        data['label'] = data['label'].replace('+', '').replace('-', '')
            for graph in variation2['left']: 
                for node, data in graph.nodes(data=True):
                    if 'label' in data:
                        data['label'] = data['label'].replace('+', '').replace('-', '')  
            

            # do a consisteny check on the fist compound
            variation1 = consistency_check(variation1,comp_graphs[0],False)
            variation2 = consistency_check(variation2,comp_graphs[0],False) 
                        
            variation1 = check_atom_number(variation1)
            variation2 = check_atom_number(variation2)
                   
            continue_check = [len(variation1['left']),len(variation1['right']),len(variation2['left']),len(variation2['right'])]
            if any(x > 100 for x in continue_check) :   
                with open('List_BigRCLASSES.txt','a') as err:
                    err.write(rxn+' ')
                continue
            
            for mo in comp_graphs:
                if len(variation1['left']) != 0 and len(variation1['right']) != 0:
                    variation1 = matchingMolecule(mo,variation1)
                    
                if len(variation2['left']) != 0 and len(variation2['right']) != 0:   
                    variation2 = matchingMolecule(mo,variation2) 

            variation1 = check_atom_number(variation1)
            variation2 = check_atom_number(variation2)
            
            # Case 1: variation1 is complete empty and variation2 works > have to be variation2
            if len(variation1['left']) == 0 or len(variation1['right']) == 0:
                if len(variation2['left']) != 0 and len(variation2['right']) != 0:
                    rule = variation2
                    # sub graph check: sometimes the differences are only on hydrogen node less, delete this one
                    rule = subgraph_check(rule)
                else:
                    # Error: non of the variation are suitable 
                    counter_notwork = counter_notwork + 1
                    rule = None
                    with open('List_ErrorRCLASS', 'a') as file:
                        file.write(rxn + ' Compound checks faild! \n')
                    continue    
            # Case 2: Variation 1 works and variation2 is empty -> have to be variation1
            elif len(variation2['left']) == 0 or len(variation2['right']) == 0:   
                    rule = variation1                 
                     # sub graph check: sometimes the differences are only on hydrogen node less, delete this one
                    rule = subgraph_check(rule)

            # Case 3: in both variations are still working         
            elif len(variation2['left']) != 0 and len(variation2['right']) != 0:  
                if len(variation1['left']) != 0 and len(variation1['right']) != 0:
                    rule = 'both' 
            else:

                # Case 4: all Variations are empty! non of the variation are suitable 
                counter_notwork = counter_notwork + 1
                rule = None
                with open('List_ErrorRCLASS', 'a') as file:
                    file.write(rxn + ' Compound checks faild!\n')        
                continue 
                
            # Networkx overwrites the label attribute with the atom information, so this must be stored again in a separate attribute  
            def create_atomAttribute(g):     
                for node, data in g.nodes(data=True):
                    if 'label' in data:
                        data['atom'] = data['label'] 
                return g

            # save rules
            new_path = output_path + rxn
            if not os.path.exists(new_path):
                os.makedirs(new_path)
                      
            if rule == variation1 or rule == variation2:
                counter_work = counter_work + 1 
                count_rules = 1

                for r1 in rule['left']: 
                    r1 = create_atomAttribute(r1)
                    filename = new_path + '/' + rxn + '_0' + str(count_rules) + "_{'left'}.gml"
                    nx.write_gml(r1, filename)
                    count_rules = count_rules+1
     
                count_rules = 1           
                for r2 in rule['right']: 
                    r2 = create_atomAttribute(r2)
                    filename = new_path + '/' + rxn + '_0' + str(count_rules) + "_{'right'}.gml"
                    nx.write_gml(r2, filename)      
                    count_rules = count_rules+1 
                            
            if rule == 'both':               
                counter_work = counter_work + 1
                
                v1 = subgraph_check(variation1)
                count_rules = 1
                for r1 in v1['left']: 
                    r1 = create_atomAttribute(r1)
                    filename = new_path + '/' + rxn + '_0' + str(count_rules) + "_{'left'}.gml"
                    nx.write_gml(r1, filename)
                    count_rules = count_rules+1
                for r1 in v1['right']: 
                    r1 = create_atomAttribute(r1)
                    filename = new_path + '/' + rxn + '_0' + str(count_rules) + "_{'right'}.gml"
                    nx.write_gml(r1, filename)
                    count_rules = count_rules+1
                  
                v2 = subgraph_check(variation2)    
                count_rules = 1           
                for r2 in v2['left']: 
                    r2 = create_atomAttribute(r2)
                    filename = new_path + '/' + rxn + '_0' + str(count_rules) + "_{'left'}.gml"
                    nx.write_gml(r2, filename)      
                    count_rules = count_rules+1 
                for r2 in v2['right']: 
                    r2 = create_atomAttribute(r2)
                    filename = new_path + '/' + rxn + '_0' + str(count_rules) + "_{'right'}.gml"
                    nx.write_gml(r2, filename)
                    count_rules = count_rules+1 
                            
print('No. of input Rules: ', counter_input)
print('No. of generierted Rules: ', counter_work)
print('No. of Errors: ', counter_notwork)                
                
