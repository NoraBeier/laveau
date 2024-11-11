from networkx.algorithms import isomorphism as iso
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit.Chem import AllChem
from rdkit import Chem
import itertools
from collections import Counter
from itertools import product
from itertools import islice
import networkx as nx
import copy
import matplotlib.pyplot as plt
import requests
import glob
import re
import os


#This function takes a nested list and returns a flat list.
def flatten(nested_list):

    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))
        else:
            flat_list.append(item)

    return flat_list

# load compound mol format form KEGG Database
def get_compound_mol(compound_id):
    url = f'http://rest.kegg.jp/get/compound:{compound_id}/mol'
    try:
        response = requests.get(url)
        response.raise_for_status()  # This throws an exception if the status code is not 200 (OK)
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Error in the request: {e}")
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
 
def node_match(a, b):
    if b['atom'] == 'R' and a['label'] != 'H':
        return True
    else:
        return a['label'] == b['atom']   

def node_matchD(a, b):
    if b['label'] == 'H':
        return True
    return a['label'] == b['label']         
 
 # find all matches of the RCLass in the Molecule  
def subgraph_isomorphism_with_attributes(M,R,node_match_opt):
  
    def edge_match(a, b):
        return str(a['label']) == str(b['label'])       

    GM = iso.GraphMatcher(M,R,node_match=node_match_opt,edge_match=edge_match)
    subgraph_matches = list(GM.subgraph_monomorphisms_iter())       
    return subgraph_matches

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
 
# adds the mapping from the RCLASS Graphs to the Reaction
def add_mapping(mappings,comp,rc,list_mapping_num,sides):
    # copy maps and atomtypes vom RCLASS to Compound Graph
    mapped_comp = [comp]
    mapping_no = 0
    right_side = False
    for mapping in mappings:
        L = copy.deepcopy(comp)    
        
        # if it is the right side, take the mappings from the left side
        if len(list_mapping_num) == 0:
            right_side = True
        for n in mapping.keys():   
            if not isinstance(rc.nodes[mapping[n]]['atomtype'],list):
                rc.nodes[mapping[n]]['atomtype'] = [rc.nodes[mapping[n]]['atomtype']]     
            L.nodes[n]['atomtype'] = rc.nodes[mapping[n]]['atomtype']           
            # add matched RCLASS comp_pair
            L.nodes[n]['comp_pair'] = sides 
            L_maps = []
            if right_side == True:
                for i in rc.nodes[mapping[n]]['map']:                    
                    L_maps.append(str(i)+'_' + str(mapping_no))
                    if mapping_no not in list_mapping_num:
                        list_mapping_num.append(mapping_no)
            else:
                if len(list_mapping_num) != 0:
                    mapping_no = list_mapping_num.pop(0)  
                    right_side = True
                for i in rc.nodes[mapping[n]]['map']:           
                        L_maps.append(str(i)+'_' + str(mapping_no))                                                
            L.nodes[n]['map'] = L_maps       
        mapping_no = mapping_no +1      
        if L not in mapped_comp:
            mapped_comp.append(L)    
     
    # combine double Mapping on one Molecule
    mapping_no_copy = mapping_no
    if len(mapped_comp) > 2:
        combinations = []
        for r in range(1, len(mapped_comp) + 1):
            max_combinations = 1000 
            combinations.extend(itertools.islice(itertools.combinations(mapped_comp, r), max_combinations))
        result = [list(comb) for comb in combinations]
        mapped_comp = []
        for comp in result:
            mapping_no = mapping_no_copy
            if len(comp) > 1:   
                H = copy.deepcopy(comp[0]) 
                for g in range(1,len(comp)):
                    G = copy.deepcopy(comp[g])
                    for n in G.nodes():
                        data = copy.deepcopy(H.nodes[n])
                        if 'atomtype' in G.nodes[n]:
                            if 'atomtype' in H.nodes[n]:
                                if not isinstance(H.nodes[n]['atomtype'],list):
                                    H.nodes[n]['atomtype'] = [H.nodes[n]['atomtype']] 
                                if isinstance(G.nodes[n]['atomtype'],list):
                                    H.nodes[n]['atomtype'].extend(G.nodes[n]['atomtype'])
                                else: 
                                    H.nodes[n]['atomtype'].append(G.nodes[n]['atomtype'])
  
                                # add new mapping_no to differ the same RCLASSes
                                if not isinstance(H.nodes[n]['map'],list):
                                    H.nodes[n]['map'] = [H.nodes[n]['map']] 
                                map_att_new = []  
                                for i in G.nodes[n]['map']:
                                    new_num = i.split('_')
                                    map_att_new.append(new_num[0] +'_' + str(mapping_no))     
                                H.nodes[n]['map'].extend(map_att_new)
                            else: 
                                if isinstance(G.nodes[n]['atomtype'],list):
                                    H.nodes[n]['atomtype'] = G.nodes[n]['atomtype']      
                                else:    
                                    H.nodes[n]['atomtype'] = [G.nodes[n]['atomtype']]    
                                if isinstance(G.nodes[n]['map'],list):
                                    H.nodes[n]['map'] = G.nodes[n]['map']      
                                else:    
                                    H.nodes[n]['map'] = [G.nodes[n]['map']] 
 
                    mapping_no = mapping_no+1
                mapped_comp.append(H)
            else:
                mapped_comp.append(comp[0])
        
    return mapped_comp,list_mapping_num
 
def create_subgraph_with_DAtoms(graph):
    subgraph_nodes = []
    for node in graph.nodes():
        data = graph.nodes[node]
        if 'atomtype' in data:
            if data['atomtype'] is not None:
                if 'd' in data['atomtype']:
                    subgraph_nodes.append(node)
    subgraph = graph.subgraph(subgraph_nodes)
    subgraph_edges = [(u, v) for u, v in subgraph.edges() if u in subgraph_nodes and v in subgraph_nodes]
    
    return subgraph.subgraph(subgraph_nodes)
 
def give_unique_map(graph, increment):
    for node in graph.nodes():
        data = graph.nodes[node]
        if 'map' in data:
            if isinstance(data['map'],list):
                for i in range(len(data['map'])):
                    data['map'][i] = data['map'][i]+increment                    
            else:
                data['map'] += increment

# Tautomerism function which in graphene transforms -OH to =O, and vice versa
def tautomerism(graph):
    nodes_to_remove = []
    check_OH = False
    # -OH to =O
    g = copy.deepcopy(graph) 
    for node in graph.nodes():
        if graph.nodes[node].get('label') == 'O':
            neighbors = list(graph.neighbors(node))         
            for neighbor in neighbors:
                if graph.nodes[neighbor].get('label') == 'H':
                    check_OH = True
                    nodes_to_remove.append(neighbor)
                    g.remove_node(neighbor)
                    remaining_neighbors = [n for n in g.neighbors(node) if n != neighbor]
                    for remaining_neighbor in remaining_neighbors:
                        g.edges[node, remaining_neighbor]['label'] = '='
    g.remove_nodes_from(nodes_to_remove)                    
    if not check_OH:
       # =O to -OH
       for node in graph.nodes():
            if graph.nodes[node].get('label') == 'O':
                neighbors = list(graph.neighbors(node))
                for neighbor in neighbors:
                    if graph.edges[node, neighbor].get('label') == '=':
                        g.edges[node, neighbor]['label'] = '-'
                        new_node_id = max(g.nodes) + 1
                        new_node_attributes = g.nodes[node].copy()
                        new_node_attributes['label'] = 'H'                  
                        g.add_node(new_node_id, **new_node_attributes)
                        g.add_edge(node, new_node_id, label='-')                    
                        new_nodes.append(new_node_id) 
    return g
                
def find_Datoms(Side,Otherside):
    new_Otherside = []
    for i in Otherside:
        graph_copy = copy.deepcopy(i)
        new_Otherside.append([graph_copy])
    Otherside = new_Otherside  
    for l in Side:
        # create subgraphs with D-Atoms which are not covered with other RCLASSes
        # it listet also diffrent versions: dominante D (all nodes with D), restricted D (all nodes only with D)
        d_list = []
        d_sub = create_subgraph_with_DAtoms(l)
        d_list.append(d_sub.copy())  
        d_taut = tautomerism(d_sub.copy())
        d_list.append(d_taut)
   
        for d in d_list:   
            # remove the Hydrogen atoms
            hydrogens = []
            for n in d.nodes():
                a = d.nodes[n]
                if a['label'] == 'H':
                    hydrogens.append(n)
            d.remove_nodes_from(hydrogens)
            
            # define forbitten Molecules
            forbidden_mol = []       
            for n in d.nodes():
                for m in d.nodes[n]['comp_pair']:
                    if m not in forbidden_mol:
                        forbidden_mol.append(m)

            # find D subgraphs
            if len(d.nodes()) > 0:
                for r in range(len(Otherside)):        
                    for graph in range(len(Otherside[r])):                  

                        # find Mappings of D_sub on Otherside              
                        mapping_datom = subgraph_isomorphism_with_attributes(Otherside[r][graph],d,node_matchD)                   
                        mapping_datom_mod = delete_sub_mappings(mapping_datom)
                           
                        if len(mapping_datom_mod) > 1000:
                            continue
                        for mapping in mapping_datom_mod:       
                                # check if the matching are R-/M-Atoms or only D-Atoms or no atomtypes at all
                                RM_true = False
                                D_true = False
                                noAT_true = False
                                AT_check = True
                                for node in mapping.keys():
                                    if 'atomtype' in Otherside[r][graph].nodes[node]:
                                        if Otherside[r][graph].nodes[node]['comp'] in forbidden_mol:
                                            AT_check = False
                                            break                                       
                                        if noAT_true == True:
                                            AT_check = False
                                            break         
                                        if Otherside[r][graph].nodes[node]['atomtype'] is not None:            
                                            if 'd' in Otherside[r][graph].nodes[node]['atomtype']:
                                                if RM_true ==False:
                                                    D_true = True
                                                else:
                                                    AT_check = False
                                                    break     
                                        if Otherside[r][graph].nodes[node]['atomtype'] == 'r' or Otherside[r][graph].nodes[node]['atomtype'] == 'm':    
                                            if D_true ==False:
                                                RM_true = True
                                            else:
                                                AT_check = False
                                                break         
                                    else:
                                        noAT_true == True                            
                                if AT_check == False:
                                    continue
      
                                # add the found D-Subgraphs in the Molecules
                                r_copy = copy.deepcopy(Otherside[r][graph])                                   
                                for node in mapping.keys():
                                    if 'map' in r_copy.nodes[node]:
                                        if d.nodes[mapping[node]]['map'] == r_copy.nodes[node]['map']:
                                            continue
                                        if not isinstance(d.nodes[mapping[node]]['map'],list):
                                            r_copy.nodes[node]['map'].append(d.nodes[mapping[node]]['map'])
                                        else:
                                            r_copy.nodes[node]['map'].extend(d.nodes[mapping[node]]['map'])
                                        
                                        if not isinstance(d.nodes[mapping[node]]['atomtype'],list):
                                            r_copy.nodes[node]['atomtype'].append(d.nodes[mapping[node]]['atomtype'])
                                        else:
                                            r_copy.nodes[node]['atomtype'].extend(d.nodes[mapping[node]]['atomtype'])                                   
                                    else:
                                        if not isinstance(d.nodes[mapping[node]]['atomtype'],list):
                                            r_copy.nodes[node]['atomtype'] = [d.nodes[mapping[node]]['atomtype']]
                                        else:
                                            r_copy.nodes[node]['atomtype'] = d.nodes[mapping[node]]['atomtype']
                                        
                                        if not isinstance(d.nodes[mapping[node]]['map'],list):                              
                                            r_copy.nodes[node]['map'] = [d.nodes[mapping[node]]['map']]
                                        else:
                                            r_copy.nodes[node]['map'] = d.nodes[mapping[node]]['map']      

                                new_Otherside[r].append(r_copy)           
        
        Otherside = new_Otherside    
    return Otherside

def find_atommaps(graph):
    r_atommap_set = set()
    m_atommap_set = set()
    for node in graph.nodes():
        attrs = graph.nodes[node]
        if 'map' in attrs:
            if isinstance(attrs['map'], list) == False:
                attrs['map'] = [attrs['map']]
            if isinstance(attrs['atomtype'], list) == False:
                attrs['atomtype'] = [attrs['atomtype']]                
            for typ in range(len(attrs['atomtype'])):
                if attrs['atomtype'][typ] == 'r' and attrs['label'] != 'H':          
                    r_atommap_set.add(attrs['map'][typ])     
                if attrs['atomtype'][typ] == 'm' and attrs['label'] != 'H': 
                    m_atommap_set.add(attrs['map'][typ])
                      
    return sorted(r_atommap_set), sorted(m_atommap_set)

# translate the nx bond names to the MOD bond names 
def update_bound_attribute(G):
    for u, v in G.edges():
        data = G.edges[u,v]
        bound_value = data['label']
        if bound_value == 1.0:
            data['label'] = '-'
        elif bound_value == 2.0:
            data['label'] = '='
        elif bound_value == 1.5:
            data['label'] = ':'
        elif bound_value == 3:
            data['label'] = '#'                      
    return G
     
# Merge the RCLASS graphs with each other    
def add_Mappings(combined_list):  

        max_combinations = 1000   # test 10000 is to much
        combinations = list(itertools.islice(itertools.product(*combined_list),max_combinations))

        if len(combinations) == max_combinations:
            toBig = True 
        overlapping_graphs = []
        Side = []
        for comb in combinations:
            transposed = zip(*comb)
            overlapping_graphs.append(list(map(list, transposed)))

        # funktion checks if a map is already in the node
        def add_overlapping_maps(n1,n2):
            add_mapList = []
            add_atomtypeList = []
            for i in range(len(n2['map'])):
                check_map = True
                for j in range(len(n1['map'])):
                    if n2['map'][i] == n1['map'][j] and n2['atomtype'][i] == n1['atomtype'][j]:
                        check_map = False
                if check_map == True:
                    add_mapList.append(n2['map'][i])
                    add_atomtypeList.append(n2['atomtype'][i])
            return add_mapList,add_atomtypeList
            
        # add maps of the other RCLASSES to the Graphs          
        for i in overlapping_graphs:
            molecule_List = []
            for mol in i:
                if mol[0] == False:
                    continue
                G1 = copy.deepcopy(mol[0])       
                for g2 in range(1,len(mol)):
                    if mol[g2] == False:
                        continue
                    G2 = copy.deepcopy(mol[g2])
                    for n in G1.nodes.keys(): # add mapping from other RCLASS          
                        # find the nodes which were mapped in Overlapping Graph an add it to the fist graph
                        if 'atomtype' in G2.nodes[n]: 
                            if 'atomtype' in G1.nodes[n]:    
                                if isinstance(G1.nodes[n]['atomtype'],list):
                                    if isinstance(G2.nodes[n]['atomtype'],list): #both are lists
                                        add_mapList,add_atomtypeList = add_overlapping_maps(G1.nodes[n],G2.nodes[n])
                                        G1.nodes[n]['atomtype'].extend(add_atomtypeList)
                                        G1.nodes[n]['map'].extend(add_mapList)
                                    else: # only G1 is a List
                                        G2.nodes[n]['atomtype'] = [G2.nodes[n]['atomtype']]
                                        if isinstance(G2.nodes[n]['map'],list):
                                            G2.nodes[n]['map'] = [G2.nodes[n]['map']]
                                        add_mapList,add_atomtypeList = add_overlapping_maps(G1.nodes[n],G2.nodes[n])
                                        G1.nodes[n]['atomtype'].extend(add_atomtypeList)
                                        G1.nodes[n]['map'].extend(add_mapList)
                                else: # G1 is a single string entry
                                    G1.nodes[n]['atomtype'] = [G1.nodes[n]['atomtype']]
                                    if isinstance(G1.nodes[n]['map'],list):
                                        G1.nodes[n]['map'] = [G1.nodes[n]['map']]
                                    if isinstance(G2.nodes[n]['atomtype'],list):
                                        add_mapList,add_atomtypeList = add_overlapping_maps(G1.nodes[n],G2.nodes[n])
                                        G1.nodes[n]['atomtype'].extend(add_atomtypeList)
                                        G1.nodes[n]['map'].extend(add_mapList)
                                    else: 
                                        G2.nodes[n]['atomtype'] = [G2.nodes[n]['atomtype']]
                                        if isinstance(G2.nodes[n]['map'],list):
                                            G2.nodes[n]['map'] = [G2.nodes[n]['map']]
                                        add_mapList,add_atomtypeList = add_overlapping_maps(G1.nodes[n],G2.nodes[n])
                                        G1.nodes[n]['atomtype'].extend(add_atomtypeList)
                                        G1.nodes[n]['map'].extend(add_mapList)
                            else: # G1 has no atomtype
                                G1.nodes[n]['atomtype'] = G2.nodes[n]['atomtype']
                                G1.nodes[n]['comp_pair'] = G2.nodes[n]['comp_pair']
                                G1.nodes[n]['map'] = G2.nodes[n]['map']                                                       

                molecule_List.append(G1)           
            Side.append(molecule_List)             
        return Side      
      
def rclass_inAllAtomtype_Versions(g1,g2):
    # NOTE if g2 == None then the case occurs where several RCLASSes come on top of each other -> There only M has dominance, otherwise the overlaps are from one RCLASS itself, there D or M can dominate
    # If R-Atom ind atomtype then always selected as R-Atom 
    def domoinant_Atomtype(g,t):    
        for node in g.nodes():
            data = g.nodes[node]
            if 'atomtype' in data:
                if isinstance(data['atomtype'],list):
                    pos = [index for index, value in enumerate(data['atomtype']) if value == t]  
                    if len(pos) > 1 and len(pos) == len(data['atomtype']) or not pos:
                        continue
                    else:
                        new_a = []
                        new_m = []
                        for p in pos:                       
                            new_a.append(data['atomtype'][p])    
                            new_m.append(data['map'][p])   
                        data['atomtype'] = new_a  
                        data['map'] = new_m

    # Dominant R-Atom
    domoinant_Atomtype(g1,'r')
    if g2 != None:
        domoinant_Atomtype(g2,'r')
                    
    #check if double assignment of atom types M and D take place at all
    case_check = False
    for node in g1.nodes():
        data = g1.nodes[node]
        if 'atomtype' in data:
            if len(data['atomtype']) > 1 and 'm' in data['atomtype'] and 'd' in data['atomtype']:
                case_check = True
    if g2 != None:        
        for node in g2.nodes():
            data = g2.nodes[node]
            if len(data['atomtype']) > 1 and 'm' in data['atomtype'] and 'd' in data['atomtype']:
                case_check = True

    if case_check == True:           
        if g2 != None:      
            g1_M = copy.deepcopy(g1)
            g1_D = copy.deepcopy(g1)
            g2_M = copy.deepcopy(g2)    
            g2_D = copy.deepcopy(g2)    
        else:
            g1_M = copy.deepcopy(g1)
        # Dominant M-Atom
        domoinant_Atomtype(g1_M,'m')
        if g2 != None:
            # Domonant D-Atom alternativ
            domoinant_Atomtype(g1_D,'d')
            domoinant_Atomtype(g2_M,'m')
            domoinant_Atomtype(g2_D,'d')    
            return [(g1_M,g2_M),(g1_D,g2_D)]
        else:
            return [g1_M]
    else:
        if g2 != None:
            return [(g1,g2)]
        else:
             return [g1]        
             
             
def find_connected_nodes(graph, nodes):
    connected_nodes = set()
    atomlabel_list = []
    for u, v in graph.edges():
        if u in nodes:
            connected_nodes.add(u)
        if v in nodes:
            connected_nodes.add(v)            
    for i in list(connected_nodes):
        atomlabel_list.append(graph.nodes[i]['label'])
    return connected_nodes,atomlabel_list
             
             
#spezial case CO2: When C03 is eliminated, the D subgraph cannot be found, as the single bond from an oxygen to carbon becomes a double bond.       
#The function finds this special situation and maps the D atoms 
def c02_case(side1,side2,marker1,marker2,node_no,dict_node,Context_Text):       
    seen_dnodes = []
    if any(data.get('comp') == 'C00011' for _, data in side1.nodes(data=True)):
        # find the missing D-Atom Subgraph without hydrogens
        d_sub = create_subgraph_with_DAtoms(side2)
        d = d_sub.copy()
        hydrogens = []
        for n in d.nodes():
            a = d.nodes[n]
            if a['label'] == 'H':
                hydrogens.append(n)
        d.remove_nodes_from(hydrogens)

        co2_nodes,co2_atomlabels = find_connected_nodes(d, d.nodes())
        if sorted(co2_atomlabels) == ['C', 'O', 'O']:
            for node_r in side1.nodes():    
                if side1.nodes[node_r]['comp'] == 'C00011':
                    for node_l in co2_nodes:
                        if str(node_l+marker2) not in seen_dnodes and str(node_r+marker1) not in seen_dnodes and side1.nodes[node_r]['label'] == d.nodes[node_l]['label']:
                            dict_node[node_no] = [node_l,node_r] #New No. = [old No. L, old No. R]
                            Context_Text.append('  node [ id '+ str(node_no) + ' label "'+ d.nodes[node_l]['label'] +'" ]\n')
                            node_no = node_no + 1
                            seen_dnodes.append(str(node_l+marker2))
                            seen_dnodes.append(str(node_r+marker1))
                            break
    return node_no,dict_node,Context_Text

# correct the elemination of molecule parts into small molecules                          
def smallMolecule_Corr(graph,S,Context_Text,Left_Text,Right_Text,node_no,dict_node):
    small_result = True
    for node in graph.nodes:
        if graph.nodes[node].get('atomtype') == ['d'] and graph.nodes[node].get('label') != 'H':
            small_result = False
            neighbors = [n for n in graph.neighbors(node) if graph.nodes[n].get('label') != 'H']
            non_atomtype_neighbors = [n for n in neighbors if 'atomtype' not in graph.nodes[n]]
            atomtype_neighbors = [n for n in neighbors if 'atomtype' in graph.nodes[n]]
            # if there is a Neighbor with is part of the reaction rule than no filling of hydrogens is necessary
            if atomtype_neighbors:
                small_result = None
            if len(non_atomtype_neighbors) == 0:
                small_result = True
                continue
            # if there are only one non-hydrogen neighbor return the node end the edge for corraction     
            if len(non_atomtype_neighbors) == 1:
                neighbor = non_atomtype_neighbors[0]
                edge_attributes = graph.get_edge_data(node, neighbor)
                small_result = [neighbor, edge_attributes,node]  
                
                # check for dopple bond oxygens in the neighborhood of the non_atomtype_neighbors
                alcohol = False
                neighbors_of_non_atomtype_neighbors = [n for n in graph.neighbors(neighbor) if graph.nodes[n].get('label') == 'O' and  'atomtype' not in graph.nodes[n]]
                if len(neighbors_of_non_atomtype_neighbors) > 0:
                    for oxygen in neighbors_of_non_atomtype_neighbors:
                        alcohols = [n for n in graph.neighbors(oxygen) if graph.nodes[oxygen].get('label') == 'H']
                        alcohol = True
           
            # add bond changing and hydrogens to rule 
            if small_result != None and small_result != False:
                # add neighbor atom
                Context_Text.append('  node [ id '+ str(node_no) + ' label "'+ graph.nodes[small_result[0]]['label'] +'" ]\n')
                no_small = str(node_no)
                node_no = node_no + 1                                                         
                if small_result[1] == {'label': 1.0}:
                    bond_label = '-'
                    h = 1    
                if small_result[1] == {'label': 1.5}:
                    bond_label = ':'
                    h = 1                                            
                if small_result[1] == {'label': 2.0}:
                    bond_label = '='
                    h = 2
                if small_result[1] == {'label': 3.0}:
                    bond_label = '#'
                    h = 3                                                       
                # add edge     
                for nr in dict_node.keys():
                    if S =='L':
                        if small_result[2] == dict_node[nr][0]:
                            Left_Text.append('  edge [ source '+ str(nr) + ' target '+ no_small + ' label "'+ bond_label + '" ]\n') 
                    if S =='R':
                        if small_result[2] == dict_node[nr][1]:                                          
                            Right_Text.append('  edge [ source '+ str(nr) + ' target '+ no_small + ' label "'+ bond_label + '" ]\n')                 
                # add hydrogens                                                            
                while h > 0 and alcohol == False:
                    if S == 'L':
                        Left_Text.append('  node [ id '+ str(node_no) + ' label "H+" ]\n')
                        Right_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')                                                             
                        Right_Text.append('  edge [ source '+ no_small + ' target '+ str(node_no) + ' label "-" ]\n')
                    if S == 'R':
                        Right_Text.append('  node [ id '+ str(node_no) + ' label "H+" ]\n')
                        Left_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')                                                              
                        Left_Text.append('  edge [ source '+ no_small + ' target '+ str(node_no) + ' label "-" ]\n')
                    node_no = node_no + 1 
                    h = h-1

    if small_result == False:
            return False,False,False,False
    else:
            return Context_Text,Left_Text,Right_Text,node_no                         
    return Context_Text,Left_Text,Right_Text,node_no 


def findRCLASS_inMOL(compound,sides,toBig,S,reaction_with_mapping):
    if compound[0] == sides[0]:   
        mapping_left = subgraph_isomorphism_with_attributes(compound[1],r[0],node_match)    
        mapping_left = delete_sub_mappings(mapping_left)                  
        if len(mapping_left) > 100:
            toBig = True 
            return False,toBig    
        if len(mapping_left) == 0:
            return False, toBig
        if reaction_with_mapping != False:
            reaction_with_mapping[S].append({'comp_graph':compound[1],'rc_graph':r[0],'mapping':mapping_left}) 
    elif compound[0] == sides[1]:                        
        mapping_left = subgraph_isomorphism_with_attributes(compound[1],r[1],node_match)       
        mapping_left = delete_sub_mappings(mapping_left) 
        if len(mapping_left) > 100:
            toBig = True
            return False,toBig
        if len(mapping_left) == 0:
            return False,toBig  
        if reaction_with_mapping != False:                   
            reaction_with_mapping[S].append({'comp_graph':compound[1],'rc_graph':r[1],'mapping':mapping_left})  
    elif reaction_with_mapping != False:
        reaction_with_mapping[S].append({'comp_graph':[compound[1]],'rc_graph':None,'mapping':None})        
    
    return reaction_with_mapping, toBig
                 
# Stats
no_reactions = 0 #No. of TOTAL Reactions
no_reactions_false = 0 #No. of Reaction were no DPO-Rule could be created      
no_reactions_true = 0 #No. of Reaction were a DPO-Rule succsesful could be created         
no_reaction_toBig = 0 #No. of Reaction were the combinatorics became too big
missing_comp = 0 #No. of missing Compounds


if not os.path.exists('./03_stats'):
    os.makedirs('./03_stats')

log_unsucc = './03_stats/UnsuccessfulRXN.log'
log_big = './03_stats/ToBigReactions.log'
log_missC = './03_stats/faild_Compounds.txt'


# load reaction data
reactions = {}
input_path = './Additional_Files/REACTION_RCLASS_DATA.txt'
with open(input_path,'r') as in_f:
    lines = in_f.readlines()
rn = None
comp = None
rc = None   
for line in lines:
    line = line.split(':')
    if line[0] == 'Reaction ID':
        rn = line[1]
    if line[0] == 'Compound IDs':
        comp = line[1]
    if line[0] == 'RCLASS':
        rc = line[1]
    if rn != None and comp !=None and rc !=None:
        if rc != '' and comp != '':
            reactions[rn.rstrip()] =  {'compounds': comp, 'rclass': rc}
            rn = None
            comp = None
            rc = None  

# start progress   
for rn in reactions.keys():
        toBig = False #markes if some cases has to be scipted because of combinatoric explosion    
        no_reactions = no_reactions+1   
        print('START REACTION', rn)
        # load RCLASS 
        rxn_split = reactions[rn]['rclass'].split("'")
        rc_list = []
        rc_comp_list = []
        rc_graphs = {}
        for rc in rxn_split:
            if rc.startswith('RC'):
                rc_list.append(rc)
            elif rc.startswith('C'):
                rc_comp_list.append(rc)
        
        #check if there are RCLASSes listed
        if len(rc_list) == 0:
            exit()
        else:
            rclass_num = len(rc_list) # is used later to check the right number of R-Atoms
            
        for rc in range(len(rc_list)):
            filepath = input_path+rc_list[rc]
            gml_files = glob.glob(filepath+'/*.gml', recursive=True)
            if len(gml_files) == 0:
                break
            graphs_l = []
            graphs_r = []    
            for g in gml_files:      
                if 'left' in g:
                    graphs_l.append(nx.read_gml(g))
                if 'right' in g:
                    graphs_r.append(nx.read_gml(g))          
            rc_graphs[rc_list[rc]] = {rc_comp_list[rc]:{'left':graphs_l,'right':graphs_r}}     

        if len(gml_files) == 0:
            exit()
            
        # List of compounds which RCLASS are discribring
        rc_comp = []
        for c in rc_comp_list:
            c = c.split('_')
            rc_comp = rc_comp + c

        # build reaction
        reaction = {}
        reaction['Left'] = []
        reaction['Right'] = []
        right_side = False
        eq_list = reactions[rn]['compounds'].split("'")
        s = None
        for eq in eq_list:
            if eq.isdigit() and right_side == True:
                s = int(eq)
            elif eq.startswith('C') and right_side == True:
                reaction['Right'].append((s,eq))
                s = None
            elif eq.isdigit() and right_side == False:
                s = int(eq)
            elif eq.startswith('C') and right_side == False:
                reaction['Left'].append((s,eq))
                s = None
            if eq == '<=>':
                right_side = True

        compund_list = []
        for entry in reaction['Left']:
            try:
                m1 = get_compound_mol(entry[1])
                if m1 is None and entry[1] in rc_comp:
                    compund_list = []
                    break
                g1 = mol_to_graph(m1)                
                # add compound name to graph
                for n in g1.nodes():
                    g1.nodes[n]['comp'] = entry[1]
                
                count = entry[0]
                while count != 0:
                    compund_list.append((entry[1],g1))
                    if count == None:
                        count = 0
                    else:
                        count = count - 1 
            except:
                if entry[1] in rc_comp:
                    compund_list = []
                    break

        if len(compund_list) == 0:
            missing_comp = missing_comp+1    
            with open(log_missC,'a') as log:
                log.write(str(rn)+'\n')
            continue
                    
        reaction['Left'] = compund_list
        
        compund_list = []
        for entry in reaction['Right']:
            try:
                m1 = get_compound_mol(entry[1])
                if m1 is None and entry[1] in rc_comp:
                    compund_list = []
                    break
                g1 = mol_to_graph(m1) 
                # add compound name to graph
                for n in g1.nodes():
                    g1.nodes[n]['comp'] = entry[1]       
                count = entry[0]
                while count != 0:
                    compund_list.append((entry[1],g1))
                    if count == None:
                        count = 0
                    else:
                        count = count - 1 
            except:
                if entry[1] in rc_comp:
                    compund_list = []
                    break

        if len(compund_list) == 0:
            missing_comp = missing_comp+1    
            with open(log_missC,'a') as log:
                log.write(str(rn)+'\n')
            continue

        reaction['Right'] = compund_list

        # Combine all the possibilities of RCLASSes
        # A list is created which shows the reactions with the respective mapped RCLASS: [(Compound Pair of RCLASS]{'Left':[{'comp','comp_graph','rc_graph','mapping'},{...}],'Right':[{...},{...}]})]   
        inc = 100 # number unique map ids
        mapping_RCLASS = {}        
        for comp_pair in rc_graphs.keys():   
            mol_withMapping = []
            sides_list = rc_graphs[comp_pair].keys()
            for k in sides_list:
                sides = k.split('_')  
                rc_combis = list(product(rc_graphs[comp_pair][k]['left'],rc_graphs[comp_pair][k]['right']))

            # build all atomytpe versions 
            rc_combis_new = []
            for rc in rc_combis:
                rc = rclass_inAllAtomtype_Versions(rc[0],rc[1])             
                rc_combis_new.append(rc)
            rc_combis = rc_combis_new  
                        
            if len(rc_combis) > 300:
                toBig = True
                count_graph = len(rc_combis)
                break            
            # try every combination of RCLASS
            for rc in rc_combis:
                for r in rc:         
                    give_unique_map(r[0], inc)
                    give_unique_map(r[1], inc)
                    reaction_with_mapping = {}
                    reaction_with_mapping['Left'] = []
                    reaction_with_mapping['Right'] = []
                    switch = False # information if rclass compunds have do be switched
 
                    # left side matching of the linked compound 
                    for compound in reaction['Left']:   

                        reaction_with_mapping,toBig = findRCLASS_inMOL(compound,sides,toBig,'Left',reaction_with_mapping)
                        if not reaction_with_mapping: 
                            reaction_with_mapping = {}
                            reaction_with_mapping['Left'] = []
                            reaction_with_mapping['Right'] = []
                            sides[0], sides[1] = sides[1], sides[0]
                            switch = True
                            for compound in reaction['Left']: 
                                reaction_with_mapping,toBig = findRCLASS_inMOL(compound,sides,toBig,'Left',reaction_with_mapping)
                                if not reaction_with_mapping:        
                                    break       
                    if not reaction_with_mapping: 
                        break
                        
                    # right side matching of the linked compound               
                    for compound in reaction['Right']:                             
                        reaction_with_mapping,toBig = findRCLASS_inMOL(compound,sides,toBig,'Right',reaction_with_mapping)
                        if not reaction_with_mapping:  
                            break   
                    if not reaction_with_mapping and switch == True:
                        break   
                    if not reaction_with_mapping and switch == False:      
                        reaction_with_mapping = {}
                        reaction_with_mapping['Left'] = []
                        reaction_with_mapping['Right'] = []  
                        sides[0], sides[1] = sides[1], sides[0]             
                        for compound in reaction['Left']:     
                            reaction_with_mapping,toBig = findRCLASS_inMOL(compound,sides,toBig,'Left',reaction_with_mapping)
                            if not reaction_with_mapping:
                                break
                        for compound in reaction['Right']:  
                            reaction_with_mapping,toBig = findRCLASS_inMOL(compound,sides,toBig,'Right',reaction_with_mapping)    
                            if not reaction_with_mapping:
                                break                                                    
                    if not reaction_with_mapping:
                        break                     

                    # overwrite mapping information from RCLASS to reaction metabolite for each mapping alternative             
                    list_mapping_num = []
                    for l in reaction_with_mapping['Left']:
                        if l['mapping'] != None:                                           
                            mapped_comp_left,list_mapping_num = add_mapping(l['mapping'],l['comp_graph'],l['rc_graph'],list_mapping_num,sides)                                    
                            for gl in mapped_comp_left:
                                for n in gl.nodes():
                                    if 'comp_pair' not in gl.nodes[n] and 'atomtype' in gl.nodes[n]:
                                        gl.nodes[n]['comp_pair'] = sides                                     
                            l['comp_graph'] = mapped_comp_left                    
                        if l['comp_graph'] == False:
                            break     
                    for r in reaction_with_mapping['Right']:
                        if r['mapping'] != None:
                            mapped_comp_right,list_mapping_num = add_mapping(r['mapping'],r['comp_graph'],r['rc_graph'],list_mapping_num,sides)                          
                            for gl in mapped_comp_right:
                                for n in gl.nodes():
                                    if 'comp_pair' not in gl.nodes[n] and 'atomtype' in gl.nodes[n]:
                                        gl.nodes[n]['comp_pair'] = sides  

                            r['comp_graph'] = mapped_comp_right
                        if r['comp_graph'] == False:
                            break

                    if len(reaction_with_mapping['Left']) != 0:
                        if len(reaction_with_mapping['Right']) != 0: 
                            if inc not in  mapping_RCLASS.keys():
                                mapping_RCLASS[inc] = []
                            mapping_RCLASS[inc].append(reaction_with_mapping)
                inc = inc + 100     
         
        # Combine all compounds on each side for each possible mapping on the side for an RCLASS
        left_side = []
        right_side = []     
        for rc_var in mapping_RCLASS.keys(): 
            for r in mapping_RCLASS[rc_var]:
                l_side = []
                r_side = []   
                # filtere mapped compounds out 
                for dict_comp in r['Left']: 
                    if isinstance(dict_comp['comp_graph'],list):
                        l_side.append(dict_comp['comp_graph']) 
                    else: 
                        l_side.append([dict_comp['comp_graph']]) 
                for dict_comp in r['Right']: 
                    if isinstance(dict_comp['comp_graph'],list):
                        r_side.append(dict_comp['comp_graph']) 
                    else: 
                        r_side.append([dict_comp['comp_graph']]) 
                left_side.append(l_side)
                right_side.append(r_side)
                
            # create all Combinations of mapped molecules
            combined_list_right = []
            combined_list_left = []
            for sublist in left_side:
                        combined_sublist = [list(p) for p in product(*sublist)]
                        combined_list_left.append(combined_sublist)    
            for sublist in right_side:
                        combined_sublist = [list(p) for p in product(*sublist)]
                        combined_list_right.append(combined_sublist)

                            
        # all feasible mappings from both reaction sides                 
        Left = add_Mappings(combined_list_left)       
        Right = add_Mappings(combined_list_right)  

        #clear of dominant Atoms (R- over M- over D-Atoms)
        new = []
        for L in Left:
                l_new = []
                for l in L:
                    l_dom = rclass_inAllAtomtype_Versions(l,None)
                    l_new.extend(l_dom)
                new.append(l_new)
        Left = new
        new = []
        for L in Right:
                l_new = []
                for l in L:
                    l_dom = rclass_inAllAtomtype_Versions(l,None)
                    l_new.extend(l_dom)
                new.append(l_new)
        Right = new                               
        if Left == False or Right == False:
                continue                      
        # create Rules 
        no = 1
        succ = False   
        # filter to big Reactions
        if len(Left)>1000:
                count_graph = len(Left)
                toBig = True
                continue
        if len(Right)>1000:
                count_graph = len(Right)
                toBig = True
                continue
       
        for L in Left:                        
                for R in Right:               
                    # filter Cases I with to much Molecules
                    L_new = []
                    for i in L:
                        if len(i)<1000:
                            L_new.append(i)
                        else:
                            count_graph = len(i)
                            toBig = True
                    R_new = []
                    for i in R:
                        if len(i)<1000:
                            R_new.append(i)  
                        else:
                            count_graph = len(i)
                            toBig = True    
                    
                    # find for all combinations of D-Atoms matching
                    Right1 = find_Datoms(L_new,R_new)                    
                    Left1 = find_Datoms(R_new,L_new)

                    # filter Cases II with to much D-Atom mapped Subgraphs
                    def count_entries(data):
                        count = 0
                        for sublist in data:
                            count += len(sublist)
                        return count
                        
                    if count_entries(Left1) > 300:
                        count_graph = len(Left1)
                        toBig = True
                        continue                     
                    if count_entries(Right1) > 300:
                        count_graph = len(Right1)
                        toBig = True
                        continue  
                                        
                    Right1 = [list(p) for p in product(*Right1)]
                    Left1 = [list(p) for p in product(*Left1)]
                                                    
                    if len(Left1) > 300:
                        toBig = True
                        count_graph = len(Left1)
                        continue
                    if len(Right1) > 300:
                        toBig = True
                        count_graph = len(Right1)
                        continue
                           
                    #try all combinations for every D-Atom Matching
                    for l1 in Left1:                        
                        # compose all molecule at left side together
                        Left_new = []
                        graph_num = 0
                        for left in l1:
                            l = left.copy()
                            l = nx.relabel_nodes(l, {node: f"{graph_num}_{node}" for node in l.nodes()})
                            graph_num = graph_num + 1       
                            Left_new.append(l)
                        if len(Left_new) == 0:
                            #print('if len(Left_new) == 0:')
                            break

                        Left_new = nx.compose_all(Left_new)     
                        
                                           
                        # Resolve unwanted nesting of map lists 
                        for n,data in Left_new.nodes(data=True):
                            if 'map' in data:
                                Left_new.nodes[n]['map'] = flatten(Left_new.nodes[n]['map'])  
                       
                       
                        # write DPO-Rules for all possibile reaction atom to atom maps     
                        # create Lists of the Atommaps for all Atom Typs        
                        r_nodes_Left, m_nodes_Left = find_atommaps(Left_new)                                              
                        for r1 in Right1:   
                            if len(r1) != 0:
                                # compose all molecule at right side together
                                Right_new = []
                                graph_num = 0       
                                for right in r1:
                                    r = right.copy()
                                    r = nx.relabel_nodes(r, {node: f"{graph_num}_{node}" for node in r.nodes()})
                                    graph_num = graph_num + 1
                                    Right_new.append(r)
                                Right_new = nx.compose_all(Right_new) 
                                if len(Right_new) == 0:
                                    continue
                           
                                # Resolve unwanted nesting of map lists 
                                for n,data in  Right_new.nodes(data=True):
                                    if 'map' in data:
                                        Right_new.nodes[n]['map'] = flatten(Right_new.nodes[n]['map'])   
                             
                                # write DPO-Rules for all possibile reaction atom to atom maps     
                                # create Lists of the Atommaps for all Atom Typs  
                                r_nodes_Right, m_nodes_Right = find_atommaps(Right_new) 
                                if len(r_nodes_Left) ==  0 or len(r_nodes_Right) == 0:
                                    continue
                          
                                if len(r_nodes_Left) < rclass_num or len(r_nodes_Right) < rclass_num :   
                                    continue
                          
                                if len(r_nodes_Left) != len(r_nodes_Right) or len(m_nodes_Left) != len(m_nodes_Right):                                          
                                    continue  
                               
                                # genarate the DPO for all possibile atom-to-atom maps
                                # Initiallise lists with the entries to be written for each part of the DPO rule
                                Left_Text = []
                                Context_Text = []
                                Right_Text = []  
                                node_no = 1 # new node number system for DPO rule 
                                dict_node = {} # dictonary for rewrite edges with new number system
                                                                                       
                                #check for right Number of R-/D- and M-Atoms
                                left_R_atoms = []
                                left_D_atoms = []
                                left_M_atoms = []
                                for n in Left_new.nodes():
                                    data = Left_new.nodes[n]
                                    if 'atomtype' in data:
                                        if data['label'] != 'H':
                                            if 'r' in data['atomtype']:
                                                left_R_atoms.append(n)
                                            elif 'd' in data['atomtype']:
                                                left_D_atoms.append(n)
                                            elif 'm' in data['atomtype']:
                                                left_M_atoms.append(n)
                                right_R_atoms = []
                                right_D_atoms = []
                                right_M_atoms = []
                                for n in Right_new.nodes():
                                    data = Right_new.nodes[n]
                                    if 'atomtype' in data:
                                        if data['label'] != 'H':
                                            if 'r' in data['atomtype']:
                                                right_R_atoms.append(n)
                                            elif 'd' in data['atomtype']:
                                                right_D_atoms.append(n)
                                            elif 'm' in data['atomtype']:
                                                right_M_atoms.append(n)          

                                if len(left_R_atoms) != len(right_R_atoms) or len(left_M_atoms) != len(right_M_atoms):                                                         
                                    continue           
                                                                                                                      
                                # write R-Atoms      
                                seen_nodes = []                
                                found_maps = []    
                                hydrogen_left = []
                                hydrogen_right = []

                                for m in r_nodes_Left:       
                                    m_o = m.split('_')
                                    mp1 = int(m_o[0])+1
                                    mm1 = int(m_o[0])-1
                                    for node_l in Left_new.nodes(): 
                                        attrs_l = Left_new.nodes[node_l]
                                        if 'map' not in attrs_l:
                                            continue      
                                        if m in attrs_l['map'] and attrs_l['label'] != 'H':
                                            atom_l = attrs_l['label']
                                            node_r = None
                                            for r in Right_new.nodes():
                                                attrs_r = Right_new.nodes[r]
                                                if 'map' in attrs_r:
                                                    map_list = []
                                                    for m in attrs_r['map']:
                                                        m_o = m.split('_')
                                                        mp2 = int(m_o[0])
                                                        map_list.append(mp2)

                                                    if (mp1 in map_list or mm1 in map_list) and attrs_r['label'] == atom_l and r not in seen_nodes and 'r' in attrs_r['atomtype']:
                                                        node_r = r
                                                        seen_nodes.append(node_r)
                                                        found_maps.append(m)
                                                        break
                                                                        
                                            if [node_l,node_r] not in dict_node.values():
                                                Context_Text.append('  node [ id '+ str(node_no) + ' label "'+ atom_l +'" ]\n')
                                                dict_node[node_no] = [node_l,node_r] #New No. = [old No. L, old No. R]
                                                node_no = node_no + 1

                                # check if all R_Atoms were found on both sides
                                if len(found_maps) != len(r_nodes_Left):
                                   continue                                
                                # write R-Hydrogen                                
                                for node in Left_new.nodes():
                                    if 'map' in Left_new.nodes[node]:
                                        if Left_new.nodes[node]['label'] == 'H' and 'r' in Left_new.nodes[node]['atomtype']:
                                            hydrogen_left.append(node)
                                            for r in range(len(Left_new.nodes[node]['atomtype'])):
                                                if Left_new.nodes[node]['atomtype'][r] == 'r':
                                                    if not isinstance(Left_new.nodes[node]['map'],list) and r == 0:
                                                        m = Left_new.nodes[node]['map'] 
                                                    else:
                                                        m = Left_new.nodes[node]['map'][r]
                                            m_o = m.split('_')
                                            mp1 = int(m_o[0])+1
                                            mp1 = str(mp1)+'_'+m_o[1]
                                            mm1 = int(m_o[0])-1
                                            mm1 = str(mm1)+'_'+m_o[1]
                                            found = False
                                            for node in Right_new.nodes():
                                                if 'map' in Right_new.nodes[node]:
                                                    if (mp1 in Right_new.nodes[node]['map'] or mm1 in Right_new.nodes[node]['map']) and 'r' in Right_new.nodes[node]['atomtype'] and Right_new.nodes[node]['label'] == 'H' and node not in seen_nodes:
                                                        found = True
                                                        seen_nodes.append(node)
                                                        hydrogen_right.append(node)   
                                                        break
                                            if found == False:
                                                hydrogen_right.append(None)    
                                         
                                # find R-Hydrogens which are only one the right side
                                for node in Right_new.nodes():
                                    if 'map' in Right_new.nodes[node] and node not in seen_nodes:
                                        if Right_new.nodes[node]['label'] == 'H' and 'r' in Right_new.nodes[node]['atomtype']:
                                            seen_nodes.append(node)
                                            hydrogen_left.append(None)  
                                            hydrogen_right.append(node)          
                             
                                # if no hydrogens got lost than write then in the contect 
                                # if Right Graph has more hydrogens then write this one in Right and as H+ in Left          
                                for i in range(len(hydrogen_left)):
                                    if [hydrogen_left[i],hydrogen_right[i]] not in dict_node.values() and hydrogen_left[i] != None and hydrogen_right[i] != None:
                                        Context_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')
                                        dict_node[node_no] = [hydrogen_left[i],hydrogen_right[i]] #New No. = [old No. L, old No. R]
                                        node_no = node_no + 1    
                                    elif [hydrogen_left[i],hydrogen_right[i]] not in dict_node.values() and hydrogen_left[i] == None and hydrogen_right[i] != None:
                                        Left_Text.append('  node [ id '+ str(node_no) + ' label "H+" ]\n')
                                        Right_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')
                                        dict_node[node_no] = [None,hydrogen_right[i]] #New No. = [old No. L, old No. R]
                                        node_no = node_no + 1                     
                                    elif [hydrogen_left[i],hydrogen_right[i]] not in dict_node.values() and hydrogen_left[i] != None and hydrogen_right[i] == None:
                                        Right_Text.append('  node [ id '+ str(node_no) + ' label "H+" ]\n')
                                        Left_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')
                                        dict_node[node_no] = [hydrogen_left[i],None] #New No. = [old No. L, old No. R]
                                        node_no = node_no + 1         
         
                                # write D-Atoms        
                                D_check = True # turn in to False if D-Atom has no Patner, and mark the side where the Partner is missing   
                                for node_l in Left_new.nodes(): 
                                    attrs = Left_new.nodes[node_l]
                                    if 'map' not in attrs:
                                        continue      
                                    if 'd' in attrs['atomtype'] and attrs['label'] != 'H' and D_check == True:        
                                        D_check = 'Right'     
                                        atom_l = attrs['label']
                                        map_l = attrs['map']
                                        for node_r in Right_new.nodes():   
                                            attrs = Right_new.nodes[node_r]
                                            if 'map' not in attrs:
                                                continue       
                                            if attrs['map'] == map_l and atom_l == attrs['label'] and node_r not in seen_nodes:
                                                D_check = True
                                                dict_node[node_no] = [node_l,node_r] #New No. = [old No. L, old No. R]
                                                Context_Text.append('  node [ id '+ str(node_no) + ' label "'+ atom_l +'" ]\n')
                                                node_no = node_no + 1
                                                seen_nodes.append(node_r)
                                                break

                                # check D-Atom Partners also for the Right Side
                                for node_r in Right_new.nodes(): 
                                    attrs = Right_new.nodes[node_r]
                                    if 'map' not in attrs:
                                        continue      
                                    if 'd' in attrs['atomtype'] and attrs['label'] != 'H' and D_check == True:        
                                        if D_check == True:
                                            D_check = 'Left'     
                                        else: 
                                            D_check = 'Both' 
                                        atom_r = attrs['label']
                                        map_r = attrs['map']
                                        for node_l in Left_new.nodes():   
                                            attrs = Left_new.nodes[node_l]
                                            if 'map' not in attrs:
                                                continue     
                                            if attrs['map'] == map_r and atom_r == attrs['label'] and node_l not in seen_nodes:
                                                D_check = True           

                                # if one D Atom has no partner -> check for spezial CO2 Case
                                if D_check != True:
                                    if D_check == 'Left':
                                        node_no,dict_node,Context_Text = c02_case(Left_new,Right_new,'L','R',node_no,dict_node,Context_Text)
                                    if D_check == 'Right':
                                        node_no,dict_node,Context_Text = c02_case(Right_new,Left_new,'R','L',node_no,dict_node,Context_Text)      
                                        #print('New Context',Context_Text)           
                                    if D_check == 'Both':
                                        node_no,dict_node,Context_Text = c02_case(Left_new,Right_new,'L','R',node_no,dict_node,Context_Text)
                                        node_no,dict_node,Context_Text = c02_case(Right_new,Left_new,'R','L',node_no,dict_node,Context_Text) 
                                    # check again if D-Atoms have missing partners
                                    for n in Left_new.nodes():                                  
                                        if 'atomtype' in Left_new.nodes[n]:
                                            if 'd' in Left_new.nodes[n]['atomtype'] and Left_new.nodes[n]['label'] != 'H':
                                                D_check = False 
                                                for i in dict_node.values():
                                                    if n in i: 
                                                        D_check = True
                                                if D_check == False:
                                                    break

                                    for n in Right_new.nodes():
                                        if 'atomtype' in Right_new.nodes[n]:
                                            if 'd' in Right_new.nodes[n]['atomtype'] and Right_new.nodes[n]['label'] != 'H':
                                                D_check = False 
                                                for i in dict_node.values():
                                                    if n in i: 
                                                        D_check = True
                                                if D_check == False:
                                                    break  
                                    # if still D_check Fails -> -> Atom to Atom Map not complete! Try Next Combi
                                    if D_check == False:
                                        continue
                                                            
                                # correced of small molecules
                                Context_Text,Left_Text,Right_Text,node_no = smallMolecule_Corr(Left_new,'L',Context_Text,Left_Text,Right_Text,node_no,dict_node)
                                #if the mapped the D-Atom to a bigger structure the mapping is probably wrong
                                if Context_Text == False:
                                    continue 
                                Context_Text,Left_Text,Right_Text,node_no = smallMolecule_Corr(Right_new,'R',Context_Text,Left_Text,Right_Text,node_no,dict_node)
                                if Context_Text == False:
                                    continue 
                                    
                                # write M-Atoms
                                #  inicialize control structures
                                m_nodes_no = 0
                                for i in Left_new.nodes():
                                    if 'map' in Left_new.nodes[i]:
                                        if Left_new.nodes[i]['label'] != 'H' and 'd' not in Left_new.nodes[i]['atomtype'] and 'r' not in Left_new.nodes[i]['atomtype']:
                                            m_nodes_no = m_nodes_no + 1 
                                            
                                seen_nodes = []                        
                                # find the matching non hydrogen M-Atoms 
                                for m in range(len(m_nodes_Left)): 
                                    # find node with Atommap ID
                                    for node_l in Left_new.nodes():      
                                        attrs1 = Left_new.nodes[node_l]  
                                        if 'map' in attrs1:                                
                                            if m_nodes_Left[m] in attrs1['map'] and 'r' not in attrs1['atomtype'] and attrs1['label'] != 'H':
                                                atom_l = attrs1['label']  
                                                # find corresponding M-Atom on the Right Side
                                                for node in Right_new.nodes():
                                                    attrs2 = Right_new.nodes[node]
                                                    if 'map' in attrs2:
                                                        if m_nodes_Right[m] in attrs2['map'] and attrs2['label'] == atom_l and node not in seen_nodes:
                                                            node_r = None
                                                            if len(attrs2['map']) != len(attrs1['map']):
                                                                continue
                                                            if len(attrs2['map']) == 1:
                                                                node_r = node 
                                                                seen_nodes.append(node)                                                  
                                                            else:
                                                                # check if all maps fit together
                                                                fit_check = True
                                                                for i in range(len(attrs2['map'])):
                                                                    if attrs2['map'][i] in m_nodes_Right and attrs1['map'][i] in m_nodes_Left:
                                                                        if m_nodes_Right.index(attrs2['map'][i]) != m_nodes_Left.index(attrs1['map'][i]):
                                                                            fit_check = False
                                                                            break
                                                                    else: 
                                                                        fit_check = False
                                                                        break
                                                                if fit_check == True:   
                                                                    # check if the neighborhood is the same
                                                                    m_neighbors_left = [n for n in Left_new.neighbors(node_l) if Left_new.nodes[n].get('label') != 'H']  
                                                                    list_neighbor_atoms_left = [Left_new.nodes[n].get('label') for n in m_neighbors_left] 
                                                                    m_neighbors_right = [n for n in Right_new.neighbors(node) if Right_new.nodes[n].get('label') != 'H']  
                                                                    list_neighbor_atoms_right = [Right_new.nodes[n].get('label') for n in m_neighbors_right] 
                                                                    if sorted(list_neighbor_atoms_left) == sorted(list_neighbor_atoms_right):
                                                                        node_r = node 
                                                                        seen_nodes.append(node)  
                                                            if node_r != None:
                                                                if [node_l,node_r] not in dict_node.values():
                                                                    Context_Text.append('  node [ id '+ str(node_no) + ' label "'+ atom_l +'" ]\n')
                                                                    dict_node[node_no] = [node_l,node_r] #New No. = [old No. L, old No. R]
                                                                    node_no = node_no + 1       
                                                                    break        
  
                                # check if all M-Atoms were found on both sides
                                if len(seen_nodes) != m_nodes_no: 
                                    continue
                                                                         
                                # write M-Hydrogens 
                                hydrogen_left = []
                                hydrogen_right = []
                                for node in Left_new.nodes():
                                    attrs = Left_new.nodes[node]
                                    if 'map' in attrs:
                                        if attrs['label'] == 'H' and 'r' not in attrs['atomtype'] and 'd' not in attrs['atomtype']:
                                            hydrogen_left.append(node)
                                            # find the corresponding atom map ID for the M-hydrogen
                                            if not isinstance(attrs['atomtype'],list):
                                                attrs['atomtype'] = [attrs['atomtype']]
                                            for r in range(len(attrs['atomtype'])):
                                                if attrs['atomtype'][r] == 'm':
                                                    if isinstance(attrs['map'],list):
                                                        m = attrs['map'][r]
                                                    else:
                                                        m = attrs['map']
                                                for n in range(len(m_nodes_Left)):
                                                    if m_nodes_Left[n] == m:
                                                        mp1 = m_nodes_Right[n]        

                                            found = False
                                            for node in Right_new.nodes():
                                                attrs = Right_new.nodes[node]
                                                if 'map' in attrs:
                                                    if mp1 in attrs['map'] and 'r' not in attrs['atomtype'] and 'd' not in attrs['atomtype'] and attrs['label'] == 'H' and node not in seen_nodes:
                                                        found = True
                                                        seen_nodes.append(node)
                                                        hydrogen_right.append(node)   
                                                        break
                                            if found == False:
                                                hydrogen_right.append(None)   
                                                        
                                for node in Right_new.nodes():
                                    attrs = Right_new.nodes[node]
                                    if 'map' in attrs and node not in seen_nodes:
                                        if attrs['label'] == 'H' and 'r' not in attrs['atomtype'] and 'd' not in attrs['atomtype']:
                                            seen_nodes.append(node)
                                            hydrogen_left.append(None)  
                                            hydrogen_right.append(node)                                       

                                # if no hydrogens got lost than write then in the contect 
                                # if Right Graph has more hydrogens then write this one in Right and as H+ in Left          
                                for i in range(len(hydrogen_left)):
                                    if [hydrogen_left[i],hydrogen_right[i]] not in dict_node.values() and hydrogen_left[i] != None and hydrogen_right[i] != None:
                                        Context_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')
                                        dict_node[node_no] = [hydrogen_left[i],hydrogen_right[i]] #New No. = [old No. L, old No. R]
                                        node_no = node_no + 1    
                                    if [hydrogen_left[i],hydrogen_right[i]] not in dict_node.values() and hydrogen_left[i] == None and hydrogen_right[i] != None:
                                        Left_Text.append('  node [ id '+ str(node_no) + ' label "H+" ]\n')
                                        Right_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')
                                        dict_node[node_no] = [None,hydrogen_right[i]] #New No. = [old No. L, old No. R]
                                        node_no = node_no + 1                     
                                    if [hydrogen_left[i],hydrogen_right[i]] not in dict_node.values() and hydrogen_left[i] != None and hydrogen_right[i] == None:
                                        Right_Text.append('  node [ id '+ str(node_no) + ' label "H+" ]\n')
                                        Left_Text.append('  node [ id '+ str(node_no) + ' label "H" ]\n')
                                        dict_node[node_no] = [hydrogen_left[i],None] #New No. = [old No. L, old No. R]
                                        node_no = node_no + 1    
                                                                   
                                # translate numeric edge labelling in symbolic labelling 
                                L1 = update_bound_attribute(copy.deepcopy(Left_new))                                
                                R1 = update_bound_attribute(copy.deepcopy(Right_new))                   
                                                                 
                                # write Edges
                                trans_Ledges = {}
                                trans_Redges = {}  
                                ori_Ledges = list(L1.edges(data=True))
                                ori_Redges = list(R1.edges(data=True))

                                # translete Edges with new node numer system
                                for edge in L1.edges():
                                    node1 = None
                                    node2 = None
                                    for i in dict_node:
                                        if edge[0] == dict_node[i][0]:
                                            node1 = i
                                        if edge[1] == dict_node[i][0]:
                                            node2 = i
                                    if node1 == None or node2 == None:
                                        L1.remove_edge(edge[0],edge[1])
                                    else:
                                        trans_Ledges[(str(node1),str(node2))] = L1.edges[edge]['label']

                                for edge in R1.edges():
                                    node1 = None
                                    node2 = None
                                    for i in dict_node:
                                        if edge[0] == dict_node[i][1]:
                                            node1 = i
                                        if edge[1] == dict_node[i][1]:
                                            node2 = i
                                    if node1 == None or node2 == None:
                                        R1.remove_edge(edge[0],edge[1])
                                        
                                    
                                    else:
                                        trans_Redges[(str(node1),str(node2))] = R1.edges[edge]['label']          
                                
                                # check for common edges
                                seen_edges = []
                                for e1,e2 in trans_Ledges.keys():
                                    # first node order e1,e2
                                    if (e1,e2) in trans_Redges.keys() and (e1,e2) not in seen_edges:
                                        seen_edges.append((e1,e2))
                                        if trans_Ledges[(e1,e2)] == trans_Redges[(e1,e2)]:
                                            Context_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                        else:
                                            Left_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                            Right_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Redges[(e1,e2)]+ '" ]\n')
                                    # second node order e2,e1
                                    elif (e2,e1) in trans_Redges.keys() and (e2,e1) not in seen_edges:
                                        seen_edges.append((e2,e1))
                                        if trans_Ledges[(e1,e2)] == trans_Redges[(e2,e1)]:
                                            Context_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                        else:
                                            Left_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                            Right_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Redges[(e2,e1)]+ '" ]\n')
                               
                                # check for edges only of one side        
                                for e1,e2 in trans_Ledges.keys():             
                                    if (e1,e2) not in seen_edges and (e2,e1) not in seen_edges:
                                        Left_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')                 
                                for e1,e2 in trans_Redges.keys():                
                                    if (e1,e2) not in seen_edges and (e2,e1) not in seen_edges:
                                        Right_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Redges[(e1,e2)] + '" ]\n')         
                                              
                                # check if are found changes on left or right side
                                if len(Left_Text) == 0 and len(Right_Text) == 0: 
                                    continue         
                                
                                # write gml file                       
                                if not os.path.exists('./03_DPO_Rules/'+rn):
                                    os.makedirs('./03_DPO_Rules/'+rn) 
                                with open('./03_DPO_Rules/'+rn+'/'+rn+'_'+str(no)+'.gml', 'w') as out: 
                                    out.write('rule [\n')
                                    out.write(' ruleID "' + rn + '"\n')
                                    out.write(' left [\n') 
                                    for i in Left_Text:
                                        out.write(i)
                                    out.write(' ]\n')
                                    out.write(' context [\n')    
                                    for i in Context_Text:
                                        out.write(i)                      
                                    out.write(' ]\n')      
                                    out.write(' right [\n')        
                                    for i in Right_Text:
                                        out.write(i)        
                                    out.write(' ]\n')    
                                    out.write(']\n')          
                                succ = True
                                no = no + 1
        
        if succ == False:
            if toBig == True:
                no_reaction_toBig = no_reaction_toBig +1
                with open(log_big,'a') as log:
                    log.write(str(rn)+' No. of graphs:'+str(count_graph)+'\n')
            else:
                no_reactions_false = no_reactions_false+1    
                with open(log_unsucc,'a') as log:
                    log.write(str(rn)+'\n')

        else:
            no_reactions_true = no_reactions_true+1        
 
print('############# Statistics #############')
print('No. of TOTAL Rections',no_reactions)
print('No. of created DPO-Rules',no_reactions_true)
print('No. of Rections were no DPO-Rule could be created',no_reactions_false)
print('No. of Reactions which to big to create', no_reaction_toBig)    
print('No. of Reaction with missing Compounds',missing_comp)
