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
import sys
<<<<<<< HEAD
from tqdm import tqdm
=======
>>>>>>> dbf44ddee9f253e25b1bdc61ff6cc036e1c13eb1


# shows the progress over the algorithm
def progress_bar(progress, total, length=40):
    percent = int(progress / total * 100)
    filled = int(length * progress // total)
    bar = "#" * filled + "-" * (length - filled)
    sys.stdout.write(f"\r[{bar}] {percent}%")
    sys.stdout.flush()


# This function takes a nested list and returns a flat list.
def flatten(nested_list):

    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))
        else:
            flat_list.append(item)

    return flat_list


# Load compound mol format form KEGG Database
def get_compound_mol(compound_id):
    url = f"http://rest.kegg.jp/get/compound:{compound_id}/mol"
    try:
        response = requests.get(url)
        response.raise_for_status()  # This throws an exception if the status code is not 200 (OK)
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Error in the request: {e}")
        return None


# Create the compound from KEGG mol data
def mol_to_graph(mol_str):
    mol = Chem.MolFromMolBlock(mol_str)
    Chem.rdmolops.AssignStereochemistry(mol, force=True, cleanIt=True)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    Chem.SanitizeMol(
        mol, Chem.SanitizeFlags.SANITIZE_ALL - Chem.SanitizeFlags.SANITIZE_PROPERTIES
    )

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
    if b["atom"] == "R" and a["label"] != "H":
        return True
    else:
        return a["label"] == b["atom"]


def node_matchD(a, b):
    if b["label"] == "H":
        return True
    return a["label"] == b["label"]


# Find all matches of the RCLASS in the molecule
def subgraph_isomorphism_with_attributes(M, R, node_match_opt):

    def edge_match(a, b):
        return str(a["label"]) == str(b["label"])

    GM = iso.GraphMatcher(M, R, node_match=node_match_opt, edge_match=edge_match)
    subgraph_matches = list(GM.subgraph_monomorphisms_iter())
    return subgraph_matches


# Check for permutaion of matching nodes (most because of hydrogens)
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


# Adds the mapping from the RCLASS graphs to the reaction
def add_mapping(mappings, comp, rc, list_mapping_num, sides):
    # Copy maps and atom types from RCLASS to compound graph
    mapped_comp = [comp]
    mapping_no = 0
    right_side = False
    for mapping in mappings:
        L = copy.deepcopy(comp)

        # If it is the right side, take the mappings from the left side
        if len(list_mapping_num) == 0:
            right_side = True
        for n in mapping.keys():
            if not isinstance(rc.nodes[mapping[n]]["atomtype"], list):
                rc.nodes[mapping[n]]["atomtype"] = [rc.nodes[mapping[n]]["atomtype"]]
            L.nodes[n]["atomtype"] = rc.nodes[mapping[n]]["atomtype"]
            L.nodes[n]["comp_pair"] = sides
            L_maps = []
            if right_side == True:
                for i in rc.nodes[mapping[n]]["map"]:
                    L_maps.append(str(i) + "_" + str(mapping_no))
                    if mapping_no not in list_mapping_num:
                        list_mapping_num.append(mapping_no)
            else:
                if len(list_mapping_num) != 0:
                    mapping_no = list_mapping_num.pop(0)
                    right_side = True
                for i in rc.nodes[mapping[n]]["map"]:
                    L_maps.append(str(i) + "_" + str(mapping_no))
            L.nodes[n]["map"] = L_maps
        mapping_no = mapping_no + 1
        if L not in mapped_comp:
            mapped_comp.append(L)

    # Combine double mapping on one molecule
    mapping_no_copy = mapping_no
    if len(mapped_comp) > 2:
        combinations = []
        for r in range(1, len(mapped_comp) + 1):
            max_combinations = 1000
            combinations.extend(
                itertools.islice(
                    itertools.combinations(mapped_comp, r), max_combinations
                )
            )
        result = [list(comb) for comb in combinations]
        mapped_comp = []
        for comp in result:
            mapping_no = mapping_no_copy
            if len(comp) > 1:
                H = copy.deepcopy(comp[0])
                for g in range(1, len(comp)):
                    G = copy.deepcopy(comp[g])
                    for n in G.nodes():
                        data = copy.deepcopy(H.nodes[n])
                        if "atomtype" in G.nodes[n]:
                            if "atomtype" in H.nodes[n]:
                                if not isinstance(H.nodes[n]["atomtype"], list):
                                    H.nodes[n]["atomtype"] = [H.nodes[n]["atomtype"]]
                                if isinstance(G.nodes[n]["atomtype"], list):
                                    H.nodes[n]["atomtype"].extend(
                                        G.nodes[n]["atomtype"]
                                    )
                                else:
                                    H.nodes[n]["atomtype"].append(
                                        G.nodes[n]["atomtype"]
                                    )

                                # Add new mapping_no to distinguish the same RCLASSes
                                if not isinstance(H.nodes[n]["map"], list):
                                    H.nodes[n]["map"] = [H.nodes[n]["map"]]
                                map_att_new = []
                                for i in G.nodes[n]["map"]:
                                    new_num = i.split("_")
                                    map_att_new.append(
                                        new_num[0] + "_" + str(mapping_no)
                                    )
                                H.nodes[n]["map"].extend(map_att_new)
                            else:
                                if isinstance(G.nodes[n]["atomtype"], list):
                                    H.nodes[n]["atomtype"] = G.nodes[n]["atomtype"]
                                else:
                                    H.nodes[n]["atomtype"] = [G.nodes[n]["atomtype"]]
                                if isinstance(G.nodes[n]["map"], list):
                                    H.nodes[n]["map"] = G.nodes[n]["map"]
                                else:
                                    H.nodes[n]["map"] = [G.nodes[n]["map"]]

                    mapping_no = mapping_no + 1
                mapped_comp.append(H)
            else:
                mapped_comp.append(comp[0])

    return mapped_comp, list_mapping_num


def create_subgraph_with_DAtoms(graph):
    subgraph_nodes = []
    for node in graph.nodes():
        data = graph.nodes[node]
        if "atomtype" in data:
            if data["atomtype"] is not None:
                if "d" in data["atomtype"]:
                    subgraph_nodes.append(node)
    subgraph = graph.subgraph(subgraph_nodes)
    subgraph_edges = [
        (u, v)
        for u, v in subgraph.edges()
        if u in subgraph_nodes and v in subgraph_nodes
    ]

    return subgraph.subgraph(subgraph_nodes)


def give_unique_map(graph, increment):
    for node in graph.nodes():
        data = graph.nodes[node]
        if "map" in data:
            if isinstance(data["map"], list):
                for i in range(len(data["map"])):
                    data["map"][i] = data["map"][i] + increment
            else:
                data["map"] += increment


# Tautomerism function which in graphene transforms -OH to =O, and vice versa
def tautomerism(graph):
    nodes_to_remove = []
    check_OH = False
    # -OH to =O
    g = copy.deepcopy(graph)
    for node in graph.nodes():
        if graph.nodes[node].get("label") == "O":
            neighbors = list(graph.neighbors(node))
            for neighbor in neighbors:
                if graph.nodes[neighbor].get("label") == "H":
                    check_OH = True
                    nodes_to_remove.append(neighbor)
                    g.remove_node(neighbor)
                    remaining_neighbors = [
                        n for n in g.neighbors(node) if n != neighbor
                    ]
                    for remaining_neighbor in remaining_neighbors:
                        g.edges[node, remaining_neighbor]["label"] = "="
    g.remove_nodes_from(nodes_to_remove)
    if not check_OH:
        # =O to -OH
        for node in graph.nodes():
            if graph.nodes[node].get("label") == "O":
                neighbors = list(graph.neighbors(node))
                for neighbor in neighbors:
                    if graph.edges[node, neighbor].get("label") == "=":
                        g.edges[node, neighbor]["label"] = "-"
                        new_node_id = max(g.nodes) + 1
                        new_node_attributes = g.nodes[node].copy()
                        new_node_attributes["label"] = "H"
                        g.add_node(new_node_id, **new_node_attributes)
                        g.add_edge(node, new_node_id, label="-")
                        new_nodes.append(new_node_id)
    return g


def find_Datoms(Side, Otherside):
    new_Otherside = []
    for i in Otherside:
        graph_copy = copy.deepcopy(i)
        new_Otherside.append([graph_copy])
    Otherside = new_Otherside
    for l in Side:
        # Create subgraphs with D-atoms which are not covered by other RCLASSes
        # It lists also different versions: Dominant D (all nodes with D), restricted D (all nodes only with D)
        d_list = []
        d_sub = create_subgraph_with_DAtoms(l)
        d_list.append(d_sub.copy())
        d_taut = tautomerism(d_sub.copy())
        d_list.append(d_taut)

        for d in d_list:
            # Remove the Hydrogen atoms
            hydrogens = []
            for n in d.nodes():
                a = d.nodes[n]
                if a["label"] == "H":
                    hydrogens.append(n)
            d.remove_nodes_from(hydrogens)

            # Define forbidden molecules
            forbidden_mol = []
            for n in d.nodes():
                for m in d.nodes[n]["comp_pair"]:
                    if m not in forbidden_mol:
                        forbidden_mol.append(m)

            # Find D subgraphs
            if len(d.nodes()) > 0:
                for r in range(len(Otherside)):
                    for graph in range(len(Otherside[r])):

                        # Find mappings of D_sub on Otherside
                        mapping_datom = subgraph_isomorphism_with_attributes(
                            Otherside[r][graph], d, node_matchD
                        )
                        mapping_datom_mod = delete_sub_mappings(mapping_datom)

                        if len(mapping_datom_mod) > 1000:
                            continue
                        for mapping in mapping_datom_mod:
                            # Check if the matching are R-/M-atoms or only D-atoms or no atom types at all
                            RM_true = False
                            D_true = False
                            noAT_true = False
                            AT_check = True
                            for node in mapping.keys():
                                if "atomtype" in Otherside[r][graph].nodes[node]:
                                    if (
                                        Otherside[r][graph].nodes[node]["comp"]
                                        in forbidden_mol
                                    ):
                                        AT_check = False
                                        break
                                    if noAT_true == True:
                                        AT_check = False
                                        break
                                    if (
                                        Otherside[r][graph].nodes[node]["atomtype"]
                                        is not None
                                    ):
                                        if (
                                            "d"
                                            in Otherside[r][graph].nodes[node][
                                                "atomtype"
                                            ]
                                        ):
                                            if RM_true == False:
                                                D_true = True
                                            else:
                                                AT_check = False
                                                break
                                    if (
                                        Otherside[r][graph].nodes[node]["atomtype"]
                                        == "r"
                                        or Otherside[r][graph].nodes[node]["atomtype"]
                                        == "m"
                                    ):
                                        if D_true == False:
                                            RM_true = True
                                        else:
                                            AT_check = False
                                            break
                                else:
                                    noAT_true == True
                            if AT_check == False:
                                continue

                            # Add the found D-subgraphs in the moleculess
                            r_copy = copy.deepcopy(Otherside[r][graph])
                            for node in mapping.keys():
                                if "map" in r_copy.nodes[node]:
                                    if (
                                        d.nodes[mapping[node]]["map"]
                                        == r_copy.nodes[node]["map"]
                                    ):
                                        continue
                                    if not isinstance(
                                        d.nodes[mapping[node]]["map"], list
                                    ):
                                        r_copy.nodes[node]["map"].append(
                                            d.nodes[mapping[node]]["map"]
                                        )
                                    else:
                                        r_copy.nodes[node]["map"].extend(
                                            d.nodes[mapping[node]]["map"]
                                        )

                                    if not isinstance(
                                        d.nodes[mapping[node]]["atomtype"], list
                                    ):
                                        r_copy.nodes[node]["atomtype"].append(
                                            d.nodes[mapping[node]]["atomtype"]
                                        )
                                    else:
                                        r_copy.nodes[node]["atomtype"].extend(
                                            d.nodes[mapping[node]]["atomtype"]
                                        )
                                else:
                                    if not isinstance(
                                        d.nodes[mapping[node]]["atomtype"], list
                                    ):
                                        r_copy.nodes[node]["atomtype"] = [
                                            d.nodes[mapping[node]]["atomtype"]
                                        ]
                                    else:
                                        r_copy.nodes[node]["atomtype"] = d.nodes[
                                            mapping[node]
                                        ]["atomtype"]

                                    if not isinstance(
                                        d.nodes[mapping[node]]["map"], list
                                    ):
                                        r_copy.nodes[node]["map"] = [
                                            d.nodes[mapping[node]]["map"]
                                        ]
                                    else:
                                        r_copy.nodes[node]["map"] = d.nodes[
                                            mapping[node]
                                        ]["map"]

                            new_Otherside[r].append(r_copy)

        Otherside = new_Otherside
    return Otherside


def find_atommaps(graph):
    r_atommap_set = set()
    m_atommap_set = set()
    for node in graph.nodes():
        attrs = graph.nodes[node]
        if "map" in attrs:
            if isinstance(attrs["map"], list) == False:
                attrs["map"] = [attrs["map"]]
            if isinstance(attrs["atomtype"], list) == False:
                attrs["atomtype"] = [attrs["atomtype"]]
            for typ in range(len(attrs["atomtype"])):
                if attrs["atomtype"][typ] == "r" and attrs["label"] != "H":
                    r_atommap_set.add(attrs["map"][typ])
                if attrs["atomtype"][typ] == "m" and attrs["label"] != "H":
                    m_atommap_set.add(attrs["map"][typ])

    return sorted(r_atommap_set), sorted(m_atommap_set)


# Translate the nx bond names to the MOD bond names
def update_bound_attribute(G):
    for u, v in G.edges():
        data = G.edges[u, v]
        bound_value = data["label"]
        if bound_value == 1.0:
            data["label"] = "-"
        elif bound_value == 2.0:
            data["label"] = "="
        elif bound_value == 1.5:
            data["label"] = ":"
        elif bound_value == 3:
            data["label"] = "#"
    return G


# Merge the RCLASS graphs with each other
def add_Mappings(combined_list):

    max_combinations = 10000
    combinations = list(
        itertools.islice(itertools.product(*combined_list), max_combinations)
    )

    if len(combinations) == max_combinations:
        toBig = True
    overlapping_graphs = []
    Side = []
    for comb in combinations:
        transposed = zip(*comb)
        overlapping_graphs.append(list(map(list, transposed)))

    # Function checks if a map is already in the node
    def add_overlapping_maps(n1, n2):
        add_mapList = []
        add_atomtypeList = []
        for i in range(len(n2["map"])):
            check_map = True
            for j in range(len(n1["map"])):
                if (
                    n2["map"][i] == n1["map"][j]
                    and n2["atomtype"][i] == n1["atomtype"][j]
                ):
                    check_map = False
            if check_map == True:
                add_mapList.append(n2["map"][i])
                add_atomtypeList.append(n2["atomtype"][i])
        return add_mapList, add_atomtypeList

    # Add maps of the other RCLASSES to the Graphs
    print("\nStart Map Adding")
    total = len(overlapping_graphs)
    status = 0
    for i in overlapping_graphs:
        status += 1
        progress_bar(status, total)
        molecule_List = []
        for mol in i:
            if mol[0] == False:
                continue
            G1 = copy.deepcopy(mol[0])
            for g2 in range(1, len(mol)):
                if mol[g2] == False:
                    continue
                G2 = copy.deepcopy(mol[g2])
                for n in G1.nodes.keys():  # Add mapping from other RCLASS
                    # Find the nodes which were mapped in Overlapping Graph an add it to the fist graph
                    if "atomtype" in G2.nodes[n]:
                        if "atomtype" in G1.nodes[n]:
                            if isinstance(G1.nodes[n]["atomtype"], list):
                                if isinstance(
                                    G2.nodes[n]["atomtype"], list
                                ):  # Both are lists
                                    add_mapList, add_atomtypeList = (
                                        add_overlapping_maps(G1.nodes[n], G2.nodes[n])
                                    )
                                    G1.nodes[n]["atomtype"].extend(add_atomtypeList)
                                    G1.nodes[n]["map"].extend(add_mapList)
                                else:  # only G1 is a list
                                    G2.nodes[n]["atomtype"] = [G2.nodes[n]["atomtype"]]
                                    if isinstance(G2.nodes[n]["map"], list):
                                        G2.nodes[n]["map"] = [G2.nodes[n]["map"]]
                                    add_mapList, add_atomtypeList = (
                                        add_overlapping_maps(G1.nodes[n], G2.nodes[n])
                                    )
                                    G1.nodes[n]["atomtype"].extend(add_atomtypeList)
                                    G1.nodes[n]["map"].extend(add_mapList)
                            else:  # G1 is a single string entry
                                G1.nodes[n]["atomtype"] = [G1.nodes[n]["atomtype"]]
                                if isinstance(G1.nodes[n]["map"], list):
                                    G1.nodes[n]["map"] = [G1.nodes[n]["map"]]
                                if isinstance(G2.nodes[n]["atomtype"], list):
                                    add_mapList, add_atomtypeList = (
                                        add_overlapping_maps(G1.nodes[n], G2.nodes[n])
                                    )
                                    G1.nodes[n]["atomtype"].extend(add_atomtypeList)
                                    G1.nodes[n]["map"].extend(add_mapList)
                                else:
                                    G2.nodes[n]["atomtype"] = [G2.nodes[n]["atomtype"]]
                                    if isinstance(G2.nodes[n]["map"], list):
                                        G2.nodes[n]["map"] = [G2.nodes[n]["map"]]
                                    add_mapList, add_atomtypeList = (
                                        add_overlapping_maps(G1.nodes[n], G2.nodes[n])
                                    )
                                    G1.nodes[n]["atomtype"].extend(add_atomtypeList)
                                    G1.nodes[n]["map"].extend(add_mapList)
                        else:  # G1 has no atom type
                            G1.nodes[n]["atomtype"] = G2.nodes[n]["atomtype"]
                            G1.nodes[n]["comp_pair"] = G2.nodes[n]["comp_pair"]
                            G1.nodes[n]["map"] = G2.nodes[n]["map"]

            molecule_List.append(G1)
        Side.append(molecule_List)
    return Side


def rclass_inAllAtomtype_Versions(g1, g2):
    """NOTE if g2 == None then the case occurs where several RCLASSes come on top of each other -> There only M has dominance,
    otherwise the overlaps are from one RCLASS itself, there D or M can dominate"""

    # If R-Atom is in atom type then always selected as R-Atom
    def domoinant_Atomtype(g, t):
        for node in g.nodes():
            data = g.nodes[node]
            if "atomtype" in data:
                if isinstance(data["atomtype"], list):
                    pos = [
                        index
                        for index, value in enumerate(data["atomtype"])
                        if value == t
                    ]
                    if len(pos) > 1 and len(pos) == len(data["atomtype"]) or not pos:
                        continue
                    else:
                        new_a = []
                        new_m = []
                        for p in pos:
                            new_a.append(data["atomtype"][p])
                            new_m.append(data["map"][p])
                        data["atomtype"] = new_a
                        data["map"] = new_m

    # Dominant R-atom
    domoinant_Atomtype(g1, "r")
    if g2 != None:
        domoinant_Atomtype(g2, "r")

    # Check if double assignments of atom types M and D take place at all
    case_check = False
    for node in g1.nodes():
        data = g1.nodes[node]
        if "atomtype" in data:
            if (
                len(data["atomtype"]) > 1
                and "m" in data["atomtype"]
                and "d" in data["atomtype"]
            ):
                case_check = True
    if g2 != None:
        for node in g2.nodes():
            data = g2.nodes[node]
            if (
                len(data["atomtype"]) > 1
                and "m" in data["atomtype"]
                and "d" in data["atomtype"]
            ):
                case_check = True

    if case_check == True:
        if g2 != None:
            g1_M = copy.deepcopy(g1)
            g1_D = copy.deepcopy(g1)
            g2_M = copy.deepcopy(g2)
            g2_D = copy.deepcopy(g2)
        else:
            g1_M = copy.deepcopy(g1)
        # Dominant M-atom
        domoinant_Atomtype(g1_M, "m")
        if g2 != None:
            # Dominant D-atom alternative
            domoinant_Atomtype(g1_D, "d")
            domoinant_Atomtype(g2_M, "m")
            domoinant_Atomtype(g2_D, "d")
            return [(g1_M, g2_M), (g1_D, g2_D)]
        else:
            return [g1_M]
    else:
        if g2 != None:
            return [(g1, g2)]
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
        atomlabel_list.append(graph.nodes[i]["label"])
    return connected_nodes, atomlabel_list


"""Spezial case CO2: When CO3 is eliminated, the D subgraph cannot be found, 
as the single bond from an oxygen to carbon becomes a double bond."""


# The function finds this special situation and maps the D-atoms
def c02_case(side1, side2, marker1, marker2, node_no, dict_node, Context_Text):
    seen_dnodes = []
    if any(data.get("comp") == "C00011" for _, data in side1.nodes(data=True)):
        # Find the missing D-atom subgraph without hydrogens
        d_sub = create_subgraph_with_DAtoms(side2)
        d = d_sub.copy()
        hydrogens = []
        for n in d.nodes():
            a = d.nodes[n]
            if a["label"] == "H":
                hydrogens.append(n)
        d.remove_nodes_from(hydrogens)

        co2_nodes, co2_atomlabels = find_connected_nodes(d, d.nodes())
        if sorted(co2_atomlabels) == ["C", "O", "O"]:
            for node_r in side1.nodes():
                if side1.nodes[node_r]["comp"] == "C00011":
                    for node_l in co2_nodes:
                        if (
                            str(node_l + marker2) not in seen_dnodes
                            and str(node_r + marker1) not in seen_dnodes
                            and side1.nodes[node_r]["label"] == d.nodes[node_l]["label"]
                        ):
                            dict_node[node_no] = [node_l, node_r]
                            Context_Text.append(
                                "  node [ id "
                                + str(node_no)
                                + ' label "'
                                + d.nodes[node_l]["label"]
                                + '" ]\n'
                            )
                            node_no = node_no + 1
                            seen_dnodes.append(str(node_l + marker2))
                            seen_dnodes.append(str(node_r + marker1))
                            break
    return node_no, dict_node, Context_Text


# Correct the elimination of molecule parts into small molecules
def smallMolecule_Corr(
    graph, S, Context_Text, Left_Text, Right_Text, node_no, dict_node
):
    small_result = True
    for node in graph.nodes:
        if (
            graph.nodes[node].get("atomtype") == ["d"]
            and graph.nodes[node].get("label") != "H"
        ):
            small_result = False
            neighbors = [
                n for n in graph.neighbors(node) if graph.nodes[n].get("label") != "H"
            ]
            non_atomtype_neighbors = [
                n for n in neighbors if "atomtype" not in graph.nodes[n]
            ]
            atomtype_neighbors = [n for n in neighbors if "atomtype" in graph.nodes[n]]
            # If there is a neighbor which is part of the reaction rule then no filling of hydrogens is necessary
            if atomtype_neighbors:
                small_result = None
            if len(non_atomtype_neighbors) == 0:
                small_result = True
                continue
            # If there is only one non-hydrogen neighbor then return the node and the edge for correction
            if len(non_atomtype_neighbors) == 1:
                neighbor = non_atomtype_neighbors[0]
                edge_attributes = graph.get_edge_data(node, neighbor)
                small_result = [neighbor, edge_attributes, node]

                # Check for double bond oxygens in the neighborhood of the non_atomtype_neighbors
                alcohol = False
                neighbors_of_non_atomtype_neighbors = [
                    n
                    for n in graph.neighbors(neighbor)
                    if graph.nodes[n].get("label") == "O"
                    and "atomtype" not in graph.nodes[n]
                ]
                if len(neighbors_of_non_atomtype_neighbors) > 0:
                    for oxygen in neighbors_of_non_atomtype_neighbors:
                        alcohols = [
                            n
                            for n in graph.neighbors(oxygen)
                            if graph.nodes[oxygen].get("label") == "H"
                        ]
                        alcohol = True

            # Add bond changing and hydrogens to rule
            if small_result != None and small_result != False:
                # Add neighbor atom
                Context_Text.append(
                    "  node [ id "
                    + str(node_no)
                    + ' label "'
                    + graph.nodes[small_result[0]]["label"]
                    + '" ]\n'
                )
                no_small = str(node_no)
                node_no = node_no + 1
                if small_result[1] == {"label": 1.0}:
                    bond_label = "-"
                    h = 1
                if small_result[1] == {"label": 1.5}:
                    bond_label = ":"
                    h = 1
                if small_result[1] == {"label": 2.0}:
                    bond_label = "="
                    h = 2
                if small_result[1] == {"label": 3.0}:
                    bond_label = "#"
                    h = 3
                # Add edge
                for nr in dict_node.keys():
                    if S == "L":
                        if small_result[2] == dict_node[nr][0]:
                            Left_Text.append(
                                "  edge [ source "
                                + str(nr)
                                + " target "
                                + no_small
                                + ' label "'
                                + bond_label
                                + '" ]\n'
                            )
                    if S == "R":
                        if small_result[2] == dict_node[nr][1]:
                            Right_Text.append(
                                "  edge [ source "
                                + str(nr)
                                + " target "
                                + no_small
                                + ' label "'
                                + bond_label
                                + '" ]\n'
                            )
                # Add hydrogens
                while h > 0 and alcohol == False:
                    if S == "L":
                        Left_Text.append(
                            "  node [ id " + str(node_no) + ' label "H+" ]\n'
                        )
                        Right_Text.append(
                            "  node [ id " + str(node_no) + ' label "H" ]\n'
                        )
                        Right_Text.append(
                            "  edge [ source "
                            + no_small
                            + " target "
                            + str(node_no)
                            + ' label "-" ]\n'
                        )
                    if S == "R":
                        Right_Text.append(
                            "  node [ id " + str(node_no) + ' label "H+" ]\n'
                        )
                        Left_Text.append(
                            "  node [ id " + str(node_no) + ' label "H" ]\n'
                        )
                        Left_Text.append(
                            "  edge [ source "
                            + no_small
                            + " target "
                            + str(node_no)
                            + ' label "-" ]\n'
                        )
                    node_no = node_no + 1
                    h = h - 1

    if small_result == False:
        return False, False, False, False
    else:
        return Context_Text, Left_Text, Right_Text, node_no
    return Context_Text, Left_Text, Right_Text, node_no


def findRCLASS_inMOL(compound, sides, toBig, S, reaction_with_mapping):
    if compound[0] == sides[0]:
        mapping_left = subgraph_isomorphism_with_attributes(
            compound[1], r[0], node_match
        )
        mapping_left = delete_sub_mappings(mapping_left)
        if len(mapping_left) > 100:
            toBig = True
            return False, toBig
        if len(mapping_left) == 0:
            return False, toBig
        if reaction_with_mapping != False:
            reaction_with_mapping[S].append(
                {"comp_graph": compound[1], "rc_graph": r[0], "mapping": mapping_left}
            )
    elif compound[0] == sides[1]:
        mapping_left = subgraph_isomorphism_with_attributes(
            compound[1], r[1], node_match
        )
        mapping_left = delete_sub_mappings(mapping_left)
        if len(mapping_left) > 100:
            toBig = True
            return False, toBig
        if len(mapping_left) == 0:
            return False, toBig
        if reaction_with_mapping != False:
            reaction_with_mapping[S].append(
                {"comp_graph": compound[1], "rc_graph": r[1], "mapping": mapping_left}
            )
    elif reaction_with_mapping != False:
        reaction_with_mapping[S].append(
            {"comp_graph": [compound[1]], "rc_graph": None, "mapping": None}
        )

    return reaction_with_mapping, toBig


# Statistics
no_reactions = 0  # Number of TOTAL Reactions
no_reactions_false = 0  # Number. of reactions were no DPO-Rule could be created
no_reactions_true = 0  # Number of reactions were a DPO-Rule successful could be created
no_reaction_toBig = 0  # Number of reactions were the combinatorics became too big
missing_comp = 0  # Number of missing compounds

input_path = sys.argv[1]
<<<<<<< HEAD
input_rxn = sys.argv[2]
output_path = os.path.dirname(input_path)

if not os.path.exists(output_path + "/03_DPO_Rules/"):
    os.makedirs(output_path + "/03_DPO_Rules/")

log_unsucc = output_path + "UnsuccessfulRXN.log"
log_big = output_path + "ToBigReactions.log"
log_missC = output_path + "faild_Compounds.txt"

output_path = output_path + "/03_DPO_Rules/"

# load reactions
rxn_list = []
with open(input_rxn, "r") as f:
    rxns = f.readlines()
for rxn in rxns:
    rxn = rxn.strip().split("\t")
    rxn_list.append(rxn[0])

# Load reaction data
reactions = {}
path = "./Additional_Files/REACTION_RCLASS_DATA.txt"
with open(path, "r") as in_f:
=======

output_path = os.path.dirname(input_path)

if not os.path.exists(output_path+'/03_DPO_Rules/'):
    os.makedirs(output_path+'/03_DPO_Rules/')
if not os.path.exists(output_path+'/03_stats'):
    os.makedirs(output_path+'/03_stats')    

log_unsucc = output_path+'UnsuccessfulRXN.log'
log_big = output_path+'ToBigReactions.log'
log_missC = output_path+'faild_Compounds.txt'

output_path = output_path+'/03_DPO_Rules/'

# Load reaction data
reactions = {}
path = './Additional_Files/REACTION_RCLASS_DATA.txt'
with open(path,'r') as in_f:
>>>>>>> dbf44ddee9f253e25b1bdc61ff6cc036e1c13eb1
    lines = in_f.readlines()
rn = None
comp = None
rc = None
for line in lines:
    line = line.split(":")
    if line[0] == "Reaction ID":
        rn = line[1]
    if line[0] == "Compound IDs":
        comp = line[1]
    if line[0] == "RCLASS":
        rc = line[1]
    if rn != None and comp != None and rc != None:
        if rc != "" and comp != "":
            reactions[rn.rstrip()] = {"compounds": comp, "rclass": rc}
            rn = None
            comp = None
            rc = None

<<<<<<< HEAD
# Start progress
for rn in rxn_list:
    toBig = False  # Marks if some cases have to be scripted because of combinatoric explosion
    no_reactions = no_reactions + 1
    print("\nSTART REACTION", rn)
    # Load RCLASS
    rxn_split = reactions[rn]["rclass"].split("'")
    rc_list = []
    rc_comp_list = []
    rc_graphs = {}
    for rc in rxn_split:
        if rc.startswith("RC"):
            rc_list.append(rc)
        elif rc.startswith("C"):
            rc_comp_list.append(rc)
=======
# Start progress  
for rn in reactions.keys():
        toBig = False # Marks if some cases have to be scripted because of combinatoric explosion   
        no_reactions = no_reactions+1   
        print('START REACTION', rn)
        # Load RCLASS  
        rxn_split = reactions[rn]['rclass'].split("'")
        rc_list = []
        rc_comp_list = []
        rc_graphs = {}
        for rc in rxn_split:
            if rc.startswith('RC'):
                rc_list.append(rc)
            elif rc.startswith('C'):
                rc_comp_list.append(rc)
        
        #Check if there are RCLASSes listed
        if len(rc_list) == 0:
            exit()
        else:
            rclass_num = len(rc_list) # Is used later to check the right number of R-atoms
            
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
>>>>>>> dbf44ddee9f253e25b1bdc61ff6cc036e1c13eb1

    # Check if there are RCLASSes listed
    if len(rc_list) == 0:
        print("No RCLASSE  found")
        exit()
    else:
        rclass_num = len(rc_list)  # Is used later to check the right number of R-atoms

    for rc in range(len(rc_list)):
        filepath = input_path + "/" + rc_list[rc]
        gml_files = glob.glob(filepath + "/*.gml", recursive=True)
        if len(gml_files) == 0:
            break
        graphs_l = []
        graphs_r = []
        for g in gml_files:
            if "left" in g:
                graphs_l.append(nx.read_gml(g))
            if "right" in g:
                graphs_r.append(nx.read_gml(g))
        rc_graphs[rc_list[rc]] = {
            rc_comp_list[rc]: {"left": graphs_l, "right": graphs_r}
        }

    if len(gml_files) == 0:
        print("No GML found")
        exit()

    # List of compounds which RCLASS are describing
    rc_comp = []
    for c in rc_comp_list:
        c = c.split("_")
        rc_comp = rc_comp + c

    # Build reaction
    reaction = {}
    reaction["Left"] = []
    reaction["Right"] = []
    right_side = False
    eq_list = reactions[rn]["compounds"].split("'")
    s = None
    for eq in eq_list:
        if eq.isdigit() and right_side == True:
            s = int(eq)
        elif eq.startswith("C") and right_side == True:
            reaction["Right"].append((s, eq))
            s = None
        elif eq.isdigit() and right_side == False:
            s = int(eq)
        elif eq.startswith("C") and right_side == False:
            reaction["Left"].append((s, eq))
            s = None
        if eq == "<=>":
            right_side = True

    compund_list = []
    for entry in reaction["Left"]:
        try:
            m1 = get_compound_mol(entry[1])
            if m1 is None and entry[1] in rc_comp:
                compund_list = []
                break
            g1 = mol_to_graph(m1)
            # Add compound name to graph
            for n in g1.nodes():
                g1.nodes[n]["comp"] = entry[1]

            count = entry[0]
            while count != 0:
                compund_list.append((entry[1], g1))
                if count == None:
                    count = 0
                else:
                    count = count - 1
        except:
            if entry[1] in rc_comp:
                compund_list = []
                break

    if len(compund_list) == 0:
        missing_comp = missing_comp + 1
        with open(log_missC, "a") as log:
            log.write(str(rn) + "\n")
        continue

    reaction["Left"] = compund_list

    compund_list = []
    for entry in reaction["Right"]:
        try:
            m1 = get_compound_mol(entry[1])
            if m1 is None and entry[1] in rc_comp:
                compund_list = []
                break
            g1 = mol_to_graph(m1)
            # Add compound name to graph
            for n in g1.nodes():
                g1.nodes[n]["comp"] = entry[1]
            count = entry[0]
            while count != 0:
                compund_list.append((entry[1], g1))
                if count == None:
                    count = 0
                else:
                    count = count - 1
        except:
            if entry[1] in rc_comp:
                compund_list = []
                break

    if len(compund_list) == 0:
        missing_comp = missing_comp + 1
        with open(log_missC, "a") as log:
            log.write(str(rn) + "\n")
        continue

    reaction["Right"] = compund_list

    # Combine all the possibilities of RCLASSes
    """ A list is created which shows the reactions with the respective mapped RCLASS: 
        [(Compound Pair of RCLASS]{'Left':[{'comp','comp_graph','rc_graph','mapping'},{...}],'Right':[{...},{...}]})] """
    inc = 100  # Number unique map IDs
    mapping_RCLASS = {}
    for comp_pair in rc_graphs.keys():
        mol_withMapping = []
        sides_list = rc_graphs[comp_pair].keys()
        for k in sides_list:
            sides = k.split("_")
            rc_combis = list(
                product(
                    rc_graphs[comp_pair][k]["left"], rc_graphs[comp_pair][k]["right"]
                )
            )

        # Build all atom type versions
        rc_combis_new = []
        for rc in rc_combis:

            rc = rclass_inAllAtomtype_Versions(rc[0], rc[1])
            rc_combis_new.append(rc)
        rc_combis = rc_combis_new

        if len(rc_combis) > 10000:
            toBig = True
            count_graph = len(rc_combis)
            break
        # Try every combination of RCLASS
        for rc in rc_combis:
            for r in rc:
                give_unique_map(r[0], inc)
                give_unique_map(r[1], inc)
                reaction_with_mapping = {}
                reaction_with_mapping["Left"] = []
                reaction_with_mapping["Right"] = []
                switch = False  # Information if RCLASS compounds have to be switched

                # Left side matching of the linked compound
                for compound in reaction["Left"]:

                    reaction_with_mapping, toBig = findRCLASS_inMOL(
                        compound, sides, toBig, "Left", reaction_with_mapping
                    )
                    if not reaction_with_mapping:
                        reaction_with_mapping = {}
                        reaction_with_mapping["Left"] = []
                        reaction_with_mapping["Right"] = []
                        sides[0], sides[1] = sides[1], sides[0]
                        switch = True
                        for compound in reaction["Left"]:
                            reaction_with_mapping, toBig = findRCLASS_inMOL(
                                compound, sides, toBig, "Left", reaction_with_mapping
                            )
                            if not reaction_with_mapping:
                                break
                if not reaction_with_mapping:
                    break

                #  Right side matching of the linked compound
                for compound in reaction["Right"]:
                    reaction_with_mapping, toBig = findRCLASS_inMOL(
                        compound, sides, toBig, "Right", reaction_with_mapping
                    )
                    if not reaction_with_mapping:
                        break
                if not reaction_with_mapping and switch == True:
                    break
                if not reaction_with_mapping and switch == False:
                    reaction_with_mapping = {}
                    reaction_with_mapping["Left"] = []
                    reaction_with_mapping["Right"] = []
                    sides[0], sides[1] = sides[1], sides[0]
                    for compound in reaction["Left"]:
                        reaction_with_mapping, toBig = findRCLASS_inMOL(
                            compound, sides, toBig, "Left", reaction_with_mapping
                        )
                        if not reaction_with_mapping:
                            break
                    for compound in reaction["Right"]:
                        reaction_with_mapping, toBig = findRCLASS_inMOL(
                            compound, sides, toBig, "Right", reaction_with_mapping
                        )
                        if not reaction_with_mapping:
                            break
                if not reaction_with_mapping:
                    break

                # Overwrite mapping information from RCLASS to reaction metabolite for each mapping alternative
                list_mapping_num = []
                for l in reaction_with_mapping["Left"]:
                    if l["mapping"] != None:
                        mapped_comp_left, list_mapping_num = add_mapping(
                            l["mapping"],
                            l["comp_graph"],
                            l["rc_graph"],
                            list_mapping_num,
                            sides,
                        )
                        for gl in mapped_comp_left:
                            for n in gl.nodes():
                                if (
                                    "comp_pair" not in gl.nodes[n]
                                    and "atomtype" in gl.nodes[n]
                                ):
                                    gl.nodes[n]["comp_pair"] = sides
                        l["comp_graph"] = mapped_comp_left
                    if l["comp_graph"] == False:
                        break
                for r in reaction_with_mapping["Right"]:
                    if r["mapping"] != None:
                        mapped_comp_right, list_mapping_num = add_mapping(
                            r["mapping"],
                            r["comp_graph"],
                            r["rc_graph"],
                            list_mapping_num,
                            sides,
                        )
                        for gl in mapped_comp_right:
                            for n in gl.nodes():
                                if (
                                    "comp_pair" not in gl.nodes[n]
                                    and "atomtype" in gl.nodes[n]
                                ):
                                    gl.nodes[n]["comp_pair"] = sides

                        r["comp_graph"] = mapped_comp_right
                    if r["comp_graph"] == False:
                        break

                if len(reaction_with_mapping["Left"]) != 0:
                    if len(reaction_with_mapping["Right"]) != 0:
                        if inc not in mapping_RCLASS.keys():
                            mapping_RCLASS[inc] = []
                        mapping_RCLASS[inc].append(reaction_with_mapping)
            inc = inc + 100

    if toBig == True and len(mapping_RCLASS.keys()) == 0:
        print("tobig")
        break
    # Combine all compounds on each side for each possible mapping on the side for an RCLASS
    left_side = []
    right_side = []
    for rc_var in mapping_RCLASS.keys():
        for r in mapping_RCLASS[rc_var]:
            l_side = []
            r_side = []
            # Filter out mapped compoundst
            for dict_comp in r["Left"]:
                if isinstance(dict_comp["comp_graph"], list):
                    l_side.append(dict_comp["comp_graph"])
                else:
                    l_side.append([dict_comp["comp_graph"]])
            for dict_comp in r["Right"]:
                if isinstance(dict_comp["comp_graph"], list):
                    r_side.append(dict_comp["comp_graph"])
                else:
                    r_side.append([dict_comp["comp_graph"]])
            left_side.append(l_side)
            right_side.append(r_side)

        # Create all combinations of mapped molecules
        combined_list_right = []
        combined_list_left = []
        for sublist in left_side:
            combined_sublist = [list(p) for p in product(*sublist)]
            combined_list_left.append(combined_sublist)
        for sublist in right_side:
            combined_sublist = [list(p) for p in product(*sublist)]
            combined_list_right.append(combined_sublist)

    # All feasible mappings from both reaction sides
    Left = add_Mappings(combined_list_left)
    Right = add_Mappings(combined_list_right)
    # Remove of dominant Atoms (R- over M- over D-Atoms)
    new = []
    for L in Left:
        l_new = []
        for l in L:
            l_dom = rclass_inAllAtomtype_Versions(l, None)
            l_new.extend(l_dom)
        new.append(l_new)
    Left = new
    new = []
    for L in Right:
        l_new = []
        for l in L:
            l_dom = rclass_inAllAtomtype_Versions(l, None)
            l_new.extend(l_dom)
        new.append(l_new)
    Right = new
    if Left == False or Right == False:
        continue
    # Create rules
    no = 1
    succ = False
    # Filter too big reactions
    if len(Left) > 10000:
        count_graph = len(Left)
        toBig = True
        continue
    if len(Right) > 10000:
        count_graph = len(Right)
        toBig = True
        continue
    total = len(Left)  # to show progress
    status = 0
    print("\nStart Combining")
    for L in Left:
        status += 1
        progress_bar(status, total)
        for R in Right:
            # Filter cases I with too many molecules
            L_new = []
            for i in L:
                if len(i) < 10000:
                    L_new.append(i)
                else:
                    count_graph = len(i)
                    toBig = True
            R_new = []
            for i in R:
                if len(i) < 10000:
                    R_new.append(i)
                else:
                    count_graph = len(i)
                    toBig = True
            # Find all combinations of D-atoms matching
            Right1 = find_Datoms(L_new, R_new)
            Left1 = find_Datoms(R_new, L_new)

            # Filter cases II with too many D-atom mapped subgraphs
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
            # Try all combinations for every D-atom matching
            for l1 in Left1:
                # Compose all molecules at left side together
                Left_new = []
                graph_num = 0
                for left in l1:
                    l = left.copy()
                    l = nx.relabel_nodes(
                        l, {node: f"{graph_num}_{node}" for node in l.nodes()}
                    )
                    graph_num = graph_num + 1
                    Left_new.append(l)
                if len(Left_new) == 0:
                    break

                Left_new = nx.compose_all(Left_new)

                # Resolve unwanted nesting of map lists
                for n, data in Left_new.nodes(data=True):
                    if "map" in data:
                        Left_new.nodes[n]["map"] = flatten(Left_new.nodes[n]["map"])
                # Write DPO-Rules for all possible reaction atom to atom maps
                r_nodes_Left, m_nodes_Left = find_atommaps(Left_new)
                for r1 in Right1:
                    if len(r1) != 0:
                        # Compose all molecules at right side together
                        Right_new = []
                        graph_num = 0
                        for right in r1:
                            r = right.copy()
                            r = nx.relabel_nodes(
                                r, {node: f"{graph_num}_{node}" for node in r.nodes()}
                            )
                            graph_num = graph_num + 1
                            Right_new.append(r)
                        Right_new = nx.compose_all(Right_new)
                        if len(Right_new) == 0:
                            continue

                        # Resolve unwanted nesting of map lists
                        for n, data in Right_new.nodes(data=True):
                            if "map" in data:
                                Right_new.nodes[n]["map"] = flatten(
                                    Right_new.nodes[n]["map"]
                                )

                        # Write DPO-Rules for all possible reaction atom to atom maps
                        # Create Lists of the atom maps for all atom types
                        r_nodes_Right, m_nodes_Right = find_atommaps(Right_new)
                        if len(r_nodes_Left) == 0 or len(r_nodes_Right) == 0:
                            continue

                        if (
                            len(r_nodes_Left) < rclass_num
                            or len(r_nodes_Right) < rclass_num
                        ):
                            continue

                        if len(r_nodes_Left) != len(r_nodes_Right) or len(
                            m_nodes_Left
                        ) != len(m_nodes_Right):
                            continue
                        # Generate the DPO for all possible atom-to-atom maps
                        # Initialize lists with the entries to be written for each part of the DPO rule
                        Left_Text = []
                        Context_Text = []
                        Right_Text = []
                        node_no = 1  # New node number system for DPO rule
                        dict_node = (
                            {}
                        )  # Dictonary to rewrite edges with new number system

                        # Check for right Number of R-/D- and M-atoms
                        left_R_atoms = []
                        left_D_atoms = []
                        left_M_atoms = []
                        for n in Left_new.nodes():
                            data = Left_new.nodes[n]
                            if "atomtype" in data:
                                if data["label"] != "H":
                                    if "r" in data["atomtype"]:
                                        left_R_atoms.append(n)
                                    elif "d" in data["atomtype"]:
                                        left_D_atoms.append(n)
                                    elif "m" in data["atomtype"]:
                                        left_M_atoms.append(n)
                        right_R_atoms = []
                        right_D_atoms = []
                        right_M_atoms = []
                        for n in Right_new.nodes():
                            data = Right_new.nodes[n]
                            if "atomtype" in data:
                                if data["label"] != "H":
                                    if "r" in data["atomtype"]:
                                        right_R_atoms.append(n)
                                    elif "d" in data["atomtype"]:
                                        right_D_atoms.append(n)
                                    elif "m" in data["atomtype"]:
                                        right_M_atoms.append(n)

                        if len(left_R_atoms) != len(right_R_atoms) or len(
                            left_M_atoms
                        ) != len(right_M_atoms):
                            continue
                        # Write R-atoms
                        seen_nodes = []
                        found_maps = []
                        hydrogen_left = []
                        hydrogen_right = []

                        for m in r_nodes_Left:
                            m_o = m.split("_")
                            mp1 = int(m_o[0]) + 1
                            mm1 = int(m_o[0]) - 1
                            for node_l in Left_new.nodes():
                                attrs_l = Left_new.nodes[node_l]
                                if "map" not in attrs_l:
                                    continue
                                if m in attrs_l["map"] and attrs_l["label"] != "H":
                                    atom_l = attrs_l["label"]
                                    node_r = None
                                    for r in Right_new.nodes():
                                        attrs_r = Right_new.nodes[r]
                                        if "map" in attrs_r:
                                            map_list = []
                                            for m in attrs_r["map"]:
                                                m_o = m.split("_")
                                                mp2 = int(m_o[0])
                                                map_list.append(mp2)

                                            if (
                                                (mp1 in map_list or mm1 in map_list)
                                                and attrs_r["label"] == atom_l
                                                and r not in seen_nodes
                                                and "r" in attrs_r["atomtype"]
                                            ):
                                                node_r = r
                                                seen_nodes.append(node_r)
                                                found_maps.append(m)
                                                break

                                    if [node_l, node_r] not in dict_node.values():
                                        Context_Text.append(
                                            "  node [ id "
                                            + str(node_no)
                                            + ' label "'
                                            + atom_l
                                            + '" ]\n'
                                        )
                                        dict_node[node_no] = [node_l, node_r]
                                        node_no = node_no + 1

                        # Check if all R-atoms were found on both sides
                        if len(found_maps) != len(r_nodes_Left):
                            continue
                        # print('R')
                        # Write R-hydrogen
                        for node in Left_new.nodes():
                            if "map" in Left_new.nodes[node]:
                                if (
                                    Left_new.nodes[node]["label"] == "H"
                                    and "r" in Left_new.nodes[node]["atomtype"]
                                ):
                                    hydrogen_left.append(node)
                                    for r in range(
                                        len(Left_new.nodes[node]["atomtype"])
                                    ):
                                        if Left_new.nodes[node]["atomtype"][r] == "r":
                                            if (
                                                not isinstance(
                                                    Left_new.nodes[node]["map"], list
                                                )
                                                and r == 0
                                            ):
                                                m = Left_new.nodes[node]["map"]
                                            else:
                                                m = Left_new.nodes[node]["map"][r]
                                    m_o = m.split("_")
                                    mp1 = int(m_o[0]) + 1
                                    mp1 = str(mp1) + "_" + m_o[1]
                                    mm1 = int(m_o[0]) - 1
                                    mm1 = str(mm1) + "_" + m_o[1]
                                    found = False
                                    for node in Right_new.nodes():
                                        if "map" in Right_new.nodes[node]:
                                            if (
                                                (
                                                    mp1 in Right_new.nodes[node]["map"]
                                                    or mm1
                                                    in Right_new.nodes[node]["map"]
                                                )
                                                and "r"
                                                in Right_new.nodes[node]["atomtype"]
                                                and Right_new.nodes[node]["label"]
                                                == "H"
                                                and node not in seen_nodes
                                            ):
                                                found = True
                                                seen_nodes.append(node)
                                                hydrogen_right.append(node)
                                                break
                                    if found == False:
                                        hydrogen_right.append(None)

                        # find R-Hydrogens which are only one the right side
                        for node in Right_new.nodes():
                            if (
                                "map" in Right_new.nodes[node]
                                and node not in seen_nodes
                            ):
                                if (
                                    Right_new.nodes[node]["label"] == "H"
                                    and "r" in Right_new.nodes[node]["atomtype"]
                                ):
                                    seen_nodes.append(node)
                                    hydrogen_left.append(None)
                                    hydrogen_right.append(node)

                        # If no hydrogens got lost then write them in the Context
                        # If right graph has more hydrogens then write this one in right and as H+ in left
                        for i in range(len(hydrogen_left)):
                            if (
                                [hydrogen_left[i], hydrogen_right[i]]
                                not in dict_node.values()
                                and hydrogen_left[i] != None
                                and hydrogen_right[i] != None
                            ):
                                Context_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H" ]\n'
                                )
                                dict_node[node_no] = [
                                    hydrogen_left[i],
                                    hydrogen_right[i],
                                ]
                                node_no = node_no + 1
                            elif (
                                [hydrogen_left[i], hydrogen_right[i]]
                                not in dict_node.values()
                                and hydrogen_left[i] == None
                                and hydrogen_right[i] != None
                            ):
                                Left_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H+" ]\n'
                                )
                                Right_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H" ]\n'
                                )
                                dict_node[node_no] = [None, hydrogen_right[i]]
                                node_no = node_no + 1
                            elif (
                                [hydrogen_left[i], hydrogen_right[i]]
                                not in dict_node.values()
                                and hydrogen_left[i] != None
                                and hydrogen_right[i] == None
                            ):
                                Right_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H+" ]\n'
                                )
                                Left_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H" ]\n'
                                )
                                dict_node[node_no] = [hydrogen_left[i], None]
                                node_no = node_no + 1
                        # write D-atoms
                        # print('D')
                        D_check = True  # Turn into False if D-atom has no partner, and mark the side where the partner is missing
                        for node_l in Left_new.nodes():
                            attrs = Left_new.nodes[node_l]
                            if "map" not in attrs:
                                continue
                            if (
                                "d" in attrs["atomtype"]
                                and attrs["label"] != "H"
                                and D_check == True
                            ):
                                D_check = "Right"
                                atom_l = attrs["label"]
                                map_l = attrs["map"]
                                for node_r in Right_new.nodes():
                                    attrs = Right_new.nodes[node_r]
                                    if "map" not in attrs:
                                        continue
                                    if (
                                        attrs["map"] == map_l
                                        and atom_l == attrs["label"]
                                        and node_r not in seen_nodes
                                    ):
                                        D_check = True
                                        dict_node[node_no] = [node_l, node_r]
                                        Context_Text.append(
                                            "  node [ id "
                                            + str(node_no)
                                            + ' label "'
                                            + atom_l
                                            + '" ]\n'
                                        )
                                        node_no = node_no + 1
                                        seen_nodes.append(node_r)
                                        break

                        # Check D-atom partners also for the right side
                        for node_r in Right_new.nodes():
                            attrs = Right_new.nodes[node_r]
                            if "map" not in attrs:
                                continue
                            if (
                                "d" in attrs["atomtype"]
                                and attrs["label"] != "H"
                                and D_check == True
                            ):
                                if D_check == True:
                                    D_check = "Left"
                                else:
                                    D_check = "Both"
                                atom_r = attrs["label"]
                                map_r = attrs["map"]
                                for node_l in Left_new.nodes():
                                    attrs = Left_new.nodes[node_l]
                                    if "map" not in attrs:
                                        continue
                                    if (
                                        attrs["map"] == map_r
                                        and atom_r == attrs["label"]
                                        and node_l not in seen_nodes
                                    ):
                                        D_check = True

                        # If one D-atom has no partner then check for special CO2 Case
                        if D_check != True:
                            if D_check == "Left":
                                node_no, dict_node, Context_Text = c02_case(
                                    Left_new,
                                    Right_new,
                                    "L",
                                    "R",
                                    node_no,
                                    dict_node,
                                    Context_Text,
                                )
                            if D_check == "Right":
                                node_no, dict_node, Context_Text = c02_case(
                                    Right_new,
                                    Left_new,
                                    "R",
                                    "L",
                                    node_no,
                                    dict_node,
                                    Context_Text,
                                )
                            if D_check == "Both":
                                node_no, dict_node, Context_Text = c02_case(
                                    Left_new,
                                    Right_new,
                                    "L",
                                    "R",
                                    node_no,
                                    dict_node,
                                    Context_Text,
                                )
                                node_no, dict_node, Context_Text = c02_case(
                                    Right_new,
                                    Left_new,
                                    "R",
                                    "L",
                                    node_no,
                                    dict_node,
                                    Context_Text,
                                )
                            # Check again if D-Atoms have missing partners
                            for n in Left_new.nodes():
                                if "atomtype" in Left_new.nodes[n]:
                                    if (
                                        "d" in Left_new.nodes[n]["atomtype"]
                                        and Left_new.nodes[n]["label"] != "H"
                                    ):
                                        D_check = False
                                        for i in dict_node.values():
                                            if n in i:
                                                D_check = True
                                        if D_check == False:
                                            break

                            for n in Right_new.nodes():
                                if "atomtype" in Right_new.nodes[n]:
                                    if (
                                        "d" in Right_new.nodes[n]["atomtype"]
                                        and Right_new.nodes[n]["label"] != "H"
                                    ):
                                        D_check = False
                                        for i in dict_node.values():
                                            if n in i:
                                                D_check = True
                                        if D_check == False:
                                            break
                            # If D_check still fails then the Atom to Atom Map is not complete! Continue with the next Combination
                            if D_check == False:
                                continue
                        # Correction of small molecules
                        Context_Text, Left_Text, Right_Text, node_no = (
                            smallMolecule_Corr(
                                Left_new,
                                "L",
                                Context_Text,
                                Left_Text,
                                Right_Text,
                                node_no,
                                dict_node,
                            )
                        )
                        if Context_Text == False:
                            continue
                        Context_Text, Left_Text, Right_Text, node_no = (
                            smallMolecule_Corr(
                                Right_new,
                                "R",
                                Context_Text,
                                Left_Text,
                                Right_Text,
                                node_no,
                                dict_node,
                            )
                        )
                        if Context_Text == False:
                            continue
                        # print('M')
                        # Write M-Atoms
                        #  Initialize control structures
                        m_nodes_no = 0
                        for i in Left_new.nodes():
                            if "map" in Left_new.nodes[i]:
                                if (
                                    Left_new.nodes[i]["label"] != "H"
                                    and "d" not in Left_new.nodes[i]["atomtype"]
                                    and "r" not in Left_new.nodes[i]["atomtype"]
                                ):
                                    m_nodes_no = m_nodes_no + 1

                        seen_nodes = []
                        # Find the matching non hydrogen M-atoms
                        for m in range(len(m_nodes_Left)):
                            # Find node with Atommap ID
                            for node_l in Left_new.nodes():
                                attrs1 = Left_new.nodes[node_l]
                                if "map" in attrs1:
                                    if (
                                        m_nodes_Left[m] in attrs1["map"]
                                        and "r" not in attrs1["atomtype"]
                                        and attrs1["label"] != "H"
                                    ):
                                        atom_l = attrs1["label"]
                                        # Find corresponding M-atom on the right side
                                        for node in Right_new.nodes():
                                            attrs2 = Right_new.nodes[node]
                                            if "map" in attrs2:
                                                if (
                                                    m_nodes_Right[m] in attrs2["map"]
                                                    and attrs2["label"] == atom_l
                                                    and node not in seen_nodes
                                                ):
                                                    node_r = None
                                                    if len(attrs2["map"]) != len(
                                                        attrs1["map"]
                                                    ):
                                                        continue
                                                    if len(attrs2["map"]) == 1:
                                                        # if attrs1['label']  == 'O':
                                                        # bindings_Right = [d['label'] for u, v, d in Right_new.edges(data=True) if node in (u, v) and 'label' in d]
                                                        # bindings_Left = [d['label'] for u, v, d in Left_new.edges(data=True) if node_l in (u, v) and 'label' in d]
                                                        # if (2.0 in bindings_Right) ^ (2.0 in bindings_Left) == False:
                                                        #    continue
                                                        node_r = node
                                                        seen_nodes.append(node)
                                                    else:
                                                        # Check if all maps fit together
                                                        fit_check = True
                                                        for i in range(
                                                            len(attrs2["map"])
                                                        ):
                                                            if (
                                                                attrs2["map"][i]
                                                                in m_nodes_Right
                                                                and attrs1["map"][i]
                                                                in m_nodes_Left
                                                            ):
                                                                if m_nodes_Right.index(
                                                                    attrs2["map"][i]
                                                                ) != m_nodes_Left.index(
                                                                    attrs1["map"][i]
                                                                ):
                                                                    fit_check = False
                                                                    break
                                                            else:
                                                                fit_check = False
                                                                break
                                                        if fit_check == True:
                                                            # Check if the neighborhood is the same
                                                            m_neighbors_left = [
                                                                n
                                                                for n in Left_new.neighbors(
                                                                    node_l
                                                                )
                                                                if Left_new.nodes[
                                                                    n
                                                                ].get("label")
                                                                != "H"
                                                            ]
                                                            list_neighbor_atoms_left = [
                                                                Left_new.nodes[n].get(
                                                                    "label"
                                                                )
                                                                for n in m_neighbors_left
                                                            ]
                                                            m_neighbors_right = [
                                                                n
                                                                for n in Right_new.neighbors(
                                                                    node
                                                                )
                                                                if Right_new.nodes[
                                                                    n
                                                                ].get("label")
                                                                != "H"
                                                            ]
                                                            list_neighbor_atoms_right = [
                                                                Right_new.nodes[n].get(
                                                                    "label"
                                                                )
                                                                for n in m_neighbors_right
                                                            ]
                                                            if sorted(
                                                                list_neighbor_atoms_left
                                                            ) == sorted(
                                                                list_neighbor_atoms_right
                                                            ):
                                                                node_r = node
                                                                seen_nodes.append(node)
                                                    if node_r != None:
                                                        if [
                                                            node_l,
                                                            node_r,
                                                        ] not in dict_node.values():
                                                            Context_Text.append(
                                                                "  node [ id "
                                                                + str(node_no)
                                                                + ' label "'
                                                                + atom_l
                                                                + '" ]\n'
                                                            )
                                                            dict_node[node_no] = [
                                                                node_l,
                                                                node_r,
                                                            ]
                                                            node_no = node_no + 1
                                                            break

                        # Check if all M-atoms were found on both sides
                        if len(seen_nodes) != m_nodes_no:
                            continue
                        # Write M-hydrogens
                        hydrogen_left = []
                        hydrogen_right = []
                        for node in Left_new.nodes():
                            attrs = Left_new.nodes[node]
                            if "map" in attrs:
                                if (
                                    attrs["label"] == "H"
                                    and "r" not in attrs["atomtype"]
                                    and "d" not in attrs["atomtype"]
                                ):
                                    hydrogen_left.append(node)
                                    # Find the corresponding atom map ID for the M-hydrogen
                                    if not isinstance(attrs["atomtype"], list):
                                        attrs["atomtype"] = [attrs["atomtype"]]
                                    for r in range(len(attrs["atomtype"])):
                                        if attrs["atomtype"][r] == "m":
                                            if isinstance(attrs["map"], list):
                                                m = attrs["map"][r]
                                            else:
                                                m = attrs["map"]
                                        for n in range(len(m_nodes_Left)):
                                            if m_nodes_Left[n] == m:
                                                mp1 = m_nodes_Right[n]

                                    found = False
                                    for node in Right_new.nodes():
                                        attrs = Right_new.nodes[node]
                                        if "map" in attrs:
                                            if (
                                                mp1 in attrs["map"]
                                                and "r" not in attrs["atomtype"]
                                                and "d" not in attrs["atomtype"]
                                                and attrs["label"] == "H"
                                                and node not in seen_nodes
                                            ):
                                                found = True
                                                seen_nodes.append(node)
                                                hydrogen_right.append(node)
                                                break
                                    if found == False:
                                        hydrogen_right.append(None)

<<<<<<< HEAD
                        for node in Right_new.nodes():
                            attrs = Right_new.nodes[node]
                            if "map" in attrs and node not in seen_nodes:
                                if (
                                    attrs["label"] == "H"
                                    and "r" not in attrs["atomtype"]
                                    and "d" not in attrs["atomtype"]
                                ):
                                    seen_nodes.append(node)
                                    hydrogen_left.append(None)
                                    hydrogen_right.append(node)
=======
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
                                
                                # Check for common edges
                                seen_edges = []
                                for e1,e2 in trans_Ledges.keys():
                                    # First node order e1,e2
                                    if (e1,e2) in trans_Redges.keys() and (e1,e2) not in seen_edges:
                                        seen_edges.append((e1,e2))
                                        if trans_Ledges[(e1,e2)] == trans_Redges[(e1,e2)]:
                                            Context_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                        else:
                                            Left_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                            Right_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Redges[(e1,e2)]+ '" ]\n')
                                    # Second node order e2,e1
                                    elif (e2,e1) in trans_Redges.keys() and (e2,e1) not in seen_edges:
                                        seen_edges.append((e2,e1))
                                        if trans_Ledges[(e1,e2)] == trans_Redges[(e2,e1)]:
                                            Context_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                        else:
                                            Left_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')
                                            Right_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Redges[(e2,e1)]+ '" ]\n')
                               
                                # Check for edges only of one side       
                                for e1,e2 in trans_Ledges.keys():             
                                    if (e1,e2) not in seen_edges and (e2,e1) not in seen_edges:
                                        Left_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Ledges[(e1,e2)] + '" ]\n')                 
                                for e1,e2 in trans_Redges.keys():                
                                    if (e1,e2) not in seen_edges and (e2,e1) not in seen_edges:
                                        Right_Text.append('  edge [ source '+ e1 + ' target '+ e2 + ' label "' + trans_Redges[(e1,e2)] + '" ]\n')         
                                              
                                # Check if changes are found on left or right side
                                if len(Left_Text) == 0 and len(Right_Text) == 0: 
                                    continue         
                                
                                # Write gml file                       
                                with open(output_path+rn+'/'+rn+'_'+str(no)+'.gml', 'w') as out: 
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
>>>>>>> dbf44ddee9f253e25b1bdc61ff6cc036e1c13eb1

                        # If no hydrogens got lost then write them in the 'Context'
                        # If right graph has more hydrogens then write this into 'Right' and as H+ in 'Left'
                        for i in range(len(hydrogen_left)):
                            if (
                                [hydrogen_left[i], hydrogen_right[i]]
                                not in dict_node.values()
                                and hydrogen_left[i] != None
                                and hydrogen_right[i] != None
                            ):
                                Context_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H" ]\n'
                                )
                                dict_node[node_no] = [
                                    hydrogen_left[i],
                                    hydrogen_right[i],
                                ]
                                node_no = node_no + 1
                            if (
                                [hydrogen_left[i], hydrogen_right[i]]
                                not in dict_node.values()
                                and hydrogen_left[i] == None
                                and hydrogen_right[i] != None
                            ):
                                Left_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H+" ]\n'
                                )
                                Right_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H" ]\n'
                                )
                                dict_node[node_no] = [None, hydrogen_right[i]]
                                node_no = node_no + 1
                            if (
                                [hydrogen_left[i], hydrogen_right[i]]
                                not in dict_node.values()
                                and hydrogen_left[i] != None
                                and hydrogen_right[i] == None
                            ):
                                Right_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H+" ]\n'
                                )
                                Left_Text.append(
                                    "  node [ id " + str(node_no) + ' label "H" ]\n'
                                )
                                dict_node[node_no] = [hydrogen_left[i], None]
                                node_no = node_no + 1

                        # Translate numeric edge labelling into symbolic labelling
                        L1 = update_bound_attribute(copy.deepcopy(Left_new))
                        R1 = update_bound_attribute(copy.deepcopy(Right_new))

                        # Write edges
                        trans_Ledges = {}
                        trans_Redges = {}
                        ori_Ledges = list(L1.edges(data=True))
                        ori_Redges = list(R1.edges(data=True))

                        # Translate edges with new node number system
                        for edge in L1.edges():
                            node1 = None
                            node2 = None
                            for i in dict_node:
                                if edge[0] == dict_node[i][0]:
                                    node1 = i
                                if edge[1] == dict_node[i][0]:
                                    node2 = i
                            if node1 == None or node2 == None:
                                L1.remove_edge(edge[0], edge[1])
                            else:
                                trans_Ledges[(str(node1), str(node2))] = L1.edges[edge][
                                    "label"
                                ]

                        for edge in R1.edges():
                            node1 = None
                            node2 = None
                            for i in dict_node:
                                if edge[0] == dict_node[i][1]:
                                    node1 = i
                                if edge[1] == dict_node[i][1]:
                                    node2 = i
                            if node1 == None or node2 == None:
                                R1.remove_edge(edge[0], edge[1])
                            else:
                                trans_Redges[(str(node1), str(node2))] = R1.edges[edge][
                                    "label"
                                ]
                        # Check for common edges
                        seen_edges = []
                        for e1, e2 in trans_Ledges.keys():
                            # First node order e1,e2
                            if (e1, e2) in trans_Redges.keys() and (
                                e1,
                                e2,
                            ) not in seen_edges:
                                seen_edges.append((e1, e2))
                                if trans_Ledges[(e1, e2)] == trans_Redges[(e1, e2)]:
                                    Context_Text.append(
                                        "  edge [ source "
                                        + e1
                                        + " target "
                                        + e2
                                        + ' label "'
                                        + trans_Ledges[(e1, e2)]
                                        + '" ]\n'
                                    )
                                else:
                                    Left_Text.append(
                                        "  edge [ source "
                                        + e1
                                        + " target "
                                        + e2
                                        + ' label "'
                                        + trans_Ledges[(e1, e2)]
                                        + '" ]\n'
                                    )
                                    Right_Text.append(
                                        "  edge [ source "
                                        + e1
                                        + " target "
                                        + e2
                                        + ' label "'
                                        + trans_Redges[(e1, e2)]
                                        + '" ]\n'
                                    )
                            # Second node order e2,e1
                            elif (e2, e1) in trans_Redges.keys() and (
                                e2,
                                e1,
                            ) not in seen_edges:
                                seen_edges.append((e2, e1))
                                if trans_Ledges[(e1, e2)] == trans_Redges[(e2, e1)]:
                                    Context_Text.append(
                                        "  edge [ source "
                                        + e1
                                        + " target "
                                        + e2
                                        + ' label "'
                                        + trans_Ledges[(e1, e2)]
                                        + '" ]\n'
                                    )
                                else:
                                    Left_Text.append(
                                        "  edge [ source "
                                        + e1
                                        + " target "
                                        + e2
                                        + ' label "'
                                        + trans_Ledges[(e1, e2)]
                                        + '" ]\n'
                                    )
                                    Right_Text.append(
                                        "  edge [ source "
                                        + e1
                                        + " target "
                                        + e2
                                        + ' label "'
                                        + trans_Redges[(e2, e1)]
                                        + '" ]\n'
                                    )

                        # Check for edges only of one side
                        for e1, e2 in trans_Ledges.keys():
                            if (e1, e2) not in seen_edges and (
                                e2,
                                e1,
                            ) not in seen_edges:
                                Left_Text.append(
                                    "  edge [ source "
                                    + e1
                                    + " target "
                                    + e2
                                    + ' label "'
                                    + trans_Ledges[(e1, e2)]
                                    + '" ]\n'
                                )
                        for e1, e2 in trans_Redges.keys():
                            if (e1, e2) not in seen_edges and (
                                e2,
                                e1,
                            ) not in seen_edges:
                                Right_Text.append(
                                    "  edge [ source "
                                    + e1
                                    + " target "
                                    + e2
                                    + ' label "'
                                    + trans_Redges[(e1, e2)]
                                    + '" ]\n'
                                )

                        # Check if changes are found on left or right side
                        if len(Left_Text) == 0 and len(Right_Text) == 0:
                            continue
                        # Write gml file
                        if not os.path.exists(output_path + rn):
                            os.makedirs(output_path + rn)

                        with open(
                            output_path + rn + "/" + rn + "_" + str(no) + ".gml", "w"
                        ) as out:
                            out.write("rule [\n")
                            out.write(' ruleID "' + rn + '"\n')
                            out.write(" left [\n")
                            for i in Left_Text:
                                out.write(i)
                            out.write(" ]\n")
                            out.write(" context [\n")
                            for i in Context_Text:
                                out.write(i)
                            out.write(" ]\n")
                            out.write(" right [\n")
                            for i in Right_Text:
                                out.write(i)
                            out.write(" ]\n")
                            out.write("]\n")
                        succ = True
                        no = no + 1

    if succ == False:
        if toBig == True:
            no_reaction_toBig = no_reaction_toBig + 1
            with open(log_big, "a") as log:
                log.write(str(rn) + " No. of graphs:" + str(count_graph) + "\n")
        else:
            no_reactions_false = no_reactions_false + 1
            with open(log_unsucc, "a") as log:
                log.write(str(rn) + "\n")

    else:
        no_reactions_true = no_reactions_true + 1

print("\n############# Statistics #############")
print("No. of TOTAL Rections", no_reactions)
print("No. of created DPO-Rules", no_reactions_true)
print("No. of Rections were no DPO-Rule could be created", no_reactions_false)
print("No. of Reactions which to big to create", no_reaction_toBig)
print("No. of Reaction with missing Compounds", missing_comp)
