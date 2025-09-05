import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import glob
import itertools
import os
import sys
from copy import copy
from networkx.algorithms import isomorphism as iso
from itertools import combinations_with_replacement
from itertools import product
from tqdm import tqdm

############################################
# code only work with python version > 3.6 #
############################################


# shows the progress over the algorithm
def progress_bar(progress, total, length=40):
    percent = int(progress / total * 100)
    filled = int(length * progress // total)
    bar = "#" * filled + "-" * (length - filled)
    sys.stdout.write(f"\r[{bar}] {percent}%")
    sys.stdout.flush()


def test_new_Graph(correct_graphs, G):
    def edge_match(a, b):
        return a["bound"] == b["bound"]

    def node_match(a, b):
        return a["atom"] == b["atom"]

    for g in correct_graphs:
        GM = iso.GraphMatcher(g, G, node_match=node_match, edge_match=edge_match)
        is_isomorphic = GM.is_isomorphic()
        if is_isomorphic == True:
            return False
    return True


def count_sameLabelLeafs(G):
    label_count = {}
    label_nodes = {}
    for n in G.nodes(data=True):
        if "label" in n[1]:
            if n[1]["label"] in label_count.keys():
                label_count[n[1]["label"]] = label_count[n[1]["label"]] + 1
                label_nodes[n[1]["label"]].append(n[0])
            else:
                label_count[n[1]["label"]] = 1
                label_nodes[n[1]["label"]] = [n[0]]
    return label_count, label_nodes


def translate_LabelNode(component, node, label_path, i):
    # Defined attributes of label
    node_name = component.nodes[node]["label"]
    side = component.nodes[node]["side"]
    aam = component.nodes[node]["map"]
    typ = component.nodes[node]["atomtype"]

    # Load label graph
    with open(label_path + node_name + ".gml", "r") as file:
        gml_code = file.read()
    subgra = nx.parse_gml(gml_code)

    for sub_node, attrs in subgra.nodes(data=True):
        # Copy information from RDM Pattern node
        new_node_id = f"subgraph_{sub_node}" + str(i)
        component.add_node(new_node_id, **attrs)
        component.nodes[new_node_id]["side"] = side
        component.nodes[new_node_id]["map"] = aam
        component.nodes[new_node_id]["atomtype"] = typ
        # conacted the Labelgraph with the RDM pattern tree
        if component.nodes[new_node_id]["status"] == "1":
            neighbors_label = list(component.neighbors(node))
            for n in neighbors_label:
                component.add_edge(n, new_node_id, kind=3)
            component.remove_node(node)

    # Add the edges of the subgraph to the component
    for sub_edge in subgra.edges(data=True):
        component.add_edge(
            f"subgraph_{sub_edge[0]}" + str(i),
            f"subgraph_{sub_edge[1]}" + str(i),
            bound=sub_edge[2]["bound"],
        )


def find_neighborhood(G, node, forbidden):

    neighbor_nodes = list(G.neighbors(node))
    neighbor_nodes.remove(forbidden)

    # remove all neighbors with edges incl Attribute 'kind':3
    for edge in G.edges(node, data=True):
        u, v, attr = edge
        if "kind" in attr:
            if u != node and u != forbidden:
                neighbor_nodes.remove(u)
            if v != node and v != forbidden:
                neighbor_nodes.remove(v)

    return neighbor_nodes


def compare_neighborhood(G, neighbor_node1, neighbor_node2, found_node, partner):

    # Create a list of the atoms with their bonding from the neighborhood
    def get_edge_attributes_for_node(graph, node, forbidden):
        edge_attributes = []
        for edge in graph.edges(node, data=True):
            u, v, attr = edge
            if "bound" in attr:
                if u == node and v != forbidden:
                    edge_attributes.append([G.nodes[v]["atom"], attr["bound"]])
                if v == node and u != forbidden:
                    edge_attributes.append([G.nodes[u]["atom"], attr["bound"]])
        return edge_attributes

    """ Check if the neighborhoods are fit together, i.e. whether they have the same atoms as neighbors with the same bonding. 
    Exception: The neighbornodes can have also rest in the place of an atom"""

    def neighborhood_check(neighborhood1, neighborhood2):
        # check for same neighbor nodes
        combined_hist = {}
        itemset = set()
        for n1 in neighborhood1:
            itemset.add((n1[0], n1[1]))
        for n2 in neighborhood2:
            itemset.add((n2[0], n2[1]))

        for item in itemset:
            combined_hist[item] = 0

        bond_hist_r1 = {"-": 0, "=": 0, ":": 0, "r": 0, "dr": 0, "#": 0}
        bond_hist_r2 = {"-": 0, "=": 0, ":": 0, "r": 0, "dr": 0, "#": 0}
        for n1 in neighborhood1:
            if n1[0] != "R":
                combined_hist[(n1[0], n1[1])] = combined_hist[(n1[0], n1[1])] + 1
            else:
                bond_hist_r1[n1[1]] = bond_hist_r1[n1[1]] + 1
        for n2 in neighborhood2:
            if n2[0] != "R":
                combined_hist[(n2[0], n2[1])] = combined_hist[(n2[0], n2[1])] - 1
            else:
                bond_hist_r2[n2[1]] = bond_hist_r2[n2[1]] + 1

        no_match_n1 = {}
        no_match_n2 = {}
        for n in combined_hist:
            if combined_hist[n] > 0:
                no_match_n1[n] = combined_hist[n]
            elif combined_hist[n] < 0:
                no_match_n2[n] = -1 * combined_hist[n]

        for n1 in no_match_n1:
            bond_hist_r2[n1[1]] = bond_hist_r2[n1[1]] - no_match_n1[n1]

        for n2 in no_match_n2:
            bond_hist_r1[n2[1]] = bond_hist_r1[n2[1]] - no_match_n2[n2]

        for b in bond_hist_r1:
            if bond_hist_r1[b] != bond_hist_r2[b]:
                return False

        return True

    found_neighborhood = get_edge_attributes_for_node(G, found_node, neighbor_node1)
    partner_neighborhood = get_edge_attributes_for_node(G, partner, neighbor_node2)
    neighbor_node1_neighborhood = get_edge_attributes_for_node(
        G, neighbor_node1, found_node
    )
    neighbor_node2_neighborhood = get_edge_attributes_for_node(
        G, neighbor_node2, partner
    )
    check_neighborhood = False

    # Case 1: N1 = atom and N2 = atom
    if G.nodes[neighbor_node2]["atom"] == G.nodes[found_node]["atom"]:
        if G.nodes[neighbor_node1]["atom"] == G.nodes[partner]["atom"]:
            if (
                neighborhood_check(found_neighborhood, neighbor_node2_neighborhood)
                == True
            ):
                if (
                    neighborhood_check(
                        partner_neighborhood, neighbor_node1_neighborhood
                    )
                    == True
                ):
                    check_neighborhood = True
            # Special case: If N1 and N2 are atoms and have no neighbors than pair is also right
            if check_neighborhood == False:
                if len(neighbor_node2_neighborhood) == 0:
                    if len(neighbor_node1_neighborhood) == 0:
                        check_neighborhood = True
    # Case 2: N1 = atom and N2 = R
    if G.nodes[neighbor_node2]["atom"] == "R":
        if G.nodes[neighbor_node1]["atom"] != "R":
            if G.nodes[neighbor_node1]["atom"] != "H":
                if (
                    neighborhood_check(
                        partner_neighborhood, neighbor_node1_neighborhood
                    )
                    == True
                ):
                    check_neighborhood = True
            # Special case
            if check_neighborhood == False:
                if len(neighbor_node2_neighborhood) == 0:
                    if len(neighbor_node1_neighborhood) == 0:
                        check_neighborhood = True
    # Case 3: N1 = R and N2 = atom
    if G.nodes[neighbor_node1]["atom"] == "R":
        if G.nodes[neighbor_node2]["atom"] != "R":
            if G.nodes[neighbor_node2]["atom"] != "H":
                if (
                    neighborhood_check(found_neighborhood, neighbor_node2_neighborhood)
                    == True
                ):
                    check_neighborhood = True
                # Special case
                if check_neighborhood == False:
                    if len(neighbor_node2_neighborhood) == 0:
                        if len(neighbor_node1_neighborhood) == 0:
                            check_neighborhood = True
    # Case 4: N1 and N2 are R
    if G.nodes[neighbor_node1]["atom"] == "R":
        if G.nodes[neighbor_node2]["atom"] == "R":
            check_neighborhood = True

    return check_neighborhood


def compare_bonds(G, node_list, found_node, partner):
    node_list_filtered = []
    for node1, node2 in node_list:
        if G.edges[found_node, node1].get("bound") == G.edges[partner, node2].get(
            "bound"
        ):
            node_list_filtered.append((node1, node2))

    return node_list_filtered


def find_neighbornodes(G, found_node, partner):
    # Find all possible neighbor nodes for connecting: combi = produkt(N1,N2)
    neighbor_nodes1 = find_neighborhood(G, found_node, partner)
    neighbor_nodes2 = find_neighborhood(G, partner, found_node)
    neighbor_combis = product(neighbor_nodes1, neighbor_nodes2)
    neighbor_combis = list(neighbor_combis)

    final_atom_to_atom = []
    # Try first to connect from atom to atom
    for x in neighbor_combis:
        atom_to_atom = []
        if G.nodes[x[0]]["atom"] != "H":
            if G.nodes[x[1]]["atom"] == G.nodes[found_node]["atom"]:
                if G.nodes[x[0]]["atom"] == G.nodes[partner]["atom"]:
                    atom_to_atom.append(x)
        # check the atom to atom connections are valid with the neighborhood
        for pair in atom_to_atom:
            check_neighborhood = compare_neighborhood(
                G, pair[0], pair[1], found_node, partner
            )
            if check_neighborhood == True:
                final_atom_to_atom.append(pair)

    # If no atom to atom bond possible then try atom to rest or rest to rest
    if len(final_atom_to_atom) == 0:
        atom_to_atom = []
        for x in neighbor_combis:
            if (
                G.nodes[x[0]]["atom"] != "H"
                and G.nodes[x[1]]["atom"] != "H"
                and (G.nodes[x[0]]["atom"] == "R" or G.nodes[x[1]]["atom"] == "R")
            ):
                atom_to_atom.append(x)

        # Check if the atom to atom connections with the neighborhood are valid
        for pair in atom_to_atom:
            check_neighborhood = compare_neighborhood(
                G, pair[0], pair[1], found_node, partner
            )
            if check_neighborhood == True:
                final_atom_to_atom.append(pair)

    # Filter of pairs with different bondings
    final_atom_to_atom = compare_bonds(G, final_atom_to_atom, found_node, partner)

    return final_atom_to_atom


def connected_kind3_Edge(G, neighbor_node1, neighbor_node2, found_node, partner):

    check_edge = True

    # Modify the graph so that edge kind 3 got fused together
    if neighbor_node1 != False and neighbor_node1 in G.nodes():
        bound1 = G.edges[found_node, neighbor_node1]["bound"]
        delete_nodes = find_neighborhood(G, neighbor_node1, found_node)
        delete_nodes.append(neighbor_node1)
        # Delete neighbor_node1 and neighbor_node2 and the neighborhood which would be unconnected afterwards
        for n in delete_nodes:
            G.remove_node(n)
    else:
        check_edge = False

    if neighbor_node2 != False and neighbor_node2 in G.nodes():
        delete_nodes = find_neighborhood(G, neighbor_node2, partner)
        delete_nodes.append(neighbor_node2)
        # Delete neighbor_node1 and neighbor_node2 and the neighborhood which would be unconnected afterwards
        for n in delete_nodes:
            G.remove_node(n)
    else:
        check_edge = False

    # Rename the kind 3 edge to a molecule edge
    if check_edge == True:
        G.edges[found_node, partner]["kind"] = "#"
        G.edges[found_node, partner]["bound"] = bound1

        # Check if the graph has loose nodes, the preceding neighborhood check ensures that these occur twice and once in a bound manner in the graph and can therefore be deleted.     if not nx.is_connected(G):
        isolated_nodes = [node for node in G.nodes if G.degree(node) == 0]
        for node in isolated_nodes:
            G.remove_node(node)

    return check_edge


def translateMolecule(G, correct_graphs, forrest_count, print_no, root, print_r):

    # Choose random edge with attribute 'kind':3 and find node and partner
    found = False
    for found_node, partner, attr in G.edges(data=True):
        if "kind" in attr and attr["kind"] == 3:
            kind3 = (found_node, partner)
            # find nodes for connecting
            atom_to_atom = find_neighbornodes(G, found_node, partner)
            found = True
            break

    if not found:
        iso_test = test_new_Graph(correct_graphs, G)

        if iso_test == True:
            correct_graphs.append(G)
        return

    check_edge = False
    # Try the different paths to connect the subgraph at the edges with attribute 'kind':3
    for neighbor_found, neighbor_partner in atom_to_atom:
        newgraph = G.copy()

        con_check = connected_kind3_Edge(
            newgraph, neighbor_found, neighbor_partner, kind3[0], kind3[1]
        )
        if not con_check:
            return
        translateMolecule(
            newgraph, correct_graphs, forrest_count, print_no, root, print_r
        )


def check_allAlternativCombis(
    G, list_alternatives, rxn, counter_work, counter_notwork, counter_multi
):

    forrest_count = 0  # Number of the subgraphs (in file name)
    sub_correct = True  # Turn into False if one subgraph could not be built
    multi = False  # Turn into True if more than one graph per subgraph could be built

    # Separate connected compounds and replace nodes with label graphs
    connected_components = list(nx.connected_components(G))
    subgraphs = [G.subgraph(component) for component in connected_components]
    translation_count = []
    print_r = 100
    print_no = 0

    for comp in subgraphs:
        forrest_count = forrest_count + 1
        round_num = 1  # None of version (in folder name)
        correct_graphs = []  # List of generated molecule structures

        # Find the root node and try all alternative label graphs
        for node in comp.nodes():
            if comp.nodes[node]["atomtype"] == "r":
                root = 0
                for alt in list_alternatives:
                    root = root + 1
                    # Copy before transforming root label
                    component = comp.copy()
                    sub_num = (
                        1  # Number is needed to rename the node IDs with unique numbers
                    )
                    translate_LabelNode(component, node, alt, sub_num)
                    sub_num = sub_num + 1
                    # Create all combinations of leave labeling
                    entry_count, label_to_nodes = count_sameLabelLeafs(component)
                    label_list = list(entry_count.keys())
                    combi_lists = {}
                    for j in entry_count:
                        combi = list(
                            combinations_with_replacement(
                                range(0, len(list_alternatives)), entry_count[j]
                            )
                        )
                        combi_lists[j] = combi
                    # Check if the root has only one leaf, otherwise find the labeling combis for all
                    if len(combi_lists) == 1:
                        for comb in combi_lists[
                            label_list[0]
                        ]:  # comb can be a vector e.g. (0,1)
                            component_alt = (
                                component.copy()
                            )  # create second copy for each component and translate the rest

                            for i in range(
                                len(comb)
                            ):  # each position in the combination
                                alt_index = comb[i]
                                translate_LabelNode(
                                    component_alt,
                                    label_to_nodes[label_list[0]][i],
                                    list_alternatives[alt_index],
                                    sub_num,
                                )
                                sub_num = sub_num + 1

                            # Try to translate this variant of subgraph
                            translateMolecule(
                                component_alt,
                                correct_graphs,
                                forrest_count,
                                print_no,
                                root,
                                print_r,
                            )
                    else:
                        # Find all label combis for all leaves
                        combis = list(
                            product(
                                combi_lists[label_list[0]], combi_lists[label_list[1]]
                            )
                        )
                        for i in range(2, len(label_list)):
                            combis = [
                                (*a, b)
                                for (a, b) in product(
                                    combis, combi_lists[label_list[i]]
                                )
                            ]

                        # Load all combinations of label graphs for the leafs
                        for (
                            comb
                        ) in (
                            combis
                        ):  # For a combination of alternative labels, e.g. ((2,), (1, 1))

                            component_alt = (
                                component.copy()
                            )  # create second copy for each component and translate the rest
                            for c, l in zip(comb, label_list):
                                for i, alt_index in enumerate(c):
                                    translate_LabelNode(
                                        component_alt,
                                        label_to_nodes[l][i],
                                        list_alternatives[alt_index],
                                        sub_num,
                                    )
                                    sub_num = sub_num + 1

                            # Try to translate this variant of subgraph
                            translateMolecule(
                                component_alt,
                                correct_graphs,
                                forrest_count,
                                print_no,
                                root,
                                print_r,
                            )
        # Save graph
        round_num = 1
        if len(correct_graphs) == 0:
            sub_correct = False
            with open("01_RCLASS_Error", "a") as file:
                file.write(
                    rxn
                    + " Rule with the Subgraph No. "
                    + str(forrest_count)
                    + " doesnt work!"
                    + "\n"
                )

        r_check = False
        for i in correct_graphs:
            # Finds all molecular subgraphs which were earlier created and load them in left and right lists
            r_count = 0
            for node, data in i.nodes(data=True):
                if data["atomtype"] == "r":
                    if data["atom"] != "H":
                        r_count += 1
            if r_count > 1:
                # Special case C6a with O6a (Carboxy-Endings)
                for node, data in i.nodes(data=True):
                    if data["atomtype"] == "r" and data["atom"] == "C":
                        neighbors = list(i.neighbors(node))
                        if len(neighbors) == 3:
                            correct_neighbors = []
                            for n in neighbors:
                                correct_neighbors.append(
                                    (i.nodes[n]["atom"], i.edges[node, n]["bound"])
                                )
                            if ("O", "-") and ("O", "=") in correct_neighbors:
                                saveRule(i, rxn, forrest_count, round_num)
                                round_num = round_num + 1
                                r_check = True
            else:
                saveRule(i, rxn, forrest_count, round_num)
                round_num = round_num + 1
                r_check = True

        if r_check == False:
            sub_correct = False

        if round_num > 2:
            multi = True

    # Counting for statistics
    if sub_correct == True:
        counter_work = counter_work + 1
    else:
        counter_notwork = counter_notwork + 1
    if multi == True:
        counter_multi = counter_multi + 1
    return counter_work, counter_notwork, counter_multi


def saveRule(G1, rxn, forrest_count, round_num):
    # Save the graphs as GML
    # Check if the 'side' attribute is the same for all nodes in the
    unique_side_values = set(G1.nodes[node]["side"] for node in G1.nodes())
    # Delete unnecessary attributes
    for node in G1.nodes():
        if "kind" in G1.nodes[node]:
            del G1.nodes[node]["kind"]
        if "status" in G1.nodes[node]:
            del G1.nodes[node]["status"]
        if "side" in G1.nodes[node]:
            del G1.nodes[node]["side"]
        for u, v, attr in G1.edges(data=True):
            if "kind" in attr:
                del G1[u][v]["kind"]
    # If all nodes have the same 'side' attribute, save the component as a GML file
    folder = folder_path.split(".gexf")
    new_path = path + "/01_RDM_Graphs/" + rxn + "_" + str(round_num) + "/"
    if os.path.exists(new_path) == False:
        os.makedirs(new_path)
    filename = (
        new_path
        + str(rxn)
        + "_"
        + str(unique_side_values)
        + "_"
        + str(forrest_count)
        + ".gml"
    )
    nx.write_gml(G1, filename)


#### Different Edge Kinds ####
# 1 = Edge between RCLASS label node with name giving Node of the label subgraph
# 2 = edges of the subgraph
# 3 = edges of the graph

# Paths for Labelgraphs
path_ori = "./Labels/01_Labelgraphs/"
path_alt1 = "./Labels/01_Labelgraphs_alternativ/"
path_alt2 = "./Labels/01_Labelgraphs_alternativ2/"
path_alt3 = "./Labels/01_Labelgraphs_alternativ3/"
path_alt4 = "./Labels/01_Labelgraphs_alternativ4/"
list_alternatives = [path_ori, path_alt1, path_alt2, path_alt3, path_alt4]

# Load RDM pattern trees
folder_path = sys.argv[1]
file_names = glob.glob(folder_path + "/*.gexf")
graphs = []
out_name = []

# create output folder
path = os.path.dirname(folder_path)
if not os.path.exists(path + "/01_RDM_Graphs/"):
    os.makedirs(path + "/01_RDM_Graphs/")

for file_name in file_names:
    rxn = file_name.split(".gexf")
    rxn = rxn[0].split("/")
    out_name.append(rxn[len(rxn) - 1])
    graph = nx.read_gexf(file_name)
    graphs.append(graph)
out_name_old = out_name

# Counter for statistics
counter_work = 0
counter_notwork = 0
counter_multi = 0

# Translate RCLASS in molecule structure
total = len(graphs)
status = 0
for graph in graphs:
    status += 1
    progress_bar(status, total)
    # Load the RDM pattern trees
    G = nx.Graph()
    node_attrList = []
    for node, attrs in graph.nodes(data=True):
        new_node_id = f"graph_{node}"
        G.add_node(new_node_id, **attrs)
        node_attrList.append(attrs["label"])
    for edge in graph.edges():
        new_edge = tuple(f"graph_{node}" for node in edge)
        G.add_edge(*new_edge, kind=3)

    rxn = out_name.pop(0)
    # rxn = rxn.split('/')
    # rxn = rxn[2]
    multi_check = False

    # Try to translate the RDM trees in molecular structure
    counter_work, counter_notwork, counter_multi = check_allAlternativCombis(
        G, list_alternatives, rxn, counter_work, counter_notwork, counter_multi
    )

print("\nNumber of TOTAL read Files:", len(file_names))
print("Number of constructed Molecule-Substructure: ", counter_work)
print("Number of constructed but multiple solutions Structures:", counter_multi)
print("not working Molecules: ", counter_notwork)
