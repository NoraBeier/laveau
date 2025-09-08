import networkx as nx
import matplotlib.pyplot as plt
import re
import os
import sys


# Function to parse the GML-like text and create a graph
def parse_gml_and_create_graphs(gml_content):
    context_nodes = {}
    context_edges = []
    left_nodes = {}
    left_edges = []
    right_nodes = {}
    right_edges = []
    constrains = False
    # Define the patterns for the nodes and edges
    node_pattern = re.compile(r'node \[ id (\d+) label "([^"]+)" \]')
    edge_pattern = re.compile(r'edge \[ source (\d+) target (\d+) label "([^"]+)" \]')

    current_section = None
    with open(gml_content, "r") as f:
        lines = f.readlines()
    # Go through the content line by line
    for line in lines:
        line = line.strip()
        if line.startswith("context ["):
            current_section = "context"
        elif line.startswith("left ["):
            current_section = "left"
        elif line.startswith("right ["):
            current_section = "right"
        elif line == "]":
            current_section = None
        if line.startswith("ruleID"):
            rule = line.split('"')[1]
        if line.startswith("constrainLabelAny"):
            constrains = True
        # Process nodes
        node_match = node_pattern.match(line)
        if node_match:
            node_id = int(node_match.group(1))
            label = node_match.group(2)
            if current_section == "context":
                context_nodes[node_id] = label
            elif current_section == "left":
                left_nodes[node_id] = label
            elif current_section == "right":
                right_nodes[node_id] = label

        # Process edges
        edge_match = edge_pattern.match(line)
        if edge_match:
            source = int(edge_match.group(1))
            target = int(edge_match.group(2))
            label = edge_match.group(3)
            if current_section == "context":
                context_edges.append((source, target, label))
            elif current_section == "left":
                left_edges.append((source, target, label))
            elif current_section == "right":
                right_edges.append((source, target, label))

    # Create the graphs
    G_left = create_graph(context_nodes, context_edges, left_nodes, left_edges, "left")
    G_right = create_graph(
        context_nodes, context_edges, right_nodes, right_edges, "right"
    )

    return G_left, G_right, rule, constrains


# Function to create a graph based on context and another area (left or right)
def create_graph(
    context_nodes, context_edges, additional_nodes, additional_edges, label
):
    G = nx.Graph()

    # Add the nodes from the context
    for node_id, node_label in context_nodes.items():
        G.add_node(node_id, label=node_label, context=True)

    # Add the edges from the context
    for source, target, edge_label in context_edges:
        G.add_edge(source, target, label=edge_label, context=True)

    # Add the additional nodes (left or right)
    for node_id, node_label in additional_nodes.items():
        G.add_node(node_id, label=node_label, context=False, rule=label)

    # Add the additional edges (left or right)
    for source, target, edge_label in additional_edges:
        G.add_edge(source, target, label=edge_label, context=False, rule=label)

    return G


# Function for comparing graphs and finding swapped neighbors
def find_swapped_neighbors(G1, G2):
    swapped_edges = []

    for u, v in G1.edges():
        if G2.has_edge(u, v):
            # Find the neighbors of u and v in both graphs
            neighbors_G1_u = set(G1.neighbors(u)) - {v}
            neighbors_G1_v = set(G1.neighbors(v)) - {u}
            neighbors_G2_u = set(G2.neighbors(u)) - {v}
            neighbors_G2_v = set(G2.neighbors(v)) - {u}

            # Check whether the neighbors are swapped
            if neighbors_G1_u == neighbors_G2_v and neighbors_G1_v == neighbors_G2_u:
                swapped_edges.append((u, v))

    return swapped_edges


def find_Node_in_swapped_edges(tupel_list, integer):
    for idx, tupel in enumerate(tupel_list):
        if integer in tupel:
            position = tupel.index(integer)
            return (tupel, position)
    return None


# Function for drawing the graph
def draw_graph(G, title):
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G, "label")
    node_labels = nx.get_node_attributes(G, "label")

    plt.figure(figsize=(10, 7))
    nx.draw(
        G,
        pos,
        with_labels=True,
        labels=node_labels,
        node_size=700,
        node_color="skyblue",
        font_size=15,
    )
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.title(title)
    plt.show()


def process_graphs(G1, G2, Left_text, Right_text, correction):
    """
    Finds nodes in G1 with the label 'O' that have only one neighbor,
    and checks the edge labels in G1 and G2.
    Depending on the label difference, new nodes are added to both graphs.
    """
    node_num = len(G1.nodes()) + 1
    for node in G1.nodes(data=True):
        # Check whether the node has the label 'O' and has only one neighbor
        if node[1].get("label") == "O" and len(list(G1.neighbors(node[0]))) == 1:
            neighbor = list(G1.neighbors(node[0]))[0]
            # Check whether the edge has a label and exists
            if G1.has_edge(node[0], neighbor) and G2.has_edge(node[0], neighbor):
                label_G1 = G1.edges[node[0], neighbor].get("label")
                label_G2 = G2.edges[node[0], neighbor].get("label")
                # If the labels are the same, do nothing
                if label_G1 == label_G2:
                    continue
                # If G1 has the edge with '-' and G2 with '=', add node 'H' in G1 and 'H+' in G2
                elif label_G1 == "-" and label_G2 == "=":
                    Left_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                    Left_text.append(
                        "  edge [ source "
                        + str(node_num)
                        + " target "
                        + str(node[0])
                        + ' label "-" ]\n'
                    )
                    Right_text.append(
                        "  node [ id " + str(node_num) + ' label "H+" ]\n'
                    )
                    node_num = node_num + 1
                # If G1 has the edge with '=' and G2 with '-', add node 'H+' in G1 and 'H' in G2
                elif label_G1 == "=" and label_G2 == "-":
                    if (
                        G2.nodes[node[0]]["label"] == "O"
                        and len(list(G2.neighbors(node[0]))) == 1
                    ):
                        Right_text.append(
                            "  node [ id " + str(node_num) + ' label "H" ]\n'
                        )
                        Right_text.append(
                            "  edge [ source "
                            + str(node_num)
                            + " target "
                            + str(node[0])
                            + ' label "-" ]\n'
                        )
                        Left_text.append(
                            "  node [ id " + str(node_num) + ' label "H+" ]\n'
                        )
                        node_num = node_num + 1
                correction = True
            elif G1.edges[node[0], neighbor].get("label") != "=":
                Left_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                Left_text.append(
                    "  edge [ source "
                    + str(node_num)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                Right_text.append("  node [ id " + str(node_num) + ' label "H+" ]\n')
                node_num = node_num + 1
                correction = True
        if node[1].get("label") == "O" and len(list(G2.neighbors(node[0]))) == 1:
            neighbor = list(G2.neighbors(node[0]))[0]
            # Check whether the edge has a label and exists
            if (
                G2.edges[node[0], neighbor].get("label") != "="
                and not G1.has_edge(node[0], neighbor)
                and not G2.has_edge(node[0], neighbor)
                or not len(list(G1.neighbors(node[0]))) == 0
            ):
                Right_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                Right_text.append(
                    "  edge [ source "
                    + str(node_num)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                Left_text.append("  node [ id " + str(node_num) + ' label "H+" ]\n')
                node_num = node_num + 1
                correction = True
        # check Case of Water Left Side
        if node[1].get("label") == "O" and len(list(G1.neighbors(node[0]))) == 0:
            if len(list(G2.neighbors(node[0]))) == 1:
                neighbor = list(G2.neighbors(node[0]))[0]
                Left_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                Left_text.append("  node [ id " + str(node_num + 1) + ' label "H" ]\n')
                Left_text.append(
                    "  edge [ source "
                    + str(node_num)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                Left_text.append(
                    "  edge [ source "
                    + str(node_num + 1)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                if G2.edges[node[0], neighbor].get("label") == "=":
                    Right_text.append(
                        "  node [ id " + str(node_num) + ' label "H+" ]\n'
                    )
                    Right_text.append(
                        "  node [ id " + str(node_num + 1) + ' label "H+" ]\n'
                    )
                    node_num += 2
                    correction = True
                if G2.edges[node[0], neighbor].get("label") == "-":
                    Right_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                    Right_text.append(
                        "  edge [ source "
                        + str(node_num)
                        + " target "
                        + str(node[0])
                        + ' label "-" ]\n'
                    )
                    Right_text.append(
                        "  node [ id " + str(node_num + 1) + ' label "H+" ]\n'
                    )
                    node_num += 2
                    correction = True
        # check Case of Water Left Side
        if node[1].get("label") == "O" and len(list(G2.neighbors(node[0]))) == 0:
            if len(list(G1.neighbors(node[0]))) == 1:
                neighbor = list(G1.neighbors(node[0]))[0]
                Left_text.append("  node [ id " + str(node_num) + ' label "H+" ]\n')
                Left_text.append("  node [ id " + str(node_num + 1) + ' label "H+" ]\n')
                Right_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                Right_text.append("  node [ id " + str(node_num + 1) + ' label "H" ]\n')
                Right_text.append(
                    "  edge [ source "
                    + str(node_num)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                Right_text.append(
                    "  edge [ source "
                    + str(node_num + 1)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                node_num += 2
                correction = True

        # check case ammonia
        if node[1].get("label") == "N" and len(list(G1.neighbors(node[0]))) == 1:
            neighbor_left = list(G1.neighbors(node[0]))[0]
            neighbor_right = list(G2.neighbors(node[0]))
            if len(list(G2.neighbors(node[0]))) == 0:
                if G1.edges[node[0], neighbor_left].get("label") == "=":
                    Right_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                    Right_text.append(
                        "  node [ id " + str(node_num + 1) + ' label "H" ]\n'
                    )
                    Right_text.append(
                        "  edge [ source "
                        + str(node_num)
                        + " target "
                        + str(node[0])
                        + ' label "-" ]\n'
                    )
                    Right_text.append(
                        "  edge [ source "
                        + str(node_num + 1)
                        + " target "
                        + str(node[0])
                        + ' label "-" ]\n'
                    )
                    node_num = node_num + 2
                    correction = True
                if G1.edges[node[0], neighbor_left].get("label") == "-":
                    Right_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                    Right_text.append(
                        "  edge [ source "
                        + str(node_num)
                        + " target "
                        + str(node[0])
                        + ' label "-" ]\n'
                    )
                    node_num = node_num + 1
                    correction = True

        # check phosphate case: if O is bonded at two atoms and got split to one of them -> Add a Hydrogen
        if (
            node[1].get("label") == "O"
            and len(list(G2.neighbors(node[0]))) == 1
            and len(list(G1.neighbors(node[0]))) == 2
        ):
            neighbor = list(G2.neighbors(node[0]))[0]
            if G2.edges[node[0], neighbor].get("label") == "-":
                Left_text.append("  node [ id " + str(node_num) + ' label "H+" ]\n')
                Right_text.append("  node [ id " + str(node_num) + ' label "H" ]\n')
                Right_text.append(
                    "  edge [ source "
                    + str(node_num)
                    + " target "
                    + str(node[0])
                    + ' label "-" ]\n'
                )
                node_num = node_num + 1
                correction = True
    # Check if the other side as degree = 2
    # if yes add Hydrogen to the Right side
    # if not do nothing

    return Left_text, Right_text, correction


def correct_GMLs(gml_content):
    correction = False
    # Create the graphs based on the GML file
    G_left, G_right, rule, constrains = parse_gml_and_create_graphs(gml_content)

    swapped_edge = find_swapped_neighbors(G_left, G_right)
    Context_text = []
    Left_text = []
    Right_text = []

    # write Nodes
    for node in G_left.nodes():
        attrs = G_left.nodes[node]
        if attrs["context"] == True:
            Context_text.append(
                "  node [ id " + str(node) + ' label "' + attrs["label"] + '" ]\n'
            )
        else:
            if attrs["rule"] == "left":
                Left_text.append(
                    "  node [ id " + str(node) + ' label "' + attrs["label"] + '" ]\n'
                )

    for node in G_right.nodes():
        attrs = G_right.nodes[node]
        if attrs["context"] != True:
            Right_text.append(
                "  node [ id " + str(node) + ' label "' + attrs["label"] + '" ]\n'
            )

    # write Edges
    for edge in G_left.edges():
        attrs = G_left.edges[edge]
        if attrs["context"] == True:
            Context_text.append(
                "  edge [ source "
                + str(edge[0])
                + " target "
                + str(edge[1])
                + ' label "'
                + attrs["label"]
                + '" ]\n'
            )
        else:
            Left_text.append(
                "  edge [ source "
                + str(edge[0])
                + " target "
                + str(edge[1])
                + ' label "'
                + attrs["label"]
                + '" ]\n'
            )

    for edge in G_right.edges():
        edge_s0 = False
        edge_s1 = False
        attrs = G_right.edges[edge]
        if attrs["context"] != True:
            swapped = find_Node_in_swapped_edges(swapped_edge, edge[0])
            if swapped != None:
                edge_s0 = swapped[0][abs(swapped[1] - 1)]
            swapped = find_Node_in_swapped_edges(swapped_edge, edge[1])
            if swapped != None:
                edge_s1 = swapped[0][abs(swapped[1] - 1)]
            if edge_s0 == False and edge_s1 != False:
                correction = True
                Right_text.append(
                    "  edge [ source "
                    + str(edge[0])
                    + " target "
                    + str(edge_s1)
                    + ' label "'
                    + attrs["label"]
                    + '" ]\n'
                )
            elif edge_s0 != False and edge_s1 == False:
                correction = True
                Right_text.append(
                    "  edge [ source "
                    + str(edge_s0)
                    + " target "
                    + str(edge[1])
                    + ' label "'
                    + attrs["label"]
                    + '" ]\n'
                )
            else:
                Right_text.append(
                    "  edge [ source "
                    + str(edge[0])
                    + " target "
                    + str(edge[1])
                    + ' label "'
                    + attrs["label"]
                    + '" ]\n'
                )

    # correction of alcohol problematic
    Left_text, Right_text, correction = process_graphs(
        G_left, G_right, Left_text, Right_text, correction
    )

    if correction == True:
        file_name = gml_content.replace(".gml", "_corr.gml")
        with open(file_name, "w") as out:
            out.write("rule [\n")
            out.write(' ruleID "' + rule + '_corr"\n')
            out.write(" left [\n")
            for i in Left_text:
                out.write(i)
            out.write(" ]\n")
            out.write(" context [\n")
            for i in Context_text:
                out.write(i)
            out.write(" ]\n")
            out.write(" right [\n")
            for i in Right_text:
                out.write(i)
            out.write(" ]\n")
            if constrains == True:
                out.write("constrainLabelAny [\n")
                out.write('    label "Rest(_R)"\n')
                out.write(
                    '    labels [label "Rest(C)" label "Rest(N)" label "Rest(O)" label "Rest(P)" label "Rest(S)"]\n'
                )
                out.write("    ]\n")
            out.write("]\n")


def process_all_subfolders(root_directory):
    """Runs through all subfolders of a directory and executes the duplicate script in each one."""
    file_paths = []
    for dirpath, dirnames, filenames in os.walk(root_directory):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            file_paths.append(file_path)
    for subdir_path in file_paths:

        if not "_corr" in subdir_path:
            try:
                correct_GMLs(subdir_path)
            except:
                print("correction faild")


path = sys.argv[1]
process_all_subfolders(path)
