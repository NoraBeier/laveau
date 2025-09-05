import mod
import glob
import os
import re
import regex
import itertools
import requests
import networkx as nx
import copy
from mod import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import openbabel


config.rule.ignoreConstraintsDuringInversion = True


# Load compound mol format form KEGG Database
def get_compound_mol(compound_id):
    url = f"http://rest.kegg.jp/get/compound:{compound_id}/mol"
    try:
        response = requests.get(url)
        response.raise_for_status()  # This throws an exception if the status code is not 200 (OK)
        return response.text
    except requests.exceptions.RequestException as e:
        return None


def mol_to_smiles(molfile_content):
    if isinstance(molfile_content, bytes):
        molfile_content = molfile_content.decode("utf-8")
    elif not isinstance(molfile_content, str):
        return Falseprint

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "smi")
    mol = openbabel.OBMol()
    success = obConversion.ReadString(mol, molfile_content)

    if not success:
        return False

    smiles = obConversion.WriteString(mol).strip()
    if smiles != "[H+]":
        smiles_without_charges = re.sub(
            r"\[\w+[\+\-]\d*\]",
            lambda x: x.group(0).replace("+", "").replace("-", ""),
            smiles,
        )
    else:
        return smiles
    return smiles_without_charges


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


def count_protons(text):
    # Define patterns for the left block and right block
    left_block_pattern = r"left\s*\[(.*?)\]\s*\]"
    right_block_pattern = r"right\s*\[(.*?)\]\s*\]"

    # Extract the left block and right block
    left_block = re.search(left_block_pattern, text, re.DOTALL)
    right_block = re.search(right_block_pattern, text, re.DOTALL)

    # Count the number of 'H+' in the blocks
    left_h_plus_count = (
        len(re.findall(r"label\s*\"H\+\"", left_block.group(1))) if left_block else 0
    )
    right_h_plus_count = (
        len(re.findall(r"label\s*\"H\+\"", right_block.group(1))) if right_block else 0
    )

    return left_h_plus_count, right_h_plus_count


def getReactionSmiles(dg, correct_edge):
    origSmiles = {}
    # Change the class labels for every atom
    for v in dg.vertices:
        s = v.graph.smilesWithIds
        s = regex.sub(":([0-9]+)\]", ":o\\1]", s)
        origSmiles[v.graph] = s
    res = {}

    # maps = DGVertexMapper(correct_edge,leftLimit=1,rightLimit=1)
    maps = DGVertexMapper(correct_edge, rightLimit=1)
    # Replace the original local IDs with global IDs
    eductSmiles = [origSmiles[g] for g in maps.left]

    for v in maps.left.vertices:
        s = eductSmiles[v.graphIndex]
        s = s.replace(":o%d]" % v.vertex.id, ":%d]" % v.id)
        eductSmiles[v.graphIndex] = s

    strs = set()
    for r in maps:
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
    all_combinations = []
    for r in range(1, len(input_list)):
        combinations_r = list(itertools.combinations(input_list, r))
        all_combinations.extend([list(comb) for comb in combinations_r])
    return all_combinations


def create_Educts_and_Products(reaction_data, path, rule, rule_var, universe):

    educ = []  # List of educts known by the database for finding correct hyperedge
    prod = []  # List of products known by the database for finding correct hyperedge
    subspace = []
    # load all compounants of the reaction
    educts = copy.deepcopy(reaction_data[rule][0])
    products = copy.deepcopy(reaction_data[rule][1])

    # Define all educts
    for e in educts.keys():
        if e == "C00080":
            continue
        while educts[e] > 0:
            if e == False:
                break
            elif isinstance(universe[e], Graph):  # and universe[e].id:
                educ.append(universe[e])
                educts[e] = educts[e] - 1
            else:
                break
        if e != False and isinstance(universe[e], Graph):  # and universe[e].id:
            subspace.append(universe[e])

    # Define all products
    for e in products.keys():
        if e == "C00080":
            continue
        while products[e] > 0:
            if e == False:
                break
            elif isinstance(universe[e], Graph):  # and universe[e].id:
                prod.append(universe[e])
                products[e] = products[e] - 1
            else:
                break
        if e != False:
            iso_check = False
            for g1 in subspace:
                if isinstance(g1, Graph) and g1.id:
                    if isinstance(universe[e], Graph):  # and universe[e].id:
                        if g1.isomorphism(universe[e]) == True:
                            iso_check = True
            if iso_check == False and isinstance(
                universe[e], Graph
            ):  # and universe[e].id:
                subspace.append(universe[e])

    # add protons to Educts and Products
    with open(path + "/" + rule + "/" + rule_var, "r") as f:
        lines = f.readlines()
    string_rule = " ".join(lines)
    left_protons, right_protons = count_protons(string_rule)
    if left_protons != 0 or right_protons != 0:
        h_mol = smiles("[H+]", name="C00080")
        subspace.append(h_mol)
        while left_protons > 0:       
            educ.append(h_mol)
            subspace.append(h_mol)
            left_protons = left_protons - 1
        while right_protons > 0:
            prod.append(h_mol)
            right_protons = right_protons - 1
    return educ, prod, subspace


def find_AAMs(dg, educ, prod, path, found, count_aam, count_aam_list):

    def write_rule(dg, global_edge, outfile, count_aam, count_aam_list, found):
        try:
            # parse graphs to SMILES Format
            res = getReactionSmiles(dg, global_edge)
        except:
            with open(path + "/WildcardCase.txt", mode="a") as outp:
                outp.write(rule + " " + str(s) + "\n")
            with open(path + "/Variations_Wildcard.txt", mode="a") as out:
                out.write(rule_var + "\n")
            count_aam += 1
            count_aam_list.append(rule)
            return True, count_aam, count_aam_list, True
        else:
            assert len(res) != 0
            for s in res:
                if found == False:
                    count_aam += 1
                    count_aam_list.append(rule)
                    found = True
                with open(path + outfile, mode="a") as outp:
                    outp.write(rule + " " + str(s) + "\n")
                with open(path + "/Variations_" + outfile, mode="a") as out:
                    out.write(rule_var + "\n")
            return count_aam, count_aam_list, found

    # search for HyperEdge which creates all the products -> global atom-to-atom maps
    try:
        global_edge = dg.findEdge(educ, prod)
    except:
        global_edge = False
    if global_edge:
        count_aam, count_aam_list, found = write_rule(
            dg, global_edge, "GlobalReactions.txt", count_aam, count_aam_list, found
        )
    else:
        # try with changed directions
        try:
            global_edge = dg.findEdge(prod, educ)
        except:
            global_edge = False
        if global_edge:
            count_aam, count_aam_list, found = write_rule(
                dg, global_edge, "GlobalReactions.txt", count_aam, count_aam_list, found
            )
        else:
            ## Filter cofaktors out
            cofactors = [
                "C00006",  # NADP
                "C00005",  # NADPH
                "C00002",  # ATP
                "C00020",  # AMP
                "C00013",  # Diphosphate
                "C00008",  # ADP
                "C00003",  # NAD+
                "C00004",  # NADH
                "C00016",  # FAD
                "C01352",  # FADH2
                "C00143",  # 5,10-Methylenetetrahydrofolate
                "C00101",  # Tetrahydrofolate
                "C06453",  # Methylcobalamin
                "C05776",  # Vitamin B12
                "C01115",  # L-Galactono-1,4-lactone
                "C00072",  # Vitamin C
                "C00010",  # CoA
                "C00024",  # Acetyl-CoA
            ]

            # Remove Cofactors from the molecule lists
            educ_f = educ
            for v in educ:
                if v.name in cofactors:
                    educ_f = [x for x in educ_f if x.name != v.name]
                if v in prod and v.name != "C00080":
                    prod_f.remove(v)
                    educ_f.remove(v)
            prod_f = prod
            for v in prod:
                if v.name in cofactors:
                    prod_f = [x for x in prod_f if x.name != v.name]
            try:
                global_edge = dg.findEdge(educ_f, prod_f)
            except:
                global_edge = False
            if global_edge:
                count_aam, count_aam_list, found = write_rule(
                    dg,
                    global_edge,
                    "GlobalReactions.txt",
                    count_aam,
                    count_aam_list,
                    found,
                )
            else:
                # try with changed directions
                try:
                    global_edge = dg.findEdge(prod_f, educ_f)
                except:
                    global_edge = False
                if global_edge:
                    count_aam, count_aam_list, found = write_rule(
                        dg,
                        global_edge,
                        "GlobalReactions.txt",
                        count_aam,
                        count_aam_list,
                        found,
                    )

    # check for partial AAMs: educts and products must be a subset of the given molecules
    num_partial = 0
    if global_edge == False:
        for hyperedge in dg.edges:
            if any(v.graph not in educ for v in hyperedge.sources):
                continue
            if any(v.graph not in prod for v in hyperedge.targets):
                continue
            num_partial += 1
            count_aam, count_aam_list, found = write_rule(
                dg,
                global_edge,
                "PartialReactions.txt",
                count_aam,
                count_aam_list,
                found,
            )
        # other direction
        for hyperedge in dg.edges:
            if any(v.graph not in prod for v in hyperedge.sources):
                continue
            if any(v.graph not in educ for v in hyperedge.targets):
                continue
            num_partial += 1
            count_aam, count_aam_list, found = write_rule(
                dg,
                global_edge,
                "PartialReactions.txt",
                count_aam,
                count_aam_list,
                found,
            )

    if num_partial != 0 or found == True:
        return True, count_aam, count_aam_list, True
    else:
        return False, count_aam, count_aam_list, False


def createMoleculeDatabase(reaction_data, universe):
    for molecule in reaction_data:
        for m in molecule.keys():
            s = mol_db[m]
            if s == "FALSE":
                universe[m] = s
            else:
                mol = smiles(s, name=m, allowAbstract=True)
                if mol not in universe:
                    universe[m] = mol

    return universe

# load molecule data form KEGG database
mol_db = {}
with open("./Additional_Files/KEGG_MoleculeDB.txt", "r") as f:
    lines = f.readlines()
for line in lines:
    line = line.split(",")
    sm = line[1].split("\n")
    mol_db[line[0]] = sm[0]

# load reaction data form KEGG database
with open("./Additional_Files/REACTION_RCLASS_DATA.txt", "r") as in_f:
    lines = in_f.readlines()

reaction_data = {}
for line in lines:
    # Extract reaction name
    if line.startswith("Re"):
        line = line.split(":")
        line = line[1].split("\n")
        rn = line[0]
        continue
    # Extract compounds and stochiometric data
    if line.startswith("C"):
        compound_ids_str = re.search(r"Compound IDs:\[(.*?)\]", line).group(1)
        compound_ids = eval(compound_ids_str)
        educts = {}
        products = {}
        delimiter = "<=>"
        left_side = True

        for i in range(len(compound_ids)):
            item = compound_ids[i]
            if item == delimiter:
                left_side = False
                continue
            if item.startswith("C"):
                if left_side:
                    count = 1
                    if i > 0 and compound_ids[i - 1].isdigit():
                        count = int(compound_ids[i - 1])
                    educts[item] = count
                else:
                    count = 1
                    if i > 0 and compound_ids[i - 1].isdigit():
                        count = int(compound_ids[i - 1])
                    products[item] = count
        reaction_data[rn] = [educts, products]

# Statistics
count_aam = 0  # counts the number of generated aams
count_miss = 0  # counts the number of reactions which could not be generated
count_aam_list = []
count_miss_list = []

# Load all Reaction rules witch AAMs should be generated
path = '../demo/03_DPO_Rules/'
path_results = '../demo/04_AAM_results/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

rules_names = os.listdir(path)
universe = {}  # List of all known Molecules by KEGG DB


for rule in rules_names:

    found = False  # marker if a aam could be generated

    # complete universe
    universe = createMoleculeDatabase(reaction_data[rule], universe)
    rule_variations = os.listdir(path + "/" + rule)  # find all rules for one reaction

    for rule_var in rule_variations:
        print(rule_var)
        var_succ = False  # marker if a variation of a rule could generated a aam

        # generated educts and products
        educ, prod, subspace = create_Educts_and_Products(
            reaction_data, path, rule, rule_var, universe
        )

        # Create Network
        dg = DG(
            labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation),
            graphDatabase=educ,
        )
        dg.build().apply(
            subspace,
            ruleGML(
                path + "/" + rule + "/" + rule_var, add=False, printStereoWarnings=False
            ),
            onlyProper=False,
        )

        # find atom-to-atom maps in Network
        found, count_aam, count_aam_list, var_succ = find_AAMs(
            dg, educ, prod, path_results, found, count_aam, count_aam_list
        )

        if found == False:
            # try the opposite direction
            # generated educts and products
            educ, prod, subspace = create_Educts_and_Products(
                reaction_data, path, rule, rule_var, universe
            )

            # Create Network
            dg = DG(
                labelSettings=LabelSettings(
                    LabelType.Term, LabelRelation.Specialisation
                ),
                graphDatabase=educ,
            )
            dg.build().apply(
                subspace,
                ruleGML(
                    path + "/" + rule + "/" + rule_var,
                    invert=True,
                    add=False,
                    printStereoWarnings=False,
                ),
                onlyProper=False,
            )
            correct_edge = False

            # find atom-to-atom maps in Network
            found, count_aam, count_aam_list, var_succ = find_AAMs(
                dg, educ, prod, path_results, found, count_aam, count_aam_list
            )

        if var_succ == False:
            with open(
                path_results + "/FalseRuleVariation.txt", mode="a"
            ) as outp:
                outp.write(rule_var + "\n")

        else:
            with open(path_results + "/TrueRuleVariation.txt", mode="a") as outp:
                outp.write(rule_var + "\n")

    if found == False:
        count_miss += 1
        count_miss_list.append(rule)
        with open(path_results + "/FalseReactions.txt", mode="a") as outp:
            outp.write(rule + "\n")

print("######################## SUMMARY #############################")
print("Total No. of Reactions: ", len(rules_names))
print("No. of generated AAMs: ", count_aam, count_aam_list)
print("No. of RXN which faild: ", count_miss, count_miss_list)
print("#####################################################")
