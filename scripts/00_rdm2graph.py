import networkx as nx
import sys
import re 
import matplotlib.pyplot as plt

not_work = []  # List of reactions which could not be built up
multiple_solutions = [] # List of reaction

# Create lists of the R-atoms, D-atoms, M-atoms for both reaction sides
def parse_equation(eq):
    left_r = []
    right_r = []
    left_d = []
    right_d = []
    left_m = [] 
    right_m = []

    eq_parts = eq.split(':')

    for parts in range(len(eq_parts)):
        side_parts = eq_parts[parts].split('-')
        left_side = side_parts[0].split('+')
        right_side = side_parts[1].split('+')
        if parts == 0:
            for i in left_side:
                left_r.append(i)
            for i in right_side:
                right_r.append(i) 
        if parts == 1:
            for i in left_side:
                if i is not '*':
                    left_d.append(i)
            for i in right_side:
                if i is not '*':
                    right_d.append(i)
        if parts == 2:
            for i in left_side:
                if i is not '*':
                    left_m.append(i)
            for i in right_side:
                if i is not '*':
                    right_m.append(i) 

    return left_r, right_r, left_d, right_d, left_m, right_m 

      
def createNewLeaf(tree,node_list,atom_type,node_num,attr,side,parent):
    for child in node_list:
        attr[node_num] = {'label': child,'atomtype': atom_type,'map': node_num,'side':side}
        tree.add_node(node_num)
        tree.add_edge(parent, node_num)
        node_num = node_num + 1

    return tree, node_num, attr

undefindedAtoms = ['Z','C0','N0','O0','S0','X']
count_correct = 0
count_undefinedAtoms = 0
rc_list = []

# Read file
path = sys.argv[1]
with open(path, 'r') as f:
    lines = f.readlines()
for line in lines:
    if line.startswith('INFO'):
        continue
        
    equations = []
    equation = line.strip().split('\t')
    
    # check if there are undefined Atoms, if yes than skip
    check_atoms = True
    for substring in undefindedAtoms:
        if substring in equation[1]:
            count_undefinedAtoms = count_undefinedAtoms + 1         
            check_atoms = False
    if check_atoms == False: 
        continue
    
    rxn = equation[0]
    equations = equation[1].split(' ')
    
    graphs_left = []
    graphs_right = []
    node_num = 0  
    # create for every R-Atom a single Graph
    for eq in equations:

            # separete the R-Atoms form the M-Atoms and the D-Atoms
            left_r, right_r, left_d, right_d, left_m, right_m = parse_equation(eq)                    
            print('START ' + rxn + ' : ' + eq)
            tree_left = nx.Graph() 
            tree_right = nx.Graph() 
                       
            # set roots           
            parent_left = node_num
            attr_l = {node_num:{'label': left_r[0],'atomtype': 'r','map': node_num,'side':'left'}}
            tree_left.add_node(parent_left)
            node_num = node_num+1
             
            parent_right = node_num
            attr_r = {node_num:{'label': right_r[0],'atomtype': 'r','map': node_num,'side':'right'}}
            tree_right.add_node(parent_right)
            node_num = node_num+1 

            # Create leafs for D-atoms
            tree_left, node_num, attr_l = createNewLeaf(tree_left,left_d,'d',node_num, attr_l,'left',parent_left)
            tree_right, node_num, attr_r = createNewLeaf(tree_right,right_d,'d',node_num, attr_r,'right',parent_right)
            
            # Create leafs for M-atoms
            tree_left, node_num, attr_l = createNewLeaf(tree_left,left_m,'m',node_num, attr_l,'left',parent_left)
            tree_right, node_num, attr_r = createNewLeaf(tree_right,right_m,'m',node_num, attr_r,'right',parent_right)
              
            # Set Attributes
            nx.set_node_attributes(tree_left, attr_l)
            nx.set_node_attributes(tree_right, attr_r)

            # For writing file, but all created graphs in a list
            graphs_left.append(tree_left.copy())
            graphs_right.append(tree_right.copy())

    # Create files with RDM patterns as graphs       
    multi_graph = nx.MultiGraph() 
    for left_graph in graphs_left:         
        mapping1 = {}  
        for u, v in left_graph.edges():
            if u not in mapping1:
                mapping1[u] = u + 1000
            if v not in mapping1:
                mapping1[v] = v + 1000
            multi_graph.add_edge(mapping1[u], mapping1[v])                
            attr = {}
            for node in left_graph.nodes(data=True):    
                attr[node[0]+1000] = node[1]
            nx.set_node_attributes(multi_graph, attr)     
    for right_graph in graphs_right:
            multi_graph.add_edges_from(right_graph.edges(data=True))
            multi_graph.add_nodes_from(right_graph.nodes(data=True))    
                 
    nx.write_gexf(multi_graph, './00_RCLASS_Graphs/' + str(rxn) + '.gexf')  
    count_correct = count_correct +1        

print('#### Number of Generated RCLASS Graphs: ',count_correct)
print('#### Number of RCLASS with undefined Atoms: ',count_undefinedAtoms)


