import networkx as nx
from itertools import permutations, combinations
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import BamHI, EcoRI, BsmBI # BamHI cuts GGATCC /// EcoRI cuts CTTAAG
import io


def graph_assembly(list_seqs: list):
    '''
    graph_assembly creates a graph with all the possibilities of assembly between the sequences on the list

    Parameters
    ----------
    list_seqs : list
        list of the sorted sequences with sticky ends for assembly (order is important)

    Returns
    -------
    _type_
        graph
    '''
    G = nx.DiGraph()
    
    for comb in permutations(list_seqs, 2):
        seq1, seq2 = comb
        try: 
            seq1 + seq2
        except:
            continue
        else:
            node1 = list_seqs.index(seq1)+1
            node2 = list_seqs.index(seq2)+1
            G.add_node(node1, dseq=seq1)
            G.add_node(node2, dseq=seq2)
            G.add_edge(node1, node2)
        
    if G.nodes:
        return G
    else:
        raise ValueError('Assembly not possible with these sequences')
    

def find_all_paths(graph):
    '''
    find_all_paths is a function that returns a list of all the paths in a graph

    Parameters
    ----------
    graph : _type_
        a graph with multiple nodes and edges (it can be directed or not)

    Returns
    -------
    all_paths : list
        a list of all the possible paths in the graph
    '''
    def find_paths(start, end, path=[]):
        '''
        find_paths  is a secondary function for recursive use in all the nodes of the graph

        Parameters
        ----------
        start : _type_
            node of a graph to start the path
        end : _type_
            node of a graph that ends the path
        path : list, optional
            list of the current path when searching, by default []

        Returns
        -------
        paths : list
            list of all the nodes included in the current path
        '''
        path = path + [start]
        if start == end and len(path) > 1:
            return [path]
        if start not in graph:
            return []
        paths = []
        for node in graph[start]:
            if node not in path or (node == end): 
                newpaths = find_paths(node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

    all_paths = []
    for start in graph:
        for end in graph:
            if start != end:
                paths = find_paths(start, end)
                all_paths.extend(paths)

        # Add circular paths
        circular_paths = find_paths(start, start)
        all_paths.extend(circular_paths)

    return all_paths


def find_paths_seqs(paths, grafo):
    '''
    find_paths_seqs is a function that returns a dictionary, the keys are all the paths in a graph and the values are the sequences assembled (for each path) 

    Parameters
    ----------
    paths : list
        list of lists of paths (sequences that can be assembled by that specific order)
    grafo : _type_
        graph that contains all the nodes and edges (sequences and possible assemblies between them)

    Returns
    -------
    sequencias : dict
        a dictionary, the keys are all the paths in a graph and the values are the sequences assembled (for each path) 
    '''
    sequencias = {}
    for path in paths:
        soma = None
        circular = False # bool
        if len(set(path)) != len(path):
            circular = True
        for i, seq in enumerate(list(set(path))): # se for circular não se pode repetir para juntar as sequencias 
            if i > len(list(set(path)))-1: # quando chega ao fim da lista
                break
            if soma is None: # quando está no inicio do path
                soma = grafo.nodes[seq]['dseq']
            else: # continua no path e ir somando à sequencia existente
                soma += grafo.nodes[seq]['dseq']

        if circular:
            soma = Dseqrecord(soma, circular=True)

        sequencias[(tuple(path))] = soma

    # complementar com calculadora seguid para eliminar sequencias repetidas (alguns paths acabam por ser iguais, ciruclares) 
    # se seguid não der:
    repetidos = [] # retirar sequencias iguais
    for comb in combinations(sequencias.items(), 2):
        if comb[0][1] == comb[1][1] and comb[1] not in repetidos:
            repetidos.append(comb[1])
    for i in repetidos:
        sequencias.pop(i[0])

    return sequencias


def to_graphviz(g):
    '''
    to_graphviz is a secondary function that prints code for using in the website GraphvizOnline for better visualization of the graph. Link: https://dreampuf.github.io/GraphvizOnline/ 

    Parameters
    ----------
    g : _type_
        graph (it can be a object of NetworkX graph)

    Returns
    -------
    _type_
        Code printed for copy+paste in the website GraphvizOnline
    '''
    with io.StringIO() as F:
        print("""digraph{
            node[shape = "circle", style = "filled"];
            """, file = F)
        for u in g:
            for v in g[u]:
                if g[u][v] is not {}:
                    print(f'{u} -> {v}[label = "{g[u][v]}"];', file = F)
                    # print(f'{u} -> {v}[label = ""];', file = F) # backup
                    # print(f"{u} -> {v};", file = F)
                else:
                    print(f'{u} -> {v}[label = "{g[u][v]}"];', file = F)
        print("}", file = F)
        return F.getvalue()
        


if __name__ == '__main__':
    print()

    a = Dseqrecord("CTTAAGatgccctaaccccGAATTC")
    b = Dseqrecord("GAATTCatgccctaaccccGAATTC")
    c = Dseqrecord("GAATTCatgcccgggggggggggccCTTAAG")

    s1,a2  = a.cut(EcoRI)
    b1,s2,b3 = b.cut(EcoRI)
    c1,s3 = c.cut(EcoRI)
    
    lista_seqs = [s1,s2,s3]

    grafo = graph_assembly(lista_seqs)
    # print(grafo)
    # print('Nodes: ', grafo.nodes)
    # print('Edges: ', grafo.edges)

    # print()
    # print(to_graphviz(grafo))

    paths = find_all_paths(grafo)
    # print()
    # for path in find_all_paths(grafo):
    #     print(path)
    # print()

    dic_paths = find_paths_seqs(paths, grafo)
    for i,x in dic_paths.items():
        print(i)
        print(x.figure())

    