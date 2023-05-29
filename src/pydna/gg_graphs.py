import networkx as nx
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio.Restriction import BamHI, EcoRI, BsmBI # BamHI cuts GGATCC /// EcoRI cuts CTTAAG
import io


def graph_assembly(list_seqs: list):
    '''
    graph_assembly creates a graph with all the possibilities of assembly between the sequences on the list

    Parameters
    ----------
    list_seqs : list
        list of sequences with sticky ends for assembly

    Returns
    -------
    _type_
        graph
    '''
    G = nx.DiGraph()

    for index, sequence in enumerate(list_seqs):
        lista = list_seqs.copy()
        lista.remove(sequence)
        for i, seq in enumerate(lista):
            try: 
                new_seq = sequence + seq
                node1 = index+1 #,sequence.figure().split('\n')[0][12:14]) # index da sequencia (na lista original) + tamanho da sequencia para referencia
                node2 = list_seqs.index(seq)+1#,seq.figure().split('\n')[0][12:14])
                G.add_node(node1, size = sequence.figure().split('\n')[0][12:14], dseq=sequence)
                G.add_node(node2, size = seq.figure().split('\n')[0][12:14], dseq=seq)
                G.add_edge(node1,node2, w=new_seq)
            except:
                pass

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
        if start == end and len(path) > 1: #dif
            return [path]
        if start not in graph:
            return []
        paths = []
        for node in graph[start]:
            if node not in path or (node == end): #dif
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
    sequencias = {}
    for path in paths:
        # print(path)
        soma = None
        circular = False
        if len(set(path)) != len(path):
            circular = True
        for i, seq in enumerate(list(set(path))):
            
            if i > len(list(set(path)))-1:
                break
            if soma is None:
                soma = grafo.nodes[seq]['dseq']
                # print('soma inicial',i, soma)
            else:
                # print()
                # print('i',i,  grafo.nodes[seq])
                # print('i', grafo.nodes[seq+1])
                soma += grafo.nodes[seq]['dseq']
        # print('SOMA', soma.figure())
        # print(soma)
        if circular:
            soma = Dseqrecord(soma, circular=True)

        sequencias[(tuple(path))] = soma

    # print(sequencias)
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
        
# pprint.pprint(g)


if __name__ == '__main__':

    a = Dseqrecord("CTTAAGatgccctaaccccGAATTC")
    b = Dseqrecord("GAATTCatgccctaaccccGAATTC")
    c = Dseqrecord("GAATTCatgcccgggggggggggccCTTAAG")
    # print(a.figure())
    # print(b.figure())
    # print(c.figure()) 

    s1,a2  = a.cut(EcoRI)
    b1,s2,b3 = b.cut(EcoRI)
    c1,s3 = c.cut(EcoRI)
    # d = Dseqrecord("CTTAAGatgccctttttttGAATTC")
    # s4,d2 = d.cut(EcoRI)

    # só funciona com Dseq (em vez de Dseqrecord):
    # print('5', s1.five_prime_end())
    # print('3', s1.three_prime_end())
    # print()

    # print(s1.figure())
    # print()
    # print(s2.figure())
    # print()
    # print(s3.figure())
    # print()

    # d = s1+s2+s3
    # e = s1+s3
    # print(d.figure())
    # print(e.figure())
    # print(repr(s1))
    # print()
    # print(repr(s2))
    # print()
    # print(repr(s3))
    # print()
    
    # lista_seqs = [s1,s2,s3]
    lista_seqs = [s1,s1]

    grafo = graph_assembly(lista_seqs)
    # print(grafo)
    # print('Nodes: ', grafo.nodes)
    # print('Edges: ', grafo.edges)

    # print()

    # print()
    # print(to_graphviz(grafo))

    # print('All the paths on the graph: ')
    paths = find_all_paths(grafo)
    # print()
    # for path in find_all_paths(grafo):
    #     print(path)
    # print()

    dic_paths = find_paths_seqs(paths, grafo)
    for i,x in dic_paths.items():
        print(i)
        print(x.figure())

    # print(unique_paths(paths, grafo))
    # print(grafo.nodes[0]['dseq']) # forma de aceder à informação toda da sequencia em questão

