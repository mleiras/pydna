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
    # print(list_seqs)

    # for index, elem in enumerate(layers):        
    #     if(index<(len(layers)-1)):    
    #         thiseLayer= elem        
    #         nexteLayer= layers[index+1]


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
                G.add_edge(node1,node2)
                # G.add_node((index,sequence.figure().split('\n')[0])) # em vez de figure tem de ser a sequencia apenas, algo do genero?
                # G.add_node((list_seqs.index(seq),seq.figure().split('\n')[0]))
                # G.add_edge(sequence.figure().split('\n')[0],seq.figure().split('\n')[0])
            except:
                pass

    print(G)
    return G


def find_all_paths(graph):
    def find_paths(start, end, path=[]):
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

    for path in all_paths:
        print(path)

    return all_paths


def to_graphviz(g):
    with io.StringIO() as F:
        print(f'\n Code to GraphvizOnline: \n') #(https://dreampuf.github.io/GraphvizOnline/)')
        print("""digraph{
            node[shape = "circle", style = "filled"];
            """, file = F)
        for u in g:
            for v in g[u]:
                if g[u][v] is not {}:
                    print(f'{u} -> {v}[label = ""];', file = F)
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


    # print(s1.figure())
    # print()
    # print(s2.figure())
    # print()
    # print(s3.figure())
    # print()

    d = s1+s2+s3
    e = s1+s3
    # print(d.figure())
    # print(e.figure())
    # print(repr(s1))
    # print()
    # print(repr(s2))
    # print()
    # print(repr(s3))
    # print()
    
    lista_seqs = [s1,s2,s3]


    print()
    print('GRAPH: \n')

    grafo = graph_assembly(lista_seqs)

    print(to_graphviz(grafo))
    find_all_paths(grafo)
    print()
    # print(grafo.nodes[0]['dseq']) # forma de aceder à informação toda da sequencia em questão




    # só funciona com Dseq (em vez de Dseqrecord)

    # print('5', s1.five_prime_end())
    # print('3', s1.three_prime_end())
    # print()

    # print('5', s2.five_prime_end())
    # print('3', s2.three_prime_end())
    # print()

    # print('5', s3.five_prime_end())
    # print('3', s3.three_prime_end())


