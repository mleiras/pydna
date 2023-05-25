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
    G = nx.Graph()
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
                # print(index, i, '1try')
                G.add_node(sequence.figure()) # em vez de figure tem de ser a sequencia apenas, algo do genero?
                G.add_node(seq.figure())
                G.add_edge(sequence.figure(),seq.figure())
                # print('new_seq', new_seq.figure())
                # print()
            except:
                pass


    
        # if seq + list_seqs[index+1]:
        #     new_seq = seq + list_seqs[index+1]
        #     print(index, new_seq.figure())

    print(G)

    return G


def to_graphviz(g):
    with io.StringIO() as F:
        print("""digraph{
            node[shape = "circle", style = "filled"];
            """, file = F)
        for u in g:
            for v in g[u]:
                if g[u][v] is not {}:
                    print(f'{u} -> {v}[label = ""];', file = F)
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


    # só funciona com Dseq (em vez de Dseqrecord)

    # print('5', s1.five_prime_end())
    # print('3', s1.three_prime_end())
    # print()

    # print('5', s2.five_prime_end())
    # print('3', s2.three_prime_end())
    # print()

    # print('5', s3.five_prime_end())
    # print('3', s3.three_prime_end())


