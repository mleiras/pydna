import networkx as nx
from itertools import permutations, combinations
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import *
import io
from pydna.seqrecord import SeqRecord
from pydna.amplify import pcr
from Bio.Restriction import *
import goldenhinges
from goldenhinges import OverhangsSelector

# selector = OverhangsSelector(
#     gc_min=0.25,
#     gc_max=0.5,
#     differences=2,
#     forbidden_overhangs=['TTCC'] #, 'ACTC']
# )
# overhangs = selector.generate_overhangs_set(n_overhangs=1)
# print (overhangs)


rb_iis = RestrictionBatch([AcuI, AlwI, BaeI, BbsI, BbvI, BccI, BceAI, BcgI, BciVI, BcoDI, BfuAI, BmrI, BpmI, BpuEI, BsaI, BsaXI, BseRI, BsgI, BsmAI, BsmBI, BsmFI, BsmI, BspCNI, BspMI, BspQI, BsrDI, BsrI, BtgZI, BtsCI, BtsI, BtsIMutI, CspCI, EarI, EciI, Esp3I, FauI, FokI, HgaI, HphI, HpyAV, MboII, MlyI, MmeI, MnlI, NmeAIII, PaqCI, PleI, SapI, SfaNI])


# DESIGN
# STEP 1 -> list all sticky ends...

def list_sticky_ends(seqs):
    sticky_ends_dict = {} # com dicionario dá para ver sticky ends de cada sequencia (ou blunt)
    sticky_ends = [] # só para registar diferentes tipos de stiky ends (sem repetir e não por ordem, nem blunts)
    for sequence in seqs:
        end_5 = sequence.seq.five_prime_end()
        end_3 = sequence.seq.three_prime_end()
        sticky_ends_dict[seqs.index(sequence)] = []
        if end_5[0] == 'blunt':
            sticky_ends_dict[seqs.index(sequence)].append((end_5[0])) # blunt
        else:
            sticky_ends_dict[seqs.index(sequence)].append((end_5[1])) # overhang
            sticky_ends.append(end_5[1].upper()) # será preciso saber a ponta 5 ou 3?
        if end_3[0] == 'blunt':
            sticky_ends_dict[seqs.index(sequence)].append((end_3[0])) # blunt
        else:
            sticky_ends_dict[seqs.index(sequence)].append((end_3[1])) # overhang
            sticky_ends.append(end_3[1].upper())

    sticky_ends = list(set(sticky_ends))

    return sticky_ends_dict, sticky_ends


# STEP 2
# compatible enzymes (enzymes that do not cut within the sequences)

def compatible_enzyme(seqs, enzs):
    
    comp_enzymes = enzs.copy()
    
    for e in enzs:
        if e not in rb_iis: # only type IIs restriction enzymes
            comp_enzymes.remove(e)
        else:
            for seq in seqs:
                if e.search(seq.seq) != []:
                    comp_enzymes.remove(e)
                    break

    if comp_enzymes == []:
        raise ValueError('''No type IIs restriction enzymes available on this list
                        that do not cut within the sequences''')
    
    return comp_enzymes


# STEP 3 

# design relevant primer if needed (add enzyme and comp sticky end)

def design_seqs(seqs, sticky_ends_dict, compatible_enzymes, circular, sticky_ends):

    def design_primers(sequence, compatible_enzymes, f_ovhg, r_ovhg, ind):
        golden_hinges = OverhangsSelector(
                    gc_min=0.25,
                    gc_max=0.5,
                    differences=2,
                    forbidden_overhangs= sticky_ends)
        
        overhangs = golden_hinges.generate_overhangs_set(n_overhangs=2)
        
        if f_ovhg == 'blunt':
            overhang_f = overhangs[0]
            sticky_ends.append(overhang_f)
        else:
            overhang_f = SeqRecord(f_ovhg).reverse_complement().seq
        if r_ovhg == 'blunt':
            overhang_r = overhangs[1]
            sticky_ends.append(overhang_r)   
        else:
            overhang_r = SeqRecord(r_ovhg).reverse_complement().seq

        fp = compatible_enzymes[0].site + 'a' + overhang_f+ sequence.forward_primer
        rp = compatible_enzymes[0].site + 'a' + overhang_r + sequence.reverse_primer

        pcr_prod = pcr(fp, rp, sequence)
        new_seq = Dseqrecord(pcr_prod)

        sticky_ends_dict[ind][0] = overhang_f
        sticky_ends_dict[ind][1] = overhang_r

        return new_seq

    new_seqs = []

    for ind in range(len(seqs)): # index da sequencia na lista (até ao penultimo)

        if ind == 0:
            if circular:
                previous_sticky_end = sticky_ends_dict[len(seqs)-1][1]
            else:
                previous_sticky_end = 'blunt' # temporario -> mudar para None
        else:
            previous_sticky_end = sticky_ends_dict[ind-1][1]

        if ind == len(seqs)-1:
            if circular:
                next_sticky_end = sticky_ends_dict[0][0]
            else:
                next_sticky_end = 'blunt' # temporario -> mudar para None
        else:
            next_sticky_end = sticky_ends_dict[ind+1][0] 


        if sticky_ends_dict[ind][0] == 'blunt' and sticky_ends_dict[ind][1] == 'blunt':
            new_seq = design_primers(seqs[ind], compatible_enzymes = compatible_enzymes, f_ovhg = previous_sticky_end, r_ovhg = next_sticky_end, ind=ind)

            new_seqs.append(new_seq)

        elif sticky_ends_dict[ind][0] != 'blunt' and sticky_ends_dict[ind][1] != 'blunt':
            new_seqs.append(seqs[ind]) # não fazer nada porque já tem sticky ends (não é preciso enzima porque o design dos outros fragmentos vão coincidir)
    
        else:
            raise ValueError('ERRO DESIGN PRIMER')

    return new_seqs
    


## ASSEMBLY
# STEP 1

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
    
# STEP 2

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

# STEP 3

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
            path.pop()
        for i, seq in enumerate(path): # se for circular não se pode repetir para juntar as sequencias 
            if i >= len(list(set(path))): # quando chega ao fim da lista
                break
            if soma is None: # quando está no inicio do path
                soma = grafo.nodes[seq]['dseq']
            else: # continua no path e ir somando à sequencia existente
                soma += grafo.nodes[seq]['dseq']

        if circular:
            soma = Dseqrecord(soma, circular=True)

        sequencias[(tuple(path))] = soma

    # complementar com calculadora seguid para eliminar sequencias repetidas (alguns paths acabam por ser iguais, ciruclares) 
    #### AQUI: retorna 4 seqs (sequencias unicas)
    dicio = sequencias.copy()

    for comb in combinations(sequencias.items(), 2):
        if comb[0][1].circular and comb[1][1].circular:
            if comb[0][1].cseguid() == comb[1][1].cseguid():
                if comb[0][0] in dicio:
                    dicio.pop(comb[0][0])
        elif comb[0][1].useguid() == comb[1][1].useguid(): 
            if comb[0][0] in dicio:
                dicio.pop(comb[0][0])
    

    return dicio




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

