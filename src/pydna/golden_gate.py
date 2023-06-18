from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from Bio.Restriction import *
from pydna.seqrecord import SeqRecord
from pydna.amplify import pcr
import goldenhinges
from goldenhinges import OverhangsSelector
import networkx as nx
from itertools import permutations, combinations
from pydna.seqrecord import SeqRecord


rb_iis = RestrictionBatch([AcuI, AlwI, BaeI, BbsI, BbvI, BccI, BceAI, BcgI, BciVI, BcoDI, BfuAI, BmrI, BpmI, BpuEI, BsaI, BsaXI, BseRI, BsgI, BsmAI, BsmBI, BsmFI, BsmI, BspCNI, BspMI, BspQI, BsrDI, BsrI, BtgZI, BtsCI, BtsI, BtsIMutI, CspCI, EarI, EciI, Esp3I, FauI, FokI, HgaI, HphI, HpyAV, MboII, MlyI, MmeI, MnlI, NmeAIII, PaqCI, PleI, SapI, SfaNI])

def _list_sticky_ends_(seqs):
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
            sticky_ends.append(end_5[1].upper()) 
        if end_3[0] == 'blunt':
            sticky_ends_dict[seqs.index(sequence)].append((end_3[0])) # blunt
        else:
            sticky_ends_dict[seqs.index(sequence)].append((end_3[1])) # overhang
            sticky_ends.append(end_3[1].upper())

    sticky_ends = list(set(sticky_ends))

    return sticky_ends_dict, sticky_ends


def _compatible_enzyme_(seqs, enzs):
    
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


def _design_seqs_(seqs, sticky_ends_dict, compatible_enzymes, circular, sticky_ends):

    def _design_primers_(sequence, compatible_enzymes, f_ovhg, r_ovhg, ind):

        if isinstance(sequence, Dseqrecord):
            sequence = primer_design(sequence)

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

        fp = compatible_enzymes[0].site + 'a' + overhang_f + sequence.forward_primer
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
                previous_sticky_end = 'blunt' # Desta forma está a criar primer com sticky end "aleatório" a não corresponder à ultima seq // Algternativa é mudar para None e depois alterar na função em que não cria um primer para essa ponta mas depois a DNA ligase poderá ligar as duas sequencias
        else:
            previous_sticky_end = sticky_ends_dict[ind-1][1]

        if ind == len(seqs)-1:
            if circular:
                next_sticky_end = sticky_ends_dict[0][0]
            else:
                next_sticky_end = 'blunt' # Desta forma está a criar primer com sticky end "aleatório" a não corresponder à primeira seq
        else:
            next_sticky_end = sticky_ends_dict[ind+1][0] 


        if sticky_ends_dict[ind][0] == 'blunt' and sticky_ends_dict[ind][1] == 'blunt':
            new_seq = _design_primers_(seqs[ind], compatible_enzymes = compatible_enzymes, f_ovhg = previous_sticky_end, r_ovhg = next_sticky_end, ind=ind)

            new_seqs.append(new_seq)

        elif sticky_ends_dict[ind][0] != 'blunt' and sticky_ends_dict[ind][1] != 'blunt':
            new_seqs.append(seqs[ind]) # não fazer nada porque já tem sticky ends (não é preciso enzima porque o design dos outros fragmentos vão coincidir)
    
        else:
            raise ValueError('ERRO DESIGN PRIMER')

    return new_seqs


def GoldenGateDesigner(seqs: list, enzs: list, circular: bool = True) -> list:
    '''
    GoldenGateDesigner is a function that receives a list of sequences and a list of enzymes and, using these, it must return a list of sequences with compatible primers or overhangs.

    Parameters
    ----------
    seqs : list
        list of sequences (order is important)
    enzs : list
        list of enzymes (order is important)
    circular : bool, optional 
        To define if intended final sequence is circular or not, by default True

    Returns
    -------
    new_seqs: list
        list with Dseqrecords and amplicons
    '''

    # STEP 1:  list all sticky ends
    sticky_ends_dict,  sticky_ends = _list_sticky_ends_(seqs)
   

    # STEP 2:  compatible enzymes (enzymes that do not cut within the sequences)
    compatible_enzymes = _compatible_enzyme_(seqs, enzs) 

    # STEP 3: for loop over sequences
    # design relevant primer if needed (add enzyme and comp sticky end)
    new_seqs = _design_seqs_(seqs, sticky_ends_dict, compatible_enzymes, circular, sticky_ends)
   
    # STEP 4: Return a list with Dseqrecords and amplicons

    return new_seqs , compatible_enzymes[0] # o user tem de saber qual a enzima compativel




def _graph_assembly_(list_seqs: list):
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
    

def _find_all_paths_(graph):
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
    def _find_paths_(start, end, path=[]):
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
                newpaths = _find_paths_(node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

    all_paths = []
    for start in graph:
        for end in graph:
            if start != end:
                paths = _find_paths_(start, end)
                all_paths.extend(paths)

        # Add circular paths
        circular_paths = _find_paths_(start, start)
        all_paths.extend(circular_paths)

    return all_paths



def _find_paths_seqs_(paths, grafo):
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




def GoldenGateAssembler(seqs: list) -> dict:
    '''
    GoldenGateAssembler is a function that receives a list of sequences (fragments) and tries all the possible assembles between them, returns a dictionary of all the possible sequences ligated (with only DNA ligase).

    Parameters
    ----------
    seqs : list
        list of sequences (fragments of DNA)

    Returns
    -------
    dict_seqs : dict
        dictionary with the possible sequences ligated (keys are the order of the sequences ligated; values are the dseqrecords of the sequences assembled)
    '''

    # STEP 1: creates a graph where each node is a sequence and each edge represents the link between sequences (each possible combination in a list of sequences)
    graph = _graph_assembly_(seqs)
    
    # STEP 2: find all paths (each combination of a possible assembled sequence)  
    paths = _find_all_paths_(graph)

    # STEP 3: using the paths as keys, returns all the assembled sequences for each of them (without repetitions) 
    dict_seqs = _find_paths_seqs_(paths, graph)
    
    return dict_seqs

    

if __name__ == '__main__':

    print()

    from pydna.dseqrecord import Dseqrecord
    from Bio.Restriction import *


    a = Dseqrecord("GGTCTCatgccctaaccccccctaacccacGAGACC")
    b = Dseqrecord("GGTCTCacccagtaaccatgtagtaaaaaaatgccGAGACC")
    c = Dseqrecord("GGTCTCtatgctcgggggggggggcagtacagtaGAGACC")
    d = Dseqrecord("GGTCTCccagttgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatgGAGACC")

    lista_enz = [EcoRI, BsaI, BsmBI, BamHI, BpmI]
    # lista_enz = [EcoRI, BamHI] # not type IIs
    
    a1,a2,a3  = a.cut(BsaI)
    b1,b2,b3 = b.cut(BsaI)
    c1,c2,c3 = c.cut(BsaI)
    d1,d2,d3 = d.cut(BsaI)

    # lista_seqs = [a2,b2,c2,d2]

    # print('Sequencias: \n')
    # for i in lista_seqs:
    #     print(i.figure())

    e = primer_design(Dseqrecord("ccagttgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")) # AMPLICON
    f = primer_design(Dseqrecord("ccagttgctaacccttccttggtcgtttcgttcgaaacttacgatg")) # AMPLICON
    g = primer_design(Dseqrecord("tatgctcgggggggaacaagatcgacgacatttcgttcgaaacttacgagggggcagtacagta")) # AMPLICON

    # lista_seqs2 = [a2,b2,c2,d2,e]
    # lista_seqs3 = [e,f, a2,b2,c2,d2,g]
    lista_seqs3 = [e,f,g]

    seq1 = Dseqrecord(
    "ccagttgctaacccttttcgatttcgaaacttacgatg")
    seq2 = primer_design(Dseqrecord(
    "ccagttgctaaccttcttcgttcgaaacttacgatg"))

    seqs = [seq1, seq2]


    gg_designer = GoldenGateDesigner(seqs, lista_enz, circular=False)

    print()
    print('GoldenGateDesigner:')
    print(gg_designer)
    print()

    nova_lista = []

    for i in gg_designer[0]:
        # print(i)
        if i.cut(gg_designer[1]): # se tiver primers BsaI
            i1, i2, i3 = i.cut(gg_designer[1])
            nova_lista.append(i2)
            # print(i2.figure())
        else:
            nova_lista.append(i)
            # print(i.figure())


    # print(gg_designer[0])
    # e1,e2,e3  = gg_designer[0].cut(BsaI)
    # print()
    # print('Sequência 1 pós corte enzimatico: ')
    # print(e2.figure())


    # nova_lista = [e2, a2,b2,c2,d2]

    print()
    print('GoldenGateAssembler:')
    gg_assembler = GoldenGateAssembler(nova_lista)
    print(gg_assembler)

    print()

    for i in gg_assembler.values():
        print(i.figure())





    
