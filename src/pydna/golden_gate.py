from pydna.amplicon import Amplicon
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from goldenhinges import OverhangsSelector
# from Bio import Restriction
from Bio.Restriction import *
from pydna.assembly import Assembly
import Bio.Seq as Seq
from gg_graphs import *


def GoldenGateDesigner(seqs, enzs):

    # STEP 0
    # list all sticky ends...
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

            sticky_ends.append(end_5[1]) # será preciso saber a ponta 5 ou 3?
        if end_3[0] == 'blunt':
            sticky_ends_dict[seqs.index(sequence)].append((end_3[0])) # blunt
        else:
            sticky_ends_dict[seqs.index(sequence)].append((end_3[1])) # overhang
            sticky_ends.append(end_3[1])

    sticky_ends = list(set(sticky_ends))
    print(sticky_ends_dict)
    # print(sticky_ends)


    # STEP 1
    # compatible enzymes (enzymes that do not cut within the sequences)
    compatible_enzymes = enzs.copy()
    
    for e in enzs:
        for seq in seqs:
            if e.search(seq.seq) != []:
                compatible_enzymes.remove(e)
                break

    # print(compatible_enzymes)

    # STEP 2
    # loop over pairs of seqs

    # STEP 3 
    # design relevant primer if needed (add enzyme and comp sticky end)
    # FUNÇÃO AUXILIAR para loop do step 2
    
    
    def design_relevant_primer(seq, ovhg = None, end=5):
        primers = primer_design(seq)
        # print('teste primer')
        if ovhg is None:
            ovhg = 'aatt'# temporario # goldenhinges para forbidden sticky ends      
        if end == 5:
            seq = Amplicon(seq, forward_primer=(compatible_enzymes[0].site+ovhg+primers.forward_primer))
            # print(seq)
        elif end == 3:
            seq = Amplicon(seq, reverse_primer=(compatible_enzymes[0].site+ovhg+primers.reverse_primer))
            # print(seq)
        else:
            print('ERRO')

        # print('primer design')
        print(seq)
        
        return seq

    new_seqs = []

    for ind in range(len(seqs)): # index da sequencia na lista (até ao penultimo)
        if ind == len(seqs)-1:
            print('teste')
            new_seqs.append(seqs[ind])

            break
        print(ind)

        # if ind == 0 and sticky_ends_dict[0][1] == 'blunt':
        #     amp = design_relevant_primer(seqs[0]) ## goldenhinges ??
        #     new_seqs.append(amp)
        #     amp2 = design_relevant_primer(seqs[ind+1], end=3) # mas depois tem de atualizar na lista na proxima iteração -- como?
        #     new_seqs.append(amp2)
        #     continue
        # elif ind == 0 and sticky_ends_dict[0][1] != 'blunt':
        #     new_seqs.append(seq[0])
        #     continue
        # print(ind, ind+1)

        print(sticky_ends_dict[ind][1], sticky_ends_dict[ind+1][0])

        if sticky_ends_dict[ind][1] == 'blunt' and sticky_ends_dict[ind+1][0] == 'blunt': ## 2 blunts (2 primers)
            print('blunt test')
            amp = design_relevant_primer(seqs[ind]) ## goldenhinges ??
            new_seqs.append(amp) ## primeira sequencia
            amp2 = design_relevant_primer(seqs[ind+1], end=3) # mas depois tem de atualizar na lista na proxima iteração -- como?
            new_seqs.append(amp2)

        elif sticky_ends_dict[ind][1] != 'blunt' and sticky_ends_dict[ind+1][0] != 'blunt': ## dseq + dseq
            print('dseq + dseq test')
            new_seqs.append(seqs[ind])

        elif sticky_ends_dict[ind][1] == 'blunt' and sticky_ends_dict[ind+1][0] != 'blunt': ## blunt + dseq
            print('blunt + dseq test')
            # print(seqs[ind])
            ovhg = sticky_ends_dict[ind+1][0]
            print('ovhg: ', ovhg)
            amp = design_relevant_primer(seqs[ind], ovhg = ovhg, end=5)
            new_seqs.append(amp)

        elif sticky_ends_dict[ind][1] != 'blunt' and sticky_ends_dict[ind+1][0] == 'blunt': ## dseq + blunt
            print('dseq + blunt test')
            # print(seqs[ind])
            ovhg = sticky_ends_dict[ind][1]
            print('ovhg: ', ovhg)
            amp = design_relevant_primer(seqs[ind+1], ovhg = ovhg, end=3)
            new_seqs.append(amp)

        else: 
            print('ERRO')


        # if isinstance(seqs[ind], Dseqrecord) and isinstance(seqs[ind+1], Dseqrecord):
        #     # 2 dseqrecords - não acontece nada as sequencias
        #     # print('dseqrecords')
        #     new_seqs.append(seqs[ind])
        # elif isinstance(seqs[ind], Amplicon) and isinstance(seqs[ind+1], Amplicon):
        #     # 2 primers
        #     # print('amplicons')
        #     amp = design_relevant_primer(seqs[ind])
        #     new_seqs.append(amp)
        # elif isinstance(seqs[ind], Dseqrecord) and isinstance(seqs[ind+1], Amplicon):
        #     new_seqs.append(seqs[ind])
        #     # fazer primer na segunda seq?
        #     design_relevant_primer(seqs[ind], seqs[ind+1])
        # else: # amp + dseq
        #     pass
            # print('teste amp + dseqrecord')

    print(new_seqs)
    
    
    
    # STEP 4
    # if two amplicons, design both relevant primers (lista sticky ends serve aqui)
    def design_relevant_primers(seq, forbidden_ends = sticky_ends): # forbidden ends -> goldenhinges package
        
        ampl = primer_design(seq)
        
        return seq
    
    
    # STEP 5
    # Return a list with Dseqrecords and amplicons

    return new_seqs



##############################################################

    # STEP 1: Divide input into Dseqrecords and amplicons
    dseqrecords = [seq for seq in seqs if isinstance(seq, Dseqrecord)]
    amplicons = [seq for seq in seqs if isinstance(seq, Amplicon)]
    
    # STEP 2: Verify that all sequences are linear (apenas dseqrecords certo?)
    for seq in seqs:
        if seq.circular:
            raise ValueError("All sequences must be linear")
    
    # STEP 3: Check type and length of overhangs on dseqrecords
    # overhangs = {}
    # print(dseqrecords)

    # for seq in dseqrecords:
    #     dseq = seq.seq
    #     length = abs(dseq.ovhg)
    #     w = dseq.watson[:length]
    #     # c = dseq.crick[:length]
    #     overhangs[w] = length
    # print('overhangs\n', overhangs)
    
    
    # STEP 4: Select compatible enzymes

    # check if type IIs restriction enzyme
    # criar um dicionario com as enzimas de type IIS ? e permitir na lista apenas essas
    # compatible_enzymes = [enzyme for enzyme in enz if (enzyme.site == BsaI.site or enzyme.site == BsmBI.site)]

    # comp_enzymes = [] # se tiver o site da enzima em alguma sequencia, é adicionado a uma lista de enzimas final
    # for enzyme in compatible_enzymes:
    #     for seq in dseqrecords: 
    #         if enzyme.search(seq) is not None:
    #             comp_enzymes.append(enzyme)
    #     for seq in amplicons: 
    #         if enzyme.search(seq) is not None:
    #             comp_enzymes.append(enzyme)                


    # # STEP 5: Check all amplicons for internal restriction sites, select Golden Gate enzymes that do not cut internally
    # for amplicon in amplicons:
    #     for enzyme in compatible_enzymes:
    #         if enzyme.search(amplicon.seq) is not None:
    #             # gg_enzyme = enzyme
    #             print(enzyme.search(amplicon.seq))
    #             comp_enzymes.append(enzyme)
    


    # STEP 6: Use goldenhinges to design assembly fragments  ### GOLDENHINGES package doesn't seem to work 

    

    # STEP 7: Add appropriate tails to PCR primers of amplicons
    # overhangs tem de corresponder à ordem que queremos ligar os fragmentos (imagem artigo P1 liga-se ao amplicon, o P2 liga-se ao P3 e o P4 liga-se ao outro lado do amplicon)


    lista = [] # lista dos fragmentos depois da enzima cortar (tem de ocorrer PCR aqui dentro?)


    # STEP 8: Return a list with Dseqrecords and amplicons
    return lista


def GoldenGateAssembler(seqs: list):
    '''
    GoldenGateAssembler is a function that receives a list of sequences (fragments) and tries all the possible assembles between them, returns a list/dictionary of all the possible sequences ligated (with only DNA ligase)

    Parameters
    ----------
    seqs : list
        fragments of DNA sequences

    Returns
    -------
    dic_paths : dict
        dictionary with the possible sequences ligated (keys are the order of the sequences; values are the dseqrecords of the sequences assembled)
    '''
    # cria um grafo em que cada nodo é uma sequencia e cada edge representa a ligação entre sequencias (cada combinação possível numa lista de sequencias)
    graph = graph_assembly(seqs)
    
    # encontra todos os paths - encontra todas as combinações de sequencias possiveis 
    paths = find_all_paths(graph)

    # usando os paths como keys, retorna as sequencias formadas em cada uma destes (sem repetições)  
    dic_paths = find_paths_seqs(paths, graph)
    
    return dic_paths

    

if __name__ == '__main__':
    print()
    # testar o codigo
    from Bio.Seq import Seq
    from pydna.dseqrecord import Dseqrecord
    from pydna.amplicon import Amplicon
    from pydna.primer import Primer
    from Bio.Restriction import BsaI, BsmBI, BamHI
    from pydna.dseq import Dseq
    from pydna.seq import Seq

    # # Create a Dseqrecord object
    # seq1 = Dseqrecord("atgccctaa")
    # seq1.add_feature()
    # print(seq1.figure())

    # #instruçoes pydna geral PCR + cut and paste cloning

    # from pydna.design import primer_design
    # amplicon = primer_design(seq1, limit=3, target_tm = 0)
    # print(amplicon.figure())
    # print(amplicon.forward_primer)

    # forward_primer = 'ccccGGATCC' + amplicon.forward_primer
    # reverse_primer = 'ttttGGATCC' + amplicon.reverse_primer
    
    # from pydna.amplify import pcr
    # pcr_product = pcr(forward_primer, reverse_primer, seq1, limit=3)
    # print(pcr_product.figure())


    # a,b,c = pcr_product.cut(BamHI)
    # print(a.figure())
    # print(b.figure())
    # print(c.figure())

    # print(seq1.seq)

    a = Dseqrecord("CTTAAGatgccctaaccccccctaacccacGAATTC")
    b = Dseqrecord("CTGGAGatgccctaaccatgtagtaaaaaaatgccGAATTC")
    b3 = Dseqrecord("CTGGAGatgccctaaccatgtagtatgccCTGGAG")
    b2 = Dseqrecord("CTGGAGatgccctaaccccaggaattagagCTGGAG")
    c = Dseqrecord("GAATTCatgcccgggggggggggcagtacagtaCTTAAG")
    d = Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    # ampl = primer_design(d)
    # ampl.figure()

    lista_enz = [EcoRI, BsaI, BsmBI, BamHI, BpmI]
    
    # print(b.cut(EcoRI))

    s1,a2  = a.cut(EcoRI)
    s2, b1 = b.cut(EcoRI)
    # b1,s2 = b2.cut(BpmI)
    # print(s2.figure())
    c1,s3 = c.cut(EcoRI)
    
    lista = [a,b,c]
    # lista_seqs = [s1,b2,s3]

    lista_seqs = [s1,s2,s3]
    for i in lista_seqs:
        print(i.figure())
        
    # ampl = primer_design(a)
    # print(ampl.figure())

    print(GoldenGateAssembler(lista_seqs))
    print()

    # print(GoldenGateDesigner(lista, lista_enz))
    print(GoldenGateDesigner(lista_seqs, lista_enz))
    lista_gg = GoldenGateDesigner(lista_seqs, lista_enz)
    print()
    for i in lista_gg:
        print(lista_gg.index(i))
        print(i)

    print()
    print(GoldenGateAssembler(lista_gg))
    print()


    




