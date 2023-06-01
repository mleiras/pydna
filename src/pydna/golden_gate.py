from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from Bio.Restriction import *
from golden_gate_auxiliar import *

from pydna.amplicon import Amplicon
from pydna.seqrecord import SeqRecord
from pydna.amplify import pcr
from goldenhinges import OverhangsSelector
from pydna.assembly import Assembly
import Bio.Seq as Seq



def GoldenGateDesigner(seqs: list, enzs: list, circular: bool = True) -> list:

    # STEP 1:  list all sticky ends
    sticky_ends_dict = list_sticky_ends(seqs)

    # STEP 2:  compatible enzymes (enzymes that do not cut within the sequences)
    compatible_enzymes = compatible_enzyme(seqs, enzs) 

    # STEP 3: for loop over sequences
    # design relevant primer if needed (add enzyme and comp sticky end)
    new_seqs = design_seqs(seqs, sticky_ends_dict, compatible_enzymes, circular)
   
    # STEP 4: Return a list with Dseqrecords and amplicons

    return new_seqs



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

    lista_seqs = [a2,b2,c2,d2]

    # print('Sequencias: \n')
    # for i in lista_seqs:
    #     print(i.figure())

    e = primer_design(Dseqrecord("ccagttgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")) # AMPLICON
    

    lista_seqs2 = [a2,b2,c2,d2,e]

    lista_seqs3 = [e, a2,b2,c2,d2]

    gg_designer = GoldenGateDesigner(lista_seqs3, lista_enz)
    print()
    print('GoldenGateDesigner:')
    print(gg_designer)
    print()
    for i in gg_designer:
        print()
        print(i.figure())

    # print(gg_designer[0])
    e1,e2,e3  = gg_designer[0].cut(BsaI)
    print()
    print(e2.figure())

    nova_lista = [e2, a2,b2,c2,d2]
    print('GoldenGateAssembler:')
    print(GoldenGateAssembler(nova_lista))





# def GoldenGateDesigner(seqs: list, enzs: list, circular: bool = True) -> list:

#     sticky_ends_dict = list_sticky_ends(seqs)

#     compatible_enzymes = compatible_enzyme(seqs, enzs) 

#     new_seqs = design_seqs(seqs, sticky_ends_dict, compatible_enzymes, circular)
   
#     return new_seqs


# def list_sticky_ends(seqs):
    # sticky_ends_dict = {} 
    # for sequence in seqs:
    #     end_5 = sequence.seq.five_prime_end()
    #     end_3 = sequence.seq.three_prime_end()
    #     sticky_ends_dict[seqs.index(sequence)] = []
    #     if end_5[0] == 'blunt':
    #         sticky_ends_dict[seqs.index(sequence)].append((end_5[0])) # blunt
    #     else:
    #         sticky_ends_dict[seqs.index(sequence)].append((end_5[1])) # overhang
    #     if end_3[0] == 'blunt':
    #         sticky_ends_dict[seqs.index(sequence)].append((end_3[0])) # blunt
    #     else:
    #         sticky_ends_dict[seqs.index(sequence)].append((end_3[1])) # overhang

    # return sticky_ends_dict


# def compatible_enzyme(seqs, enzs):

#     type_IIS_enzymes = [AcuI, AlwI, BaeI, BbsI, BbvI, BccI, BceAI, BcgI, BciVI, BcoDI, BfuAI, BmrI, BpmI, BpuEI, BsaI, BsaXI, BseRI, BsgI, BsmAI, BsmBI, BsmFI, BsmI, BspCNI, BspMI, BspQI, BsrDI, BsrI, BtgZI, BtsCI, BtsI, BtsIMutI, CspCI, EarI, EciI, Esp3I, FauI, FokI, HgaI, HphI, HpyAV, MboII, MlyI, MmeI, MnlI, NmeAIII, PaqCI, PleI, SapI, SfaNI]
    
#     comp_enzymes = enzs.copy()
    
#     for e in enzs:
#         if e not in type_IIS_enzymes: # only type IIs restriction enzymes
#             comp_enzymes.remove(e)
#         else:
#             for seq in seqs:
#                 if e.search(seq.seq) != []:
#                     comp_enzymes.remove(e)
#                     break

#     if comp_enzymes == []:
#         raise ValueError('''No type IIs restriction enzymes available on this list
#                         that do not cut within the sequences''')
    
#     return comp_enzymes


# def design_seqs(seqs, sticky_ends_dict, compatible_enzymes, circular):

#     def design_primers(sequence, compatible_enzymes, f_ovhg, r_ovhg, ind):
#         if f_ovhg == 'blunt':
#             pass
#             # goldenhinge
#             # overhang_f = 
#         else:
#             overhang_f = SeqRecord(f_ovhg).reverse_complement().seq
#         if r_ovhg == 'blunt':
#             pass
#             # goldenhinge
#             # overhang_r = 
#         else:
#             overhang_r = SeqRecord(r_ovhg).reverse_complement().seq

#         fp = compatible_enzymes[0].site + 'a' + overhang_f+ sequence.forward_primer
#         rp = compatible_enzymes[0].site + 'a' + overhang_r + sequence.reverse_primer

#         pcr_prod = pcr(fp, rp, e)
#         new_seq = Dseqrecord(pcr_prod)

#         sticky_ends_dict[ind][0] = overhang_f
#         sticky_ends_dict[ind][1] = overhang_r

#         return new_seq

#     new_seqs = []

#     for ind in range(len(seqs)-1): 

#         if ind == 0:
#             if circular:
#                 previous_sticky_end = sticky_ends_dict[len(seqs)-1][1]
#             else:
#                 previous_sticky_end = 'blunt' # temporario -> mudar para None
#         else:
#             previous_sticky_end = sticky_ends_dict[ind-1][1]

#         if ind == len(seqs)-1:
#             if circular:
#                 next_sticky_end = sticky_ends_dict[0][0]
#             else:
#                 next_sticky_end = 'blunt' # temporario -> mudar para None
#         else:
#             next_sticky_end = sticky_ends_dict[ind+1][0] 

#         if sticky_ends_dict[ind][0] == 'blunt' and sticky_ends_dict[ind][1] == 'blunt':
#             new_seq = design_primers(seqs[ind], compatible_enzymes = compatible_enzymes,
#                                      f_ovhg = previous_sticky_end, r_ovhg = next_sticky_end, ind=ind)
#             new_seqs.append(new_seq)

#         elif sticky_ends_dict[ind][0] != 'blunt' and sticky_ends_dict[ind][1] != 'blunt':
#             new_seqs.append(seqs[ind]) 

#         else:
#             raise ValueError('ERRO DESIGN PRIMER')

#     return new_seqs
    

