from pydna.amplicon import Amplicon
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from goldenhinges import OverhangsSelector
from Bio.Restriction import RestrictionBatch
from pydna.assembly import Assembly
import Bio.Seq as Seq
from gg_graphs import *

def GoldenGateDesigner(seqs, enz):
    # STEP 1: Divide input into Dseqrecords and amplicons
    dseqrecords = [seq for seq in seqs if isinstance(seq, Dseqrecord)]
    amplicons = [seq for seq in seqs if isinstance(seq, Amplicon)]
    
    # STEP 2: Verify that all sequences are linear (apenas dseqrecords certo?)
    for seq in dseqrecords:
        if seq.circular:
            raise ValueError("All sequences must be linear")
    
    # STEP 3: Check type and length of overhangs on dseqrecords
    overhangs = {}

    for seq in dseqrecords:
        dseq = seq.seq
        length = abs(dseq.ovhg)
        w = dseq.watson[:length]
        # c = dseq.crick[:length]
        overhangs[w] = length
    print('overhangs\n', overhangs)
    
    
    # STEP 4: Select compatible enzymes

    # check if type IIs restriction enzyme
    # criar um dicionario com as enzimas de type IIS ? e permitir na lista apenas essas
    compatible_enzymes = [enzyme for enzyme in enz if (enzyme.site == BsaI.site or enzyme.site == BsmBI.site)]

    comp_enzymes = [] # se tiver o site da enzima em alguma sequencia, é adicionado a uma lista de enzimas final
    for enzyme in compatible_enzymes:
        for seq in dseqrecords: 
            if enzyme.search(seq) is not None:
                comp_enzymes.append(enzyme)
        for seq in amplicons: 
            if enzyme.search(seq) is not None:
                comp_enzymes.append(enzyme)                


    # STEP 5: Check all amplicons for internal restriction sites, select Golden Gate enzymes that do not cut internally
    for amplicon in amplicons:
        for enzyme in compatible_enzymes:
            if enzyme.search(amplicon.seq) is not None:
                # gg_enzyme = enzyme
                print(enzyme.search(amplicon.seq))
                comp_enzymes.append(enzyme)
    


    # STEP 6: Use goldenhinges to design assembly fragments  ### GOLDENHINGES package doesn't seem to work 

    

    # STEP 7: Add appropriate tails to PCR primers of amplicons
    # overhangs tem de corresponder à ordem que queremos ligar os fragmentos (imagem artigo P1 liga-se ao amplicon, o P2 liga-se ao P3 e o P4 liga-se ao outro lado do amplicon)


    lista = [] # lista dos fragmentos depois da enzima cortar (tem de ocorrer PCR aqui dentro?)


    # STEP 8: Return a list with Dseqrecords and amplicons
    return lista


def GoldenGateAssembler(seqs, enz):

    # queremos todas as possibilidades de ligar os fragmentos (porque podem haver problemas no design dos fragmentos e este package serve para isso)

    # Step 1: Cut all Amplicons
    # amplicons, dseqrecords = seqs
    amplicons = seqs
    # passar por todas as sequencias com todas as enzimas até cortar?
    for amp in amplicons:
        for enzyme in enz:
            amp1 = amp.cut(enzyme) # tem de ser com o try (não vai cortar com qualquer enzima)
    
            seqs.extend(amp1)
    
    
    # Step 2: Try to add all sequences in order
    graph = graph_assembly(seqs)
    paths = find_all_paths(graph)
    dic_paths = find_paths_seqs(paths, graph)
    for i,x in dic_paths.items():
        print(i)
        print(x.figure())


    # Step 4: return Dseqrecord # ou vários?
    # return result


    

if __name__ == '__main__':
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

    a = Dseqrecord("CTTAAGatgccctaaccccGAATTC")
    b = Dseqrecord("GAATTCatgccctaaccccGAATTC")
    c = Dseqrecord("GAATTCatgcccgggggggggggccCTTAAG")

    s1,a2  = a.cut(EcoRI)
    b1,s2,b3 = b.cut(EcoRI)
    c1,s3 = c.cut(EcoRI)
    
    lista_seqs = [a,b,c]

    GoldenGateAssembler(lista_seqs, [EcoRI])

