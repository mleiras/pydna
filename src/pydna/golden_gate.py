from pydna.amplicon import Amplicon
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from goldenhinges import OverhangsSelector
from Bio.Restriction import RestrictionBatch
from pydna.assembly import Assembly
import Bio.Seq as Seq

def GoldenGateDesigner(seqs, enz):
    # Step 1: Divide input into Dseqrecords and amplicons
    dseqrecords = [seq for seq in seqs if isinstance(seq, Dseqrecord)]
    amplicons = [seq for seq in seqs if isinstance(seq, Amplicon)]
    
    # Step 2: Verify that all sequences are linear
    for seq in dseqrecords:
        if seq.circular:
            raise ValueError("All sequences must be linear")
    
    # Step 3: Check type and length of overhangs on dseqrecords
    # type?


    
    # Step 4: Select compatible enzymes
    compatible_enzymes = [enzyme for enzyme in enz if (enzyme.site == BsaI.site or enzyme.site == BsmBI.site)]
    print(compatible_enzymes)
    
    # Step 5: Check all amplicons for internal restriction sites, select Golden Gate enzymes that do not cut internally
    for amplicon in amplicons:
        for enzyme in compatible_enzymes:
            if enzyme.search(amplicon.seq):
                print(enzyme.search(amplicon.seq))
                # compatible_enzymes.remove(enzyme)
    
    ### GOLDENHINGES package doesn't seem to work 
    # Step 6: Use goldenhinges to design assembly fragments
    # selector = OverhangsSelector()
    # overhangs = selector.select_overhangs(dseqrecords)

    # Step 7: Add appropriate tails to PCR primers of amplicons
    for amplicon in amplicons:
        # for overhang in overhangs:
        #     if overhang[0] in amplicon.seq:
        #         amplicon.forward_primer = Seq.Seq(overhang[1] + str(amplicon.forward_primer.seq))
        #     if overhang[2] in amplicon.seq:
        #         amplicon.reverse_primer = Seq.Seq(overhang[3] + str(amplicon.reverse_primer.seq))
        if amplicon.forward_primer is not None:
            amplicon.forward_primer = Seq(str(amplicon.forward_primer.seq))
        if amplicon.reverse_primer is not None:
            amplicon.reverse_primer = Seq(str(amplicon.reverse_primer.seq))
    
    # Step 8: Return a list with Dseqrecords and amplicons
    return dseqrecords + amplicons


def GoldenGateAssembler(seqs, enz):

    # Step 1: Cut all Amplicons
    cut_seqs = []
    for seq in seqs:
        for enzyme in enz:
            fragments = enzyme.catalyze(seq.seq if isinstance(seq, Dseqrecord) else seq)
            cut_seqs.extend([Dseqrecord(frag) for frag in fragments])
    print(cut_seqs)
    
    # Step 2: Try to add all sequences in order
    assembly = Assembly(cut_seqs)
    print(assembly)
    try:
        result = assembly.assemble_linear()
        # result = assembly.assemble_circular()
    except IndexError:
        # Step 3: If assembly not possible raise exception
        raise ValueError("Assembly not possible with provided sequences and enzymes")
    
    # Step 4: return Dseqrecord
    return result



def check_overhangs(dseqrecord):
    """
    This function checks the overhangs of a Dseqrecord object.
    It assumes that the overhangs are located at the ends of the sequence.
    """
    # Get the Dseq object from the Dseqrecord
    dseq = dseqrecord.seq
    
    # Check for overhangs
    print(dseq.watson, dseq.crick)
    if len(dseq.watson) > len(dseq.crick):
        print("5' overhang on top strand:", dseq.watson[len(dseq.crick):])
    elif len(dseq.watson) < len(dseq.crick):
        print("5' overhang on bottom strand:", dseq.crick[len(dseq.watson):])

    return True



    

if __name__ == '__main__':
    # testar o codigo
    from Bio.Seq import Seq
    from pydna.dseqrecord import Dseqrecord
    from pydna.amplicon import Amplicon
    from pydna.primer import Primer
    from Bio.Restriction import BsaI, BsmBI, BamHI, EcoRI, HindIII
    from pydna.dseq import Dseq
    from pydna.seq import Seq

    # Create a Dseqrecord object
    seq1 = Dseqrecord("atgccctaa")
    seq1.add_feature()
    print(seq1.figure())

    #instruçoes pydna geral PCR + cut and paste cloning

    from pydna.design import primer_design
    amplicon = primer_design(seq1, limit=3, target_tm = 0)
    print(amplicon.figure())
    print(amplicon.forward_primer)

    forward_primer = 'ccccGGATCC' + amplicon.forward_primer
    reverse_primer = 'ttttGGATCC' + amplicon.reverse_primer
    
    from pydna.amplify import pcr
    pcr_product = pcr(forward_primer, reverse_primer, seq1, limit=3)
    print(pcr_product.figure())


    a,b,c = pcr_product.cut(BamHI)
    print(a.figure())
    print(b.figure())
    print(c.figure())

    # print(seq1.seq)

    # Check the overhangs
    # five_prime_overhang, three_prime_overhang = check_overhangs(dseqrec)

    # print("5' overhang:", five_prime_overhang)
    # print("3' overhang:", three_prime_overhang)





# if __name__ == '__main__':
#     # testar o codigo
#     from Bio.Seq import Seq
#     from pydna.dseqrecord import Dseqrecord
#     from pydna.amplicon import Amplicon
#     from pydna.primer import Primer
#     from Bio.Restriction import BsaI, BsmBI, EcoRI, HindIII
#     from pydna.dseq import Dseq

    # # Create some test sequences
    # seq1 = Dseqrecord("ATGCGTACGTACGATCGTAGCTACGTACG")
    # seq2 = Dseqrecord("ATGCGTACGTACGATCGTAGCTACGTACG")
    # seq3 = Dseqrecord("ATGCGTACGTACGATCGTAGCTACGTACG")

    # # Create some test amplicons
    # primer1 = Primer("ATGCGTACGTACGATCGTAGCTACGTACG")
    # primer2 = Primer("CGTACGTACGATCGTAGCTACGTACGCAT")
    # amplicon1 = Amplicon(seq1, primer1, primer2)
    # amplicon2 = Amplicon(seq2, primer1, primer2)

    # # Create a list of sequences and enzymes
    # seqs = [seq1, seq2, seq3, amplicon1, amplicon2]
    # enz = [BsaI, BsmBI]

    # result = GoldenGateDesigner(seqs, enz)

    # # for res in result:
    # #     print(res)

    # print('\n GOLDENGATE ASSEMBLY')


    # # Create some test sequences
    

    # # Test the function
    # result2 = GoldenGateAssembler(result, [EcoRI, HindIII])
    # print(result2)


    ############################

    # Create the fragments sequences
    # fragment1 = Dseqrecord("ATGGTCTCAAGCTTACCGGTAAGCTTGGTCTCATG")
    # fragment2 = Dseqrecord("ATGGTCTCAAGCTTACCGGTAAGCTTGGTCTCATG")

    # # Create the destination vector
    # amplicon = Amplicon('ATGCGTCTCAAGCTTACCGGTAGCTTGGTCTCGATGCGTCTCAAGCTTACCGGTAGCTTGGTCTCGATG')

    # # Create a list of sequences and enzymes
    # seqs = [amplicon, fragment1, fragment2]
    # enz = [BsaI,BsmBI]

    # result = GoldenGateDesigner(seqs, enz)

    # print(result)

    # for res in result:
    #     print(res)

    # print('\n GOLDENGATE ASSEMBLY')


    # # Create some test sequences
    

    # # Test the function
    # result2 = GoldenGateAssembler(result, [EcoRI, HindIII])
    # print(result2,'\n')

    #############

    


