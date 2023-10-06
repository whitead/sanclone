#this script takes a input sequence and outputs a sequence ready for restriction cloning 
#input - gene sequence, restriction sites to avoid, flanking restriction sites, maximum size of DNA fragment to synthesize
#output - language on what was done, final gene sequence, gene fragments to order

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

#this function cleans up ambiguous characters in DNA
def dna_clean(sequence = ""):
    
    #make all uppercase
    sequence = sequence.upper()
    
    #convert RNA to DNA
    sequence = sequence.replace("U","T")
    
    return sequence


#this function returns the reverse complement of a DNA sequence
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(dna)])

#this outputs a random codon for each amino acid for the purposes of removing restriction sites
def codonchoice(aminoacid):
    aminoacid = aminoacid.upper()
    codon_dict = {'A':["GCA","GCC","GCG","GCT"],
    'C':["TGC","TGT"],
    'D':["GAC","GAT"],
    'E':["GAA","GAG"],
    'F':["TTC","TTT"],
    'G':["GGA","GGC","GGG","GGT"],
    'H':["CAC","CAT"],
    'I':["ATA","ATC","ATT"],
    'K':["AAA","AAG"],
    'L':["CTA","CTC","CTG","CTT","TTA","TTG"],
    'M':["ATG"],
    'N':["AAC","AAT"],
    'P':["CCA","CCC","CCG","CCT"],
    'Q':["CAA","CAG"],
    'R':["AGA","AGG","CGA","CGC","CGG","CGT"],
    'S':["AGC","AGT","TCA","TCC","TCG","TCT"],
    'T':["ACA","ACC","ACG","ACT"],
    'V':["GTA","GTC","GTG","GTT"],
    'W':["TGG"],
    'Y':["TAC","TAT"]}
    
    return random.choice(codon_dict[aminoacid])


#this function removes restriction sites (input is name of restriction sites) from the gene sequence by replacing codons with 
def remove_restriction(restriction_sites = [], sequence = "", ):
    
    #default restriction enzymes to remove - these include common golden gate enzymes
    default_restriction = ['BSAI', 'BBSI', 'BSMBI']
    
    #convert restriction sites to uppercase
    restriction_sites = [x.upper() for x in restriction_sites]
    
    all_restriction_sites = default_restriction + restriction_sites
    
    #file of restriction enzymes
    restriction_csv = 'restriction_enzymes.csv'
    restriction_df = pd.read_csv(restriction_csv,sep=',')
    
    new_nucleotide_sequence = sequence.upper()
    
    #convert restriction enzyme names to sequences
    restriction_sites_seqs = []
    for restriction_site in all_restrition_sites:
        site_row = restriction_df.loc[restriction_df['Name'] == restriction_site]
        if len(site_row.shape[0]) > 0:
            this_seq = site_row.iloc[0]['Sequence_cleaned']
            restriction_sites_seqs.append(this_seq)
            restriction_sites.seqs.append(reverse_complement(this_seq))
    
    
    #iterate through the restriction site sequences and remove them
    for restriction_site in restriction_sites_seqs:
        restriction_site = restriction_site.upper()
        
        while True:
            #find all coordinates of restriction site
            coordinates = [m.start() for m in re.finditer(restriction_site, new_nucleotide_sequence)]
            for coordinate in coordinates:

                #find the index of the start of the codon
                remainder = coordinate % 3
                codon_start = coordinate - remainder

                while True:      
                    #determine the codon
                    codon = new_nucleotide_sequence[codon_start:codon_start + 3]

                    #translate the codon
                    amino_acid = str(Seq(codon).translate())

                    #if it's an ATG or stop codon, then move on to the next codon instead
                    if amino_acid == 'M' or amino_acid == '*' or amino_acid =='W':
                        skip = random.choice([3])
                        codon_start = codon_start + skip                        
                    else:
                        break
                        
                #comeup with new codon
                while True:
                    new_codon = codonchoice(amino_acid)
                    if new_codon != codon:
                        break

                #adjust the sequence 
                new_nucleotide_sequence = new_nucleotide_sequence[0:codon_start] + new_codon + new_nucleotide_sequence[codon_start+3:]

            #check if all sites are successfully removed
            coordinates = [m.start() for m in re.finditer(restriction_site, new_nucleotide_sequence)]
            if len(coordinates) == 0:
                break
            else:
                break
      
    return new_nucleotide_sequence  

#this function breaks the insert into chunks for DNA synthesis
def eblockerizer(sequence, max_fragment_length = 500, min_fragment_length = 300, min_gibson_overlap = 25, min_overlap_Tm = 55):
    
    full_sequence = sequence.upper()
    
    current_part = 0
    sequences = []
    remaining_sequence = full_sequence
    
    while True:
        current_part = current_part + 1
        
        #is sequence already < max_fragment_length?
        if len(remaining_sequence) <= max_fragment_length:
            sequences.append(["{}_{}".format(seq_name,current_part),remaining_sequence])
            break
        else:
            #can the rest be split into two equal pieces?
            if len(remaining_sequence) <= 2*max_fragment_length-70:
                target_length = len(remaining_sequence)//2
            else:
                target_length = max_fragment_length
            
            #check if Tm of overlap is greater than 55C and extend overlap if not
            i=0
            while True:
                this_overlap = remaining_sequence[target_length-min_gibson_overlap-i:target_length]
                this_overlap_Tm = mt.Tm_NN(Seq(this_overlap))
                
                if this_overlap_Tm > min_overlap_Tm:
                    #output this fragment and update the remaining sequence
                    sequences.append(["{}_{}".format(seq_name,current_part),remaining_sequence[0:target_length]])
                    remaining_sequence = remaining_sequence[target_length-min_gibson_overlap-i:]
                    break
                else:
                    #if needed extent the overlap
                    i = i+1
    return sequences


#this function actually does the full task
def do_the_thing(sequence = '', restriction_sites = [], cloning_sites = [], max_fragment_length = 500):
    sequence = dna_clean(sequence)
    new_sequence = remove_restriction(restriction_sites = restriction_sites[],sequence = sequence)
    
    fragments = eblockerizer(sequence = new_sequence)
    #add the new sequence to the State
    
    return "The following restriction sites were removed:" + restriction_sites
    
    