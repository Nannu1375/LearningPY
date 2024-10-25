#tally ACGT count
def nuclCount(input_seq): 
    #turn all bases to upper case: cuts of case issues
    cap_seq = input_seq.upper() 

    #create an empty directory to store the bases and their counts
    count_dir = {} 
    
    #count nucleotides
    for base in cap_seq:
        if base in count_dir:

            #adds 1 to existing count
            count_dir[base] += 1 
        else:

            #allots 1 count to the base
            count_dir[base] = 1  

    #retrieves the values of each key (bases here) in a list as integers
    base_tally = [count_dir.get(nuc, 0) for nuc in ['A', 'C', 'G', 'T']]

    #converts the list of integers into a string of the key values, separated by space
    result = ' '.join(map(str, base_tally))

    return result

def transcription(input_seq):
    #convert to upper case
    cap_seq = input_seq.upper()
    
    #DNA to RNA conversion
    rna_seq = cap_seq.replace('T', 'U')

    return rna_seq

def reverse_comp(input_seq):
    #convert to upper case
    cap_seq = input_seq.upper()

    #create a complement table
    comp_dir = cap_seq.maketrans('ATCG', 'TAGC')

    #create the complement of the strand
    comp = cap_seq.translate(comp_dir)

    #reverse the complement
    rev_comp = comp[::-1]

    return rev_comp

def modFibonacci(months, pairs):
    #throws 1 if months are not enough
    if months < 3:
        return 1
        
    #1st and 2nd month pairs; genOld = F(n-2), genNew = F(n-1)
    genOld, genNew = 1,1 
    
    #create a loop to iterate through the months
    for iteration in range(3, months+1):

        #3 pairs from each older generation, offspring = Fn
        offspring = genNew + pairs * genOld
        
        #replace the new set of pairs
        genOld, genNew = genNew, offspring

    #return final number of offsprings
    return offspring

def gcContent(sequence):
    """Calculate GC content percentage in a DNA sequence."""
    gc_count = (sequence.upper().count("G") + sequence.upper().count("C")) / len(sequence) * 100
    return gc_count

def parse_fasta(file_path):
    """Parse FASTA file and return a dictionary of sequences with their IDs."""
    sequences = {}
    with open(file_path, 'r') as f:
        seq_id = None
        seq_data = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:  # If seq_id already exists, store the previous sequence
                    sequences[seq_id] = ''.join(seq_data)
                seq_id = line[1:]  # New sequence ID, remove ">"
                seq_data = []  # Reset sequence data
            else:
                seq_data.append(line)  # Collect sequence lines
        if seq_id:
            sequences[seq_id] = ''.join(seq_data)  # Store the last sequence
    return sequences

def find_highest_gc_content(file_path):
    """Find the sequence with the highest GC content from a FASTA file."""
    sequences = parse_fasta(file_path)
    max_gc_id = None
    max_gc_content = 0
    for seq_id, sequence in sequences.items():
        gc_content = gcContent(sequence)
        if gc_content > max_gc_content:
            max_gc_content = gc_content
            max_gc_id = seq_id
    return max_gc_id, max_gc_content

