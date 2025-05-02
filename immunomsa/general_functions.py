import sys
def split_list(input_list, num_sublists):
    """
    Split lists into a specified number of (num_sublists) list of lists
    input_list: list
    num_sublists: int
    """
    # Calculate the size of each sublist
    sublist_size = len(input_list) // num_sublists
    remainder = len(input_list) % num_sublists

    # Initialize variables
    start = 0
    end = 0
    sublists = []

    # Create sublists
    for i in range(num_sublists):
        start = end
        end += sublist_size + (1 if i < remainder else 0)
        sublists.append(input_list[start:end])

    return sublists



def read_accs_and_sequences_from_fasta(infile):
        """
        Input: readfile: Fasta file. 
        Outputs: List of tuples. Containing accs and sequences, e.g. [(acc, aTHNtem..)..()]. 
        """
        
        if not infile.is_file():
            sys.exit(f"The input file was invalid: {infile}")

        accs = list()
        sequences = list()
        seq = ""

        read_acc = False    
        infile = open(infile, "r")
        readfile = infile.readlines()
        
        infile.close()

        for line in readfile:
            line = line.strip()
            if line.startswith(">"):
                acc = line.split(">")[1]
                if read_acc:
                    accs.append(acc)
                    sequences.append(seq)
                    #reset sequence string
                    seq = ""
                #catch first accesion.
                else:
                    accs.append(acc)
            else:
                seq += line
                read_acc = True

        #get last sequence
        sequences.append(seq)

        if accs == False or sequences == False:
            sys.exit(f"No accessions or sequences found in fasta file. Please check file: {infile}")

     

        return accs, sequences