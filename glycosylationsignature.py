#Question 5.
#1. Translate the DNA sequence to protein sequences.
#2. check for Glycosylation signatures.
#3. Convert the signature region to an upper case. 

def glyco_signat(x):
        
    def translate(a):
        aa = ""
        codontable = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
            }    
        
        length = 0
        if len(a)%3 == 0:
            length = len(a)
        else:
            length = len(a) - (len(a)%3)
            
        for n in range(0,length,3):
            
            aa = aa + codontable[a[n:n+3]]
        global protein 
        protein = aa.lower()
        
        
    
    import sys
    import os.path
    import re
    
        
    os.path.dirname(x)
        
    m = open(x,'r')            #open file
    firstline = m.readline()    #read the first line of FASTA file
    format_checker = "CLUSTAL multiple sequence alignment by MUSCLE (3.8)"   #check whether the file format is consistent with Clustal
    lines = m.readlines()
    seq = {}
    seq_orig = {}
    matches = ''  #we will add the star sequences to this string.
    b = '' #blank string
    #seqid =        #sequence i.d.
    
    if firstline[:8] == format_checker[:8]:    #condition to perform question 2
        proceed = "The format allows us to proceed."
        regx_star = re.compile("^\*")      # if line starts with '*'
        regx_blank = re.compile("")        # if the line is a blank line
        
        counter = 0
        for line in lines:
            sequence = line.split()   #what is proposed to be the sequence based 
            
            if sequence != []:         # is it a blank line?
                o=regx_star.match(sequence[0])   #o = the value of the match of * reg expression
                n=regx_blank.match(sequence[0])   #n = the value of the match of [] reg expression
                
                if counter == 0:        #for the first set of blank lines and sequences, we create our dictionary
                    if o == None:
                        seq[sequence[0]] = sequence[1]
                    else:
                        matches = "".join(sequence)   # now we add the line that contains the stars to the match string.
                        counter = 1      # change counter to 1 to show that now we are no longer creating a dictionary but modifying values of the dictionary
                
                elif counter == 1:    # now we are modifying values of the dictionary
                    if o == None:      # if the result of 'o' is not a star sequences and is not a blank sequences, therefore it is a sequence that contains nucleotides
                        seq[sequence[0]] = seq[sequence[0]] + sequence[1]
                    else:               #at this stage we want to count the number of stars therefore we add the star sequences to match and remove the spaces.
                        matches = matches + "".join(sequence)
        
        seqno = 0
        seq_orig_len = 0
        average_seqlen = 0
        null_proteins = {}
        for k in seq:
            seq_orig[k] = seq[k].replace('-','')     #seq is the sequence after alignment i.e. without any dashes introduced
            seq_orig_len = seq_orig_len + len(seq_orig[k])   # seq_orig is the sequence before alignement i.e. with dashes introduced
            seq_len = len(seq[k])
            seqno += 1
        
        print "\n",
        for nucleo in seq_orig:
            translate(seq_orig[nucleo]) 
            regx_glyco = re.compile("n[^p][st]")
            result = regx_glyco.findall(protein)
            
            new_protein = ''
            if result != []:
                a = ''
                new_protein = protein
                for i in result:
                    a = i.upper()
                    new_protein = new_protein.replace(i,a)
                print "Glycosylation Signatures found in sequence %s: \n\n" %(nucleo), "'%s'"%(new_protein) 
                
            else:
                null_proteins[nucleo] = nucleo
        
        print " \n\nThere are no signatures in sequences:", 
        
        spacer = '' 
        for i in null_proteins:
            print spacer, null_proteins[i], 
            if i != []:
                spacer = 'and'
        print "\n"
            
           


    
