import sys
import os.path
import re
def export_to_fasta(x):
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
        return aa
            
    def get_fasta(a,b):
        d = b[0]
        for i in range(1,len(b)):
            width = 60
            if (i+1)%width == 0:
                d = d + b[i] + '\n'
            else:
                d = d + b[i]
        
        c = ">" + a + "\n" + d + '\n'       
        return c
    
    # #######################################    

    
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
        print "%-16s %s"%("\nFile name: ",x), '\n'
        
        option = raw_input('If you want to create Multiple files(one per sequence) press M \n or press S for a single file that includes all the sequences \n:')
        DNAorProtein = raw_input('Do you want to save the DNA sequences(input D) or the proteins (input P): ')
        
       # ######################################################################################################     
        protein = ''
        if DNAorProtein == 'P' or DNAorProtein == 'p':
            # #################
            for k in seq:
                seq_orig[k] = seq[k].replace('-','')     #seq is the sequence after alignment i.e. without any dashes introduced
                seq_orig_len = seq_orig_len + len(seq_orig[k])   # seq_orig is the sequence before alignement i.e. with dashes introduced
                seq_len = len(seq[k])
                seqno += 1
                
            if option == 'M' or option == 'm':          
                for nucleo in seq_orig:
                    protein = translate(seq_orig[nucleo])
                    get_fasta(nucleo,protein)
                    newfasta = open('%s.fasta'%(nucleo),'w+')
                    newfasta.writelines(get_fasta(nucleo,protein))
                    newfasta.close()
                
                        
            elif option == 'S' or option == 's':
                filename = raw_input("\n Please input the name of the file to save: ")
                newfasta = open (filename, 'w+')
                for nucleo in seq_orig:
                    protein = translate(seq_orig[nucleo])
                    get_fasta(nucleo,protein)
                    newfasta.writelines(get_fasta(nucleo,protein))
                newfasta.close()
                
                conclusion = "\n File %s.fasta successfully saved.\n" %(filename)
                print conclusion
            
            else:
                option = raw_input('CAUTION!!! Press S/s or M/m ONLY: ')
                
       # ####################################################################################         
        elif DNAorProtein == 'D' or DNAorProtein == 'd':
            
            for k in seq:
                seq_orig[k] = seq[k].replace('-','')     #seq is the sequence after alignment i.e. without any dashes introduced
                seq_orig_len = seq_orig_len + len(seq_orig[k])   # seq_orig is the sequence before alignement i.e. with dashes introduced
                seq_len = len(seq[k])
                seqno += 1
            if option == 'M' or option == 'm':
                for k in seq:
                    get_fasta(k,seq_orig[k])    
        
                    newfasta = open('%s.fasta'%(k),'w+')
                    newfasta.writelines(get_fasta(k,seq_orig[k]))
                    newfasta.close()
                
            
            elif option == 'S' or option == 's':
                filename = raw_input("\n Please input the name of the file to save: ")
                           
                newfasta = open(filename, 'w+')
                for k in seq:
                    get_fasta(k,seq_orig[k])
                    newfasta.writelines(get_fasta(k,seq_orig[k]))
                newfasta.close()
                conclusion = "\n File %s.fasta successfully saved.\n" %(filename)
                print conclusion 
                            
                
            else:
                option = raw_input('CAUTION!!! Press S/s or M/m ONLY: ')
        else:
            DNAorProtein = raw_input('CAUTION!!! Press D/d or P/p ONLY): ')
            