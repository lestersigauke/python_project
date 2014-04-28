#1. User to indicate path for file to be analyzed
#2. Manage errors to the given path (e.g. File not found)
#3. Data should be loaded into memory
#4. Error messages should be displayed in case the file does not follow the ClustalW Format. 

##############################################
# Question1 and Question 2
import sys
import os.path
import re
    
def open_alignment(x):    
    
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
        print "%-16s %s"%("File name: ",x)
        print '\n'
        print "%-16s" %("Sequences: "),
        for k in seq:
            seq_orig[k] = seq[k].replace('-','')     #seq is the sequence after alignment i.e. without any dashes introduced
            seq_orig_len = seq_orig_len + len(seq_orig[k])   # seq_orig is the sequence before alignement i.e. with dashes introduced
            seq_len = len(seq[k])
            seqno += 1
            #print seq_orig[k]
            print k,
        print '\n'
        res_collection = ()
        for i in range(0,seq_len):
            res_tuple = ""
            
            for k in seq:
                res_tuple = res_tuple + seq[k][i]
            #print res_tuple
            res_collection = res_collection + (res_tuple,)
        A_count = 0
        G_count = 0
        T_count = 0
        C_count = 0
        n = {}
        
        for l in range(0,seqno+1):
            n[l]= 0 
            
        for k in range(0,seq_len):
            res_collection[k]
            A_count = res_collection[k].count('A',0,seqno)
            C_count = res_collection[k].count('C',0,seqno)
            G_count = res_collection[k].count('G',0,seqno)
            T_count = res_collection[k].count('T',0,seqno)
        
            for l in range(0,seqno+1): 
                if A_count == l:
                    n[l] = n[l] + 1
                elif C_count == l:
                    n[l] = n[l] + 1
                elif G_count == l:
                    n[l] = n[l] + 1   
                elif T_count == l:
                    n[l] = n[l] + 1
        
        print "%-16s"%("Number of matches: ") 
        print '\n',
        for l in range(2,seqno+1):
            print "\t\t[%02i]="%(l), n[l]
                    
        print '\n',                
        #print '\n',res_collection[5]
        
        print "%-16s %i"%("Seq_origi_len: ", seq_orig_len)
        print '\n',
        average_seqlen = float(seq_orig_len/seqno)
        print "%-16s %.2f"%('Average length: ', average_seqlen)
        print '\n',    
        number_of_matches = len(matches)
        match_percent = (float(number_of_matches)/seq_len)*100
        print "%-16s %.2f%%" %("% of Matches: ",match_percent) 
    else:
        proceed = "There is an error with the format of your input file" 
        print proceed

#open_alignment('example.aln')