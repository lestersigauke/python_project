#Question 4
#1. Count the number of bases per sequence.

import sys
import os.path
import re

def seq_info(x,yu):
    
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
        print '\n',
        print "%-6s %s"%("ID: ",yu)
        print '\n',
        
        for k in seq:
            seq_orig[k] = seq[k].replace('-','')     #seq is the sequence after alignment i.e. without any dashes introduced
            seq_orig_len = seq_orig_len + len(seq_orig[k])   # seq_orig is the sequence before alignement i.e. with dashes introduced
            seq_len = len(seq[k])
            seqno += 1
        
        print 'Length:', len(seq_orig[yu])
        print "\n",
        print "Frequency per base: "
        count = {}
        string_id = ('AGCT')
        
        count[string_id[0]] = seq_orig[yu].count('A')
        count[string_id[1]] = seq_orig[yu].count('G')
        count[string_id[2]] = seq_orig[yu].count('C')
        count[string_id[3]] = seq_orig[yu].count('T')
        
        print "\n",
        for i in count:
            print '%-19s[%s]:'%("",i),count[i], '\n'
        
        print "Sequence: \n"
        new_seq = ' '
        for i in seq_orig[yu]:
            if len(new_seq)%60 == 0:
                new_seq = new_seq + i + '\n\n '
                
            else:
                new_seq = new_seq + i
        print new_seq
                
                
    