# ###Question 2###
#a. extract the sequences into separate strings
#b. give them indexes
#c. print them out in Clustal format

import sys
import os.path
import re

def seq_segment(start_seq,start,end):
    start_seq = 'example.aln'
    m = open(start_seq,'r')            #open file
    firstline = m.readline()    #read the first line of FASTA file
    format_checker = "CLUSTAL multiple sequence alignment by MUSCLE (3.8)"   #check whether the file format is consistent with Clustal
    lines = m.readlines()
    seq = {}
    segment = {}
    
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
                    
    
                    segment[line[:line.index(sequence[1])]] = line[line.index(sequence[1]):-1] 
                    if o == None:
                        pass
                    else:
                        counter = 1      # change counter to 1 to show that now we are no longer creating a dictionary but modifying values of the dictionary
                        
                
                elif counter == 1:    # now we are modifying values of the dictionary
                    
                    segment[line[:line.index(sequence[1])]] = segment[line[:line.index(sequence[1])]] + line[line.index(sequence[1]):]
        
        
        print "\n"
        print "CLUSTAL segment [%i-%i] of the %s alignment."%(start,end,start_seq)
        print "\n"
        
        if start < end:
            for k in segment:
                segment[k] = segment[k].replace('\n','')
                 
                c = segment[k][start-1:end]
                string = k
                switch = k
                switch = ''
                stitch = 0
                
                for i in range(1,len(c)+1):
                    width = 60
                    if stitch <= 60:
                        string = string + c[stitch]
                        stitch = stitch + 1 
                        
                print string
                
            print '\n\n',
            if len(c) > 60:
                for k in segment:
                    segment[k] = segment[k].replace('\n','')
                     
                    c = segment[k][start-1:end]
                    strong = k
                    switch = k
                    switch = ''
                    
                    
                    for i in range(61,len(c)):
                        
                        if 61 <= stitch:
                            strong = strong + c[i]
                            stitch = stitch + 1 
                    
                    print strong
