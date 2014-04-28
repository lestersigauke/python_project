import openmultiplealignment
import Seq_in_clustal_format
import sequenceinfo
import glycosylationsignature
import exporttofasta

#Question1 and 2 ################################################
input_sequence = raw_input('Please enter the path: ')
openmultiplealignment.open_alignment(input_sequence)

#Question3       ################################################
y = input("Enter the start of segment: ")
z = input("Enter the end of segment: ")
Seq_in_clustal_format.seq_segment(input_sequence,y,z)

#Question4       ################################################

dictionary = openmultiplealignment.open_alignment
seqid = raw_input("Please enter the sequence ID: ")
sequenceinfo.seq_info(input_sequence,seqid)

#Question5        ###############################################

glycosylationsignature.glyco_signat(input_sequence)

#Question6        ###############################################

exporttofasta.export_to_fasta(input_sequence)


