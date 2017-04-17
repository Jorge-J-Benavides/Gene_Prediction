#Take in seq of states and return the postion where the states switch, because it could be a possible codon


seq = 'GGGGNNNGGG'

for i in range(0,len(seq)-2):
    
    if seq[i] != seq[i+1]:
        print i