# Project 
# Erik Holbrook Jorge Benavides Cosima Jackson

from collections import defaultdict
import argparse
from os import path
import itertools


tags = None
tag_counts = None
################################################################################################################################
#    References:
#        
################################################################################################################################
################################################################################################################################
#    Usage function to help the user type in the correct command to run the file
#    (e.g. python jorgebenavides_hw1.py -F <filename> 
################################################################################################################################
def usage():
    
    print ('\npython viterbi.py -D <DNA_seq_filename> -S <DNA_state_filename> \n')
    print ("where -D is the name of a file to be read and analyzed")
    print ("where -S is the name of a file to be read and analyzed")

################################################################################################################################
#    Takes takes the dna_seq filepath and the dna_state file path and converts the data into a new file called seq_state that 
#    has the correct format of the combined data ->nucleotide    state that the viterbi functions takes in
################################################################################################################################
def get_dna_into_correct_format(dna_seq_filepath, dna_state_filepath):
    
    #open dna_state file and parse it into dna_state_array 
    with open(dna_state_filepath,'r') as fl:
        dna_state_array = [l.strip() for l in fl.readlines()]
        dna_state_array = "".join(dna_state_array)
        
    #open dna_seq file and parse it into dna_seq_array 
    with open(dna_seq_filepath,'r') as fl:
        dna_seq_array = [l.strip() for l in fl.readlines()]
        dna_seq_array = "".join(dna_seq_array)
        
    #name for new state and seq combined in correct format file 
    seq_state = "seq_state.txt"
    
    #create new seq_state file in directory with correct data format for viterbi function
    with open(seq_state,"w") as fl:
        for nucleotide, state in zip(dna_seq_array, dna_state_array):
        #for nucleotide in seq array and state in state array combine both with tab divider
            fl.write(nucleotide + '\t' + state + '\n')

    # return the name of the file created with the seq and state in the correct format
    return seq_state, dna_state_array 

################################################################################################################################

# each of these takes two lists of "G" and "N"
# the first is always the correct list
# the second is the output of viterbi
# note that they should be the same length!

def accuracy(correct,guesses):
    # accuracy = number of correct / number of total
    num_correct = 0
    for i in range(len(correct)):
        if guesses[i]==correct[i]: num_correct = num_correct+1
    return num_correct*1.0 / len(correct)

def sensitivity(correct,guesses):
    # sensitivity = number of true positives caught by guesses
    num_correct = 0
    for i in range(len(correct)):
        if correct[i]=="G" and guesses[i]=="G": num_correct = num_correct+1
    return num_correct * 1.0 / len([c for c in correct if c=="G"])

def specificity(correct,guesses):
    # exact same thing as above except with "N"
    num_correct = 0
    for i in range(len(correct)):
        if correct[i]=="N" and guesses[i]=="N": num_correct = num_correct+1
    return num_correct * 1.0 / len([c for c in correct if c=="N"])

'''
def precision(correct,guesses):
    # precision = number of true positives / number of total positives
    true_positives = 0
    for i in range(len(correct)):
        if correct[i]=="G" and guesses[i]=="G": true_positives = true_positives+1
    return true_positives *1.0 / len([g for g in guesses if g=="G"])
'''


################################################################################################################################
# the actual algorithm
def viterbi(s,trans,emits,counts):
    # insert '' at begining and end for start and end tags
    s.insert(0,'')
    s.append('')
    global tags
    # initialize the trellis to zeros except for start
    trellis = [{t:0 for t in tags} for _ in range(len(s))]
    trellis[0]['s'] = 1
    # initialize the backpointer
    backpointer = [{t:'' for t in tags} for _ in range(len(s)-1)]
    # do the actual viterbi loop
    global tag_counts
    for i in range(1,len(s)):
        for t in tags:
            # calculate all probabilities for this word/tag (makes the backpointer easier)
            # note that prob =  max([tag emits word * tag2 -> tag * tag2]) where tag is the current tag and tag2 is previous
            probs = [(t2,emits[t][s[i]]/tag_counts[t] * trans[t2][t]/tag_counts[t2] * trellis[i-1][t2]) for t2 in tags]
            # store only max in trellis (for forward algorithm, this would be sum)
            trellis[i][t] = max(list(zip(*probs))[1])
            # now store most probable t2 from previous step
            backpointer[i-1][t] = probs[list(zip(*probs))[1].index(trellis[i][t])][0]
    # now we can start from backpointer[len(s)-2]['\s'] and work our way back
    i = len(s)-2
    t = '\s'
    taglist = []
    while i > 0:
        taglist.insert(0,backpointer[i][t])
        t = backpointer[i][t]
        i = i-1
    return taglist



################################################################################################################################
################################################################################################################################
def main():
    #call function to get seq and states into one file in the correct format
    seq_state_filename, correct_state_array = get_dna_into_correct_format(dna_seq_filepath ,dna_state_filepath)

    #read in all lines
    with open(seq_state_filename,'r') as fl:
        lines = [l.strip() for l in fl.readlines()]
        
    # add nucleotide sequence as lists of lines, split by newlines
    nucleotide_state = [[]]

    for line in lines:
        if line == '':
            nucleotide_state.append([])
        else:
            nucleotide_state[-1].append(line)

    # convert nucleotide_state lists to (nucleotide, state) tupples
    word_tags = [[(l.split('\t')[0], l.split('\t')[1]) for l in tuple] for tuple in nucleotide_state]
    
    #insert artificial start and end tags
    for sentence in word_tags:
        sentence.insert(0,('','s'))
        sentence.append(('','\s'))

    # get a list of all tags, appending the start/end because its easier to loop over lines and not nucleotide sequence
    global tags
    tags = list(set([l.split('\t')[1] for l in lines if l!= '']))
    tags.append('s')
    tags.append('\s')

    # get counts for all tags
    global tag_counts
    tag_counts = defaultdict(int)
    for s in word_tags:
        for w in s:
            tag_counts[w[1]] += 1

    # now calculate transition probabilities and emmission probabilities
    # first, initialize counts to 0 via defaultdict
    transcounts = {t:defaultdict(int) for t in tags}
    emitcounts = {t:defaultdict(int) for t in tags}

    # loop over all nucleotide sequence and increment the correct transmission and emmision count
    for s in word_tags:
        for i in range(len(s)):
            # indexing works like: transcounts[current tag][next tag]
            # the if statement avoids this for '\s', since it doens't transition to anything
            if i < len(s) - 1: transcounts[s[i][1]][s[i+1][1]] += 1
            # indexing works like: emitcounts[current tag][current word]
            emitcounts[s[i][1]][s[i][0]] += 1
    
    # now run the nucleotide sequence as specified in the description
    # note that my implimentation required the use of the tag counts
    
    #using loop to generate a guess of states to plug into vitrbi that will be the correct length
    guesses_states_array = []
    for i in range(len(correct_state_array)):
        n = 'N'
        guesses_states_array.append(','.join(n))
    

    viterbi_returned_states = (viterbi(guesses_states_array, transcounts, emitcounts,tag_counts))
    
    print ("Guessed States Plugged into Viterbi" , "".join(guesses_states_array))
    print ("Viterbi Output States              " , "".join(viterbi_returned_states))
    print ("Correct States                     " , correct_state_array)
    print ("Accuracy: ", accuracy(correct_state_array, viterbi_returned_states))
    print ("Sensitivity: ", sensitivity(correct_state_array, viterbi_returned_states))
    print ("Specificity: ", specificity(correct_state_array, viterbi_returned_states))
    '''
    print (precision(correct_state_array, viterbi_returned_states))
    '''
        
################################################################################################################################
# Have to use this line of code because our program needs to run with command line arguments
################################################################################################################################
if __name__ == "__main__":
    #create an instance of argument parser
    parse = argparse.ArgumentParser()
    
    #add arguments to parser so if command line argumnet matches -f or -F or --filename it will take the next argument
    #and put it into the --filename. See reference for more details on argparse.
    parse.add_argument('-D','-d', '--dna_filename')
    parse.add_argument('-S','-s', '--state_filename')
    #create instance of parse_args
    args = parse.parse_args()
    
    #set the args to variable names for easy reference
    dna_seq_filename = args.dna_filename
    dna_state_filename = args.state_filename
    
    #error handling if any variables are empty something went wrong call usage function to inform user
    if(dna_seq_filename == None or dna_state_filename == None):
        usage()
        exit(2)
        
    #get absolute file path for filename    
    dna_seq_filepath = path.abspath(dna_seq_filename)
    dna_state_filepath = path.abspath(dna_state_filename)
    
    seq_state_filename = get_dna_into_correct_format(dna_seq_filepath ,dna_state_filepath)

    main()
