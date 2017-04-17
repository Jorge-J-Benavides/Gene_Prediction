import random
import numpy.random

# function I stole from http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
# this makes a random choice with weights
def weighted_choice(choices):
   total = sum(w for c, w in choices.items())
   r = random.uniform(0, total)
   upto = 0
   for c, w in choices.items():
      if upto + w >= r:
         return c
      upto += w
   assert False, "Shouldn't get here"


# various variables should be set
filename = "test_viterbi.txt"
# lenght of the file you want
length = 10
# probability of transitioning between states
gene_to_notgene = 0.5
notgene_to_gene = 0.5

# probability of emiting any particular letter in the different states
gene_emits = {"A":0.3, "C":0.2,"T":0.3,"G":0.2}
notgene_emits = {"A":0.2, "C":0.3,"T":0.2,"G":0.3}

# randomly initialize state, G=gene reigon, N=nongene 
states = ["G","N"]
current_state = random.choice(states)

# loop to generate the file
with open(filename,"w") as fl:
    for _ in range(length):
        # randomly choose a letter based on current state:
        if current_state == "G":
            letter = weighted_choice(gene_emits)
        else: # state is nongene
            letter = weighted_choice(notgene_emits)
        # write the state and letter to file
        fl.write(letter+'\t'+current_state+'\n')
        # randomly change states
        if current_state=="G": 
            if numpy.random.binomial(1,gene_to_notgene)==1: current_state="N"
        else:
            if numpy.random.binomial(1,notgene_to_gene)==1: current_state="G"
