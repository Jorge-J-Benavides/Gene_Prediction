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

def precision(correct,guesses):
    # precision = number of true positives / number of total positives
    true_positives = 0
    for i in range(len(correct)):
        if correct[i]=="G" and guesses[i]=="G": true_positives = true_positives+1
    return true_positives *1.0 / len([g for g in guesses if g=="G"])

