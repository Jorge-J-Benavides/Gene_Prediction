from collections import defaultdict
import itertools

#read in all lines
with open('test_viterbi.txt','r') as fl:
    lines = [l.strip() for l in fl.readlines()]

# add sentences as lists of lines, split by newlines
sentences = [[]]

for line in lines:
    if line == '':
        sentences.append([])
    else:
        sentences[-1].append(line)

# convert sentence lists to (word,tag) tupples
word_tags = [[(l.split('\t')[0], l.split('\t')[1]) for l in sentence] for sentence in sentences]

#insert artificial start and end tags
for sentence in word_tags:
    sentence.insert(0,('','s'))
    sentence.append(('','\s'))

# get a list of all tags, appending the start/end because its easier to loop over lines and not sentences
tags = list(set([l.split('\t')[1] for l in lines if l!= '']))
tags.append('s')
tags.append('\s')

# get counts for all tags
tag_counts = defaultdict(int)
for s in word_tags:
    for w in s:
        tag_counts[w[1]] += 1

# now calculate transition probabilities and emmission probabilities
# first, initialize counts to 0 via defaultdict
transcounts = {t:defaultdict(int) for t in tags}
emitcounts = {t:defaultdict(int) for t in tags}

# loop over all sentences and increment the correct transmission and emmision count
for s in word_tags:
    for i in range(len(s)):
        # indexing works like: transcounts[current tag][next tag]
        # the if statement avoids this for '\s', since it doens't transition to anything
        if i < len(s) - 1: transcounts[s[i][1]][s[i+1][1]] += 1
        # indexing works like: emitcounts[current tag][current word]
        emitcounts[s[i][1]][s[i][0]] += 1

# the actual algorithm
def viterbi(s,trans,emits,counts):
    # insert '' at begining and end for start and end tags
    s.insert(0,'')
    s.append('')
    # initialize the trellis to zeros except for start
    trellis = [{t:0 for t in tags} for _ in range(len(s))]
    trellis[0]['s'] = 1
    # initialize the backpointer
    backpointer = [{t:'' for t in tags} for _ in range(len(s)-1)]
    # do the actual viterbi loop
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

# now run the sentences as specified in the description
# note that my implimentation required the use of the tag counts
#print(viterbi(['This','is','a','sentence','.'], transcounts, emitcounts,tag_counts))
#print(viterbi(['This','might','produce','a','result','if','the','system','works','well','.'], transcounts, emitcounts,tag_counts))
print(viterbi(['N','N','N','N','G','G','G'], transcounts, emitcounts,tag_counts))
#print(viterbi(['Can','a','can','move','a','can','?'], transcounts, emitcounts,tag_counts))
#print(viterbi(['Can','you','walk','the','walk','and','talk','the','talk','?'], transcounts, emitcounts,tag_counts))
