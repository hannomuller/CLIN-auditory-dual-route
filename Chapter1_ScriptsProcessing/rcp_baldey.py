#! usr/bin/python
# rcp_baldey.py
# H. Muller
# 2019-04-29

# Input: Baldey, sep=' '
# Output: Baldey with additional column for root completion point (rcp)

#-------------------------------------------------
# Import modules
import argparse, sys, re, os

#-------------------------------------------------
# Input / Output
def inputfile(infile):
    with open(infile, 'r') as myfile:
        text = myfile.read().splitlines()
    return(text)

def outputfile(outfile, text):
    with open(outfile, 'w') as myfile:
        myfile.write(text)
    return

#-------------------------------------------------
# Parse arguments
def parseargs():

    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile", required=True,
	help="specify the filename")
    ap.add_argument("-o", "--outdir", required=True,
	help="specify the outdir")

    infile = vars(ap.parse_args())["infile"]
    outdir = vars(ap.parse_args())["outdir"]

    return(infile, outdir)

#-------------------------------------------------
# Determine rcp
def determine_rcp(baldey):

    # store results in rcp (root completion pint)
    rcp = []

    # get headers
    header = baldey[0].split(' ')

    # determine column 'stem' and 'word'
    for i,j in enumerate(header):
        if j == 'stem':
            stem_index = i
        if j == 'word':
            word_index = i

    # get stems and words
    stems = [item.split(' ')[stem_index] for item in baldey]
    stems = stems[1:]
    words = [item.split(' ')[word_index] for item in baldey]
    words = words[1:]

    # remove '+' in stems
    stems = [re.sub(r'\+', '', word) for word in stems]

    # for every word
    for i in range(len(words)):

        # preprocess in order to smooth vowels, 'v', 'f', 's', 'z'
        stem_simple = re.sub('s', 'z', stems[i])
        stem_simple = re.sub('f', 'v', stem_simple)
        stem_simple = re.sub(r'([aeiou])\1+', r'\1', stem_simple)

        word_simple = re.sub('s', 'z', words[i])
        word_simple = re.sub('f', 'v', word_simple)
        word_simple = re.sub(r'([aeiou])\1+', r'\1', word_simple)

        # create search pattern like prefix(es) + stem (+ suffix)
        my_regex = str('.+(' + stem_simple + '){1}')
        r = re.compile(my_regex)

        # if stem not in word, word is opaque
        if stem_simple not in word_simple:
            stem_new = 'opaque'
            rcp.append(stem_new)

        # else check for prefixes
        elif r.match(word_simple) != None:

            # get prefix(es)
            my_regex = str('.+(?=' + stem_simple + ')')
            r = re.compile(my_regex)
            prefix = r.findall(word_simple)[0]

            # combine prefix(es) + stem to full stem
            stem_new = str(prefix + stems[i])

            # append full stem (incl. prefix) to rcp
            rcp.append(stem_new)

        # else append normal stem
        else:
            rcp.append(stems[i])

    return(rcp)

#-------------------------------------------------
# Main caller
def main():

    # parse arguments
    infile, outdir = parseargs()

    # read files
    baldey = inputfile(infile)

    # create outfilename
    filename = infile.split('.')[0]
    outname = str(filename + '_only_rcp.txt')
    outfile = os.path.join(outdir, outname)

    # determine rcp
    rcp = determine_rcp(baldey)

    # write new lexicon
    text = '\n'.join(rcp)
    text = text + '\n'
    outputfile(outfile, text)

    print('rcp written to: ' + outfile)

    return

#-------------------------------------------------

main()

#--- EOF -----------------------------------------
