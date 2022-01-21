#! usr/bin/python
# rcp_baldey_integrate.py
# H. Muller
# 2019-05-07

# Input: baldey txt file (sep = space), baldey rcp file (sep = tab)
# Output: baldey txt file with additional column for rcp

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
	help="specify the filename (baldey)")  
    ap.add_argument("-r", "--rcpfile", required=True,
	help="specify the rcp filename") 
    ap.add_argument("-o", "--outfile", required=True,
	help="specify the outfile")

    infile = vars(ap.parse_args())["infile"]
    rcpfile = vars(ap.parse_args())["rcpfile"]
    outfile = vars(ap.parse_args())["outfile"]

    return(infile, rcpfile, outfile) 

#-------------------------------------------------
# Combine baldey and rcp file
def transform_rcp(baldey, rcp):

    # get headers
    header = baldey[0].split(' ')

    # determine column transcription, lip and fip
    for i,j in enumerate(header):
        if j == 'transcription':
            transcription_index = i
        if j == 'lip':
            lip_index = i
        if j == 'fip':
            fip_index = i

    # get transcriptions, lip and fip
    transcription = [item.split(' ')[transcription_index] for item in baldey]
    transcription = transcription[1:]
    lip = [item.split(' ')[lip_index] for item in baldey]
    lip = lip[1:]
    fip = [item.split(' ')[fip_index] for item in baldey]
    fip = fip[1:]

    # get rcp and its transcription
    rcp_ort = [item.split('\t')[0] for item in rcp]
    rcp_raw = [item.split('\t')[1] for item in rcp]
    rcp_raw = [re.sub(' ', '', item) for item in rcp_raw]  

    # make rcp_phon and transcription all lowercase
    rcp_phon_low = [item.lower() for item in rcp_raw]
    trans_low = [item.lower() for item in transcription]

    # shift some characters like g -> x
    rcp_phon_low = [re.sub('g', 'x', item) for item in rcp_phon_low] 
    rcp_phon_low = [re.sub('b', 'p', item) for item in rcp_phon_low] 
    rcp_phon_low = [re.sub('s', 'z', item) for item in rcp_phon_low] 
    rcp_phon_low = [re.sub('t', 'd', item) for item in rcp_phon_low] 
    rcp_phon_low = [re.sub('f', 'v', item) for item in rcp_phon_low] 

    trans_low = [re.sub('g', 'x', item) for item in trans_low] 
    trans_low = [re.sub('b', 'p', item) for item in trans_low] 
    trans_low = [re.sub('s', 'z', item) for item in trans_low] 
    trans_low = [re.sub('t', 'd', item) for item in trans_low] 
    trans_low = [re.sub('f', 'v', item) for item in trans_low] 

    # 1st round: for every word in list, optimize transcription
    for i, word in enumerate(rcp_phon_low):

        # no transcription for opaque words
        if word == 'OpAk':
            pass

        # if word in transcription, keep word
        elif word in trans_low[i]:
            pass

        # else correct transcription based on pattern search
        else:
            correction = pattern(word, trans_low[i])
            rcp_phon_low[i] = correction
    my_dic = {}
    # 2nd round: for every word in list, optimize transcription based on frequent substitutions
    for i, word in enumerate(rcp_phon_low):

        # apply substitutions
        if word not in trans_low[i] and word != 'opak':
            word = re.sub('ui', 'l', word)
            word = re.sub('au', 'mw', word)
            word = re.sub('ei', 'k', word)
            word = re.sub('dz', 'x', word)
            word = re.sub(r'@n\b', '@', word)

            # correct last words manually
            if word not in trans_low[i]:
                if trans_low[i] == 'drinxpar':
                    word = 'drinx'
                elif trans_low[i] == 'merxpar':
                    word = 'merx'
                elif trans_low[i] == 'wimpar':
                    word = 'win'
                elif trans_low[i] == 'zpaxpar':
                    word = 'zpe'
                elif trans_low[i] == 'dr}kz':
                    word = 'dr}kz'
                elif trans_low[i] == 'prexpar':
                    word = 'prex'
                elif trans_low[i] == 'p@mox@':
                    word = 'p@mok'
                elif trans_low[i] == 'v@rm}rw':
                    word = 'v@rm}rw'
                elif trans_low[i] == 'kl}pvodlz':
                    word = 'kl}pvodl'
                elif trans_low[i] == 'pywd@':
                    word = 'pyw'
                elif trans_low[i] == 'kakm':
                    word = 'kakm'
                elif trans_low[i] == 'plmkzxk':
                    word = 'plmkzxk'
                elif trans_low[i] == 'ep@rz':
                    word = 'ep@'
                elif trans_low[i] == 'ymprin':
                    word = 'ymprin'
                elif trans_low[i] == 'werxpar':
                    word = 'werx'
                elif trans_low[i] == 'nld@':
                    word = 'nld@'
                elif trans_low[i] == 'dompar':
                    word = 'dom'
                elif trans_low[i] == 'ondhix':
                    word = 'ondhix'
                elif trans_low[i] == 'maxpar':
                    word = 'max'
                elif trans_low[i] == 'drmp@r@x':
                    word = 'drmp@'
                elif trans_low[i] == 'mokaz)':
                    word = 'mokaz)'
                elif trans_low[i] == 'rexpar':
                    word = 'rex'
                elif trans_low[i] == 'denxpar':
                    word = 'denx'

            # apply changes to rcp_phon_low
            rcp_phon_low[i] = word

    # backtransform substitutions and recreate cases
    rcp_phon = rcp_raw[:]
    for i, word in enumerate(rcp_phon_low):
        
        # exclude marked words
        if word == 'opak':
            rcp_phon[i] = 'opaque'

        # get as many chars from transcription as word is long
        else:        
            chars = len(word)
            word = transcription[i][:chars]
            rcp_phon[i] = word

    # insert headers
    rcp_phon.insert(0, 'rcp')
    rcp_raw.insert(0, 'rcp_raw')
    rcp_ort.insert(0, 'rcp_orthography')

    # combine new columns
    new_cols = []
    for i in range(len(rcp_phon)):
        line = str(rcp_ort[i] + ' ' + rcp_raw[i] + ' ' + rcp_phon[i])
        new_cols.append(line)
 
    return(new_cols)

#-------------------------------------------------
# search the longest matching pattern between word and transcription strating from the end
def pattern(word, transcript):

    match = []
    # get end of word iteratively
    for j in range(len(word)):

        # iteratively increase sequence
        j += 1
        my_regex = str('.+(' + word[-j:] + ').*')
        r = re.compile(my_regex)

        # increase pattern until it doesn't match anymore
        if r.match(transcript) != None:
            match = r.findall(transcript)[-1]

        else:
            break

    # if a pattern could be found, add to rcp_phon_opt
    if match != []:
        my_regex = str('.*' + match)
        r = re.compile(my_regex)
        target = r.findall(transcript)[0]
   
    # else add 'wrong' transcription; will be filtered later
    else:
        target = word

    return(target)

#-------------------------------------------------
# Combine BALDEY and RCP data
def combine(baldey, new_cols):

    # get headers
    header = baldey[0].split(' ')

    # determine fip.ms index
    for i,j in enumerate(header):
        if j == 'fip.ms':
            fipms_index = i

    # insert new columns after fip.ms
    baldey_rcp = []
    for i, line in enumerate(baldey):
        first = line.split(' ')[:fipms_index+1]
        first = ' '.join(first)
        second = line.split(' ')[fipms_index+1:]
        second = ' '.join(second)
        newline = str(first + ' ' + new_cols[i] + ' ' + second)
        baldey_rcp.append(newline)

    return(baldey_rcp)

#-------------------------------------------------
# Main caller
def main():

    # parse arguments
    infile, rcpfile, outfile = parseargs()

    # read files
    baldey = inputfile(infile)
    rcp = inputfile(rcpfile)

    # optimize rcp transcription
    new_cols = transform_rcp(baldey, rcp)

    # make into one table
    baldey_rcp = combine(baldey, new_cols)

    # write new lexicon
    text = '\n'.join(baldey_rcp)
    text = text + '\n'
    outputfile(outfile, text)

    print('rcp written to: ' + outfile)

    return

#-------------------------------------------------

main()

#--- EOF -----------------------------------------
