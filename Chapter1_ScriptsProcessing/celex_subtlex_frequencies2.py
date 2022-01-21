#! usr/bin/python
# celex_subtlex_frequencies.py
# H. Muller
# 2020-11-04

# Input: celex lemma and form file and subtlex
# Output: lemma, stem and form frequency of every word in both databases

#-------------------------------------------------
# Import modules
import argparse, sys, re, os
#import collections
import pandas as pd

#-------------------------------------------------
# Input / Output
def inputfile(infile):
    with open(infile, 'r') as myfile:
        text = myfile.read().splitlines()
    return(text)

#-------------------------------------------------
# Parse arguments
def parseargs():

    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--baldey", required=True,
	help="specify the BALDEY filename") 
    ap.add_argument("-l", "--l_celex", required=True,
	help="specify the CELEX orthographic lemma filename") 
    ap.add_argument("-f", "--f_celex", required=True,
	help="specify the CELEX orthographic form filename")
    ap.add_argument("-i", "--lf_celex", required=True,
	help="specify the CELEX frequency lemma filename") 
    ap.add_argument("-j", "--ff_celex", required=True,
	help="specify the CELEX frequency form filename")
    ap.add_argument("-p", "--p_celex", required=True,
	help="specify the CELEX phonetic lemma filename") 
    ap.add_argument("-q", "--q_celex", required=True,
	help="specify the CELEX phonetic form filename")
    ap.add_argument("-s", "--subtlex", required=True,
	help="specify the CELEX morphemic lemma filename") 
    ap.add_argument("-o", "--outdir", required=True,
	help="specify the outdir")

    baldey = vars(ap.parse_args())["baldey"]
    lemmas_ort = vars(ap.parse_args())["l_celex"]
    form_ort = vars(ap.parse_args())["f_celex"]
    lemmas_freq = vars(ap.parse_args())["lf_celex"]
    form_freq = vars(ap.parse_args())["ff_celex"]
    lemmas_phon = vars(ap.parse_args())["p_celex"]
    form_phon = vars(ap.parse_args())["q_celex"]
    subtlex = vars(ap.parse_args())["subtlex"]
    outdir = vars(ap.parse_args())["outdir"]

    return(baldey, lemmas_ort, form_ort, lemmas_freq, form_freq, lemmas_phon, form_phon, subtlex, outdir) 

#-------------------------------------------------
# Preprocess CELEX
def preprocess(lemmas_ort, form_ort, lemmas_freq, form_freq, lemmas_phon, form_phon, subtlex):

    # make pandas from relevant CELEX lemma columns
    form = [item.split('\\')[1] for item in form_phon]
    form_lemma = [item.split('\\')[3] for item in form_phon]
    form_trans = [item.split('\\')[4] for item in form_phon]
    form_trans = [item.replace("'", '').replace('-', '') for item in form_trans]
    form_freq = [item.split('\\')[3] for item in form_freq]
    form_db = pd.DataFrame(list(zip(form, form_lemma, form_trans, form_freq)), columns =['word', 'index', 'transcription', 'CELEXform_freq']) 

    # make pandas from relevant CELEX lemma columns
    index = [item.split('\\')[0] for item in lemmas_ort]
    word = [item.split('\\')[1] for item in lemmas_ort]
    freq = [item.split('\\')[2] for item in lemmas_freq]
    #composition = [item.split('\\')[8] for item in lemmas_ort]
    #morphs = [item.split('\\')[12] for item in lemmas_ort]
    trans = [item.split('\\')[3] for item in lemmas_phon]
    trans = [item.replace("'", '').replace('-', '') for item in trans]
    lemma_db = pd.DataFrame(list(zip(index, word, trans, freq)), columns =['index', 'lemma', 'lemma_trans', 'CELEXlemma_freq']) 

    # use for faster development
    #form_db = form_db.head(10000)
    #lemma_db = lemma_db.head(1000)

    ### get stems and lemma frequencies
    form_db = pd.merge(form_db, lemma_db, how='left', on=['index'])

    print('CELEX preprocessed')
    return(form_db)

#-------------------------------------------------
# Retrieve frequency for every wordform
def stem_to_form(my_stem, direction):

    old_stem = my_stem
    alternation = 0

    # form to stem
    if direction == 'from_form':
        if my_stem[-1] == 'G':
            my_stem = my_stem[:-1] + 'x'
        if my_stem[-1] == 'v':
            my_stem = my_stem[:-1] + 'f'
        if my_stem[-1] == 'd':
            my_stem = my_stem[:-1] + 't'
        if my_stem[-1] == 'b':
            my_stem = my_stem[:-1] + 'p'
        if my_stem[-1] == 'z':
            my_stem = my_stem[:-1] + 's'
        if my_stem[-1] == 'g':
            my_stem = my_stem[:-1] + 'k'

    # stem to form
    if direction == 'from_stem':
        if my_stem[-1] == 'x':
            my_stem = my_stem[:-1] + 'G'
        if my_stem[-1] == 'f':
            my_stem = my_stem[:-1] + 'v'
        if my_stem[-1] == 't':
            my_stem = my_stem[:-1] + 'd'
        if my_stem[-1] == 'p':
            my_stem = my_stem[:-1] + 'b'
        if my_stem[-1] == 's':
            my_stem = my_stem[:-1] + 'z'
        if my_stem[-1] == 'k':
            my_stem = my_stem[:-1] + 'g'

    # save whether alternation takes place
    if old_stem != my_stem:
        alternation = 1

    return(my_stem, alternation)

#-------------------------------------------------
# Retrieve frequency for every wordform
def forms_retrieve(baldey, ort, trans, my_db):

    altcol = ort + '-alternation'
    baldey[altcol] = 0
    baldey[ort] = 0

    # find words with frequencies in CELEX
    # preprocess my_baldey
    my_affix = baldey.loc[baldey['one_compound'] == 0,]
    my_affix = my_affix.loc[my_affix['one_prefix'] == 0,]
    my_affix = my_affix.loc[my_affix['only_suffix'] == ort, ['transcription', 'rcp']]
    my_affix = my_affix.drop_duplicates()
    my_affix = my_affix.reset_index()

    for i in range(len(my_affix)):
        my_form = my_affix.loc[i, 'transcription']
        my_stem = my_affix.loc[i, 'rcp']
        my_stem, alternation = stem_to_form(my_stem, 'from_form')

        # mark row
        baldey.loc[baldey['transcription'] == my_form, ort] = 1

        # find stem frequency
        if my_stem in my_db.transcription.unique():
            CELEXstem_freq = my_db.loc[my_db['transcription'] == my_stem, 'CELEXform_freq']
            CELEXstem_freq = float(CELEXstem_freq.sum())#; print(CELEXstem_freq)
        else:
            CELEXstem_freq = 0
            #print('form processing: ' + my_form + ' ' + my_stem)
            #print(my_stem + ' not found')

    
        # add frequencies to baldey
        baldey.loc[baldey['transcription'] == my_form, 'CELEXstem_freq'] = CELEXstem_freq
        baldey.loc[baldey['transcription'] == my_form, altcol] = alternation

    return(baldey)

#-------------------------------------------------
# Derive forms of stem and retrieve frequencies
def stems_retrieve(baldey, ort, trans, my_db):

    altcol = ort + '-alternation'
    newcol = 'CELEX' + ort + '_form_freq'
    baldey[newcol] = 0

    # find monomorphemic words
    # preprocess my_baldey
    my_stems = baldey.loc[baldey['one_compound'] == 0,]
    my_stems = my_stems.loc[my_stems['one_prefix'] == 0,]
    my_stems = my_stems.loc[baldey['transcription'] == baldey['rcp'], :]
    my_stems = my_stems.rcp.unique()

    for my_stem in my_stems:

        if ort not in ['ig', 's', 't']:
            my_form, alternation = stem_to_form(my_stem, 'from_stem')
            my_form = my_form + trans
        else:
            alternation = 0
            my_form = my_stem + trans

        # find lemma and affixed form frequency
        if my_form in my_db.transcription.unique():
            CELEXform_freq = my_db.loc[my_db['transcription'] == my_form, 'CELEXform_freq']
            CELEXform_freq = float(CELEXform_freq.sum())
            #print('stem-form-pattern: ' + my_stem + ' ' + my_form)
            baldey.loc[baldey['transcription'] == my_stem, ort] = 1

        else:
            CELEXform_freq = 0
    
        # add frequencies to baldey
        baldey.loc[baldey['transcription'] == my_stem, newcol] = CELEXform_freq
        baldey.loc[baldey['transcription'] == my_form, altcol] = alternation

    return(baldey)

#-------------------------------------------------
# Main caller
def general_retrieve(baldey, suffixes, my_db):

    # get stems
    my_stems = baldey.loc[baldey['one_compound'] == 0,]
    my_stems = my_stems.loc[my_stems['one_prefix'] == 0,]
    my_stems = my_stems.loc[baldey['transcription'] == baldey['rcp'], :]
    my_stems = my_stems.transcription.unique()

    # find stem frequency
    for my_stem in my_stems:
        print(my_stem)
        if my_stem in my_db.transcription.unique():
            CELEXstem_freq = my_db.loc[my_db['transcription'] == my_stem, 'CELEXform_freq']
            CELEXstem_freq = float(CELEXstem_freq.sum())
            CELEXform_freq = CELEXstem_freq
            CELEXlemma_freq = my_db.loc[my_db['transcription'] == my_stem, 'CELEXlemma_freq']
            CELEXlemma_freq = float(CELEXlemma_freq.sum())
        else:
            CELEXstem_freq = 'NA'
            CELEXform_freq = 'NA'
            CELEXlemma_freq = 'NA'
        baldey.loc[baldey['transcription'] == my_stem, 'CELEXstem_freq'] = CELEXstem_freq
        baldey.loc[baldey['transcription'] == my_stem, 'CELEXform_freq'] = CELEXform_freq
        baldey.loc[baldey['transcription'] == my_stem, 'CELEXlemma_freq'] = CELEXlemma_freq

    # get suffixed words
    my_suffixes = list(suffixes.keys())
    my_forms = baldey.loc[baldey['one_compound'] == 0,]
    my_forms = my_forms.loc[my_forms['one_prefix'] == 0,]
    my_forms = my_forms.loc[my_forms['only_suffix'].isin(my_suffixes), :]
    my_forms = my_forms.transcription.unique()

    # get form frequencies
    for my_form in my_forms:
        print(my_form)
        if my_form in my_db.transcription.unique():
            CELEXform_freq = my_db.loc[my_db['transcription'] == my_form, 'CELEXform_freq']
            CELEXform_freq = float(CELEXform_freq.sum())
            CELEXlemma_freq = my_db.loc[my_db['transcription'] == my_form, 'CELEXlemma_freq']
            CELEXlemma_freq = float(CELEXlemma_freq.sum())
        else:
            CELEXform_freq = 'NA'
            CELEXlemma_freq = 'NA'
        baldey.loc[baldey['transcription'] == my_form, 'CELEXform_freq'] = CELEXform_freq
        baldey.loc[baldey['transcription'] == my_form, 'CELEXlemma_freq'] = CELEXlemma_freq

    return(baldey)

#-------------------------------------------------
# Main caller
def main():

    # parse arguments
    baldey, lemmas_ort, form_ort, lemmas_freq, form_freq, lemmas_phon, form_phon, subtlex, outdir = parseargs()

    # read files
    baldey = pd.read_csv(baldey, sep='\t', header=0)
    #baldey = baldey[:100]
    lemmas_ort = inputfile(lemmas_ort)
    form_ort = inputfile(form_ort)
    lemmas_freq = inputfile(lemmas_freq)
    form_freq = inputfile(form_freq)
    lemmas_phon = inputfile(lemmas_phon)
    form_phon = inputfile(form_phon)
    subtlex = pd.read_csv(subtlex, sep='\t')
    subtlex = subtlex.drop_duplicates()

    # create outfilename
    outfile = 'my_baldey_frequencies.txt'

    # preprocess databases
    my_db = preprocess(lemmas_ort, form_ort, lemmas_freq, form_freq, lemmas_phon, form_phon, subtlex)

    # retrieve frequencies
    suffixes = baldey[['only_suffix', 'only_suffix_trans']]
    suffixes = suffixes.drop_duplicates()
    suffixes = suffixes.groupby('only_suffix')['only_suffix_trans'].apply(list).to_dict()
    suffixes.pop('-')

    # retrieve stem and form frequencies for monomorphemic words and forms
    baldey = general_retrieve(baldey, suffixes, my_db)
    print('stems\' frequencies retrieved')
    
    # find forms buildt of stem + suffix and retrieve frequencies
    for suffix in suffixes:
        ort = suffix
        trans = suffixes[suffix][0]
        print('*'*30); print('processing\t\'' + ort + '\'\t/' + trans +'/')
        baldey = forms_retrieve(baldey, ort, trans, my_db)
        baldey = stems_retrieve(baldey, ort, trans, my_db)
        n_stems = baldey.loc[baldey['transcription'] == baldey['rcp'], :]
        n_stems = n_stems.loc[baldey[ort] == 1, :]
        n_stems = len(n_stems)
        n_forms = baldey.loc[baldey['only_suffix'] == ort, :]
        n_forms = baldey.loc[baldey[ort] == 1, :]
        n_forms = len(n_forms)
        n_total = n_stems + n_forms
        feedback = 'total: ' + str(n_total) + '\tstems: ' + str(n_stems) + '\tforms: ' + str(n_forms)
        print(feedback); print('*'*30)

    # write new database
    baldey.to_csv(outfile, sep='\t')

    print('Lemma, stem and form frequencies written to: ' + outfile)

    return

#-------------------------------------------------

main()

#--- EOF -----------------------------------------

