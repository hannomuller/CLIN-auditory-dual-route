#! usr/bin/python
# baldey_morphology.py
# H. Muller
# 2021-03-14

# Input: celex and baldey
# Output: baldey enriched with morhological information

#-------------------------------------------------
# Import modules
import argparse, sys, re
import pandas as pd
import numpy as np

#-------------------------------------------------
# Parse arguments
def parseargs():

    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--baldey", required=True,
	help="specify the BALDEY file") 
    ap.add_argument("-c", "--celex", required=True,
	help="specify the CELEX database") 
    ap.add_argument("-o", "--outdir", required=True,
	help="specify the outdir")

    baldey = vars(ap.parse_args())["baldey"]
    celex = vars(ap.parse_args())["celex"]
    outdir = vars(ap.parse_args())["outdir"]

    return(baldey, celex, outdir) 

#-------------------------------------------------
# compute relative frequency for suffixed items in baldey
def get_relfreq(baldey, celex):

    # celex remove duplicate entries
    celex = celex.sort_values(by=['word', 'cob'], ascending=True)
    celex = celex.drop_duplicates(subset=['word'], keep='first')

    # merge words with celex information
    words = baldey[['word', 'stem']]
    words = words.merge(celex[['word', 'composition', 'cob', 'stem_cob']], on='word', how='left')
    words = words.drop_duplicates()

    # merge frequency and morphological information with baldey again
    baldey = baldey.merge(words[['word', 'cob', 'stem_cob']], on='word', how='left')

    return(baldey)        

#-------------------------------------------------
# compute family size for items in baldey
def get_family(baldey, celex):

    celex['transcription'] = celex['transcription'].str.replace("'","").str.replace("-","")

    ### CAN PROBABLY BE ADAPTED SO THAT ONLY WORDS ABOVE CERTAIN FREQUENCY THRESHOLD ARE INCLUDED ###

    ### whole family
    families = celex['stem'].value_counts()
    families = families.to_frame()
    families.rename(columns = {'stem':'family_size'}, inplace=True)
    families['stem'] = families.index

    # merge family information with baldey
    baldey = pd.merge(baldey, families, on='stem', how='left')

    ### onset-aligned family
    # get CELEX stems' transcription
    forms = celex[['word', 'transcription', 'stem', 'morph1', 'morph2', 'morph3', 'morph4', 'morph5', 'morph6']]
    #forms = forms.drop_duplicates(['word'])

    # aggregate possible transcriptions in list
    #forms = celex[['word', 'transcription']].drop_duplicates(['word', 'transcription'])
    #forms = forms.groupby(['word']).agg({'transcription': lambda x: x.tolist()}).reset_index()

    # get baldey items' transcription
    stems = forms[['word', 'transcription']]
    stems.rename(columns = {'word':'stem', 'transcription':'stem_transcription'}, inplace=True)
    my_baldey = baldey.merge(stems, on='stem', how='left')
    
    # get onset aligned family members and continuation forms
    family_oa = pd.DataFrame(columns=['word', 'family_size_oa', 'morph_continuations', 'continuations'])
    items = my_baldey[['word', 'transcription', 'stem']]
    items = items.merge(forms[['word', 'morph1', 'morph2', 'morph3', 'morph4', 'morph5', 'morph6']], on='word')
    items = items.drop_duplicates(['word'])
    items = items[items['word'].notna()]

    for i, row in items.iterrows():
        my_word = row['word']
        my_trans = row['transcription']
        my_stem = row['stem']

        # continuations
        continuations = forms.loc[(forms['transcription'].str.startswith(my_trans, na=False)), ['word', 'transcription']]
        continuations = continuations.drop_duplicates('word')

        # morphological related continuations
        morph_continuations = forms.loc[(forms['stem'] == my_stem) & (forms['transcription'].str.startswith(my_trans, na=False)), ['word', 'transcription']]
        morph_continuations = morph_continuations.drop_duplicates('word')

        # onset aligned family members  
        if (row['morph1'] == my_stem) or (row['morph1']==False):
            oa_members = forms.loc[(forms['morph1'] == my_stem), ['word', 'transcription']]
            oa_members = oa_members.drop_duplicates('word')
        elif row['morph2'] == my_stem:
            oa_members = forms.loc[((forms['morph2'] == my_stem) & (forms['morph1'] == row['morph1'])), ['word', 'transcription']]
            oa_members = oa_members.drop_duplicates('word')
        else:
            oa_members = []   
    
        ### NOTE: 'pious' has member 'pieties' since shared first stem regardeless of pronounciation. Should be improved


        ''' use if variabilty in transcription should be reflected in family size
        my_family = pd.DataFrame()
        for j in my_trans:
            match = celex.loc[(celex['stem'] == my_stem) & (celex['transcription'].str.startswith(j, na=False)), 'word']
            my_family = my_family.append(match)
        '''

        continu = (len(continuations))
        morph_continu = (len(morph_continuations))
        family_size_oa = (len(oa_members))
        family_oa = family_oa.append({'word': my_word, 'family_size_oa': family_size_oa, 'continuations':continu, 'morph_continuations':morph_continu}, ignore_index=True)

    baldey = baldey.merge(family_oa, on='word', how='left')

    ### ONSET-ALIGNED FAMILY SIZE MIGHT NOT BE CORRECT FOR IRREGULAR WORDS, OR WORDS WITH SEVERAL SUFFIXES ###

    return(baldey)

#-------------------------------------------------
# Main caller
def main():

    # parse arguments
    baldey, celex, outdir = parseargs()

    # read files
    baldey = pd.read_csv(baldey, sep='\t', header=0)
    #baldey = baldey[:20] # use for faster development
    #baldey = baldey.drop(['family_size', 'family_size_oa'],axis=1)
    celex = pd.read_csv(celex, sep='\t', header=0)

    # compute relative frequency
    baldey = get_relfreq(baldey, celex)

    # compute family size
    #baldey = get_family(baldey, celex)  

    # save file to disk
    #baldey = baldey[['word', 'family_size', 'family_size_oa', 'continuations', 'morph_continuations', 'cob', 'stem_cob']]
    baldey = baldey[['word', 'cob', 'stem_cob']]
    baldey.fillna('-', inplace=True)
    outfile = 'my_baldey_morphology.txt'
    baldey.to_csv(outfile, sep='\t')

    print('Family size and stem frequencies written to: ' + outfile)

    return

#-------------------------------------------------
main()

#--- EOF -----------------------------------------

