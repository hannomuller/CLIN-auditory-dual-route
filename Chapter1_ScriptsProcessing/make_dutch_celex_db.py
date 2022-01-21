#! usr/bin/python
# make_dutch_celex_db.py
# H. Muller
# 2021-03-14

# Input: celex (DUTCH)
# Output: celex db

### NOTE: In Dutch CELEX, lemmas are not correctly coded for prefixed words. In addition, compounds' lemmas coding is questionable too

#-------------------------------------------------
# Import modules
import argparse, sys, re
import pandas as pd
import numpy as np

#-------------------------------------------------
# Parse arguments
def parseargs():

	ap = argparse.ArgumentParser()
	ap.add_argument("-l", "--celex_lemma", required=True,
	help="specify the CELEX lemma file")
	ap.add_argument("-w", "--celex_form", required=True,
	help="specify the CELEX wordform file")
	ap.add_argument("-f", "--celex_freq", required=True,
	help="specify the CELEX frequency count file")
	ap.add_argument("-p", "--celex_phonetics", required=True,
	help="specify the CELEX phonetic transcription file")
	ap.add_argument("-o", "--outdir", required=True,
	help="specify the outdir")

	lemma = vars(ap.parse_args())["celex_lemma"]
	form = vars(ap.parse_args())["celex_form"]
	freq = vars(ap.parse_args())["celex_freq"]
	phon = vars(ap.parse_args())["celex_phonetics"]
	outdir = vars(ap.parse_args())["outdir"]

	return(lemma, form, freq, phon, outdir)

#-------------------------------------------------
# read a file line by line
def inputfile(filename):
	with open(filename, "r") as text:
		lines = text.readlines()
	return(lines)

#-------------------------------------------------
# explain stuff here
def preprocess_celex(lemma, form, freq, phon):

	# for lemma file
	my_lemma = []
	for line in lemma:
		lemma_index = line.split('\\')[0]
		lemma = line.split('\\')[1]
		morphs = line.split('\\')[8]
		composition = line.split('\\')[12]
		my_line = lemma_index + '\t' + lemma + '\t' + morphs + '\t' + composition
		my_lemma.append(my_line)
	my_lemma = pd.DataFrame([line.split('\t') for line in my_lemma])
	my_lemma.columns = ['lemma_index', 'word', 'morphs', 'composition']

	# for form file
	my_form = []
	for line in form:
		form_index = line.split('\\')[0]
		lemma_index = line.split('\\')[3]
		word = line.split('\\')[1]
		my_line = form_index + '\t' + lemma_index + '\t' + word
		my_form.append(my_line)
	my_form = pd.DataFrame([line.split('\t') for line in my_form])
	my_form.columns = ['form_index', 'lemma_index', 'word']

	# for frequency file
	my_freq = []
	for line in freq:
		index = line.split('\\')[0]
		word = line.split('\\')[1]
		cob = line.split('\\')[3]
		my_line = index + '\t' + word + '\t' + cob
		my_freq.append(my_line)
	my_freq = pd.DataFrame([line.split('\t') for line in my_freq])
	my_freq.columns = ['form_index', 'word', 'cob']

	# for phonetics file
	my_phon = []
	for line in phon:
		index = line.split('\\')[0]
		word = line.split('\\')[1]
		trans = line.split('\\')[4]
		my_line = index + '\t' + word + '\t' + trans
		my_phon.append(my_line)
	my_phon = pd.DataFrame([line.split('\t') for line in my_phon])
	my_phon.columns = ['form_index', 'word', 'transcription']

	# merge lemmas and freq
	my_lemma = pd.merge(my_lemma, my_freq, on='word', how='left')
	my_lemma = my_lemma.rename(columns = {'cob':'lemma_cob', 'word':'lemma'})
	my_lemma.pop('form_index')
	my_lemma.drop_duplicates()

	# merge forms and freq and trans
	my_form = pd.merge(my_form, my_freq, on=['form_index', 'word'], how='left')
	my_form.drop_duplicates()
	my_form = pd.merge(my_form, my_phon, on=['form_index', 'word'], how='left')
	my_form.drop_duplicates()

	# merge lemmas and form
	my_db = pd.merge(my_form, my_lemma, on='lemma_index', how='left')
	my_db[['form_index', 'lemma_index', 'cob', 'lemma_cob']] = my_db[['form_index', 'lemma_index', 'cob', 'lemma_cob']].apply(pd.to_numeric, errors='coerce')

	# remove duplicate entries
	my_db.sort_values(by=['form_index', 'lemma_cob'], inplace=True, ascending=True)
	my_db = my_db.drop_duplicates(subset=['form_index', 'lemma_index', 'word', 'cob', 'lemma', 'morphs', 'composition'], keep='first') # keep last might be plausible too. Where do the differences in freq come from?

	return(my_db)

#-------------------------------------------------
# correct stem
def true_stem(my_db):

	# get seperable lemmas
	mask = my_db['composition'].str.contains(r',', na=False)
	multimorph = my_db.loc[mask,]
	splits = multimorph.composition.str.split(',', expand=True)
	multimorph = pd.concat([multimorph, splits], axis=1)

	# get classes
	multimorph['class1'] = multimorph[0].str.extract(r'\[(.*?)\]')
	multimorph['class2'] = multimorph[1].str.extract(r'\[(.*?)\]')
	multimorph['class3'] = multimorph[2].str.extract(r'\[(.*?)\]')
	multimorph['class4'] = multimorph[3].str.extract(r'\[(.*?)\]')
	multimorph['class5'] = multimorph[4].str.extract(r'\[(.*?)\]')
	multimorph['class6'] = multimorph[5].str.extract(r'\[(.*?)\]')

	# get morphs
	multimorph['morph1'] = multimorph[0].str.extract(r'\(*(.*?)\)')
	multimorph['morph2'] = multimorph[1].str.extract(r'\(*(.*?)\)')
	multimorph['morph3'] = multimorph[2].str.extract(r'\(*(.*?)\)')
	multimorph['morph4'] = multimorph[3].str.extract(r'\(*(.*?)\)')
	multimorph['morph5'] = multimorph[4].str.extract(r'\(*(.*?)\)')
	multimorph['morph6'] = multimorph[5].str.extract(r'\(*(.*?)\)')

	# check if component is monomorphemic word
	multimorph['mask1'] = multimorph['class1'].str.contains(r'\|', na=True)
	multimorph['mask2'] = multimorph['class2'].str.contains(r'\|', na=True)
	multimorph['mask3'] = multimorph['class3'].str.contains(r'\|', na=True)
	multimorph['mask4'] = multimorph['class4'].str.contains(r'\|', na=True)
	multimorph['mask5'] = multimorph['class5'].str.contains(r'\|', na=True)
	multimorph['mask6'] = multimorph['class6'].str.contains(r'\|', na=True)

	# determine true stem
	multimorph['stem'] = np.select([(multimorph.mask6==False), (multimorph.mask5==False), (multimorph.mask4==False), (multimorph.mask3==False), (multimorph.mask2==False)], [multimorph.morph6, multimorph.morph5, multimorph.morph4, multimorph.morph3, multimorph.morph2], default=multimorph.morph1)

	#print(multimorph[['form_index','stem', 'morph1', 'mask1', 'morph2', 'mask2', 'morph3', 'mask3', 'morph4', 'mask4', 'morph5', 'mask5', 'morph6', 'mask6']])

	# merge true stems with words having correct stems
	my_db = pd.merge(my_db, multimorph[['form_index','stem', 'morph1', 'morph2', 'morph3', 'morph4', 'morph5', 'morph6']], on='form_index', how='left')
	my_db.stem.fillna(my_db.lemma, inplace=True)

	# get stem frequency
	freqs = my_db[['word', 'cob']]
	freqs = freqs.rename(columns = {'word':'stem', 'cob':'stem_cob'})
	freqs.sort_values(by=['stem', 'stem_cob'], inplace=True, ascending=True)
	freqs = freqs.drop_duplicates(subset=['stem'], keep='first')
	my_db = pd.merge(my_db, freqs, on='stem', how='left')

	return(my_db)

#-------------------------------------------------
# Main caller
def main():

	# parse arguments
	lemma, form, freq, phon, outdir = parseargs()

	# read files
	lemma = inputfile(lemma)
	form = inputfile(form)
	freq = inputfile(freq)
	phon = inputfile(phon)

	# preprocess CELEX
	celex = preprocess_celex(lemma, form, freq, phon)

	# correct stems
	celex = true_stem(celex)

	# save file to disk
	outfile = 'celex_dutch_db.txt'
	celex.to_csv(outfile, sep='\t')

	print('CELEX DB written to: ' + outfile)

	return

#-------------------------------------------------
main()

#--- EOF -----------------------------------------
