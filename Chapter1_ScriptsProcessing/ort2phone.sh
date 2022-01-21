#!/usr/bin/env bash

t_path=/vol/tensusers2/hmueller/ORT2PHO/
KALDIbin=/vol/tensusers2/eyilmaz/local/bin

OOVin=$1
Lexout=$2

export PATH=$KALDIbin:$PATH

N=1

"${KALDIbin}/phonetisaurus-apply" --model "${t_path}model.fst" --word_list $OOVin -n $N > $Lexout
