#!/bin/bash

while getopts a:t:q:o: flag
do
    case "${flag}" in
        a) alignment=${OPTARG};;
        t) template_dir=${OPTARG};;
        q) query_dir=${OPTARG};;
        o) output_file=${OPTARG};;
    esac
done

echo "Template,Query,Aligned length,RMSD,SeqID,TMScore" >> ${output_file}

for t in ${template_dir}/*
do
    echo $t
    for q in ${query_dir}/*
    do
        echo $q
        if [ $alignment = 'tm' ]
        then
            TMalign $t $q >> temp.txt
        elif [ $alignment = 'mm' ]
        then
            ./MMalign $t $q >> temp.txt
        fi
        v1=${t%.pdb}
        v1=${v1##*/}
        v2=${q%.pdb}
        v2=${v2##*/}
        v1+=",${v2},"
        v1+=$(awk '/^Aligned length/ {print $3,$5,$7}' temp.txt)
        v1+=","
        if [ $alignment = 'tm' ]
        then
            v1+=$(awk '/by length of Chain_1/ {print $2}' temp.txt)
        elif [ $alignment = 'mm' ]
        then
            v1+=$(awk '/by length of Structure_1/ {print $2}' temp.txt)
        fi
        echo $v1 | tr -d ' ' >> ${output_file}
        rm temp.txt
    done
done