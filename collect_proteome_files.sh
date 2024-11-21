#!/bin/env bash
#Wed  9 Oct 13:30:37 BST 2024

#constants and functions
fastamainpath="$UP/proteomes/fasta"

#to get filepath for fasta file given a proteomeid
function proteomeid2path {
  proteome_id=$1;
  if [ ${#proteome_id} -le 2 ]; then
    path=$(echo "$proteome_id" | sed 's/.$//' | sed 's/./&\//g')
  path="$proteome_id/"else
    path=$(echo "$proteome_id" | sed 's/..$//' | sed 's/./&\//g')
  fi
  fasta_file="${fastamainpath}/${path}proteome_${proteome_id}.fa";
  if [ -s $fasta_file ]; then
    echo $fasta_file
  else
    echo "no"
  fi
}

#main
if [ "X$1" == "X" ]; then
  >&2 echo "You need to provide as first argument a file with a list of proteome_ids"
  >&2 echo "example:"
  >&2 echo "   $0 ecoli_proteomeids.out"
  exit 22
fi
if [ "X$2" == "X" ]; then
  >&2 echo "You need to provide as second argument a directory name where fasta files will be linked to"
  >&2 echo "example:"
  >&2 echo "   $0 ecoli_proteomeids.out ecoli"
  exit 22
fi
proteomeslist=$1
outdir=$2
set -e
set -u

if [ ! -d "$outdir" ]; then
  >&2 echo "NOTICE: creating fastafiles dir $outdir"
  mkdir -p $outdir
fi
if [ ! -s "$proteomeslist" ]; then
  >&2 echo "ERROR: no such file $proteomeslist"
  exit 1
fi
workingdir="$(pwd)"
proteomeslist="$workingdir/$proteomeslist"
cd $outdir

for proteome_id in $(cat "$proteomeslist"); do
  file=$(proteomeid2path $proteome_id)
  if [ $file == "no" ]; then
    >&2 echo "WARNING: no fasta file for proteome_id $proteome_id"
  else
    ln -s $file .
    echo "${outdir}/$file" >>"${workingdir}/${outdir}.tsv" #also create .tsv file (useful if too many)
  fi
done

cd -
