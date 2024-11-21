#!/bin/env bash
#Sun 17 Nov 12:52:18 GMT 2024

if [ "X$1" == "X" ]; then
  >&2 echo "You need to provide as first argument a file containing a list of tsv filepaths"
  >&2 echo "examples:"
  >&2 echo "   $0 ecoli.tsv"
  exit 22
fi
fasta_list=$1

set -e
set -u

if [ ! -s "$fasta_list" ]; then
  >&2 echo "ERROR: empty or missing $fasta_list"
  exit 1
fi

for filepath in $(cat $fasta_list); do
  #note: requires ffdb.py
  indexer.py -e '>' -i '^(.+\|.*)$' -f <(cat $filepath; echo '>') > ${filepath}.idx
done
