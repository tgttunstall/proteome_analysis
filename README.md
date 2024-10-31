# Task is to add the proteom ID to the file containsing protein clusters so we know for which proteome the protein belongs to
## Overview of the files:
### proteome fasta
For the pig data, there are about 14 files in the proteome dir, each containing about 50-60k proteins beginning with ">":

PR1.fa: 50k
>ABC1|7
ABCDEFG
>ABC2|10
EFGHIJKLMN

PR2.fa: 60k
>PQR1|6
XYZABC
>PQR2|12
RSTUVGHJKKVV

### protein cluster.tsv
This is the results file containing the protein cluster analysis. In the pig dataset, there are 680633 lines in the format:

ABC1|7	ABC1|7
ABC1|7	PQR1|6
ABC2|10	PQR1|6
ABC2|10	PQR2|12