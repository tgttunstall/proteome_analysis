#!/usr/bin/env python3


import venn
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
#atb_ids
#pedro_ids
#insana_ids

with open('/home/pub/Work/data_arise_proteome/spneumo_dataset/IDs.txt', 'r') as file:
    atb_ids = file.readlines()
    print(len(atb_ids))
    atb_ids_S = set(atb_ids)
    print(len(atb_ids_S))


with open('/home/pub/Work/data_arise_proteome/spneumo_dataset/pedroIDs.txt', 'r') as file:
    pedro_ids = file.readlines()
    #print(f"Pedro Ids: {len(pedro_ids)}")
    print(len(pedro_ids))
    pedro_ids_S  = set(pedro_ids)
    print(len(pedro_ids_S))


with open('/home/pub/Work/data_arise_proteome/spneumo_dataset/insanaIDs.txt', 'r') as file:
    insana_ids = file.readlines()
    print(len(insana_ids))
    insana_ids_S = set(insana_ids)
    print(len(insana_ids_S))

#venn.venn3([atb_ids_S, pedro_ids_S, insana_ids_S], set_labels=('ATBs', 'P','G'))
#venn3([atb_ids_S, pedro_ids_S, insana_ids_S])
#plt.show()

print("printing difference")
difference1 = pedro_ids_S-insana_ids_S
print("Only in Pedro", len(difference1))
print(difference1)

difference2p = atb_ids_S - pedro_ids_S
print("Only in ATB, not in Pedro:", len(difference2p))
#print(difference2p)

difference3i= atb_ids_S-insana_ids_S
print("Only in ATB, not in Insana:", len(difference3i))
#print(difference3i)

difference4 = insana_ids_S - atb_ids_S
print("Only in Insana, not ATB:", len(difference4))
#print(difference4)

common_elements = atb_ids_S.intersection(insana_ids_S)
print("Common b/w ATB and Insana:", len(common_elements))

common_elements2 = atb_ids_S.intersection(pedro_ids_S)
print("Common b/w ATB and Pedro:", len(common_elements2))
