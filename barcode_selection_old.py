#Define the path where fastq file is stored.
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator 
from time import localtime, strftime
from collections import defaultdict
from collections import Counter

path=os.path.join(os.path.expanduser("~"), "Desktop", "SVA1_S1_R2_001.fastq")
stats = counter()
#create a new file handle to write the file
path_new=os.path.join(os.path.expanduser("~"), "Desktop", "trimmedfastq.fastq")
path_ust_barcode=os.path.join(os.path.expanduser("~"), "Desktop", "UST40MBarcodes.txt")
write_file=open(path_new, "w")
#open the text file containg the correct pool of barcodes
now = strftime("%H:%M:%S", localtime())
print ("#%s\tLoading barcodes")
barcodes = set()
barcode_tree = defaultdict(set)
with open(path_ust_barcode, "r") as fh:
   for bc in fh:
       barcodes.add(bc.rstrip())
       barcode_tree[bc[:11]].add(bc[11:])
now = strftime("%H:%M:%S", localtime())
print ("#%s\tBarcodes loaded\n#%s\tMatching reads" % (now, now)))
handle=open(path, "r")  
output = []
#open the fastq file
for title, sequence, qual in FastqGeneralIterator(handle):
    stats["reads"] += 1
    if sequence in barcodes:
        stats["perfect_match"] += 1
        output.append("@%s\n%s+\n%s\n" % (title, sequence, qual))
        continue
    if sequence[:11] in barcode_tree:
        corr_barcode = sequence[:11] + find_best_match(sequence[11:], barcode_tree[sequence[:11]])
        stats["First 11 match"] += 1
    else:
        stats["Full search"] += 1
        corr_barcode = find_best_match(sequence, barcodes)
    output.append("@%s\n%s+\n%s\n" % (title, corr_barcode, qual))
    if stats["reads"] % 100 == 0:
        now = strftime("%H:%M:%S", localtime())
        print ("#%s\t%i reads processed, %i perfect matches, %i first 11 matches, %i full searches" % 
               (now, stats["reads"], stats["perfect_match"], stats["First 11 match"], stats["Full search"]))
now = strftime("%H:%M:%S", localtime())
print ("#%s\tAll Reads Processed\n#%s\tOutputting Results" % (now, now)))
with open(path_new, "a") as wf:
    for line in output:
        wf.write(line)
handle.close()
now = strftime("%H:%M:%S", localtime())
print("#%s\tDone" % now)
