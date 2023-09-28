# goal of this script is to take `.fastq.gz` files as input. Analyze their barcode which is in the read name, and correct the barcode against the whitelist and an edit distance of my chosing. The edit distance will apply to each of the three barcodes.

import pandas
import gzip
import os
import gzip
import Bio.SeqIO.QualityIO

# Return the Hamming distance between string1 and string2.
def hamming_distance(string1, string2): 
    # Start with a distance of zero, and count up
    distance = 0
    # Loop over the indices of the string
    L = len(string1)
    for i in range(L):
        # Add 1 to the distance if these two characters are not equal
        if string1[i] != string2[i]:
            distance += 1
    # Return the final count of differences
    return distance

def find_best_match(sequence, barcode_list, min_distance):
    '''Scans the barcodes in barcode list and returns
    the barcode with the lowest hamming distance, and that
    distance. Use min_distance if you have already scanned
    a different list of barcodes.'''
    # if a perfect match exists, just output as is
    if sequence in barcode_list:
        return({'input':sequence, 'correction':sequence, 'distance':0})
    for read in barcode_list:
        h_distance = hamming_distance(sequence, read)
        # since if there was a perfect match we won't be here, 
        # if distance is 1, we'll never do better
        if h_distance == 1:
            return({'input':sequence, 'correction':read, 'distance':1})
            break
        if h_distance <= min_distance:
            return({'input':sequence, 'correction':read, 'distance':h_distance})
            break
    else: 
        return({'input':sequence, 'correction':'NNNNNNNN', 'distance':'no match'})

whiteList = pandas.read_table(snakemake.params['barcodeList'], sep =',', skiprows =4)
barcodeOne = whiteList['sequence'][0:24]
barcodeTwo = whiteList['sequence'][25:121]
barcodeThree = whiteList['sequence'][123:218]

# identify the minimum inherent edit distance within each barcode set
def min_hamming_set(inputList):
    '''Determines the minimum Hamming distance of a list of barcode inputs'''
    inputList = list(inputList)
    output = []
    for i in range(len(inputList)-1):
        for bc in inputList[i+1:]:
            output.append(hamming_distance(inputList[i], bc))
    return(min(output))

minDistOne = min_hamming_set(barcodeOne) # minimum Hamming distance is 2
minDistTwo = min_hamming_set(barcodeTwo) # minimum Hamming distance is 4
minDistThree = min_hamming_set(barcodeThree) # minimum Hamming distance is 4

# for now read in one specific file and operate on it. in the future this will be an argument from snakemake, i.e. {input}.
matchesOne = []
matchesTwo = []
matchesThree = []
output = []
seqGz = gzip.open(snakemake.input[0], 'rt')
for title, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(seqGz):
    barcode = title.split(sep = '_')[1]
    bc1 = barcode[0:8]
    bc2 = barcode[8:16]
    bc3 = barcode[16:24]
    tempOne = find_best_match(bc1, barcodeOne, minDistOne)
    tempTwo = find_best_match(bc2, barcodeTwo, minDistTwo)
    tempThree = find_best_match(bc3, barcodeThree, minDistThree)
    if((tempOne['correction'] != 'NNNNNNNN') & (tempTwo['correction'] != 'NNNNNNNN') & (tempThree['correction'] != 'NNNNNNNN')):
        titleOut = title.split(sep = '_')[0] + '_' + tempOne['correction'] + tempTwo['correction'] + tempThree['correction'] + '_' + title.split(sep = '_')[2]
        output.append("@%s\n%s\n+\n%s\n" % (titleOut, sequence, qual))
    matchesOne.append(tempOne)
    matchesTwo.append(tempTwo)
    matchesThree.append(tempThree)

matchesOne = pandas.DataFrame(matchesOne)
matchesTwo = pandas.DataFrame(matchesTwo)
matchesThree = pandas.DataFrame(matchesThree)

# need to correct existing barcodes and write out log files
with gzip.open(snakemake.output[0], 'wt') as f:
    for line in output:
        f.write(line)

with open(snakemake.log[0], 'w') as f:
    f.write('Reads analyzed: ' + str(matchesOne.shape[0]))
    f.write('\nReads passing barcode filter: ' + str(len(output))) 
    f.write('\nBarcode one distance distribution (min dist ' + str(minDistOne) + '):\n')
    f.write('dist count\n' + matchesOne['distance'].value_counts().to_string())
    f.write('\nBarcode two distance distribution (min dist ' + str(minDistTwo) + '):\n')
    f.write('dist count\n' + matchesTwo['distance'].value_counts().to_string())
    f.write('\nBarcode three distance distribution (min dist ' + str(minDistThree) + '):\n')
    f.write('dist count\n' + matchesThree['distance'].value_counts().to_string())
    f.write('\nBarcode one barcode distribution (min dist ' + str(minDistOne) + '):\n')
    f.write('dist count\n' + matchesOne['correction'].value_counts().to_string())
    f.write('\nBarcode two barcode distribution (min dist ' + str(minDistTwo) + '):\n')
    f.write('dist count\n' + matchesTwo['correction'].value_counts().to_string())
    f.write('\nBarcode three barcode distribution (min dist ' + str(minDistThree) + '):\n')
    f.write('dist count\n' + matchesThree['correction'].value_counts().to_string())
