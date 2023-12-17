# Parsebio pipeline

## Quick-start

Copy the pipeline file `sn-rna-seq-parsebio.smk`, the `scripts` folder and the `barcodes` folder into the base directory of the sequencing run.

The script also expects 
    1. a folder starting with six digits (the date). Inside this folder should be the Sample Sheet (`SampleSheet.csv`) with relevant names and barcodes. For more information on what the Sample Sheet should look like, see XXX. 
    2. If basecalled already, add all `.fastq.gz` files into a folder called `01_basecalled`.

## FAQ

### Why do I have to copy the files?
Stuff about version control etc...


## Kits and associated barcodes

| kit     | chem | nwells | bc1     | bc2 | bc3 | ktype   |
| ---     | ---  | ---    | ---     | --- | --- | ---     |
| custom  | NA   | 1      | NA      | NA  | NA  | special |
| WT_mini | v1   | 12     | n24_v4  | v1  | v1  | normal  |
| WT_mini | v2   | 12     | n24_v4  | v1  | v1  | normal  |
| WT      | v1   | 48     | v2      | v1  | v1  | normal  |
| WT      | v2   | 48     | n96_v4  | v1  | v1  | normal  |
| WT_mega | v1   | 96     | n192_v4 | v1  | v1  | normal  |
| WT_mega | v2   | 96     | n192_v4 | v1  | v1  | normal  |
