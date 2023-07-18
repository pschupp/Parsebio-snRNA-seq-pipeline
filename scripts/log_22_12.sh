# goal is process latest run which was performed on a NextSeq 550 High throughput. Used Midi scale kit from parse
# important priors
# number of reads from next seq 550 - 400 million reads
# number of cells from midi scale - 6000 cells as input
# what is the run architecture - r1:74, index 7: 6, r2: 86
# library structure
# r1 is read, r2 is library
# 5k + 500 cells
# 2 n7 indeces
# run took 21.5 hours!
# first issue is that there are many undetermined reads, possibly due to the fact that basecallin wsa done with insufficient edit distance. This is actually unlikely as this should be the minority of reads...
#   - alot of these reads give perfectly good sequence, but they don't show what their barcode is. therefore difficult to determine reason for exclusion, some of these certainly have all N reads or are shorter, but most reads that are undetermined seem perfectly fine.
#   - it's actually estimated about 57 million reads out of an expected 400 million, so while not  great, it's an acceptable loss.
#   - might be good try to see how basecalling was done, if there is a log...
#       - it appears its doing trimming as well, with no input from me



# 1) get data from run into folder - done
cd ~/@patrick/parsebio/dec_run/221129_NS500257_0186_AH7T7KBGXN
# 2) download and install parse pipeline 
# 3) do base calling
/opt/illumina/bcl2fastq/bin/bcl2fastq \
                --input-dir ~/@patrick/parsebio/dec_run/221129_NS500257_0186_AH7T7KBGXN/Data/Intensities/BaseCalls/  \
                --output-dir ~/@patrick/parsebio/dec_run/00_basecalled_fastq  \
                --runfolder-dir . \
                --no-lane-splitting \
                --sample-sheet SampleSheet.csv \
                --with-failed-reads \
                --ignore-missing-bcls \
                --barcode-mismatches 2 \
                --fastq-compression-level 6 \
                --loading-threads 15 \
                --processing-threads 15 \
                --writing-threads 15
# stats from report: 
#
# Clusters (Raw)	Clusters(PF)	Yield (MBases)
# 428,991,768	     400,761,426	64,122
# 99% perfect barcode
# 93.3% cluster PF
# mean quality score: 34
# 1 - CAGATC
# Clusters (Raw)	Clusters(PF)	Yield (MBases)
# 28,598,047	    28,598,047	    4,576
# 2 - ACTTGA
# Clusters (Raw)	Clusters(PF)	Yield (MBases)
# 316,042,386	    316,042,386	    50,567
# Undetermined
# Clusters (Raw)	Clusters(PF)	Yield (MBases)
# 84,351,335	    56,120,993	    8,979
# most common unkown barcodes are GGGGGG (33.9M), NNNNNN (21M), and AATAAA (6.3M)
# max is 2 mismatches in this case

# 3) run parse pipeline
# 4) set up special genome dir - DONE
split-pipe \
    --mode mkref \
    --genome_name hg38 \
    --fasta /home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa \
    --genes /home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf \
    --output_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome
# 5) run parse pipeline - 
split-pipe \
    --mode all \
    --kit WT_mini \
    --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
    --fq1 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/2_S2_R1_001.fastq.gz \
    --fq2 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/2_S2_R2_001.fastq.gz \
    --output_dir /home/patrick/@patrick/parsebio/dec_run/01_parse_pipe_proc/ \
    --sample 500run A1  \
    --nthreads 15

split-pipe \
    --mode all \
    --kit WT_mini \
    --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
    --fq1 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/1/1_S1_R1_001.fastq.gz \
    --fq2 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/1/1_S1_R2_001.fastq.gz \
    --output_dir /home/patrick/@patrick/parsebio/dec_run/01_parse_pipe_proc/ \
    --nthreads 15

split-pipe \
    --mode all \
    --kit WT_mini \
    --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
    --fq1 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/2/2_S2_R1_001.fastq.gz \
    --fq2 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/2/2_S2_R2_001.fastq.gz \
    --output_dir /home/patrick/@patrick/parsebio/dec_run/02_parse_pipe_proc/ \
    --nthreads 15

split-pipe \
    --mode all \
    --kit WT_mini \
    --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
    --fq1 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/2/2_S2_R1_001.fastq.gz \
    --fq2 /home/patrick/@patrick/parsebio/dec_run/00_basecalled_fastq/2/2_S2_R2_001.fastq.gz \
    --output_dir /home/patrick/@patrick/parsebio/dec_run/03_parse_pipe_proc/ \
    --sample s9 A1-A3  \
    --sample s44 A10-A12  \
    --sample s76 A4-A6  \
    --sample s118 A7-A9  \
    --nthreads 15
