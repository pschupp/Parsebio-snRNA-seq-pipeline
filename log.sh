# steps to get the data processed
# 1) get data from run into folder - DONE
# 2) download and install parse pipeline - DONE
# 3) determine input of parse pipeline - FASTQ files... do not see them, need to base call - DONE
cd ~/@patrick/parsebio/oct_run/221014_M02564_0812_000000000-DH5B5
/opt/illumina/bcl2fastq/bin/bcl2fastq -i ~/@patrick/parsebio/oct_run/221014_M02564_0812_000000000-DH5B5/Data/Intensities/BaseCalls/ --runfolder=. -p 15 --output-dir ~/@patrick/parsebio/oct_run/00_basecalled_fastq --no-lane-splitting --sample-sheet=SampleSheet.csv --with-failed-reads --ignore-missing-bcls


cd ~/@patrick/parsebio/nov_run/
/opt/illumina/bcl2fastq/bin/bcl2fastq -i ~/@patrick/parsebio/nov_run/221106_M00179_0379_000000000-DH6KP/Data/Intensities/BaseCalls/ --runfolder=. -p 15 --output-dir ~/@patrick/parsebio/nov_run/00_basecalled_fastq --no-lane-splitting --sample-sheet=SampleSheet.csv --with-failed-reads --ignore-missing-bcls
# 4) set up special genome dir - DONE
split-pipe \
    --mode mkref \
    --genome_name hg38 \
    --fasta /home/shared/hg_align_db/GRCh38_gencode_primary/GRCh38.primary_assembly.genome.fa \
    --genes /home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf \
    --output_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome
# 5) run parse pipeline - DONE
split-pipe \
    --mode all \
    --kit WT_mini \
    --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
    --fq1 /home/patrick/@patrick/parsebio/nov_run/00_basecalled_fastq/Parse01_S1_R1_001.fastq.gz \
    --fq2 /home/patrick/@patrick/parsebio/nov_run/00_basecalled_fastq/Parse01_S1_R2_001.fastq.gz \
    --output_dir /home/patrick/@patrick/parsebio/nov_run/01_parse_pipe_proc/ \
    --sample Parse01 A1 
