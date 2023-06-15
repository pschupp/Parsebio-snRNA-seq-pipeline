# basecalling has been done with only 10%
# # 1) get data from run into folder - done
cd /home/patrick/@patrick/parsebio/jan_run
# 3) run parse pipeline
# 5) run parse pipeline - 
done
for i in {1..5}
do
done

for i in {1..8}
do
    echo "${i}_S${i}_L001_R1_001.fastq.gz"
    echo "${i}_S${i}_L001_R2_001.fastq.gz"
    split-pipe \
        --mode all \
        --kit WT \
        --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
        --fq1 "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/Data/Intensities/BaseCalls/${i}_S${i}_L001_R1_001.fastq.gz" \
        --fq2 "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/Data/Intensities/BaseCalls/${i}_S${i}_L001_R2_001.fastq.gz" \
        --output_dir "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/lib_${i}" \
        --nthreads 15
done

for i in {1..8}
do
    echo "${i}_S${i}_L001_R1_001.fastq.gz"
    echo "${i}_S${i}_L001_R2_001.fastq.gz"
    split-pipe \
        --mode all \
        --kit WT \
        --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
        --fq1 "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/Data/Intensities/BaseCalls/${i}_S${i}_L001_R1_001.fastq.gz" \
        --fq2 "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/Data/Intensities/BaseCalls/${i}_S${i}_L001_R2_001.fastq.gz" \
        --tscp_use 20 \
        --output_dir "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/lib_20_cutoff_${i}" \
        --nthreads 15
done


split-pipe \
        --mode comb \
        --kit WT \
        --genome_dir /home/shared/hg_align_db/GRCh38_gencode_primary/parse_split_pipe_genome \
        --sublibraries ./lib_1/ ./lib_5/ ./lib_3/ ./lib_6/ ./lib_7/ ./lib_8/ \
        --output_dir "~/@patrick/parsebio/jan_run/Parse_100K_230210_M00179_0395_000000000-G6B9V/combined" \
        --nthreads 15
