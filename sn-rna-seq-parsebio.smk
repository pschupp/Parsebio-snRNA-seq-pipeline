# start with expectation that input is folder from sequencing machine in the form:
# [YY][MM][DD]_M00179_0434_000000000-DL3KG
# also includes the following usefull files:
# RTA Logs folder - Contains log files that describe each step performed by RTA (real-time analysis) for each read
# InterOp folder - Contains binary files used by Sequencing Analysis Viewer (SAV) to summarize various primary analysis metric such as cluster density, intensities, quality scores, and overall run quality
# Logs folder - Contains log files that describe eavery step performed by the intrument for each cycle
# RunInfo.xml - Contains high-level run information, such as the number of reads and cycles in the sequencing run.
# runParameters.xml - Contains a summary of run parameters and information about run components sucah as the RFID of the flow cell and reagents associated with the run.
# SampleSheet.csv - Provides parameters for teh run and subsequent analysis.
# Images and Thumbnail_Images folders - Contains focus images and thumbnail images. 

threads = workflow.cores
import os
import numpy
import pandas
import re
import math
cwd = os.getcwd()

folder = ''.join([fol for fol in os.listdir() if re.search('^[0-9]{6}_.*', fol)])

# find when '[Data]' begins
i=0
with open(folder + '/SampleSheet.csv', 'r') as f:
    for line in f:
        if 'Data' in line:
            inLine = i + 1
            break
        i=i+1

# read in the sample table under '[Data]'
sample=pandas.read_table(folder + '/SampleSheet.csv', header=0, skiprows=inLine, sep=',', dtype = 'str')

# open file until '[Data]' and keep this header
with open(folder + '/SampleSheet.csv', 'r') as f:
    i=0
    lines = []
    for line in f:
        i = i + 1
        lines.append(line)
        if i == inLine + 0:
            break

# write only header to file
with open(folder + '/SampleSheet.csv', 'w') as f:
    for line in lines:
        _ = f.write(line)

# determine which names need to be replaced
replaceIndeces = [str(x) == 'nan' for x in sample['Sample_Name']]
# create names for unnamed librararies
sample.loc[replaceIndeces,'Sample_Name'] = ['unamed-sublibrary-'+str(i) for i in numpy.where(replaceIndeces)[0]]
# write out new sample table with names
sample.to_csv(path_or_buf = folder + '/SampleSheet.csv', sep = ',', header = True, index = False, mode = 'a') 
# create sample names and file names (with read 1/2 label)
sampleN = sample['Sample_Name'] + '_S' + sample['Sample_ID']
rn = [['_R1']*len(sample['Sample_Name']), ['_R2']*len(sample['Sample_Name'])]
rn = [x for y in rn for x in y]
fileN = sample['Sample_Name'] + '_S' + sample['Sample_ID']
fileN = (pandas.concat([fileN, fileN]) + rn).sort_values().tolist()

rule all:
    input: "09_multiqc/multiqc_report.html"
            
rule bcl_2_fastq:
    output: "01-basecalled/basecalling.log.txt"
    params: runfolder = folder
    shell:
        """
        /opt/illumina/bcl2fastq/bin/bcl2fastq \\
        --runfolder={params.runfolder} \\
        --output-dir=01-basecalled \\
        --sample-sheet={params.runfolder}/SampleSheet.csv \\
        --minimum-trimmed-read-length=30 \\
        --with-failed-reads \\
        --barcode-mismatches=1 \\
        --no-lane-splitting \\
        --loading-threads=5 \\
        --processing-threads=15 \\
        --writing-threads=5 \\
        2> 01-basecalled/basecalling.log.txt
        """
# --minimum-trimmed-read-length arg (=35)         minimum read length after adapter trimming
# --with-failed-reads                             include non-PF clusters
# --no-lane-splitting                             do not split fastq files by lane.
# --barcode-mismatches arg (=1)                   number of allowed mismatches per index. Multiple, comma delimited, entries allowed. Each entry is applied to the corresponding index; last entry applies to all remaining indices. Accepted values: 0, 1, 2.
# --loading-threads Number of threads to load BCL data
# --processing-threads Number of threads to process demultiplexing data.
# --writing-threads Number of threads to write FASTQ data. This number must be lower than number of samples.

#"02-fastqc/{file_names}.std_out.txt"
rule fastqc:
    input: "01-basecalled/basecalling.log.txt"
    output: "02-fastqc/{fileFull}_001_fastqc.html"
    params: "01-basecalled/{fileFull}_001.fastq.gz"
    log: 
            stdOut = "02-fastqc/{fileFull}.std_out.txt",
            stdErr = "02-fastqc/{fileFull}.std_err.txt"
        #         stdOut = expand("02-fastqc/{fileFull}.std_out.txt", fileFull=fileN),
        #         stdErr = expand("02-fastqc/{fileFull}.std_err.txt", fileFull=fileN)
    threads: 3
    shell:
        """
        fastqc {params} \\
        --outdir 02-fastqc/ \\
        --noextract \\
        --threads {threads} \\
        > {log.stdOut} \\
        2> {log.stdErr} 
        """

# Use read 2 to deconvolute barcode and extract umi. 
# TODO: is there any hamming code correction? 
rule umi_tools:
    input: "01-basecalled/basecalling.log.txt"
    output: "03-umi-extract/{file_names}_R1_001.fastq.gz"
    log: "03-umi-extract/{file_names}.log.txt"
    params:
        bcRegex=lambda wildcards: '"(?P<umi_1>.{10})(?P<cell_3>.{8})(?P<discard_1>GTGGCCGATGTTTCGCATCGGCGTACGACT){s<=30}(?P<cell_2>.{8})(?P<discard_2>ATCCACGTGCTTGAGACTGTGG){s<=22}(?P<cell_1>.{8}).*"',
        files = "01-basecalled/{file_names}_R1_001.fastq.gz",
        whitelist = "."
    threads: 1
    shell:
        """
        umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern={params.bcRegex} \\
        --read2-in={params.files} \\
        --read2-out={output} \\
        --stdin=01-basecalled/{wildcards.file_names}_R2_001.fastq.gz \\
        --stdout=03-umi-extract/{wildcards.file_names}_R2_001.fastq.gz \\
         --log={log}
        """
#         --whitelist={params.whitelist} 

rule multiqc:
    input: 
        inDir=".",
        fastqcReport = expand("02-fastqc/{fileFull}.std_err.txt", fileFull=fileN),
        processedReads = expand("03-umi-extract/{file_names}_R1_001.fastq.gz", file_names = sampleN)
    output:
        outFile="09_multiqc/multiqc_report.html"
    shell:
        """
        multiqc {input.inDir} \\
        -o 09_multiqc
        """
# 
# rule trimmomatic:
#     input:"01_basecalled/{file_names}.umi.fastq.gz"
#     output:"02_trimmed/{file_names}.trimmed.fastq"
#     params:
#         headcrop=4,
#         clipfile="/usr/share/trimmomatic/TruSeq3-SE.polyA.fa"
#     threads: 1
#     shell:
#         """
#         TrimmomaticSE {input} {output} \\
#         -trimlog {output}.trim_full_log.txt \\
#         -threads {threads} \\
#         ILLUMINACLIP:{params.clipfile}:2:30:7 \\
#         LEADING:10 \\
#         TRAILING:10 \\
#         SLIDINGWINDOW:4:15 \\
#         MINLEN:30 \\
#         HEADCROP:{params.headcrop} \\
#         > 02_trimmed/{wildcards.file_names}.std_out.txt \\
#         2> 02_trimmed/{wildcards.file_names}.std_err.txt
#         """
# 
# rule star_align:
#     input: 
#         dummy="03_posttrim_fastqc/{file_names}.std_out.txt",
#         files="03_trimmed/{file_names}.trimmed.fastq"
#     output: "04_aligned/{file_names}.Aligned.out.bam"
#     threads: 15
#     params:
#         genomeDir=config['genomeDir']
#     shell:
#         """
#         /opt/STAR_old/STAR-2.6.1e/bin/Linux_x86_64/STAR --genomeDir config['genomeDir'] \\
#         --readFilesIn {input.files} --readFilesCommand cat \\
#         --outFileNamePrefix 04_aligned/{wildcards.file_names}. \\
#         --outSAMtype BAM Unsorted \\
#         --outFilterMultimapNmax 4 \\
#         --outFilterType BySJout \\
#         --outSAMmultNmax 1  \\
#         --alignSJoverhangMin 5 \\
#         --alignSJDBoverhangMin 1 \\
#         --outFilterMismatchNmax 999 \\
#         --outFilterMismatchNoverReadLmax 0.3 \\
#         --outFilterScoreMinOverLread 0.3 \\
#         --outFilterMatchNminOverLread 0.3 \\
#         --alignIntronMin 20 \\
#         --alignIntronMax 1000000 \\
#         --genomeLoad LoadAndKeep --limitBAMsortRAM 322122382273  \\
#         --runThreadN {threads} \\
#         > 04_aligned/{wildcards.file_names}.std_out.txt
#         2> 04_aligned/{wildcards.file_names}.std_err.txt
#         """
# # already done by trimmomatic
# #        --clip3pAdapterSeq AAAAAAAA --clip3pAdapterMMp 0.1 \\
# 
# rule create_field_OL:
#     input: "04_aligned/{file_names}.Aligned.out.bam"
#     output: "04_aligned/{file_names}.fixed_tag.bam"
#     threads: 1
#     shell:
#         """
#         set +o pipefail;
#         samtools view -h -@ {threads} {input} \\
#         | sed -E 's/([0-9]{{3,}})(:)([0-9]{{3,}})(:)([0-9]{{3,}})(_)([A-Z]{{6}})(.*)$/\\1:\\3:\\5:\\7\\8\tOL:Z:\\7/g' \\
#         | samtools view -O BAM -b -@ {threads} -o {output} - 
#         """
# 
# rule sort_query:
#     input: "04_aligned/{file_names}.fixed_tag.bam"
#     output: "05_sorted/{file_names}.bam"
#     threads: 5
#     shell:
#         """
#         PicardCommandLine SortSam \\
#         I={input} \\
#         O={output} \\
#         SORT_ORDER=queryname \\
#         > 05_sorted/{wildcards.file_names}.sort.std_out.txt \\
#         2> 05_sorted/{wildcards.file_names}.sort.std_err.txt
#         """
# 
# rule deduplicate:
#     input: "05_sorted/{file_names}.bam"
#     output: "06_deduplicated/{file_names}.deduplicated.bam"
#     threads: 5
#     shell:
#         """
#         PicardCommandLine MarkDuplicates \\
#         I={input} \\
#         O={output} \\
#         TAGGING_POLICY=All \\
#         ASSUME_SORT_ORDER=queryname \\
#         BARCODE_TAG=OL \\
#         MOLECULAR_IDENTIFIER_TAG=OL \\
#         REMOVE_DUPLICATES=TRUE \\
#         READ_NAME_REGEX=null \\
#         SORTING_COLLECTION_SIZE_RATIO=.8 \\
#         MAX_FILE_HANDLES=1000 \\
#         METRICS_FILE={output}.metrics \\
#         > 06_deduplicated/{wildcards.file_names}.markdup.std_out.txt \\
#         2> 06_deduplicated/{wildcards.file_names}.markdup.std_err.txt
#         """
# 
# rule feature_counts:
#     input: "06_deduplicated/{file_names}.deduplicated.bam"
#     output: "07_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
#     threads: 15
#     params:
#         gtfFile=config['gtfFile']
#     shell:
#         """
#         featureCounts \\
#         {input} \\
#         -a {params.gtfFile} \\
#         -o 08_feature_count/{wildcards.file_names}.count_matrix.tsv \\
#         -R BAM \\
#         -g gene_id \\
#         -M -O \\
#         -t gene \\
#         --fraction \\
#         -T {threads} \\
#         > 08_feature_count/{wildcards.file_names}.std_out.txt \\
#         2> 08_feature_count/{wildcards.file_names}.std_err.txt
#         """
# # options:
# # -a gtf file input
# # -o output file count matrix
# # -R input file format
# # -g ninth column of gtf file. the attribute type by which to group features. look to gtf file for correct naming. could also be transcript_id or exon_id
# # -t third column of gtf file. indicates what type of features should be considered. When using -g gene_id should use gene. If using transcript_id or exon_id, should use processed transcript.
# # -M multi-mapping reads will counted. Uses NH tag for multimappers.
# # -O assign reads to all overlaping features
# # --fraction Assign fractional counts to features
# 
# rule bam_sort:
#     input: "07_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
#     output: "08_sort_index/{file_names}.sort.bam"
#     shell:
#         """
#         samtools sort \\
#         {input} \\
#         -o {output}
#         > 08_sort_index/{wildcards.file_names}.sort.std_out.txt \\
#         2> 08_sort_index/{wildcards.file_names}.sort.std_err.txt
#         """
# 
# rule bam_index:
#     input: "08_sort_index/{file_names}.sort.bam"
#     output: "08_sort_index/{file_names}.sort.bam.bai"
#     shell:
#         """
#         samtools index \\
#         {input} \\
#         > 08_sort_index/{wildcards.file_names}.index.std_out.txt \\
#         2> 08_sort_index/{wildcards.file_names}.index.std_err.txt
#         """
# 
# rule genebody_cov:
#     input: 
#         real="08_sort_index/{file_names}.sort.bam",
#         dummy="08_sort_index/{file_names}.sort.bam.bai"
#     output: "08_sort_index/{file_names}.geneBodyCoverage.curves.pdf"
#     params:
#         bedFile=config['bedFile']
#     shell:
#         """
#         geneBody_coverage.py \\
#         -i {input.real} \\
#         -r {params.bedFile} \\
#         -o 11_sort_index/{wildcards.file_names} \\
#         -f pdf
#         """
# 
# rule multiqc:
#     input: 
#         inDir=".",
# #        bodycov=expand("09_sort_index/{file_names}.sort.bam.bai", file_names=sample),
#         extra=expand("08_sort_index/{file_names}.geneBodyCoverage.curves.pdf", file_names=sampleN)
# #         dummy=expand("08_feature_count/{file_names}.deduplicated.bam.featureCounts.bam", file_names=sample)
#     output:
#         outFile="09_multiqc/multiqc_report.html"
#     shell:
#         """
#         multiqc {input.inDir} \\
#         -o 09_multiqc; \\
#         multiqc 02_pretrim_fastqc -o 02_pretrim_fastqc
#         """
