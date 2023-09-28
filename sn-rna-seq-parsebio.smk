# start with expectation that input is folder from sequencing machine in the form:
# [YY][MM][DD]_M00179_0434_000000000-DL3KG
# {{{ expalanation of useful files in the run folder
# RTA Logs folder - Contains log files that describe each step performed by RTA (real-time analysis) for each read
# InterOp folder - Contains binary files used by Sequencing Analysis Viewer (SAV) to summarize various primary analysis metric such as cluster density, intensities, quality scores, and overall run quality
# Logs folder - Contains log files that describe eavery step performed by the intrument for each cycle
# RunInfo.xml - Contains high-level run information, such as the number of reads and cycles in the sequencing run.
# runParameters.xml - Contains a summary of run parameters and information about run components sucah as the RFID of the flow cell and reagents associated with the run.
# SampleSheet.csv - Provides parameters for teh run and subsequent analysis.
# Images and Thumbnail_Images folders - Contains focus images and thumbnail images. 
# }}}

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
# {{{ documentation for options in bcl2fastq
# --minimum-trimmed-read-length arg (=35)         minimum read length after adapter trimming
# --with-failed-reads                             include non-PF clusters
# --no-lane-splitting                             do not split fastq files by lane.
# --barcode-mismatches arg (=1)                   number of allowed mismatches per index. Multiple, comma delimited, entries allowed. Each entry is applied to the corresponding index; last entry applies to all remaining indices. Accepted values: 0, 1, 2.
# --loading-threads Number of threads to load BCL data
# --processing-threads Number of threads to process demultiplexing data.
# --writing-threads Number of threads to write FASTQ data. This number must be lower than number of samples.
# }}}

rule fastqc:
    input: "01-basecalled/basecalling.log.txt"
    output: "02-fastqc/{fileFull}_001_fastqc.html"
    params: "01-basecalled/{fileFull}_001.fastq.gz"
    log: 
            stdOut = "02-fastqc/{fileFull}.std_out.txt",
            stdErr = "02-fastqc/{fileFull}.std_err.txt"
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

rule umi_tools_extract:
    input: "01-basecalled/basecalling.log.txt"
    output: "02-barcode-umi-extract/{file_names}_R1_001.fastq.gz"
    log: "02-barcode-umi-extract/{file_names}.log.txt"
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
        --stdout=02-barcode-umi-extract/{wildcards.file_names}_R2_001.fastq.gz \\
        --log={log}
        """

rule custom_barcode_correction_script:
    input: "02-barcode-umi-extract/{file_names}_R1_001.fastq.gz"
    output: "03-barcode-filter/{file_names}_R1_001.fastq.gz"
    log: "03-barcode-filter/{file_names}.log.txt"
    params:
        barcodeList = 'barcode_data.csv'
    threads: 1
    script:
        "scripts/barcode_selection_correction.py"

rule trimmomatic:
    input: "03-barcode-filter/{file_names}_R1_001.fastq.gz"
    output: "04-trim-reads/{file_names}_R1_001.fastq.gz"
    params:
        headcrop=4,
        clipfile="/usr/share/trimmomatic/TruSeq3-SE.polyA.fa"
    threads: 1
    shell:
        """
        TrimmomaticSE {input} {output} \\
        -trimlog {output}.trim_full_log.txt \\
        -threads {threads} \\
        ILLUMINACLIP:{params.clipfile}:2:30:7 \\
        LEADING:10 \\
        TRAILING:10 \\
        SLIDINGWINDOW:4:15 \\
        MINLEN:30 \\
        HEADCROP:{params.headcrop} \\
        > 04-trim-reads/{wildcards.file_names}.std_out.txt \\
        2> 04-trim-reads/{wildcards.file_names}.std_err.txt
        """

rule star_align:
    input: "04-trim-reads/{file_names}_R1_001.fastq.gz"
    output: "05_aligned/{file_names}.Aligned.out.bam"
    threads: 15
    params:
        genomeDir=config['genomeDir']
    shell:
        """
        /opt/STAR_old/STAR-2.6.1e/bin/Linux_x86_64/STAR \\
        --genomeDir {params.genomeDir} \\
        --readFilesIn {input} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix 05_aligned/{wildcards.file_names}. \\
        --outSAMtype BAM Unsorted \\
        --outFilterMultimapNmax 4 \\
        --outFilterType BySJout \\
        --outSAMmultNmax 1  \\
        --alignSJoverhangMin 5 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.3 \\
        --outFilterScoreMinOverLread 0.3 \\
        --outFilterMatchNminOverLread 0.3 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --genomeLoad LoadAndKeep \\
        --limitBAMsortRAM 322122382273  \\
        --runThreadN {threads} \\
        > 05_aligned/{wildcards.file_names}.std_out.txt \\
        2> 05_aligned/{wildcards.file_names}.std_err.txt
        """

rule create_field_CB_and_UB:
    input: "05_aligned/{file_names}.Aligned.out.bam"
    output: "05_aligned/{file_names}.fixed_tag.bam"
    threads: 1
    shell:
        """
        set +o pipefail;
        samtools view -h -@ {threads} {input} \\
                | sed -E 's/_([A-Z]{{24}})_([A-Z]{{10}})(.*)$/_\\1_\\2\\3\tCB:Z:\\1\\2\tUB:Z:\\2\tBC:Z:\\1/g' \\
        | samtools view -O BAM -b -@ {threads} -o {output} - 
        """
# {{{ documentation for command        
# create field CB (cell barcode) and UB (unique UMI barcode)
# }}}

rule sort_query:
   input: "05_aligned/{file_names}.fixed_tag.bam"
   output: "06_sorted/{file_names}.bam"
   threads: 5
   shell:
       """
       PicardCommandLine SortSam \\
       I={input} \\
       O={output} \\
       SORT_ORDER=coordinate\\
       > 06_sorted/{wildcards.file_names}.sort.std_out.txt \\
       2> 06_sorted/{wildcards.file_names}.sort.std_err.txt
       """

rule deduplicate:
    input: "06_sorted/{file_names}.bam"
    output: "07_deduplicated/{file_names}.deduplicated.bam"
    threads: 5
    shell:
        """
        PicardCommandLine MarkDuplicates \\
        INPUT={input} \\
        METRICS_FILE={output}.barcode-metrics.txt \\
        OUTPUT={output} \\
        ASSUME_SORT_ORDER=coordinate \\
        BARCODE_TAG=CB \\
        CLEAR_DT=FALSE \\
        DUPLEX_UMI=FALSE \\
        DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \\
        MAX_FILE_HANDLES=1000 \\
        MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
        REMOVE_DUPLICATES=TRUE \\
        SORTING_COLLECTION_SIZE_RATIO=.8 \\
        TAGGING_POLICY=All \\
        > 07_deduplicated/{wildcards.file_names}.markdup.std_out.txt \\
        2> 07_deduplicated/{wildcards.file_names}.markdup.std_err.txt
        """
# {{{ options for deduplicate tool,
# note that UmiAwareMarkDuplicatesWithMateCigar does not seem to work for me with this data.
# need to deduplicate based on sublibrary (SM/CB tag) and UMI (OL/UB tag)
# --INPUT,-I:String             One or more input SAM or BAM files to analyze. Must be coordinate sorted. This argument must be specified at least once. Required .
# --METRICS_FILE,-M:File        File to write duplication metrics to  Required.
# --OUTPUT,-O:File              The output file to write marked records to  Required. 
# --arguments_file:File read one or more arguments files and add them to the command line  This argument may be specified 0 or more times. Default value: null.
# --ASSUME_SORT_ORDER,-ASO:SortOrder  If not null, assume that the input file has this order even if the header says otherwise. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate,  unknown}  Cannot be used in conjuction with argument(s) ASSUME_SORTED (AS)  
# --ASSUME_SORTED,-AS:Boolean If true, assume that the input file is coordinate sorted even if the header says otherwise. Deprecated, used ASSUME_SORT_ORDER=coordinate instead.  Default value: false. Possible values: {true, false}  Cannot be used in conjuction with argument(s)  ASSUME_SORT_ORDER (ASO) 
# --BARCODE_TAG:String  Barcode SAM tag (ex. BC for 10X Genomics)  Default value: null. 
# --CLEAR_DT:Boolean  Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this tag.  Default true  Default value: true. Possible values: {true, false} 
# --COMMENT,-CO:String  Comment(s) to include in the output file's header.  This argument may be specified 0 or more times. Default value: null.  
# --DUPLEX_UMI:Boolean  Treat UMIs as being duplex stranded.  This option requires that the UMI consist of two equal length strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered duplicates if, in addition to standard definition, have identical normalized UMIs.  A UMI from the 'bottom' strand is normalized by swapping its content around the hyphen (eg. ATC-GTC becomes GTC-ATC).  A UMI from the 'top' strand is already normalized as it is. Both reads from a read pair considered top strand if the read 1 unclipped 5' coordinate is less than the read 2 unclipped 5' coordinate. All chimeric reads and read fragments are treated as having come from the top strand. With this option is it required that the  BARCODE_TAG hold non-normalized UMIs. Default false.  Default value: false. Possible values: {true, false}  
# --DUPLICATE_SCORING_STRATEGY,-DS:ScoringStrategy  The scoring strategy for choosing the non-duplicate among candidates.  Default value: 
#   SUM_OF_BASE_QUALITIES. Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}  
# --help,-h:Boolean display the help message  Default value: false. Possible values: {true, false}  
# --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP,-MAX_FILE_HANDLES:Integer Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system.  Default value: 8000. 
# --MAX_OPTICAL_DUPLICATE_SET_SIZE:Long  This number is the maximum size of a set of duplicate reads for which we will attempt to  determine which are optical duplicates.  Please be aware that if you raise this value too high and do encounter a very large set of duplicate reads, it will severely affect the  runtime of this tool.  To completely disable this check, set the value to -1.  Default  value: 300000.  
# --MOLECULAR_IDENTIFIER_TAG:String SAM tag to uniquely identify the molecule from which a read was derived.  Use of this option requires that the BARCODE_TAG option be set to a non null value.  Default null.  Default value: null. 
# --OPTICAL_DUPLICATE_PIXEL_DISTANCE:Integer The maximum offset between two duplicate clusters in order to consider them optical   duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and  models, users should experiment to find what works best.  Default value: 100. 
# --PROGRAM_GROUP_COMMAND_LINE,-PG_COMMAND:String Value of CL tag of PG record to be created. If not supplied the command line will be   detected automatically.  Default value: null. 
# --PROGRAM_GROUP_NAME,-PG_NAME:String  Value of PN tag of PG record to be created.  Default value: MarkDuplicates. 
# --PROGRAM_GROUP_VERSION,-PG_VERSION:String  Value of VN tag of PG record to be created. If not specified, the version will be detected  automatically.  Default value: null.  
# --PROGRAM_RECORD_ID,-PG:StringThe program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation.  This string may have a suffix appended to avoid collision  with other program record IDs.  Default value: MarkDuplicates.  
# --READ_NAME_REGEX:String  MarkDuplicates can use the tile and cluster positions to estimate the rate of optical duplication in addition to the dominant source of duplication, PCR, to provide a more accurate estimation of library size. By default (with no READ_NAME_REGEX specified),  MarkDuplicates will attempt to extract coordinates using a split on ':' (see Note below).  Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. Note that without optical duplicate counts, library size estimation will be less accurate. If the read name  does not follow a standard Illumina colon-separation convention, but does contain tile and  x,y coordinates, a regular expression can be specified to extract three variables:  tile/region, x coordinate and y coordinate from a read name. The regular expression must  contain three capture groups for the three variables, in order. It must match the entire  read name. e.g. if field names were separated by semi-colon (';') this example regex  could be specified  (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ Note that if no  READ_NAME_REGEX is specified, the read name is split on ':'. For 5 element names, the  3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values.  Default value: <optimized capture of last three ':' separated fields as numeric values>.  
# --READ_ONE_BARCODE_TAG:String Read one barcode SAM tag (ex. BX for 10X Genomics)  Default value: null.  
# --READ_TWO_BARCODE_TAG:String Read two barcode SAM tag (ex. BX for 10X Genomics)  Default value: null.  
# --REMOVE_DUPLICATES:Boolean If true do not write duplicates to the output file instead of writing them with appropriate flags set.  Default value: false. Possible values: {true, false}  
# --REMOVE_SEQUENCING_DUPLICATES:Boolean  If true remove 'optical' duplicates and other duplicates that appear to have arisen from  the sequencing process instead of the library preparation process, even if  REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and  this option is ignored.  Default value: false. Possible values: {true, false} 
# --SORTING_COLLECTION_SIZE_RATIO:Double  This number, plus the maximum RAM available to the JVM, determine the memory footprint   used by some of the sorting collections.  If you are running out of memory, try reducing  this number.  Default value: 0.25.  
# --TAG_DUPLICATE_SET_MEMBERS:Boolean  If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), indicates the size of the duplicate set. The smallest possible DS value is 2 which   occurs when two reads map to the same portion of the reference only one of which is marked  as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier  for the duplicate set to which the record belongs. This identifier is the index-in-file of  the representative read that was selected out of the duplicate set.  Default value: false.  Possible values: {true, false}  
# --TAGGING_POLICY:DuplicateTaggingPolicy Determines how duplicate types are recorded in the DT optional attribute.  Default value: DontTag. Possible values: {DontTag, OpticalOnly, All} 
# --version:Boolean display the version number for this tool  Default value: false. Possible values: {true, false}
# }}}

rule bam_index:
    input: "07_deduplicated/{file_names}.deduplicated.bam"
    output: "07_deduplicated/{file_names}.deduplicated.bam.bai"
    shell:
        """
        samtools index \\
        {input} \\
        > 07_deduplicated/{wildcards.file_names}.index.std_out.txt \\
        2> 07_deduplicated/{wildcards.file_names}.index.std_err.txt
        """
# featureCounts does not include options to account umi and barcodes
# can still use feature counts if we splitup the file into the component single-nucleus libraries
# otherwise will use htseq-counts which is single-threaded but fine for this case 

# option 1: split each bam file into a bam file for each single-nucleus library
# option 2: use htseq-counts, with barcoded, position-sorted bams
# for now try option 2
# cb tag is barcode, ub tag is umi
# will need to investigage options 1, options 2 is way too slow

rule hts_seq_counts:
    input: 
        file="07_deduplicated/{file_names}.deduplicated.bam",
        index="07_deduplicated/{file_names}.deduplicated.bam.bai"
    output: "08_feature_count/{file_names}.deduplicated.bam.featureCounts.tsv"
    threads: 1
    params:
        gtfFile=config['gtfFile']
    shell:
        """
        htseq-count-barcodes {input.file} {params.gtfFile} \\
        --stranded yes \\
        --minaqual 0 \\
        --type gene \\
        --idattr gene_id \\
        --mode union \\
        --nonunique all \\
        --secondary-alignments ignore \\
        --supplementary-alignments ignore \\
        --delimiter "\t" \\
        --counts_output {output} \\
        --cell-barcode BC \\
        --UMI UB \\
        > 08_feature_count/{wildcards.file_names}.std_out.txt \\
        2> 08_feature_count/{wildcards.file_names}.std_err.txt
        """
# {{{ options documentation
# This script takes one alignment file in SAM/BAM format and a feature file in GFF format and calculates for each feature the numberof reads mapping to it, accounting for barcodes. See http://htseq.readthedocs.io/en/master/count.html for details. positional arguments:
# - samfilename Path to the SAM/BAM file containing the barcoded, mapped reads. If '-' is selected, read from standard input
# - featuresfilename Path to the GTF file containing the features
# optional arguments:
# -h, --help show this help message and exit
# -f {sam,bam,auto}, --format {sam,bam,auto} Type of <alignment_file> data. DEPRECATED: file format is detected automatically. This option is ignored.
# -r {pos,name}, --order {pos,name} 'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing data must be sorted either by position or by read name, and the sorting order must be specified. Ignored for single-end data.
# --max-reads-in-buffer MAX_BUFFER_SIZE When <alignment_file> is paired end sorted by position, allow only so many reads to stay in memory until the mates are found (raising this number will use more memory). Has no effect for single end or paired end sorted by name
# -s {yes,no,reverse}, --stranded {yes,no,reverse} Whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: yes). 'reverse' means 'yes' with reversed strand interpretation
# -a MINAQUAL, --minaqual MINAQUAL Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads.
# -t FEATURETYPE, --type FEATURETYPE Feature type (3rd column in GTF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)
# -i IDATTR, --idattr IDATTR GTF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id)
# --additional-attr ADDITIONAL_ATTR Additional feature attributes (default: none, suitable for Ensembl GTF files: gene_name). Use multiple times for each different attribute
# -m {union,intersection-strict,intersection-nonempty}, --mode {union,intersection-strict,intersection-nonempty} Mode to handle reads overlapping more than one feature (choices: union, intersection-strict, intersection- nonempty; default: union)
# --nonunique {none,all} Whether to score reads that are not uniquely aligned or ambiguously assigned to features
# --secondary-alignments {score,ignore} Whether to score secondary alignments (0x100 flag)
# --supplementary-alignments {score,ignore} Whether to score supplementary alignments (0x800 flag)
# -o SAMOUT, --samout SAMOUT Write out all SAM alignment records into aSAM/BAM file, annotating each line with its feature assignment (as an optional field with tag 'XF'). See the -p option to use BAM instead of SAM.
# -p {SAM,BAM,sam,bam}, --samout-format {SAM,BAM,sam,bam} Format to use with the --samout option.
# -d OUTPUT_DELIMITER, --delimiter OUTPUT_DELIMITER Column delimiter in output (default: TAB).
# -c OUTPUT_FILENAME, --counts_output OUTPUT_FILENAME TSV/CSV filename to output the counts to instead of stdout.
# --cell-barcode CB_TAG BAM tag used for the cell barcode (default compatible with 10X Genomics Chromium is CB).
# --UMI UB_TAG BAM tag used for the unique molecular identifier, also known as molecular barcode (default compatible with
#10X Genomics Chromium is UB).
# -q, --quiet Suppress progress report
# --version Show software version and exit
# }}}

# {{{ rule feature_counts:
#     input: "07_deduplicated/{file_names}.deduplicated.bam"
#     output: "08_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
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

# options:
# -a gtf file input
# -o output file count matrix
# -R input file format
# -g ninth column of gtf file. the attribute type by which to group features. look to gtf file for correct naming. could also be transcript_id or exon_id
# -t third column of gtf file. indicates what type of features should be considered. When using -g gene_id should use gene. If using transcript_id or exon_id, should use processed transcript.
# -M multi-mapping reads will counted. Uses NH tag for multimappers.
# -O assign reads to all overlaping features
# --fraction Assign fractional counts to features
# }}}

rule multiqc:
    input: 
        inDir=".",
        fastqcReport = expand("02-fastqc/{fileFull}.std_err.txt", fileFull=fileN),
        processedReads = expand("08_feature_count/{file_names}.deduplicated.bam.featureCounts.tsv", file_names = sampleN)
    output:
        outFile="09_multiqc/multiqc_report.html"
    shell:
        """
        multiqc {input.inDir} \\
        -o 09_multiqc
        """
# {{{ extra, TODO: delete 
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
# }}}
