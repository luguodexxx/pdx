import os
import glob
hINDEX = "/mnt/raid61/Personal_data/zhouran/reference/WES/hg19_mm10_bwa/hg19_mm10"
hGENOME = "/mnt/raid61/Microwell/hg19/refdata-cellranger-hg19-3.0.0/fasta/genome.fa"
hg38bed = '/mnt/raid61/Personal_data/zhouran/reference/WES/hg19_mm10_bwa/hg19.bed'
# mINDEX = "/mnt/data3/zhouran/mm10/genome"
# mGENOME = "/mnt/data3/zhouran/mm10/genome.fa"
# bwa = '/home/zhouran/data/soft/bwa/bwa-master/bwa'
INTERVAL = '/mnt/raid61/Personal_data/zhouran/reference/WES/hg19_mm10_bwa/S07604514_Covered.noChr.interval'
SITE1 = "/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/1000G_phase1.indels.b37.vcf.gz"
SITE2 = "/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/dbsnp_138.b37.vcf.gz"
SITE3 = "/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
Preprocessed_interval = "/mnt/raid61/Personal_data/zhouran/reference/WES/hg19_mm10_bwa/S07604514_Covered.noChr.preprocessed.interval.list"

hmap = '/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/hapmap_3.3.b37.vcf.gz'
omin = '/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/1000G_omni2.5.b37.vcf.gz'
G1000 = '/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/1000G_phase1.snps.high_confidence.b37.vcf.gz'
dbsnp = '/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/dbsnp_138.b37.vcf.gz'

germline = '/mnt/raid61/Personal_data/zhouran/reference/public_set/WES/af-only-gnomad.raw.sites.b37.vcf.gz'

TARGET = "/mnt/raid61/Personal_data/zhouran/reference/WES/hg19_mm10_bwa/S07604514_Covered.noChr.bed"

picard = "/mnt/raid61/Personal_data/zhouran/soft"
# labels = []
TMP = "--java-options \"-Djava.io.tmpdir=tmp\""
# for i in glob.glob('/mnt/data5/zhouran/WES/fq/first/*'):
#     labels.append(os.path.split(i)[1].split('.')[0])

NORMAL_tissue = "hNPC_WT_C"

labels = ["hNPC_TNP_C",
			"hNPC_V2TC_C",
			"hNPC_WT_C",
			"T37_TNP_S_2M",
			"T37_TNP_T_2M",
			"T46_TNP_T_2M",
			"T88_TNP_CC_2M",
			"T88_TNP_S_2M",
			"T88_TNP_T_2M",
			"T89_TNP_CC_2M",
			"T89_TNP_S_2M",
			"T89_TNP_T_2M",
			"T93_TNP_CC_2M",
			"T93_TNP_S_2M",
			"T93_TNP_T_2M",
			"T94_TNP_CC_2M",
			"T94_TNP_S_2M",
			"T94_TNP_T_2M",
			"WXF100_TNP_S_EM",
			"WXF100_TNP_T_EM",
			"WXF48_TNP_T_EM",
			"WXF49_TNP_T_EM",
			"WXF50_TNP_T_EM",
			"WXF83_TNP_T_EM",
			"WXF96_TNP_T_EM"
			]
# labels = [
# 	"T37_TNP_S_2M",
# 	"T37_TNP_T_2M",
# 	"T46_TNP_T_2M",
# 	"T88_TNP_CC_2M",
# 	"T88_TNP_S_2M",
# 	"T88_TNP_T_2M",
# 	"T89_TNP_CC_2M",
# 	"T89_TNP_S_2M",
# 	"T89_TNP_T_2M",
# 	"T93_TNP_CC_2M",
# 	"T93_TNP_S_2M",
# 	"T93_TNP_T_2M",
# 	"T94_TNP_CC_2M",
# 	"T94_TNP_S_2M",
# 	"T94_TNP_T_2M"
# ]

# labels = ["hNPC_TNP_C",
# 			"hNPC_V2TC_C",
			# "hNPC_WT_C"]
# 

reps = ['r1','r2']

ruleorder: bwa_aln>bwa_sampe>run_sambamba>sambamba_index>bam_validate>fetch_target>sortByName>FilterAndFix>FetchUnPair>FetchPair>PairMarkdup>baserecall>applyBQSR>gvcf
rule all:
	input:
		# expand('data/bwa_aln/{label}.{rep}.fq.sai',label=labels,rep=reps),
		# "data/callsnp/combined/all.vcf.gz",
		expand("data/varscan/{label}.varscan.snp",label=labels),
		expand("data/CNV/ModelSegments/{label}.cr.called.seg",label=labels)
		# expand("data/mutect2/hNPC_WT_C.{label}.mutect2.vcf",label=labels),
		# expand("data/callsnp/gvcf/{label}.g.vcf.gz",label=labels),
		# expand("data/bwa_aln/{label}.hsmetrics.hg19.txt",label=labels),
		# 'data/callsnp/combined/snps.threeWT.VQSR.vcf.gz',
		# 'data/callsnp/combined/snps.threeWT.VQSR.indel.vcf.gz'

rule bwa_aln:
	input:'/mnt/data5/zhouran/WES/fq/first/{label}.{rep}.fq.gz'
	output:temp('data/bwa_aln/{label}.{rep}.fq.sai')
	# log:
	# 	err = 'log/bwa_aln/{label}.{rep}_bwa_aln.err'
	resources:load=10
	shell:
		"bwa aln -t 6 {hINDEX} {input}>{output}"



# bwa sampe
rule bwa_sampe:
	input:
		'/mnt/data5/zhouran/WES/fq/first/{label}.r1.fq.gz','/mnt/data5/zhouran/WES/fq/first/{label}.r2.fq.gz',
		'data/bwa_aln/{label}.r1.fq.sai','data/bwa_aln/{label}.r2.fq.sai'
	output:
		temp('data/bwa_aln/{label}.tmp.bam')
	log:
		out = "log/bwa_aln/{label}_bwa_sampe.log",
		err = "log/bwa_aln/{label}_bwa_sampe.err"
	resources:load=1
	params:
		sample = '{label}'
	shell:
		"""
		bwa sampe -P -r '@RG\\tID:0\\tSM:{params.sample}\\tPL:Illumina' \
		{hINDEX} {input[2]} {input[3]} {input[0]} {input[1]} | \
		samtools view -F 256 -bu - | samtools fixmate - - | \
		sambamba sort -l 1 -t 8 --tmpdir ./tmp -o {output} /dev/stdin 2> {log.err} 1> {log.out}
		"""

# # sambamba
rule run_sambamba:
	input:
		bam = 'data/bwa_aln/{label}.tmp.bam'
	output:
		bam = 'data/bwa_aln/{label}.tmp.markdups.bam'
	log:
		out = "log/bwa_aln/{label}_run_sambamba.log",
		err = "log/bwa_aln/{label}_run_sambamba.err"
	resources:load=30
	shell:
		"""
		sambamba markdup -t 2 -p --tmpdir=./tmp {input.bam} {output.bam} 2> {log.err} 1> {log.out}
		"""

# # sambamba_index
rule sambamba_index:
	input:
		bam = 'data/bwa_aln/{label}.tmp.markdups.bam'
	output:
		bai = 'data/bwa_aln/{label}.tmp.markdups.bam.bai'
	log:
		out = "log/bwa_aln/{label}_sambamba_index.log",
		err = "log/bwa_aln/{label}_sambamba_index.err"
	resources:load=1
	shell:
		"""
		sambamba index -t 2 {input.bam} 2> {log.err} 1> {log.out}
		"""

# # validate bam
rule bam_validate:
	input:
		bai = 'data/bwa_aln/{label}.tmp.markdups.bam.bai',
		bam = 'data/bwa_aln/{label}.tmp.markdups.bam'
	output:
		txt = 'data/bwa_aln/{label}_validate.txt'
	resources:load=1
	shell:
		"""
		bam validate --in {input.bam} --verbose > {output.txt} 2> {output.txt}
		"""

rule fetch_target:
	input:
		bam = 'data/bwa_aln/{label}.tmp.markdups.bam'
	output:
		bam = temp("data/bwa_aln/{label}.header_correct_hg19.bam")
	log:
		err = "log/bwa_aln/{label}_fetch_target.err"
	resources:load=1
	shell:
		"""
		samtools view -hL {hg38bed} {input.bam} | sed "/@SQ\tSN:mm10*/d" |\
		samtools view -bt {hGENOME} - > {output.bam} 2> {log.err}
		"""

rule sortByName:
	input:
		bam = "data/bwa_aln/{label}.header_correct_hg19.bam"
	output:
		bam = temp("data/bwa_aln/{label}.hg19_sorted_by_name.bam")
	log:
		out = "log/bwa_aln/{label}.sortByName.log",
		err = "log/bwa_aln/{label}.sortByName.err"
	resources:load=1
	shell:
		"""
		sambamba sort -l 0 -t 2 -n --tmpdir ./tmp -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule FilterAndFix:
	input:
		bam = "data/bwa_aln/{label}.hg19_sorted_by_name.bam"
	output:
		bam = temp("data/bwa_aln/{label}.hg19_sorted.bam")
	log:
		out = "data/bwa_aln/{label}.FilterAndFix.log",
		err = "data/bwa_aln/{label}.FilterAndFix.err"
	resources:load=1
	shell:
		"""
		samtools view -F 256 -Sbh {input.bam} | \
		samtools fixmate - - | \
		sambamba sort -l 0 -t 2 --tmpdir ./tmp -o {output.bam} /dev/stdin 2> {log.err} 1> {log.out}
		"""

rule FetchUnPair:
	input:
		bam = "data/bwa_aln/{label}.hg19_sorted.bam"
	output:
		bam = "data/bwa_aln/{label}.hybrid_unpaired_hg19.bam"
	log:
		err = "data/bwa_aln/{label}.FetchUnPair.err"
	resources:load=1
	shell:
		"""
		samtools view -F 1 -Sbh {input.bam} > {output.bam} 2> {log.err}
		"""

rule FetchPair:
	input:
		bam = "data/bwa_aln/{label}.hg19_sorted.bam"
	output:
		bam = "data/bwa_aln/{label}.hybrid_paired_hg19.bam"
	log:
		err = "data/bwa_aln/{label}.FetchPair.err"
	resources:load=1
	shell:
		"""
		samtools view -f 1 -Sbh {input.bam} > {output.bam} 2> {log.err}
		"""

rule PairMarkdup:
	input:
		bam = "data/bwa_aln/{label}.hybrid_paired_hg19.bam"
	output:
		bam = "data/bwa_aln/{label}.final.hybrid_hg19.bam"
	log:
		out = "data/bwa_aln/{label}.PairMarkdup.log",
		err = "data/bwa_aln/{label}.PairMarkdup.err"
	resources:load=1
	shell:
		"""
		sambamba markdup -t 2 -p --tmpdir ./tmp {input.bam} {output.bam} 2> {log.err} 1> {log.out} 
		"""

# align stats
rule alignstats:
	input:
		"data/bwa_aln/{label}.final.hybrid_hg19.bam"
	output:
		# "data/alignstats/{label}.hybrid_hg19.txt"
	log:
		# out = "log/alignstats/{label}_realign_stats.log",
		# err = "log/alignstats/{label}_realign_stats.err"
	params:
	shell:
		"""
		# alignstats -v -i {input} -o {output} -t {TARGET} -C -W 2> {log.err} 1> {log.out}
		"""

rule HsMetrics:
    input:
    	bam = "data/bwa_aln/{label}.final.hybrid_hg19.bam"
    output:
    	metrics="data/bwa_aln/{label}.hsmetrics.hg19.txt"
    shell:
        "gatk CollectHsMetrics -BI {INTERVAL} -I {input.bam} -O {output.metrics} -TI {INTERVAL} --VALIDATION_STRINGENCY=LENIENT"

rule baserecall:
    input:
        bam = "data/bwa_aln/{label}.final.hybrid_hg19.bam"
    output:
        table = "data/callsnp/{label}.final.hybrid_hg19.table"
    resources:load=1
    # params:
    #     site1 = "{SITE1}",
    #     site2 = "{SITE2}",
    #     site3 = "{SITE3}"
    shell:
        "gatk BaseRecalibrator -R {hGENOME} \
        -I {input.bam}\
        --known-sites {SITE1} \
        --known-sites {SITE2} \
        --known-sites {SITE3} \
        -O {output.table}\
        "

rule applyBQSR:
    input:
        table = "data/callsnp/{label}.final.hybrid_hg19.table",
        bam = "data/bwa_aln/{label}.final.hybrid_hg19.bam"
    output:
        bam = "data/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
    resources:load=1
    shell:
        "gatk ApplyBQSR \
        --bqsr-recal-file {input.table} \
        -I {input.bam} \
        -O {output.bam}"

# """
# normal_pileup=\"/home/zhouran/data/miniconda3/bin/samtools mpileup -f {hGENOME} {params.normal} -\" && \
# tumor_pileup_{params.label_}=\"/home/zhouran/data/miniconda3/bin/samtools mpileup -f {hGENOME} {input} -\" && \
# /home/zhouran/data/miniconda3/bin/varscan somatic <(/bin/bash -c \"$normal_pileup\") <(/bin/bash -c \"$tumor_pileup_{params.label_}\") {params.prefix} 
# """

# normal_pileup=\"/home/zhouran/data/miniconda3/bin/samtools mpileup -f {hGENOME} {params.normal} -\" && \
# tumor_pileup_{params.label_}=\"/home/zhouran/data/miniconda3/bin/samtools mpileup -f {hGENOME} {input} -\" && \
# --min-coverage-normal 10 --min-coverage-tumor 14 --min-var-freq 0.07
# filter criterion: Spatial intratumoral heterogeneity and temporal clonal evolution in esophageal squamous cell carcinoma

rule gvcf:
    input:"data/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
    output:"data/callsnp/gvcf/{label}.g.vcf.gz"
   	resources:load=1
    shell:
        "gatk HaplotypeCaller -R {hGENOME} \
        -I {input} -O {output} -ERC GVCF --java-options \'-DGATK_STACKTRACE_ON_USER_EXCEPTION=true\'"

rule combined:
	input:expand("data/callsnp/gvcf/{label}.g.vcf.gz",label=labels)
	output:"data/callsnp/combined/threeWT.g.vcf.gz"
	# params:labels = lambda wildcards: wildcards.label
	run:
		inputs = ' -V '.join(input)
		shell("gatk CombineGVCFs -R {hGENOME} -V {inputs} -O {output}")


rule GenotypeGVCFs:
	input:"data/callsnp/combined/threeWT.g.vcf.gz"
	output:"data/callsnp/combined/threeWT.vcf.gz"
	shell:
		"gatk GenotypeGVCFs -R {hGENOME} -V {input} -O {output}"


######### here was for population analysis, for now we only care about the three normal cells
#############################################
# 											# 
#     here was for three normal cells		# 
# 											# 
#############################################

rule SNP_mod:
	input:
		"data/callsnp/combined/threeWT.vcf.gz"
	output:
		recalfile  = "data/callsnp/combined/threeWT.snps.recal",
		tranchefile = "data/callsnp/combined/threeWT.snps.tranches"
	params:
		rscriptfile = "data/callsnp/combined/threeWT.snps.plots.R"
	shell:
		"""
		gatk VariantRecalibrator \
		-R {hGENOME} \
		-V {input} \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hmap} \
		-resource:omin,known=false,training=true,truth=false,prior=12.0 {omin} \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 {G1000} \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
		-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
		-mode SNP \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
		--rscript-file {params.rscriptfile} \
		--tranches-file {output.tranchefile} \
		-O {output.recalfile}
		"""

rule SNP_apply:
	input:
		vcf = "data/callsnp/combined/threeWT.vcf.gz",
		recalfile  = "data/callsnp/combined/threeWT.snps.recal",
		tranchefile = "data/callsnp/combined/threeWT.snps.tranches"
	output:
		'data/callsnp/combined/snps.threeWT.VQSR.vcf.gz'
	shell:
		"""
		gatk ApplyVQSR \
		-R {hGENOME} \
		-V {input.vcf} \
		--recal-file {input.recalfile} \
		--tranches-file {input.tranchefile} \
		-mode SNP \
		-O {output}
		"""
rule indel_mod:
	input:
		'data/callsnp/combined/snps.threeWT.VQSR.vcf.gz'
	output:
		recalfile  = "data/callsnp/combined/snps.threeWT.indel.recal",
		tranchefile = "data/callsnp/combined/snps.threeWT.indel.tranches"
	params:
		rscriptfile = "data/callsnp/combined/snps.threeWT.indel.plots.R"
	shell:
		"""
		gatk VariantRecalibrator \
		-R {hGENOME} \
		-V {input} \
		-resource:mills,known=true,training=true,truth=true,prior=12.0 {SITE3} \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
		-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
		-mode INDEL \
		--max-gaussians 6 \
		--rscript-file  {params.rscriptfile} \
		--tranches-file {output.tranchefile} \
		-O {output.recalfile}
		"""
rule indel_apply:
	input:
		vcf = 'data/callsnp/combined/snps.threeWT.VQSR.vcf.gz',
		recalfile  = "data/callsnp/combined/snps.threeWT.indel.recal",
		tranchefile = "data/callsnp/combined/snps.threeWT.indel.tranches"
	output:
		'data/callsnp/combined/snps.threeWT.VQSR.indel.vcf.gz'
	shell:
		"""
		gatk ApplyVQSR \
		-R {hGENOME} \
		-V {input.vcf} \
		-ts-filter-level 99.0 \
		--recal-file {input.recalfile} \
		--tranches-file {input.tranchefile} \
		-mode INDEL \
		-O {output}
		"""

#############################################
# 											# 
#          here was for mutect2 		    # 
# 											# 
#############################################

rule mutect2:
	input:
		bam = "data/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
	output:
		vcf = "data/mutect2/hNPC_WT_C.{label}.mutect2.vcf",
		bam = "data/mutect2/hNPC_WT_C.{label}.bam"
	params:
		normal_bam = "data/finalhgbam/hNPC_WT_C.final.hybrid_hg19.BQSR.bam",
		normal_label = "hNPC_WT_C",
		tumor_label = "{label}"
	shell:
		"""
		gatk Mutect2 -R {hGENOME} \
		-I {params.normal_bam} -normal {params.normal_bam} \
		-I {input.bam} -tumor {params.tumor_label} \
		--germline-resource {germline} \
		--af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
		--bam-output {output.bam} -O {output.vcf} --native-pair-hmm-threads 4
		"""

#############################################
# 											# 
#          here was for CNV anlysis		    # 
# 											# 
#############################################

rule CNV_CollectReadCounts:
	input:
		bam = "data/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
	output:
		hdf5 = "data/CNV/CollectReadCounts/{label}.CollectReadCounts.hdf5"
	shell:
		"""
		gatk CollectReadCounts \
		-I {input.bam} \
		-L {Preprocessed_interval} \
		-R {hGENOME} \
		--format HDF5 \
		--interval-merging-rule OVERLAPPING_ONLY \
		--output {output.hdf5}
		"""
rule CNV_CollectAllelicCounts:
	input:
		bam = "data/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
	output:
		"data/CNV/CollectReadCounts/{label}.CollectAllelicCounts.txt"
	shell:
		"""
		gatk CollectAllelicCounts \
		-I {input.bam} \
		-L {Preprocessed_interval} \
		-R {hGENOME}
		--minimum-base-quality 20 \
		--output {output}
		"""

rule CNV_CreateReadCountPanelOfNormals:
	input:
		expand("data/CNV/CollectReadCounts/{label}.CollectAllelicCounts.txt",label=labels),
		expand("data/CNV/CollectReadCounts/{label}.CollectReadCounts.hdf5",label=labels)
	output:
		"data/CNV/PoN.hdf5"
	params:
		normal = lambda wildcards: "data/CNV/CollectReadCounts/{wildcards.NORMAL_tissue}.CollectReadCounts.hdf5"
	shell:
		"""
		gatk CreateReadCountPanelOfNormals --input {params.normal} --minimum-interval-median-percentile  5.0  \
		--output {output}
		"""

rule CNV_DenoiseReadCounts:
	input:
		hdf5 = "data/CNV/CollectReadCounts/{label}.CollectReadCounts.hdf5",
		PoN = "data/CNV/PoN.hdf5"
	output:
		standardizedCR = "data/CNV/DenoiseReadCounts/{label}.standardizedCR.txt",
		denoisedCR = "data/CNV/DenoiseReadCounts/{label}.denoisedCR.txt"
	shell:
		"""
		gatk DenoiseReadCounts --input {input.hdf5} \
		--count-panel-of-normals {input.PoN} \
		--standardized-copy-ratios {output.standardizedCR} \
		--denoised-copy-ratios {output.denoisedCR}
		"""
rule CNV_ModelSegments:
	input:
		standardizedCR = "data/CNV/DenoiseReadCounts/{label}.standardizedCR.txt",
		denoisedCR = "data/CNV/DenoiseReadCounts/{label}.denoisedCR.txt",
		allelicCounts = "data/CNV/CollectReadCounts/{label}.CollectAllelicCounts.txt"
	output:
		"data/CNV/ModelSegments/{label}/{label}.cr.seg"
	params:
		normal = lambda wildcards: "data/CNV/CollectReadCounts/{wildcards.NORMAL_tissue}.CollectAllelicCounts.txt",
		outputdir = "data/CNV/ModelSegments/{label}",
		prefix = "{label}"
	shell:
		"""
		gatk ModelSegments --denoised-copy-ratios {input.denoisedCR} \
		--allelic-counts {input.allelicCounts} \
		--normal-allelic-counts {params.normal} \
		--output {params.outputdir} --output-prefix {params.prefix} \
		
		"""

# gatk --java-options "-Xmx5000m" ModelSegments --denoised-copy-ratios samplename.denoisedCR.tsv  \
# --allelic-counts samplename_T.allelic_counts --normal-allelic-counts samplename_N.allelic_counts --output segments --output-prefix samplename \
###segments 需要存在，否则报错 A USER ERROR has occurred: Output directory segments does not exist.

rule CNV_CallCopyRatioSegments:
	input:
		"data/CNV/ModelSegments/{label}/{label}.cr.seg"
	output:
		"data/CNV/ModelSegments/{label}.cr.called.seg"
	shell:
		"""
		gatk CallCopyRatioSegments -I {input} -O {output}
		"""

# time gatk --java-options "-Xmx5000m"  CallCopyRatioSegments -I segments/samplename.cr.seg -O segments/samplename.called.seg


#############################################
# 											# 
#          one line for varscan  		    # 
# 											# 
#############################################

rule bam_mpileup:
	input:
		"data/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
	output:
		"data/varscan/{label}.varscan.snp"
	params:
		normal = "data/finalhgbam/hNPC_WT_C.final.hybrid_hg19.BQSR.bam",
		prefix = "data/varscan/{label}.varscan",
	shell:
		"""
		samtools mpileup -f {hGENOME} -q 1 -B {params.normal} {input}| varscan somatic - {params.prefix} --mpileup 1
		"""

