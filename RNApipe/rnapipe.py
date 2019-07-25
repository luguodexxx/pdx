import os
import glob
configfile: "config.allinone.yaml"

# hINDEX = "/mnt/data3/zhouran/hg19_mm10_star/"
# hGENOME = "/mnt/data3/zhouran/hg19_mm10/hg19_mm10.fa"
# hg19bed = '/mnt/data3/zhouran/hg19_mm10/hg19.bed'
# mINDEX = "/mnt/data3/zhouran/mm10/genome"
# mGENOME = "/mnt/data3/zhouran/mm10/genome.fa"

INTERVAL = '/mnt/data3/zhouran/hg38/hg38.exon_interval.bed'
fqdir = config['fq']
labels = []
TMP = "--java-options \"-Djava.io.tmpdir=tmp\""


# SPECIES = ['hg19','mm10']
SPECIES = ['hg19']
for i in glob.glob(f'{fqdir}/*'):
    labels.append(os.path.split(i)[1].split('.')[0])
# labels = labels[0]
# labels = labels[:4]
# labels = ['WXF144_V2TC_T_EM_b2','WXF121_V2TC_T_EM_b1','WXF150_V2TC_T_EM_b1','T7_V2TC_T_1M_b1']
rule all:
	input:expand("data/{species}/{label}.star_{species}.sort.fixmate.paired.bam",label = labels,species=SPECIES),
		expand("data/{species}/{label}.r2.fq.gz",label=labels,species=SPECIES),
		expand("data/STAR/{label}.star_hybrid.bam",label=labels),
		expand("data/final/{species}/{label}.bam",species=SPECIES,label=labels),
		expand("data/rsem/{species}/{label}.transcript.bam",species=SPECIES,label=labels)

rule star_align:
	input:
		fq1 = config['fq'] + '/{label}.r1.fq.gz',
		fq2 = config['fq'] + '/{label}.r2.fq.gz'
	output:
		bam = "data/STAR/{label}.star_hybrid.bam"
	log:
		out = "log/STAR/star_align.{label}.log",
		err = "log/STAR/star_align.{label}.err"
	params:
		outprefix = "data/STAR/{label}.",
		indexdir = config['INDEX']['hybird']
	shell:  
		"""
		STAR \
		--genomeDir {params.indexdir} \
		--outFileNamePrefix {params.outprefix} \
		--readFilesIn {input.fq1} {input.fq2} \
		--readFilesCommand zcat \
		--runThreadN 6 \
		--outFilterScoreMinOverLread 0.33 \
		--outFilterMatchNminOverLread 0.33 \
		--outStd SAM \
		--genomeLoad LoadAndKeep \
		--outSAMunmapped Within | samtools view - -b -S -o {output.bam} 2> {log.err} 1> {log.out}
		"""
		# --chimSegmentMin 18 \
		# --chimScoreMin 12 \
		# --outFilterMultimapNmax 20 \
		# --outFilterMismatchNoverLmax 0.04 \
		# --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		# --alignIntronMax 200000 \
		# --outSAMstrandField intronMotif \

rule fetch_species:
	input:
		bam = "data/STAR/{label}.star_hybrid.bam"
	output:
		bam = "data/{species}/{label}.star_{species}.bam"
		
	log:
		err = "log/{species}/{label}.fetch_{species}.err"
	params:
		genome = lambda wildcards: config['GENOME'][wildcards.species],
		bed = lambda wildcards: config['BED'][wildcards.species]
	shell:
		"""
		samtools view -hL {params.bed} {input.bam}| \
		samtools view -bt {params.genome} - > {output.bam} 2> {log.err}
		"""
	
		# if wildcards.species == 'hg19':
		# 	shell("""
		# 	samtools view -hL {params.bed} {input.bam}| \
		# 	samtools view -bt {params.genome} - > {output.bam} 2> {log.err}
		# 	""")
		# else:
		# 	shell("""
		# 	samtools view -hL {params.bed} {input.bam} | \
		# 	samtools view -bt {params.genome} - > {output.bam} 2> {log.err}
		# 	""")



rule sortByname:
	input:
		bam = "data/{species}/{label}.star_{species}.bam"
	output:
		bam = temp("data/{species}/{label}.star_{species}.sortByname.bam")
	log:
		out = "log/{species}/{label}.sortByname.log",
		err = "log/{species}/{label}.sortByname.err"

	shell:
		"""
		sambamba sort -l 0 -t 8 -n --tmpdir ./tmp -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule fixmate:
	input:
		bam = "data/{species}/{label}star_{species}.sortByname.bam"
	output:
		bam = temp("data/{species}/{label}star_{species}.sortByname.fixmate.bam")
	log:
		out = "log/{species}/{label}.fixmate.log",
		err = "log/{species}/{label}.fixmate.err"
	shell:
		"""
		samtools view -F 256 -bh {input.bam} | samtools fixmate - {output.bam} 2> {log.err} 1> {log.out}
		"""

rule sort:
	input:
		bam = "data/{species}/{label}.star_{species}.sortByname.fixmate.bam"
	output:
		bam = temp("data/{species}/{label}.star_{species}.sort.fixmate.bam")
	log:
		out = "log/{species}/sort.{label}.log",
		err = "log/{species}/sort.{label}.err"

	shell:
		"""
		sambamba sort -l 0 -t 8 --tmpdir ./tmp -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""
rule fetch_paired:
	input:
		bam = "data/{species}/{label}.star_{species}.sort.fixmate.bam"
	output:
		bam1 = "data/{species}/{label}.star_{species}.sort.fixmate.unpaired.bam",
		bam2 = "data/{species}/{label}.star_{species}.sort.fixmate.paired.bam"
	log:
		err = "log/{species}/{label}.fetch_paired.err"
	shell:
		"""
		samtools view -F 1 -bh {input.bam} > {output.bam1} 2> {log.err}
		samtools view -f 1 -bh {input.bam} > {output.bam2} 2> {log.err}
		"""

# rule markdup_for_paired:
# 	input:
# 		bam = "data/{species}/{label}.star_hg19.sort.fixmate.paired.bam"
# 	output:
# 		bam = "data/{species}/{label}.star_hg19_final.bam"
# 	log:
# 		out = "log/{species}/{label}.process_bam7.log",
# 		err = "log/{species}/{label}.process_bam7.err"
# 	shell:
# 		"""
# 		sambamba markdup -t 8 -p --tmpdir ./tmp {input.bam} {output.bam} 2> {log.err} 1> {log.out}
# 		"""

rule indexbam:
	input:
		bam = "data/{species}/{label}.star_{species}.sort.fixmate.paired.bam"
	output:
		bai = "data/{species}/{label}.star_{species}.sort.fixmate.paired.bam.bai"
	log:
		out = "log/{species}/indexbam.{label}.log",
		err = "log/{species}/indexbam.{label}.err"
	shell:
		"""
		samtools index {input.bam} 2> {log.err} 1> {log.out}
		"""

rule bam2fastq:
	input:
		bam = "data/{species}/{label}.star_{species}.sort.fixmate.paired.bam",
		bai = "data/{species}/{label}.star_{species}.sort.fixmate.paired.bam.bai"
	output:
		fq1 = "data/{species}/{label}.r1.fq.gz",
		fq2 = "data/{species}/{label}.r2.fq.gz"
	log:
		out = "log/{species}/bam2fastq.{label}.log",
		err = "log/{species}/bam2fastq.{label}.err"
	shell:
		"""
		java -Xmx4g -jar /home/zhouran/data/soft/picard/picard.jar SamToFastq INPUT={input.bam} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2} 2> {log.err} 1> {log.out}
		"""

rule star_realign:
	input:
		fq1 = "data/{species}/{label}.r1.fq.gz",
		fq2 = "data/{species}/{label}.r2.fq.gz"
	output:
		bam = "data/final/{species}/{label}.bam",
		isobam = "data/final/{species}/{label}.Aligned.toTranscriptome.out.bam"
	log:
		out = "log/final/{species}/star_realign.{label}.log",
		err = "log/final/{species}/star_realign.{label}.err"
	params:
		# prefix = 'data/final/{species}/{label}.',
		# hgindex = lambda wildcards: config["INDEX"][wildcards.species]
		prefix = 'data/final/{species}/{label}.',
		hgindex = lambda wildcards: config["INDEX"][wildcards.species]
	shell:
		"""
		STAR \
		--genomeDir {params.hgindex} \
		--outFileNamePrefix {params.prefix} \
		--readFilesIn {input.fq1} {input.fq2} \
		--quantMode TranscriptomeSAM GeneCounts\
		--readFilesCommand zcat \
		--chimSegmentMin 18 \
		--chimScoreMin 12 \
		--runThreadN 6 \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNoverLmax 0.04 \
		--outFilterScoreMinOverLread 0.33 \
		--outFilterMatchNminOverLread 0.33 \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--alignIntronMax 200000 \
		--outSAMstrandField intronMotif \
		--chimJunctionOverhangMin 12 \
		--outSJfilterOverhangMin 12 12 12 12 \
		--outSJfilterCountUniqueMin 1 1 1 1 \
		--outSJfilterCountTotalMin 1 1 1 1 \
		--outStd SAM \
		--outSAMunmapped Within | samtools view - -b -S -o {output.bam} 2> {log.err} 1> {log.out}
		"""
		# --genomeLoad LoadAndKeep \

rule rsem:
	input:
		isobam = "data/final/{species}/{label}.Aligned.toTranscriptome.out.bam"
	output:
		bam  = "data/rsem/{species}/{label}.transcript.bam"
	params:
		prefix = 'data/rsem/{species}/{label}',
		rsemindex = lambda wildcards: config["RSEM"][wildcards.species]
	log:
		out = "log/rsem/{species}/rsem.{label}.log",
		err = "log/rsem/{species}/rsem.{label}.err"

	shell:
		"""
		rsem-calculate-expression --strandedness reverse \
		--paired-end --alignments --calc-pme --calc-ci --sort-bam-by-coordinate \
		-p 4 {input.isobam} \
		{params.rsemindex} {params.prefix} 2> {log.err} 1> {log.out}
		"""





