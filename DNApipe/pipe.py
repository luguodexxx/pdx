import os
import glob
hINDEX = "/mnt/data3/zhouran/hg19_mm10/hg19_mm10"
hGENOME = "/mnt/raid61/Microwell/hg19/refdata-cellranger-hg19-3.0.0/fasta/genome.fa"
hg38bed = '/mnt/data3/zhouran/hg19_mm10/hg19.bed'
# mINDEX = "/mnt/data3/zhouran/mm10/genome"
# mGENOME = "/mnt/data3/zhouran/mm10/genome.fa"

INTERVAL = '/mnt/data3/zhouran/hg38/hg38.exon_interval.bed'
SITE1 = "/mnt/data3/zhouran/publicbase/1000G_phase1.indels.b37.vcf.gz"
SITE2 = "/mnt/data3/zhouran/publicbase/dbsnp_138.b37.vcf.gz"
SITE3 = "/mnt/data3/zhouran/publicbase/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

labels = []
TMP = "--java-options \"-Djava.io.tmpdir=tmp\""
for i in glob.glob('/mnt/data3/zhouran/projtmp/fq/*'):
    labels.append(os.path.split(i)[1].split('.')[0])

labels = labels[:2]
reps = ['r1','r2']
rule all:
	input:
		expand('testnewmethod/bwa_aln/{label}.{rep}.fq.sai',label=labels,rep=reps),
		expand("testnewmethod/bwa_aln/{label}.final.hybrid_hg19.bam",label=labels),
		expand("testnewmethod/callsnp/gvcf/{label}.g.vcf.gz",label=labels)

rule bwa_aln:
    input:'/mnt/data3/zhouran/projtmp/fq/{label}.{rep}.fq.gz'
    output:'testnewmethod/bwa_aln/{label}.{rep}.fq.sai'
    log:'testnewmethod/bwa_aln/{label}.{rep}.log'
    shell:
        "bwa aln -t 8 {hINDEX} {input}>{output} 2> {log}"



# bwa sampe
rule bwa_sampe:
	input:
		'/mnt/data3/zhouran/projtmp/fq/{label}.r1.fq.gz','/mnt/data3/zhouran/projtmp/fq/{label}.r2.fq.gz',
		'testnewmethod/bwa_aln/{label}.r1.fq.sai','testnewmethod/bwa_aln/{label}.r2.fq.sai'
	output:
		'testnewmethod/bwa_aln/{label}.tmp.bam'
	log:
		out = "testnewmethod/bwa_aln/{label}_bwa_sampe.log",
		err = "testnewmethod/bwa_aln/{label}_bwa_sampe.err"
	params:
		sample = '{label}'
	# threads: 2
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
		bam = 'testnewmethod/bwa_aln/{label}.tmp.bam'
	output:
		bam = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam'
	log:
		out = "testnewmethod/bwa_aln/{label}_run_sambamba.log",
		err = "testnewmethod/bwa_aln/{label}_run_sambamba.err"

	shell:
		"""
		sambamba markdup -t 2 -p --tmpdir=./tmp {input.bam} {output.bam} 2> {log.err} 1> {log.out}
		"""

# # sambamba_index
rule sambamba_index:
	input:
		bam = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam'
	output:
		bai = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam.bai'
	log:
		out = "testnewmethod/bwa_aln/{label}_sambamba_index.log",
		err = "testnewmethod/bwa_aln/{label}_sambamba_index.err"

	shell:
		"""
		sambamba index -t 8 {input.bam} 2> {log.err} 1> {log.out}
		"""

# # validate bam
rule bam_validate:
	input:
		bai = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam.bai',
		bam = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam'
	output:
		txt = 'testnewmethod/bwa_aln/{label}_validate.txt'

	shell:
		"""
		bam validate --in {input.bam} --verbose > {output.txt} 2> {output.txt}
		"""

# align stats
# rule realign_stats:
# 	input:
# 		bam = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam'
# 	output:
# 		realignstats = bam = "testnewmethod/bwa_aln/{label}_hybrid_alignstats.txt"
# 	log:
# 		out = "testnewmethod/bwa_aln/{label}_realign_stats.log",
# 		err = "testnewmethod/bwa_aln/{label}_realign_stats.err"
# 	params:
# 		outdir = "testnewmethod/bwa_aln/"
# 		realignstats = config['dirs']['outdirs']['realignstatsdir'] + "{file}" + "/",
# 		vcrome_bed = config['data']['bed']['vcrome_bed']
# 	threads: 2
# 	shell:
# 		"""
# 		mkdir -p {params.realignstats}

# 		{params.alignstats} -v -i {input.bam} -o {output.realignstats} -t {params.vcrome_bed} -C -W 2> {log.err} 1> {log.out}
# 		"""

rule process_bam1:
	input:
		bam = 'testnewmethod/bwa_aln/{label}.tmp.markdups.bam'
	output:
		bam = "testnewmethod/bwa_aln/{label}.header_correct_hg19.bam"
	log:
		err = "testnewmethod/bwa_aln/{label}_process_bam1.err"

	shell:
		"""
		samtools view -hL {hg38bed} {input.bam} | sed "/@SQ\tSN:mm10*/d" |\
		samtools view -bt {hGENOME} - > {output.bam} 2> {log.err}
		"""

rule process_bam3:
	input:
		bam = "testnewmethod/bwa_aln/{label}.header_correct_hg19.bam"
	output:
		bam = "testnewmethod/bwa_aln/{label}.hg19_sorted_by_name.bam"
	log:
		out = "testnewmethod/bwa_aln/{label}.process_bam3.log",
		err = "testnewmethod/bwa_aln/{label}.process_bam3.err"

	shell:
		"""
		sambamba sort -l 0 -t 8 -n --tmpdir ./tmp -o {output.bam} {input.bam} 2> {log.err} 1> {log.out}
		"""

rule process_bam4:
	input:
		bam = "testnewmethod/bwa_aln/{label}.hg19_sorted_by_name.bam"
	output:
		bam = "testnewmethod/bwa_aln/{label}.hg19_sorted.bam"
	log:
		out = "testnewmethod/bwa_aln/{label}.process_bam4.log",
		err = "testnewmethod/bwa_aln/{label}.process_bam4.err"
	shell:
		"""
		samtools view -F 256 -Sbh {input.bam} | \
		samtools fixmate - - | \
		sambamba sort -l 0 -t 8 --tmpdir ./tmp -o {output.bam} /dev/stdin 2> {log.err} 1> {log.out}
		"""

rule process_bam5:
	input:
		bam = "testnewmethod/bwa_aln/{label}.hg19_sorted.bam"
	output:
		bam = "testnewmethod/bwa_aln/{label}.hybrid_unpaired_hg19.bam"
	log:
		err = "testnewmethod/bwa_aln/{label}.process_bam5.err"
	shell:
		"""
		samtools view -F 1 -Sbh {input.bam} > {output.bam} 2> {log.err}
		"""

rule process_bam6:
	input:
		bam = "testnewmethod/bwa_aln/{label}.hg19_sorted.bam"
	output:
		bam = "testnewmethod/bwa_aln/{label}.hybrid_paired_hg19.bam"
	log:
		err = "testnewmethod/bwa_aln/{label}.process_bam6.err"
	shell:
		"""
		samtools view -f 1 -Sbh {input.bam} > {output.bam} 2> {log.err}
		"""

rule process_bam7:
	input:
		bam = "testnewmethod/bwa_aln/{label}.hybrid_paired_hg19.bam"
	output:
		bam = "testnewmethod/bwa_aln/{label}.final.hybrid_hg19.bam"
	log:
		out = "testnewmethod/bwa_aln/{label}.process_bam7.log",
		err = "testnewmethod/bwa_aln/{label}.process_bam7.err"
	shell:
		"""
		sambamba markdup -t 8 -p --tmpdir ./tmp {input.bam} {output.bam} 2> {log.err} 1> {log.out} 
		"""

# rule HsMetrics:
#     input:
#     	bam = "testnewmethod/bwa_aln/{label}.final.hybrid_hg19.bam"
#     output:
#     	metrics="testnewmethod/bwa_aln/{label}.hsmetrics.hg19.txt"
#     shell:
#         "gatk CollectHsMetrics -BI {INTERVAL} -I {input.bam} -O {output.metrics} -TI {INTERVAL}"

rule baserecall:
    input:
        bam = "testnewmethod/bwa_aln/{label}.final.hybrid_hg19.bam"
    output:
        table = "testnewmethod/callsnp/{label}.final.hybrid_hg19.table"
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
        table = "testnewmethod/callsnp/{label}.final.hybrid_hg19.table",
        bam = "testnewmethod/bwa_aln/{label}.final.hybrid_hg19.bam"
    output:
        bam = "testnewmethod/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
    shell:
        "gatk ApplyBQSR \
        --bqsr-recal-file {input.table} \
        -I {input.bam} \
        -O {output.bam}"

rule gvcf:
    input:"testnewmethod/finalhgbam/{label}.final.hybrid_hg19.BQSR.bam"
    output:"testnewmethod/callsnp/gvcf/{label}.g.vcf.gz"
    shell:
        "gatk HaplotypeCaller -R {hGENOME} \
        -I {input} -O {output} -ERC GVCF --java-options \'-DGATK_STACKTRACE_ON_USER_EXCEPTION=true\'"
