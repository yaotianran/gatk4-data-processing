configfile: "config.yaml"

sample = config['sample']

rule all:
    input:
        expand('{sample}/7_recal.bam', sample = sample)


rule map:
    input:
        R1 = config['fastq']['R1'],
        R2 = config['fastq']['R2'],
        reference = config['support']['reference']
    params:
        rg = expand(r'@RG\tID:{sample}\tSM:{sample}', sample = sample)
    output:
        expand('{sample}/2_raw.sam', sample = sample)
    shell:
        "bwa mem -K 10000000 -P -o {output} -R '{params.rg}' {input.R1} {input.R2} {input.reference}"

rule sort_queryname:
    input:
        expand('{sample}/2_raw.sam', sample = sample)
    output:
        expand('{sample}/3_queryname_sorted.bam', sample = sample)
    shell:
        "samtools sort -n -o {output} {input}"


rule mark_duplicates:
    input:
        expand('{sample}/3_queryname_sorted.bam', sample = sample)
    output:
        expand('{sample}/4_dedup.bam', sample = sample)
    shell:
        "gatk MarkDuplicates -I {input} -O {output} --METRICS_FILE {sample}/MarkDuplicates.metrics.txt --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER 'queryname'"


rule sort_position:
    input:
        expand('{sample}/4_dedup.bam', sample = sample)
    output:
        expand('{sample}/5_position_sorted.bam', sample = sample)
    shell:
        "gatk SortSam -I {input} -O {output} --SORT_ORDER 'coordinate' --CREATE_INDEX false --CREATE_MD5_FILE false"


rule fix_tags:
    input:
        bam = expand('{sample}/5_position_sorted.bam', sample = sample), 
        reference = config['support']['reference']
    output:
        expand('{sample}/6_fixed.bam', sample = sample)
    shell:
        "gatk SetNmMdAndUqTags -I {input.bam} -O {output} --REFERENCE_SEQUENCE {input.reference} --CREATE_INDEX true"


rule calculate_recalibration:
    input:
        bam = expand('{sample}/6_fixed.bam', sample = sample), 
        reference = config['support']['reference'], 
        indel = config['support']['known_sites_indel'], 
        snp = config['support']['known_sites_snp']
    output:
        expand('{sample}/7_recal.table', sample = sample)
    shell:
        "gatk BaseRecalibrator --use-original-qualities -I {input.bam} -O {output} -R {input.reference} --known-sites {input.indel} --known-sites {input.snp}"


rule apply_recalibration:
    input:
        bam = expand('{sample}/6_fixed.bam', sample = sample), 
        reference = config['support']['reference'], 
        table = expand('{sample}/7_recal.table', sample = sample)
    output:
        expand('{sample}/7_recal.bam', sample = sample)
    shell:
        "gatk ApplyBQSR -bqsr {input.table} -I {input.bam} -O {output} -R {input.reference} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --use-original-qualities --emit-original-quals true"
