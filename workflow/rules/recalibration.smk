rule baseRecalibrator:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.bam",
        dbsnp=dbsnp,
        genome=genome,
    output:
        table=wrkdir / "metrics" / "{sample}_recal_data.table",
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.cram",
        bai=wrkdir / "alignments" / "{sample}_dedup.recall.cram.crai",
        analyse_covariates=wrkdir / "metrics" / "{sample}_covariates.pdf",
    conda:
        "../envs/gatk.yaml"
    threads: 8
    resources:
        mem_mb=50000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "gatk/{sample}_recal.log",
    message:
        "Recalibrating with GATK BaseRecalibrator"
    shell:
        " ( "
        " gatk BaseRecalibrator -I {input.bam} -R {input.genome} "
        " --known-sites {input.dbsnp} "
        " -O {output.table} && "
        " gatk ApplyBQSR -I {input.bam} -R {genome} --bqsr-recal-file {output.table} -O {output.bam} && "
        " gatk AnalyzeCovariates "
        " -bqsr {output.table} "
        " -plots {output.analyse_covariates} && "
        " samtools index {output.bam}"
        " ) &> {log} "
        # Please refer to issue https://github.com/broadinstitute/gatk/issues/5299 with regards to mv command
        # in-short gatk can't create a crai extension and only support bai
        # further note samtools index command has to be used to because renaming bai to crai doesnt work for mosdepth
