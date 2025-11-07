rule flagstatt:
    input:
        bam=wrkdir / "alignments" / "{sample}_{ext}.cram",
    output:
        wrkdir / "metrics" / "{sample}_{ext}.flagstat",
    conda:
        "../envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "sambamba/{sample}_{ext}.log",
    message:
        "Running Flagstat"
    shell:
        "(samtools flagstat {input.bam} > {output}) &> {log}"
