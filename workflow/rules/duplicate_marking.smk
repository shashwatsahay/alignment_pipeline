rule duplicates:
    input:
        bam=wrkdir / "alignments" / "{sample}.merged.sorted.bam",
        bai=wrkdir / "alignments" / "{sample}.merged.sorted.bam.bai",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}_dedup.bam"),
    conda:
        "../envs/sambamba.yaml"
    threads: 8
    log:
        logdir / "gatk/{sample}_dedup.log",
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    message:
        "Marking duplicates"
    shell:
        " ( "
        " sambamba markdup -t {threads} --tmpdir {resources.tmpdir} {input.bam} {output.bam}"
        " ) &> {log} "
