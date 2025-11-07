rule sort_index:
    input:
        bam=wrkdir / "alignments" / "{sample}.bam",
    output:
        bam=temp(wrkdir / "alignments" / "{sample}.sorted.bam"),
        bai=temp(wrkdir / "alignments" / "{sample}.sorted.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    threads: 8
    params:
        mem_thread=8000,
    resources:
        mem_mb=9 * 8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "samtools/{sample}_sort.log",
    message:
        "Sorting and indexing recalibrated bam file"
    shell:
        " ( "
        " samtools sort --threads 8 -m{params.mem_thread}m -o {output.bam}##idx##{output.bai} {input.bam} -T {resources.tmpdir} --write-index"
        " ) &> {log} "
