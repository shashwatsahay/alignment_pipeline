# check if Mosdepth is run with in Exome/Panel or WGS mode

if seq_type in ["Panel", "WES"]:

    rule mosdepth:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.cram",
            bai=wrkdir / "alignments" / "{sample}_dedup.recall.cram.crai",
            target_regions=target_regions,
            genome=genome,
        params:
            temp_idx=str(wrkdir / "alignments" / "{sample}_dedup.recall.crai"),
            prefix=str(wrkdir / "metrics" / "{sample}"),
        output:
            out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
            out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/mosdepth.yaml"
        log:
            logdir / "mosdepth/{sample}.log",
        message:
            "Running mosdepth for WES/panel data"
        shell:
            "cp {input.bai} {params.temp_idx} && "
            "mosdepth -f {input.genome} --by {input.target_regions} -n {params.prefix} {input.bam} &> {log}"

else:

    rule mosdepth:
        input:
            bam=wrkdir / "alignments" / "{sample}_dedup.recall.cram",
            genome=genome,
        params:
            prefix=str(wrkdir / "metrics" / "{sample}"),
        output:
            out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
            out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
        threads: 1
        resources:
            mem_mb=8000,
            runtime=24 * 60,
            nodes=1,
            tmpdir=scratch_dir,
        conda:
            "../envs/mosdepth.yaml"
        message:
            "Running mosdepth for WGS data"
        log:
            logdir / "mosdepth/{sample}.log",
        shell:
            "mosdepth -n {params.prefix} {input.bam} &> {log}"
