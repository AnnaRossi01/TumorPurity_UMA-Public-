import pandas as pd 

# ---- - Configurations ---- #
configfile: "config/config.yaml"

# ---- Manifest parsing ---- #
manifest_df = pd.read_csv("config/manifest.csv")
manifest_df.set_index("sample_id", inplace=True)
print("Manifest loaded successfully. Available samples:")
print(manifest_df.index.tolist())

SAMPLES = manifest_df.index.tolist()

# ---- VARIABLES ---- #
BAM = manifest_df["bam"].to_dict()
VCF = manifest_df["vcf"].to_dict()
ICHOR = manifest_df["ichor_res"].to_dict()

PURECN_SIF = f"{config['singularity_image']}/purecn_2.10.0-amd64.sif"
ASCAT_SIF  = f"{config['singularity_image']}/ascatsc_env_v3.sif"

# ---- - Main rule ---- #
rule all:
    input:
        "results/purity_comparison/purity_comparison.csv",
        expand("results/ascat_sc/{sample}/ascat.done", sample=SAMPLES),
        expand("results/purecn/{sample}/{sample}.csv", sample=SAMPLES)


# ----ASCAT.sc rule ---- #
rule ascat_sc:
    input: 
        tumor = lambda wildcards: BAM[wildcards.sample],
        normal = config['normal']
    output:
        out_dir = directory("results/ascat_sc/{sample}"),
        done    = touch("results/ascat_sc/{sample}/ascat.done")   # .done is a common convention to indicate the completion of a task
    params: 
        bed = config['bed_file']
    log: 
        "logs/ascat_sc/{sample}.log"
    resources:
        mem_mb = 10000,
        cpus_per_task = 4
    singularity: 
        ASCAT_SIF
    shell: 
        """
        mkdir -p {output.out_dir}
        Rscript workflow/scripts/ascat/run_ASCAT.R \
            {wildcards.sample} \
            {input.tumor} \
            {input.normal} \
            {output.out_dir} \
            {params.bed} 2> {log}
        """ 



# ---- PureCN rule ---- #

# Phase 1: Generate intervals for PureCN
rule purecn_intervals:
    input:
        baits = config['baits_clean'],
        ref   = config['refGen']
    output:
        intervals = "resources/purecn_assets/baits_optimized_hg19.txt"
    params: 
        mappability = config ['mappability']
    singularity:
        PURECN_SIF
    log:
        "logs/purecn/intervals_generation.log"
    shell: 
        """
        Rscript $PURECN/IntervalFile.R \
            --in-file {input.baits} \
            --fasta {input.ref} \
            --genome hg19 \
            --mappability {params.mappability} \
            --off-target \
            --out-file {output.intervals} > {log} 2>&1
        """


# Phase 2: Coverage and NormalDB
rule purecn_normalcoverage:
    input:
        intervals = "resources/purecn_assets/baits_optimized_hg19.txt" 
    output:
        out_dir = directory("resources/purecn_assets/coverage"), 
        done    = touch("resources/purecn_assets/coverage/normal_coverage.done")
    params: 
        normals_dir = "/mnt/drive02/PureCN/input_file/bam/normals",
        bam_list = "resources/purecn_assets/normals_bam.list"
    log:
        "logs/purecn/pon_coverage_batch.log"
    resources:
        mem_mb = 32000,     
        cpus_per_task = 8  
    singularity:
        PURECN_SIF
    shell: 
        """
        mkdir -p {output.out_dir}
        
        # 1. Dynamically recreate the BAM list as in your bash script
        find {params.normals_dir} -name "*.bam" > {params.bam_list}
        
        # 2. Exeution
        Rscript $PURECN/Coverage.R \
            --out-dir {output.out_dir} \
            --bam {params.bam_list} \
            --intervals {input.intervals} \
            --cores {resources.cpus_per_task} > {log} 2>&1
        """

# Phase 3: NormalDB 
rule purecn_normaldb:
    input:
        coverage_normal = "resources/purecn_assets/coverage/normal_coverage.done"
    output:
        normaldb = "resources/purecn_assets/normalDB_panel_aio_hg19.rds"
    params: 
        # N.B. Assicurati che questa cartella sia esattamente quella dove la regola precedente ha salvato i .txt.gz
        cov_dir = "resources/purecn_assets/coverage", 
        out_dir = "resources/purecn_assets/",
        cov_list = "resources/purecn_assets/coverage_list.txt"
    log:
        "logs/purecn/normaldb.log"
    resources:
        mem_mb = 16000       
    singularity:
        PURECN_SIF
    shell: 
        """
        # 1. create the list file for NormalDB.R
        find {params.cov_dir} -name "*_coverage.txt*" > {params.cov_list}
        
        # 2. execute NormalDB.R
        Rscript $PURECN/NormalDB.R \
            --out-dir {params.out_dir} \
            --coverage-files {params.cov_list} \
            --genome hg19 \
            --assay panel_aio > {log} 2>&1
        """


# Phase 4: Tumor coverage for each sample
rule purecn_tumor_coverage:
    input:
        tumor_bam = lambda wildcards: BAM[wildcards.sample],
        intervals = "resources/purecn_assets/baits_optimized_hg19.txt"
    output:
        tumor_cov = "results/purecn/coverage/{sample}/{sample}_coverage.txt.gz"
    log:
        "logs/purecn/{sample}_coverage.log"
    resources:
        mem_mb = 10000
    singularity:
        PURECN_SIF
    shell:
        """
        OUTDIR=$(dirname {output.tumor_cov})
        mkdir -p $OUTDIR
        
        Rscript $PURECN/Coverage.R \
            --out-dir $OUTDIR \
            --bam {input.tumor_bam} \
            --intervals {input.intervals} > {log} 2>&1
            
        # rename the output file to match the expected output
        BAM_BASENAME=$(basename {input.tumor_bam} .bam)
        mv $OUTDIR/${{BAM_BASENAME}}_coverage.txt.gz {output.tumor_cov}
        """

# Phase 5: Execute PureCN for each sample
rule purecn_run:
    input:
        tumor_cov = "results/purecn/coverage/{sample}/{sample}_coverage.txt.gz",
        vcf       = lambda wildcards: VCF[wildcards.sample],
        intervals = "resources/purecn_assets/baits_optimized_hg19.txt",
        normaldb  = "resources/purecn_assets/normalDB_panel_aio_hg19.rds",
    output:
        out_dir = directory("results/purecn/{sample}"),
        csv_res = "results/purecn/{sample}/{sample}.csv"
    params:
        sample_id = "{sample}",
        genome    = "hg19"
    log:
        "logs/purecn/{sample}_run.log"
    resources:
        mem_mb = 16000,
        cpus_per_task = 4
    singularity:
        PURECN_SIF
    shell:
        """
        # FOLDER
        mkdir -p {output.out_dir}

        # Eseguiamo PureCN e forziamo Bash a ignorare il suo exit code con '|| true'
        Rscript $PURECN/PureCN.R \
            --out {output.out_dir}/{params.sample_id} \
            --tumor {input.tumor_cov} \
            --sampleid {params.sample_id} \
            --vcf {input.vcf} \
            --normaldb {input.normaldb} \
            --intervals {input.intervals} \
            --genome {params.genome} \
            --model betabin \
            --fun-segmentation PSCBS \
            --postoptimize \
            --force \
            --seed 123 \
            --cores {resources.cpus_per_task} >> {log} 2>&1 || true

        # Ora controlliamo solo se il file CSV di output esiste
        if [ ! -f {output.csv_res} ]; then
            echo "[FALLBACK] PureCN non ha prodotto il file CSV (non-convergenza o crash)." >> {log}
            echo "[FALLBACK] Creazione del file dummy con valori NA per {params.sample_id}..." >> {log}
            
            # Creiamo il file
            echo "Sample.ID,purity,ploidy" > {output.csv_res}
            echo "{params.sample_id},NA,NA" >> {output.csv_res}
            
            echo "[FALLBACK] File dummy creato correttamente." >> {log}
        else
            echo "[SUCCESS] PureCN ha generato correttamente il file CSV." >> {log}
        fi
        
        # Obblighiamo Bash a chiudere il blocco con successo, in modo che Snakemake 
        # creda che il job sia andato a buon fine (e quindi conservi il file CSV dummy)
        exit 0
        """

# tumor purity comaprison 
rule purity_comparison: 
    input: 
        ascat_done = expand("results/ascat_sc/{sample}/ascat.done", sample=SAMPLES),
        purecn_csv = expand("results/purecn/{sample}/{sample}.csv", sample=SAMPLES), 
        manifest = "config/manifest.csv",
        script = "workflow/scripts/purity_comparisonOOP.py"
    output:
        report = "results/purity_comparison/purity_comparison.csv"
    params: 
        flow_db = "/mnt/drive02/Tumor_purity_pipeline/flow_cytometry_db.xlsx",
        ascat_dir = "results/ascat_sc/",
        purecn_dir = "results/purecn/"
    log: 
        "logs/purity_comparison.log"
    conda: 
        "workflow/envs/python_env.yaml" 
    shell: 
        """
        python workflow/scripts/purity_comparisonOOP.py \
            --manifest {input.manifest} \
            --flow_db {params.flow_db} \
            --ascat_dir {params.ascat_dir} \
            --purecn_dir {params.purecn_dir} \
            --output {output.report} > {log} 2>&1
        """

