import pandas as pd 
import os
import yaml

## CONFIG -----------------------------------------------------------------------------
# Load configuration from config.yaml
config = yaml.safe_load(open("config.yaml"))

# Extract configuration parameters
LD_blocks_file = config["ld_blocks_file"]
ped_file = config["ped_file"]
POPULATIONS = config["populations"]
output_bcf_dir = config["output_bcf_dir"]
output_tsv_dir = config["output_tsv_dir"]

# Extract VCF files using pattern from config
vcf_pattern = os.path.join(config["input_dir"], config["vcf_pattern"])
CHROMS, = glob_wildcards(vcf_pattern)

CHROMS = list(set(CHROMS) - set(["chrM", "chrY", "others"]))  # Exclude some chromosomes

# Quality control parameters
MIN_MAF = config["min_maf"]
MIN_MAC = config["min_mac"]
MAX_MISSING = config["max_missing"]

# Environment paths
# CONDA = config["conda_path"]
# MAMBA = config["mamba_path"]
BCFTOOLS_PLUGINS = config["bcftools_plugins"]

# -------------------------------------------------------------------------------------
if not os.path.isfile(LD_blocks_file):
    IOError(f"LD blocks file {LD_blocks_file} does not exist")


LD_blocks = pd.read_table(LD_blocks_file)
print(LD_blocks.head())
LD_blocks["block"] = LD_blocks["chr"].astype(str) + "_" + LD_blocks["start"].astype(str) + "_" + LD_blocks["end"].astype(str)

LD_blocks["bcf"] = output_bcf_dir + "/" + LD_blocks["chr"].astype(str) + ".recalibrated_variants.bcf"


rule all:
    input:
        expand(os.path.join(output_bcf_dir, "{CHROM}.recalibrated_variants.bcf"), CHROM = CHROMS),
        expand(os.path.join(output_bcf_dir, "{POP}.samples"), POP = POPULATIONS),
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "filtered.bcf"), BLOCK = LD_blocks.block.values),
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "filtered_{POP}.bcf"), BLOCK = LD_blocks.block.values, POP = POPULATIONS),
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "filtered.bed"), BLOCK = LD_blocks.block.values),
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "filtered_{POP}.bed"), BLOCK = LD_blocks.block.values, POP = POPULATIONS),
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "variant_list.tsv"), BLOCK = LD_blocks.block.values),
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "variant_list_{POP}.tsv"), BLOCK = LD_blocks.block.values, POP = POPULATIONS),
        expand(os.path.join(output_tsv_dir, "variant_list_{POP}.tsv"), POP = POPULATIONS)

rule bcf:
    input:
        vcf = os.path.join(config["input_dir"], config["vcf_pattern"].replace("{ID}", "{CHROM}"))
    output:
        bcf = os.path.join(output_bcf_dir, "{CHROM}.recalibrated_variants.bcf")
    threads: 2 
    resources:
        mem = "2G",
        partition = "shared",
        time = "15:30:00"
    shell:
        """
        # module load bcftools
        bcftools view -m2 -M2 -v snps -f PASS -i "MAC >= {MIN_MAC} && F_MISSING < {MAX_MISSING}" -O b {input.vcf} > {output.bcf}
        bcftools index {output.bcf}
        """

rule pop_samples:
    input:
        ped = ped_file
    output:
        os.path.join(output_bcf_dir, "{POP}.samples")
    threads: 1
    resources:
        mem = "1G",
        partition = "express",
        time = "01:30:00"
    shell:
        """
        tail -n+2 {input.ped} | grep {wildcards.POP} | cut -f2 -d' ' > {output}
        """

rule filter:
    input:
        expand(os.path.join(output_bcf_dir, "{CHROM}.recalibrated_variants.bcf"), CHROM = CHROMS)
    output:
        bcf = os.path.join(output_bcf_dir, "{BLOCK}", "filtered.bcf")
    threads: 2 
    resources:
        mem = "2G",
        partition = "shared",
        time = "05:30:00"
    params:
        CHROM = lambda wildcards: LD_blocks.chr[LD_blocks.block == wildcards.BLOCK].values[0],
        START = lambda wildcards: LD_blocks.start[LD_blocks.block == wildcards.BLOCK].values[0],
        END = lambda wildcards: LD_blocks.end[LD_blocks.block == wildcards.BLOCK].values[0],
        BCF = lambda wildcards: LD_blocks.bcf[LD_blocks.block == wildcards.BLOCK].values[0]
    shell:
        """
        mkdir -p $(dirname {output.bcf})
        # module load bcftools
        bcftools view -r {params.CHROM}:{params.START}-{params.END} -O u {params.BCF} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b > {output.bcf}
        bcftools index {output.bcf}
        """

rule filter_pop:
    input:
        bcf = expand(os.path.join(output_bcf_dir, "{CHROM}.recalibrated_variants.bcf"), CHROM = CHROMS),
        samples = os.path.join(output_bcf_dir, "{POP}.samples")
    output:
        bcf = os.path.join(output_bcf_dir, "{BLOCK}", "filtered_{POP}.bcf")
    threads: 2 
    resources:
        mem = "2G",
        partition = "shared",
        time = "05:30:00"
    params:
        CHROM = lambda wildcards: LD_blocks.chr[LD_blocks.block == wildcards.BLOCK].values[0],
        START = lambda wildcards: LD_blocks.start[LD_blocks.block == wildcards.BLOCK].values[0],
        END = lambda wildcards: LD_blocks.end[LD_blocks.block == wildcards.BLOCK].values[0],
        BCF = lambda wildcards: LD_blocks.bcf[LD_blocks.block == wildcards.BLOCK].values[0]
    shell:
        """
        mkdir -p $(dirname {output.bcf})

        # export BCFTOOLS_PLUGINS={BCFTOOLS_PLUGINS}

        bcftools view -r {params.CHROM}:{params.START}-{params.END} -S {input.samples} -O u {params.BCF} | bcftools +fill-tags -O u | bcftools view -i "MAF > {MIN_MAF}" -O u | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b > {output.bcf}
        bcftools index {output.bcf}
        """
        
rule convert_to_plink:
    input:
        bcf = os.path.join(output_bcf_dir, "{BLOCK}", "filtered.bcf")
    output:
        variant_file = os.path.join(output_bcf_dir, "{BLOCK}", "variant_list.tsv"),
        bed = os.path.join(output_bcf_dir, "{BLOCK}", "filtered.bed")
    threads: 2 
    params:
        prefix = lambda wildcards: f"{output_bcf_dir}/{wildcards.BLOCK}/filtered"
    resources:
        mem = "2G",
        partition = "shared",
        time = "05:30:00"
    shell:
        """
        # module load bcftools
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input.bcf} > {output.variant_file}
        # plink2 does not reorder alleles, unlike plink1.9
        plink2 --bcf {input.bcf} --make-bed --out {params.prefix}
        """

rule convert_to_plink_pop:
    input:
        bcf = os.path.join(output_bcf_dir, "{BLOCK}", "filtered_{POP}.bcf")
    output:
        variant_file = os.path.join(output_bcf_dir, "{BLOCK}", "variant_list_{POP}.tsv"),
        bed = os.path.join(output_bcf_dir, "{BLOCK}", "filtered_{POP}.bed")
    threads: 2 
    params:
        prefix = lambda wildcards: f"{output_bcf_dir}/{wildcards.BLOCK}/filtered_{wildcards.POP}"
    resources:
        mem = "2G",
        partition = "shared",
        time = "05:30:00"
    shell:
        """
        # module load bcftools
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input.bcf} > {output.variant_file}
        # plink2 does not reorder alleles, unlike plink1.9
        plink2 --bcf {input.bcf} --make-bed --out {params.prefix}
        """

rule combine_variant_lists:
    input:
        variant_files = expand(os.path.join(output_bcf_dir, "{BLOCK}", "variant_list_{POP}.tsv"), BLOCK = LD_blocks.block.values, POP = POPULATIONS)
    output:
        variant_list = os.path.join(output_tsv_dir, "variant_list_{POP}.tsv")
    params:
        dirname = output_tsv_dir
    threads: 1 
    resources:
        mem = "1G",
        time = "01:30:00"
    shell:
        """
        mkdir -p {params.dirname}
        echo -e "CHROM\tPOS\tID\tREF\tALT" > {output.variant_list}
        for file in {input.variant_files}; do
            tail -n +2 "$file" >> {output.variant_list}
        done
        """
