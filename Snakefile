import pandas as pd 
import os

## CONFIG -----------------------------------------------------------------------------
# one line per region
LD_blocks_file = "/data/abattle4/april/hi_julia/LDblocks_GRCh38/data/pyrho_EUR_LD_blocks.bed"

# No X chromsome!
CHROMS, = glob_wildcards("../input/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{ID}.recalibrated_variants.vcf.gz")

ped_file = "../input/20130606_g1k_3202_samples_ped_population.txt"

# POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
POPULATIONS = ["EUR"] # only using EUR LD blocks for now

output_bcf_dir = "../output/bcf"

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
        expand(os.path.join(output_bcf_dir, "{BLOCK}", "variant_list_{POP}.tsv"), BLOCK = LD_blocks.block.values, POP = POPULATIONS)

rule bcf:
    input:
        vcf = "../input/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{ID}.recalibrated_variants.vcf.gz"
    output:
        bcf = os.path.join(output_bcf_dir, "{ID}.recalibrated_variants.bcf")
    threads: 2 
    resources:
        mem = "2G",
        partition = "shared",
        time = "15:30:00"
    shell:
        """
        module load bcftools
        bcftools view -m2 -M2 -v snps -f PASS -i "MAC >= 5 && F_MISSING < 0.15" -O b {input.vcf} > {output.bcf}
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
        module load bcftools
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
        module load bcftools
        bcftools view -r {params.CHROM}:{params.START}-{params.END} -S {input.samples} -O u {params.BCF} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b > {output.bcf}
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
        module load bcftools
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
        module load bcftools
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input.bcf} > {output.variant_file}
        # plink2 does not reorder alleles, unlike plink1.9
        plink2 --bcf {input.bcf} --make-bed --out {params.prefix}
        """

