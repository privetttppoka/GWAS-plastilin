rule all:
    input:
        "results/gwas/gwas_GLM.assoc.linear",
        "results/plots/manhattan.png",
        "results/plots/qqplot.png",
        "results/tree/phylo_tree.pdf",
        "results/pca/pca.eigenvec",
        "results/plots/pca.png"

# -----------------------
# 1. Filter VCF
# -----------------------

rule filter_vcf:
    input: "data/test_data.vcf"
    output: "results/clean_data.vcf"
    shell:
        """
        bcftools view -e 'ALT=REF' {input} -Oz -o {output}
        """

# -----------------------
# 2. PLINK: VCF → BED
# -----------------------

rule vcf_to_plink:
    input: "results/clean_data.vcf"
    output:
        bed="results/data0.bed",
        bim="results/data0.bim",
        fam="results/data0.fam"
    shell:
        """
        plink --vcf {input} --double-id --make-bed --out results/data0
        """

# -----------------------
# 3. QC (MAF + GENO)
# -----------------------

rule qc:
    input:
        bed="results/data0.bed",
        bim="results/data0.bim",
        fam="results/data0.fam"
    output:
        bed="results/data_qc.bed",
        bim="results/data_qc.bim",
        fam="results/data_qc.fam"
    shell:
        """
        plink --bfile results/data0 \
              --maf 0.05 --geno 0.1 \
              --make-bed --out results/data_qc
        """

# -----------------------
# 4. LD pruning
# -----------------------

rule prune:
    input: "results/data_qc.bed"
    output:
        prunin="results/pruned.prune.in",
        prunout="results/pruned.prune.out"
    shell:
        """
        plink --bfile results/data_qc \
              --indep-pairwise 200 50 0.2 \
              --out results/pruned
        """

rule prune_extract:
    input:
        bed="results/data_qc.bed",
        prunin="results/pruned.prune.in"
    output:
        bed="results/data_pruned.bed",
        bim="results/data_pruned.bim",
        fam="results/data_pruned.fam"
    shell:
        """
        plink --bfile results/data_qc \
              --extract results/pruned.prune.in \
              --make-bed --out results/data_pruned
        """

# -----------------------
# 5. PCA
# -----------------------

rule pca:
    input: "results/data_pruned.bed"
    output:
        eigenvec="results/pca/pca.eigenvec",
        eigenval="results/pca/pca.eigenval"
    shell:
        """
        mkdir -p results/pca
        plink --bfile results/data_pruned \
              --pca 10 --out results/pca/pca
        """

# -----------------------
# PCA plot
# -----------------------

rule pca_plot:
    input: "results/pca/pca.eigenvec"
    output: "results/plots/pca.png"
    shell:
        """
        mkdir -p results/plots
        python scripts/plot_pca.py {input} {output}
        """

# -----------------------
# 6. Build Q matrix
# -----------------------

rule q_matrix:
    input: "results/pca/pca.eigenvec"
    output: "results/Q_matrix.txt"
    shell:
        """
        awk 'BEGIN{{print "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10"}} \
             {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' \
             {input} > {output}
        """

# -----------------------
# 7. Prepare phenotype file
# -----------------------

rule pheno:
    input: "data/test_data.tsv"
    output: "results/phenotypes_plink.tsv"
    shell:
        """
        mkdir -p results
        echo -e "FID\tIID\tyield" > {output}
        tail -n +2 {input} | awk '{{print $1"\\t"$1"\\t"$2}}' >> {output}
        """

# -----------------------
# 8. GWAS (PLINK linear)
# -----------------------

rule gwas:
    input:
        bed="results/data_qc.bed",
        pheno="results/phenotypes_plink.tsv",
        covar="results/Q_matrix.txt"
    output: "results/gwas/gwas_GLM.assoc.linear"
    shell:
        """
        mkdir -p results/gwas
        plink --bfile results/data_qc \
              --pheno results/phenotypes_plink.tsv --pheno-name yield \
              --covar results/Q_matrix.txt \
              --linear hide-covar \
              --allow-no-sex \
              --out results/gwas/gwas_GLM
        """

# -----------------------
# 9. Manhattan plot (Python)
# -----------------------

rule manhattan:
    input: "results/gwas/gwas_GLM.assoc.linear"
    output: "results/plots/manhattan.png"
    shell:
        """
        mkdir -p results/plots
        python scripts/plot_manhattan.py {input} {output}
        """

# -----------------------
# 10. QQ plot (R)
# -----------------------

rule qqplot:
    input: "results/gwas/gwas_GLM.assoc.linear"
    output: "results/plots/qqplot.png"
    shell:
        """
        mkdir -p results/plots
        Rscript scripts/plot_qq.R {input} {output}
        """

# -----------------------
# 11. Distance matrix
# -----------------------

rule dist:
    input: "results/data_pruned.bed"
    output:
        matrix="results/dist/dist.mdist",
        ids="results/dist/dist.mdist.id"
    shell:
        """
        mkdir -p results/dist
        plink --bfile results/data_pruned \
              --distance square 1-ibs \
              --out results/dist/dist
        """

# -----------------------
# 12. Phylogenetic tree
# -----------------------

rule tree:
    input:
        matrix="results/dist/dist.mdist",
        ids="results/dist/dist.mdist.id"
    output: "results/tree/phylo_tree.pdf"
    shell:
        """
        mkdir -p results/tree
        Rscript scripts/plot_tree.R {input.matrix} {input.ids} {output}
        """
