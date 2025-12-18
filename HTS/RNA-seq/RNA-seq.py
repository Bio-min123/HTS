import subprocess
import yaml
from pathlib import Path
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# -----------------------
# Logging setup
# -----------------------
logging.basicConfig(
    filename="logs/pipeline.log",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s"
)

# -----------------------
# Load config
# -----------------------
with open("config.yaml") as f:
    config = yaml.safe_load(f)

FASTQ_DIR = Path(config["paths"]["fastq"])
TRIM_DIR = Path(config["paths"]["trimmed"])
QC_DIR = Path(config["paths"]["qc"])
STAR_DIR = Path(config["paths"]["star"])
COUNTS_DIR = Path(config["paths"]["counts"])
EDGER_DIR = Path(config["paths"]["edgeR"])
PLOTS_DIR = Path(config["paths"]["plots"])
GENOME_DIR = config["reference"]["genome_dir"]
GTF_FILE = config["reference"]["gtf"]

for d in [TRIM_DIR, QC_DIR, STAR_DIR, COUNTS_DIR, EDGER_DIR, PLOTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# -----------------------
# Step 1: FastQC + MultiQC
# -----------------------
def run_fastqc():
    logging.info("Starting FastQC...")
    for fq in FASTQ_DIR.glob("*.fastq.gz"):
        subprocess.run(["fastqc", fq, "--outdir", QC_DIR], check=True)
    subprocess.run(["multiqc", QC_DIR, "-o", QC_DIR], check=True)
    logging.info("FastQC + MultiQC finished.")


# -----------------------
# Step 2: Trim Galore
# -----------------------
def trim_reads():
    logging.info("Starting Trim Galore...")
    samples = set(f.name.split("_R")[0] for f in FASTQ_DIR.glob("*_R1.fastq.gz"))
    for sample in samples:
        r1 = FASTQ_DIR / f"{sample}_R1.fastq.gz"
        r2 = FASTQ_DIR / f"{sample}_R2.fastq.gz"
        subprocess.run([
            "trim_galore", "--paired",
            "--output_dir", TRIM_DIR,
            r1, r2
        ], check=True)
    logging.info("Trimming finished.")


# -----------------------
# Step 3: STAR alignment (BAM)
# -----------------------
def run_star():
    logging.info("Starting STAR alignment...")
    samples = set(f.name.split("_R")[0] for f in FASTQ_DIR.glob("*_R1.fastq.gz"))
    for sample in samples:
        r1 = TRIM_DIR / f"{sample}_R1_val_1.fq.gz"
        r2 = TRIM_DIR / f"{sample}_R2_val_2.fq.gz"
        out_prefix = STAR_DIR / sample
        out_prefix.mkdir(exist_ok=True)
        subprocess.run([
            "STAR",
            "--runThreadN", str(config["threads"]),
            "--genomeDir", GENOME_DIR,
            "--readFilesIn", str(r1), str(r2),
            "--readFilesCommand", "zcat",
            "--outFileNamePrefix", str(out_prefix) + "/",
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--quantMode", "GeneCounts"
        ], check=True)
    logging.info("STAR alignment finished.")


# -----------------------
# Step 4: Generate counts using featureCounts
# -----------------------
def run_featurecounts():
    logging.info("Running featureCounts on BAM files...")
    bam_files = [str(f / "Aligned.sortedByCoord.out.bam") for f in STAR_DIR.iterdir() if f.is_dir()]
    out_file = COUNTS_DIR / "gene_counts.txt"
    subprocess.run([
                       "featureCounts",
                       "-T", str(config["threads"]),
                       "-a", GTF_FILE,
                       "-o", str(out_file),
                       "-p",  # paired-end
                       "-B",  # only count properly paired reads
                       "-C",  # skip chimeric fragments
                   ] + bam_files, check=True)
    logging.info(f"Counts saved to {out_file}")
    return out_file


# -----------------------
# Step 5: Differential expression using edgeR (R script)
# -----------------------
def run_edger(counts_file):
    logging.info("Starting edgeR differential expression...")
    subprocess.run([
        "Rscript", "scripts/run_edger.R",
        str(counts_file),
        "data/metadata.csv",
        str(EDGER_DIR / "edgeR_results.csv")
    ], check=True)
    logging.info("edgeR finished.")


# -----------------------
# Step 6: Visualization (PCA, heatmap, volcano)
# -----------------------
def plot_results():
    logging.info("Plotting results...")
    res = pd.read_csv(EDGER_DIR / "edgeR_results.csv", index_col=0)
    counts = pd.read_csv(COUNTS_DIR / "gene_counts.txt", sep="\t", comment="#", index_col=0)

    # Keep only counts columns (skip first columns STAR adds)
    sample_cols = res.columns.intersection(counts.columns)
    counts_matrix = counts[sample_cols].T

    # PCA
    pca = PCA(n_components=2)
    X = pca.fit_transform(counts_matrix)
    meta = pd.read_csv("data/metadata.csv", index_col=0)
    df_pca = pd.DataFrame(X, columns=["PC1", "PC2"], index=counts_matrix.index)
    df_pca = df_pca.merge(meta, left_index=True, right_index=True)

    plt.figure()
    sns.scatterplot(data=df_pca, x="PC1", y="PC2", hue="condition", s=100)
    plt.title("PCA of samples")
    plt.savefig(PLOTS_DIR / "pca.png", bbox_inches="tight")

    # Heatmap of top 50 DE genes
    top_genes = res.sort_values("FDR").head(50).index
    sns.clustermap(counts.loc[top_genes, sample_cols], z_score=0, cmap="vlag")
    plt.savefig(PLOTS_DIR / "heatmap_top50.png", bbox_inches="tight")

    # Volcano plot
    plt.figure()
    sns.scatterplot(
        x=res["logFC"],
        y=-np.log10(res["FDR"]),
        hue=(res["FDR"] < 0.05) & (abs(res["logFC"]) > 1),
        palette={True: "red", False: "grey"}
    )
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 FDR")
    plt.title("Volcano plot")
    plt.savefig(PLOTS_DIR / "volcano.png", bbox_inches="tight")

    logging.info("Plots saved.")


# -----------------------
# Main
# -----------------------
if __name__ == "__main__":
    run_fastqc()
    trim_reads()
    run_star()
    counts_file = run_featurecounts()
    run_edger(counts_file)
    plot_results()
