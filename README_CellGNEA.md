# CellGNEA

This repository contains the toy dataset and R code used in the analyses for the paper:

**"Inferring Cell Line-Specific Gene Network Enrichment Patterns Along Continuous Phenotypes"**

---

## 🔬 Overview

CellGNEA is a framework designed to identify **pathway-level enrichment patterns** from **cell line–specific gene regulatory networks** along a continuous phenotype (e.g., drug response, molecular gradient).

The method integrates:
- Cell-specific network inference
- Network topology features
- Gene-level scoring
- Enrichment analysis with permutation testing

---

## 📊 Graphical Workflow

```
Input Data
 ├── Gene Expression Matrix (EXP)
 ├── Phenotype / Modulator (Modulator)
 └── Pathway Gene Set (PW_genes)

        ↓

[Step 1] Network Inference
→ Construct cell-specific gene regulatory networks (EdgeW)

        ↓

[Step 2] Feature Extraction
→ Clustering Coefficient (C*)
→ PageRank
→ Regulatory Effect (RE)

        ↓

[Step 3] Gene Scoring
→ Integrate network features into gene-level scores

        ↓

[Step 4] Association Analysis
→ Correlate gene scores with phenotype

        ↓

[Step 5] Enrichment Analysis
→ Compute enrichment score (ES)
→ Permutation test → p-value

        ↓

Output
 ├── Enrichment Score (ES)
 └── Statistical Significance (p-value)
```

---

## 📁 Repository Structure

```
ToyDATA_CellGNEA/
 ├── EXP.csv
 ├── PathwayGENES.csv
 ├── Modulator.csv
 ├── BETA_Sample1.csv
 ├── BETA_Sample2.csv
 └── ...
```

---

## ⚙️ Requirements

```r
library(data.table)
library(igraph)
```

---

## 🚀 Step-by-Step Usage

### Step 1. Load Data

```r
EXP <- read.table("ToyDATA_CellGNEA/EXP.csv", sep=",")
PW_genes <- read.table("ToyDATA_CellGNEA/PathwayGENES.csv", sep=",")
Modulator <- read.table("ToyDATA_CellGNEA/Modulator.csv", sep=",")
```

---

### Step 2. Network-Based Gene Scoring

For each cell line:
- Load network (BETA_Sample i)
- Compute clustering coefficient, PageRank, and regulatory effect
- Combine into gene scores

---

### Step 3. Correlation with Phenotype

```r
CORR <- apply(SCORE, 2, function(x) cor(x, Modulator))
```

---

### Step 4. Enrichment Analysis

```r
NOpm <- 1001
```

---

### Step 5. Visualization

```r
plot(erSCORE_ORG[,6], type="l", main="Enrichment Score")
abline(h=0, lty=2)
abline(v=which.max(abs(erSCORE_ORG[,6])), col="red", lwd=2)
```

---

## ✅ Expected Output

- Enrichment Score (ES)
- Normalized ES
- p-value
- Enrichment plot

---

## 📌 Notes

- Number of permutations = 1000
- First permutation = observed pathway
- Others = random gene sets

---

## 📬 Contact

Please open an issue for questions or feedback.
