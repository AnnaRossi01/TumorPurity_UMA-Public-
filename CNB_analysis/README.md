# Copy Number Burden Workflow

## Objective

This analysis evaluates the genomic instability of each tumor sample through Copy Number Burden (CNB), with the goal of understanding and justifying the behavior of tumor purity estimation algorithms.

## Workflow 

### Phase 1: Orthogonal Algorithm

The .cns files obtained from CNVkit are converted to .seg format. This step allows obtaining an independent representation of the actual genomic state of the sample, without introducing bias from the tools being evaluated.

### Phase 2: Flat vs Altered Classification

The Python script processes the .seg files and calculates the Copy Number Burden, defined as the percentage of altered genome within the target panel. If the burden is below a predefined threshold, the sample is classified as "Flat"; otherwise, it is classified as "Altered".

### Phase 3: Stratification

In the main database containing tumor samples, a column is added with the "Flat" or "Altered" classification for each ID.

### Phase 4: Critical Analysis

The database is divided into the two subgroups and the main visualizations are re-executed separately for the Flat cohort and the Altered cohort, including correlation matrices, MAE, and RMSE.

---

## Thesis Notes: Copy Number Burden (CNB) Analysis

### 1. Logical Basis

CNB analysis was introduced to investigate the collapse of performance in tumor purity estimation algorithms (ASCAT.sc, IchorCNA, PureCN) in the Multiple Myeloma patient cohort.

Copy Number-based tools do not directly quantify tumor cells but estimate their fraction by analyzing genomic instability. The central hypothesis is that their failure is not due to a limitation of the methodology itself, but rather to biological and technical constraints specific to this context.

- Biological limit: MM is heterogeneous and not all patients show a strongly unstable genome.
- Technical limit: targeted sequencing covers a very small portion of the genome (0.46 MB). If alterations do not fall within the targeted regions, the panel appears virtually flat despite the actual presence of genomic instability.

### 2. Methodology and Orthogonal Tool

To avoid internal validation bias, the segmentation files produced by the tools under examination were not used. CNVkit was chosen as an orthogonal and independent approach for genome segmentation.

### 3. Parameters Evaluated

The Python script extracted the following parameters from CNVkit .seg files:

- seg.mean (Log2 Ratio): measures the deviation of copy number from the diploid baseline.
- num.mark: indicates the number of markers or probes included in each segment.
- Copy Number Burden (CNB): fraction or percentage of the target panel affected by copy number alterations.

Formula:

$$
CNB = \frac{\text{sum of num.mark from altered segments}}{\text{total sum of num.mark in the panel}}
$$

### 4. Cut-offs and Thresholds Used

Two analytical thresholds were defined:

- Segment alteration threshold: Log2 Ratio values greater than 0.2 or less than -0.2.
- Genome classification threshold: CNB equal to 15% (0.15).

Classification:

- Flat: CNB < 0.15
- Altered: CNB >= 0.15

### 5. Goal of Stratification

The division into two groups allows recalculating validation metrics such as Spearman correlation, MAE, and RMSE separately for Flat and Altered samples.

- If performance improves in the Altered group, the utility of the algorithms is confirmed, but limited to patients with high genomic burden.
- If performance remains unacceptable in both groups, this suggests a structural limitation of targeted sequencing at 0.46 MB that cannot be overcome.

SAMPLES -> 156 (98 altered, 58 flat)
