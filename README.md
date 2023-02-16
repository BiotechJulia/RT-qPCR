# RT-qPCR
Simple RT-qPCR data analysis and visualization.

This script lets you calculate delta CT, assuming perfect primer efficiency. If your primer efficiency is not between 90-100%, calculate an efficiency score by dividing the percent score by 100 and adding 1. For example, if the primer efficiency was 85%, the efficiency score would be: 85/100+1=1.85. In that case, replace the value "2" in lines 36 and 38 with the real efficiency scores for the reference and target genes.

Next, an analysis of variance (ANOVA) is performed with Tukey's HSD (honestly significant difference) post-hoc test. Values belonging to the same homogeneous group are denoted with the same letter (p<0.05).

Finally, mean and standard deviation are calculated for relative expression values. These values are plotted using bar plots with the letters denoting statistical significance above.

The sample data includes average CT values obtained from three technical replicates; three biological replicates were analyzed. The expression of three genes (two target genes and one reference gene) for two lines (wild type and mutant) was analyzed.

![Test RT-qPCR2](https://user-images.githubusercontent.com/125211875/219329001-67fc810f-81cb-410f-ae5e-5dba7b77a3b1.png)
