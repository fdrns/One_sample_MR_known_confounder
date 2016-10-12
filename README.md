# One sample MR with known confounder

This script allows the simulation of a one sample MR with a known confounder.

Two correlated measures are generated and 20 SNPs. 10 of the SNPs have an effect in both traits violating the MR assumptions. Using a basic implementation of a genetic algorithm we use different criteria to find the optimal subset of SNPs to be used in MR of the first measure.

The problem can be approached as a feature selection problem.  We can either use a suitable criterion to identify the SNPs that need to be removed (example R square for the confounder) or a combination (example difference of R squares for the two measures) criterion that can provide a parsimonious solution while removing pleiotropic SNPs.

When individual SNPs are used, the difference between the R squares provide a good selection while when we use a gene score the ratio of the F-statistics can also provide an acceptable, but less efficient, solutions.  The use of the F-statistic in individual SNPs performs less well than when this is applied in the unweighted score of SNPs.
