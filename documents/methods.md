Detailed methods
===============

Preparing RNA values
--------------------
The sequence count data set consists of 11 measurements for each of 6579 genes: three technical replicates for each of the three in vivo conditions, and two technical replicates for the in vitro input condition. We first find the logarithm of the sequence count plus one. We use a logarithm because the sequence counts are approximately lognormally distributed, and we use the addition of one because some of the counts are zero, particularly in the in vitro experiments where the values are generally much smaller.
We then normalize the replicates against their mean, to remove technical and experimental differences, and then find the mean of the two or three replicates for each gene within each experimental group. We normalize each group against the in vitro input as a control. Finally we adjust the standard deviation to 0.1. The resulting values are appoximately normally distributed around one, but with a spike at one. 

Combining gene expressions for the metabolic model
---------------------------------------------------
We use a metabolic network with Gene-Protein-Reaction mappings from the literature. To combine this with the RNA data, we map each gene in the RNA dataset to the equivalent gene in the metabolic model. Where multiple genes correspond to a particular reaction, we combine them via a continuous extrapolation of the boolean function encoded in the Gene-Protein_Reaction mappings. Specifically, if *both* of the genes or gene sets are required (i.e. an `AND` relation), we take the minimum activation of the two; if *one* of the genes or gene sets are required (i.e. an `OR` relation), we take the maximum activation of the two. Where we do not have information about a gene, we assume that is abundant by default, leaving the other gene in the expression to dictate the output level. 

Finding fluxes based on gene expression
---------------------------------------
We first use standard flux balance analysis to find the base flux levels without any gene expression restrictions. We then calculate a new flux "target" for each group by multiplying the base flux level by the activation, and set new flux bounds as 10% of this target on either side. Based on these new flux bounds, we conduct our final round of flux balance analysis to find the actual predicted fluxes.

