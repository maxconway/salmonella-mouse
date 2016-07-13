Firstly, to try to better highlight important reactions, we modified the filtering procedure. Now, after the first optimization round, we do a second minimization round that attempts to reduce all fluxes as much as possible while maintaining the same total result. This removes two sources of noise: loops that don't affect total production of metabolites, and variation in fluxes that doesn't affect the result but results from the internal workings of the optimization library.
The flux suppression resulting from this also allowed us to increase the effect of gene expression on flux bounds.

With four contrasts are of particular interest: Input vs Group1, Input vs Group2, Group1 vs Group2 and Group2 vs Group3. Originally I tried to group discussion into these contrasts, but since many reactions display differences across multiple contrasts, a reaction - foremost view was more effective.

Group 3 unfortunately has not produced a viable model under the techniques that we've applied. This is interesting, and quite possibly biologically meaningful, because on many metrics it appears between groups 1 and 2. This implies that rather than being due to just a lack of data, the model does not cover some pathway that is required in group 3.

Eight reactions display strong variation between Input, Group 1 and Group 2:
- The glycolate oxidases, GLYCTO2 and GLYCTO4
- The 2 3 diketo 5 methylthio 1 phosphopentane degradation reactions, DKMPPD and DKMPPD3
- sulfite reductases SULRi and SO3R
- Acetyl CoA ACP transacylase, ACOATA
- beta ketoacyl ACP synthases, KAS14 and KAS15
- ATP maintenance allowance

In terms of glycolate oxidases, we find that input uses GLYTO4, whereas the two groups use GLYTO2. This means that Input oxidises Glycolate using Demethylmenaquinone 8, whereas Groups 1 and two use quinone 8. In order of rate, Group2 > Input > Group1, but they're all within 10% of each other.

In 2 3 diketo 5 methylthio 1 phosphopentane degradation, group 2 and input use Nicotinamide adenine dinucleotide to break down dkmpp, while group 1 does not. Group 1's reaction is about 10% slower, so this may be a restriction, rather than an advantage.

In sulfite redutases, we once again see group 1 behave differently. It uses nadh to reduce so3, but Input and group 3 use nadph. Group 1 is also slower, which indicates it may be restricted. This could indicate that it requires nad more than nadp, or that it has excess nadh vs nadph. This ties in with the previous reaction to suggest that group 1 is restricted in it's ability to convert nadh to nad, or somehow else faces an excess of nadh and a dearth of nad.

In the pseudo reaction ATPM, which converts ATP to ADP without any work getting done, as a proxy for mischelanous processes that require ATP, Group 1 operates at about half the speed of Group 2 and Input (0.376 vs 0.791 and 0.774 respectively). This may be linked to a general shortage of phosphate which interpheres with the previous reactions as well.

The interaction between ACOATA, KAS 14 and KAS 15 is best described together. Input uses ACOATA and KAS14 together in sequence to perform the same reaction as KAS15. Because KAS14 and ACOATA are at the same rate, there is no net difference between these strategies, so this change doesn't have any knock on effects.