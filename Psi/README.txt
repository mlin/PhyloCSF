This directory contains files for calculating PhyloCSF-Psi from the PhyloCSF score and region length, as described in Section 4.2 of Lin et al. 2011 (https://academic.oup.com/bioinformatics/article/27/13/i275/178183). PhyloCSF-Psi is a length-dependent score that uses the empirical distributions of PhyloCSF scores in coding and non-coding regions of various lengths to increase accuracy, set thresholds, and calculate a p-value.

Each of the coefficients files contains six numbers, which are the values of the six constants needed to evaluate the formual for Psi in Section 4.2, namely:
	MUcoding Acoding Bcoding MUnoncoding Anoncoding Bnoncoding.
The name of each coefficients file is of the form 
    PsiCoefs.{Alignment Set}.{PhyloCSF Parameters}.{PhyloCSF strategy}.txt, 
where:
- {Alignment Set} defines the whole-genome multispecies alignment from which the local alignment supplied to PhyloCSF was extracted. The Psi score will not be meaningful unless the local alignment was extracted from this whole-genome alignment and the only species with aligned sequence in the whole-genome alignment that were removed are those not included in the Alignment Set. Alignment Sets are defined here: https://data.broadinstitute.org/compbio1/cav.php?Alnsets.
- {PhyloCSF Parameters} is the name of the parameter set that was supplied as the first argument to PhyloCSF (for example 12flies).
- {PhyloCSF strategy} is the strategy supplied in the --strategy argument to PhyloCSF, which is either mle (default), fixed, or omega.

More about the motivation and use of PhyloCSF-Psi can be found on the PhyloCSF wiki here: https://github.com/mlin/PhyloCSF/wiki.

