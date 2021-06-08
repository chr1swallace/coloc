# coloc 5.1.0
* release to CRAN of coloc.susie
* remove nref now not recommended

# coloc 5.0.0.9002
* in runsusie, allow nref to be overridden if user passes z_ld_weight

# coloc 5.0.0.9001
* deprecate coloc.signals for multiple causal variants
* introduce coloc.susie for multiple causal variants
* new vignette 02_data giving more details on structuring your data properly

# coloc 4.0-5
* update error checking

# coloc 4.0-3
* BUGFIX
- missing sdY for type="quant" would error in finemap.signals or coloc.signals, now it is estimated if missing.
- warning for factors used as snp names added (issue #29)

# coloc 4.0-2
* BUGFIX
- snps in the second dataset might not have been masked as intended. This would result in effect using the "single" option for trait 2 when masking was expected.

# coloc 4.0-0
*new functions: coloc.signals, finemap.abf*
- analogues of coloc.abf and finemap.abf that allow for multiple causal variants.  See vignette on conditioning/masking
*new function: sensitivity*
- post-hoc, determine the sensitivity of coloc results to changes in the prior.  See vignette on sensitivity

2019-06-25

# coloc 3.2-1
* BUGFIX: finemap.abf()
- in low power situations, the posterior for H0 was previously too low.  This will only affect datasets where the minimum p value was > 1e-7 - ie where the posterior for H0 would be expected to be much above 0.

# coloc 3.2 
- fix bug in process.datasets which suggested MAF was needed for cc data when beta/varbeta also present
- added pkgdown
	
2018-09-28

# coloc 3.1
*new function: finemap.abf*
- improved clarity of error messages in coloc.abf() sub functions
- added finemap.abf to fine map a single trait
- fixed warnings caused by CRAN disliking the BioConductor devel branch
	
2018-02-22
	
# coloc 2.4 
- added stratification to proportional testing approaches
	
2017-03-21 

# coloc 2.3-2 
- tidied code relating mainly to proportional
- colocalisation testing methods, making more methods confirm to S4.
- pcs.prepare now imputates missing genotypes by default
2013-05-22 

# coloc 2.3 
## coloc BUGFIX
- Introduced a function to estimate trait variance from supplied coefficients and standard errors.  This is used within the approach implemented in coloc.abf(), and replaces the earlier version which implicity assumed that var(Y)=1 for quantitative traits, which could lead to incorrect inference when var(Y) was far from 1.

2013-09-25  

# coloc 2.2
- Merged coloc.abf and coloc.abf.imputed(), so that datasets for wheich beta, var(beta) are available can be matched to datasets with only p values and maf.2 This means the arguments to coloc.abf() have been changed!  Please
check ?coloc.abf for the new function.

2013-19-06  

# coloc 2.1
- Bug fix for coloc.abf() function, which used p12 instead of log(p12) to calculate L4.
- New function coloc.abf.imputed() to make better use of fuller information on imputed data.

2013-03-06  

# coloc 2.0
- New function, coloc.abf(), to implement the colocalisation approach described by Giambartolomei et. al.
- Changes in the coloc.test() and coloc.bma() functions to make them consistent with regards arguments and output.

2013-22-05  

# coloc 1.11
- added principal components functions pcs.prepare(), pcs.model().
- Restructed the coloc objects to separate Bayesian and non-Bayesian inference.

2012-12-10  

# coloc 1.10
- added Credible Interval calculation to coloc.test().

2012-11-28  

# coloc 1.09
- updated to return u and Var(u) in addition to chisq statistic.

2012-07-12  

# coloc 1.08
- fixed error in documentation, added MASS to Depends.

2012-04-12  

# coloc 1.06
- some tidying.

2012-01-05  

# coloc 1.05
- moved to S4 methods.

2012-01-03  

# coloc 1.04
- made the means of generating plots more flexible.

2011-12-27  

