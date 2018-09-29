# coloc 3.2 
*new function: colocqq*
- added colocqq for quantitative traits measured on the same individuals
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

