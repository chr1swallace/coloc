# Prior sensitivity for coloc

Shows how prior and posterior per-hypothesis probabilities change as a
function of p12

## Usage

``` r
sensitivity(
  obj,
  rule = "",
  dataset1 = NULL,
  dataset2 = NULL,
  npoints = 100,
  doplot = TRUE,
  plot.manhattans = TRUE,
  preserve.par = FALSE,
  row = 1
)
```

## Arguments

- obj:

  output of coloc.detail or coloc.process

- rule:

  a decision rule. This states what values of posterior probabilities
  "pass" some threshold. This is a string which will be parsed and
  evaluated, better explained by examples. "H4 \> 0.5" says post prob of
  H4 \> 0.5 is a pass. "H4 \> 0.9 & H4/H3 \> 3" says post prob of H4
  must be \> 0.9 AND it must be at least 3 times the post prob of H3."

- dataset1:

  optional the dataset1 used to run SuSiE. This will be used to make a
  Manhattan plot if plot.manhattans=TRUE.

- dataset2:

  optional the dataset2 used to run SuSiE. This will be used to make a
  Manhattan plot if plot.manhattans=TRUE.

- npoints:

  the number of points over which to evaluate the prior values for p12,
  equally spaced on a log scale between p1\*p2 and min(p1,p2) - these
  are logical limits on p12, but not scientifically sensible values.

- doplot:

  draw the plot. set to FALSE if you want to just evaluate the prior and
  posterior matrices and work with them yourself

- plot.manhattans:

  if TRUE, show Manhattans of input data

- preserve.par:

  if TRUE, do not change par() of current graphics device - this is to
  allow sensitivity plots to be incoporated into a larger set of plots,
  or to be plot one per page on a pdf, for example

- row:

  when coloc.signals() has been used and multiple rows are returned in
  the coloc summary, which row to plot

## Value

list of 3: prior matrix, posterior matrix, and a pass/fail indicator
(returned invisibly)

## Details

Function is called mainly for plotting side effect. It draws two plots,
showing how prior and posterior probabilities of each coloc hypothesis
change with changing p12. A decision rule sets the values of the
posterior probabilities considered acceptable, and is used to shade in
green the region of the plot for which the p12 prior would give and
acceptable result. The user is encouraged to consider carefully whether
some prior values shown within the green shaded region are sensible
before accepting the hypothesis. If no shading is shown, then no priors
give rise to an accepted result.

## Author

Chris Wallace
