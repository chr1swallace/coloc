# plot a coloc dataset

Plot a coloc structured dataset

## Usage

``` r
plot_dataset(
  d,
  susie_obj = NULL,
  highlight_list = NULL,
  alty = NULL,
  ylab = "-log10(p)",
  show_legend = TRUE,
  color = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "gold1", "skyblue2",
    "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
    "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
    "yellow3", "darkorange4", "brown"),
  ...
)
```

## Arguments

- d:

  a coloc dataset

- susie_obj:

  optional, the output of a call to runsusie()

- highlight_list:

  optional, a list of character vectors. any snp in the character vector
  will be highlighted, using a different colour for each list.

- alty:

  default is to plot a standard manhattan. If you wish to plot a
  different y value, pass it here. You may also want to change ylab to
  describe what you are plotting.

- ylab:

  label for y axis, default is -log10(p) and assumes you are plotting a
  manhattan

- show_legend:

  optional, show the legend or not. default is TRUE

- color:

  optional, specify the colours to use for each credible set when
  susie_obj is supplied. Default is shamelessly copied from
  susieR::susie_plot() so that colours will match

- ...:

  other arguments passed to the base graphics plot() function

## Author

Chris Wallace
