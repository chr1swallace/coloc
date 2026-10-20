# logdiff

Internal function, logdiff

## Usage

``` r
logdiff(x, y)
```

## Arguments

- x:

  numeric

- y:

  numeric

## Value

max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))

## Details

This function calculates the log of the difference of the exponentiated
logs taking out the max, i.e. insuring that the difference is not
negative

## Author

Chris Wallace
