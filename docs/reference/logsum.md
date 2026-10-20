# logsum

Internal function, logsum

## Usage

``` r
logsum(x)
```

## Arguments

- x:

  numeric vector

## Value

max(x) + log(sum(exp(x - max(x))))

## Details

This function calculates the log of the sum of the exponentiated logs
taking out the max, i.e. insuring that the sum is not Inf

## Author

Claudia Giambartolomei
