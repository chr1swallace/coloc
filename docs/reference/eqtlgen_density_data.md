# eQTLGen estimated distance density

eQTLGen estimated distance density

## Usage

``` r
data(eqtlgen_density_data)
```

## Format

A data.frame with two columns, "rel_dist" the relative distance to the
TSS and "dens_value" the estimated value of the density. This density is
esimated from the eQTLGen data. For details of how the density was
estiamted see Pullin and Wallace (2024+).

## Examples

``` r
data(eqtlgen_density_data)
names(eqtlgen_density_data)
#> [1] "rel_dist"   "dens_value"
```
