[![DOI](https://zenodo.org/badge/235801923.svg)](https://zenodo.org/badge/latestdoi/235801923)

# ErrViewLib

Package for the statistical analysis and plotting of error sets from
computational chemistry methods.

Now augmented with UQ validation methods.

## Install

```r
install.packages("devtools")
devtools::install_github("ppernot/ErrViewLib")
```

## Usage

`ErrViewLib` is used by the [ErrView](https://github.com/ppernot/ErrView) 
shiny code, but can be used in base R scripts. For instance:

```r
library(ErrViewLib)

# Get data
data = read.csv(
  file = 'myData.csv',
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract error sets and misc. info.
systems <- data[, 1]                   # System names
rownames(data) = systems
Ref  = data[, 2]                       # Reference values
Data = data[, -c(1, 2), drop = FALSE]  # Data
methList = colnames(Data)              # Methods names
Errors   = Ref - Data                  # Errors

# Estimate statistics with bootstrap uncertainties 
# (MUE, Q95 and GMCF)
bs = ErrViewLib::estBS1(
  Errors,
  props = c('mue','q95hd','gimc'),
  do.sip = FALSE
)

# Generate and print simple results table
df = ErrViewLib::genTabStat(bs,comp = FALSE)
print(knitr::kable(df))

# Plot Lorenz curves
ErrViewLib::plotLorenz(
  Errors,
  gPars = ErrViewLib::setgPars( # Redefine plot params
    type = 'plot',
    gPars = list(
      lwd = 1,
      cols = rainbow(length(methList))
    )
  )
)
```
