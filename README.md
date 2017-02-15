# MIRMMR (murmur)
Microsatellite Instability Regression using Methylation and Mutations in R
---
# TODO
1. Complete README.md
  1. Usage statements
  2. Modules
    1. Penalized
    2. Stepwise
    3. Univariate
    4. Predict
    5. Compare
2. Clean up code / small improvements
3. Automate output of paper stats
4. Improve optimal lambda search (smooth lambda/error curve?)  
---
# Usage
```
Rscript murmur.R -m module -f data.frame -i msi.status -c first.data.column -o output.prefix -d output.directory [options]
```
---
# Modules
### Penalized
```
Rscript murmur.R -m penalized -f data.frame -i msi.status -c first.data.column -o output.prefix -d output.directory [options]
```
### Stepwise
```
Rscript murmur.R -m penalized -f data.frame -i msi.status -c first.data.column -o output.prefix -d output.directory [options]
```
### Univariate
### Predict
### Compare
