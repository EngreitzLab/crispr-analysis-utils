# R API

The full R reference is published under `/r/` via pkgdown.

Open the R docs:

- [R site home](/crispr-analysis-utils/r/)
- [R reference index](/crispr-analysis-utils/r/reference/)

Build it locally after the MkDocs site:

```r
roxygen2::roxygenise()
pkgdown::build_site(new_process = FALSE, override = list(destination = "site/r"))
```
