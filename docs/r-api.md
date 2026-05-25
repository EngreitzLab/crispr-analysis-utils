# R API

The full R reference is published under `/r/` via pkgdown.

Build it locally after the MkDocs site:

```r
roxygen2::roxygenise()
pkgdown::build_site(new_process = FALSE, override = list(destination = "site/r"))
```
