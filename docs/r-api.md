# R API

The full R reference is published under `/r/` via pkgdown.

Reference index: `/r/reference/`

Build it locally after the MkDocs site:

```r
roxygen2::roxygenise()
pkgdown::build_site(new_process = FALSE, override = list(destination = "site/r"))
```
