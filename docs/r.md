# R API

The published GitHub Pages site replaces this placeholder with a pkgdown site
for the R package.

Build it locally after the MkDocs site:

```r
roxygen2::roxygenise()
pkgdown::build_site(new_process = FALSE, override = list(destination = "site/r"))
```
