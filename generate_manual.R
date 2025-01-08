devtools::document()
devtools::check(args="--as-cran")
devtools::build_manual(path=".")
