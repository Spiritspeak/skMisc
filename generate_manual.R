devtools::document()
devtools::build_manual(path=".")
devtools::check(args="--as-cran")
