library(lintr)

mylinters <- linters_with_tags(tags=NULL,exclude_tags=c("deprecated","readability"))
mylinters <- linters_with_defaults(defaults=mylinters,
                                   return_linter=NULL,
                                   trailing_whitespace_linter=NULL,
                                   trailing_blank_lines_linter=NULL,
                                   condition_call_linter=NULL,
                                   object_name_linter=NULL,
                                   implicit_integer_linter=NULL)
lint_package(path=".",linters=mylinters)




