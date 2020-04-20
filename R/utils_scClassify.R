
# A function coded by Hadley Wickham to suppress print/cat output
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}
