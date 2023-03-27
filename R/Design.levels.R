Design.levels <- function(df, at) {
    ac <- at$assume.code
    for (nn in names(df)) {
        j <- match(nn, at$name, 0)
        if (j > 0) {
            if ((ac[j] == 5 | ac[j] == 8) & length(lev <- at$parms[[nn]])) {
                df[[nn]] <- factor(df[[nn]], lev)
            }
        }
    }
    df
}
