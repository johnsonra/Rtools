# power.prop.test2.R
# Calculate power as in power.prop.test(), but allow for differing sample sizes
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created January 24, 2007
# Last Modified January 24, 2007

power.prop.test2 <- function (n1 = NULL, n2 = NULL, p1 = NULL, p2 = NULL,
                              sig.level = 0.05, power = NULL,
                              alternative = c("two.sided", "one.sided"), strict = FALSE) 
{
    if (sum(sapply(list(n1, n2, p1, p2, power, sig.level), is.null)) != 
        1) 
        stop("exactly one of 'n', 'p1', 'p2', 'power', and 'sig.level' must be NULL")
    
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
        sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    p.body <- quote(pnorm(((sqrt(n) * abs(p1 - p2) - (qnorm(sig.level/tside, 
        lower = FALSE) * sqrt((p1 + p2) * (1 - (p1 + p2)/2))))/sqrt(p1 * 
        (1 - p1) + p2 * (1 - p2)))))
    if (strict & tside == 2) 
        p.body <- quote({
            qu <- qnorm(sig.level/tside, lower = FALSE)
            d <- abs(p1 - p2)
            q1 <- 1 - p1
            q2 <- 1 - p2
            pbar <- (p1 + p2)/2
            qbar <- 1 - pbar
            v1 <- p1 * q1
            v2 <- p2 * q2
            vbar <- pbar * qbar
            pnorm((sqrt(n) * d - qu * sqrt(2 * vbar))/sqrt(v1 + 
                v2)) + pnorm((sqrt(n) * d + qu * sqrt(2 * vbar))/sqrt(v1 + 
                v2), lower = FALSE)
        })
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(1, 1e+07))$root
    else if (is.null(p1)) 
        p1 <- uniroot(function(p1) eval(p.body) - power, c(0, 
            p2))$root
    else if (is.null(p2)) 
        p2 <- uniroot(function(p2) eval(p.body) - power, c(p1, 
            1))$root
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10))$root
    else stop("internal error")
    NOTE <- "n is number in *each* group"
    METHOD <- "Two-sample comparison of proportions power calculation"
    structure(list(n = n, p1 = p1, p2 = p2, sig.level = sig.level, 
        power = power, alternative = alternative, note = NOTE, 
        method = METHOD), class = "power.htest")
}
