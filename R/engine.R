## ---------------------------------------------------------------------------
## The colocalisation engine.
##
## coloc.bf_bf() reorganised around one idea: every degenerate case sets
## $status instead of returning a differently shaped object. The arithmetic is
## unchanged. combine.abf()'s five-hypothesis enumeration, the SNP.PP.H4
## softmax and logsum/logdiff all carry over exactly as they are.
##
## The engine's contract is that bf1 and bf2 are the ONLY trait-specific
## information it ever sees. It knows nothing about datasets, LD, susieR or
## how the numbers were made, and that ignorance is what lets per-trait
## methods mix freely.
##
## CITATIONS. A file:line reference in this branch's own code always means
## 6.0.3, at the pinned commit 8f20f0b, which never moves and so can be cited
## precisely forever. Live code in this tree is named by symbol instead, never
## by line, because a line number against our own files rots the next time
## anyone inserts above it.
## ---------------------------------------------------------------------------

## The defaults live in ONE place: the constructors in R/prior.R. Re-declaring
## them here is how 6.0.3 ended up with p12 = 1e-5 in coloc.abf() and 5e-6 in
## coloc.bf_bf() (allowlist AL01), so the engine derives them instead.
##
## The OLD doors keep the literals in their own signatures on purpose: those
## are 6.0.3's published API, frozen until the hard break, and they are always
## passed down explicitly rather than falling through to these.
.default_prior   <- function() unclass(coloc_prior())
.default_control <- function() unclass(coloc_control())

.coloc_engine <- function(bf1, bf2,
                          prior = .default_prior(),
                          control = .default_control()) {

    msgs <- character(0)
    ## keep emitting warning() as well as recording $messages: users rely on
    ## the warnings and the harness compares conditions, while the structured
    ## copy lets someone running ten thousand pairs ask which were trimmed
    ## without parsing console output
    note <- function(...) {
        m <- paste0(...)
        warning(m, call. = FALSE)
        msgs <<- c(msgs, m)
    }

    p1 <- prior$p1; p2 <- prior$p2; p12 <- prior$p12
    w1 <- prior$weights1; w2 <- prior$weights2

    ## ---- 0. invariants the adapters are supposed to have established ------
    stopifnot(is.matrix(bf1), is.matrix(bf2),
              !is.null(colnames(bf1)), !is.null(colnames(bf2)))
    if (!is.null(w1) && ncol(bf1) != length(w1))
        stop("Length of prior_weights1 must match size of dataset 1")
    if (!is.null(w2) && ncol(bf2) != length(w2))
        stop("Length of prior_weights2 must match size of dataset 2")

    ## ---- 1. an adapter returned no signals --------------------------------
    ## Newly reachable: 6.0.3 could not produce a zero-row matrix, and its
    ## 1:nrow(bf1) would yield c(1, 0) and a bogus pair grid.
    if (nrow(bf1) == 0L || nrow(bf2) == 0L) {
        note("at least one trait has no signals, nothing to colocalise")
        return(.new_coloc_result(status = "no_signal", messages = msgs, prior = prior))
    }

    ## ---- 2. snp intersection ----------------------------------------------
    isnps <- setdiff(intersect(colnames(bf1), colnames(bf2)), "null")
    if (!length(isnps)) {
        note("no snps in common between the two traits")
        return(.new_coloc_result(status = "no_overlap", messages = msgs, prior = prior))
    }

    ## ---- 3. null column convention, tested per matrix ---------------------
    nulllast <- function(m) m[, c(setdiff(colnames(m), "null"), "null"), drop = FALSE]
    if ("null" %in% colnames(bf1))
        bf1 <- nulllast(bf1) - matrix(bf1[, "null"], nrow(bf1), ncol(bf1))
    if ("null" %in% colnames(bf2))
        bf2 <- nulllast(bf2) - matrix(bf2[, "null"], nrow(bf2), ncol(bf2))

    ## ---- 4. posteriors and coverage, on the UNRESTRICTED matrices ---------
    pp1 <- logbf_to_pp(bf1, p1, last_is_null = "null" %in% colnames(bf1))
    pp2 <- logbf_to_pp(bf2, p2, last_is_null = "null" %in% colnames(bf2))
    ph0 <- function(pp) if ("null" %in% colnames(pp)) pp[, "null"] else 1 - rowSums(pp)
    prop <- function(pp) rowSums(pp[, isnps, drop = FALSE]) /
                         rowSums(pp[, setdiff(colnames(pp), "null"), drop = FALSE])
    ph0.1 <- ph0(pp1); ph0.2 <- ph0(pp2)
    prop1 <- prop(pp1); prop2 <- prop(pp2)

    ## ---- 5. the pair grid, and trimming -----------------------------------
    todo <- as.data.table(expand.grid(i = seq_len(nrow(bf1)), j = seq_len(nrow(bf2))))
    n_in <- nrow(todo)
    drop <- if (isTRUE(control$trim_by_posterior))
                prop1[todo$i] < control$overlap.min | prop2[todo$j] < control$overlap.min
            else rep(FALSE, n_in)
    ## Coverage is unassessable when prop is NA, and that takes only ONE NA
    ## anywhere in a bf row, because prop is a ratio of rowSums. Let the loop
    ## mark such a pair degenerate rather than trying to trim on it.
    ##
    ## 6.0.3 has no such line and errors on ANY NA in drop, via whichever of
    ## its two guards catches it: all() short-circuits on FALSE and any() on
    ## TRUE, so a lone NA trips if(all(drop)) at R/susie.R:285 while
    ## c(NA, FALSE) gets past that one and trips if(any(drop)) at :297.
    ## Either way its own all-NA branch at :350 is unreachable under the
    ## default trim_by_posterior = TRUE.
    drop[is.na(drop)] <- FALSE

    if (all(drop)) {
        note("snp overlap too small between datasets: too few snps with high ",
             "posterior in one trait represented in other")
        ## 6.0.3 names the hit from the POSTERIOR here and from the BF on the
        ## ok path below. Those disagree whenever a trait's best snp is not in
        ## the overlap, which is exactly this situation. Preserved on purpose:
        ## unifying them would silently move a reported column.
        pk <- function(pp, ix) colnames(pp)[apply(pp, 1, which.max)][ix]
        skel <- data.table(nsnps = length(isnps),
                           hit1 = pk(pp1, todo$i), hit2 = pk(pp2, todo$j),
                           PP.H0.abf = pmin(ph0.1[todo$i], ph0.2[todo$j]),
                           PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
                           PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
                           idx1 = todo$i, idx2 = todo$j)
        return(.new_coloc_result(summary = skel, status = "all_trimmed",
                                 messages = msgs, prior = prior,
                                 meta = list(n_pairs_in = n_in, n_pairs_kept = 0L,
                                             dropped = "overlap.min")))
    }
    if (any(drop))
        todo <- todo[!drop]

    ## ---- 6. restrict to the intersection, carrying per-snp priors ---------
    if (length(p12) > 1 && (length(isnps) != ncol(bf1) || length(isnps) != ncol(bf2)))
        stop("p12 must have length 1 or equal length to p1 and p2, but different numbers of snps between datasets")
    restrict <- function(bf, p, w) {
        if (length(isnps) == ncol(bf)) return(list(bf = bf, p = p, w = w))
        keep <- match(isnps, colnames(bf))
        list(bf = bf[, keep, drop = FALSE],
             p = if (length(p) > 1) p[keep] else p,
             w = w[keep])
    }
    r1 <- restrict(bf1, p1, w1); bf1 <- r1$bf; p1 <- r1$p; w1 <- r1$w
    r2 <- restrict(bf2, p2, w2); bf2 <- r2$bf; p2 <- r2$p; w2 <- r2$w

    ## Which combining function. combine.abf() is scalar-only: given a per-snp
    ## p1 it recycles rather than collapses, which is how 6.0.3 turned length-8
    ## priors into a 33-hypothesis summary (allowlist AL10).
    ## combine_abf_weighted() collapses correctly with logsum(), but it
    ## reassociates the H3 term, so it is NOT bit-identical to combine.abf()
    ## on the scalar path even though it is algebraically equal. D6 says the
    ## scalar path keeps the function it has always used.
    user_weighted <- !is.null(w1) || !is.null(w2)
    per_snp_prior <- length(p1) > 1 || length(p2) > 1
    use_weighted  <- user_weighted || per_snp_prior

    if (per_snp_prior) {
        if (length(p12) > 1)
            stop("p12 must be a single value when p1 or p2 is per-snp")
        ## re-express a per-snp prior as (scalar, weights), which is exact:
        ## combine_abf_weighted() builds Q * p * w/sum(w), and with w = p and
        ## scalar mean(p) that is p itself
        as_weights <- function(p, w) {
            if (length(p) > 1) list(p = mean(p), w = p) else list(p = p, w = w)
        }
        c1 <- as_weights(p1, w1); p1 <- c1$p; w1 <- c1$w
        c2 <- as_weights(p2, w2); p2 <- c2$p; w2 <- c2$w
    }

    ## ---- 7. one signal pair at a time -------------------------------------
    hit_of <- function(bf, i) {
        w <- which.max(bf[i, ])
        if (!length(w)) NA_character_ else colnames(bf)[w]
    }
    pps  <- matrix(NA_real_, nrow(todo), 5L)
    snpp <- vector("list", nrow(todo))
    h1 <- h2 <- character(nrow(todo))
    degenerate <- FALSE

    for (k in seq_len(nrow(todo))) {
        l1 <- bf1[todo$i[k], isnps]
        l2 <- bf2[todo$j[k], isnps]
        lsum <- l1 + l2
        sp <- exp(lsum - logsum(lsum))
        pp <- if (use_weighted)
                  combine_abf_weighted(l1, l2, p1, p2, p12, w1, w2, quiet = TRUE)
              else
                  combine.abf(l1, l2, p1, p2, p12, quiet = TRUE)
        h1[k] <- hit_of(bf1, todo$i[k])
        h2[k] <- hit_of(bf2, todo$j[k])
        ## 6.0.3 rewrote these to certainties: an all-NA pair became
        ## PP.H0 = 1, then immediately PP.H2 = 1 because which.max on an
        ## all-NA row returns nothing. Under the governing rule, a posterior
        ## that could not be computed is NA.
        if (all(is.na(sp)) || is.na(h1[k]) || is.na(h2[k])) {
            degenerate <- TRUE
            sp <- rep(NA_real_, length(isnps))
            pp <- rep(NA_real_, 5L)
        }
        snpp[[k]] <- sp
        pps[k, ] <- unname(pp)
    }

    ## ---- 8. one shape out --------------------------------------------------
    skel <- data.table(nsnps = length(isnps), hit1 = h1, hit2 = h2,
                       PP.H0.abf = pps[, 1], PP.H1.abf = pps[, 2],
                       PP.H2.abf = pps[, 3], PP.H3.abf = pps[, 4],
                       PP.H4.abf = pps[, 5],
                       idx1 = todo$i, idx2 = todo$j)
    PP <- as.data.table(do.call(cbind, snpp))
    colnames(PP) <- if (ncol(PP) > 1) paste0("SNP.PP.H4.row", seq_len(ncol(PP)))
                    else "SNP.PP.H4.abf"

    .new_coloc_result(summary = skel,
                      results = cbind(data.table(snp = isnps), PP),
                      prior = prior,
                      weights = if (user_weighted) list(w1 = w1, w2 = w2),
                      status = if (degenerate) "degenerate" else "ok",
                      messages = msgs,
                      meta = list(n_pairs_in = n_in, n_pairs_kept = nrow(todo),
                                  dropped = if (any(drop)) "overlap.min" else character(0)))
}
