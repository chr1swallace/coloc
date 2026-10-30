globalVariables(c("grp", "lab", "prob", "id", "spl", "fac1", "fac2", "maxPP", "n_support", "shape", "show")) # data.table columns in plot_query()

##' Default colour palette for coloc group plots
##'
##' Paul Tol's "muted" qualitative palette, without its two palest colours (sand and cyan), which are hard to read as
##' text on white. It is designed to stay distinguishable under the common forms of colour blindness. See
##' <https://personal.sron.nl/~pault/>.
##'
##' @return character vector of 7 colours
##' @export
coloc_palette <- function() c("#332288", "#117733", "#44AA99", "#999933", "#CC6677", "#882255", "#AA4499")

## Resolve palette (unnamed vector, named vector, or function(n)) to one colour per level, named by level
.group_palette <- function(palette, lev) {
    if (is.function(palette)) palette <- palette(length(lev))
    if (is.null(names(palette))) return(stats::setNames(rep_len(palette, length(lev)), lev)) # recycled: neighbouring groups still differ
    miss <- setdiff(lev, names(palette))
    if (length(miss)) stop("palette has no colour for group(s): ", paste(miss, collapse = ", "), call. = FALSE)
    palette[lev]
}

## Abbreviations for point labels, applied in order (regular expressions, case-insensitive; longer phrases first)
.label_abbrev <- c(
    "\\s*\\(UKB data field [0-9]+\\)" = "", # identifiers carry no meaning on a plot
    "\\s*\\((?:Lipidomics ID|GC-TOF|PheCode|agilent_)[^)]*\\)" = "",
    "\\s*\\(X[0-9.]+\\)" = "", # protein assay ids
    "^ICD10 [A-Z0-9.]+:\\s*" = "",
    ",?\\s*inv-norm transformed" = "",
    "\\s*\\(Sex, Age[^)]*\\)" = "", # environment lists in interaction traits
    "commonly-tested environments" = "env.",
    "all environments" = "all env.",
    "\\binteraction\\b" = "int.",
    "\\bgene-based burden\\b" = "burden",
    "\\bbrain rsfmri measurement \\(functional connectivity\\)" = "rsFMRI connectivity",
    "\\s+protein levels?\\b" = " protein",
    "\\s+levels?\\b" = "",
    "\\bcerebrospinal fluid\\b" = "CSF",
    "\\bchronic kidney disease\\b" = "CKD",
    "\\bcoronary artery disease\\b" = "CAD",
    "\\binflammatory bowel disease\\b" = "IBD",
    "\\btype 1 diabetes\\b" = "T1D",
    "\\btype 2 diabetes\\b" = "T2D",
    "\\bbody mass index\\b" = "BMI",
    "\\bsystolic blood pressure\\b" = "SBP",
    "\\bdiastolic blood pressure\\b" = "DBP",
    "\\bblood pressure\\b" = "BP",
    "\\bred blood cell\\b" = "RBC",
    "\\bwhite blood cell\\b" = "WBC",
    "\\bvery low density lipoprotein\\b" = "VLDL",
    "\\blow density lipoprotein cholesterol\\b" = "LDL-C",
    "\\bhigh density lipoprotein cholesterol\\b" = "HDL-C",
    "\\blow density lipoprotein\\b" = "LDL",
    "\\bhigh density lipoprotein\\b" = "HDL",
    "\\bcholesteryl esters\\b" = "CE",
    "\\btriglycerides\\b" = "TG",
    "\\bcholesterol\\b" = "chol.",
    "\\bpercentage\\b" = "%",
    "\\bunidentified\\b" = "unid.",
    "\\bseropositivity for\\b" = "Seropos.",
    "\\babundance\\b" = "abund.",
    "\\brelative\\b" = "rel.",
    "\\bmeasurement\\b" = "meas.",
    "\\bfunctional connectivity\\b" = "func. conn.",
    "\\bhemisphere\\b" = "hemi.",
    "\\bmetabolite\\b" = "metab.")

## Shorten point labels: abbreviate (see .label_abbrev), then cut to width characters, at a word break where possible
.shorten_label <- function(x, width = 40) {
    for (i in seq_along(.label_abbrev)) x <- gsub(names(.label_abbrev)[i], .label_abbrev[[i]], x, ignore.case = TRUE, perl = TRUE)
    x <- trimws(gsub("\\s+", " ", x))
    if (is.null(width) || !is.finite(width)) return(x)
    long <- !is.na(x) & nchar(x) > width
    cut <- substr(x[long], 1, width - 1)
    brk <- regexpr("\\s\\S*$", cut) # last word break
    cut <- ifelse(brk > width / 2, substr(cut, 1, brk - 1), cut) # break at a word unless that loses over half the label
    x[long] <- sub("[[:space:][:punct:]]+$", "", cut)
    x
}

##' PheWAS-style plot of colocalisation results by group
##'
##' One point per gene and group (eg tissue or cell type): its height is the maximum PP.H4 over all results for that
##' gene in that group, and its size the number of results with PP.H4 above `threshold`. Points are hollow when no
##' result passes `threshold`. Groups are coloured, and the x axis labels are coloured to match.
##'
##' @param x data.frame of coloc results, eg from [coloc_query()] with a grouping column added
##' @param group column of x giving the group by which to group genes. default "group label" 
##'     are sorted alphabetically
##' @param label column of x used to label points. By default "label" if x has that column (as added by
##'     [annotate_ot()] and [annotate_ecat()]: gene symbol, or trait for GWAS), otherwise "symbol"
##' @param label_n label at most this many points per group in each panel: those with the highest pp. Unlabelled
##'     points are still plotted, and labels are placed to avoid them. Use Inf to label every point, or 0 for none
##' @param label_width point labels longer than this many characters are shortened: first by abbreviating common
##'     words and phrases in trait names (eg "cerebrospinal fluid" to "CSF", dropping identifiers such as "(UKB data
##'     field 21002)"), then by cutting at a word break and adding an ellipsis. NULL or Inf to only abbreviate
##' @param pp column of x containing the posterior probability to plot. default "PP.H4"
##' @param split_by column of x giving separate plots, or NULL for one plot. By default "query_cs" if x has that
##'     column (as output by [coloc_query_many()]), so there is one plot per index credible set
##' @param facet_by up to two columns of x by which to split each plot into panels, or NULL for a single panel. With
##'     one column, panels are side by side (or wrapped, if there are more than three); with two, the first gives rows
##'     and the second columns. Each panel shows only its own groups on the x axis, and is as wide as it needs to be. By default "study_class" if x has that column (as added by [annotate_ot()]), so
##'     GWAS and QTL results are side by side
##' @param count_by optional column of x (eg "dataset_id" or "studyId"). If given, size counts distinct values of this
##'     column above threshold, so that one dataset seen through several quantification methods counts once. If NULL,
##'     rows are counted
##' @param threshold results with pp above this count towards point size
##' @param min_pp genes whose maximum pp in a group is below this are not shown
##' @param palette colours for the groups: an unnamed vector (recycled over the groups in order), a vector named by
##'     group (use this to keep colours the same across several figures), or a function taking the number of groups,
##'     eg [grDevices::hcl.colors]. Default [coloc_palette()]
##' @param seed seed for the horizontal jitter, so the plot is reproducible
##' @return a ggplot object or, if split_by is set, a named list of ggplot objects, one per value of split_by
##'     (in order of first appearance in x), titled with that value. Values with nothing to plot are left out, with a
##'     message. Group colours and point size scales are the same across the list, so the plots can be compared or
##'     combined, eg with `patchwork::wrap_plots()`
##' @export
##' @author Chris Wallace
##' @examples
##' \dontrun{
##' plot_query(res, group = "group_label", count_by = "dataset_id")
##'
##' ## many credible sets: one plot each
##' res <- annotate_ot(coloc_query_many(cs_list, query_ot_many(cs_list)))
##' plots <- plot_query(res, count_by = "studyId")
##' plots[["BACH2"]]
##' patchwork::wrap_plots(plots[1:4], guides = "collect")
##' }
plot_query <- function(x, group = "group_label", label = if ("label" %in% names(x)) "label" else "symbol", pp = "PP.H4",
                              split_by = if ("query_cs" %in% names(x)) "query_cs", facet_by = intersect("study_class", names(x)),
                              label_n = 1, label_width = 40, count_by = NULL, threshold = 0.8, min_pp = 0.5, palette = coloc_palette(), seed = 1) {
    for (p in c("ggrepel", "ggtext"))
        if (!requireNamespace(p, quietly = TRUE)) stop("plot_query() needs the ", p, " package: install.packages(\"", p, "\")", call. = FALSE)
    need <- c(group, label, pp, count_by, split_by, facet_by)
    if (!is.data.frame(x) || !all(need %in% names(x))) stop("x must be a data.frame with columns: ", paste(need, collapse = ", "), call. = FALSE)
    if (length(split_by) > 1) stop("split_by must be a single column, or NULL", call. = FALSE)
    if (length(facet_by) > 2) stop("facet_by can have at most two columns", call. = FALSE)
    d <- data.table(grp = x[[group]], lab = x[[label]], prob = x[[pp]], id = if (is.null(count_by)) seq_len(nrow(x)) else x[[count_by]],
                    spl = if (is.null(split_by)) "" else as.character(x[[split_by]]),
                    fac1 = if (length(facet_by) >= 1) x[[facet_by[1]]] else "", fac2 = if (length(facet_by) == 2) x[[facet_by[2]]] else "")
    d <- d[!is.na(grp) & !is.na(lab) & !is.na(prob)]
    agg <- d[, list(maxPP = max(prob), n_support = uniqueN(id[prob > threshold])), by = c("spl", "grp", "lab", "fac1", "fac2")][maxPP >= min_pp]
    if (!nrow(agg)) stop("no ", label, " has ", pp, " >= min_pp (", min_pp, ")", call. = FALSE)

    ## shared across all plots, so colours and point sizes mean the same in each
    lev <- if (is.factor(x[[group]])) intersect(levels(x[[group]]), as.character(agg$grp)) else sort(unique(as.character(agg$grp)))
    agg[, grp := factor(as.character(grp), levels = lev)]
    agg[, shape := ifelse(n_support > 0, 16, 1)] # hollow: best result is below threshold
    agg[, lab := .shorten_label(as.character(lab), label_width)] # after aggregating, so shortening never merges different labels
    agg[, show := frank(-maxPP, ties.method = "first") <= label_n, by = c("spl", "grp", "fac1", "fac2")]
    agg[, lab := ifelse(show, lab, "")] # "" is not drawn, but ggrepel keeps other labels clear of the point
    pal <- .group_palette(palette, lev)
    esc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE), fixed = TRUE)
    nmax <- max(1, agg$n_support)
    br <- unique(c(0, pretty(c(1, nmax), n = 3)))
    br <- br[br >= 0 & br <= nmax & br == round(br)]
    pos <- position_jitter(width = 0.3, height = 0, seed = seed) # shared by points and labels so they stay together
    facets <- switch(length(facet_by) + 1, NULL,
                     if (uniqueN(agg$fac1) <= 3) facet_grid(~fac1, scales = "free_x", space = "free_x") # few panels (eg gwas/qtl): side by side, width by number of groups
                     else facet_wrap(~fac1, scales = "free_x"),
                     facet_grid(fac1 ~ fac2, scales = "free_x", space = "free_x"))

    one_plot <- function(a, title = NULL)
        ggplot(a, aes(x = .data$grp, y = .data$maxPP)) +
            geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey60") +
            geom_point(aes(size = .data$n_support, colour = .data$grp, shape = .data$shape), alpha = 0.8, position = pos) +
            ggrepel::geom_text_repel(aes(label = .data$lab), position = pos, size = 3, max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey70", colour = "grey20",
                                     direction = "y", force = 5, force_pull = 0.02, box.padding = 0.4, # spread labels up and down the panel, not just next to their points
                                     max.time = 1, max.iter = 1e5, seed = seed) +
            scale_shape_identity() +
            scale_colour_manual(values = pal, guide = "none") + # x axis labels already name the groups
            scale_size(range = c(1.5, 7), limits = c(0, nmax), breaks = br) +
            guides(size = guide_legend(override.aes = list(shape = ifelse(br > 0, 16, 1)))) +
            scale_x_discrete(labels = function(z) sprintf("<span style='color:%s'>%s</span>", pal[z], esc(z))) +
            scale_y_continuous(limits = c(min_pp, 1)) +
            labs(x = NULL, y = paste("max", pp), size = paste0("n with ", pp, "\n> ", threshold), title = title) +
            theme_bw() +
            theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1), panel.grid.major.x = element_blank()) +
            facets # NULL for a single panel: adding NULL to a ggplot does nothing

    if (is.null(split_by))
        return(one_plot(agg))
    keep <- unique(d$spl) # in order of first appearance in x
    empty <- setdiff(keep, agg$spl)
    if (length(empty))
        message(length(empty), " ", split_by, " value(s) with no ", label, " reaching min_pp (", min_pp, ") are not plotted: ", paste(utils::head(empty, 5), collapse = ", "), if (length(empty) > 5) ", ...")
    keep <- setdiff(keep, empty)
    stats::setNames(lapply(keep, function(k) one_plot(agg[agg$spl == k], title = k)), keep)
}
