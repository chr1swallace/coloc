#' @importFrom graphics points
globalVariables(c("variable","pp","position","value"))

                                        # ds1: dataset1
# ds2: dataset2
# both ds1 and ds2 should contain the same snps in the same order

# Pos: positions of all snps in ds1 or in ds2
# chr: Chromosome
# pos.start: lower bound of positions
# pos.end: upper bound of positions
# trait1: name of trait 1
# trait2: name of trait 2   
ymin <- NULL
##'
##' @title plot a coloc dataset
##' @param d a coloc dataset
##' @param susie_obj optional, the output of a call to runsusie()
##' @param highlight_list optional, a list of character vectors. any snp in the
##'   character vector will be highlighted, using a different colour for each
##'   list.
##' @param alty default is to plot a standard manhattan. If you wish to plot a
##'   different y value, pass it here. You may also want to change ylab to
##'   describe what you are plotting.
##' @param ylab label for y axis, default is -log10(p) and assumes you are
##'   plotting a manhattan
##' @param show_legend optional, show the legend or not. default is TRUE
##' @param color optional, specify the colours to use for each credible set when
##'   susie_obj is supplied. Default is shamelessly copied from
##'   susieR::susie_plot() so that colours will match
##' @param ... other arguments passed to the base graphics plot() function
##' @author Chris Wallace
##' @export
plot_dataset <- function(d,
                         susie_obj=NULL,
                         highlight_list=NULL,
                         alty=NULL,ylab="-log10(p)",
                         show_legend=TRUE,
                         color = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                                   "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
                                   "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
                                   "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
                                   "yellow3", "darkorange4", "brown"),
                         ...) {
  if(!("position" %in% names(d)))
    stop("no position element given")
  if(is.null(alty)) {
    y= - ( pnorm(-abs(d$beta)/sqrt(d$varbeta),log.p=TRUE) + log(2) ) / log(10)
  } else {
    y=alty
  }
  plot(d$position,y,xlab="Position",ylab=ylab,pch=16,col="grey",...)
  if(!is.null(susie_obj))
    highlight_list=lapply(susie_obj$sets$cs, names)
  if(!is.null(highlight_list)) {
    for(i in 1:length(highlight_list)) {
      w=which(d$snp %in% highlight_[[i]])
      points(d$position[w], y[w], col=color[i], cex=2)
    }
    if(show_legend)
      legend("topright",col=color[1:length(highlight_list)],
             pch=rep(1,length(highlight_list)),legend=1:length(highlight_list))
  }
}

ymax <- NULL

##' Plot a coloc structured dataset
##'
##' @title plot a coloc dataset
##' @param d a coloc dataset
##' @param susie_obj optional, the output of a call to runsusie()
##' @param highlight_list optional, a list of character vectors. any snp in the
##'   character vector will be highlighted, using a different colour for each
##'   list.
##' @param alty default is to plot a standard manhattan. If you wish to plot a
##'   different y value, pass it here. You may also want to change ylab to
##'   describe what you are plotting.
##' @param ylab label for y axis, default is -log10(p) and assumes you are
##'   plotting a manhattan
##' @param show_legend optional, show the legend or not. default is TRUE
##' @param color optional, specify the colours to use for each credible set when
##'   susie_obj is supplied. Default is shamelessly copied from
##'   susieR::susie_plot() so that colours will match
##' @param ... other arguments passed to the base graphics plot() function
##' @author Chris Wallace
##' @export
plot_dataset <- function(d,
                         susie_obj=NULL,
                         highlight_list=NULL,
                         alty=NULL,ylab="-log10(p)",
                         show_legend=TRUE,
                         color = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                                   "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
                                   "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
                                   "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
                                   "yellow3", "darkorange4", "brown"),
                         ...) {
  if(!("position" %in% names(d)))
    stop("no position element given")
  if(is.null(alty)) {
    y= - ( pnorm(-abs(d$beta)/sqrt(d$varbeta),log.p=TRUE) + log(2) ) / log(10)
  } else {
    y=alty
  }
  plot(d$position,y,xlab="Position",ylab=ylab,pch=16,col="grey",...)
  if(!is.null(susie_obj))
    highlight_list=lapply(susie_obj$sets$cs, names)
  if(!is.null(highlight_list)) {
    for(i in 1:length(highlight_list)) {
      w=which(d$snp %in% highlight_list[[i]])
      points(d$position[w], y[w], col=color[i], cex=2)
    }
    if(show_legend)
      legend("topright",col=color[1:length(highlight_list)],
             pch=rep(1,length(highlight_list)),legend=1:length(highlight_list))
  }
}

##' Plot a pair of coloc datasets in an extended format
##' 
##' @title Draw extended plot of summary statistics for two coloc datasets
##'
##' @description Draw Manhattan-style locus plots for p-values in each dataset, gene annotations, 
##' a scatter plot of z-scores, and a table of coloc summary statistics.
##'
##' @details This function expects that the first two elements of the coloc 
##' dataset list d contain summary statistics. Columns labelled 'chromosome',
##' 'position', and 'p' (p-values) are expected in each dataset. The packages \code{locuszoomr}
##' and one of \code{EnsDb.Hsapiens.v86} or \code{EnsDb.Hsapiens.v75} are required
##' for this function. One of the \code{EnsDb.Hsapiens} libraries containing the Ensembl database
##'  specified by ens_db must be loaded prior to  use (see the example). EnsDb.Hsapiens.v86
##' is the default database (GRCh38/hg19); for GRCh37/hg19, use EnsDb.Hsapiens.v75.
##'
##' @param d a coloc dataset
##' @param x object of class \code{coloc_abf} returned by \code{coloc.abf()}
##' @param first_highlight_list character vector of snps to highlight in the first dataset; 
##' \code{'index'} can be passed to highlight the most significant SNP
##' @param second_highlight_list character vector of snps to highlight in the second dataset;
##' \code{'index'} can be passed to highlight the most significant SNP
##' @param first_index_snp snp to designate as index snp in the first dataset 
##' for the purpose of LD visualisation
##' @param second_index_snp snp to designate as index snp in the second dataset
##' for the purpose of LD visualisation
##' @param first_trait label for the first trait
##' @param second_trait label for the second trait
##' @param snp_label label for the snp column
##' @param ld_label label for the LD column
##' @param show_ld logical, whether to show LD in the locus plots
##' @param locus_plot_ylim numeric vector of length 2 specifying the y-axis 
##' limits for the locus plots on the -log10 scale
##' @param ens_db character string specifying Ensembl database package from which to get gene positions
##' @return a gtable object
##'
##' @importFrom ggplot2 ggplot geom_point geom_vline geom_hline coord_fixed xlim ylim
##' @importFrom data.table data.table
##' @importFrom gridExtra arrangeGrob tableGrob ttheme_minimal
##' @export
##' @author Tom Willis
##' @examples
##' \dontrun{
##' library(coloc)
##' library(EnsDb.Hsapiens.v86)
##' plot(plot_extended_dataset(list(first_dataset, second_dataset), coloc,
##' first_highlight_list = c("rs123", "rs456"), 
##' second_highlight_list = c("rs789", "rs1011"), 
##' ens_db = "EnsDb.Hsapiens.v86"))
##' }
plot_extended_dataset <- function(d, 
                                  x, 
                                  first_highlight_list = NULL, 
                                  second_highlight_list = NULL, 
                                  first_index_snp = NULL,
                                  second_index_snp = NULL,
                                  first_trait = "Trait 1",
                                  second_trait = "Trait 2",
                                  snp_label = "snp",
                                  ld_label = NULL,
                                  show_ld = FALSE,
                                  locus_plot_ylim = NULL,
                                  ens_db = "EnsDb.Hsapiens.v86") {
  if(!requireNamespace("locuszoomr", quietly = TRUE)) {
    stop("locuszoomr package is required for this function.")
  } 
 
  if(ens_db == "EnsDb.Hsapiens.v86") {
    if(!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      stop("EnsDb.Hsapiens.v86 package is required for this function.")
    }
  } else if(ens_db == "EnsDb.Hsapiens.v75") {
    if(!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) {
      stop("EnsDb.Hsapiens.v75 package is required for this function.")
    }
  } else {
    stop("Invalid Ensembl database specified.")
  }
  
  # From locuszoomr::locus
  if (!ens_db %in% (.packages())) {
    stop("Ensembl database not loaded. Try: library(", ens_db, ")",
          call. = FALSE)
  }
  
  # See https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html;
  # tl;dr it stops the 'undefined global variables' warning in R CMD check
  p <- z.x <- z.y <- pretty_value <- NULL  
  if(!("chromosome" %in% names(d[[1]]) & "chromosome" %in% names(d[[1]]))) {
    stop("A \"chromosome\" column is required in both datasets.")
  }

  if(!("position" %in% names(d[[1]]) & "position" %in% names(d[[1]]))) {
    stop("A \"position\" column is required in both datasets.")
  }

  if(!("p" %in% names(d[[1]]) & "p" %in% names(d[[1]]))) {
    stop("A \"p\" column is required in both datasets.")
  }

  if(length(unique(d[[1]]$chromosome)) != 1 | 
     length(unique(d[[2]]$chromosome)) != 1 | 
     unique(d[[1]]$chromosome) != unique(d[[2]]$chromosome)) {
    print(unique(d[[1]]$chromosome))
    print(unique(d[[2]]$chromosome))
    stop("Dataset should contain summary statistics for variants on only one chromosome, 
    which should be the same in both datasets.")
  } else {
    chromosome <- unique(d[[1]]$chromosome)
  }
  
  if(show_ld & is.null(ld_label)) {
      stop("If show_ld is TRUE, ld_label must be specified.")
   }

  first_dataset <- as.data.table(d[[1]])
  first_dataset[, p := as.numeric(p)]
  first_dataset <- first_dataset[p > 0 & !is.na(p)]
  first_dataset[, z := beta/sqrt(varbeta)]
  
  second_dataset <- as.data.table(d[[2]])
  second_dataset[, p := as.numeric(p)]
  second_dataset <- second_dataset[p > 0 & !is.na(p)]
  second_dataset[, z := beta/sqrt(varbeta)]

  min_pos <- min(first_dataset[, min(position)], second_dataset[, min(position)])
  max_pos <- max(first_dataset[, max(position)], second_dataset[, max(position)])

  min_z <- min(first_dataset[, min(z)], second_dataset[, min(z)])
  max_z <- max(first_dataset[, max(z)], second_dataset[, max(z)])

  first_loc  <- locuszoomr::locus(data = first_dataset, chrom = "chromosome", 
    pos = "position", p = "p", labs = snp_label, seqname = chromosome, 
    xrange = c(min_pos, max_pos), LD = ld_label, index_snp = first_index_snp,
    ens_db = ens_db)
  second_loc <- locuszoomr::locus(data = second_dataset, chrom = "chromosome", 
    pos = "position", p = "p", labs = snp_label, seqname = chromosome, 
    xrange = c(min_pos, max_pos), LD = ld_label, index_snp = second_index_snp,
    ens_db = ens_db)
  
  if(is.null(locus_plot_ylim)) {
    # If we have any highlights we may have to make space for them by adding a 
    # few log units to the y axis
    first_locus_plot_ylim <- c(0, max(first_loc$data$logP, na.rm = T))
    second_locus_plot_ylim <- c(0, max(second_loc$data$logP, na.rm = T))
   if(!is.null(first_highlight_list)) {
     first_locus_plot_ylim[2] <- first_locus_plot_ylim[2] + 8 
   }
   if(!is.null(second_highlight_list)) {
     second_locus_plot_ylim[2] <- second_locus_plot_ylim[2] + 8 
   }
  } else {
    first_locus_plot_ylim <- locus_plot_ylim
    second_locus_plot_ylim <- locus_plot_ylim
  }

  pls <- list()

  pls[[1]] <- locuszoomr::gg_scatter(first_loc, min.segment.length = 0,
   nudge_y = 5, labels = first_highlight_list, showLD = show_ld,
   legend_pos = 'topright', ylim = first_locus_plot_ylim)+
    ggtitle(first_trait)
  pls[[2]] <- locuszoomr::gg_scatter(second_loc, min.segment.length = 0,
   nudge_y = 5, labels = second_highlight_list, showLD = show_ld,
   legend_pos = 'topright', ylim = second_locus_plot_ylim)+
    ggtitle(second_trait)
  pls[[3]] <- locuszoomr::gg_genetracks(first_loc)

  merged <- merge(first_dataset[, .(snp, z), env = list(snp = snp_label)], 
    second_dataset[, .(snp, z), env = list(snp = snp_label)], 
    by = snp_label, suffixes = c(".x", ".y"))

  z_cor <- merged[, cor(z.x, z.y, use = "pairwise.complete.obs")]

  ggplot(merged)+
    geom_point(aes(x = z.x, y = z.y))+
    xlab(first_trait)+
    ylab(second_trait)+
    geom_vline(xintercept = qnorm(2.5e-8), linetype = "dashed", col = "blue")+
    geom_vline(xintercept = qnorm(2.5e-8, lower.tail = F), linetype = "dashed", col = "blue")+
    geom_hline(yintercept = qnorm(2.5e-8), linetype = "dashed", col = "blue")+
    geom_hline(yintercept = qnorm(2.5e-8, lower.tail = F), linetype = "dashed", col = "blue")+
    coord_fixed(ratio = 1)+
    xlim(c(min_z, max_z))+
    ylim(c(min_z, max_z))+
    ggtitle(sprintf('Pearson correlation: %.3f', z_cor))+
    theme_bw() -> pls[[4]]
    
  coloc_summary <- data.table(variable = names(x$summary), value = x$summary)
  coloc_summary[variable == "nsnps",  pretty_value := format(value, big.mark = ",")]
  coloc_summary[variable != "nsnps",  pretty_value := signif(value, digits =  2)]

  pls[[5]] <- tableGrob(coloc_summary[, .(Variable = variable, Value = pretty_value)],
   theme = ttheme_minimal(), rows = NULL)

  arrangeGrob(grobs = pls, layout_matrix = cbind(1:3, c(4,4,5)))
}

##' Print summary of a coloc.abf run
##'
##' @title print.coloc_abf
##' @param x object of class \code{coloc_abf} returned by coloc.abf()
##'   or coloc.signals()
##' @param ... optional arguments: "trait1" name of trait 1, "trait2"
##'   name of trait 2
##' @return x, invisibly
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname print-methods
print.coloc_abf <- function(x,...) {
  trait1="trait 1"
  trait2="trait 2"
  args <- list(...)
  if("trait1" %in% names(args))
    trait1=args$trait1
  if("trait2" %in% names(args))
    trait2=args$trait2
  
    message("Coloc analysis of ",trait1,", ",trait2)
    message("\nSNP Priors")
    print(x$priors)
    message("\nHypothesis Priors")
    ns <- if(is.data.table(x$summary)) { x$summary$nsnps[1] } else { x$summary["nsnps"] }
    hprior <- prior.snp2hyp(nsnp=ns,
                        p1=x$priors["p1"],
                        p2=x$priors["p2"],
                        p12=x$priors["p12"])
    rownames(hprior) <- ""
    print(hprior)
    message("\nPosterior")
    if(is.data.table(x$summary)) {
        summ <- copy(x$summary)
        setnames(summ, gsub("PP.|.abf","",names(summ)))
        print(summ[,list(nsnps,hit1,hit2,H0,H1,H2,H3,H4)])
    } else {
        summ <- x$summary
        names(summ) <- gsub("PP.|.abf","",names(summ))
        print(summ)
    }
    invisible(x)
}

##' plot a coloc_abf object
##'
##' @title plot a coloc_abf object
##' @param x coloc_abf object to be plotted
##' @param ... other arguments
##' @return ggplot object
##' @docType methods
##' @export
##' @rdname plot-methods
##' @author Chris Wallace
plot.coloc_abf <- function(x,...) {
  x=as.data.table(x$results)
  if(!("position" %in% names(x)))
    x$position <- 1:nrow(x)
  m= melt(x,id.vars=c("snp","position"), measure.vars=grep("z.df",names(x),value=TRUE))
  m2=melt(x,id.vars=c("snp","position"), measure.vars=grep("PP",names(x),value=TRUE))
  setnames(m2,"value","pp")
  if(length(grep("row",m$variable))) {
    m[,c("z","df","row"):=tstrsplit(variable,"\\.")]
    m2[,row:=sub(".*\\.","",variable)]
    m <- merge(m,m2[,.(snp,row,pp)],by=c("snp","row"))
    ggplot(m, aes(x=position,y=abs(value),col=pp,size=pp)) +
      geom_point() +
      facet_grid(df ~ row)
  } else {
    m[,c("z","df"):=tstrsplit(variable,"\\.")]
    m <- merge(m,m2[,.(snp,pp)],by=c("snp"))
    ggplot(m, aes(x=position,y=abs(value),col=pp,size=pp)) +
      geom_point() +
      facet_grid(df ~ .)
  }
}

