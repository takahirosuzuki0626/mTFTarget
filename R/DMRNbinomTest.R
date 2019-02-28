#' negative binomial test
#'
#' This function performs negative binomial test for methyl-regulating-TF target
#' @param target_positionsList list of lists of motif positions
#' @param random_positionsList list of lost of random positions
#' @param outname name for output file
#'
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom fitdistrplus fitdist
#' @importFrom stats dnbinom pnbinom
#' @importFrom rlang enquos quo_is_missing
#'
#' @return result of negative binomial exact test (a negative binomial test fitting plot pdf is also exported)

#'
#' @keywords gff3 metadata
#' @export

DMRNbinomTest <- function (target_positionsList,random_positionsList, outname){
    ##To avoid note in devtools::heck()
    value <- NULL
    key <- NULL
    num <- NULL
    nb <- NULL

    ## negative binomial fitting
    target_motif_num_each <- data.frame(target=sapply(target_positionsList[[1]], length))
    random_motif_num_each <- data.frame(random=sapply(random_positionsList[[1]], length))

    #library(fitdistrplus)
    fit <- fitdist(unlist(random_motif_num_each), lower=c(0,0), distr="nbinom")
    random_simulation <- data.frame(num = 0:10, nb = dnbinom(0:10, size=fit$estimate["size"], mu=fit$estimate["mu"])*length(random_motif_num_each[[1]]))

    ## Plot of motif count distribution and negative binomial distribution model
    #library(ggplot2)
    df1 <- tidyr::gather(target_motif_num_each)
    df2 <- tidyr::gather(random_motif_num_each)
    df3 <- tidyr::gather(random_simulation)
    df <- rbind(df1,df2)

    theme <- theme(panel.background = element_blank(),    # initialization
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                plot.title = element_text(size = 11,hjust = 0.5),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12),
                axis.ticks = element_line(size=0.5),
                axis.line = element_line(size=0.5),
                plot.margin = unit(c(1,1,1,1),"line"),
                legend.title = element_blank())

    g <- ggplot(df)
    g <- g + geom_histogram(aes(x=value, color=key, fill=key), alpha=0.2, binwidth = 1, position = "identity")
    g <- g + scale_fill_manual(values = c("black", "blue"))
    g <- g + scale_color_manual(values = c("gray40", "blue4"))
    g <- g + geom_line(data=random_simulation, aes(x=num, y=nb), color="black", alpha=0.8)
    g <- g + ggtitle("Disributioin of motif number and Negative Binomial distribution Fitting")
    g <- g + xlab("Number of motif")
    g <- g + theme
    g <- g + theme(legend.position = c(0.8, 0.8))

    outfile_nbinom <- paste0(outname, "_nbinom_model.pdf")
    ggsave(filename=outfile_nbinom, plot=g)

    ## Computation of p-value based on exact test of negative binomial distribution model
    nbinom_pval <- sapply(target_motif_num_each, function(x){pnbinom(x, size=fit$estimate["size"], mu=fit$estimate["mu"], lower.tail = FALSE, log.p = FALSE)})
    rownames(nbinom_pval) <- rownames(target_motif_num_each)
    return(nbinom_pval)
}