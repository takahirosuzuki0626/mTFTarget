#' Number of significant targets
#' 
#' This function compute the number of significant targets of methyl-regulating TFs in differentranges.
#' Number of motifs in a region is dependent of the rength of the region. 
#' To potimize the range of the region, this function return the number of significatn targets changing the range.
#' Becaue fitting function is not robust to outliers, if a fitting return the error, repret the fitting by 100 times.
#' output figures consist of histgram of motif number at DMR and at random regions, fitting model line plot for each results.
#' and return tabel of the result

#' @param ranges a vector. ranegs to be tested.
#' @param eg_connection_list a list of connection data between enhancers of the methyl-regulating TF targets and thir downstream gene. the data have to have score at 3rd column.
#' @param sig_genes_list  a list of vectors of the methyl-regulating TF target promoters
#' @param outname name for output file
#'
#' @import GenomicRanges
#' @importFrom igraph graph.data.frame simplify V layout_with_fr V<- E E<- layout.kamada.kawai layout.reingold.tilford
#' @importFrom graphics legend plot
#' @importFrom grDevices adjustcolor dev.off pdf
#'
#' @return a dataframe of connections (and export a network plot as a pdf)
#'
#' @keywords binominal fitting optimization
#' @export


ranges <- seq(100, 2000, length=20)




function ( ranges, )

totalvstaraget <- NULL
for(i in ranges){
    fiterror <- TRUE
    fitreps <- 0
    while(fiterror == TRUE){    # repeat if the fitting is error
        ## if while loop repeated more than 100 time. break the loop
        if (fitreps > 100){
            cat("fitting failed!\n")
            totalvstaraget <- rbind(totalvstaraget, rep("NA", 6))
            next
        }

        ## parametors setting
        outname = paste0("RUNX1OE_demet_", i)
        seq_range = c(-i, i)

        ## count the motifs in the indicated range
        motif_positionsList <- MotPosiList(infile=infile,
                                        motifDBList=motifDBList,
                                        ControlColnum = ControlColnum,
                                        TreatmentColnum = TreatmentColnum,
                                        MethylDemethyl=MethylDemethyl,
                                        version = version,
                                        seq_range = seq_range,
                                        outname = outname)

        target_positionsList <- motif_positionsList[[1]]
        random_positionsList <- motif_positionsList[[2]]

        target_motif_num_each <- data.frame(target=sapply(target_positionsList[[1]], length))
        random_motif_num_each <- data.frame(random=sapply(random_positionsList[[1]], length))

        ## fitting 
        library(fitdistrplus)
        fit <- try(fitdist(unlist(random_motif_num_each), lower=c(0,0), distr="nbinom"))
        cat("OK\n")
        if(class(fit) == "try-error"){    # fitting error check
            cat ("fitting error! repeat the fitting\n")
            fiterror <- TRUE
            fitreps <- fitreps + 1
        }else{
            fiterror <- FALSE    # to break the while loop
            size=fit$estimate["size"]    # siz  e of fitted distribution
            mu=fit$estimate["mu"]    # mu of fitted distribution
            zero_p <- pnbinom(0, size=size, mu=mu, lower.tail = FALSE, log.p = FALSE)    #p-value of zero
            nbinom_pval <- DMRNbinomTest(target_positionsList,random_positionsList, outname)    # computation of p-value
            sig_targets <- nbinom_pval[nbinom_pval <= nbiom_cutoff,]    # significant targets
            nsig_targets <- length(sig_targets)    # number of significant targets
        
            totalvstaraget <- rbind(totalvstaraget, c(i, size=size, mu=mu, zero_p, length(motif_positionsList[[1]][[1]]), nsig_targets))
        }
    }
}
colnames(totalvstaraget) <- c("range", "size", "mu", "zero_pval", "total", "target")