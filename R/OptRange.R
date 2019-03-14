#' Number of significant targets
#' 
#' This function compute the number of significant targets of methyl-regulating TFs in differentranges.
#' Number of motifs in a region is dependent of the rength of the region. 
#' To potimize the range of the region, this function return the number of significatn targets changing the range.
#' Becaue fitting function is not robust to outliers, if a fitting return the error, repret the fitting by 100 times.
#' output figures consist of histgram of motif number at DMR and at random regions, fitting model line plot for each results.
#' and return tabel of the result

#' @param ranges a vector. ranegs to be tested.
#' @param infile input m-value file
#' @param motifDBList list of motif PWM
#' @param ControlColnum control column index(s)
#' @param TreatmentColnum treatment column index(s)
#' @param MethylDemethyl Metylation or Demethylation analysis
#' @param version version of methylation array 450 or EPIC(850)
#' @param nbiom_cutoff cutoff p-value for negative binomial test
#' @param outname output file name
#'
#' @import GenomicRanges
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pnbinom
#'
#' @return a matrix of results of fitting and negative binominal test
#'
#' @keywords binominal fitting optimization
#' @export

OptRange <- function (infile,
                      motifDBList,
                      ControlColnum,
                      TreatmentColnum,
                      MethylDemethyl,
                      version,
                      outname,
                      nbiom_cutoff,
                      ranges){
    totalvstaraget <- NULL
    for(i in ranges){
        cat(paste0("---------- Range:",i," ----------\n"))
        fiterror <- TRUE
        fitreps <- 0
        while (fiterror == TRUE){    # repeat if the fitting is error
            ## if while loop repeated more than 100 time. break the loop
            if(fitreps > 100){
                cat("fitting failed!\n")
                totalvstaraget <- rbind(totalvstaraget, rep("NA", 6))
                next
            }

            ## parametors setting
            outname = paste0(outname, i)
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
            #library(fitdistrplus)
            fit <- try(fitdist(unlist(random_motif_num_each), lower=c(0,0), distr="nbinom"))
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
                cat("nbinominal exact test Done...\n")
            }
        }
    }
    colnames(totalvstaraget) <- c("range", "size", "mu", "zero_pval", "total", "target")
    return(totalvstaraget)
}