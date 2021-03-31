#' Motif position list at DMRs
#'
#' This function performes differentially methylation analysis and motif enrichment analysis at the DMR
#' @param infile input m-value file
#' @param motifDBList list of motif PWM
#' @param ControlColnum control column index(s)
#' @param TreatmentColnum treatment column index(s)
#' @param p.cutoff cutoff for multi-sample comarison (p-value)
#' @param cutoff cutoff for single sample comparison (delta M)
#' @param MethylDemethyl Metylation or Demethylation analysis
#' @param version version of methylation array 450 or EPIC(850)
#' @param sampling number of sampling. if FALSE, use all data
#' @param seq_range range from CpG to be analyzed
#' @param outname output file name
#' @param nbiom_cutoff cutoff of negative binomial model test
#' @param min.score he minimum score for counting a match. Can be given as a character string containing a percentage (e.g. "85%") of the highest possible score or as a single number. 
#'
#' @import ggplot2
#' @import InfiniumDiffMetMotR
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom utils read.table
#' @importFrom stats runif na.omit
#'
#' @return a List of lists of motif positions in DMR from CpG sites and of that in random regions
#'
#' @keywords DMR Infinium motif TF
#' @export

MotPosiList <- function(infile="sel_processed_Mval.txt",
                        motifDBList,
                        ControlColnum,
                        TreatmentColnum,
                        p.cutoff = 0.001,
                        cutoff = 2,
                        MethylDemethyl="Demethyl",
                        version = "850",
                        sampling = FALSE,
                        seq_range = c(-200, 200),
                        outname = "motif_num_dist",
                        nbiom_cutoff = 0.05,
                        min.score = "90%"
){
    #library("InfiniumDiffMetMotR")
    selDataMatrix <- read.table(infile)    #Reading of M-value matrix

    DMP_IDs <- DmpId(selDataMatrix=selDataMatrix, ControlColnum = ControlColnum, TreatmentColnum = TreatmentColnum, p.cutoff=p.cutoff, cutoff= cutoff, MethylDemethyl=MethylDemethyl)

    if(!sampling == FALSE & length(DMP_IDs) >= sampling){ ##fast option: sampling of 300 DMPs
        ##fast option
        cat(paste0("Analysis will be run with ", sampling, "randomly selected DMPs."))
        cat("\n")
        DMP_IDs <- DMP_IDs[floor(runif(sampling, 1, length(DMP_IDs)))]
    }

    if(version == "450"){
        probe_annotation <- InfiniumDiffMetMotR::Methyl450anno
    } else if((version == "EPIC") || (version == "850")){
        probe_annotation <- InfiniumDiffMetMotR::EPICanno
    }
    
    target_position <- na.omit(probeID2position(probe_IDs = DMP_IDs, anno_info = probe_annotation)) #conversion of DMP IDs to position
    randomProbe_IDs <- stratSampling(target_IDs = DMP_IDs, anno_info = probe_annotation) #stratified sampling for negative control
    random_position <- na.omit(probeID2position(probe_IDs = randomProbe_IDs, anno_info = probe_annotation)) #conversion of NC probe IDs to position
    positionsList <- list("target" = target_position, "random" = random_position) #integrate DMP positoins and NC positions

    ## sequence extraction
    ## Read human hg19 genomic sequence
    library("BSgenome.Hsapiens.UCSC.hg19")
    tmp <- ls(paste("package", "BSgenome.Hsapiens.UCSC.hg19", sep=":"))
    genome <- eval(parse(text=tmp))
    cat("Retreave the sequences...\n")
    seq_range <- seq_range #range from the CpG position to be extracted
    sequences <- lapply(positionsList , function(x){seqExtract(positions = x, genome = genome, seq_range)})

    ## make a templrary output directory
    tempDir <- paste(outname, "_temp", sep="")
    dir.create(tempDir)
    ## writing the sequences to splitted files
    seqs <- sequences$target
    target_all_filenames <- writeSplitSeq (seqs=seqs, split_num = 2500, tempDir=tempDir, output_file="target" )
    seqs <- sequences$random
    random_all_filenames <- writeSplitSeq (seqs=seqs, split_num = 2500, tempDir=tempDir, output_file="random" )
    rm(sequences)
    rm(seqs)
    invisible(replicate(3, gc()))

    ## Motif search
    cat(paste("motif serch: Total ", length(motifDBList), " motifs", sep=""))
    cat("\n\tTarget regions...\n")
    ## ((multi-fasta file(multi-seqs) x motif) x [motif number])) x [multi-fasta file number]
    target_positionsList <- splitSeqMotDist(filenames=target_all_filenames,  motif_list=motifDBList, min.score = min.score)
    file.remove(target_all_filenames)
    cat("\tbackground regions...\n")
    random_positionsList <- splitSeqMotDist(filenames=random_all_filenames,  motif_list=motifDBList, min.score = min.score)
    file.remove(random_all_filenames)
    gc()
    file.remove(tempDir)
    motif_positionsList <- list(target=target_positionsList, random=random_positionsList)
    return(motif_positionsList)
}