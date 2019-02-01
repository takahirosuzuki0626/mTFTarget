#' Extraction of DMR associated enhancer
#'
#' This function extract a Enhancer horboring DMR as a GRanges object
#' @param sig_targets_gr GRanges object
#'
#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits
#'
#' @return a GRanges object of DMR horboring enhancers
#'
#' @keywords DMR Enhancer GRanges
#' @export

DMREnhancer <- function(sig_targets_gr){
#library("GenomecRanges")

## read GeneHancer data
#gh_bed <- read.table("/osc-fs_home/t-suzuki/GeneHancer/GeneHancer_version_4_4_2_hg19.bed", header=F, stringsAsFactors=F)
gh_bed <- gh_bed
gh_bed <- gh_bed[grep("gl000",gh_bed[,1], invert=T),]    #removing GL000 seq
gh_gr <- bed2granges(gh_bed)

sig_enhancer_gr <- gh_gr[subjectHits(findOverlaps(sig_targets_gr, gh_gr))]    # intersect between enhancer and DMR

return(sig_enhancer_gr)
}