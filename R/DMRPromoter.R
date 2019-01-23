#' Extraction of DMR associated promoter
#'
#' This function extract a promoter horboring DMR
#' Gene model used in this function is Gencode_hg19(gene only)
#' Promoter is defined as 1.5kbp upstream and 500 bp downstream from TSS
#' @param sig_targets_gr a Granges object
#'
#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits
#'
#' @return a GRanges object of DMR horboring genes at their promoter
#'
#' @keywords DMR Promoter GRanges
#' @export

DMRPromoter <- function(sig_targets_gr){
    #library("GenomicRanges")
    ## Read the gencode gff3 file
    #file <- "/osc-fs_home/t-suzuki/Gencode/human/gencode.v19.gene_only.gff3"
    #gencode_gr <- gff32granges(file)
    gencode_gr <- gencode_gr
    gencode_prom_gr <- promoters(gencode_gr, upstream = 1500, downstream = 500)

    ## overlap
    sig_prom_gr <- unique(gencode_prom_gr[subjectHits(findOverlaps(sig_targets_gr,gencode_prom_gr))])
    meta_id <- elementMetadata(sig_prom_gr)[["id"]]
    sig_gencode_gr <- gencode_gr[elementMetadata(gencode_gr)[["id"]] %in% meta_id]
    return(sig_gencode_gr)
}