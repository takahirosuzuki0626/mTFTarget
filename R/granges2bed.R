#' GRanges to bed
#'
#' This function write a bed format file from a GRanges object
#' @param gr GeRanges object
#' @param file file name of output bed file
#' @param names names for bed fromat
#' @param TSS logical if TRUE use start and  start+1
#'
#' @import GenomicRanges
#' @importFrom utils write.table
#'
#' @return a data frame of bed-like format (and export as a bed file)
#'
#' @keywords GRanges bed
#' @export

granges2bed <- function(gr, file= "mybedfile.bed", names=NULL, TSS=FALSE){
    if(length(names) < length(gr) || is.null(names)){
        cat("name is NULL or less than element! Use '.' as name.\n")
        names <- rep('.', length(gr))
    }else{
        if(TSS==TRUE){
            cat("Export as TSS bed.\n")
            ends=start(gr)+1
        }else if(TSS==FALSE){
            ends=end(gr)
        }
        bed_df <- data.frame(seqnames=seqnames(gr),
                                starts=start(gr),
                                ends=ends,
                                names=names,
                                scores=c(rep(".", length(gr))),
                                strands=strand(gr))
    }

    write.table(bed_df, file=file, quote=F, sep="\t", row.names=F, col.names=F)
    return(invisible(bed_df))
}