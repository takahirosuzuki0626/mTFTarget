#' gff3 to GRanges
#'
#' This function loads gff3 format file and store it as a GRanges object
#' @param gff3file location of a gff3 file
#'
#' @import GenomicRanges
#' @importFrom utils read.table
#'
#' @return a GRange object
#' 
#' @keywords gff3 GRanges
#' @export

gff32granges <- function(gff3file){
    df <- read.table(gff3file, header=F, stringsAsFactors=F)
    if(length(df) != 9){
        stop("File is not gff3 format")
    }

    header <- c('chr','source','typr','start','end','score','strand','phase','attributes')
    names(df) <- header[1:length(names(df))]
    if('strand' %in% colnames(df)){
        df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
    }

    #library("GenomicRanges")
    gr <- with(df, GRanges(chr, IRanges(start, end), id=attributes, score=score, strand=strand))
    return(gr)
}
