#' position to bed3
#'
#' This function converts position informatin to bed3 format, ranging the indicated range centering the coordination
#' @param positions a dataframe or matrix of genomic coordinate
#' @param seq_range vector indicating lower and upper position from the position
#'
#' @return a dataframe of bed3-like format
#' 
#' @keywords DMR Promoter GRanges
#' @export

position2bed3 <- function (positions, seq_range){
    positions <- positions[(positions[,1] != ""),]
    chr <- positions[,1]
    start <- as.numeric(positions[,2])+(seq_range[1])	# lower range 
    end <- as.numeric(positions[,2])+(seq_range[2])		# upper range
    bed3 <- data.frame(chr=chr, start=start, end=end)
    return(bed3)
}
