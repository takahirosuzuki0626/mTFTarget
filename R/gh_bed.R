#' GeneHancer bed-like format object
#'
#' This data comes from GeneHancer database,
#' \url{https://academic.oup.com/database/article/doi/10.1093/database/bax028/3737828}. 
#' The data are integrated human enhancer data set and gene connection of each enhancer.
#' The original data is downloaded from
#' \url{https://genecardsdata.blob.core.windows.net/public/GeneHancer_version_4_4.xlsx}
#' The header was removed and stored as tab-deleted gff3-like format file.
#' The gff3-like format file was forther converted to a  bed format file, then liftover from hg38 to hg19.
#' the hg19 bed format file was loaded by R as gh_bed.
#'
#' @docType data
#' @name gh_bed
#' @usage gh_bed
#' @format A bed-like object
#' @examples
#' gh_bed
NULL