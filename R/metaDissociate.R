#' extraction of metadata infomation
#'
#' This function dissociate metadata of GRanges class of gff3 files
#' @param metadata semi-colon deleted metadata

#' @return a data frame of links between enhancers and connected genes

#' @keywords gff3 metadata
#' @export

metaDissociate <- function(metadata){
    meta_elements <- unlist(strsplit(metadata, ";"))
    con_gene_index <- grep("connected_gene=", meta_elements)
    con_gene <- NULL
    score <- NULL
    for(i in con_gene_index){
        con_gene <- c(con_gene, gsub("connected_gene=", "", meta_elements[i]))
        score <- c(score, as.numeric(gsub("score=", "", meta_elements[i+1])))
    }
    gh_id <- gsub("genehancer_id=", "",meta_elements[1])
    enhancer_connection <- data.frame(gh_id = rep(gh_id, length(con_gene_index)), connected_gene=con_gene, score=score)
    return(enhancer_connection)
}