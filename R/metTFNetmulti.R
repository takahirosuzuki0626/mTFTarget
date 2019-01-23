#' Target network of multiple methyl-regulationg TF
#' 
#' This function draws a DNA methylation regulatory network conposed of methyl-regulating TF, target promoter, target enhancer, and enhancer connected genes
#' This fanction is designed for multiple methyl-regulating TF input.
#' output "TF_methyl_network.pdf" and connection dataframe

#' @param MethylDemethyl a vector. methylation or demethylation mediated by the methyl-regulating TF. "Methyl" or "Demethyl"
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
#' @keywords network enhancer promoter methylation
#' @export

metTFNet.multi <- function(MethylDemethyl="Demethyl",eg_connection_list=eg_connection_list, sig_genes_list=sig_genes_list, outname=""){
    ## select top 50 hoghly scored enhancers-gene connections
    TFs <- sapply(strsplit(names(eg_connection_list), "_"), function(x){x[1]})
    df <- NULL
    for (i in  seq(length(eg_connection_list))){
        TF <- TFs[i]
        enhancer_gene_connection <- eg_connection_list[[i]]
        sig_enhancer_gene_connection <- enhancer_gene_connection[order(enhancer_gene_connection[,3], decreasing=TRUE),]    # order the enhancer gene connection based on score
        #eg_net_df <- sig_enhancer_gene_connection[1:100,c(1,2)]    #extract top 50
        eg_net_df <- data.frame(sig_enhancer_gene_connection[sig_enhancer_gene_connection[,3] >=10, c(1,2),], type="eg")    #extract score >= 10

        ## Demethyl TF-enhancer connection
        net_enhancer <- as.vector(unique(eg_net_df[,1]))
        te_net_df <- data.frame(gh_id = rep(TF, length(net_enhancer)), connected_gene=net_enhancer, type="te")

        ## Demethyl TF-prom connection
        sig_genes <- sig_genes_list[[i]]
        tp_net_df <- data.frame(gh_id = rep(TF, length(sig_genes)), connected_gene=sig_genes, type="tp")
        ## conbine henhancer-gene and demethyl TF enhancer connections
        df <- rbind(df, tp_net_df, te_net_df, eg_net_df)
    }

    ## prepare network graph data
    #library(igraph)
    g <- graph.data.frame(df, directed=T)
    g <- simplify(g, remove.multiple = T, remove.loops = F)

    ## set node type
    ntype <- NULL
    ntype[V(g)$name %in% TFs] <- "TF"
    ntype[V(g)$name %in% df[df[,"type"]=="tp","connected_gene"]] <- "prom"
    ntype[V(g)$name %in% df[df[,"type"]=="eg","gh_id"]] <- "enh"
    ntype[V(g)$name %in% df[df[,"type"]=="eg","connected_gene"]] <- "gene"
    ntype[V(g)$name %in% df[df[,"type"]=="eg","connected_gene"] & V(g)$name %in% df[df[,"type"]=="tp","connect_gene"]] <- "pro_enh"

    V(g)$Faction <- ntype

    ## color of nodes
    color <- adjustcolor( c("darkmagenta", "seagreen3", "cyan", "green4", "green"), alpha.f=.6)
    col <- NULL
    col[V(g)$Faction == "TF"] <- color[1]
    col[V(g)$Faction == "prom"] <- color[2]
    col[V(g)$Faction == "enh"] <- color[3]
    col[V(g)$Faction == "gene"] <- color[4]
    col[V(g)$Faction == "pro_enh"] <- color[5]

    ## shape of nodes
    sha <- NULL
    sha[V(g)$Faction == "TF"] <- "circle"
    sha[V(g)$Faction == "prom"] <- "circle"
    sha[V(g)$Faction == "enh"] <- "square"
    sha[V(g)$Faction == "gene"] <- "circle"
    sha[V(g)$Faction == "pro_enh"] <- "circle"

    ## size of nodes
    siz <- NULL
    siz[V(g)$Faction == "TF"] <- 5
    siz[V(g)$Faction == "prom"] <- 3
    siz[V(g)$Faction == "enh"] <- 3
    siz[V(g)$Faction == "gene"] <- 2
    siz[V(g)$Faction == "pro_enh"] <- 4

    ## plot the network
    l <- layout.kamada.kawai(g)
    #l <- layout.reingold.tilford(g)

    num.of.v <- length(V(g))
    num.of.e <- length(E(g))
    V(g)$size  <- siz
    V(g)$color <- col
    V(g)$shape <- sha
    V(g)$label <- names(as.list(V(g)))
    V(g)$label.cex   <- rep(0.1, num.of.v)
    V(g)$label.color <- rep("black", num.of.v)
    V(g)$frame.color <- "transparent"
    E(g)$width  <- rep(0.2, num.of.e)
    E(g)$arrow.size <- rep(0.1, num.of.e)

    ## title
    title <-  paste0(outname,"_network")
    
    ##plot
    outfile <- paste0(outname, "_network.pdf")
    pdf(outfile)
    plot(g, main=title, layout=l, edge.curved=.1)

    ## legend
    legend(x=-1,
           y=-1,
           legend=c("Demethyl TF", "promoter", "enhancer", "enhancer connected gene", "promoter and enhancer connected gene"),
           col = color,
           bty = "n",
           pch=c(16,16,15,16,16),
           pt.cex = 2,
           cex = 1,
           text.col="black",
           horiz = F)
    dev.off()
    return(invisible(df))
}