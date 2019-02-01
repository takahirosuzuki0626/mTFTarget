#' reads enrichment plots
#' 
#' This function draws a line plot foraverage enrichment of  ngs reads and 
#' heatmaps for ngs reads enrichment (k-means and hierarchical clusterings).
#'
#' @param bam.files a vector of bam.file locations
#' @param peaks a GRanges object for regions to be analyzed.
#' @param bin.num number of bin for score matrix
#' @param sample_names  sample names of bam files
#' @param outname name for output file
#' @param range a raneg to be analyzed
#'
#' @import GenomicRanges
#' @importFrom genomation scale ScoreMatrixList plotMeta multiHeatMatrix
#' @importFrom dendextend cutree 
#' @importFrom fastcluster hclust 
#' @importFrom stats dist kmeans
#' @importFrom graphics legend plot
#' @importFrom grDevices adjustcolor dev.off pdf
#' @importFrom GiNA brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @return a scorMatrix object and export a line plot, heatmaps(k-mean and hierarchical clusterings) as a pdf)
#'
#' @keywords ngsplot bam ChIP-seq ngs heatmap enrichment
#' @export

ngsplotr <- function(bam.files, peaks, bin.num = 100, sample_names, range=c(-5000,5000), outname="my_enrichment"){
    
    peaks <- reduce(peaks)
    peaks = resize(peaks, width = 5000, fix = "center")

    sml = ScoreMatrixList(bam.files, peaks, strand.aware=TRUE, bin.num=bin.num, type = "bam")
    names(sml) <- sample_names

    sml.scaled = scaleScoreMatrixList(sml)

    smlcolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(sml))

    outfile <- paste0(outname, "_enrichment_plots.pdf")
    pdf (outfile)
    plotMeta(sml.scaled,
            line.col=smlcolors,
            xcoords = range,
            centralTend="mean",
            dispersion="se",
            winsorize=c(0,99),
            lwd=2,
            smoothfun=function(x) stats::lowess(x, f = 1/5))
    legend("bottomright", 
        names(sml.scaled),
        lty=c(1,1), 
        lwd=c(2.5,2.5),
        col=smlcolors,
        bty = "n")

     plotMeta(sml,
            line.col=smlcolors,
            xcoords = range,
            centralTend="mean",
            dispersion="se",
            winsorize=c(0,99),
            lwd=2,
            smoothfun=function(x) stats::lowess(x, f = 1/5))
    legend("bottomright", 
        names(sml.scaled),
        lty=c(1,1), 
        lwd=c(2.5,2.5),
        col=smlcolors,
        bty = "n")


    clust_num <- 30
    multiHeatMatrix(sml.scaled,
                    xcoords = range,
                    clustfun=function(x)kmeans(x, centers=clust_num)$cluster,
                    common.scale=T,
                    cex.axis=0.5,
                    cex.lab=0.5,
                    cex.legend=0.5,
                    col = c("lightgray", "blue"))

    cl2 <- function(x) cutree(hclust(dist(x), method="complete"), k=clust_num)
    multiHeatMatrix(sml.scaled,
                    xcoords = range,
                    clustfun=cl2,
                    common.scale=T,
                    cex.axis=0.5,
                    cex.lab=0.5,
                    cex.legend=0.5,
                    col = c("lightgray", "blue"))

    multiHeatMatrix(sml,
                    xcoords = range,
                    clustfun=function(x)kmeans(x, centers=clust_num)$cluster,
                    common.scale=T,
                    cex.axis=0.5,
                    cex.lab=0.5,
                    cex.legend=0.5,
                    col = c("lightgray", "blue"),
                    winsorize=c(0,95))

    multiHeatMatrix(sml,
                    xcoords = range,
                    clustfun=cl2,
                    common.scale=T,
                    cex.axis=0.5,
                    cex.lab=0.5,
                    cex.legend=0.5,
                    col = c("lightgray", "blue"),
                    winsorize=c(0,95))
    dev.off()
    return(invisible(sml.scaled))
}
