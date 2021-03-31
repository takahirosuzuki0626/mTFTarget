mTFTarget
===================
Version: 0.5

Description: This is a R package to analyze methylation-regulating TF target.

Last Update: 2019-01-22  

Depends: R (>= 3.5.0), Biobase (>= 2.5.5)  

Author: Takahiro Suzuki  

Updated by: takahirosuzuki1980@gmail.com

Install
-------
Install ofmTFTarget from github
```r
install.packages("devtools")
require(devtools)
install_github("takahirosuzuki0626/mTFTarget")
```
---
---
**Example 1: Single methylation-regulating TF analysis**
-------
## methylation-regulating TF target identification

#### 1. Load of InfiniumDiffMetMotR
```r
library(mTFTarget)
```

#### 2. Making Motif database
```r
library(InfiniumDiffMetMotR)

## Target Motif
target_TF = "GATA6"
motifDBList <- IMAGE_PWMlist
motif_names <- target_TF
motifDBList <- motifDBList[grep(motif_names,names(motifDBList))]
```

#### 3. Set parametors
```r
## Parametors
infile="sel_processed_Mval.txt"
ControlColnum = 7
TreatmentColnum = 8
p.cutoff = 0.001
cutoff = 2
MethylDemethyl="Demethyl"
version = "850"
sampling = FALSE
outname = "Hep_00_07_demet_GATA6"
nbiom_cutoff = 0.05
seq_range = c(-200, 200)
```

#### 4. Otimization of range
```r
ranges <- seq(100, 2000, length=20)    # ranges to be teste
nsig_target_by_range <-OptRange(infile=infile,
motifDBList=motifDBList,
TreatmentColnum = TreatmentColnum,
ControlColnum = ControlColnum,
MethylDemethyl=MethylDemethyl,
version = version,
ranges=ranges,
nbiom_cutoff = nbiom_cutoff,
outname = outname)
```

#### 5. Identification of DMR
```r
motif_positionsList <- MotPosiList(infile=infile,
                                   motifDBList=motifDBList,
                                   ControlColnum = ControlColnum,
                                   TreatmentColnum = TreatmentColnum,
                                   MethylDemethyl=MethylDemethyl,
                                   version = version,
                                   seq_range = seq_range,
                                   outname = outname)
```

#### 6. Exact test for nagative binomial model
```r
target_positionsList <- motif_positionsList[[1]]
random_positionsList <- motif_positionsList[[2]]
nbinom_pval <- DMRNbinomTest(target_positionsList,random_positionsList, outname)
```

#### 7. Extraction of significant target regions
```r
sig_targets <- nbinom_pval[nbinom_pval <= nbiom_cutoff,]
if(version=="450"){
    sig_targets_posi <- na.omit(probeID2position(probe_IDs=names(sig_targets),anno_info=Methyl450anno))	#conversion of DMP IDs to position
}else if ((version=="EPIC")||(version=="850")){
    sig_targets_posi <- na.omit(probeID2position(probe_IDs=names(sig_targets), anno_info=EPICanno))	#conversion of DMP IDs to position
}
sig_targets_bed <- position2bed3(sig_targets_posi, c(0,1))
sig_targets_gr <- bed2granges(sig_targets_bed)
```

#### 8. DMR associated promoter and enhancer identification
```r
sig_gene_gr <- DMRPromoter(sig_targets_gr = sig_targets_gr)    # promoter
sig_enhancer_gr <- DMREnhancer(sig_targets_gr = sig_targets_gr)    # enhancer
```

#### 9. Drawing of methyl-TF target Network
```r
## enhancer targed gene matrix
enhancer_gene_connection <- NULL
for(i in seq(length(elementMetadata(sig_enhancer_gr)[["id"]]))){
    enhancer_gene_connection <- rbind(enhancer_gene_connection, metaDissociate(elementMetadata(sig_enhancer_gr)[["id"]][i]))
}

## extract gene name
sig_genes <- sapply(strsplit(elementMetadata(sig_gene_gr)[["id"]], ";"), function(x){gsub("gene_name=", "", x[6])})

## drawing network
DMRnetwork <- metTFNet (TF=target_TF, MethylDemethyl=MethylDemethyl, enhancer_gene_connection=enhancer_gene_connection, sig_genes=sig_genes, outname=outname)
```


---
## Functional analysis
#### 1. GREAT analysis
```r
# background data
library(dplyr)
library(InfiniumDiffMetMtR)
anno_info <- EPICanno
CHR37 <- unlist(anno_info %>% select("CHR"))
CHR37 <- paste("chr", CHR37, sep="")
CPG37 <- anno_info %>% select("MAPINFO")
bg_position <- na.omit(cbind(CHR37, CPG37))
bg_bed <- position2bed3(bg_position, c(200,200))
bg_bed_gr <- bed2granges(sig_targets_bed)

## GREAT analysis
library(rGREAT)
job <- submitGreatJob(sig_targets_gr, bg = bg_bed_gr)
go_tb = getEnrichmentTables(job)
panther_tb = getEnrichmentTables(job, ontology = "PANTHER Pathway")
```

#### 2.  GO analysis by Enrichr for promoter overlapped genes
```r
## Enrichr analysis
library(enrichR)
dbs <- listEnrichrDbs()
#### all database ####
enriched_all <- enrichr(as.vector(DMRnetwork[,2]), dbs[[1]])
#### only GO ####
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched_GO <- enrichr(as.vector(DMRnetwork[,2]), dbs)
```
For web Enrichr
```r
write.table(as.vector(DMRnetwork[,2]), file=paste0(outname, "_target_enrichr.txt"), sep="\t", quote=FALSE, row.names = FALSE,col.names = FALSE)
```


---
## comparison with expression (CAGE data)
## Option 1: ngsplot
### 1. bed export
```r
## folder generation
if(!("ngsplot" %in% list.files())){
    dir.create("ngsplot")
}
## All targets
outfile_DMRbed <- paste0("ngsplot/", outname, "_DMR.bed")
write.table(sig_targets_bed, file=outfile_DMRbed, quote=F, sep="\t", row.names=F, col.names=F)
sig_targets_gr <- bed2granges(sig_targets_bed)

## target promoters TSS
sig_genecode_genes <- sapply(strsplit(elementMetadata(sig_gene_gr)[["id"]], ";"), function(x){gsub("gene_name=", "", x[6])})
outfile_TSS <- paste0("ngsplot/", outname, "_DMR_target_TSS.bed")
granges2bed(gr=sig_gene_gr, file=outfile_TSS, names=sig_genecode_genes, TSS=TRUE)

## target enhancers
sig_enhancer_id <- sapply(strsplit(elementMetadata(sig_enhancer_gr)[["id"]], ";"), function(x){gsub("genehancer_id=", "", x[1])})
outfile_enh <- paste0("ngsplot/", outname, "_DMR_target_enhancer.bed")
granges2bed(gr=sig_enhancer_gr, file=outfile_enh, names=sig_enhancer_id, TSS=FALSE)
```

#### 2. ngsplot
In `shell`
Move to ngsplot directory `cd ngsplot`
Checking of numer of region
```bash
all_nregion=`wc -l *DMR.bed|cut -f1 -d" "`
gene_nregion=`wc -l *DMR_target_TSS.bed|cut -f1 -d" "`
enh_nregion=`wc -l *DMR_target_enhancer.bed|cut -f1 -d" "`
## heatmap height index
all_height=`echo "scale=5; $all_nregion / 3000 * 30" | bc`
gene_height=`echo "scale=5; $gene_nregion / 3000 * 30" | bc`
enh_height=`echo "scale=5; $enh_nregion / 3000 * 30" | bc`
```
plot by ngsplot
```bash
DMR=`ls *DMR.bed`
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_00_1.bam ${DMR} \"00\"">list_all.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_07_1.bam ${DMR} \"07\"">>list_all.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_14_1.bam ${DMR} \"14\"">>list_all.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_21_1.bam ${DMR} \"21\"">>list_all.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_28_1.bam ${DMR} \"28\"">>list_all.txt

ngs.plot.r -G hg19 -R bed -C list_all.txt -O Hep_00_07_demethyl_all -T CAGE_tags -L 5000 -FL 300 -RR $all_height

DMR_target_TSS=`ls *DMR_target_TSS.bed`
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_00_1.bam ${DMR_target_TSS} \"00\"">list_gene.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_07_1.bam ${DMR_target_TSS} \"07\"">>list_gene.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_14_1.bam ${DMR_target_TSS} \"14\"">>list_gene.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_21_1.bam ${DMR_target_TSS} \"21\"">>list_gene.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_28_1.bam ${DMR_target_TSS} \"28\"">>list_gene.txt

ngs.plot.r -G hg19 -R bed -C list_gene.txt -O Hep_00_07_demethyl_gene -T CAGE_tags -L 5000 -FL 300 -RR $gene_height

DMR_enhancer=`ls *DMR_target_enhancer.bed`
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_00_1.bam ${DMR_enhancer} \"00\"">list_enh.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_07_1.bam ${DMR_enhancer} \"07\"">>list_enh.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_14_1.bam ${DMR_enhancer} \"14\"">>list_enh.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_21_1.bam ${DMR_enhancer} \"21\"">>list_enh.txt
echo "/home/t-suzuki/osc-fs/CAGE_Hep_SKM_diff/bam/HEP/Hep_28_1.bam ${DMR_enhancer} \"28\"">>list_enh.txt

ngs.plot.r -G hg19 -R bed -C list_enh.txt -O Hep_00_07_demethyl_enhancer -T CAGE_tags -L 2000 -FL 5000 -RR $enh_height
```


---
## Option 2: Genomation library
```r
bam.files_full <- list.files("/osc-fs_home/t-suzuki/CAGE_Hep_SKM_diff/bam/HEP", full.names=TRUE, pattern="1.bam$")
bam.files <- sapply(strsplit(bam.files_full, "/"), function(x){x[7]})
sample_names <- sapply(strsplit(bam.files, "_1.b"), function(x){x[1]})

## all targets
sig_targets_scoreMatrix <- ngsplotr(bam.files=bam.files_full, peaks=sig_targets_gr, bin.num = 50, sample_names, range=c(-5000,5000), outname="GATA6_Hep00_07_demet_sig_targets")

## promoters
sig_gene_scoreMatrix <- ngsplotr(bam.files=bam.files_full, peaks=sig_gene_gr, bin.num = 50, sample_names, range=c(-5000,5000), outname="GATA6_Hep00_07_demet_sig_gene")

## enhancers
sig_enhancer_scoreMatrix <- ngsplotr(bam.files=bam.files_full, peaks=sig_enhancer_gr, bin.num = 50, sample_names, range=c(-5000,5000), outname="GATA6_Hep00_07_demet_sig_enhancer")
```

**Example 2: multiple methylation-regulating TF analysis**
-------
#### 1. Load of InfiniumDiffMetMotR
```r
library(mTFTarget)
```

#### 2. Making Motif database
```r
library(InfiniumDiffMetMotR)

## Target Motif
target_TF = "multi"
motifDBList <- IMAGE_PWMlist
motif_names <- c("GATA6","GATA4")
motifDBList <- motifDBList[sapply(strsplit(names(motifDBList), "_"), function(x){x[1] %in% motif_names})]
```

#### 3. Set parametors
```r
## Parametors
infile="sel_processed_Mval.txt"
ControlColnum = 7
TreatmentColnum = 8
p.cutoff = 0.001
cutoff = 2
MethylDemethyl="Demethyl"
version = "850"
sampling = FALSE
outname = "multi_test"
nbiom_cutoff = 0.05
seq_range = c(-200, 200)
```

#### 4. Identification of DMR
```r
motif_positionsList <- MotPosiList(infile=infile,
                                   motifDBList=motifDBList,
                                   ControlColnum = ControlColnum,
                                   TreatmentColnum = TreatmentColnum,
                                   MethylDemethyl=MethylDemethyl,
                                   version = version,
                                   seq_range = seq_range,
                                   outname = outname)
```

#### 5. Exact test for nagative binomial model
```r
target_positionsList <- motif_positionsList[[1]]
random_positionsList <- motif_positionsList[[2]]

eg_connection_list <- as.list(NA)
sig_genes_list <-  as.list(NA)
for (i in seq(length(target_positionsList))){
    mot_name <- names(target_positionsList)[i]
    
    target_positionsList. <- target_positionsList[i]
    random_positionsList. <- random_positionsList[i]

    nbinom_pval <- DMRNbinomTest(target_positionsList=target_positionsList.,
                                 random_positionsList=random_positionsList.,
                                outname=mot_name)

    #### 6. Extraction of significant target regions
    sig_targets <- nbinom_pval[nbinom_pval <= nbiom_cutoff,]

    if(version=="450"){
        sig_targets_posi <- na.omit(probeID2position(probe_IDs=names(sig_targets),anno_info=Methyl450anno))	#conversion of DMP IDs to position
    }else if ((version=="EPIC")||(version=="850")){
        sig_targets_posi <- na.omit(probeID2position(probe_IDs=names(sig_targets), anno_info=EPICanno))	#conversion of DMP IDs to position
    }
    sig_targets_bed <- position2bed3(sig_targets_posi, c(0,1))
    sig_targets_gr <- bed2granges(sig_targets_bed)

    #### 7. DMR associated promoter and enhancer identification
    sig_gene_gr <- DMRPromoter(sig_targets_gr = sig_targets_gr)    # promoter
    sig_enhancer_gr <- DMREnhancer(sig_targets_gr = sig_targets_gr)    # enhancer

    #### 8. Drawing of methyl-TF target Network
    ## enhancer targed gene matrix
    enhancer_gene_connection <- NULL
    for(j in seq(length(elementMetadata(sig_enhancer_gr)[["id"]]))){
        enhancer_gene_connection <- rbind(enhancer_gene_connection, metaDissociate(elementMetadata(sig_enhancer_gr)[["id"]][j]))
    }

    ## extract gene name
    sig_genes <- sapply(strsplit(elementMetadata(sig_gene_gr)[["id"]], ";"), function(x){gsub("gene_name=", "", x[6])})
    
    eg_connection_list[[i]] <- enhancer_gene_connection
    names(eg_connection_list)[i] <- mot_name
    sig_genes_list[[i]] <- sig_genes
    names(sig_genes_list)[i] <- mot_name
}
```

#### drawing network
```
DMRmultinetwork <- metTFNet.multi(MethylDemethyl="Demethyl",eg_connection_list=eg_connection_list, sig_genes_list=sig_genes_list, outname=outname)
```
