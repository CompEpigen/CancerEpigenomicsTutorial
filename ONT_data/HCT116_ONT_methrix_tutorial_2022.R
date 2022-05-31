#
# Analysis of MinION sequencing of HCT116 WT and DKO samples using methrix
#

#### Prerequisites: package installations

install.packages("tidyverse")
BiocManager::install("annotatr")
BiocManager::install("methrix")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("org.Hs.eg.db")

library(tidyverse)
library(annotatr)
library(methrix)


BASE_DIR<-"~/Documents/teaching/DKFZ/MajorCanBiol/B370_practical/2022/"

DATA_DIR<-file.path(BASE_DIR, "data", "methylation")

OUT_DIR<-file.path(BASE_DIR, "output")

## FLONGLE_RES_DIR<-"/storage/data2/practicals_29052021_HCT116_nanopore_sequencing/original_nanopype_data//practicals_290523021/ANALYSIS/"


bdg_files <- list.files(path = DATA_DIR, pattern = "*bg$", full.names = TRUE)
hg19_cpgs <- suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg19"))


sample_info<-data.frame("Sample_Name"=gsub("-", ".", gsub("_chr1.7_metCpG.bg", "", basename(bdg_files))), "Cell_line"="HCT116", "Sample_Group"=c("WT", "DKO"))
rownames(sample_info)<-sample_info$Sample_Name

meth <- methrix::read_bedgraphs(
        files = bdg_files, 
        ref_cpgs = hg19_cpgs, 
        coldata=sample_info,
        chr_idx = 1, start_idx = 2, strand_idx=NULL, beta_idx=NULL, cov_idx = 5, M_idx=6,
        stranded = FALSE, h5=FALSE, 
        h5_dir=file.path(BASE_DIR, "methrix_hdfs/"))


####################################################################################################################################################
## Step 1: QC and filtering
####################################################################################################################################################


#meth<-load_HDF5_methrix(file.path(OUT_DIR, "methrix_objects", "methrix_15_08_19"))

# 1.1: QC report
methrix::methrix_report(meth = meth, output_dir = file.path(OUT_DIR, "methrix_report_redone"))


# 1.2: mask low coverage CpGs
meth <- methrix::mask_methrix(m = meth, low_count=5, high_quantile = 0.99)


# 1.3: filter out uncovered CpGs
meth <- methrix::remove_uncovered(m = meth)

save_HDF5_methrix(meth, dir=file.path(BASE_DIR, "methrix_objects", "methrix_filtered_06_05_2021"))


# 1.4: QC report on coverage filtered data 
methrix::methrix_report(meth = meth, recal_stats=TRUE, output_dir = file.path(OUT_DIR, "methrix_report_filtered_redone"))


####################################################################################################################################################
## Step 2: Exploratory analysis
####################################################################################################################################################


# 2.1 Global summaries
### extract interesting regions from Annotatr built-in annotation tables

data('annotations', package = 'annotatr')

ASSEMBLY<-"hg19"

types<-builtin_annotations()
types<-types[grep(ASSEMBLY, types)][c(1:15, 124)]

assembly_annotations<-build_annotations(genome=ASSEMBLY, types)

granges_list<-lapply(types, function (tt) assembly_annotations[assembly_annotations$type==tt])
names(granges_list)<-types

#granges_list<-readRDS("/C010-datasets/External/MethylomeEngram/WGBS/annotatr_regions.RDS")

### summarize methrix objects 

summary_data_tables<-list()
for(i in seq(granges_list)){
    summary_data_tables[[i]]<-get_region_summary(m = meth, regions = granges_list[[i]], type = "M", how = "mean")
}
names(summary_data_tables)<-names(granges_list)

## convert everything to a large "long" data frame
summary_dfs<-lapply(names(summary_data_tables), function(nn) {df<-as.data.frame(summary_data_tables[[nn]]); df$Region_type<-gsub("hg19_", "", nn); df})
names(summary_dfs)<-names(summary_data_tables)
lapply(summary_dfs, dim)

sum_df<-do.call("rbind", summary_dfs)
saveRDS(sum_df, file=file.path(OUT_DIR, "ont_region_summaries_updated.RDS"))


#### Plot summarized methylation value distributions
library(tidyverse)

PLOT_DIR<-file.path(BASE_DIR, "plots")

sum_df<-readRDS(file=file.path(OUT_DIR, "ont_region_summaries.RDS"))
### or on the Orchestra machine
download.file("ftp://ftp.dkfz-heidelberg.de/outgoing/ont_region_summaries_b370/ont_region_summaries_updated.RDS", "~/ont_region_summaries_updated.RDS")
sum_df<-readRDS("~/ont_region_summaries_updated.RDS")


## convert the data frame to "long" format

sum_df_long<-gather(sum_df, Sample, Mean_methylation, HCT116_DKO:HCT116_WT)


## plot distribution violins
p<-ggplot(sum_df_long, aes(x=Sample, y=Mean_methylation)) + geom_violin(aes(fill=Sample))

p<-p + facet_wrap("Region_type", ncol=5)

pdf(file.path(PLOT_DIR, "violins_by_region.pdf"), width=12, height=6)
print(p)
dev.off()

####################################################################################################################################################
## Step 3. Differential methylation analysis
####################################################################################################################################################
                     
##Try to find any interesting differences

sum_df$Delta<-sum_df$HCT116_DKO-sum_df$HCT116_WT

p<-ggplot(sum_df, aes(x=1, y=Delta)) + geom_violin()

p<-p + facet_wrap("Region", ncol=5)

pdf(file.path(PLOT_DIR, "violins_delta_by_region.pdf"), width=12, height=6)
print(p)
dev.off()

#### remove all regions where delta is NA
sum_df_nna<-sum_df[!is.na(sum_df$Delta),]

rgn_types<-unique(sum_df_nna$Region)
dmrs<-as.list(rep(list(NULL), length(rgn_types)))
names(dmrs)<-rgn_types


## call DMRs based on a simple difference threshold
dmr_counts<-list()


for(rgn in rgn_types){
    sum_df_rgn<-sum_df_nna[sum_df_nna$Region==rgn,]
    dmrs[[rgn]]$hyper<-sum_df_rgn[sum_df_rgn$Delta  > 0.2,]
    dmrs[[rgn]]$hypo<-sum_df_rgn[sum_df_rgn$Delta < -0.2,]

    dmr_counts[[length(dmr_counts)+1]]<-data.frame(Region_type=rgn, Direction="hyper", Count=nrow(dmrs[[rgn]]$hyper))
    dmr_counts[[length(dmr_counts)+1]]<-data.frame(Region_type=rgn, Direction="hypo", Count=nrow(dmrs[[rgn]]$hypo))
}

dmr_counts<-do.call("rbind", dmr_counts)

## plot DMR counts

p<-ggplot(dmr_counts, aes(x=Direction, y=Count)) + geom_bar(stat="identity", aes(fill=Direction))

p<-p + facet_wrap("Region_type", ncol=5)

pdf(file.path(PLOT_DIR, "DMR_counts_region.pdf"), width=12, height=6)
print(p)
dev.off()

#library(GenomicRanges)
#dmrs_annot<-as.list(rep(list(NULL), length(rgn_types)))
#for(rgn in rgn_types){
#    for(dir in c("hypo", "hyper")){
#        dmrs_gr<-makeGRangesFromDataFrame(dmrs[[rgn]][[dir]])
#        
#    }
#}

### explore the regions using IGV and GREAT (http://great.stanford.edu/public/html/)

head(dmrs[["genes_promoters"]][["hypo"]])


### END
