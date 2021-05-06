##### preliminaries
library(RnBeads)

DATA_DIR<-"~/data"
#SAMPLE_SHEET_FILE<-"/ngs_share/scratch/rnbeads_tutorial/Sample_sheet_complete_clean.csv"
OUT_DIR<-"/ngs_share/scratch/rnbeads_tutorial/analysis"
dir.create(OUT_DIR)
#GEO_ACCESSION<-"GSE29290"
#GEO_URL<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29290"

###### load data from GEO

rnb.options(disk.dump.big.matrices = FALSE)
##rnb.set<-rnb.read.geo("GSE29290")

## or if the connection does not work
rnb.set<-rnb.read.geo("/ngs_share/scratch/rnbeads_tutorial/GSE29290_series_matrix.txt.gz")

## leave out non-HCT116 samples
non_hct116<-which(!pheno(rnb.set)[["source_name_ch1"]] %in% "Colon")
rnb.set<-remove.samples(rnb.set, non_hct116)

## clean up the phenotypic table

head(pheno(rnb.set))
rnb.set@pheno<-pheno(rnb.set)[,sapply(pheno(rnb.set), function(x) length(unique(x)))>1]
rnb.set@pheno<-pheno(rnb.set)[,c(1:3)]
colnames(rnb.set@pheno)<-c("Sample_Title", "GEO_accession", "Group")
rnb.set@pheno$Replicate<-sapply(strsplit(pheno(rnb.set)[["Sample_Title"]], split="\\(|\\)"), "[", 2)
rnb.set@pheno$Sample_Group<-rep("WT", nrow(pheno(rnb.set)))
rnb.set@pheno$Sample_Group[grep("DKO", pheno(rnb.set)[["Sample_Title"]])]<-"DKO"
rnb.set@pheno$Sample_ID<-paste(sep="_",  pheno(rnb.set)[["Sample_Group"]], pheno(rnb.set)[["Replicate"]])
rnb.set@pheno<-pheno(rnb.set)[c("GEO_accession", "Sample_ID", "Sample_Group", "Replicate")]
pheno(rnb.set)


####### full analysis
#options(fftempdir="/mnt/data-ssd/tmp/")
#options(fftempdir="/ngs_share/tmp/")


rnb.options(analysis.name = "Re-analysis of HCT116 450k data from Dedeurwaerder et al.")
#rnb.options(logging = TRUE)
rnb.options(email = "p.lutsik@dkfz.de")
#rnb.options(assembly = "hg19")
#rnb.options(analyze.sites = TRUE)
#rnb.options(region.types = NULL)
#rnb.options(region.aggregation = "mean")
#rnb.options(region.subsegments = 0)
#rnb.options(region.subsegments.types = NULL)
rnb.options(identifiers.column = "Sample_ID")
#rnb.options(points.category = c(16,17,3,15,4,7,8))
#rnb.options(colors.category = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"))
#rnb.options(colors.gradient = c("#132B43","#56B1F7"))
#rnb.options(min.group.size = 2)
#rnb.options(max.group.count = NULL)
#rnb.options(gz.large.files = FALSE)
#rnb.options(strand.specific = FALSE)
#rnb.options(replicate.id.column = NULL)
#rnb.options(import = TRUE)
#rnb.options(import.default.data.type = "idat.dir")
#rnb.options(import.table.separator = ";")
#rnb.options(import.bed.style = "BisSNP")
#rnb.options(import.bed.columns=1:5)
#rnb.options(import.bed.frame.shift = 1)
#rnb.options(import.bed.test = TRUE)
#rnb.options(import.bed.test.only = FALSE)
#rnb.options(import.idat.chunk.size = NULL)
#rnb.options(preprocessing = TRUE)
#rnb.options(qc = TRUE)
#rnb.options(qc.boxplots = TRUE)
#rnb.options(qc.barplots = TRUE)
#rnb.options(qc.negative.boxplot = TRUE)
#rnb.options(qc.snp.heatmap = TRUE)
#rnb.options(qc.snp.distances = TRUE)
#rnb.options(qc.snp.boxplot = FALSE)
#rnb.options(qc.snp.barplot = FALSE)
#rnb.options(qc.coverage.plots = FALSE)
#rnb.options(qc.coverage.threshold.plot=1:10)
#rnb.options(qc.coverage.histograms = TRUE)
#rnb.options(qc.coverage.violins = TRUE)
#rnb.options(qc.sample.batch.size = 500)
#rnb.options(normalization = NULL)
rnb.options(normalization.method = "bmiq")
rnb.options(normalization.background.method = "none")
#rnb.options(normalization.plot.shifts = TRUE)
#rnb.options(filtering.whitelist = NULL)
#rnb.options(filtering.blacklist = NULL)
#rnb.options(filtering.context.removal = c("CC","CAG","CAH","CTG","CTH","Other"))
#rnb.options(filtering.snp = "3")
#rnb.options(filtering.sex.chromosomes.removal = FALSE)
#rnb.options(filtering.missing.value.quantile = 1)
rnb.options(filtering.coverage.threshold = 3)
#rnb.options(filtering.low.coverage.masking = FALSE)
#rnb.options(filtering.high.coverage.outliers = FALSE)
#rnb.options(filtering.greedycut = TRUE)
#rnb.options(filtering.greedycut.pvalue.threshold = 0.05)
#rnb.options(filtering.greedycut.rc.ties = "row")
#rnb.options(filtering.deviation.threshold = 0)
#rnb.options(inference = TRUE)
#rnb.options(inference.targets.sva= c("Sample_Type","Replicate","Donor","Gender","Age","FEV1","FEV1_FVC","BMI"))
#rnb.options(inference.reference.methylome.column = "cell.type")
#rnb.options(inference.max.cell.type.markers = 50000)
#rnb.options(inference.top.cell.type.markers = 500)
#rnb.options(inference.sva.num.method = "leek")
#rnb.options(exploratory = TRUE)
#rnb.options(exploratory.columns = NULL)
#rnb.options(exploratory.top.dimensions = 0)
#rnb.options(exploratory.principal.components = 8)
#rnb.options(exploratory.correlation.pvalue.threshold = 0.01)
#rnb.options(exploratory.correlation.permutations = 10000)
#rnb.options(exploratory.correlation.qc = TRUE)
#rnb.options(exploratory.beta.distribution = TRUE)
rnb.options(exploratory.intersample = FALSE)
rnb.options(exploratory.deviation.plots = FALSE)
#rnb.options(exploratory.clustering = "all")
#rnb.options(exploratory.clustering.top.sites = 1000)
#rnb.options(exploratory.region.profiles = c("genes","promoters","cpgislands"))
#rnb.options(differential = TRUE)
#rnb.options(differential.site.test.method = "limma")
#rnb.options(differential.permutations = 0)
rnb.options(differential.comparison.columns = c("Sample_Group"))
rnb.options(columns.pairing = c("Sample_Group"="Replicate"))
#rnb.options(differential.comparison.columns.all.pairwise = c("diagnosis.ad"))
#rnb.options(covariate.adjustment.columns = NULL)
#rnb.options(differential.adjustment.sva = TRUE)
#rnb.options(differential.adjustment.celltype = TRUE)
rnb.options(differential.enrichment.go = TRUE)
#rnb.options(export.to.bed = TRUE)
#rnb.options(export.to.trackhub = c("bigBed","bigWig"))
rnb.options(export.to.csv = TRUE)
rnb.options(export.to.ewasher = TRUE)
#rnb.options(export.types = "sites")
#rnb.options(colors.meth = c("#AD0021","#909090","#39278C"))
#rnb.options(colors.3.gradient = c("#832424","#FFFFFF","#3A3A98"))
#rnb.options(logging.memory = TRUE)
rnb.options(logging.disk = TRUE)
#rnb.options(logging.exit.on.error = FALSE)
#rnb.options(distribution.subsample = 1000000)
rnb.options(disk.dump.big.matrices = FALSE)
rnb.options(enforce.memory.management = FALSE)
rnb.options(enforce.destroy.disk.dumps = FALSE)


rnb.run.analysis(
        #data.source=list(DATA_DIR, file.path(DATA_DIR, SAMPLE_SHEET_FILE)), 
        #data.type="idat.dir",
        data.source=rnb.set, 
        data.type="rnb.set",
        save.rdata=TRUE,
        dir.reports=file.path(OUT_DIR, "rnb_run_pilot"))
