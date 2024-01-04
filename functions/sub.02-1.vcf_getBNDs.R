#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).vcf", call.=FALSE)
} else if (length(args)>0) {
	data <- suppressWarnings(VariantAnnotation::readVcf(args[1]))
}

if (dim(data)[1] == 0) {
	write.table(t(c("source","source.loc","target","target.loc","reads")), file="", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	quit()
}

if ( length(args) == 2 ) {
	vector_elems <- read.table(args[2], header=FALSE)$V1
} else {
	vector_elems <- unique(as.character(levels(GenomeInfoDb::seqnames(data))))
	vector_elems <- vector_elems[!vector_elems %in% c(seq(1,99,1), "X", "Y", "MT")]
}


data <- data[data@info$SVTYPE=="BND",]
data@info$alt <- gsub(":.*","",gsub("\\[N|\\]N|N\\]|N\\[|\\[|\\]","",data@fixed$ALT))
data@info$alt.range = gsub(".*:","",gsub("\\[N|\\]N|N\\]|N\\[|\\[|\\]","",data@fixed$ALT))
data <- subset(data, subset = (as.logical(as.character(GenomeInfoDb::seqnames(data)) %in% vector_elems) & as.logical(!data@info$alt %in% vector_elems)) | (as.logical(!as.character(GenomeInfoDb::seqnames(data)) %in% vector_elems) & as.logical(data@info$alt %in% vector_elems)))

df <- data.frame(source=as.character(data@rowRanges@seqnames), source.loc=as.numeric(data@rowRanges@ranges@start), target=as.character(data@info$alt), 
	target.loc=as.numeric(data@info$alt.range), reads=sapply(as.list(data@info$RNAMES),function(lst) paste(lst,collapse=",")))

## now we need to swap out any cases where source is vector and target is host
if(sum(df$source %in% vector_elems)>0) {
	df[df$source %in% vector_elems,] <- df[df$source %in% vector_elems, c(3,4,1,2,5)]
}

## now separate comma-delimited names into distinct rows
df <- tidyr::separate_rows(df, "reads", sep=",")

## and write new files
write.table(df, file="", sep="\t", quote=FALSE, row.names=FALSE)
#write.table(df[,c(1,2,5)], file=file.path(dirname(args[1]), gsub("\\..*","_host.tsv",basename(args[1]))), sep="\t")
#write.table(df[,c(3,4,5)], file=file.path(dirname(args[1]), gsub("\\..*","_vir.tsv",basename(args[1]))), sep="\t")

