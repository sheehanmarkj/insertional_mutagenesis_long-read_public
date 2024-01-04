#!/usr/bin/env Rscript

## expects up to three arguments: SV.tsv, input.bam, and vector_elems.txt (last is list of fastq names for vector elements)
## given bam and SV.tsv, pull reads with identified insertion sites
## construct bed from aligned portions, where alignment is to one of the implicated contigs from SV res
## from this, create a string describing the observed read:
## for unaligned portions, [nbase:S]
## for aligned portions, [nbase_query:chrom:start_ref-end_ref]
## output a tsv with read\tIS_read
### Update from original:
### Original was made for pbmm2 alignment
### This one should be more tolerant of different styles of cigar string
### Also needs to adjust qwidth and isize (hardclipping breaks qwidth; isize missing from bwa mem)

## argument handling:
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("At least two arguments must be supplied (alignment file).bam, (vector elements).txt; optionally, additional argument (SV reads).tsv", call.=FALSE)
} else if (length(args)==3) {
	SV_path <- file.path(args[grep("\\.tsv$",args)])
	bam_path <- file.path(args[grep("\\.bam$",args)])
	vector_elems <- read.table(args[grep("\\.txt$",args)], header=FALSE)$V1
	SV_input <- TRUE
} else if (length(args)==2) {
	bam_path <- file.path(args[grep("\\.bam$",args)])
        vector_elems <- read.table(args[grep("\\.txt$",args)], header=FALSE)$V1
        SV_input <- FALSE
}


library(magrittr)

bam <- Rsamtools::scanBam(bam_path) 
bam <- data.frame(bam[[1]])


if ( SV_input ) {
	SVs <- read.table(SV_path, header = TRUE)
	if (dim(SVs)[1] == 0) {
		write.table(t(c("reads", "qwidths", "IS_reads", "IS_reads_simplified")), "", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
		quit()
	}
	reads <- unique(SVs$reads); reads <- as.character(unlist(reads))
} else {
	reads <- unique(bam$qname); reads <- as.character(unlist(reads))
}


recalculate_qwidth_isize=(sum(grepl("H",bam$cigar)) > 0) ## recalculate these metrics per alignment, if hardclipping is found in the bam - e.g. bwa aligner

cigar_getMbed <- function(bam_line) {
	if ( bam_line$strand == "-" ) {
		bam_line$cigar <- paste0(rev(strsplit(bam_line$cigar,"(?<=[X|=|M|D|I|S|H])", perl=TRUE)[[1]]), collapse="")
	}
	softclips <- as.numeric(gsub("[S|H]","",stringr::str_extract_all(bam_line$cigar,"[0-9]*[S|H]")[[1]]))
	if ( length(softclips) != 0 ) {
		if( length(softclips)==1 ) {
			if( substr(gsub("[0-9]*","",bam_line$cigar),1,1) %in% c("H","S") ) {
				softclips <- c(softclips, 0)
			} else {
				softclips <- c(0, softclips)
			}
		}
	} else {
		softclips <- c(0,0)
	}
	bed <- c(as.character(bam_line$rname), softclips[1]+1, bam_line$qwidth-softclips[2])
	if ( bam_line$strand == "+" ) { 
		bed <- c(bed,bam_line$pos, bam_line$pos+bam_line$isize)
	} else {
		bed <- c(bed,bam_line$pos+bam_line$isize, bam_line$pos)
	}	
	return(bed)
}

qwidths <- numeric()
IS_reads <- character()
IS_reads_simplified <- character()
IS_reads_mapq <- character()
for( read in reads ) {
	bam_sub <- bam %>% dplyr::filter(qname==read)
	if ( recalculate_qwidth_isize ) {
		bam_sub <- bam_sub %>% dplyr::mutate(qwidth = max(qwidth)) %>% # primary alignment won't have hardclipped bases, and qwidth will be accurate
			dplyr::rowwise() %>% dplyr::mutate(isize = sum(as.numeric(gsub("M|D","",grep("M|D",strsplit(cigar,"(?<=[X|=|M|D|I|S|H])", perl=TRUE)[[1]], value = TRUE))))) %>%
			dplyr::ungroup()
	}
	if ( SV_input ) {
		SVs_sub <- SVs[grep(read,SVs$reads),]
		bam_sub <- bam_sub[bam_sub$rname %in% c(SVs_sub$source, SVs_sub$target),]
	}
	bed <- data.frame(rname=NULL,start=NULL,end=NULL,ref.start=NULL,ref.end=NULL) 
	for ( i in 1:nrow(bam_sub) ) {
		bed <- rbind(bed, c(cigar_getMbed(bam_sub[i,]), bam_sub$mapq[i]))
	}
	colnames(bed) <- c("rname","start","end","ref.start","ref.end","mapq")
	for ( i in 2:ncol(bed) ) {
		bed[,i] <- as.numeric(as.character(bed[,i]))
	}
	bed <- bed[order(bed$start),]

	qwidths <- c(qwidths, bam_sub$qwidth[1])

	# construct read
	IS_read <- character()
	if ( bed$start[1] != 1 ) { IS_read <- paste0("[",bed$start[1]-1,"S]") }
	IS_read <- paste0(IS_read, paste0("[",bed$end[1]-bed$start[1],":",bed$rname[1],":",bed$ref.start[1],"-",bed$ref.end[1],"]"), collapse="")
	i=1
	if ( nrow(bed) > 1 ) {
		for ( i in 2:nrow(bed) ) {
			if ( bed$start[i] > (bed$end[i-1]+1) ) { 
				IS_read <- paste0(IS_read, paste0("[",bed$start[i]-bed$end[i-1]-1,"S]"), collapse="")
			}
			IS_read <- paste0(IS_read, paste0("[",bed$end[i]-bed$start[i],":",bed$rname[i],":",bed$ref.start[i],"-",bed$ref.end[i],"]"), collapse="")
		}
	}
	if ( bed$end[i] < bam_sub$qwidth[1] ) {
		IS_read <- paste0(IS_read, "[",bam_sub$qwidth[1]-bed$end[i],"S]")
	}
	IS_reads <- c(IS_reads, IS_read)
	
	if ( SV_input) {
		# construct read, simplified (host, chrom)
		if ( bed$rname[1] %in% vector_elems ) {
			bed$rname[1] <- "vector"
		} else {
			bed$rname[1] <- "host"
		}
		temp_line <- bed[1,]
		IS_read_simplified <- character()
		for ( i in 2:nrow(bed) ) {
			if ( bed$rname[i] %in% vector_elems ) {
				bed$rname[i] <- "vector"
			} else {
				bed$rname[i] <- "host"
			}
			if ( bed$rname[i] == temp_line$rname ) {
				temp_line$end <- bed$end[i]
			} else {
				IS_read_simplified <- paste0(IS_read_simplified, paste0("[",temp_line$rname,":",temp_line$start,"-",temp_line$end,"]"), collapse="")
				temp_line <- bed[i,]
			}
		
		}
		if (temp_line$rname == bed$rname[i]) {
			IS_read_simplified <- paste0(IS_read_simplified, paste0("[",temp_line$rname,":",temp_line$start,"-",temp_line$end,"]"), collapse="")
		}
		IS_reads_simplified <- c(IS_reads_simplified, IS_read_simplified)
	
}
	# get component mapqs
	IS_read_mapq <- paste0(bed$mapq, collapse=",")
	IS_reads_mapq <- c(IS_reads_mapq, IS_read_mapq)
	
}

if ( SV_input ) {
	output <- cbind(reads, qwidths, IS_reads, IS_reads_simplified, IS_reads_mapq)
} else {
	output <- cbind(reads,qwidths, IS_reads, IS_reads_mapq)
}

write.table(output, "", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

