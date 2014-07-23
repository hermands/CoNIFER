# R script
#Run CN anal to make plot

library(rhdf5)
library(sqldf)
library(RColorBrewer)
library(ggplot2); library(grid)
library(reshape2)

#Input
rpkm_file <- "/home/data/Conifer-Coseq21/no_dup_bed/analysis.hdf5"
db_file <- "/home/genetics/Genomes/UCSC/UCSC.db";

"%win%" <- function(x,y) {
			sapply(x, function(x_i) {
				x_i >= y[1] && x_i <= y[2]
			})
}
"%any_in%" <- function(x,y) {
		any(sapply(x, function(x_i) { 
			x_i >= y[1] && x_i <= y[2]
		}))
}

"%near%" <- function(x,y,dist=1000) {
	sapply(x, function(x_i) { any(abs(x_i-y) < dist) })
}
#==================
load_rpkm <- function(file=rpkm_file) 
# Load in rpkm file to 'dat', 'samples', 'probes'
{

	info <- h5ls(rpkm_file,recursive=F)
	info_names <- info$name[-which(info$name %in% c('samples','probes'))]

	probes <- h5read(rpkm_file,name='probes')
	names(probes) <- substring(names(probes),8)

	samples <- h5read(rpkm_file,name='samples')
	samples <- do.call('c',samples)[[1]]
	samples <- substring(samples,1,nchar(samples)-5)
	
	dat <- lapply(info_names, function(x) {
		datx <- h5read(rpkm_file,name=x)
		daty <- do.call('cbind',lapply(datx, function(x) { x[,2] }) )
		colnames(daty) <- samples
		rownames(daty) <- datx[[1]][,1]
		daty
	})
	names(dat) <- info_names
	
	# Add gene annotations
	db <- dbConnect(SQLite(), dbname=db_file)
	refGene <- dbReadTable(db, "refGene")[,c("name","name2","chrom","txStart","txEnd")]
	chrs <- names(probes)

	probes <- lapply(names(probes), function(chr) {
					probe_poss <- sort(unique(c(probes[[chr]]$start, probes[[chr]]$stop)))
					refGene_chr <- unique(refGene[which(refGene$chrom == chr & (refGene$txStart %near% probe_poss | refGene$txEnd %near% probe_poss) ),c('name2','txStart','txEnd')])
					refGene_chr <- do.call('rbind',lapply(unique(sort(refGene_chr$name2)), function(x) { 
										data.frame(name2=x,
											txStart = min(refGene_chr$txStart[which(refGene_chr$name2 == x)]), 
											txEnd = max(refGene_chr$txEnd[which(refGene_chr$name2 == x)])) 
										}) )
										
					if ( !is.null(refGene_chr) ) {
						genes <- do.call('cbind',lapply(1:nrow(refGene_chr), function(i) {
							probes[[chr]]$start %win% refGene_chr[i,c('txStart','txEnd')] | probes[[chr]]$stop %win% refGene_chr[i,c('txStart','txEnd')]
						}))
						probes[[chr]]$genes <- apply(genes, 1, function(x) {
							paste(refGene_chr$name2[which(x)],collapse='|')
						})
					} else { probes[[chr]]$genes <- '' }
					probes[[chr]]
	})
	names(probes) <- chrs
	
	list(dat=dat, samples=samples, probes=probes)
}
#===============================
plot_dat <- function(a,name) {
# TODO: figure out why not plotting color for majority of tables
# TODO: fix sample ordering
		par(xpd=TRUE)
		for (x in names(a$dat)) {
			pdf( paste(name,x,'image.pdf',sep='.') ,20,5)
			sample_dist <- dist(t(a$dat[[x]]))
			sample_order <- hclust(sample_dist)$order
			# Sample ordering not working -- want to order by similarity

			long_data <- melt(t(a$dat[[x]])[sample_order,])
			names(long_data) <- c('probe','sample','ZRPKM')

			#Truncate ZRPKM's
			long_data$ZRPKM[which(long_data$ZRPKM > 3)] <- 3
			long_data$ZRPKM[which(long_data$ZRPKM < -3)] <- -3

			p <- ggplot(long_data, aes(x=sample,y=probe)) + geom_tile(aes(fill=ZRPKM), colour='white') + scale_fill_gradient(low='green',high='red')
			p <- p + labs(x='',y='') + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y=element_text(size=5,face='bold') )
			p <- p + geom_text(x=-8,y=-5,label=x)

			# genes
			genes <- unique(a$probes[[x]]$genes); genes <- genes[which(genes != '')]
			if(!is.null(genes) & length(genes) > 0) {
				gene_table <- do.call('rbind',lapply(genes, function(g) { 
						data.frame(gene=g, start=min(which(a$probes[[x]]$genes == g))+0.5,end=max(which(a$probes[[x]]$genes == g))+0.5)
					}))
				gene_table$ymin <- -4
				gene_table$textx <- rowMeans(gene_table[,c('start','end')])
				p <- p + geom_text(data=gene_table[which(!grepl(pattern='\\|',x=gene_table$gene)),],mapping=aes(x=textx,y=-2,label=gene), angle=45,hjust=1,vjust=0.5,fontface=2,size=6 )
				p <- p + scale_x_continuous(breaks=unique(c(gene_table$start,gene_table$end)) )
				p <- p + theme(plot.margin= unit(c(0.5,0,1.5,0.5),"cm") )
			}
			gt <- ggplot_gtable(ggplot_build(p))
			gt$layout$clip[gt$layout$name == "panel"] <- "off"
			grid.draw(gt)
			dev.off()
		}
}

# TODO
# label chromosomes, genes (from BED-file?), ask Sheena to remove subjects [look at criteria]
# confirm that the hdf5 has the info that I think it does