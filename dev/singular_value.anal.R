# Conifer -- pick number of vectors to exclude
library(ggplot2)

dir <- '/home/data/Conifer-Coseq21/no_dup_bed'
file <- paste(dir,"singular_values.txt",sep='/')
name <- 'Coseq21'

n_vecs_display <- 20

#read data
a <- read.table(file,sep="\t",row.names=1,header=F)
names(a) <- paste('V',1:ncol(a),sep='')
a_median <- apply(a,2,median)
a_dat <- data.frame(vector=1:length(a_median),sv=a_median)

#Calculate derivatives
a_d1 <- do.call('rbind',lapply(1:nrow(a), function(i) {
			a[i,2:ncol(a)] - a[i,1:(ncol(a)-1)]
	}) )
a_d1_median <- apply(a_d1,2,median)
a_d1_dat <- data.frame(vector=(0.5 + 1:length(a_d1_median)),sv=a_d1_median)
	
a_d2 <- do.call('rbind',lapply(1:nrow(a_d1), function(i) {
			a_d1[i,2:ncol(a_d1)] - a_d1[i,1:(ncol(a_d1)-1)]
	}) )
a_d2_median <- apply(a_d2,2,median)
a_d2_dat <- data.frame(vector=(1 + 1:length(a_d2_median)),sv=a_d2_median)

#plot values
pdf(paste(name,'singular_value_plots.pdf',sep='.'),8,8)
p <- ggplot(a_dat[1:n_vecs_display,], aes(x=vector,y=sv) )
p <- p + geom_line(colour='black') + geom_point(colour='black')
p <- p + geom_line(data=a_d1_dat[1:n_vecs_display,], colour='blue') + geom_point(data=a_d1_dat[1:n_vecs_display,], colour='blue')
p <- p + geom_line(data=a_d2_dat[1:n_vecs_display,], colour='red') + geom_point(data=a_d2_dat[1:n_vecs_display,], colour='red')
print(p)
dev.off()

# Percent change 
a_dat$percent_change <- c(NA,- a_d1_dat[,2]/a_dat[2:nrow(a_dat),2]*100)
pdf(paste(name,'singular_value_percent.pdf',sep='.'),8,8)
p <- ggplot(a_dat[1:n_vecs_display,], aes(x=vector,y=percent_change) )
p <- p + geom_line(colour='black') + geom_point(colour='black')
print(p)
dev.off()

change_val <- max(which(a_dat$percent_change > 5))
print(change_val)

#Plan to trial excluding (change_val-1):(change_val+1)

