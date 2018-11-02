# the script detects all *.gff(3) files in ./genomes/, simplify them, and save simplified files as ./genomes/*simGff(3).

# parseF9() <- simplifyGFF()
parseF9 <- function(F9) {
	trimTail <- function(x) gsub('[,|;].*','',x,perl=T)
	getVal <- function(key='.*gene_biotype=',F9) {
		y=gsub(key,'',F9,perl=T)
		y=trimTail(y)
		y[grepl('=',y)] <- NA
		y
	}
	gene_type <- getVal('.*gene_biotype=',F9)#gsub('.*gene_biotype=','',F9,perl=T)
	gene_type[!is.na(gene_type) & match(tolower(gene_type),tolower(c('lncRNA','lnc_RNA','protein_coding','pseudogene')),nomatch=0)<1] <- 'nonInterested'	
	name <- getVal('.*gene=',F9)#gsub('.*Name=','',F9,perl=T)
	cbind(geneType=gene_type,name=name)
}
simplifyGFF <- function(gffFile,chrFile='chr2NC.txt') {
	chrMap <- read.delim(chrFile,header=F,as.is=T)
	chrMap <- unique(cbind(chr=paste0(chrMap[,1],chrMap[,2]),gsub('\\.\\d+','',chrMap[,3],perl=T))) #NC=chrMap[,3])
	gff <- read.delim(gffFile,header=F,comment.char='#')
	colnames(gff)[c(1,3,4,5,9)] <- c('chrNC','fType','chrStart','chrEnd','F9')
	gff <- subset(gff,match(tolower(fType),tolower(c('lncRNA','lnc_RNA','pseudogene','gene')),nomatch=0)>0)[,c(1,4,5,3,9)] # Drop a lot of rows by filtering condition
	gff[,1] <- gsub('\\.\\d+','',gff[,1],perl=T)
	C2 <- parseF9(gff$F9)
	gff <- cbind(gff[,c(1,2,3)],C2,gff[,4:5])
	simGff <- merge(chrMap,gff,by.x=2,by.y=1,all.y=T)
	simGff <- data.frame(simGff[,1:4],geneType=as.character(simGff[,5]),name=simGff[,6],fType=as.character(simGff[,7]),stringsAsFactors=FALSE)
	if (any(is.na(simGff[,2]))) {
		cat(sum(is.na(simGff[,2])),' rows out of total ',nrow(simGff),'are DROPPED due to irregular chr name\n')
		simGff <- simGff[!is.na(simGff[,2]),] # Drop some rows by irregular chr.
	}
	if (any(is.na(simGff[,5]))) {
		if ( any(  is.na(simGff[,5])&!is.na(simGff[,7]) )) {
			n.replace <- sum(is.na(simGff[,5])&!is.na(simGff[,7]))
			simGff[is.na(simGff[,5])&!is.na(simGff[,7]),5] <- simGff[is.na(simGff[,5])&!is.na(simGff[,7]),7]
			cat(n.replace,' rows out of total ', nrow(simGff), 'have NA geneType replaced with fType\n')
		} else {
			stop('No gene_biotype and no F3 fType!!!')
		}
	}
	cat(sum(simGff[,5]=='gene'),'rows out of total',nrow(simGff),'are imputed with DEFAULT PROTEIN_CODING\n')
	simGff[tolower(simGff[,5])=='gene',5] <- 'protein_coding'
	cat(sum(simGff[,5]=='nonInterested'),'rows out of total',nrow(simGff),'are DROPPED due to uninterested gene_biotype\n')
	simGff <- simGff[simGff[,5]!='nonInterested',]
	simGff[,5] <- gsub('lnc_RNA','lncRNA',simGff[,5])
	outFile <- gsub('.gz','',gsub('gff','simGff',gffFile),fixed=T)
	cat('FINALLY,',nrow(simGff),'rows remained in',outFile,'\n')
	simGff <- unique(simGff)
	simGff <- simGff[order(as.character(gsub('chr','',simGff[,2])),simGff[,3]),]
	write.table(simGff[,c(2:6)],outFile,row.names=F,sep='\t',quote=F)
	simGff
}
gffFiles <- dir('./genomes',pattern='gff',full.names=T) # UPDATE 10/5/2018
simGFFs <- vector('list',length(gffFiles))
names(simGFFs) <- gsub('.*/','',gffFiles) # UPDATE 10/5/2018
for(i in gffFiles) {
	cat(date(),'Simplifying GTF file ',i,'\n')
	simGFFs[[i]] <- simplifyGFF(i)
}
#save.image('simplifyGFF.RData')

