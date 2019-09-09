# Rscript binomEnrich.R joined='example.out' gn='hg38' intMode='s' arg.t=12 b=1 breaks='break123' arg.d=',' region=NULL
##### OBJECTIVE: 1.to calculate Binomial probability of feature enrichment within tracts; 2.visualize the enrichment (increased risk) in barplot & pieplot; 3.visualize overall enrichP height in background of 100-feature landscape.
##### INPUT 1: prior python output file (-o argument) 
##### INPUT 2: genome (-g argument)
##### INPUT 3: mode of overlap or intersection, singleton (s) or multiplex (m). (-I argument)
##### INPUT 4: tract modality (-t argument), or JUNCTION modality (-J argument)
##### INPUT 5: boundary (-b argument)
##### INPUT 6: hinge/junction (-j argument)
##### INPUT 7: delimiter (-d argument) 
##### INPUT 8: genomic region constraint (-r argument)
##### NOTE: Nine core assemblies take the total number of nucleotides from global variable "gnNucs". Other assemblies look for this statistics from genomes.meta.
##### NOTE: Genome-specific tract stat files must be saved as .../tracts/polytract_GN_stats.txt.
##### NOTE: Precomputed enrichment P values for nearly 100 features must be saved as a file at tracts/polytrap_P_landscape.csv.
gnNucs <- c(hg38=3088269832+16569,hg19=3095677412+16571,mm10=2730871774,mm9=2716965481,rn6=2870184193, galgal5=1230258557, danrer10=1371719383, dm6=143725995, saccer3=12157105)
# groupTri(): map the 60 tri-nucleotides to 10 broad groups as given in triG10 below (in lowercase).
groupTri <- function() {
	triG10 <- c('AAC','AAG','AAT','ACC','GAC','ACT','CAG','AGG','ATC','CGG')
	triG10.sep <- strsplit(triG10,'')
	# from1to6(): input a 3-letter vector, output 6 3-letter strings.
	from1to6 <- function(letter3) {
		shift0 <- letter3
		shift1 <- letter3[c(2,3,1)]
		shift2 <- letter3[c(3,1,2)]
		mapChart <- c('A','T','C','G')
		names(mapChart) <- c('T','A','G','C')
		lst3 <- list(shift0,shift1,shift2)
		revCmp3 <- lapply(lst3, function(x,mapChart){
					mapChart[x[c(3,2,1)]]
					},
		mapChart
		)
		res <- c(lst3,revCmp3)
		res <- sapply(res,paste,collapse='')
		res
	}
	triGroups <- lapply(triG10.sep,from1to6)
	names(triGroups) <- tolower(triG10)
	triGroups
}
tract1 <- c('A','T','C','G')
tract2 <- tolower(c('TA','GA','CA','GT','CT','GC')) # Use lower-case to denote 6 dinuc types.
tract3 <- names(groupTri())
tract12 <- c(tract1,tract2)
tract23 <- c(tract2,tract3)
tract13 <- c(tract1,tract3)
tract123 <- c(tract1,tract2,tract3)
##### enrichStats_perType(): Given the four numbers involved in the classic problem of black/white balls in urn, calculate the p value using Binom distribution.
##### INPUT nEdits_inTract: number of features overlapped with a tract. The feature can be a single-nuc event or multi-nuc event.
##### INPUT nTotTracts: total nucleotides occupupied by tracts.
##### INPUT nEdits: total nucleotides in genome occupied by the interested feature.
##### INPUT nTotal: total nucleotides in genome. A fixed number dependent on genome choice (say, HG19 or HG38).
enrichStats_perType <- function(nEdits_inTract,nTotTracts,nEdits,nTotal) {
	prob <- nEdits/nTotal
	if (nEdits_inTract!=0) {
		pEnrich <- pbinom(q=nEdits_inTract,size=nTotTracts,prob=prob,lower.tail=F)
	} else {
		pEnrich = NA
	}
  c(p=pEnrich,x=nEdits_inTract,size=nTotTracts,whiteBalls=nEdits,totalBalls=nTotal,obsRate=nEdits_inTract/nTotTracts,expRate=prob)
}
### binomEnrich(): Examine the enrichment of a genomic feature within tracts, both in the whole and in various subtypes.
##### INPUT joined: Output of tract-joined features (-o file). Each row contains four fields: chr, start, end, joining info. 
##### INPUT gn: controlled vocabulary for covered genomes. It can correspond to new genomes as long as user has defined it in generating the new genome beforehand.
##### INPUT tracts: Perform subtype enrichment analysis on only the given subtypes. {tract1,tract2,tract3,...,'break1','break2','break3'}
##### INPUT b: boundary/neighboring extension.
##### INPUT arg.d: file separator. comma or tab.
##### INPUT intMode: mode of overlap or intersection, singleton (s) or multiplex (m).
##### INPUT region: {"protein","pseudo","lncRNA","NA"}
##### INPUT headed: file has header? Python-precossed joined file always has no header.
##### NOTE: this function reads in files on tract statistics, e.g., tracts/polytract_hg38_stats.txt
##### NOTE: this function may takes value of global variable gnNucs (above).
binomEnrich <- function(joined,gn=c('HG38','HG19','MM9','MM10')[1],tracts=tract123,b=c(1,2,3,0)[1],arg.d=',',intMode=c('s','m')[1],region=NULL,headed=F) {
	#### Read in previous output file (joined)
	if (arg.d==',') {
		joined <- unique(read.csv(paste('output',joined,sep='/'),header=headed,as.is=T))
	} else {
		joined <- unique(read.delim(paste('output',joined,sep='/'),header=headed,as.is=T))
	}
	joined[,ncol(joined)] <- as.character(joined[,ncol(joined)])
	n_hgEdit <- switch(intMode,
		s=nrow(joined),
		m=sum(joined[,3]-joined[,2]+1)
	)
	###### Find total genome nuc number ####
	if ( is.element( tolower(gn),names(gnNucs) ) ) {
		nTotal=gnNucs[tolower(gn)]
	} else {
		tmp <- read.csv('./new/genomes.meta',as.is=T)[,1:2]
		gnNucs <- tmp[,2]
		names(gnNucs) <- tolower(tmp[,1])
		nTotal <- gnNucs[tolower(gn)]
	}
	###### Work on overlapped input items, find numbers or overlapped nucs #####
	if (any(joined[,ncol(joined)]!='0')) {
		joined <- joined[joined[,ncol(joined)]!='0',]
		#### Obtain statistics for each tract subtype from the stored static file ####
		statFile <- paste0('tracts/polytract_',tolower(gn),'_stats.txt')
		if (file.exists(statFile)) {
			statTbl <- read.delim(statFile,row.names=1)
		} else {
			stop('No tract stat file for the given genome!!!')
		}
		statTbl <- statTbl[rownames(statTbl)%in%tracts,,drop=F]
		nucTract <- statTbl[,gsub('^\\.','',paste(region,'nucTracts',sep='.'))]
		nTract <- statTbl[,gsub('^\\.','',paste(region,'nTracts',sep='.'))]
		b <- rep(as.numeric(b),length(nucTract))
		names(b) <- names(nucTract) <- rownames(statTbl)
		b[grepl('break',names(b))] <- 0 # for break "tracts", fix so-called boudary to 0.
		nucTract <- nucTract+2*b*nTract
		nucTract <- nucTract[tracts]
		nTotTracts <- sum(nucTract)
		#### Find which tracts are present in joined, possibly extend to a full scope
    joined[,ncol(joined)] <- gsub('.*\\-','break',joined[,ncol(joined)],perl=T) # assume break values are like ca_A-3, A_T-1
		tractInstances <- strsplit(joined[,ncol(joined)],';')
		tracts.incidence <- lapply( tracts,function(x,instances) sapply(instances,function(y,x) sum(y==x),x), tractInstances )
		nEdits_inTract <- switch(intMode,
			s=sapply(tracts.incidence,function(x) sum(x>0)),
			m=sapply(tracts.incidence,sum)
		)
		names(nEdits_inTract) <- tracts
		#### For subtype-rows and overall-row generate four numbers and expRate,obsRate, and p.binom. 
		stats.perType <- mapply(enrichStats_perType,nEdits_inTract,nucTract,n_hgEdit,nTotal)
		nE_inTract_total <- switch(intMode,
			s=nrow(joined),
			m=sum(sapply(tractInstances,length))
		)
		stat.whole <- enrichStats_perType(nE_inTract_total,nTotTracts,n_hgEdit,nTotal)
		enrichTbl <- cbind(stats.perType,stat.whole)
		enrichTbl <- t(enrichTbl)
		rownames(enrichTbl)[-nrow(enrichTbl)] <- tracts
		rownames(enrichTbl)[nrow(enrichTbl)] <- 'Overall'
	} else {
		enrichTbl <- matrix(NA,nr=0,nc=7)
	}
  enrichTbl <- enrichTbl[,c(2:(ncol(enrichTbl)-2),1,(ncol(enrichTbl)-1):ncol(enrichTbl)),drop=F]
  colnames(enrichTbl) <- c('nFeatures_intract','nucTract','nFeatures_ingenome','nucGenome','pEnrich','obsRate','expRate')
  enrichTbl
}

# enrichVis(): process the enrich table to plot 1) RR barplot & 2) RR pieplot.
### INPUT enrichTbl: table out of binomEnrich()
### INPUT pth: RR bars with p<=pth will be marked with an asterisk.
### INPUT joined: stem of output file name. Generated pictures will be saved as "joined."tif 
enrichVis <- function(enrichTbl,pth=0.01,joined='output') {
	dat <- data.frame(enrichTbl)
	types <- rownames(dat)
	n_type <- length(types) # subtypes plus "Overall"
	RR <- dat$obsRate/dat$expRate
	p <- dat$pEnrich
	p[is.na(p)] <- 1
	sigIx <- p<=pth
	names(p) <- names(RR) <- types
	tiff(paste0('output/',joined,'.tif'),width=1792,height=1792)
	layout(matrix(1:2,nr=2)) 
	############# barplot ##################
	par(cex.lab=3,mar=c(5,8,4,2))
	col=rep('gray',n_type)
	if (rownames(dat)[nrow(dat)]=='Overall')
		col[length(col)] <- 'black'
	mp=barplot(RR,ylim=c(0,max(RR)+1),las=1,border=NA, col=col,main='Relative Risk of tract-trapping vs genome baseline',
			cex.lab=3,cex.axis=2,cex.names=2.5,cex.main=3)
	title(ylab='Relative Risk (RR)',line=4)
	abline(h=1,lty='dashed',col='darkgray',lwd=4)
	points(mp[sigIx],RR[sigIx]+max(RR)/50,pch='*',cex=3)
	############# pieplot ###################
	if (rownames(dat)[nrow(dat)]=='Overall') {
			RR <- RR[-nrow(dat)]
			p <- p[-nrow(dat)]
			sigIx <- sigIx[-nrow(dat)]
			dat <- dat[-nrow(dat),]
	}# trim Overall entitity from each related R object.
	cols <- c(colorRampPalette(c('white','blue'))(148))
	p[p==0] <- 1e-147 # Extremely small p is flatted to 1e-147.
	p[p<1e-27] <- 1e-147
	MP = -log10(p)
	MP = floor(1+MP)
	MPpseudo <- cols[MP] # p, MP, and MPpseudo are created for showing p-proportional color pallete. Color is absolute scale from 1 to 147. 
	pieLabels <- names(MPpseudo) <- rownames(dat)
	pieLabels[sigIx] <- paste0(pieLabels[sigIx],'*')
	slice=RR
	pie(slice[RR!=0],labels=pieLabels[RR!=0],col=MPpseudo[RR!=0],cex=2.5,main='Normalized abundance (slice) and enrichment significance (shade) of events in tracts',cex.main=3)
	sink('output/null')
        dev.off()
	file.remove('output/null')
}
# signifBarplot(): plot a landscape barplot as a background for user's Overall enrichment p-value. 
### INPUT enrichTbl: table out of binomEnrich(). 
### INPUT joined: stem of output file name. Generated picture will be saved as "joined."landscape.tif
signifBarplot <- function(enrichTbl,joined='output',landscape='tracts/polytrap_P_landscape.csv') {
	uInvP <- -log10(enrichTbl['Overall','pEnrich'])
	if (!is.na(uInvP)) {
		landscape <- read.csv(landscape,check.names=F,as.is=T)
		rownames(landscape) <- landscape$shownName
		landscape[landscape==0] = 1e-314
		invP <- -log10(landscape[,'Tri-'])
		invP <- c(invP,uInvP)
		names(invP) <- c(rownames(landscape),'USER_DATA')
		invP <- sort(invP,decreasing=T)
		cols <- rep('lightgray',length(invP))
		cols[names(invP)=='USER_DATA'] <- 'red'
		tiff(paste0('output/',joined,'.landscape.tif'),width=1792,height=1024)
		par(mar=c(14,8,8,2),cex.axis=1.5)
		coord = barplot(invP,names.arg=names(invP),las=2,ylab='',col=cols, ylim=c(0,314),
			main='New Data in Human TNR Enrichment landscape',cex.main=2,cex.lab=1.7,cex.names=1.5) #cex.axis=1.5,cex.names=2,
		mtext(text='Association significance\n( -log10(P) )',side=2,line=4,cex=2)
		axis(1,at=coord[invP>=4],labels=names(invP)[invP>=4],col.axis='blue',cex=1.5,las=2)
		axis(1,at=coord[names(invP)=='USER_DATA'],labels='USER_DATA',col.axis='red',cex=1.5,las=2)
		lines(c(min(coord),max(coord)+1),c(4,4),lty='dashed',lwd=2,col='blue')
		text(max(coord)-5,10,'p=1e-4',cex=2,col='blue')
		dev.off()
	}
}
######################################################
## Process argument and invoke subfunctions ##########
######################################################
args=commandArgs(T)
for (i in 1:length(args)) {
	eval(parse(text=args[i]))
}
if (exists('BREAKS')) {
	tracts=BREAKS	
	enrichTbl <- binomEnrich(joined,gn,tracts,0,arg.d,intMode,region)
} else {
	tractsCmmd=paste0('tracts=tract',arg.t)
	eval(parse(text=tractsCmmd))
	if (!any(is.na(breaks))) {
		tracts=c(tracts,breaks)
	}
	enrichTbl <- binomEnrich(joined,gn,tracts,b,arg.d,intMode,region)
}
write.table(enrichTbl,paste0('output/',joined,'.enrich'),sep='\t',quote=F,col.names=NA)
if (nrow(enrichTbl)>0) {
	enrichVis(enrichTbl,0.01,joined)
	signifBarplot(enrichTbl,joined)
}
