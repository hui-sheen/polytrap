# Rscript binomEnrich.R joined='example.out' gn='hg38' arg.t=12 b=1 breaks='break123' arg.d=',' region=NULL # args.output, args.genome, args.tract, str(args.boundary),  args.delimiter, args.header
##### OBJECTIVE: to calculate Binomial probability of feature enrichment within tracts, and visualize the enrichment (increased risk) in barplot & pieplot.
##### INPUT 1: prior python output file (-o argument) 
##### INPUT 2: genome (-g argument)
##### INPUT 3: tract modality (-t argument), or JUNCTION modality (-J argument)
##### INPUT 4: boundary (-b argument)
##### INPUT 5: hinge/junction (-j argument)
##### INPUT 6: delimiter (-d argument) 
##### INPUT 7: genomic region constraint (-r argument)
##### NOTE 1: Nine core assemblies take the total number of nucleotides from global variable "gnNucs". Other assemblies look for this statistics from genomes.meta.
##### NOTE 2: Genome-specific stat files must be saved as .../tracts/polytract_GN_stats.txt.
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
##### INPUT arg.d: file separator. comma or tab.
##### INPUT headed: file has header?
##### INPUT b: boundary specification (adjacent nucleotide extension).
##### INPUT tracts: Perform subtype enrichment analysis on only the given subtypes. {tract12,tract23,tract13,tract123,'break1','break2','break3'}  
##### NOTE: this function requires static files on tract statistics, e.g., tracts/polytract_hg38_stats.txt
binomEnrich <- function(joined,gn=c('HG38','HG19','MM9','MM10')[1],tracts=tract123,b=c(1,2,3,0)[1],arg.d=',',region=NULL,headed=F) {#,nTotal=gnNucs[tolower(gn)]) {
	#### Read in previous output file (joined)
	if (arg.d==',') {
		joined <- unique(read.csv(paste('output',joined,sep='/'),header=headed,as.is=T))
	} else {
		joined <- unique(read.delim(paste('output',joined,sep='/'),header=headed,as.is=T))
	}
	joined[,ncol(joined)] <- as.character(joined[,ncol(joined)])
	n_hgEdit <- nrow(joined)
	if ( is.element( tolower(gn),names(gnNucs) ) ) {
		nTotal=gnNucs[tolower(gn)]
	} else {
		tmp <- read.csv('./new/genomes.meta',as.is=T)[,1:2]
		gnNucs <- tmp[,2]
		names(gnNucs) <- tolower(tmp[,1])
		nTotal <- gnNucs[tolower(gn)]
	}
	joined <- joined[joined[,ncol(joined)]!='0',]
	if (any(joined[,ncol(joined)]!='0')) {
		#### Obtain full tracts statistics from the static file ####
		statFile <- paste0('tracts/polytract_',tolower(gn),'_stats.txt')
		if (file.exists(statFile)) {
			statTbl <- read.delim(statFile,row.names=1)
		} else {
			stop('No tract stat file for the concerned genome!!!')
		}
		statTbl <- statTbl[rownames(statTbl)%in%tracts,,drop=F]
		nucTract <- statTbl[,gsub('^\\.','',paste(region,'nucTracts',sep='.'))]
		nTract <- statTbl[,gsub('^\\.','',paste(region,'nTracts',sep='.'))]
		b <- rep(as.numeric(b),length(nucTract))
		names(b) <- names(nucTract) <- rownames(statTbl)
		b[grepl('break',names(b))] <- 0 # for subset of break tracts, fix so-called boudary to 0.
		nucTract <- nucTract+2*b*nTract
		nucTract <- nucTract[tracts]
		nTotTracts <- sum(nucTract)
		joined[,ncol(joined)] <- gsub('.*\\-','break',joined[,ncol(joined)],perl=T) # assume break values are like ca_A-3, A_T-1
		#### Find which tracts are present in joined, possibly extend to a full scope
		tractInstances <- strsplit(joined[,ncol(joined)],';') ### in future, potentially use semicolon to deliver nested cell.
		tracts.idx <- lapply( tracts,function(x,instances) sapply(instances,function(y,x) is.element(x,y),x), tractInstances )
		nEdits_inTract <- sapply(tracts.idx,sum)
		names(nEdits_inTract) <- tracts
		#### Obtain nucTract & nTract information from the static stat file or derive them from tracts files
		stats.perType <- mapply(enrichStats_perType,nEdits_inTract,nucTract,n_hgEdit,nTotal)
		stat.whole <- enrichStats_perType(nrow(joined),nTotTracts,n_hgEdit,nTotal)
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

# enrichVis(): process the enrich table to plot 1) barplot & 2) pieplot.
enrichVis <- function(enrichTbl,pth=0.01,joined='output') {
	dat <- data.frame(enrichTbl)
	types <- rownames(dat)#gsub('all.','',grep('all.',rownames(dat),value=T),fixed=T)
		n_type <- length(types) # subtypes plus "Overall"
		RR <- dat$obsRate/dat$expRate
		p <- dat$pEnrich
		p[is.na(p)] <- 1
		sigIx <- p<=pth
		names(p) <- names(RR) <- types
		tiff(paste0('output/',joined,'.tif'),width=1792,height=1792)
		layout(matrix(1:2,nr=2)) 
		# barplot
		par(cex.lab=3,mar=c(5,8,4,2))
		col=rep('gray',n_type)
		if (rownames(dat)[nrow(dat)]=='Overall')
			col[length(col)] <- 'black'
		mp=barplot(RR,ylim=c(0,max(RR)+1),las=1,border=NA, col=col,main='Relative Risk of tract-trapping vs genome baseline',
				cex.lab=3,cex.axis=2,cex.names=2.5,cex.main=3)
		title(ylab='Relative Risk (RR)',line=4)
		abline(h=1,lty='dashed',col='darkgray',lwd=4)
		points(mp[sigIx],RR[sigIx]+max(RR)/50,pch='*',cex=3)
		# pieplot
		if (rownames(dat)[nrow(dat)]=='Overall') {
				RR <- RR[-nrow(dat)]
				p <- p[-nrow(dat)]
				sigIx <- sigIx[-nrow(dat)]
				dat <- dat[-nrow(dat),]
		}
		cols <- c(colorRampPalette(c('white','blue'))(148))
		p[p==0] <- 1e-147
		p[p<1e-27] <- 1e-147
		MP = -log10(p)
		MP = floor(1+MP)#((MP-min(MP,na.rm=T))/max(MP,na.rm=T))*100)
		MPpseudo <- cols[MP] 
		pieLabels <- names(MPpseudo) <- rownames(dat)
		pieLabels[sigIx] <- paste0(pieLabels[sigIx],'*')
		slice=RR #dat$nFeatures_intract
		pie(slice[RR!=0],labels=pieLabels[RR!=0],col=MPpseudo[RR!=0],cex=2.5,main='Normalized abundance (slice) and enrichment significance (shade) of events in tracts',cex.main=3)
	sink('output/null')
        dev.off()
	file.remove('output/null')
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
	enrichTbl <- binomEnrich(joined,gn,tracts,0,arg.d,region)
} else {
	tractsCmmd=paste0('tracts=tract',arg.t)
	eval(parse(text=tractsCmmd))
	if (!any(is.na(breaks))) {
		tracts=c(tracts,breaks)
	}
	enrichTbl <- binomEnrich(joined,gn,tracts,b,arg.d,region)
}
#enrichTbl <- binomEnrich(joined,gn,tracts,b,commaSep)
write.table(enrichTbl,paste0('output/',joined,'.enrich'),sep='\t',quote=F,col.names=NA)
if (nrow(enrichTbl)>0)
	enrichVis(enrichTbl,0.01,joined)
