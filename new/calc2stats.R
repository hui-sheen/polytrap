# UPDATE 1/15/20 to reflect new series of meta-types, particularly two Single-types & two DI-types.
# groupTri() is utilized mainly to indicate the order of the 10 tri-types.
# groupTri() <- splitByType().
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
# splitByType() <- get2statByList()
# INPUT tract0: a data.frame with read-in tracts. Can be stacked tracts of diverse types (single, di, tri, and breaks).
splitByType <- function(tract0) {
	tract1 <- c('A/T','C/G') #c('A','T','C','G')
	tract2 <- c('ta','ct/ga','ca/gt','gc') #tolower(c('TA','GA','CA','GT','CT','GC')) # Use lower-case to denote 6 dinuc types.
	tract3 <- names(groupTri())
	tract.types <- c(tract1,tract2,tract3)
	Type <- factor(tract0[,4],levels=tract.types,ordered=T)#levels=intersect(tract.types,unique(tract0[,4])),ordered=T) # UPDATE 9/28
	tract.lst <- by(tract0,Type,function(x) x)
	tract.lst
}
#get2statByList() <- getRegionStat(),calc2stats()
get2statByList <- function(tracts) {
        tracts.lst <- splitByType(tracts)
        nTracts <- sapply(tracts.lst,nrow)
	if (nrow(tracts)>0) {
	        nucTracts <- sapply(tracts.lst, function(x) sum(x[,3]-(x[,2]-1))        )
	       	stat <- cbind(nTracts=nTracts,nucTracts=nucTracts)
		stat[is.null(stat)] <- 0
	} else {
		stat <- matrix(0,nr=length(tracts.lst),nc=2)
	}
        rownames(stat) <- names(tracts.lst)
	stat
}
#getRegionStat() <- calc2stats()
getRegionStat <- function(tracts,isBreak=F) {
	colnames(tracts) <- c('chr','start','end','pattern','region')
	regions <- c('protein_coding','pseudogene','lncRNA')#unique(tracts[,5])
	regions.label <- c('protein','pseudo','lncRNA','NA')
	for (i in 1:(length(regions)+1)) {
		if (i<=length(regions)) {
			tracts.i <- subset(tracts,region==regions[i])
		} else {
			tracts.i <- subset(tracts,is.na(region))
		}
		if (i==1) {
			if (!isBreak) {
				stat <- get2statByList(tracts.i)
			} else {
				stat <- c(nrow(tracts.i),sum(tracts.i[,3]-tracts.i[,2]+1))
			}
		} else {
			if (!isBreak) {
				stat <- cbind(stat,get2statByList(tracts.i))
			} else {
				stat <- c(stat,c(nrow(tracts.i),sum(tracts.i[,3]-tracts.i[,2]+1)))
			}
		}
	}
	stat <- matrix(stat,nc=(length(regions)+1)*2)
	colnames(stat) <- paste(rep(regions.label,each=2),rep(c('nTracts','nucTracts'),4),sep='.')
	stat
}
calc2stats <- function(gn='MM10',dir='./newtracts') {
	tract3 <- c('single','di','tri')
	break3 <- c(1,2,3)
	tract3files <- paste0(dir,'/polytract_',tolower(gn),'_',tract3,'.csv')
	break3files <- paste0(dir,'/polytract_breaks_',tolower(gn),'_',break3,'.csv')
	#files <- c(tract3files)#,break3files) # Not considering break files in the same way.
	######## Handling tract files ##############
	for (i in 1:length(tract3files)) {
		dati <- read.delim(tract3files[i],header=F,sep=',',as.is=T)
		if (i==1) {
			tracts <- dati
		} else {
			tracts <- rbind(tracts,dati)
		}
	}	
	stat <- get2statByList(tracts)	
	if (ncol(tracts)>4) {
		stat.regions <- getRegionStat(tracts) 
		stat <- cbind(stat,stat.regions)	
	}
	############# Handling break files #############
	for (i in 1:3) {
		dat <- read.delim(break3files[i],header=F,sep=',',as.is=T)
		break.stat.i <- matrix(c(nrow(dat),sum(dat[,3]-dat[,2]+1)),nr=1)
		if (ncol(dat)>4) {
			break.stat.regions <- getRegionStat(dat,isBreak=T)	
			break.stat.i <- cbind(break.stat.i,break.stat.regions)
		}
		if (i==1) {
			break.stat <- break.stat.i
		} else {
			break.stat <- rbind(break.stat,break.stat.i)
		}
	}
	stat <- rbind(stat,break.stat)
	rownames(stat)[(nrow(stat)-2):nrow(stat)] <- paste0('break',1:3)
	############ Save to stat file ###################
	write.table(stat,paste0(dir,'/polytract_',tolower(gn),'_stats.txt'),col.names=NA,sep='\t',quote=F)
}

gns <- commandArgs(TRUE)
#gns <- c(rn6 galGal5 danRer10 dm6 sacCer3)
#gns <- c('hg38','hg19','mm10','mm9')
for (gn in gns) {
	cat('Calculating tract statistics for genome',gn,'...\n')
	calc2stats(gn,dir='./newtracts')
}
