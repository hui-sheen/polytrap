# UPDATE 1/15/2020 Muted the two reduceOverlapTract() in succession to findPolyTract(). Within findPolyTract(), employed reduce() across sequentially similar patterns (say ACT/CTA/TAC). DNRs & TNRs are left as original tract instances without being reduced across raw types. Nominally, un-reduced consecutive tracts are tagged with a same meta-type.
# UPDATE 1/15/2020 merge_di() is changed to reflect four meta-groups, including ca/gt and ct/ga meta-groups.
# UPDATE 1/15/2020 merge_single() to perform meta-grouping of 4 Single groups to 2 meta-groups: A/T, C/G.
# UPDATE 9/27/2018 relative to workable: Combine the task of find_tract_break.R & tract_anno_type.pl in the same current directory.
#source('~/polytract/R/groupTri.R')
#source('~/polytract/R/reduceOverlapTract.R')
#source('~/polytract/R/gapBreaks.R') # UPDATE 9/27: make 4th column contain both pattern and breakLen (e.g. "AT_T-3")
# RELIEs on gn2bsgenome.map file in the same working directory.
      #genomePkg can be 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
	################## 'BSgenome.Mmusculus.UCSC.mm9', 'BSgenome.Mmusculus.UCSC.mm10', 
	################## 'Rnorvegicus.UCSC.rn6','BSgenome.Dmelanogaster.UCSC.dm6', etc.
      #matchedPatterns include 4 single, 6 dinuc, and 10 trinuc patterns.
gn2bsgenome <- read.csv('genomes.meta',as.is=T) # UPDATE 10052018
gns <- commandArgs(TRUE)
sortChar_in_array <- function(charArray) sapply(lapply(strsplit(charArray,''),function(x) sort(x)),function(x) paste(x,collapse=''))
revChar_in_array <- function(charArray) sapply(lapply(strsplit(charArray,''),function(x) x[length(x):1]),function(x) paste(x,collapse=''))
# patternRun(): to generate a mapping chart for raw 4/12/60 patterns to 4/6/20 meta-patterns with re-ordered patterns combined. 
patternRun <- function() {
	######## Single ######
	chart.single <- c('A','T','G','C')
	names(chart.single) <- chart.single
	######## DI ##########
	di6 <- c('TA','GA','CA','GT','CT','GC')
	di6.2 <- revChar_in_array(di6)
	chart.di <- rep(di6,2)
	names(chart.di) <- c(di6,di6.2)
	######## TRI ########
	order231 <- function(charArray) sapply(lapply(strsplit(charArray,''),function(x) x[c(2,3,1)]),function(x) paste(x,collapse=''))
	order312 <- function(charArray) sapply(lapply(strsplit(charArray,''),function(x) x[c(3,1,2)]),function(x) paste(x,collapse=''))
	tri10 <- c('AAC','AAG','AAT','ACC','GAC','ACT','CAG','AGG','ATC','CGG')
	tri10.rev <- c('GTT','CTT','ATT','GGT','GTC','AGT','CTG','CCT','GAT','CCG')
	chart.tri <- c(rep(tri10,3),rep(tri10.rev,3))
	names(chart.tri) <- c(tri10,order231(tri10),order312(tri10),tri10.rev,order231(tri10.rev),order312(tri10.rev))
	chart <- c(chart.single,chart.di,chart.tri) 
}
pattern.chart <- patternRun()
findPolyTract<-function(genomePkg,matchedPatterns,chart=pattern.chart) {  
      require(stringr)
      if (!require(genomePkg,character.only=T)) {
             source('http://www.bioconductor.org/biocLite.R')
             biocLite(genomePkg)
      }
      require(genomePkg,character.only=T)
      Genome=get(genomePkg)
      #resultOutAll<-NULL
			gr.all <- NULL
      for (chrName in standardChromosomes(Genome)) {
             print(chrName)
             chrSeq <- Genome[[chrName]]
             chrSeq<-as.character(chrSeq)
             for (j in 1:length(matchedPatterns)) {
                   matchedPattern<-matchedPatterns[j]
                   matchedPosition1<-str_locate_all(chrSeq, matchedPattern)[[1]]
                   if (nrow(matchedPosition1)==0) {
                          next;
                   }
                   gr<-data.frame(chr=chrName,matchedPosition1,Pattern=matchedPattern,stringsAsFactors=FALSE) #,Pattern=matchedPattern
                   gr<-makeGRangesFromDataFrame(gr,keep.extra.columns=TRUE)
                   gr.all<-c(gr.all,gr)
             }
      }
			gr.all <- unlist(GRangesList(gr.all))
      Pattern1 <- sapply(strsplit(mcols(gr.all)$Pattern,'\\(|\\)|\\|'),function(x) x[2])
			Pattern1 <- chart[Pattern1]#Pattern1 <- sortChar_in_array(Pattern1)	
			mcols(gr.all)$Pattern <- Pattern1
			############ reduce sequentially reordered patterns, not across revComp boudary  ####################
			grl <- split(gr.all,factor(mcols(gr.all)$Pattern))
			grl <- lapply(grl,reduce)
			resultOutAll <- NULL
			for (pa in names(grl)) {
				gr <- grl[[pa]]
				resultOut <- data.frame(chr=gsub('chr|Chr','',as.character(seqnames(gr))),
																	start=start(gr),
																	end=end(gr),
																	Pattern=pa,
																	stringsAsFactors=FALSE
				)
				resultOutAll <- rbind(resultOutAll,resultOut)
			}
      return(resultOutAll)      
}
# INPUT rPattern: (raw Pattern) 12 dinuc combinations.
# OUTPUT nPattern: (new Pattern) 6 dinuc combinations. ta for TA/AT; ga for GA/AG; ca for CA/AC; gt for GT/TG; ct for CT/TC; gc for GC/CG.
genTriple <- function() {
  temp1=expand.grid(x1 = c('A','T','C','G'), x2 = c('A','T','C','G'), x3 = c('A','T','C','G'))
  temp2=apply(temp1,1,function(x) length(unique(x)))
  temp1=temp1[which(temp2>1),] #remove all identical (such as GGG)
  matchedPatterns=as.character(apply(temp1,1,function(x) paste0(x,collapse='')))
  matchedPatterns=paste0('(',matchedPatterns,'){3,}')
  matchedPatterns
}
# merge_di(): Merge 12 raw dinuc patterns to 4 meta-types: ta,ca/gt,ct/ga,gc
merge_di <- function(rPattern) {
  type4 <- c('ta','ca/gt','ct/ga','gc')
	rawS <- sapply(type4,function(x) strsplit(toupper(x),'/'))
	revRawS <- lapply(rawS,revChar_in_array)
	oTypeS <- mapply(c,rawS,revRawS)
	chart <- rep(names(oTypeS),sapply(oTypeS,length))
	names(chart) <- unlist(oTypeS)
  nPattern <- chart[as.character(rPattern)]
  nPattern
}
merge_single <- function(rPattern) {
	type2 <- c('A/T','C/G')
	oTypeS <- sapply(type2,function(x) strsplit(toupper(x),'/'))
  chart <- rep(names(oTypeS),sapply(oTypeS,length))
  names(chart) <- unlist(oTypeS)
  nPattern <- chart[as.character(rPattern)]
  nPattern
}
# merge_di_pair0(): obsolete from merge_di_pair() [Sep 2018]
merge_di_pair0 <- function(rPattern) {
	type6 <- c('ta','ga','ca','gt','ct','gc')
	TYPE6 <- toupper(type6)	
	TYPE6rev <- sapply(lapply(strsplit(TYPE6,''),function(x) x[2:1]),function(x) paste(x,collapse=''))
	type6.repeated <- rep(type6,2)
	names(type6.repeated) <- c(TYPE6,TYPE6rev)
	nPattern <- type6.repeated[as.character(rPattern)]
  nPattern
}
# groupTri(): generate the 60-to-10 mapping of trinuc types.
#### OUTPUT: a list of 10 components, each consisting of three frameshit runs and three complementary frameshit runs.
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

#########################################################################################
######### sub-functions invoked by gapBreaks() ##########################################
#########################################################################################
### UPDATE 9/27/2018 relative to a workable version: Make the tracts have only 4 columns (rather than 5 columns), where the 4th is like ca_A-3, A_T-1, etc..
# gapBreaks(): main function.
gapBreaks <- function(tract,b=c(1,2,3)[1]) {
#       tract.single <- read.csv(paste0('polytract_',hg,'_single.csv'),header=F,as.is=T)
#       tract.di <- read.csv(paste0('polytract_',hg,'_di.csv'),header=F,as.is=T)
#       tract <- rbind(tract.single,tract.di)
        colnames(tract) <- c('chr','start','end','type')
        chr.fct <- factor(tract$chr)
        gaps.chrLst <- by(tract[,1:3],chr.fct,findGaps,b=b,simplify=F)
        tStarts.chrLst <- by(tract[,c(2,4)],chr.fct,function(x) {y=x[,2];names(y)=x[,1];return(y)})
        tEnds.chrLst <- by(tract[,c(3,4)],chr.fct,function(x) {y=x[,2];names(y)=x[,1];return(y)})
        breaks <- mapply(labelBreaks,gaps.chrLst,tStarts.chrLst,tEnds.chrLst,SIMPLIFY=F)
        breaks <- do.call(rbind,breaks)
        breaks <- breaks[order(as.numeric(breaks[,1]),breaks[,2]),]
        breaks
}
# findGaps() <- gapBreaks()
findGaps <- function(chrIntv,b=c(1,2,3)[1]) { #f1:chr,f2:start,f3:end
        library(GenomicRanges)
        chrIntv <- data.frame(chrIntv)
        gr <- makeGRangesFromDataFrame(chrIntv)
        gaps.gr <- gaps(gr)
        gaps <- data.frame(
                chr=as.character(seqnames(gaps.gr)),
                start=start(gaps.gr),
                end=end(gaps.gr),
                len=end(gaps.gr)-(start(gaps.gr)-1),
                stringsAsFactors=F
        )
#       cat('All gaps numbered:',nrow(gaps),'\t')
        gaps <- subset(gaps,len==b)
#       cat('Length-qualified gaps numbered:',nrow(gaps),'\n')
        gaps
}
# labelBreaks() <- gapBreaks()
labelBreaks <- function(breaks,tStarts,tEnds) {# Work on tracts of one same chr
        colnames(breaks) <- c('chr','start','end','len')
        labels <- paste(tEnds[as.character(breaks[,'start']-1)],tStarts[as.character(breaks[,'end']+1)],sep='_')
        breaks <- data.frame(breaks[,-4],label=paste(labels,breaks[,4],sep='-'))#label=labels,len=breaks[,4]) # UPDATE 9/27/2018
        breaks
}
#########################################################################################
##########  Generate 3 basic tract files for input genomes  #############################
#########################################################################################
singlePatterns <- c('(A){6,}','(T){6,}','(C){6,}','(G){6,}')
doublePatterns <- c('(AT){3,}','(AG){3,}','(AC){3,}','(TG){3,}','(TC){3,}','(CG){3,}','(TA){3,}','(GA){3,}','(CA){3,}','(GT){3,}','(CT){3,}','(GC){3,}')
triplePatterns <- genTriple()
bsgenome <- gn2bsgenome[,'BSGENOME']
names(bsgenome) <- gn2bsgenome[,'GENOME']
genomePkgs <- bsgenome[gns]#c(mm10='BSgenome.Mmusculus.UCSC.mm10',mm9='BSgenome.Mmusculus.UCSC.mm9')
patterns <- list(single=singlePatterns,di=doublePatterns,tri=triplePatterns) #list(tri=triplePatterns)#
tracts.lst <- vector('list',length(genomePkgs))
names(tracts.lst) <- names(genomePkgs)
if (!file.exists('newtracts'))
	dir.create('newtracts')
for (gn in names(genomePkgs)) {
    cat(gn,'...\n')
    tracts.lst[[gn]] <- vector('list',length(patterns))
    names(tracts.lst[[gn]]) <- names(patterns)
    for (modality in names(patterns)) {
        cat(date(),'\t',modality,'...\n')
        tracts <- tracts0 <- findPolyTract(genomePkgs[gn],patterns[[modality]])
				if (modality%in%c('single','di')) {
					nPatterns = switch(modality,
						single=merge_single(tracts0$Pattern),
						di=merge_di(tracts0$Pattern)
					)
					tracts <- cbind(tracts0[,1:3],nPatterns)
				} else {
					 mega10 <- groupTri()
					 tri60_tri10 <- cbind(unlist(mega10),rep(names(mega10),each=6)) #Each mega-group conists of 6 types.
					 tracts <- merge(tracts0,tri60_tri10,by.x=4,by.y=1)[,c(2:5)]
					 #tracts <- reduceOverlapTract(merged)
				}
        tracts <- tracts[order(as.numeric(tracts[,1]),as.character(tracts[,4]),tracts[,2]),]
        write.table(tracts,paste0('newtracts/polytract_',tolower(gn),'_',modality,'.csv'),sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)
        tracts.lst[[gn]][[modality]] <- tracts
    }
}
##########################################################################################
########### Generate 3 breaks files based on single & di tracts ##########################
##########################################################################################
bs <- c(1,2,3)
breaks.lst <- vector('list',length(gns))
names(breaks.lst) <- gns
for (gn in gns) {
        cat('\n',gn,':\n')
        breaks.lst[[gn]] <- vector('list',length(bs))
        names(breaks.lst[[gn]]) <- bs
        tract.single <- read.csv(paste0('newtracts/polytract_',tolower(gn),'_single.csv'),header=F,as.is=T)
        tract.di <- read.csv(paste0('newtracts/polytract_',tolower(gn),'_di.csv'),header=F,as.is=T)
        tract <- rbind(tract.single,tract.di)
        for (b in bs) {
                cat('breaks',b,'... ')
                breaks.lst[[gn]][[b]] <- breaks <- gapBreaks(tract,b)
                write.table(breaks,paste0('newtracts/polytract_breaks_',tolower(gn),'_',b,'.csv'),row.names=F,col.names=F,quote=F,sep=',')
        }
}

#save.image(paste0(paste(gns,collapse='_'),'.tract_break.RData'))
