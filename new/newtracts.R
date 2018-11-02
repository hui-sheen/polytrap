source('~/polytract/R/groupTri.R')
source('~/polytract/R/reduceOverlapTract.R')
source('~/polytract/R/gapBreaks.R') # UPDATE 9/27: make 4th column contain both pattern and breakLen (e.g. "AT_T-3")
# UPDATE 9/27/2018 relative to workable: Combine the task of find_tract_break.R & tract_anno_type.pl in the same current directory.
# RELIEs on gn2bsgenome.map file in the same working directory.
      #genomePkg can be 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
	################## 'BSgenome.Mmusculus.UCSC.mm9', 'BSgenome.Mmusculus.UCSC.mm10', 
	################## 'Rnorvegicus.UCSC.rn6','BSgenome.Dmelanogaster.UCSC.dm6', etc.
      #matchedPatterns include 4 single, 6 dinuc, and 10 trinuc patterns.
gn2bsgenome <- read.csv('genomes.meta',as.is=T) # UPDATE 10052018
gns <- commandArgs(TRUE)
findPolyTract<-function(genomePkg,matchedPatterns) {  
      require(stringr)
      if (!require(genomePkg,character.only=T)) {
             source('http://www.bioconductor.org/biocLite.R')
             biocLite(genomePkg)
      }
      require(genomePkg,character.only=T)
      Genome=get(genomePkg)
      resultOutAll<-NULL
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
                   resultOut1<-data.frame(chr=chrName,matchedPosition1,stringsAsFactors=FALSE) #,Pattern=matchedPattern
                   resultOut1<-makeGRangesFromDataFrame(resultOut1,keep.extra.columns=TRUE)
                   gr<-reduce(c(resultOut1))
                   resultOut <- data.frame(chr=gsub('chr|Chr','',as.character(seqnames(gr))),
                                start=start(gr),
                                end=end(gr),
                                Pattern=matchedPattern,
                                stringsAsFactors=FALSE)
                   resultOutAll<-rbind(resultOutAll,resultOut)
             }
      }
      Pattern1 <- sapply(strsplit(resultOutAll$Pattern,'\\(|\\)|\\|'),function(x) x[2])
      resultOutAll <- cbind(resultOutAll[,c('chr','start','end')],Pattern=Pattern1)
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
merge_di_pair <- function(rPattern) {
	type6 <- c('ta','ga','ca','gt','ct','gc')
	TYPE6 <- toupper(type6)	
	TYPE6rev <- sapply(lapply(strsplit(TYPE6,''),function(x) x[2:1]),function(x) paste(x,collapse=''))
	type6.repeated <- rep(type6,2)
	names(type6.repeated) <- c(TYPE6,TYPE6rev)
	nPattern <- type6.repeated[as.character(rPattern)]
  nPattern
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
        if (modality != 'single') {
          if (modality=='di') {
             nPatterns = merge_di_pair(tracts0$Pattern)
             tmp <- cbind(tracts0[,1:3],nPatterns)
             tracts <- reduceOverlapTract(tmp)             
          }
          if (modality=='tri') {
             mega10 <- groupTri()
             tri60_tri10 <- cbind(unlist(mega10),rep(names(mega10),each=6)) #Each mega-group conists of 6 types.
             merged <- merge(tracts0,tri60_tri10,by.x=4,by.y=1)[,c(2:5)]
             tracts <- reduceOverlapTract(merged)
          }
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
