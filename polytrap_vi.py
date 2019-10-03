#!/usr/bin/env python3
# python polytrap polytrap.py -d , -i example.input -o example
import gzip
import argparse
import os
from argparse import RawTextHelpFormatter
from subprocess import call
#error messages
em1="##########ERROR: OPENING POLYTRACT FILE failed, please seek the bulk TRACTS/HINGE files and PLACE them under polytrap/tracts. please see more detailed error message below:"
em2="##########ERROR: OPENING INPUT FILE, input file may not be supplied correctly, please see more detailed error message below:"

#arguments parsing
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--input", help = "your input file in BED format, (required)", required=True)
parser.add_argument("-o", "--output", help = "output file for the results (required)", required=True)
parser.add_argument("-g", "--genome", help = "hg19 = Genome Reference Consortium Human Build 37 (GRCh37)\n" +
                                             "hg38 = Genome Reference Consortium Human Build 37 (GRCh38)\n" +
											 "rhemac8 = Rhesus monkey assembly rheMac8\n" +	
                                             "mm9 = mouse genome assembly MGSCv37\n" +
                                             "mm10 = mouse genome assembly GRCm38\n" +
											 "rn6 = rat genome assembly Rnor_6.0\n" +
											 "canfam3 = dog genome assembly canFam3\n" 
											 "galal5 = chicken genome assembly Gallus_gallus-5.0\n" +
                                             "danrer10 = zebrafish genome assembly GRCz10\n" +
                                             "dm6 = fruitfly genome assembly BDGP6\n" +
                                             "saccer3 = yeast (Saccharomyces cerevisiae) genome assembly R64\n",
                                             default="hg38", choices=["hg19", "hg38", "rhemac8","mm9", "mm10","rn6","canfam3","galgal5","danrer10", "dm6", "saccer3"], type = str.lower) 
parser.add_argument("-t", "--tract", help = "1 = single nucleotide tract\n" +
                                            "2 = dinucleotide tract\n3 = trinucleotide tract\n" +
                                            "12 = single + dinucleotide tracts\n" +
                                            "13 = single + trinucleotide tracts\n" +
                                            "23 = dinucleotide + trinucleotide tracts\n" +
                                            "123 = single + dinucleotide + trinucleotide tracts\n",
                                            choices = ["1", "2", "3", "12", "23", "13", "123"], default="12")
parser.add_argument("-b", "--boundary", help = "Number of  up and downstream nucleotides from the polytract to include in the analysis", type = int, choices=range(0,11),default=1)
parser.add_argument("-j", "--junction", help = "Map input to the junctions (hinges) between tracts, it is used in conjunction with -t.", type = int, choices=range(1,4))
parser.add_argument("-J", "--JUNCTION", help = "Map input to the junctions (hinges)  between tracts only.", type = int, choices=range(1,4))
parser.add_argument("-r", "--region", help = "Consider tracts located only in specified genomic region: protein-coding, pseudogene, lncRNA", choices=["protein","pseudo","lncRNA","NA"])
parser.add_argument("-H", "--header", help = "1 = yes, 0 = no, default is 1", default="1", choices=["0", "1"])

parser.add_argument("-d", "--delimiter", help = "delimiter for input file, choice are ',' and tab, default is tab", default="t", choices=[",", "t"])
parser.add_argument("-I","--intsect",help="intersection mode: s for singleton, m for multiplex",default="s",choices=["s","m"])
args = parser.parse_args()
if args.delimiter == "t":
 d = "\t"
else:
 d = ","

boundary = args.boundary
#if boundary > 100:
# print "##########WARNING: boundary too big may cause overlap polytracts and generating unrelaible results"
 
print "starting polytract enrichment analysis with the following arguments:"
print "\tgenome = " + args.genome
print "\tdelimiter = " + args.delimiter
print "\tintersection = " + args.intsect
print "\tinput = " + args.input
print "\toutput = " + args.output
print "\ttract = " + args.tract
print "\theader = " + args.header
print "\tboundary = " + str(args.boundary)
if args.region is not None:
 print "\tregion = " + args.region + " [Applicable only to latest genome assemblies (e.g., HG38 but not HG19)]"
if args.junction is not  None:
 print "\tjunction (hinge) = " + str(args.junction)
if args.JUNCTION is not  None:
 print "\tJUNCTION (hinge) = " + str(args.JUNCTION) + " If JUNCTION is specified, -t (tract) option is ignored"


#loading input files
try:
 input = open(args.input, "r")
except (OSError, IOError) as e:
 print em2
 raise Exception(e)

#loading polytract files
def make_tractlist(opt_g, opt_J, opt_t, opt_j):
        tract_list=[]
        option_j = {1:["tracts/polytract_breaks_"+opt_g.lower()+"_1.csv.gz"],
        2:["tracts/polytract_breaks_"+opt_g.lower()+"_1.csv.gz", "tracts/polytract_breaks_"+opt_g.lower()+"_2.csv.gz"],
        3:["tracts/polytract_breaks_"+opt_g.lower()+"_1.csv.gz", "tracts/polytract_breaks_"+opt_g.lower()+"_2.csv.gz", "tracts/polytract_breaks_"+opt_g.lower()+"_3.csv.gz"]}
        if opt_J is None:
                if "1" in opt_t:
                        tract_list.append("tracts/polytract_"+opt_g.lower()+"_single.csv.gz")
                if "2" in opt_t:
                        tract_list.append("tracts/polytract_"+opt_g.lower()+"_di.csv.gz")
                if "3" in opt_t:
                        tract_list.append("tracts/polytract_"+opt_g.lower()+"_tri.csv.gz")
                if opt_j is not None:
										if opt_j in option_j:
                        	tract_list.extend(option_j[opt_j])
        elif opt_J in option_j:
                tract_list.extend(option_j[opt_J])
        return tract_list

tract_list = make_tractlist(args.genome,args.JUNCTION,args.tract,args.junction)

#creating dictionary
D={}
for f in tract_list:
 print "loading polytract reference file " + f + " and creating dictionary, this may take a few minutes"
 try:
  file = gzip.open(f, "r")
 except (OSError, IOError) as e:
  print em1
  raise Exception(e)

 #filling dictionary with poly tracts
 for line in file:
  words=line.split(",")
  start=int(words[1])
  end=int(words[2])+1
  if args.region is not None and len(words)>4:
			R={"protein":"protein_coding",
					"pseudo":"pseudogene",
					"lncRNA":"lncRNA",
				"NA":"NA"}
			if words[4].strip("\n") != R[args.region]:
				continue
  for i in range(start, end):
   	 D[words[0]+"_"+str(i)]=words[3]
if len(words)<=4 and args.region is not None:
 print "WARNING: Due to lack of region annotation in the concerned genome, --region option was ignored and total tracts were considered."
 args.region = None
print "Identify tracts overlapping " + args.input
#open an output file
output = open("output/"+args.output, "w")

#open an invalid input file
invalid_input=open("output/"+args.output+".invalid_input.csv", "w")
if args.header=="1":
 next(input)
#merge input against dictionary
for line in input:
 words=line.split(d)
 #check to see if exmpty line exist in the input file
 l=len(words)
 if not line.strip():
  print "##########ERROR: The input file " + args.input + " contains empty lines in the middle or the end of the file, please fix this and rerun the program."
  print "Terminating the program!"
  exit(1)
 #check to see if input line contains at least 3 fields
 if l < 3:
  print "##########WARNING: Record "+line.strip("\n").strip('\r')+ " is an invalid input! Three fields (chr start end) are needed. Check if correct delimiter is specified with -d option"
  invalid_input.write(line)
  keyflag=2
  continue
 #compute boundary and check if chromosome positions are integers
 try:
  start=int(words[1])-boundary
 except ValueError:
  print "##########WARNING: Record "+line.strip("\n").strip('\r')+ " is an invalid input! Starting position is not an integer"
  invalid_input.write(line)
  continue

 try:
  end=int(words[2])+boundary+1
 except ValueError:
  print "##########WARNING: Record "+line.strip("\n").strip('\r')+ " is an invalid input! Ending position is not an integer"
  invalid_input.write(line)
  continue
 #check if start is before end
 diff=end-start
 if diff < 1:
  print "##########WARNING: Record "+line.strip("\n").strip('\r')+ " is an invalid input! Ending position is before starting position"
  invalid_input.write(line)
  continue
 #merging starts here
 #keyflag=0
 #print str(start)+ str(end)+str(diff)
 if diff>=1:
  tract_labels="" # initial empty string
  for i in range(start, end):
   key=words[0]+"_"+str(i)
   if key in D:
     tract_labels +=  D[key].strip("\n")+";"
  if tract_labels=="": 
   output.write(line.strip('\n').strip('\r') + d + "0\n")
  else:
   if args.intsect=="s":
    output.write(line.strip('\n').strip('\r') + d + tract_labels.split(';')[0]+"\n")	
   else:
    output.write(line.strip('\n').strip('\r') + d + tract_labels.strip(";")+"\n")
print "Overlapping examination of "+args.input + " is successful. Output is saved to output/"+args.output
print "If any invalid input exist, they are saved to output/"+args.output+".invalid_input.csv"

#close file handles
file.close()
input.close()
output.close()
invalid_input.close()
if os.stat("output/"+args.output+".invalid_input.csv").st_size==0:
 os.remove("output/"+args.output+".invalid_input.csv")

print "Conduct enrichment analysis on " + "output/" + args.output
print "If the R script completes normally, enrichment results are saved to " + "output/" + args.output + ".enrich and " + "output/" + args.output + ".tif"
#calling Rscript to perform statistical analysis
if args.junction is None:
	junCommand="breaks=NA" 
else:
	junCommand="breaks=paste0('break',1:"+str(args.junction)+")"
if args.region is not None:
	regionCommand="region='"+args.region+"'"
else:
	regionCommand="region=NULL" 
if args.JUNCTION is None:
	call(["Rscript","binomEnrich.R","joined='"+args.output+"'","gn='"+args.genome+"'","intMode='"+args.intsect+"'","arg.t="+args.tract,"b="+str(args.boundary),junCommand,"arg.d='"+ args.delimiter+"'",regionCommand])
else:
	JUNCommand="BREAKS=paste0('break',1:"+str(args.JUNCTION)+")"
	call( ["Rscript","binomEnrich.R","joined='"+args.output+"'","gn='"+args.genome+"'","intMode='"+args.intsect+"'",JUNCommand,"arg.d='"+ args.delimiter+"'",regionCommand] )
print "ALL DONE"
