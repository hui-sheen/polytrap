## Introduction
*Polytrap: screening genomic features trapped by polytracts (single, di-, tri-nucleotide tandem repeats)*

Polytrap is useful for screening genomic features over-represented (“trapped”) by single, di-, or tri-nucleotide tandem repeats (in our terminology: polytracts). A polytract is identified if the tandem repeat meets the total length requirements, i.e., 6nt for single and di-nucleotide tracts and 9nt for tri-nucleotide tracts. Polytrap stores pre-identified polytract ranges, against which the possible enrichment of a specific genomic feature is assessed through the Binomial probability model. The genomic feature must be given with BED-like locations. I.e., the input file must contain such rows as "X,93712222,93712224". Currently Polytrap can handle 11 genomes of 9 species (human, macaque, mouse, rat, dog, chicken, fugu, fruitfly, and yeast).

## Download & Deploy
Download the polytrap package from github:

	git clone https://github.com/hui-sheen/polytrap/

Because the tracts files for 9 species exceeded repository quota set by GitHub, we put all tract files elsewhere (https://drive.google.com/open?id=16kpdgronJbrJjuRyS2AsXIvGPCogP0ns). You MUST download the tract files, and place them (*.csv.gz files) at polytrap/tracts/.

## Quickstart
	cd polytrap # the root directory of polytrap
	python polytrap.py -i in.bed -o out 

## Input
You must prepare an input file containing genomic ranges with the first three columns for chromosome, start position, and end position. We suggest putting this input file at the first-level directory (*polytrap/*). Example input files "in.bed" and "in.csv" have been included in the tarball at the first-level directory. Each row in an input file corresponds to a genomic range of an instance of the interested genomic feature. In just a few minutes*, the program will generate the following 4 output files at *polytrap/output*. 

## Output
1) out 

This is the foremost output file indicating the overlapping situation between user-given ranges and curated tracts. It has the same number of rows as in.bed. The beginning part of each row is exactly the same as in in.bed. The last cell in each row has the value "0" if the range does not overlap a tract, otherwise the name of an overlapping tract (e.g., "gc").

2) out.enrich
 
This contains the enrichment analysis statistics. Each row is for a particular type of polytracts, and there might be a bottom row ("Overall") for all groups combined. 

Column Name | Meaning
------------|--------
*nFeatures_intract* | number of ranges overlapping with a tract
*nucTract* | total number of nucleotides occupied by a polytract type
*nFeatures_ingenome* | number of user-supplied ranges (i.e., number of rows in input)
*nucGenome* | total number of nucleotides in the particular genome (default HG38)
*pEnrich* | P value calculated under Binomial distribution model
*obsRate* | *nFeatures_intract/nucTract*
*expRate* | *nFeatures_ingenome/nucGenome*


3) out.tif

This figure file translates enrichment statistics to visual display. It has two vertically stacked panels. The top panel is a barplot for *RR* (*obsRate/expRate*) values of each polytract type as well as the combined polytract set; the bottome panel is a piechart with slices representing the individual tract types. The slice size is proportional to *RR* values and the color scale is proportional to enrichment *p*. An asterisk indicates statistical significance (p<0.01).
![out.tif](/output/out.jpg)

4) out.landscape.tif

The enrichment p value for the Overall row of output file *out.enrich* is taken to draw a red bar in a landscape barplot, where pre-calculated p values for nearly 100 genomic features are depicted as a reference background. In this landscape barplot, p=1e-4 is indicated as a bonferroni corrected significance threshold and all gnomic features exceeding this threshold are labelled in blue text. 
![out.tif](/output/out.landscape.jpg)

## Arguments
Two mandatory arguments are -i (--input) and -o (--output). So users must prepare one input file containing genomic ranges. This file should contain three columns, denoting chromosome, start position, and end position, respectively. The field separator can be comma or tab. Please refer to the example input file included in the package (in.bed and in.csv). Besides, the user must also indicate the file name for the output, given as argument -o (--output).

Other optional arguments pertain to genome (-g), tract type (-t), extension or boundary (-b), hinge or junction (-j and -J), genomic region constraint (-r), intersection mode (-I), and input file specification (-H and -d). Type the following command for a comprehensive help on these options.

	python polytrap.py --help

Let us explain a bit more about hinge (-j) and intersection mode (-I). In a genome, two tandem stretches of polytracts may be separated by only one (or 1~3) nucleotide, in which case we term the separating nucleotide(s) a *hinge* site. We suggest considering 1-nt hinges, although this can be tunable within [1,3]. Intersection mode designates how an overlap between a user-given genomic interval and a polytract is defined and quantified. Under the default *singleton* mode (*-I s*), we regard a genomic interval as a singleton unit, so whenever a spatial overlap appears, we count it as ONE overlap. Under the alternative *multiplex* mode (*-I m*), we regard a genomic interval as a union of its constituent nucleotides, so when a spatial overlap appears, we take into account the number of overlapping nucleotides.  

## Extension to new genomes
<details>
	<summary>Polytrap can be extended to uncovered genomes. Click to find detailed instruction!</summary>

Assume Polytrap does not cover canFam3 (dog), and you want to incorporate this new genome into Polytrap. You can first move away all \*canfam3\* files at polytrap/tracts, and follow the instructions to test if those files are being generated.

#### Prerequisites
	<keyname> canFam3 or canfam3 (lower-case)
	<genomes.meta> need to be updated
	[optional] <chr2NC.txt> need to be updated
	[optional] <polytrap/genomes/> need to be populated with NCBI's GFF file for the new genome
#### Expected output: under polytrap/new/newtracts/, the following files will be generated
	polytrap_canfam3_single.csv
	polytrap_canfam3_di.csv
	polytrap_canfam3_tri.csv
	polytrap_breaks_canfam3_1.csv
	polytrap_breaks_canfam3_2.csv
	polytrap_breaks_canfam3_3.csv
	polytract_canfam3_stats.txt
	
#### [optional] Step 0: Prepare well-formatted GFF file for annotating genomic region for tracts.
###### 0.1 Download GFF file from the species-specific page at NCBI. Make sure the file name contains "gff".

For canFam3, go to https://www.ncbi.nlm.nih.gov/genome/?term=txid9615[Organism:noexp]. You can see a hyperlinked word "GFF", which links to ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1/GCF_000002285.3_CanFam3.1_genomic.gff.gz. 

###### 0.2 Simplify the raw GFF file. A simplified GFF file named \*simGff\* will appear under polytrap/new/genomes/.
	cd polytrap/new/
	Rscript simplifyGFF.R # this takes ~1m for a GFF file of ~30M size.
###### 0.3 Augment file polytrap/new/chr2NC.txt with chromosome mapping for the new genome. Each new line must contain three tab-separated cells. Make sure "chr" are all lower-case letters.
Chromosome mapping information is available at the bottom of NCBI's species-specific page. This information maps NC_ numbers to more intuitive chromosome names. For example, NC_006583.3 maps to canFam3's chr1. 

#### Step 1: Append genomes.meta file with one genome-specific row.
Add one more row to the bottom of file polytrap/new/genomes.meta. The comma-separated values correspond to the following sequentially: *genome keyname* ("canFam3"), *Genome nucleotide total* (2410976875), *BSgenome package name* ("BSgenome.Cfamiliaris.UCSC.canFam3"), *simplified GFF file name* ("ref_CanFam3.1_top_level.simGff3", the file generated in Step 0). If no simplified GFF is available, put "NA" instead.

Genome nucleotide total can be Ensembl’s golden path length. For canFam3, this number is available at https://uswest.ensembl.org/Canis_familiaris/Info/Annotation. BSgenome package should be included in the list http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData. For canFam3, the package name is "BSgenome.Cfamiliaris.UCSC.canFam3"  

#### Step 2: Generate tract and hinge files for the new genome. Genomic region information will be annotated for the each tract and hinge if Step 0 has been done beforehand.
	cd polytrap/new/
	./newtracts.sh canFam3 &>log/newtracts.log& #canFam3 can be replaced with the keyname of your specific new genome.
For canFam3, a genome of ~2.4G nucleotides, it took ~75 minutes* to generate tract files and ~75 minutes* to annotate genomic region for these tracts.

#### [optional] Step: Test if polytrap works well with the newly incorporated genome (canFam3) 
	cd polytrap # the root directory of polytrap
	python polytrap.py -i in.bed -o canfam3.test -g canFam3

*Time estimation was made on such a Ubuntu server: Intel(R) Xeon(R) CPU E5-2650 v4 @ 2.20GHz, 128G memory, 4 hard disks of 4TB WDC WD40EZRZ-75G 
</details>
