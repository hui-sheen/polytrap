#!/usr/bin/env perl
# tract_anno_type.pl hg38 ../genomes/GCF_000001405.38_GRCh38.p12_genomic.simGff
# $ARGV[0]: hg keyword, say hg38,hg19,mm10.
# $ARGV[1]: genome GTF file (simplified)
# OUTPUT will be the same file name as tract file, just with csv changed to tsv.
# UPDATE 9/21 relative to 9/18 *.dir.pl: Building site hash on each individual tracts, loop over Gff sites per interval.
# Originally named as tract_anno_type.tractHash.pl
use List::MoreUtils qw(uniq);
opendir(DIR,'./newtracts'); # UPDATE 10/5/2018
print STDERR "Repeat geomic annotation for all CSV files on $ARGV[0] in $ENV{'PWD'}/newtracts/\n";
@files=grep(/.*$ARGV[0].*\.csv$/,readdir(DIR));
foreach $fi (@files) {
	open(TRACTin,"./newtracts/$fi") || die "Reading tract file $fi fails!!!\n";
	%D = ();
	while (<TRACTin>) {
		chomp($_);
		@fields=split(",",$_);
		for ($i=$fields[1];$i<=$fields[2];$i++) {
			$D{$fields[0]."_".$i}=$_;
		};
	};
	%revD = ();
	push @{$revD{$D{$_}}},$_ for keys %D;
	close TRACTin;
    open(TRACT,">./newtracts/$fi") || die "Cannot write annotated tracts to ./newtracts/$fi!!!\n"; # To replace the file just read in.
	open(GTF,"$ARGV[1]") || die "Reading genome GTF file fails!!!\n";
    $datestring=localtime();
	print STDERR "$datestring\t Annotating genomic region for tract file $fi...\n";
	while (<GTF>) {
		next if $.<2;
		if (/^chr(\S+)\s  (\d+)\s  (\d+)\s  (\S+)    /x) {
			$chr=$1;
			$start=$2;
			$end=$3;
			$fType=$4;
			for ($j=$start;$j<=$end;$j++) {
				my $key=$chr."_".$j;
				if ( exists($D{$key}) ) {
					print TRACT "$D{$key},$fType\n";
					@synKeys = @{$revD{$D{$key}}};
					#$synonyms = join '\t',@synKeys;
					delete @D{@synKeys};
				};
			};
		};
	};
	close GTF;
	if ( scalar(keys(%D))>0 ) {
		@remaining = uniq(values(%D));
		foreach $tract (@remaining) {
			print TRACT "$tract,NA\n";
		};
	};
	close TRACT;
};
closedir(DIR);
print STDERR "Job done! newtracts/*$ARGV[0]*.csv got updated with genomic region annotation!!!\n";
