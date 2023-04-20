#!/usr/bin/perl
#
# Script name: rm_short_seq.pl
# Version 1.0 (Nov 14, 2011)
# Author: Laura Eme
# Goal:
#
# Arguments:
###################################################################################

use lib '/local/one/SOFTWARE/bin/perllib';
use lauralib ;
use Bio::SearchIO;
use Bio::SeqIO::fasta;
use Bio::SeqIO;

$fasta_file = $ARGV[0];
@c60_file = read_file("/local/one/Dropbox/C60freq.txt");

### Run Ed's script ###
$command = "fasta2renamedForPhylip.pl ".$fasta_file." > ".$fasta_file.".renamed";
system($command);

$command = "fasta2phylip.pl ".$fasta_file.".renamed ".$fasta_file.".renamed.phy";
system($command);

$command = "/local/one/SOFTWARE/minmax-chisq/minmax-chisq -s33 -n1000000 -l4 -u4 -obinfile.check < ".$fasta_file.".renamed.phy";
system($command);

open(MINMAX, "binfile.check")||die;
while(<MINMAX>){
	chomp($_);

	open(OUT, ">C60recoded.nex")||die;
	print OUT  "#nexus\n";
	print OUT "begin models;\n";

	%h_cat = ();

	@tab = split(/ /,$_);
	$h_cat{"A"} = $tab[0];
	$h_cat{"R"} = $tab[1];
	$h_cat{"N"} = $tab[2];
	$h_cat{"D"} = $tab[3];
	$h_cat{"C"} = $tab[4];
	$h_cat{"Q"} = $tab[5];
	$h_cat{"E"} = $tab[6];
	$h_cat{"G"} = $tab[7];
	$h_cat{"H"} = $tab[8];
	$h_cat{"I"} = $tab[9];
	$h_cat{"L"} = $tab[10];
	$h_cat{"K"} = $tab[11];
	$h_cat{"M"} = $tab[12];
	$h_cat{"F"} = $tab[13];
	$h_cat{"P"} = $tab[14];
	$h_cat{"S"} = $tab[15];
	$h_cat{"T"} = $tab[16];
	$h_cat{"W"} = $tab[17];
	$h_cat{"Y"} = $tab[18];
	$h_cat{"V"} = $tab[19];

	$alphabet = 'AGNPSTCHWYDEKQRFILMV';
	$recoding = "";
	foreach $elem(@tab){
		$recoding = $recoding.$elem;
	}
	print $recoding,"\n";

	$_ = $recoding;
	eval "tr/0123/ATCG/, 1" or die $@;
	$recoding = $_;

	my $fasta_in = Bio::SeqIO->new(-file => $ARGV[0],-format => 'fasta');
	while (my $seq = $fasta_in->next_seq){
		$s = $seq->seq;
		$name = $seq->display_id();
	    #print "$s\n";
	    $_ = $s;
	    eval "tr/$alphabet/$recoding/, 1" or die $@;
	    $s = $_;

		$fasta_out = $fasta_file;
		$fasta_out =~ s/\.fas.*/_cr.fasta/;
		open(FASTA_REC,">>$fasta_out") || die;
	    print FASTA_REC ">",$name,"\n";
		print FASTA_REC $s, "\n";
	}

	$iter_cat = 0;
	foreach $line (@c60_file){
		#print $line;
		chomp($line);

		$iter_cat++;
		foreach $cat (sort(values %h_cat)){
			$h_addedfreq{$cat} = 0;
		}
		$sum_addedfreq = 0;

		%f_aa = ();

		@tab = split(/\t/,$line);
		$c60cat = $tab[0];
		$f_aa{"A"} = $tab[1];
		$f_aa{"R"} = $tab[2];
		$f_aa{"N"} = $tab[3];
		$f_aa{"D"} = $tab[4];
		$f_aa{"C"} = $tab[5];
		$f_aa{"Q"} = $tab[6];
		$f_aa{"E"} = $tab[7];
		$f_aa{"G"} = $tab[8];
		$f_aa{"H"} = $tab[9];
		$f_aa{"I"} = $tab[10];
		$f_aa{"L"} = $tab[11];
		$f_aa{"K"} = $tab[12];
		$f_aa{"M"} = $tab[13];
		$f_aa{"F"} = $tab[14];
		$f_aa{"P"} = $tab[15];
		$f_aa{"S"} = $tab[16];
		$f_aa{"T"} = $tab[17];
		$f_aa{"W"} = $tab[18];
		$f_aa{"Y"} = $tab[19];
		$f_aa{"V"} = $tab[20];

		foreach $aa (keys(%f_aa)) {
			$cat_aa = $h_cat{$aa};
			$h_addedfreq{$cat_aa} = $h_addedfreq{$cat_aa} + $f_aa{$aa};
			$sum_addedfreq = $sum_addedfreq + $f_aa{$aa};
			#print $c60cat, ":", $aa, ':', $cat_aa, ":",$h_addedfreq{$cat_aa} , "\n";
		}

		#print "total", $sum_addedfreq, "\n" ;
		#print $c60cat, "\n";
		print OUT "frequency C60NT".$iter_cat."=";
		foreach $cat (sort(keys (%h_addedfreq))){
			#print $cat, ":", $h_addedfreq{$cat}, "\n";
			print OUT " ".$h_addedfreq{$cat};
		}
		print OUT ";\n\n";
	}

	print OUT "model C60cust4=FMIX{C60NT1,C60NT2,C60NT3,C60NT4,C60NT5,C60NT6,C60NT7,C60NT8,C60NT9,C60NT10,C60NT11,C60NT12,C60NT13,C60NT14,C60NT15,C60NT16,C60NT17,C60NT18,C60NT19,C60NT20,C60NT21,C60NT22,C60NT23,C60NT24,C60NT25,C60NT26,C60NT27,C60NT28,C60NT29,C60NT30,C60NT31,C60NT32,C60NT33,C60NT34,C60NT35,C60NT36,C60NT37,C60NT38,C60NT39,C60NT40,C60NT41,C60NT42,C60NT43,C60NT44,C60NT45,C60NT46,C60NT47,C60NT48,C60NT49,C60NT50,C60NT51,C60NT52,C60NT53,C60NT54,C60NT55,C60NT56,C60NT57,C60NT58,C60NT59,C60NT60};";
	print OUT "\n\nend;";
	close(OUT);
}

close(MINMAX);

print "Now you can run IQtree: iqtree2-linked -m \"GTR+C60cust4+G4\" -mdef C60recoded.nex -link-exchangeabilities\n\n";
exit;
