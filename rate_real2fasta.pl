#!/usr/bin/perl
#
# Version 1.0
# Author: Laura Eme
#
###################################################################################

use lib '/local/one/SOFTWARE/bin/perllib';
use lauralib ;
use Bio::SeqIO::fasta;
use Bio::SearchIO;
use Bio::DB::Taxonomy;
use Scalar::Util qw(blessed);

$USAGE = "\n [USAGE] : rate2fasta.pl [Alignment FASTA format] [IQtree rate file] [Number of rate categories to create] [highest rate cat to keep]\n\nExample: rate2fasta.pl myfile.linsi myfile.rate 10 3" ;

# Print USAGE
unless( @ARGV )
{
   print $USAGE ;
    exit ;
}

$fasta_in = $ARGV[0];
@rates = read_file($ARGV[1]);
delete $rates[0];

$num_rate_cat = $ARGV[2];
$maxrate = $ARGV[3];

$searchio = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_in);

loop_A:while( my $seqobj = $searchio->next_seq ) {

    # sequence itself
    $sequence = $seqobj->seq();

    for $rate_line (@rates) {
        @tab = split(/\t/,$rate_line);
        $index = $tab[0]-1;
        $rate = $tab[1];
        $h_ind_rate{$index} = $rate;
        $h_ind_state{$index} = substr($sequence, $index, 1);
        #print "position ",$index," has a rate of ",$rate," and corresponds to the letter ",$h_ind_state{$index} ," in the first sequence \n";
   }
   last loop_A; #only need to do that on the first seq
}

# Calculate the number of site per category and how many sites to keep
$len_ali = length($sequence);
$num_site_per_cat = int(length($sequence)/$num_rate_cat);
#print "The alignment will be divided into",$num_rate_cat," rate categories containing ",$num_site_per_cat," sites each\n";
$num_sites_tokeep = $num_site_per_cat*$maxrate;
#print $num_sites_tokeep," sites will be kept\n";

$count = 0;
foreach my $key (sort { $h_ind_rate{$a} <=> $h_ind_rate{$b} } keys %h_ind_rate) {
     $count++;
     #print "count: ", $count, "\n";
     if ($count <= $num_sites_tokeep){
          $h_index_to_keep{$key} = 1;
          #print $key,"\n";
          #print $h_ind_rate{$key}," ",$key," ",$h_ind_state{$index}, "\n";
     }
}

#
$searchio = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_in);

loop_B:while( my $seqobj = $searchio->next_seq ) {

    # name of the sequence
    $id = $seqobj->display_id();
    print ">",$id,"\n";

    # sequence itself
    $sequence = $seqobj->seq();

    foreach $index (sort (keys %h_index_to_keep)){
         #print $index,"\n";
         print substr($sequence, $index, 1);
    }
    print "\n";
}
