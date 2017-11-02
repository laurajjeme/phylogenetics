#!/usr/bin/perl
#
# Script name: tip2tip.pl
# Version 1.0 (Nov 14th, 2011)
# Author: Laura Eme
###################################################################################

use Bio::TreeIO;
use Statistics::Descriptive;

$stat = Statistics::Descriptive::Full->new();

$intree = $ARGV[0];

# parse in newick format
my $input = new Bio::TreeIO(-file   => $intree,
                            -format => "newick");
my $tree = $input->next_tree;

my @taxa = $tree->get_leaf_nodes;


#foreach $node1 (@taxa){ 
#	print $node1->id(),"\n"; 
#}

print " ------- Tree ongoing: ",$intree," -------\n";

#@matrix_distances = ();
@mean_distances = ();
@seq_names = ();

foreach $node1 (@taxa){ 

	$seq_name = $node1->id();
	push(@seq_names, $seq_name);

#	print $node1->id(),"\n";

	# list of distances is reset
	@distances = ();

	# mean long distance is reset
	$mean_dist = 0;

	foreach $node2 (@taxa){
		my $distance = $tree->distance(-nodes => [$node1,$node2]);
#		print $node1->id()," - ",$node2->id()," : ", $distance,"\n";
		push(@distances, $distance);
	}

	@distances = sort(@distances);

	# calculation of the average distance between a tip and the 5 other most distant tips
	for $i ($#distances-5 .. $#distances){
		$mean_dist = $mean_dist + ($distances[$i]);
	}
	$mean_dist = $mean_dist/5;
#	print "Dist moyenne : ", $mean_dist, "\n";
	push(@mean_distances,$mean_dist);

#	push(@matrix_distances,[@distances]);
}

$stat->add_data(@mean_distances);
$x = $stat->percentile(25);
$y = $stat->percentile(75);
$fs = $y - $x;
print "25 perc: ", $x , " - 75 perc: ", $y ,"\n";

$outlier_cutoff = 1.5*$fs;
$ext_outlier_cutoff = 3*$fs;

#print "Fs = ", $fs, " - 1.5 Fs = ", $outlier_cutoff , " - 3 Fs = ", $ext_outlier_cutoff, "\n\n";

for $i (0 .. $#mean_distances){
	$mean_dist = $mean_distances[$i];
	$seq_name = $seq_names[$i];

	if ($mean_dist > $y+(3*$fs)){
		print "Extreme outlier : ", $seq_name, " - Length: " ,$mean_dist,"\n";
	}
	elsif ($mean_dist > $y+(1.5*$fs)){
		print "Outlier : ", $seq_name, " - Length: ", $mean_dist,"\n";
	}
}

print "\n";
exit;


