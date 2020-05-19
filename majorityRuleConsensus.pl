#!/usr/bin/perl -w

#
use strict;
use warnings;
use Data::Dumper;
#program version
my $VERSION="0.1";

use Bio::AlignIO;
use Bio::Seq;
use Bio::SeqIO;


my $files = "/Users/justinvaughn/Desktop/historicUSRice/repeatsHistoric/clusterAll/files.txt";
my $inDir = "/Users/justinvaughn/Desktop/historicUSRice/repeatsHistoric/clusterAll/alignments";
my $outDir = "/Users/justinvaughn/Desktop/historicUSRice/repeatsHistoric/clusterAll/consensus";

my $missingThreshold = .2; #each seq in subalignment, 70 percent must be non missing 

open IN, "$files"; 
while (<IN>) {
	chomp $_;
	my $file = $_;
	my $io = Bio::AlignIO->new(-file => "$inDir/$file", -format => "fasta" );
	my $seqCon = '';
	my $aln = $io->next_aln();
	#print($aln->id()."\n");
	my $start = 1;
	my $length = $aln->length();
	warn ("############ alignment $length ##################");
	while ($start <= $length) {
		my $aln2 = $aln->slice($start, $start, 1); #1 keeps gaps
		my %count;
		foreach my $seq ( $aln2->each_seq() ) {
		  $count{$seq->seq()}++;
		}
		my @sort = sort {$count{$a} <= $count{$b}} keys %count;
		#warn Dumper(%count, @sort);
		my $top = shift @sort;
		$seqCon .= $top unless ($top eq '-');
		$start++;
	}
	my $o = Bio::SeqIO->new(-file => ">$outDir/$file", -format => 'fasta');
	my $newSeq = Bio::Seq->new(-id => "$file\_consensus", -seq => $seqCon);
	$o->write_seq($newSeq);
}

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 SEE ALSO

	
=head1 AUTHOR


justin vaughn (jvaughn7@utk.edu)

=cut
