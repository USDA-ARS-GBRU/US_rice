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


my $files = shift;
my $inDir = shift;
my $outDir = shift;
my $ref = shift;
my $outgroup = shift;
#my $outDir = shift;

#my $missingThreshold = .2; #each seq in subalignment, 70 percent must be non missing 

print("chrom\tindex\tstart\tend\tlength\trefStart\trefEnd\trefLength\ttype\tgapAlnCoverage\tconsensus\n");

my $buffer = 50;
my $lThresh = 50;
my $identThesh = 90;
my @lines = ('nip','cg',$outgroup);

open IN, "$files"; 
while (<IN>) {
	chomp $_;
	my $file = $_;
	my $chrom = $file;
	$chrom =~ s/\.xmfa//;
	my $refID = $ref.$chrom.".fa";
	my $index = 1;
	my $io = Bio::AlignIO->new(-file => "$inDir/$file", -format => "xmfa" );
	#my $io = Bio::AlignIO->new(-file => "aln2/Chr1_3367_3672_alt.aln", -format => "ClustalW" );
	while (my $aln = $io->next_aln()) {
		#print($aln->id()."\n");
		my $alnIndex = 1;
		my $alnLength = $aln->length();
		my $str = $aln->gap_line();
		while ($str =~ /\.(-{$lThresh,})\./g) {
			my $length = length($1);
			#next unless $length > 50;
			my $end = pos($str) - 1; 
			my $start = $end - $length + 1;
			my $lBuff = $buffer;
			my $rBuff = $buffer;
			#$lBuff = $start unless $start > $buffer;
			next unless $start > $buffer;;
			#$rBuff = $alnLength - $end unless $buffer < $alnLength - $end;
			next unless $buffer < $alnLength - $end;
			my $alnL = $aln->slice($start - 1 - $lBuff, $start - 1, 1); #1 keeps gaps
			next unless $alnL->num_sequences() == scalar(@lines);
			next unless $alnL->overall_percentage_identity() > $identThesh;
			my $alnR = $aln->slice($end + 1, $end + 1 + $rBuff, 1); #1 keeps gaps
			next unless $alnR->num_sequences() == scalar(@lines);
			next unless $alnR->overall_percentage_identity() > $identThesh;
			#print $alnL->overall_percentage_identity()."\t".$alnR->overall_percentage_identity()."\n";
			#warn "$refID";
			my $aln3 = $aln->slice($start, $end, 1);
			next unless $aln3->num_sequences() == scalar(@lines);
			my $consensus = $aln3->consensus_string();
			next if $consensus =~ /[Nn]{5,}/; 
			my %length;
			my %seq;
			foreach my $l (@lines) {
				my $ID = $l.$chrom.".fa";
				my $s = $aln3->get_seq_by_id($ID);
				#warn $s;
				$length{$l} = getLength($s);
				$seq{$l} = $s->seq();
			}
			my $ref = $aln3->get_seq_by_id($refID);
			my $refStart = $ref->start();
			my $refEnd = $ref->end();
			my $refLength = getLength($ref);
			my $type = 'noEventAssign';
			my $cov = "NA";
			if (getMatch($seq{'cg'}, $seq{'nip'}) > 95) {
				$cov = $length{"$outgroup"} / $aln3->length();
				$type = 'outgroup';
			} elsif (getMatch($seq{'cg'}, $seq{$outgroup}) > 95) {
				$cov = $length{'nip'} / $aln3->length();
				$type = "nip";# if  $cov > .95;
				#$type = 'nip_loss' if $cov < .05;
			} elsif (getMatch($seq{'nip'}, $seq{$outgroup}) > 95) {
				$cov = $length{'cg'} / $aln3->length();
				$type = "cg"; #if $cov > .95;
				#$type = 'cg_loss' if $cov < .05;
			}
			print("$chrom\t$alnIndex\_$index\t$start\t$end\t$length\t$refStart\t$refEnd\t$refLength\t$type\t$cov\t$consensus\n");
			#die Dumper(%length, %seq);
			if ($type =~ /cg|nip|noEventAssign/) {
				my $aln2 = $aln->slice($start - 1 - $lBuff, $end + 1 + $rBuff, 1);
				next unless $aln2->num_sequences() == scalar(@lines);
				my $out = Bio::AlignIO->new(-file => ">$outDir/$chrom\_$alnIndex\_$index.aln", -format => "ClustalW" );
				$out->write_aln($aln2);
			}
			$index++;
		}
		#die "test";
		$alnIndex++;
	}
	#die "test";
	
}

sub getLength {
	my $locSeq = shift;
	my $refStart = $locSeq->start();
	my $refEnd = $locSeq->end();
	my $refLength = 0;
	#print $ref->seq();
	#gaps have strange behavior; if seq is entirely gap, before and after are given
	#if there is some sequence then 1-index is used
	if ($locSeq->seq() =~ /\w/) {
		$refLength = abs($refEnd - $refStart) - 1 if $refStart > $refEnd; 
		$refLength = abs($refEnd - $refStart) + 1 if $refStart <= $refEnd;
	} 
	return $refLength;
}

sub getMatch { 
	my $a = shift;
	my $b = shift;
	$a =~ tr/acgt/ACGT/;
	$b =~ tr/acgt/ACGT/;
	my $count = 0;
	my $start = 0;
	while ($start < length($a)) {
		$count++ if substr($a, $start, 1) eq substr($b, $start, 1);
		$start++;
	}
	return ($count / length($a)) * 100	
	
}	

	# my $aln = shift;
# 	my $alnPair = $aln->select_noncont_by_name("$a$chrom.fa", "$b$chrom.fa");
# 	my $start = 1;
# 	my $length = $alnPair->length();
#
# 		my $aln2 = $aln->slice($start, $start, 1); #1 keeps gaps
# 		my %count;
# 		foreach my $seq ( $aln2->each_seq() ) {
# 		  $count{$seq->seq()}++;
# 		}
# 		my @sort = sort {$count{$a} <= $count{$b}} keys %count;
# 		#warn Dumper(%count, @sort);
# 		my $top = shift @sort;
# 		$seqCon .= $top unless ($top eq '-');
# 		$start++;
# 	}
# 	return $alnPair->overall_percent_identity();

			# my $skip = 0;
# 			foreach my $i (0 .. @lines - 1) {
# 				foreach my $j ($i + 1 .. @lines - 1) {
# 					my $alnPair = $aln3->select_noncont_by_name("$lines[$i]$chrom.fa", "$lines[$j]$chrom.fa");
# 					my $consensus = $alnPair->consensus_string();
# 					$skip = 1 if $consensus =~ /[Nn]{5,}/;
# 					print "$lines[$i]X$lines[$j]\t$consensus\n";
# 				}
# 			}
# 			next if $skip;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 SEE ALSO

	
=head1 AUTHOR


justin vaughn (jvaughn7@utk.edu)

=cut
