#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

my $SV_file = shift;
my $dir = shift;

# my $coverageThresh = 95;
# my $interveningDNAThreshParam = .12;
# my $interveningDNAThresh = 2;
# my $fusionThreshParam = .12;

my %events;
my %length;

open IN, "$SV_file";
while (my $line = <IN>) { # read blast file
    chomp $line;
    my @s = split /\t/, $line;
    next if $s[8] eq 'outgroup';
    next if $s[8] eq 'noEventAssign';
	next if ($s[9] < .51 && $s[9] > .49); #gets rid of weird gap opening mauve does 
	$events{"$s[0]\_$s[1]"} = "$s[8]\_$s[9]";
	#$length{"$s[0]\_$s[1]"} = $s[4] * $s[9]; #length of sequenc
	$length{"$s[0]\_$s[1]"} = $s[4]; #length of event
}
close IN;

opendir DIR, $dir;
my @files = readdir(DIR);
closedir DIR;

foreach my $file (@files) {
  next unless $file =~ /\.blast/;
  next if $file =~ /^\./;
  my ($id, $ext) = split /\./, $file;
  #get the first two (two longest hits);
  next unless $length{$id};
  my $l = $length{$id}; #rough approximation for check below
  my $end;
  my $begin;
  my @diag;
  open IN, "$dir/$file";
  my $i = 1;
  while ($i < 3) {
    my $line = <IN>;
    chomp $line;
    next if $line =~ /^\#/;  # ignore comments
    my($Query_id, $Subject_id, $identity, $alignment_length, $mismatches, $gap_openings, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split /\t+/, $line;
    last if $alignment_length < (.8 * $l); #alignment buffer used +- indel length, so each main diagonal must be nearly the same size as indel 
    $diag[$i] = [$qstart, $sstart, $send];
    $i++;
  }
  next if $i < 3; #catches length check
  my ($topD, $botD);
  if ($diag[1]->[0] < $diag[2]->[0]) {
	$topD = $diag[1];
	$botD = $diag[2];
    # $end = $diag[1]->[2];
#     $begin = $diag[2]->[1];
  } else {
  	$topD = $diag[2];
  	$botD = $diag[1];  
    # $end = $diag[2]->[2];
    # $begin = $diag[1]->[1];;
  }
  my $q_Bot_start = $botD->[0];
  my $q_Top_start = $topD->[0];
  my $s_Bot_start = $botD->[1];
  my $s_Top_start = $topD->[1];
  $l = ($q_Bot_start - $s_Bot_start + $s_Top_start) - $q_Top_start; #base it directly on Trace
  # unless ($begin && $end) {
  #   warn "problem: begin equals $begin, end equals $end for $id\n";
  #   next;
  # }
  my $s_Top_end = $topD->[2];
  my $d = $s_Top_end - $s_Bot_start - 1;
  #my $d = $begin - $end - 1; #-1 because you want the overlap 
  my @temp = split /\_/, $events{$id};
  print "$id\t$temp[0]\t$temp[1]\t$l\t$d\n"
}

#author: Justin Vaughn (jnvaughn@uga.edu)