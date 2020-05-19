# US_rice scripts/programs

Although anyone is free to modify and use this software in this repository, it is not meant for general application.  It has been deposited here to explicitly define key analytical components of the manuscript, "Gene disruption by structural mutations drives selection in US rice breeding over the last century".  To that end, many of the input files are hard-coded and are distributed as supplementary information associated with the publication.  It there is community interest in a particular script, please create an "Issues" thread and we can discuss expanding the tool. 

## inferInsertionOrDeletion.pl

Using whole-chromosome, multiple sequence alignments (in xmfa format) between CarGold, Nipponbare, and O. glaberrima, the script extracts any sequence (>50) of columns in which any row contains a gap.  The script then assumes that the flanking region and the deleted sequence (if applicable) pass a certain identity threshold.  The script also produces subalignments in a ClustalW format of the relevant region.

## DStatFromBlst2Seq.pl

Produces the d metric from blast2seq results (now produced using blastp) between derived and ancestral sequence.  In our case, derived and "ancestral" sequences will be either CarGold and Nipponabare or vice-versus.

## bedCovCommandForSegIndels.sh

Shell commands used to define coverage around pre-ascertained indel positions.  If the event w a deletion relative to the reference, then start and end coordinates were used.  Alternatively, if the event was an insertion relative to the reference, then the position of the insert +/- 5 bp was used.  This pre-processing step helped address the noise present when estimating coverage over a 1 bp interval.

## populationClustering.R

R program that takes a genetic relatedness matrix (in this case centered-IBS), clusters varieties/strains based on the hypothesis of speciation (or population divergence), and compiles the results in a sliding window of 25 varieties/strains through time.

## majorityRuleConsensus.R

Perl program to calculate the majority rule consensus from a clustalw alignment (although other alignment formats are easily imported via bioperl).  The approach is different than others in that it accounts for gap consensus and compression as well. The script takes a list of alignments to process and the alignement directory as input and creates a fasta file for every alignment in the specified output directory.  

