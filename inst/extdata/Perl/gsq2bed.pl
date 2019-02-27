#!/usr/bin/perl
#
use strict;
use warnings;
use diagnostics;
use Getopt::Long;


my ($gsqfile,$outfile);

GetOptions (    "gsq=s"  => \$gsqfile,		# GeneSeqer-formatted transcript spliced alignment file
                "out=s"  => \$outfile);		# Output file name [if not specified, output will have .bed extension on input filename]

my $usage = "\nUSAGE:  gsq2bed.pl -gsq gsqfile [-out outfile]\n\n";

if (defined $gsqfile) { open (INFH, $gsqfile) or die ("Could not open $gsqfile file!!"); }
else { print $usage; exit; }


my (@oftmp,$outfstem);
if (defined $outfile) { }
else { 
  @oftmp = split (/\//,$gsqfile);
  $outfstem = $oftmp[$#oftmp];
  $outfile = $outfstem . ".bed";  
}
open (OUTFH, "> $outfile") or die ("Could not open $outfile!!");


my ($start,$est_id,$gseq_id,$id_string);
my (@tss,%est_ID,%gseq_ID,%gsq_orientation,%exon_beg,%exon_end) = ();

while(<INFH>) {
# ... read lines like the following:
#	Sequence    1:   4, from 1 to 300000, both strands analyzed.
 if (/^Sequence.*:\t*\s*(.+),\s*from.*$/) {
  $gseq_id = $1;
 }
# ... read lines like the following:
#	EST sequence      8 +strand  3727 n (File: PdomTSAr1.2-008727+)
#	EST sequence     12 -strand   942 n (File: 42517803-)
 if (/^EST sequence.*strand.*\(File:\s*(.*)([+-]+)\)$/) {
  $est_id = $1; 
 }

# ... read in all the 5'-first exons from the GeneSeqer output file:
#
 if (/^\s*Exon\s+1\s+(\d+)\s+(\d+)\s+\(\s*(\d+)\s+n\);\s+cDNA\s+(\d+)\s+(\d+)\s+\(\s*(\d+)\s+n\);\s+score:\s+(\d*\.\d+)\s*\n$/) {
#CHECK:
#	print "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
  
  if ($1 == $2) {next;}	# ... ignoring the edge-case of a 1-nt 1st exon ...
  $start = $1-1;
  $id_string = $start."_#_".$est_id;
   $gseq_ID{$id_string}          = $gseq_id;
   $exon_beg{$id_string}         = $start;
   $exon_end{$id_string}         = $2;
   $est_ID{$id_string}           = $est_id;
  if ($1 < $2) {
    $gsq_orientation{$id_string} = "+";
  }
  elsif ($1 > $2) {
    $gsq_orientation{$id_string} = "-"; 
  }
  push(@tss,$id_string);
 }
}


for(my $k=0;$k<=$#tss;$k++) {
  printf OUTFH "%s\t%d\t%d\t%s\t.\t%s\n",
    $gseq_ID{$tss[$k]}, $exon_beg{$tss[$k]}, $exon_beg{$tss[$k]}+1, $est_ID{$tss[$k]}, $gsq_orientation{$tss[$k]};
} 


close INFH;
close OUTFH;
