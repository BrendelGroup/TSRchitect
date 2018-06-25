#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;

#-------------------------------------------------------------------------------
## Copyright (c) 2012, Krishnakumar Sridharan, Iowa State University
## Copyright (c) 2018, Brendel Group, Indiana University
##
## Permission to use, copy, modify, and/or distribute this software for any
## purpose with or without fee is hereby granted, provided that the above
## copyright notice and this permission notice appear in all copies.
##
## THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
## WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
## MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
## ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
## WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
## ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
## OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
##-------------------------------------------------------------------------------

### TITLE --- gsq_gff_extractor.pl 
### AUTHOR --- Krishnakumar Sridharan 09/16/2012
### PURPOSE --- To associate 5'-Exon Start positions from GeneSeqer output with genes from GFF file. 
### USAGE --- perl gsq_gff_extractor.pl -gsq <GeneSeqer_output_file> -gff <GFF_annotation_file> -dist [Distance_from_annotated_gene_start_to_be_used_for_capturing_relevant_transcripts] -out [Output_file]
### PRE-REQUISITES --- "Getopt::Long" package available at CPAN

my ($infile_gsq,$infile_gth,$infile_gff,$dist,$outfile);

GetOptions (    "gsq=s" => \$infile_gsq,                # GeneSeqer output file. 
		"gth=s" => \$infile_gth,                # GenomeThreader output file.
                "gff=s" => \$infile_gff,               	# GFF format annotation file. Optional
		"dist=i" => \$dist,			# Distance from annotated gene start to be used for getting relevant transcripts. Optional
                "out=s" => \$outfile);             	# Output file name. Optional

my $usage = "\nUSAGE: perl gsq_gff_extractor.pl -gsq [GeneSeqer_output_file] -gth [GenomeThreader_output_file] -gff [GFF_annotation_file] -dist [Distance_from_annotated_gene_start_to_be_used_for_capturing_relevant_transcripts] -out [Output_file]\n
\t<>-Mandatory parameter\t\t[]-Optional parameter\n
NOTE: Please define all mandatory paramters and install Getopt::Long perl module from CPAN before proceeding\n";

if(defined $infile_gsq || defined $infile_gth) {} else { print $usage; exit; }
#if(defined $infile_gff) {} else { print $usage; exit; }   ##### Show usage if cmd-line parameters are undefined
if(defined $dist) {} else { $dist=500; }   ####### Default distance to be used for annotation = [-500,+500]

if(defined $infile_gsq) { open (INFILE_TRANS, $infile_gsq) or die ("Could not open $infile_gsq file!!"); }
elsif(defined $infile_gth) { open (INFILE_TRANS, $infile_gth) or die ("Could not open $infile_gth file!!"); }
if (defined $infile_gff) { open (INFILE_GFF, $infile_gff) or die ("Could not open $infile_gff file!!"); }

####### Setting defaults for cmd-line parameters ################
my (@outname_temp,$outlabel);
if (defined $outfile) { }
else 
{ 
 @outname_temp = split (/\//,$infile_gsq);
 $outlabel = $outname_temp[$#outname_temp];
 $outfile = $outlabel . ".bed";  
}
open (OUTFILE, "> $outfile") or die ("Could not open $outfile!!");

my ($est_id,$est_orientation,$chr_id,$exon_start,$exon_start2, $exon_end,$exon_orientation,$species_chr_id,$gene_orientation,$gene_start,$gene_end,$gene_id,$gene_chr);
my (@transcriptSS,@exon_positions,@gff_coords,%est_ID,%est_orientation,%est_chrID,%gene_ID,%id_to_allstarts,%gsq_orientation,%gene_orientation,%exon_end,%gene_chrID) = ();

while(<INFILE_TRANS>)
{
#Sequence    1:   4, from 1 to 300000, both strands analyzed.
 if(/^Sequence.*:\t*\s*(.+),\s*from.*$/)
 {
  #print "$1\n";
  $chr_id = $1;
 }
# EST sequence     12 -strand   942 n (File: 42517803-)
# EST sequence      8 +strand  3727 n (File: PdomTSAr1.2-008727+)
 if(/^EST sequence.*strand.*\(File:\s*(.*)([+-]+)\)$/)
 {
  #print "$1\n";
  $est_id = $1; 
  $est_orientation = $2;
 }
 if(/^\s*Exon\s+1\s+(\d+)\s+(\d+)\s+\(\s*(\d+)\s+n\);\s+cDNA\s+(\d+)\s+(\d+)\s+\(\s*(\d+)\s+n\);\s+score:\s+(\d*\.\d+)\s*\n$/)  ##### Read in all the 5'-First Exons from the GeneSeqer output
 {
  #print "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
  if($1>$2)
  {
   $exon_start = $1-1;
   $exon_start2 = $1;
   $exon_end = $2;
   $exon_orientation = "-";
   my $id_string = $exon_start."_CustomSep_".$exon_start2."_CustomSep_".$est_id;	 #### Multiple ESTs can match at one position, or at multiple positions as well. This complicates the hash key.
   $gsq_orientation{$id_string}=$exon_orientation; 
   $exon_end{$id_string} = $exon_end;
   $est_ID{$id_string} = $est_id;
   $est_orientation{$id_string} = $est_orientation;
   $est_chrID{$id_string} = $chr_id;
   push(@exon_positions, $exon_start);
   push(@transcriptSS,$id_string);
  }
  elsif ($2>$1)
  {
   $exon_start = $1-1;
   $exon_start2= $1;
   $exon_end = $2;
   $exon_orientation = "+";
   my $id_string = $exon_start."_CustomSep_".$exon_start2."_CustomSep_".$est_id;	 #### Multiple ESTs can match at one position, or at multiple positions as well. This complicates the hash key.
   $gsq_orientation{$id_string} = $exon_orientation; 
   $exon_end{$id_string} = $exon_end;
   $est_ID{$id_string} = $est_id;
   $est_orientation{$id_string} = $est_orientation;
   $est_chrID{$id_string} = $chr_id;
   push(@exon_positions, $exon_start);
   push(@transcriptSS,$id_string);
  }  
 }
}


if (defined $infile_gff)
{ 
 while(<INFILE_GFF>)
 {
  if(/^[a-zA-z_]*(\d+)\s+([A-Z_a-z0-9]+)\s+(gene)\s+(\d+)\s+(\d+)\s+\.\s+([+-]*).*ID=([A-Za-z0-9_]+);.*$/)  #### GFF file format, check for compatibility 
  {												#### Read in all the gene start, end and IDs from gff file
    #print "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
    $gene_chr = $1;
    $gene_orientation = $6;
    $gene_start = $4;
    $gene_end = $5;
    $gene_id = $7;
    push(@gff_coords,$gene_start);
    push(@gff_coords,$gene_end);
    $gene_orientation{$gene_start}=$gene_orientation;
    $gene_chrID{$gene_start} = $gene_chr;
    $gene_ID{$gene_start} = $gene_id;
  }
 }

 for(my $i=0;$i<=$#transcriptSS;$i++)      ####### Associate 5'-transcript ends with gene IDs based on distance from 5'-annotated gene start and orientation
 {
  for (my $j=0;$j<=$#gff_coords;$j+=2)
  {
   if($gene_orientation{$gff_coords[$j]} eq '+')
   {
       my @start_and_id = split(/_CustomSep_/, $transcriptSS[$i]);
       #this is where I need to make the shift
    if($start_and_id[0]>($gff_coords[$j]-$dist-1) && $start_and_id[0]<($gff_coords[$j]+$dist+1) && $est_chrID{$transcriptSS[$i]} eq $gene_chrID{$gff_coords[$j]})
    {
     print OUTFILE "$est_chrID{$transcriptSS[$i]}\t$start_and_id[0]\t$start_and_id[0]\tgi|$est_ID{$transcriptSS[$i]}\t.\t$est_orientation{$transcriptSS[$i]}\n";
     last;
    }
   }
   elsif($gene_orientation{$gff_coords[$j]} eq '-')
   {
    my @start_and_id = split(/_CustomSep_/, $transcriptSS[$i]);
    if($start_and_id[0]>($gff_coords[$j+1]-$dist-1) && $start_and_id[0]<($gff_coords[$j+1]+$dist+1) && $est_chrID{$transcriptSS[$i]} eq $gene_chrID{$gff_coords[$j]})
    {
     print OUTFILE "$est_chrID{$transcriptSS[$i]}\t$start_and_id[0]\t$start_and_id[0]\tgi|$est_ID{$transcriptSS[$i]}\t.\t$est_orientation{$transcriptSS[$i]}\n";
     last;
    }
   }
  }
 }
} 		## if defined infile_gff loop
else
{
 for(my $k=0;$k<=$#transcriptSS;$k++)
 {
  my @start_and_id = split(/_CustomSep_/, $transcriptSS[$k]);
  #print "$start_and_id[0]\t$start_and_id[1]\n"; 
  print OUTFILE "$est_chrID{$transcriptSS[$k]}\t$start_and_id[0]\t$start_and_id[1]\tgi|$est_ID{$transcriptSS[$k]}\t.\t$gsq_orientation{$transcriptSS[$k]}\n";
# print OUTFILE "$est_chrID{$transcriptSS[$k]}\t$start_and_id[0]\t$start_and_id[1]\tgi|$est_ID{$transcriptSS[$k]}\t.\t$est_orientation{$transcriptSS[$k]}\n";
 } 
}
close INFILE_TRANS;
if (defined $infile_gff) { close INFILE_GFF; }
close OUTFILE;
