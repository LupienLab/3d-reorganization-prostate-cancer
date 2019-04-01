#!/usr/bin/perl
use strict;

## the script takes the all file and fasta file. It takes 55 reads for every restriction cut sites in both up an down stream. 55 reads because YT did 55 reads ##

## initialize variables
my ($fasta_file, $filename) = @ARGV;
if ((not defined $fasta_file) || (not defined $filename))
    { die ("Usage: ./fastq_from_loc.pl <fasta file> <all*.F_GC file>\n");}

open (FILE,"<$filename");
my @cuts = <FILE>; close FILE; chomp @cuts;
my $out = $filename.'_seqs'; open(FILE,">$out");

open (FASTA,"<$fasta_file"); my %genome; my $chro;
while (my $line = <FASTA>)
{ chomp $line;
  
  if (substr($line,0,1) eq '>')
  {$chro = $line; $chro =~ s/>//g;}
  else
  {$genome{$chro} .= $line;}
}
close FASTA;

foreach (@cuts)
{  chomp $_;

   my @liner = split('\t',$_); 
   my $query;
   my $chr   = $liner[2];
   my $str   = $liner[1];
   my $pos   = $liner[3];
   my $geno = $genome{$chr}; 

   my $flag = 0;
   if ($str eq "+") 
   { $query  = substr($geno,$pos-1,549);}
   elsif ($str eq "-") 
   { if ($pos > 600) {$query   = substr($geno,$pos - 549,549);}
     else            {$query   = substr($geno,0,$pos); $flag = 1;}
   }

   my $end = 0;
   if ($flag == 0) {$end = 496;}
   else            {$end = $pos;}

   my $string = $query; $string = uc($string); my $qual   = '>' x 36;
   my @new_array; $geno = (); 
 
   if ($str eq "+") 
   {
     for (my $i = 1; $i < 496; $i += 9) 
     {  my $seq  = substr($string,$i,36); my $lpos = $i + $pos; 
        print FILE "\@test_fastq:$chr:$lpos:$pos\n$seq\n+\n$qual\n";
     } 
   }
    elsif ($str eq "-") 
    { my $lstring = reverse($string); 
      for (my $i = 1; $i < $end; $i += 9) 
      { my $seq  = reverse(substr($lstring,$i,36)); my $lpos = $pos - 35 - $i;
        if (length($seq) < 36) {last;}
        print FILE "\@test_fastq:$chr:$lpos:$pos\n$seq\n+\n$qual\n";   
      }
    }
}
close FILE;
