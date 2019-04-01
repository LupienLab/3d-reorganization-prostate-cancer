#!/usr/bin/perl
use strict;

### makes the final all file.
###

if (scalar @ARGV != 3)
{ print "Please provide <all> <seqs> and <sam>\n"; exit;}

open (ALL,"<$ARGV[0]");
open (SEQ,"<$ARGV[1]");
open (SAM,"<$ARGV[2]");

my %hash_seq;
while (my $line = <SEQ>)
{ chomp $line;   
  next, if (substr($line,0,1) ne '@');

  my @liner = split(':',$line);  
  my $dan = $liner[1].'_'.$liner[3];   

  $hash_seq{$dan}++;
}
close SEQ;

my %hash_sam;
while (my $line = <SAM>) 
{ chomp $line; 
  my @liner = split('\t',$line);

  if ($liner[4] > 10)
  { my @danner = split(':',$liner[0]); 
    my $dan = $danner[1].'_'.$danner[3];

    $hash_sam{$dan}++;   
  }
}
close SAM;

while (my $line = <ALL>)
{ chomp $line; 
  my @liner = split('\t',$line);
  my $dan = $liner[2].'_'.$liner[3];

  my $obs = 0;
  $obs = $hash_sam{$dan}, if (defined $hash_sam{$dan});
  my $exp = $hash_seq{$dan};

  my $rat = 0;
  $rat = $obs/$exp, if ($obs > 0);
  print "$line\t$obs\t$exp\t$rat\n";
}
close ALL;  
