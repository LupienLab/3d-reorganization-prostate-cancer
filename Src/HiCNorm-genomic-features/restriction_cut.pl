#!/usr/bin/perl  
use strict;

### only for 6 base cutter, since it hard-codes many things.... check code carefully
### dirty code but works fine. very small modification to yaffe's call.

## get file
if (scalar @ARGV != 3)
{print "Please give the <ref_genome> <fai file> <6base cut sequence>\n"; exit;}
my $file  = $ARGV[0]; chomp $file; open (FILE,"<$file");
my $flag = 5; my $seq; my $chr;
##

## read the fasta index file and get chromosome ends
my %hash_ends;
open (FIL,"<$ARGV[1]");
while (my $line = <FIL>)
{ chomp $line;
  my @liner = split('\t',$line);
  $hash_ends{$liner[0]} = $liner[1];
}
close FIL; 
##

## read the fasta file and search for restriction cut sites
while (my $line = <FILE>)
{ chomp $line; 

  if (substr($line,0,1) eq '>') 
  { 
    if ($flag == 5) {$flag = 0; $chr = $line;} # for chr1
    else {
           $seq =~ s/\n//g;
           $seq = uc($seq);
           my @position = (); my $i = 0;           
           while ($seq =~ m/$ARGV[2]/g)
           { my $dan = pos($seq)-4; ##sid
             #print "$chr\t$dan\n"; 
             $position[$i] = $dan; $i++;
             pos($seq) = pos($seq)-5;
           }
           
           my $j = 0;
           ## gets the gc and frag length
           foreach ($i = 0; $i < scalar @position; $i++)
           { $j++; my $c = 3; my $woo = $chr; $woo =~ s/>//g;  
             if (($i == 0) and ((scalar @position) == 1)) 
             {   if ($position[$i] >= 3)
		 { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
                   my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
                   print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]+$c; print "\t$ratio\n"; $j++;  
                 } 
  		   my $string = substr($seq,$position[$i],200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//); 
    	           my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
  		   print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $hash_ends{$woo}-$position[$i]; print "\t$ratio\n";    
             } # clear
             elsif ($i == 0)
             { if ($position[$i] >= 3)
               { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//); 
		 my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
                 print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]+$c; print "\t$ratio\n"; $j++; 
               }
  	         my $string = substr($seq,$position[$i],200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
                 my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
                 print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $position[$i+1]-$position[$i]+$c; print "\t$ratio\n"; 
             }
             elsif (($i > 0) and ($i < (scalar @position) - 1)) 
             { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
	       my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
               print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]-$position[$i-1]+$c; print "\t$ratio\n"; $j++; 

               $string = substr($seq,$position[$i],200); $c_a = ($string =~ tr/A//); $c_t = ($string =~ tr/T//); $c_g = ($string =~ tr/G//); $c_c = ($string =~ tr/C//);
               $den = ($c_c + $c_t + $c_a + $c_g); $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
               print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $position[$i+1]-$position[$i]+$c; print "\t$ratio\n";
             }
             else
             { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
               my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
               print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]-$position[$i-1]+$c;  print "\t$ratio\n"; $j++;

	       $string = substr($seq,$position[$i],200); $c_a = ($string =~ tr/A//); $c_t = ($string =~ tr/T//); $c_g = ($string =~ tr/G//); $c_c = ($string =~ tr/C//);
               $den = ($c_c + $c_t + $c_a + $c_g); $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
	       print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $hash_ends{$woo}-$position[$i]; print "\t$ratio\n";
             }
           }
           $seq = (); $chr = $line;
         }
  }
  else { $seq .= $line;}
}
close FILE;

$seq =~ s/\n//g;
$seq = uc($seq);  
my @position = (); my $i = 0; 
while ($seq =~ m/$ARGV[2]/g)
{ my $dan = pos($seq)-4; 
  $position[$i] = $dan; $i++;
  pos($seq) = pos($seq)-5;
}
           
my $j = 0;
foreach ($i = 0; $i < scalar @position; $i++)
{ $j++; my $c = 3; my $woo = $chr; $woo =~ s/>//g;  
             if (($i == 0) and ((scalar @position) == 1)) 
             {   if ($position[$i] >= 3)
		 { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
                   my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
                   print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]+$c; print "\t$ratio\n"; $j++;  
	       } 
		 my $string = substr($seq,$position[$i],200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//); 
		 my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
		 print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $hash_ends{$woo}-$position[$i]; print "\t$ratio\n";    
             } # clear
             elsif ($i == 0)
             { if ($position[$i] >= 3)
               { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//); 
		 my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
                 print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]+$c; print "\t$ratio\n"; $j++; 
	     }
	       my $string = substr($seq,$position[$i],200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
	       my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
	       print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $position[$i+1]-$position[$i]+$c; print "\t$ratio\n"; 
	   }
             elsif (($i > 0) and ($i < (scalar @position) - 1)) 
             { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
	       my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
               print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]-$position[$i-1]+$c; print "\t$ratio\n"; $j++; 

               $string = substr($seq,$position[$i],200); $c_a = ($string =~ tr/A//); $c_t = ($string =~ tr/T//); $c_g = ($string =~ tr/G//); $c_c = ($string =~ tr/C//);
               $den = ($c_c + $c_t + $c_a + $c_g); $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
               print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $position[$i+1]-$position[$i]+$c; print "\t$ratio\n";
	   }
             else
             { my $string = substr($seq,$position[$i] + $c - 200,200); my $c_a = ($string =~ tr/A//); my $c_t = ($string =~ tr/T//); my $c_g = ($string =~ tr/G//); my $c_c = ($string =~ tr/C//);
               my $den = ($c_c + $c_t + $c_a + $c_g); my $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
               print "$j\t-\t$woo\t"; print $position[$i]+$c; print "\t"; print $position[$i]-$position[$i-1]+$c;  print "\t$ratio\n"; $j++;

	       $string = substr($seq,$position[$i],200); $c_a = ($string =~ tr/A//); $c_t = ($string =~ tr/T//); $c_g = ($string =~ tr/G//); $c_c = ($string =~ tr/C//);
               $den = ($c_c + $c_t + $c_a + $c_g); $ratio;  if ($den > 0) {$ratio = ($c_c + $c_g)/$den; } else {$ratio = 0;} 
	       print "$j\t+\t$woo\t"; print $position[$i];    print "\t"; print $hash_ends{$woo}-$position[$i]; print "\t$ratio\n";
	   }
}
