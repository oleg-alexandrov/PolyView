#!/usr/bin/perl
use strict;        # insist that all variables be declared
use diagnostics;   # expand the cryptic warnings
undef $/; # undefines the separator. Can read one whole file in one scalar.
MAIN:{

  open(FILE, "<../License.txt");
  my $license = join("", <FILE>);
  close(FILE);
  $license =~ s/\s*$//g;
  $license =~ s/(^|\n)/$1\/\/ /g;
  $license .= "\n";
  print "$license\n";
  
  my @files = (<*.cpp>, <*.h>);
  foreach my $file (@files){
    open(FILE, "<$file");
    print "$file\n";
    my $text = join("", <FILE>);
    $text = $license . $text;
    open(FILE, ">$file");
    print FILE $text;
    close(FILE);
  }
  
}
