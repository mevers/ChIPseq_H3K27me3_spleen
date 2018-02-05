#!/usr/bin/perl

# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Last changed: 26-07-2017

use warnings;
use strict;

my $version = 0.9;

usage() if (scalar(@ARGV) != 2);

my $dir = $ARGV[0];
my $filter = $ARGV[1];
opendir(my $dh, $dir) or doe $!;

my %ht = ();
while (my $file = readdir($dh)) {
    if ($file =~ /$filter/) {
	my $id = $file;
	$file = join("/", $dir, $file);
	open(my $fh, $file) or die $!;
	while (my $line = <$fh>) {
	    if ($line =~ /PERCENT_DUPLICATION/) {
		my @arr = split("\t", <$fh>);
		$ht{$id} = [@arr];
	    }
	}
	close($fh);
    }
}

my @cn = ("Sample", "Read pairs", "Duplicate pairs", "Percentage dupes");
printf("%s\n", join("\t", @cn));
foreach my $key (sort {$a cmp $b} keys %ht) {
    my @arr = @{$ht{$key}};
    printf("%s\n", join(
        "\t",
        $key,
        commify($arr[2]),
        commify($arr[6]),
        sprintf("%3.2f", $arr[8] * 100)));

}


sub usage {
    printf("summarise_DuplicationMetrics version %s by Maurits Evers (maurits.evers\@anu.edu.au)\n", $version);
    printf("Summarise picard-tools MarkDuplicates files.\n");
    printf("Usage:\n");
    printf("  summaris_DuplicationMetrics.pl <dir> <filter>\n");
    printf("\n");
    printf("  <dir>     Directory that contains picard-tools files.\n");
    printf("  <filter>  Regular expression to filter files.\n");
   exit;
}


# From the Perl Cookbook (Recipe 2.17)
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
