#!/usr/bin/perl

# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Last changed: 15-08-2017

use warnings;
use strict;

my $version = 0.95;

usage() if (scalar(@ARGV) != 2);

my $dir = $ARGV[0];
my $filter = $ARGV[1];
opendir(my $dh, $dir) or doe $!;

my %ht = ();
my $isPE = 0;
while (my $file = readdir($dh)) {
    if ($file =~ /\Q$filter\E/) {
	my $id = $file;
	$file = join("/", $dir, $file);
	open(my $fh, $file) or die $!;
	my $nTotal = 0;
	my $nMapped = 0;
	my $nStep = 0;
	while (my $line = <$fh>) {
        $isPE = 1 if ($line =~ /were paired/);
	    if ($line =~ /(\d+) reads; of these:/) {
		    $nTotal = $1;
		    $nStep++;
	    }
        if ($line =~ /(\d+)\s.*aligned exactly 1 time/) {
            $nMapped += $1;
            $nStep++;
        }
        if ($line =~ /(\d+)\s.*aligned >1 time/) {
            $nMapped += $1;
            $nStep++;
        }
	    if ($line =~ /(\d+)\s.*aligned concordantly exactly 1 time/) {
		    $nMapped += $1;
		    $nStep++;
	    }
	    if ($line =~ /(\d+)\s.*aligned concordantly >1 time/) {
		    $nMapped += $1;
		    $nStep++;
	    }
	    my @arr = ($nTotal, $nMapped);
	    $ht{$id} = [@arr];
	    last if ($nStep == 3);
	}
	close($fh);
    }
}

my $sep = "\t";
my @cn = (
    "Sample",
    "Number of reads",
    "Number of mapped reads",
    "Percentage mapped");
printf("%s\n", join($sep, @cn));
foreach my $key (sort {$a cmp $b} keys %ht) {
    my @arr = @{$ht{$key}};
    if ($isPE == 1) {
        $arr[0] = $arr[0] * 2;
        $arr[1] = $arr[1] * 2;
    }
    printf("%s\n", join(
        $sep,
        $key,
        commify($arr[0]),
        commify($arr[1]),
        sprintf("%3.2f", $arr[1] / $arr[0] * 100.0)));
}


sub usage {
    printf("summarise_bowtie2Alignment version %s by Maurits Evers (maurits.evers\@anu.edu.au)\n", $version);
    printf("Summarise bowtie2 alignment stats based on log files.\n");
    printf("Usage:\n");
    printf("  summarise_bowtie2Alignment.pl <dir> <filter>\n");
    printf("\n");
    printf("  <dir>     Directory that contains bowtie2 summary stats files.\n");
    printf("  <filter>  Regular expression to filter files.\n");
    exit;
}


# From the Perl Cookbook (Recipe 2.17)
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
