#!/usr/bin/perl

# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Last changed: 26-07-2017

use warnings;
use strict;
use IO::Zlib;
use File::Basename;

my $version = 0.9;

usage() if (scalar(@ARGV) != 1);

my $file = $ARGV[0];
die(sprintf("[ERROR] File %s does not exist.\n", $file)) if (! -e $file);

my $fh;
if ($file =~ /gz$/) {
    $fh = new IO::Zlib;
    $fh -> open($file, "rb");
} else {
    open($fh, $file);
}

my $ln = 0;
my $nreads = 0;
my %stats = ();
while (my $line = <$fh>) {
    if ($ln % 4 == 0) {
        if ($line =~ /^@(\w+):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):([YN]):(\d+):(\w+)/) {
            my $instrument  = $1;
            my $run_number  = $2;
            my $flowcell_id = $3;
            my $lane        = $4;
            my $tile        = $5;
            my $x_pos       = $6;
            my $y_pos       = $7;
            my $read        = $8;
            my $is_filtered = $9;
            my $control_no  = $10;
            my $sample_id   = $11;
            # Store stats
            my $arr = 0;
            my $min_x_pos = undef;
            my $max_x_pos = undef;
            my $min_y_pos = undef;
            my $max_y_pos = undef;
            if (exists($stats{$sample_id})) {
                $nreads++;
                $arr = $stats{$sample_id};
                $min_x_pos = $arr->[5];
                $max_x_pos = $arr->[6];
                $min_y_pos = $arr->[7];
                $max_y_pos = $arr->[8];
                $min_x_pos = $x_pos if ($x_pos < $min_x_pos);
                $max_x_pos = $x_pos if ($x_pos > $max_x_pos);
                $min_y_pos = $y_pos if ($y_pos < $min_y_pos);
                $max_y_pos = $y_pos if ($y_pos > $max_y_pos);
            } else {
                $nreads = 1;
                $min_x_pos = $x_pos;
                $max_x_pos = $x_pos;
                $min_y_pos = $y_pos;
                $max_y_pos = $y_pos;
            }
            $arr = [
                $nreads,
                $instrument,
                $run_number,
                $flowcell_id,
                $lane,
                $tile,
                $min_x_pos,
                $max_x_pos,
                $min_y_pos,
                $max_y_pos,
                $read,
                $is_filtered,
                $control_no];
            $stats{$sample_id} = $arr;
        }
    }
    $ln++;
#    last if $ln > 1000;
}

if ($file =~ /gz$/) {
    $fh -> close();
} else {
    close($fh);
}


my $id = [
    "Number of reads",
    "Experiment",
    "Run number",
    "Flowcell ID",
    "Lane",
    "Tile",
    "min_x_pos",
    "max_x_pos",
    "min_y_pos",
    "max_y_pos",
    "Read",
    "Is filtered",
    "Control number"];
printf("Summary for file %s\n", basename($file));
foreach my $key (sort keys %stats) {
    my $arr = $stats{$key};
    printf("Sample ID = %s\n", $key);
    for (my $i = 0; $i < scalar(@{$arr}); $i++) {
        printf("%20s    = %s\n", $id->[$i], $arr->[$i]);
    }
}



sub usage {
    printf("parse_fastq_header version %s by Maurits Evers (maurits.evers\@anu.edu.au)\n", $version);
    printf("Extract run information from FASTQ header.\n");
    printf("Usage:\n");
    printf("  summarise_bowtie2Alignment.pl <dir> <filter>\n");
    printf("\n");
    printf("  <file>    FASTQ file.\n");
    exit;
}
