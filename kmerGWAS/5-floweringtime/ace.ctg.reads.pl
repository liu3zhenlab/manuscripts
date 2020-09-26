#!/usr/bin/perl -w
use strict;
use warnings;

# CO Contig1 39 12 10 U
# AAAAACGCGAGATTGAAGATCTTCAGAAACAGCTGGCTA
#
# BQ
# 10 10 10 25 30 35 40 45 50 55 65 65 65 65 70 70 70 70 70 70 70 70 70 70 70 65 60 55 50 45 40 35 30 10 10 10 10 10 10 
#
# AF k105 U 1

my ($ctg, $ctg_len);

open(IN, $ARGV[0]) || die;
while(<IN>) {
	if (/^CO (Contig\d+) (\d+) /) {
		$ctg = $1;
		$ctg_len = $2;
	}

	if (/^AF (\S+)/) {
		print "$ctg\t$ctg_len\t$1\n";
	}

}
