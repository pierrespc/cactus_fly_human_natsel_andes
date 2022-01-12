#!/bin/perl

use strict;
use warnings;


my $treeXML=shift or die "treeXML??";
my $fam=shift or die "fam??";


open(FAM,$fam)|| die("can't open $fam\n");

my %hash;
foreach my $line(<FAM>){
	chomp $line;
	my @split=split(/\s+/,$line);
	$hash{$split[1]}=$split[0];
}
close(FAM);

#foreach my $key (keys %hash){
#	print $key." ".$hash{$key}."\n";
#}


open(XML,$treeXML) || die("can't open $treeXML\n");
open(OUT,">$treeXML.Parsed.xml") || die("can't write $treeXML.Parsed.xml\n");

foreach my $line(<XML>){
        chomp $line;
	if($line =~ /<Tree>/){
		print $line."\n";
		foreach my $key (keys %hash){
			$line=~s/\($key:/\($key-$hash{$key}:/g;
			$line=~s/,$key:/,$key-$hash{$key}:/g;
		}
		print $line."\n";
	}
	print OUT $line."\n";
}
close(XML);
close(OUT);
