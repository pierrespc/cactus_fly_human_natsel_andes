#!/bin/perl

use strict;
use warnings;

open(SAMP,"../../Gnecchio/Reformate/Outputs/IHS/Andean.sample") || die("can't open ../../Gnecchio/Reformate/Outputs/IHS/Andean.sample\n");

my @samp=<SAMP>;
close(SAMP);

open(RM,"Outputs/RemoveUros.RM")|| die("can't open Outputs/RemoveUros.RM\n");

my %hash;
foreach my $line (<RM>){
	chomp $line;
	my @split=split(/\s+/,$line);
	$hash{$split[1]}="";
}
close(RM);


system("mkdir Outputs/IHS");
open(OUT,">Outputs/IHS/Andean.sample") || die("can't write Outputs/IHS/Andean.sample\n");
print OUT @samp[(0,1)];

@samp=@samp[(2..$#samp)];
my $count=5;
my @listPOS=(0..4);
foreach my $line (@samp){
	chomp $line;
	my @split=split(/\s+/,$line);
	if(! exists $hash{$split[1]}){
		print OUT $line."\n";
		push(@listPOS,($count,$count+1));
	}
	$count=$count +2;
}
close(OUT);
print "@listPOS\n";

open(OUT,">Outputs/IHS/Andean.haps")|| die("can t write Outputs/IHS/Andean.haps\n");
open(HAP,"../../Gnecchio/Reformate/Outputs/IHS/Andean.haps") || die("can,t open ../../Gnecchio/Reformate/Outputs/IHS/Andean.haps\n");
foreach my $line (<HAP>){
	chomp $line;
	my @split=split(/\s+/,$line);
	print OUT "@split[@listPOS]\n";
}
close(OUT);



