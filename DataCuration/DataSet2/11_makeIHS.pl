#!/bin/perl

use strict;
use warnings;

open(FAM,"Outputs/6pops.KinshipFilter.AdmFilter.CPFS.fam") || die("can't open Outputs/6pops.KinshipFilter.AdmFilter.CPFS.fam\n");

my %hash;
foreach my $line (<FAM>){
	chomp $line;
	my @split=split(/\s+/,$line);
	if($split[0] eq "Andean"){
		if(exists $hash{$split[1]}){
			die($split[1]." duplicated in fam\n");
		}
		$hash{$split[1]}="";
	}
}
close(FAM);

open(SAMP,"Outputs/fineStructureOnlyNat/shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr22_alignedRef_phased.sample") || die("can't open Outputs/fineStructureOnlyNat/shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr22_alignedRef_phased.sample \n");

my @samp=<SAMP>;
close(SAMP);

system("mkdir Outputs/IHS");
open(OUT,">Outputs/IHS/Andean.sample") || die("can't write Outputs/IHS/Andean.sample\n");
print OUT @samp[(0,1)];

@samp=@samp[(2..$#samp)];
my $count=5;
my @listPOS=(0..4);
foreach my $line (@samp){
	chomp $line;
	my @split=split(/\s+/,$line);
	if(exists $hash{$split[1]}){
		if($hash{$split[1]} ne ""){
			 die($split[1]." duplicated in sample\n");
		}
		$split[0]="Andean";
		print OUT join(" ",@split)."\n";
		$hash{$split[1]}=$count;
		push(@listPOS,($count,$count+1));
	}
	$count=$count +2;
}
close(OUT);
print "@listPOS\n";

open(OUT,">Outputs/IHS/Andean.haps")|| die("can t write Outputs/IHS/Andean.haps\n");
for(my $chr=1;$chr<23; $chr ++){
	print $chr."\n";
	open(HAP,"Outputs/fineStructureOnlyNat/shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr".$chr."_alignedRef_phased.haps") || die("can't open Outputs/fineStructureOnlyNat/shapeITphased/6pops.KinshipFilter.AdmFilter.MAF0.0000001.GENO0.02.MIND0.05.chr".$chr."_alignedRef_phased.haps\n"); 
	foreach my $line (<HAP>){
		chomp $line;
		my @split=split(/\s+/,$line);
		print OUT "@split[@listPOS]\n";
	}
	close(HAP);
}
close(OUT);



