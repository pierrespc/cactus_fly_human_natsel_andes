#!/bin/perl


use strict;
use warnings;
open(MAP,"Outputs/6pops.KinshipFilter.AdmFilter.CPFS.bim") || die();
my $chr=0;


sub flip {
	if(! $_[0]){
		return("");
	}elsif($_[0] eq "A"){
		return("T");
	}elsif($_[0] eq "T"){
                return("A");
        }elsif($_[0] eq "G"){
                return("C");
        }elsif($_[0] eq "C"){
                return("G");
	}elsif($_[0] eq "-"){
                return("-");
        }elsif($_[0] eq "0" ){
                return("0");
        }else{
		die("unknown alleles ".$_[0]."\n");
	}
}


my %hash;
foreach my $line (<MAP>){
	chomp $line;
	my @split=split(/\s+/,$line);
	if(exists $hash{$split[0].".".$split[3]}){
		die("duplicated ".$split[0].".".$split[3]." in bim\n");
	}
	$hash{$split[0].".".$split[3]}=[@split[(4,5)]];
}

#foreach my $key (keys %hash){
#	print "$key -> @{$hash{$key}}\n";
#}
#die();
open(OUT,">Outputs/6pops.KinshipFilter.AdmFilter.CPFS.ancestral")||die("write\n"); 
open(WARN,">Outputs/6pops.KinshipFilter.AdmFilter.CPFS.ancestral.WARN")||die("write\n"); 

for(my $chr=1;$chr<23;$chr ++){
	print $chr."\n";
	open(ANC,"/Volumes/MARIOLOLO/Data/ancestral_hg19_1KG/1KG.hg19.chr$chr.ancestral") || die("anc\n");
	foreach my $line (<ANC>){
	        chomp $line;
        	my @split=split(/\s+/,$line);
		#$split[1]=$split[1]+1;
		if(! exists $hash{$split[0].".".$split[1]}){
			next;
		}
		my @der=split(/,/,$split[3]);
		my $anc=$split[2];
		my $minor=${$hash{$split[0].".".$split[1]}}[0];
        	my $major=${$hash{$split[0].".".$split[1]}}[1];

		if(grep(/^$minor/,(@der,$anc,"0")) ||  grep(/^$major/,(@der,$anc))){
			if(@der == 1){
				print OUT $line."\n";
			}else{
				my $ff=$der[0];
				foreach my $a (@der){
					if($a eq $minor || $a eq $major){
						$ff=$a;
						last;
					}
				}
				print OUT join("\t",(@split[(0,1)],$anc,$ff))."\n";
			}
		}else{
			#print "$anc @der\n";
			#print $line."\n";
			$anc=flip($anc);
			for(my $ff=0;$ff<@der;$ff ++){
				$der[$ff]=flip($der[$ff]);
			}
                        print "$anc @der\n";

			if(grep(/^$minor/,(@der,$anc,"0")) ||  grep(/^$major/,(@der,$anc))){
                 	       if(@der == 1){
                        	        print OUT $line."\n";
                        	}else{
                                	my $ff=$der[0];
                                	foreach my $a (@der){
                                        	if($a eq $minor || $a eq $major){
                                                	$ff=$a;
                                                	last;
                                        	}
                               	 	}
                                	print OUT join("\t",(@split[(0,1)],$anc,$ff))."\n";
                        	}
                	}else{				
				print WARN $line." vs $minor $major\n";
			}
		}
		close(ANC);
	}
}
close(OUT);
close(WARN);	 
	

