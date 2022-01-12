#!/aplic/perl/bin/perl

use strict;
use warnings;
use Cwd;
my $inFile=shift or die "plink infile (bed format) ???\n";
my $KeepFile=shift or die "file for --keep option in plink??\n";
my $thresh=shift or die "kinship threshold???\n";
my $maf=shift or die "maf filter applied???\n";
my $outPref=shift or die "outFile Prefixe (with path if needed)\n";
print "remember:
an estimated kinship coefficient range:
>0.354 <=> duplicate/MZ twin
[0.177, 0.354] <=> 1st-degree
[0.0884, 0.177] <=> 2nd-degree
[0.0442, 0.0884] <=> 3rd degree\n";


if( -e $outPref.".RemoveKin".$thresh){
	print $outPref.".RemoveKin".$thresh." already exists\n";
}else{
	####running King within Pop
	if( -e $outPref.".Relationship.kin0"){
		print "kin0 already exists!\n";
	}else{
		print "let's generate kin0 file\n";
		if( !-e $outPref.".bed" && $outPref.".fam"){
			die $outPref.".bed does not exists\n";
		}
	        system("/Users/pierrespc/Documents/PostDoc/scripts/Tools/King/Mac-king-single-thread -b ".$outPref.".bed --kinship --prefix ".$outPref.".Relationship");
#		system("rm ".$outPref.".bed");
#                system("rm ".$outPref.".log");
#                system("rm ".$outPref.".nosex");
#                system("rm ".$outPref.".bim");
	}
	if(! -e $outPref.".Relationship.ListPerInd.kin".$thresh){
	       print "make listInd\n";
		my %listInds;
		foreach my $suf ("kin0","kin"){
			my $fid1;
			my $fid2;
			my $iid1;
			my $iid2;
			my $score;
			if($suf eq "kin"){
				$fid1=0;
				$fid2=0;
				$iid1=1;
				$iid2=2;
				$score=8;
			}else{
				$fid1=0;
                                $fid2=2;
                                $iid1=1;
                                $iid2=3;
                                $score=7;
			}
			open(KIN,$outPref.".Relationship.".$suf) || die ("can t read ".$outPref.".Relationship.".$suf."\n");
			my @table=<KIN>;
			close(KIN);
			shift @table;
			foreach my $line (@table){
				chomp $line;
				my @splitted=split(/\s+/,$line);
				
				print $splitted[$score]."\n";
				if( $splitted[$score] > $thresh){				
					print $line."\n";
					my $ind1=$splitted[$iid1];
					my $ind2=$splitted[$iid2];
					print $splitted[$score]." ".$ind1." ".$ind2."\n";
					if( ! exists $listInds{$ind1}){
						@{$listInds{$ind1}}=(1,$ind2);
					}else{
						${$listInds{$ind1}}[0] ++;
						push(@{$listInds{$ind1}},$ind2);	
					}
					if( ! exists $listInds{$ind2}){
						@{$listInds{$ind2}}=(1,$ind1);
					}else{
						${$listInds{$ind2}}[0] ++;
						push(@{$listInds{$ind2}},$ind1);	
					}
				}
		 	}
		}
	
                open(LIST,">".$outPref.".Relationship.ListPerInd.kin".$thresh) || die ("can't write ".$outPref.".Relationship.ListPerInd.kin".$thresh."\n");

		foreach my $key (keys %listInds){
			print LIST $key."\t".join("\t",@{$listInds{$key}})."\n";
		}
		close(LIST);
	}
	system("Rscript /Users/pierrespc/Documents/PostDoc/scripts/Tools/King//getListIndToRm.R  ".$outPref.".Relationship.ListPerInd.kin".$thresh." ".$outPref.".fam ".$outPref.".RemoveKin".$thresh);
}

