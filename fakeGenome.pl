#!/usr/bin/perl
#============================================================================
# Name        		: fakeGenome.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Tue Dec 13 12:44:46 2011
# Last Modified By	: 
# Last Modified On	: Tue Dec 13 12:44:46 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use FindBin qw($Bin);
use lib "$Bin";
use Config::File;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Data::Dumper;
use File::Basename;
use File::Path;
use Math::Random qw(random_poisson random_normal random_exponential);
use common;

my $version = '1.0.0';
my ($in,$out);
if(@ARGV<1) { usage(); }

## for test 
my $srandSeed = int(rand (2**31));
srand ($srandSeed);
printf "srand seed:\t%s\n",$srandSeed;
## 

my $command = shift(@ARGV);
my %func = (
			germlineSite=>\&germlineSite,
			somaticSite=>\&somaticSite,
			somaticCNV=>\&somaticCNV,
			genFasta=>\&genFasta,
			);
die("Unknown command \"$command\".\n") if (!defined($func{$command}));

&{$func{$command}};

############################################################################
sub genFasta
{
	my ($opt_h,$opt_a,$outVCF,$inSite,$outSite,$outFasta,$out_coordinateMapping,$g_sampleID);
	GetOptions("h"	=>\$opt_h,"s=s"=>\$g_sampleID,"a"=>\$opt_a);
	if(@ARGV<6 || $opt_h)
	{
		die(qq/
Usage:   fakeGenome.pl genFasta  [options] <germline dir> <somatic dir> <cnv dir> <fa dir> <region dir> <out dir> 

Options:
		-s		sampleID
		-a		fasta with gap ['-']
		-h		print this help info
		\n/);
	}
	my $germlineDir=shift @ARGV;
	my $somaticDir=shift @ARGV;
	my $cnvDir=shift @ARGV;
	my $faDir=shift @ARGV;
	my $regionDir=shift @ARGV;
	my $outDir=shift @ARGV;
	mkpath $outDir;
	if(!defined($g_sampleID)) { $g_sampleID="TSomaticGenome01"; }
	printf "======================start:fakeGenome.pl genFasta(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
	my $finalVCFFile=$outDir."/$g_sampleID.all.vcf";
	open $outVCF,">",$finalVCFFile or die "Cann't open file $finalVCFFile ($!)\n";
	outputVCFHeader($outVCF,$g_sampleID);
	close $outVCF;
	open $outVCF,"|sort -k 1,1 -k 2g,2 >>$finalVCFFile" or die "($!)\n";
	for my $_chr (1..22,"X","Y")
	{
		my $g_chr="chr$_chr";
		my $regionFile = $regionDir."/$g_chr.bed";
		if(! -f $regionFile)
		{
			printf STDERR "skip $g_chr, no file: $regionFile\n";
			next;
		}

		my $outCoordinateMappingFile = $outDir."/$g_chr.coord.txt";
		open $out_coordinateMapping,">",$outCoordinateMappingFile or  die "Cann't open file $outCoordinateMappingFile($!)\n"; 
		my %newCoordinateMapTable_info=();
		my %newCoordinateMapTable_search=();

		my $outFastaFile = $outDir."/$g_chr.fa";
		open $outFasta,">",$outFastaFile or die "Cann't open file $outFastaFile($!)\n";

		my $inFastaFile = $faDir."/$g_chr.fa";
		my $pSeq=readSeq($g_chr,$inFastaFile);
		my $seqHap1=$pSeq->{'seq'};
		my $seqHap2=$pSeq->{'seq'};
		my $g_len=length($pSeq->{'seq'});
		undef $pSeq;

		my $germlineVarFile=$germlineDir."/germline.$g_chr.vcf";
		my $somaticVarFile=$somaticDir."/somatic.$g_chr.vcf";
		my $CNVFile=$cnvDir."/somaticCNV.$g_chr.vcf";
		my %VarTable=();
		my %CNVTable=();
		readVar($germlineVarFile,\%VarTable,\%CNVTable);
		readVar($somaticVarFile,\%VarTable,\%CNVTable);
		readVar($CNVFile,\%VarTable,\%CNVTable);

		##ensure no overlap among indels, SNPs
		my @_varSite=sort { $a<=>$b } keys %VarTable;
		for(my $i=0;$i<@_varSite;$i++)
		{
			if($VarTable{$_varSite[$i]}->{'type'} eq "Del")
			{
				my $_maxDelLen;
				if($i<$#_varSite)
				{
					#no cover netx var site
					$_maxDelLen=$_varSite[$i+1]-$_varSite[$i]-1;
					if(!($VarTable{$_varSite[$i]}->{'l'}<=$_maxDelLen))
					{
						printf STDERR "Deletion occur in $g_chr:%s (len:%i) is in conflict with next var site\n",$_varSite[$i],$VarTable{$_varSite[$i]}->{'l'};
						if($_maxDelLen==0) { delete $VarTable{$_varSite[$i]}; }
						else { $VarTable{$_varSite[$i]}->{'l'}=$_maxDelLen; }
					}
				}elsif($i==$#_varSite)
				{
					$_maxDelLen=$g_len-$_varSite[$i];
					if(!($VarTable{$_varSite[$i]}->{'l'}<=$_maxDelLen))
					{
						printf STDERR "Deletion occur in $g_chr:%s (len:%i) is in conflict with genome length($g_len)\n",$_varSite[$i],$VarTable{$_varSite[$i]}->{'l'};
						if($_maxDelLen==0) { delete $VarTable{$_varSite[$i]}; }
						else { $VarTable{$_varSite[$i]}->{'l'}=$_maxDelLen; }
					}
				}
			}
		}
		###
		# effect of CNV
		###
		foreach my $_breakpoint (sort { $a<=>$b } keys %CNVTable)
		{
			my $_CN=$CNVTable{$_breakpoint}->{'CN'};
			my $_CNVLen=$CNVTable{$_breakpoint}->{'l'};
			my $_CNVHaplotype=$CNVTable{$_breakpoint}->{'haplotype'};
			#get the variant in the CNV region
			my @_varInCNVR=();
			#need to optimize
			foreach my $_varSite (sort { $a<=>$b } keys %VarTable)
			{
				if(0<=$_varSite-$_breakpoint && $_varSite-$_breakpoint<=$_CNVLen) { push @_varInCNVR,$_varSite; }
			}
			## different Copy Number
			if($_CN==0)
			{
				#remove all SNP&SNV and INDEL in the CNV region
				foreach my $_varSite (@_varInCNVR)
				{
					my $_VC=$VarTable{$_varSite}->{'VC'};
					printf STDERR "remove %s($_VC) in CNV region(breakpoint:$g_chr:%s,len:$_CNVLen,CN:$_CN)\t$g_chr:%s\n",$VarTable{$_varSite}->{'type'},$_breakpoint,$_varSite;
					delete $VarTable{$_varSite};
				}
				#change the sequences;
				substr($seqHap1,$_breakpoint,$_CNVLen,"-"x$_CNVLen);
				substr($seqHap2,$_breakpoint,$_CNVLen,"-"x$_CNVLen);
			}elsif($_CN==1)
			{
				#remove all INDEL in the CNV region, remain one allele in the SNP&SNV site
				foreach my $_varSite (@_varInCNVR)
				{
					my $_VC=$VarTable{$_varSite}->{'VC'};
					if($VarTable{$_varSite}->{'type'}=~/(Ins|Del)/)
					{
						printf STDERR "remove %s($_VC) in CNV region(breakpoint:$g_chr:%s,len:$_CNVLen,CN:$_CN)\t$g_chr:%s\n",$VarTable{$_varSite}->{'type'},$_breakpoint,$_varSite;
						delete $VarTable{$_varSite};
					}else
					{
						if($VarTable{$_varSite}->{'z'} eq "1/1") 
						{
							$VarTable{$_varSite}->{'CN'}=$_CN;
							$VarTable{$_varSite}->{'CNVHaplotype'}=$_CNVHaplotype;
							next; 
						}
						if($_CNVHaplotype==$VarTable{$_varSite}->{'haplotype'})
						{
							#remain ref allele, so remove SNP & SNV
							printf STDERR "remove %s($_VC) in CNV region(breakpoint:$g_chr:%s,len:$_CNVLen,CN:$_CN)\t$g_chr:%s\n",$VarTable{$_varSite}->{'type'},$_breakpoint,$_varSite;
							delete $VarTable{$_varSite};
						}else
						{	
							#loss the ref allele
							$VarTable{$_varSite}->{'z'}="1/1";
							$VarTable{$_varSite}->{'CN'}=$_CN;
							$VarTable{$_varSite}->{'CNVHaplotype'}=$_CNVHaplotype;
							printf STDERR "loss ref allele ($_VC %s) in CNV region(breakpoint:$g_chr:%s,len:$_CNVLen,CN:$_CN)\t$g_chr:%s\n",$VarTable{$_varSite}->{'type'},$_breakpoint,$_varSite;
						}
					}
				}
				#change the sequences
				if($_CNVHaplotype==1)
				{
					substr($seqHap1,$_breakpoint,$_CNVLen,"-"x$_CNVLen);
				}elsif($_CNVHaplotype==2)
				{
					substr($seqHap2,$_breakpoint,$_CNVLen,"-"x$_CNVLen);
				}
			}else
			{
				#no effect when CN>2, only produce "new chromosome"
				foreach my $_varSite (@_varInCNVR)
				{
					$VarTable{$_varSite}->{'CN'}=$_CN;
					$VarTable{$_varSite}->{'CNVHaplotype'}=$_CNVHaplotype;
				}
			}
		}
		###
		# chang the sequence
		###
		foreach my $_refPos (sort { $a<=>$b } keys %VarTable)
		{
			#SNP & SNV
			if($VarTable{$_refPos}->{'type'} ne "SNP") { next; }
			my $_refBase=$VarTable{$_refPos}->{'ref'};
			my $_altBase=$VarTable{$_refPos}->{'alt'};
			my $_z=$VarTable{$_refPos}->{'z'};
			my $_haplotype=$VarTable{$_refPos}->{'haplotype'};
			if($_z eq "1/1")
			{
				if($VarTable{$_refPos}->{'CN'} >= 2)
				{
					substr($seqHap1,$_refPos-1,1,$_altBase);
					substr($seqHap2,$_refPos-1,1,$_altBase);
				}else
				{
					if(substr($seqHap1,$_refPos-1,1) eq "-") { substr($seqHap2,$_refPos-1,1,$_altBase); }
					else{ substr($seqHap1,$_refPos-1,1,$_altBase); }
				}
			}else
			{
				if($_haplotype==1)
				{
					substr($seqHap1,$_refPos-1,1,$_altBase);
				}else
				{
					substr($seqHap2,$_refPos-1,1,$_altBase);
				}
			}
		}
		my $vOffset=0;
		my ($_mapBeg,$_mapEnd,$_mapOp,$_mapBegV,$_mapEndV,$_mapZ)=(1,0,".",1,0,".");
		foreach my $_refPos (sort { $a<=>$b } keys %VarTable)
		{
			#INDEL
			if($VarTable{$_refPos}->{'type'}!~/(Ins|Del)/) { next; }
			my $_vPos=$_refPos+$vOffset;
			my $_refBase=$VarTable{$_refPos}->{'ref'};
			my $_altBase=$VarTable{$_refPos}->{'alt'};
			my $_z=$VarTable{$_refPos}->{'z'};
			my $_haplotype=$VarTable{$_refPos}->{'haplotype'};
			my $_op=$VarTable{$_refPos}->{'op'};
			my $_l=$VarTable{$_refPos}->{'l'};
			if($_op eq "+")
			{
				if($_z eq '0/1')
				{
					#het ins
					if($_haplotype==1) 
					{ 
						substr($seqHap1,$_vPos,0,$_altBase); 
						substr($seqHap2,$_vPos,0,"-"x$_l); 
					}else 
					{ 
						substr($seqHap1,$_vPos,0,"-"x$_l); 
						substr($seqHap2,$_vPos,0,$_altBase); 
					}
				}else
				{
					#hom ins, nothing to do
					substr($seqHap1,$_vPos,0,$_altBase);
					substr($seqHap2,$_vPos,0,$_altBase);
				}
				$vOffset+=$_l;
			}else
			{
				if($_z eq '1/1')
				{
					#hom del,
					substr($seqHap1,$_vPos,$_l,"");
					substr($seqHap2,$_vPos,$_l,"");
					$vOffset-=$_l;
				}else
				{
					#het del,modify the portation
					if($_haplotype==1) { substr($seqHap1,$_vPos,$_l,"-"x$_l); }
					else { substr($seqHap2,$_vPos,$_l,"-"x$_l); }
				}
			}
			printf $out_coordinateMapping "$g_chr\t$_mapBeg\t$_refPos\t$_mapOp\t$_mapZ\t%i\t%i\n",$_mapBegV,$_vPos;
			if(!defined($newCoordinateMapTable_search{$g_chr})) { $newCoordinateMapTable_search{$g_chr}=[]; }
			push @{$newCoordinateMapTable_search{$g_chr}},$_mapBeg,$_refPos;
			if(!defined($newCoordinateMapTable_info{$g_chr})) { $newCoordinateMapTable_info{$g_chr}={}; }
			$newCoordinateMapTable_info{$g_chr}->{"$_mapBeg-$_refPos"}={'nBeg'=>$_mapBegV,'nEnd'=>$_vPos,'op'=>$_mapOp,'z'=>$_mapZ};
			$_mapBeg=$_refPos+1;
			$_mapBegV=$_vPos+1;
			$_mapZ=$_z;
			$_mapOp=sprintf "$_op$_l";
		}
		printf $out_coordinateMapping "$g_chr\t$_mapBeg\t$g_len\t$_mapOp\t$_mapZ\t%i\t%i\n",$_mapBegV,$g_len+$vOffset;
		if(!defined($newCoordinateMapTable_search{$g_chr})) { $newCoordinateMapTable_search{$g_chr}=[]; }
		push @{$newCoordinateMapTable_search{$g_chr}},$_mapBeg,$g_len;
		if(!defined($newCoordinateMapTable_info{$g_chr})) { $newCoordinateMapTable_info{$g_chr}={}; }
		$newCoordinateMapTable_info{$g_chr}->{"$_mapBeg-$g_len"}={'nBeg'=>$_mapBegV,'nEnd'=>$g_len+$vOffset,'op'=>$_mapOp,'z'=>$_mapZ};
		#output "exome" fasta
		##$seqHap1=~s/N/G/g;
		##$seqHap2=~s/N/G/g;
		open $in,$regionFile or die "Cann't open file $regionFile ($!)\n";
		while(<$in>)
		{
			chomp;
			my @field=split /\t/;
			my ($chr,$beg,$end)=@field[0,1,2];
			$beg++;
			my $_newBeg=newCoordinate(\%newCoordinateMapTable_info,$g_chr,search(\%newCoordinateMapTable_search,$g_chr,$beg),$beg);
			my $_newEnd=newCoordinate(\%newCoordinateMapTable_info,$g_chr,search(\%newCoordinateMapTable_search,$g_chr,$end),$end);
			
			my $newSeq = substr($seqHap1,$_newBeg-1,$_newEnd-$_newBeg+1);
			my $_SEQID=sprintf "SEQ_HAP1\@($g_chr:%s-%s)#",$beg,$end;
			if(!$opt_a) { $newSeq=~s/-//g; }
			write_fasta_seq(\$newSeq,$_SEQID,$outFasta);

			$newSeq = substr($seqHap2,$_newBeg-1,$_newEnd-$_newBeg+1);
			$_SEQID=sprintf "SEQ_HAP2\@($g_chr:%s-%s)#",$beg,$end;
			if(!$opt_a) { $newSeq=~s/-//g; }
			write_fasta_seq(\$newSeq,$_SEQID,$outFasta);

		}
		#produce new chromosomes caused by CNV>2
		foreach my $_breakpoint (sort { $a<=>$b } keys %CNVTable)
		{	
			my $_CN=$CNVTable{$_breakpoint}->{'CN'};
			if($_CN<=2) { next; }
			my $_CNVLen=$CNVTable{$_breakpoint}->{'l'};
			my $_CNVHaplotype=$CNVTable{$_breakpoint}->{'haplotype'};
			#trasform coordiante
			my $_endpoint=$_breakpoint+$_CNVLen;
			$_breakpoint++;
			my $_oBreakpoint=$_breakpoint;
			$_breakpoint=newCoordinate(\%newCoordinateMapTable_info,$g_chr,search(\%newCoordinateMapTable_search,$g_chr,$_breakpoint),$_breakpoint);
			$_endpoint=newCoordinate(\%newCoordinateMapTable_info,$g_chr,search(\%newCoordinateMapTable_search,$g_chr,$_endpoint),$_endpoint);
			## for "capture"
			#$_breakpoint-=250;
			#$_endpoint+=250;
			##
			#get sequence
			my $_isGapHaplotype;
			my @_newSeq=();
			my $newCNVSeq="";
			if($_CNVHaplotype==1)
			{
				$newCNVSeq = substr($seqHap1,$_breakpoint-1,$_endpoint-$_breakpoint+1);
			}else
			{
				$newCNVSeq = substr($seqHap2,$_breakpoint-1,$_endpoint-$_breakpoint+1);
			}
			for my $i (3..$_CN)
			{
				my $_CNVID=sprintf "CNV$i\@($g_chr:%s-%s)#",$_oBreakpoint,$_oBreakpoint+$_CNVLen-1;
				if(!$opt_a) { $newCNVSeq=~s/-//g; }
				write_fasta_seq(\$newCNVSeq,$_CNVID,$outFasta);
			}
		}
		## output final vcf
		foreach my $_pos (sort { $a<=>$b } keys %VarTable)
		{
			my ($ref,$alt,$info);
			if($VarTable{$_pos}->{'type'} eq "Ins") { $ref=$VarTable{$_pos}->{'ref'}; $alt=$VarTable{$_pos}->{'ref'}.$VarTable{$_pos}->{'alt'}; }
			elsif($VarTable{$_pos}->{'type'} eq "Del") { $ref=$VarTable{$_pos}->{'ref'}.$VarTable{$_pos}->{'alt'}; $alt=$VarTable{$_pos}->{'ref'}; }
			else { $ref=$VarTable{$_pos}->{'ref'}; $alt=$VarTable{$_pos}->{'alt'}; }
			$info=sprintf "VC=%s;HAP=%s;",$VarTable{$_pos}->{'VC'},$VarTable{$_pos}->{'haplotype'};
			if($VarTable{$_pos}->{'CN'} !=2) { $info.=sprintf "CNVHAP=%s;",$VarTable{$_pos}->{'CNVHaplotype'}; }
			printf $outVCF "$g_chr\t$_pos\t%s\t%s\t%s\t.\t.\t$info\tGT:CN\t%s:%s\n",$VarTable{$_pos}->{'rid'},$ref,$alt,$VarTable{$_pos}->{'z'},$VarTable{$_pos}->{'CN'};
		}
		foreach my $_pos (sort { $a<=>$b } keys %CNVTable)
		{
			my ($ref,$alt)=($CNVTable{$_pos}->{'ref'},$CNVTable{$_pos}->{'alt'});
			my $info=sprintf "VC=%s;HAP=%s;SVLEN=%s;END=%s",$CNVTable{$_pos}->{'VC'},$CNVTable{$_pos}->{'haplotype'},$CNVTable{$_pos}->{'l'},$CNVTable{$_pos}->{'end'};
			printf $outVCF "$g_chr\t$_pos\t%s\t%s\t%s\t.\t.\t$info\tGT:CN\t./.:%s\n",$CNVTable{$_pos}->{'rid'},$ref,$alt,$CNVTable{$_pos}->{'CN'};
		}
	}
	printf "========================end:fakeGenome.pl genFasta(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
}

sub somaticCNV
{
	my ($opt_h,$outVCF,$inSite,$outSite,$out,$g_sampleID);
	GetOptions("h"	=>\$opt_h,"s=s"=>\$g_sampleID);
	if(@ARGV<3 || $opt_h)
	{
		die(qq/
Usage:   fakeGenome.pl somaticCNV [options] <config file> <region dir> <out dir>

Options: 
		-s		sampleID
		-h		print this help info
		\n/);
	}
	my $configuration_file=shift @ARGV;
	my $regionDir=shift @ARGV;
	my $outDir=shift @ARGV;
	if(!defined($g_sampleID)) { $g_sampleID="TSomaticGenome01"; }
	
	printf "======================start:fakeGenome.pl somaticCNV(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
	mkpath $outDir;
	my $config_hash = Config::File::read_config_file($configuration_file);
	my $CNV_rate=eval $config_hash->{'CNV_rate'};
	for my $_chr (1..22,"X","Y")
	{
		my $g_chr="chr$_chr";
		my $regionFile=$regionDir."/$g_chr.bed";
		#my $regionFile = $config_hash->{'dir_region'}."/$g_chr.bed";
		if(! -f $regionFile)
		{
			printf STDERR "skip $g_chr, no file: $regionFile\n";
			next;
		}
		my ($pRegion,$totRegion) = readList($regionFile);
		my $nCNV = int($CNV_rate*$totRegion);
		printf "$g_chr\tnCNV: $nCNV\n";
		my $pCNVRegion = selection_sample($pRegion,$nCNV,1);
		open $out,">","$outDir/somaticCNV.$g_chr.vcf" or die "$outDir/somaticCNV.$g_chr.vcf\n";
		outputVCFHeader($out,$g_sampleID);
		
		my $inFastaFile = $config_hash->{'file_ref'}."/$g_chr.fa";
		my $pSeq=readSeq($g_chr,$inFastaFile);
		my $seq=$pSeq->{'seq'};
		my $g_len=length($pSeq->{'seq'});
		undef $pSeq;
		
		foreach (@$pCNVRegion)
		{
			my @field=split /\t/;
			my ($chr,$beg,$end)=@field[0,1,2];	### bed format
			##Copy Number
			my $_CN;
			while(1)
			{
				$_CN=int(random_normal(1,0,2)+2);
				$_CN=$_CN<0?0:$_CN;
				if($_CN!=2) { last; }
			}
			my $_haplotype;
			if($_CN==0) { $_haplotype=3; }
			else
			{
				$_haplotype=rand(1)<0.5?1:2;
			}
			#output CNV 
			my $_refBase=substr($seq,$beg,1);
			printf $out "$chr\t$beg\t.\t$_refBase\t<CNV>\t.\t.\tVC=SomaticCNV;HAP=$_haplotype;END=$end;SVLEN=%i\tGT:CN\t./.:$_CN\n",$end-$beg;
		}
	}
	printf "========================end:fakeGenome.pl somaticCNV(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
}

sub somaticSite
{
	my ($opt_h,$outVCF,$inSite,$outSite,$out,$g_sampleID);
	GetOptions("h"	=>\$opt_h,"s=s"=>\$g_sampleID);
	if(@ARGV<4 || $opt_h)
	{
		die(qq/
Usage:   fakeGenome.pl somaticSite [options] <config file> <germline dir> <region dir> <out dir> 

Options: 
		-s		sampleID [default TSomaticGenome01]
		-h		print this help info
		\n/);
	}
	my $configuration_file=shift @ARGV;
	my $germlineDir=shift @ARGV;
	my $regionDir=shift @ARGV;
	my $outDir=shift @ARGV;
	mkpath $outDir;
	if(!defined($g_sampleID)){ $g_sampleID="TSomaticGenome01"; }

	printf "======================start:fakeGenome.pl somaticSite(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);

	my $config_hash = Config::File::read_config_file($configuration_file);
	my $SNV_rate=eval $config_hash->{'SNV_rate'};
	my $indel_rate=eval $config_hash->{'indel_rate'};
	my $CNV_rate=eval $config_hash->{'CNV_rate'};
	for my $_chr (1..22,"X","Y")
	{
		my $g_chr="chr$_chr";
		my $siteFile = "$germlineDir/germline.$g_chr.remainSite.txt.gz";
		my $regionFile=$regionDir."/$g_chr.bed";
		if(! -f $regionFile)
		{
			printf STDERR "skip $g_chr, no file: $regionFile\n";
			next;
		}
		my @siteArray=();
		if($siteFile=~/\.gz$/) { $inSite=new IO::Uncompress::Gunzip $siteFile or die "gunzip failed: $GunzipError\n"; }
		else { open $inSite, $siteFile or die "Cann't open file $siteFile ($!)\n"; }
		#my $lineNum=0;
		while(<$inSite>)
		{
			chomp;
			#$lineNum++;
			#if(/^$/){ printf "$lineNum\n";next;}
			push @siteArray,$_;
		}
		my ($pNULL,$totSite)=readRegionFile($regionFile,1);
		my $nSNV = int($totSite*$SNV_rate);
		my $nINDEL = int($totSite*$indel_rate);
		printf "$g_chr\ttotal site: $totSite\tnSNV: $nSNV\tnINDEL: $nINDEL\n";
		my $pSNVSite = selection_sample(\@siteArray,$nSNV+$nINDEL,0);
		my $pINDELSite = selection_sample($pSNVSite,$nINDEL,0);
		
		my $inFastaFile = $config_hash->{'file_ref'}."/$g_chr.fa";
		my $pSeq=readSeq($g_chr,$inFastaFile);
		my $seq=$pSeq->{'seq'};
		my $g_len=length($pSeq->{'seq'});
		undef $pSeq;

		my $het_hom_ratio = $config_hash->{'het_hom_ratio'};

		### output remaining site (mutatable sites)
		my $outSite = new IO::Compress::Gzip "$outDir/somatic.$g_chr.remainSite.txt.gz" or die "IO::Compress::Gzip failed: $GzipError\n";
		#open $outSite,"| gzip -c > $outDir/somatic.$g_chr.remainSite.txt.gz" or die "Cann't open file $outDir/somatic.$g_chr.remainSite.txt.gz ($!)\n";	
		foreach (@siteArray)
		{
			print $outSite "$_\n";
		}
		### output variants
		open $outVCF,">","$outDir/somatic.$g_chr.vcf" or die "Cann't open file $outDir/somatic.$g_chr.vcf\n";
		outputVCFHeader($outVCF,$g_sampleID);
		foreach (@$pSNVSite)
		{
			my ($chr,$pos)=split /:/;
			my ($ref,$alt,$_z,$_gt,$_haplotype);
			$ref=substr($seq,$pos-1,1);
			if($ref eq 'N') { next; }
			$alt=$OTHER_ALLELE{$ref}->[rand(scalar @{$OTHER_ALLELE{$ref}})];
			if(rand($het_hom_ratio+1)<1)
			{
				$_z="1/1";
				$_gt=$Abbrev{"$alt$alt"};
				$_haplotype=3;
			}else
			{
				$_z="0/1";
				$_gt=$Abbrev{"$ref$alt"};
				$_haplotype=rand(1)<0.5?1:2;
			}
			printf $outVCF "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tVC=SomaticSNV;HAP=$_haplotype;IUPAC=$_gt\tGT\t$_z\n";
		}
		foreach (@$pINDELSite)
		{
			my ($chr,$pos)=split /:/;
			my ($ref,$alt,$_z,$_gt,$_haplotype,$_l,$_s,$isIns);
			$ref=substr($seq,$pos-1,1);
			if(rand(1)<0.5)
			{
				$_z="0/1";
				$_haplotype=rand(1)<0.5?1:2;
			}else
			{
				$_z="1/1";
				$_haplotype=3;
			}
			$isIns=rand(1)<0.5?1:0;
			$_l=int(random_exponential())+1;
			if($isIns == 1)
			{
				#insertion (after the site specified in $_)
				$_s=randomSequence($_l);
				$_s=join("",@$_s);
				printf $outVCF "$chr\t$pos\t.\t$ref\t$ref$_s\t.\t.\tVC=SomaticINDEL;HAP=$_haplotype\tGT\t$_z\n";
			}else
			{
				#deletion (after the site specified in $_)
				$_s=substr($seq,$pos,$_l);
				printf $outVCF "$chr\t$pos\t.\t$ref$_s\t$ref\t.\t.\tVC=SomaticINDEL;HAP=$_haplotype\tGT\t$_z\n";
			}
		}
	}
	printf "========================end:fakeGenome.pl somaticSite(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
}



sub germlineSite
{
	my ($opt_h,$outVCF,$outSite,$g_sampleID);
	GetOptions("h"	=>\$opt_h,"s=s"=>\$g_sampleID);
	if(@ARGV<3 || $opt_h)
	{
		die(qq/
Usage:   fakeGenome.pl germlineSite [options] <config file> <region dir> <out dir> 

Options: 
		-s		sampleID
		-h		print this help info
		\n/);
	}
	my $configuration_file=shift @ARGV;
	my $regionDir=shift @ARGV;
	my $outDir=shift @ARGV;
	mkpath $outDir;
	if(!defined($g_sampleID)) { $g_sampleID="GermlineGenome01"; }

	printf "======================start:fakeGenome.pl germlineSite(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
	my $config_hash = Config::File::read_config_file($configuration_file);
	my $snp_random_mode=$config_hash->{'snp_random_mode'};
	my $het_hom_ratio=$config_hash->{'het_hom_ratio'};
	my $ti_tv_ratio_known=$config_hash->{'ti_tv_ratio_known'};
	my $ti_tv_ratio_novel=$config_hash->{'ti_tv_ratio_novel'};
	my $ti_tv_ratio=$config_hash->{'ti_tv_ratio'};
	my $indel_random_mode=$config_hash->{'indel_random_mode'};
	my $random_mode=$config_hash->{'random_mode'};

	my $dbSNP_SNP_rate=$config_hash->{'dbSNP_SNP_rate'};
	my $dbSNP_INDEL_rate=$config_hash->{'dbSNP_INDEL_rate'};
	my $snp_rate=eval $config_hash->{'snp_rate'};
	my $indel_rate=eval $config_hash->{'indel_rate'};
	my $CNV_rate=eval $config_hash->{'CNV_rate'};

	for my $_chr (1..22,"X","Y")
	{
		my $g_chr="chr$_chr";
		my $snpFile = $config_hash->{'dir_SNP_dbSNP'}."/$g_chr.vcf.gz";
		my $indelFile = $config_hash->{'dir_INDEL_dbSNP'}."/$g_chr.vcf.gz";
		my $regionFile = $regionDir."/$g_chr.bed";
		if(! -f $regionFile)
		{
			printf STDERR "skip $g_chr, no file: $regionFile\n";
			next;
		}
		my ($pRegionSite,$totSite) = readRegionFile($regionFile);
		my $nSNP=int($totSite*$snp_rate);
		my $nINDEL=int($totSite*$indel_rate);
		my $nCNV=int($totSite*$CNV_rate);
		
		open $outVCF,">","$outDir/germline.$g_chr.vcf" or die "Cann't open file $outDir/germline.$g_chr.vcf ($!)\n";
		outputVCFHeader($outVCF,$g_sampleID);
		#### SNP site ####
		my $nSNP_known=int($dbSNP_SNP_rate*$nSNP);
		my $nSNP_novel=$nSNP-$nSNP_known;
		my ($pKnownSNPSite,$_kSNPNum) = reservoir_sample($snpFile,$nSNP_known); 
		if($_kSNPNum<$nSNP_known) { $nSNP_known = $_kSNPNum; $nSNP_novel=$nSNP-$nSNP_known; }

		my $inFastaFile = $config_hash->{'file_ref'}."/$g_chr.fa";
		my $pSeq=readSeq($g_chr,$inFastaFile);
		my $seq=$pSeq->{'seq'};
		my $g_len=length($pSeq->{'seq'});
		undef $pSeq;

		my $het_hom_ratio = $config_hash->{'het_hom_ratio'};

		foreach (@$pKnownSNPSite)
		{
			my @field=split /\t/;
			my ($chr,$pos,$rid,$ref,$alt)=@field[0,1,2,3,4];
			$alt=(split /,/,$alt)[0];
			if(length($alt)>1 || $alt=~/[^ATCGN]/i || $ref=~/[^ATCGN]/i)
			{
				printf STDERR "strang SNP: $ref/$alt\t($_)\n";
				next;
			}
			if($ref eq 'N') { next; }
			delete $pRegionSite->{"$chr:$pos"};
			### determin genotype
			my ($_z,$_gt,$_haplotype);
			if(rand($het_hom_ratio+1)<1)
			{
				$_z="1/1";
				$_gt=$Abbrev{"$alt$alt"};
				$_haplotype=3;
			}else
			{
				$_z="0/1";
				$_gt=$Abbrev{"$ref$alt"};
				$_haplotype=rand(1)<0.5?1:2;
			}
			printf $outVCF "$chr\t$pos\t$rid\t$ref\t$alt\t.\t.\tVC=GermlineSNP;HAP=$_haplotype;IUPAC=$_gt\tGT\t$_z\n";
		}
		my @aSite=keys %$pRegionSite;
		my $pNovelSNPSite=selection_sample(\@aSite,$nSNP_novel,0);

		foreach (@$pNovelSNPSite)
		{
			my ($chr,$pos)=split /:/;
			delete $pRegionSite->{"$chr:$pos"};
			### determin genotyp
			my ($ref,$alt,$_z,$_gt,$_haplotype);
			$ref=substr($seq,$pos-1,1);
			if($ref eq 'N') { next; }
			$alt=$OTHER_ALLELE{$ref}->[rand(scalar @{$OTHER_ALLELE{$ref}})];
			if(rand($het_hom_ratio+1)<1)
			{
				$_z="1/1";
				$_gt=$Abbrev{"$alt$alt"};
				$_haplotype=3;
			}else
			{
				$_z="0/1";
				$_gt=$Abbrev{"$ref$alt"};
				$_haplotype=rand(1)<0.5?1:2;
			}
			printf $outVCF "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tVC=GermlineSNP;HAP=$_haplotype;IUPAC=$_gt\tGT\t$_z\n";
		}
		undef @aSite;
		undef @$pNovelSNPSite;
		undef @$pKnownSNPSite;
		#### INDEL site ####
		my $nINDEL_known=int($dbSNP_INDEL_rate*$nINDEL);
		my $nINDEL_novel=$nINDEL-$nINDEL_known;
		my ($pKnownINDELSite,$_kINDELNum) = reservoir_sample($indelFile,$nINDEL_known);
		if($_kINDELNum<$nINDEL_known) { $nINDEL_known = $_kINDELNum; $nINDEL_novel=$nINDEL-$nINDEL_known; }
		foreach (@$pKnownINDELSite)
		{
			my @field=split /\t/;
			my ($chr,$pos,$rid,$ref,$alt)=@field[0,1,2,3,4];
			delete $pRegionSite->{"$chr:$pos"};
			my ($_z,$_haplotype);
			if(rand(1)<0.5)
			{
				$_z="0/1";
				$_haplotype=rand(1)<0.5?1:2;
			}else
			{
				$_z="1/1";
				$_haplotype=3;
			}
			printf $outVCF "$chr\t$pos\t$rid\t$ref\t$alt\t.\t.\tVC=GermlineINDEL;HAP=$_haplotype\tGT\t$_z\n";
		}
		@aSite=keys %$pRegionSite;
		my $pNovelINDELSite=selection_sample(\@aSite,$nINDEL_novel,0);
		foreach (@$pNovelINDELSite)
		{
			my ($chr,$pos)=split /:/;
			delete $pRegionSite->{"$chr:$pos"};
			my ($ref,$alt,$_l,$_s,$isIns);
			$ref=substr($seq,$pos-1,1);
			my ($_gt,$_z,$_haplotype);
			if(rand(1)<0.5)
			{
				$_z="0/1";
				$_haplotype=rand(1)<0.5?1:2;
			}else
			{
				$_z="1/1";
				$_haplotype=3;
			}
			$isIns=rand(1)<0.5?1:0;
			$_l=int(random_exponential())+1;
			if($isIns == 1)
			{
				#insertion (after the site specified in $_)
				$_s=randomSequence($_l);
				$_s=join("",@$_s);
				printf $outVCF "$chr\t$pos\t.\t$ref\t$ref$_s\t.\t.\tVC=GermlineINDEL;HAP=$_haplotype\tGT\t$_z\n";
			}else
			{
				#deletion (after the site specified in $_)
				$_s=substr($seq,$pos,$_l);
				printf $outVCF "$chr\t$pos\t.\t$ref$_s\t$ref\t.\t.\tVC=GermlineINDEL;HAP=$_haplotype\tGT\t$_z\n";
			}
		}
		undef @aSite;
		undef @$pNovelINDELSite;
		undef @$pKnownINDELSite;
		my $outSite = new IO::Compress::Gzip "$outDir/germline.$g_chr.remainSite.txt.gz" or die "IO::Compress::Gzip failed: $GzipError\n";
		#open $outSite,"| gzip -c > $outDir/germline.$g_chr.remainSite.txt.gz" or die "Cann't open file $outDir/germline.$g_chr.remainSite.txt.gz ($!)\n";	
		foreach (keys %$pRegionSite)
		{
			printf $outSite "$_\n";
		}
		undef %$pRegionSite;
		printf "$g_chr\ttotal site: $totSite\tnSNP: $nSNP\tnINDEL: $nINDEL\n";
	}
	printf "========================end:fakeGenome.pl germlineSite(%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
}

sub readVar
{
	my ($infile,$pVarTable,$pCNVTable)=@_;
	#my $pVarTable={};
	my ($in);
	open $in,"<",$infile or die "Cann't open file $infile($!)\n";
	while(<$in>)
	{
		chomp;
		if(/^#/) { next; }
		my @field=split /\t/;
		my ($chr,$pos,$rid,$ref,$alt,$info,$z)=@field[0,1,2,3,4,7,9];
		$alt=(split /,/,$alt)[0];
		my ($VC)=/VC=(.+?)[;\t]/;
		my ($haplotype)=/HAP=(.+?)[;\t]/;
		if($alt =~ /[^ATCGN]/i && $info!~/CNV/)
		{
			printf "strange variant: alt($alt)\t$_\n";
			next;
		}
		#chrX    135566218       rs104894773     T       A       .       .       VC=GermlineSNP;IUPAC=W  GT      0/1
		#chrX    24771901        rs11573435      CT      C,CTCTC .       .       VC=GermlineINDEL        GT      1/1
		#chr22   22978075        .       C       T       .       .       VC=SomaticSNV;IUPAC=T   GT      1/1
		#chr22   35494153        .       T       TG      .       .       VC=SomaticINDEL GT      0/1
		#chr22   31223930        .       G       <CNV>   .       .       VC=SomaticCNV;END=31225018;SVLEN=1088   GT:CN   ./.:0
		if($info=~/SNP|SNV/)
		{
			$pVarTable->{$pos}={'ref'=>$ref,'alt'=>$alt,'z'=>$z,'haplotype'=>$haplotype,'type'=>'SNP','VC'=>$VC,'CN'=>2,'rid'=>$rid};
		}elsif($info=~/INDEL/)
		{
			if(length($alt)>length($ref))
			{
				## ins
				$alt=~/.(.+)/;
				my $_l=length($1);
				$pVarTable->{$pos}={'ref'=>$ref,'alt'=>$1,'z'=>$z,'haplotype'=>$haplotype,'op'=>'+','l'=>$_l,'type'=>'Ins','VC'=>$VC,'CN'=>2,'rid'=>$rid};
			}else
			{
				## del
				$ref=~/.(.+)/;
				my $_l=length($1);
				$pVarTable->{$pos}={'ref'=>$alt,'alt'=>$1,'z'=>$z,'haplotype'=>$haplotype,'op'=>'-','l'=>$_l,'type'=>'Del','VC'=>$VC,'CN'=>2,'rid'=>$rid};

			}
		}elsif($info=~/CNV/)
		{
			my ($_end,$_len)= $info=~/END=(\d+?);SVLEN=(\d+)/;
			my $CN=(split /:/,$z)[1];
			$pCNVTable->{$pos}={'ref'=>$ref,'alt'=>$alt,'haplotype'=>$haplotype,'l'=>$_len,'end'=>$_end,'CN'=>$CN,'type'=>'CNV','VC'=>$VC,'rid'=>$rid};
		}
	}
	#return $pVarTable;
}
sub readRegionFile
{
	my ($infile,$TotOnly)=@_;
	my ($in);
	my %returnHash=();
	my $returnTot=0;
	open $in,$infile or die "Cann't open file $infile ($!)\n";
	while(<$in>)
	{
		chomp;
		my @field=split /\t/;
		my ($chr,$beg,$end)=@field;
		$beg++;		## bed format
		for my $_i ($beg..$end)
		{
			if(!$TotOnly) { $returnHash{"$chr:$_i"}=1; }
			$returnTot++;
		}
	}
	return (\%returnHash,$returnTot);
}

sub readList
{
	my ($infile)=@_;
	my ($in);
	my @returnArray=();
	my $returnTot=0;
	open $in,$infile or die "Cann't open file $infile ($!)\n";
	while(<$in>)
	{
		chomp;
		push @returnArray,$_;
		$returnTot++;
	}
	return (\@returnArray,$returnTot);
}

sub randomSequence
{
	my ($len)=@_;
	#my @ALPHABET=('A','T','C','G','M','R','W','S','Y','K');
	my @ALPHABET=('A','T','C','G');
	my @result=();
	for my $i (1..$len)
	{
		push @result,$ALPHABET[rand(4)];
	}
	return \@result;
	#return sprintf "%s",join("",@result);
}
sub sample
{
	my ($num,$beg,$end)=@_;
	#$num=int($num);
	die "Too few elements (".($end-$beg+1)." to select $num from\n" unless $num<($end-$beg+1);
	die "\$end($end) must larger than \$beg($beg)\n" unless $beg<$end;
	my @result=();
	my $pos=$beg;
	while(@result<$num)
	{
		while (rand($end-$pos+1)>($num-@result))
		{
			$pos++;
		}
		push @result,$pos++;
	}
	return \@result;
}

sub selection_sample 
{
	my ($array,$num,$copy)=@_;
	return undef unless defined($copy);
	die "Too few elements (".scalar(@$array).") to select $num from\n" unless $num<@$array;
	my @result;
	my $pos=0;
	while (@result<$num) 
	{
		$pos++ while (rand(@$array-$pos)>($num-@result));
		push @result,$array->[$pos];
		if($copy == 1)
		{
			$pos++;
		}elsif($copy == 0)
		{
			splice @$array,$pos,1;
		}
	}
	return \@result;
}

sub reservoir_sample 
{
	my ($file,$num,$remain)=@_;
	my @buffer;
	my $in;
	if($file=~/\.gz$/) { open $in,"bgzip -cd $file |" or die "Cann't open file $file($!)\n"; }
	else { open $in,$file or die "Cann't open file $file($!)\n"; }
	while(<$in>) 
	{
		chomp;
		if(/^#/) { next; }
		push @buffer,$_;
		last if @buffer==$num;
	}
	if(@buffer<$num) 
	{ 
		#printf STDERR "Insufficient records\n"; 
		return (\@buffer,scalar @buffer);
	}
	my $pos=@buffer;
	while(<$in>) 
	{
		chomp;
		if(/^#/) { next; }
		$pos++;
		my $rand=rand($pos);
		if ($rand<@buffer) {
			if($remain)
			{
				printf $remain "%s\n",$buffer[$rand];
			}
			$buffer[$rand]=$_;
		}
	}
	return (\@buffer,$num);
}

sub outputVCFHeader
{
	my ($out_vcf,$_sampleID)=@_;
	printf $out_vcf "##fileformat=VCFv4.1\n";
	printf $out_vcf "##INFO=<ID=HAP,Number=1,Type=Integer,Description=\"variant haplotype(0--no variant,1--in haplotype 1,2--in haplotype2,3-- in both)\">\n";
	printf $out_vcf "##INFO=<ID=CNVHAP,Number=1,Type=Integer,Description=\"haplotype of CNV\">\n";
	printf $out_vcf "##INFO=<ID=IUPAC,Number=1,Type=String,Description=\"IUPAC genotype\">\n";
	printf $out_vcf "##INFO=<ID=VC,Number=1,Type=String,Description=\"Variation Class\">\n";
	printf $out_vcf "##INFO=<ID=Mutation,Number=1,Type=String,Description=\"Mutation\">\n";
	printf $out_vcf "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
	printf $out_vcf "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
	printf $out_vcf "##ALT=<ID=CNV,Description=\"Copy number variable region\">\n";
	printf $out_vcf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	printf $out_vcf "##FORMAT=<ID=CN,Number=1,Type=String,Description=\"CopyNumber\">\n";
	printf $out_vcf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",$_sampleID;
}
sub write_fasta_seq
{
	my $seq_p=shift;
	my $id=shift;
	my $out=shift;
	my $num_line=(@_) ? shift : 60; ##set the number of charcters in each line
	my $disp;
	printf $out ">$id\n";
	if(ref($seq_p) eq "ARRAY")
	{
		my $i=0;
		for (;$i+$num_line<=@$seq_p;$i+=$num_line)
		{
			my $t=$i+$num_line-1;
			printf $out "%s\n",join("",@{$seq_p}[$i..$t]);
		}
		printf $out "%s\n",join("",@{$seq_p}[$i..@$seq_p-1]);
	}elsif(ref($seq_p) eq "SCALAR")
	{
		for (my $i=0; $i<length($$seq_p); $i+=$num_line) 
		{
			printf $out "%s",substr($$seq_p,$i,$num_line)."\n";
		}
	}
}
sub readSeq
{
	#parameter:chr,fileName
	#read the sequence of chromosome (chr) from fasta file(fileName), return a refHash.
	my ($chr,$fileName)=@_;
	open INFILE,$fileName or die "Cann't open the file:$fileName($!)\n";
	my $reading=0;
	#reset the refHash
	my %hash=('chr'=>"",'seq'=>"");
	#read the data and set the refHash
	while(<INFILE>)
	{
		chomp;
		if(/^>$chr$/ && $reading==0) {$reading=1;next;}
		if(/^>.+/ && $reading==1) {last;}
		if($reading) { $hash{'seq'}.="\U$_"; }
		
	}
	$hash{'chr'}=$chr;
	close INFILE;
	return \%hash;
}
sub usage
{
	die(qq/
Program: fakeGenome.pl
Version: $version
Contact: Zheng Liangtao <tao2013\@gmail.com>\n
Usage:   fakeGenome.pl <command> [<arguments>]\n
Command: 
         germlineSite		general germline variant sites in vcf
         somaticSite		general somatic variant sites in vcf
         somaticCNV		general somatic CNV regions in vcf
         genFasta		generate fasta file in hapoid
\n/);
}
