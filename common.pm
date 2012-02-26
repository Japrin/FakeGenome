#============================================================================
# Name        		: common.pm
# Author      		: zhenglt (email:tao2013@gmail.com)
# Version     		: v1.00
# Created On  		: 
# Last Modified By	: zhenglt
# Last Modified On	: Fri Jun 11 00:53:31 2011
# Copyright   		: Copyright (C) 2010
# Description 		: 
#============================================================================
package  common;
require      Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(%OTHER_ALLELE %OTHER_ALLELE_TiTv %Abbrev %alphabet %alphaMinus readMapTable originalCoordinate newCoordinate search complement_reverse);      
our $VERSION   = 1.00;         

%OTHER_ALLELE=('A'=>['T','C','G'],
				  'T'=>['A','C','G'],
				  'C'=>['A','T','G'],
				  'G'=>['A','T','C'],
				  'N'=>['A','T','C','G'],
				  );
%OTHER_ALLELE_TiTv=( 'A'=>{'Ti'=>['G'],'Tv'=>['T','C']},
						'T'=>{'Ti'=>['C'],'Tv'=>['A','G']},
				  		'C'=>{'Ti'=>['T'],'Tv'=>['A','G']},
				  		'G'=>{'Ti'=>['A'],'Tv'=>['T','C']},
					);
%Abbrev=('AA'=>'A',
			'CC'=>'C',
			'GG'=>'G',
			'TT'=>'T',
			'AC'=>'M',
			'CA'=>'M',
			'AG'=>'R',
			'GA'=>'R',
			'AT'=>'W',
			'TA'=>'W',
			'CG'=>'S',
			'GC'=>'S',
			'CT'=>'Y',
			'TC'=>'Y',
			'GT'=>'K',
			'TG'=>'K',
			'A-'=>'B',
			'-A'=>'B',
			'C-'=>'D',
			'-C'=>'D',
			'G-'=>'H',
			'-G'=>'H',
			'T-'=>'V',
			'-T'=>'V',
			);
%alphabet=(  'A'=>['A','A'],
				'T'=>['T','T'],
				'C'=>['C','C'],
				'G'=>['G','G'],
				'M'=>['A','C'],
				'R'=>['A','G'],
				'W'=>['A','T'],
				'S'=>['C','G'],
				'Y'=>['C','T'],
				'K'=>['G','T'],
				'B'=>['A','-'],
				'D'=>['C','-'],
				'H'=>['G','-'],
				'V'=>['T','-'],
				'-'=>['-','-'],
			);
%alphaMinus=(	 'M'=>{'A'=>'C','C'=>'A'},
				 'R'=>{'A'=>'G','G'=>'A'},
				 'W'=>{'A'=>'T','T'=>'A'},
				 'S'=>{'C'=>'G','G'=>'C'},
				 'Y'=>{'C'=>'T','T'=>'C'},
				 'K'=>{'G'=>'T','T'=>'G'},
				 'A'=>{'A'=>'B','-'=>'A'},
				 'C'=>{'C'=>'D','-'=>'C'},
				 'G'=>{'G'=>'H','-'=>'G'},
				 'T'=>{'T'=>'V','-'=>'T'},
				 '-'=>{'-'=>'-','-'=>'-'},
			);
#######################################
##### used for coordinate mapping #####
#######################################
sub readMapTable
{
	my ($infile)=@_;
	my $pTable1={};
	my $pTable2={};
	my ($in);
	open $in,"<",$infile or die "Cann't open file $infile ($!)\n";
	while(<$in>)
	{
		chomp;
		my @field=split /\t/;
		#test_400000_450000_chr22        2217892 2236857 +1      1/1     2217998 2236964
		my ($chr,$oBeg,$oEnd,$op,$z,$vBeg,$vEnd)=@field;
		#table used for searching
		if(!defined($pTable1->{$chr})) { $pTable1->{$chr}=[]; }
		push @{$pTable1->{$chr}},$vBeg,$vEnd;
		#table used for coordinate mapping
		if(!defined($pTable2->{$chr})) { $pTable2->{$chr}={}; }
		$pTable2->{$chr}->{"$vBeg-$vEnd"}={oBeg=>$oBeg,oEnd=>$oEnd,op=>$op,z=>$z};
	}
	return ($pTable1,$pTable2);
}
sub originalCoordinate
{
	my ($mapTable,$chr,$beg,$end,$pos)=@_;
	if(!($beg<=$pos && $pos<=$end)) { warn "pos($pos) not in range($chr:$beg-$end)\n"; return; }
	my $pInfo=$mapTable->{$chr}->{"$beg-$end"};
	my ($_oBeg,$_oEnd,$_op,$_z)=($pInfo->{'oBeg'},$pInfo->{'oEnd'},$pInfo->{'op'},$pInfo->{'z'});
	#first range, directly map
	if($_op eq ".") { return $pos; }
	my $_dist=$pos-$beg;
	my $_T=abs($_op);
	if($_op>0)
	{
		#ins
		if($_dist<=$_T) { return $_oBeg; }
		else { return $_dist-$_T+$_oBeg; }

	}else
	{
		#del
		if($_z eq "1/1")
		{
			#hom
			return $_oBeg+$_dist+$_T;
		}elsif($_z eq "0/1")
		{
			#het
			return $_oBeg+$_dist;
		}
	}
}
sub newCoordinate
{
	my ($mapTable,$chr,$beg,$end,$pos)=@_;
	if(!($beg<=$pos && $pos<=$end)) { warn "pos($pos) not in range($chr:$beg-$end)\n"; return; }
	my $pInfo=$mapTable->{$chr}->{"$beg-$end"};
	my ($_nBeg,$_nEnd,$_op,$_z)=($pInfo->{'nBeg'},$pInfo->{'nEnd'},$pInfo->{'op'},$pInfo->{'z'});
	#first range, directly map
	if($_op eq ".") { return $pos; }
	my $_dist=$pos-$beg;
	my $_T=abs($_op);
	if($_op>0)
	{
		#ins
		return $_nBeg+$_dist+$_T;
	}else
	{
		#del
		if($_z eq "1/1")
		{
			#hom
			if($_dist<$_T) { return $_nBeg; }
			else { return $_nBeg+$_dist-$_T; }
		}elsif($_z eq "0/1")
		{
			#het
			return $_nBeg+$_dist;
		}
	}
}

sub search
{
	my ($list,$chr,$pos)=@_;
	my $pAry=$list->{$chr};
	if(!$pAry) { return 0; }
	my ($iBegin,$iEnd)=(0,@$pAry/2-1);
	my $iP=int(($iBegin+$iEnd)/2);
	while($iBegin<=$iEnd && $iBegin<= $iP && $iP<=$iEnd)
	{
		if($pAry->[$iP*2]<=$pos && $pos<=$pAry->[$iP*2+1])
		{
			#printf "$pos in [$pAry->[$iP*2]: $pAry->[$iP*2+1]]\n";
			return ($pAry->[$iP*2],$pAry->[$iP*2+1]);
		}elsif($pos < $pAry->[$iP*2])
		{
			$iEnd=$iP-1;
			$iP=int(($iBegin+$iEnd)/2);
		}elsif($pos > $pAry->[$iP*2+1])
		{
			$iBegin=$iP+1;
			$iP=int(($iBegin+$iEnd)/2);
		}
	}
	printf STDERR "Range not found contain pos:\t$chr\t$pos\n";
	return 0;
}
##########################################################################################################
##usage: Complement_Reverse(\$seq);
##########################################################################################################
sub complement_reverse
{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') 
	{ 	##deal with sequence
		$$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;
		$$seq_p=reverse($$seq_p);  
	}
}

1;
