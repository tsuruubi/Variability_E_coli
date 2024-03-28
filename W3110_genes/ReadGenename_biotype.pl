# File information
# /W3110_genes
# perl ReadGenename.pl

#use strict;
#use warnings;


#############################################
#Files
$reffile="AP009048.gff3"; #from genbank AP012306.gff3
$outfile="W3110BiotypeList.txt";
#$outfile_temp="MDS42CDSList_temp.txt";
#############################################

open (REF,"<$reffile") or die "$!";
my @Alltemp=();
my %Biotype=();
my %Position=();
$line=<REF>;
$line=<REF>;
$line=<REF>;
while($line=<REF>)
{
	chomp $line;
	@b=split(/\t/,$line);
	if ($b[0]=~/\#\#/) # ##FASTA
	{
		last;
	}
	if ($b[2]eq"gene")
	{
#		my $locus_tag=$b[8];
#		$locus_tag=~/locus_tag=(.*?);/;
#		$locus_tag=$1;
		my $name_1=$b[8];
		$name_1=~/ID=(.*?);/;
		$name_1=$1;
		push(@Alltemp,$name_1);
		my $name_2=$b[8];
		$name_2=~/gene_biotype=(.*?)$/;
		$name_2=$1;
		$Biotype{$name_1}=$name_2;
		$Position{$name_1}=join($b[3],"_",$b[4]);
#		my $product=$b[8];
#		$product=~/product=(.*?);/;
#		$product=$1;
#		my $strand="C";
#		if ($b[6]eq"+")
#		{
#			$strand="D";
#		}
#		my $parent=$b[8];
#		$parent=~/Parent=gene\-(.*?);/;
#		$parent=$1;
#		print "$parent\n";
	}
}
my @Alluniq=grep{!$count{$_}++} @Alltemp;
my (%tmp, %tmp2);
my @Alldup= grep { $tmp2{$_}++ < 1; } (grep { $tmp{$_}++ >= 1; } (@Alltemp));

close REF;



################################################

open (REF,"<$reffile") or die "$!";
open (OUT,">$outfile") or die "$!";
print OUT "Start\tEnd\tBiotype\tStrand\tName\tProduct\tBnb\tJWnb\tECKnb\n";
$line=<REF>;
$line=<REF>;
$line=<REF>;
my $nb=0;



while($line=<REF>)
{
	chomp $line;
	@b=split(/\t/,$line);
	if ($b[0]=~/\#\#/) # ##FASTA
	{
		last;
	}
	if ($b[2]ne"gene" && $b[2]ne"pseudogene")
	{
#		my $locus_tag=$b[8];
#		$locus_tag=~/locus_tag=(.*?);/;
#		$locus_tag=$1;
		my $parent=$b[8];
		$parent=~/Parent=(.*?);/;
		$parent=$1;
#		if (grep {$_ eq $name} @Alldup)
#		{
#			print "$parent\n";
#		}
		$pos=join($b[3],"_",$b[4]);
		unless (grep {$_ eq $parent} @Alldup)
		{
			if($Biotype{$parent}ne"" && $Position{$parent}eq$pos){
				my $gene=$parent;
				$gene=~s/gene-//;
				my $gene_biotype=$b[8];
				$gene_biotype=~/gene_biotype=(.*?);/;
				$gene_biotype=$1;
				my $product=$b[8];
				$product=~/product=(.*?);/;
				$product=$1;
				my $strand="C";
				if ($b[6]eq"+")
				{
					$strand="D";
				}
		#		print "$parent\n";
				my $gene_synonym=$b[8];
				$gene_synonym=~/Note=(.*?)[;]/;
				$gene_synonym=$1;
				my @BNumb=();
				while ($gene_synonym=~/b(\d{4})/g)
				{
					push (@BNumb,$&);
				}
				$BNumb_join=join('+',@BNumb);
				my @JWNumb=();
				while ($gene_synonym=~/JW(R?)(\d{4})/g)
				{
					push (@JWNumb,$&);
				}
				$JWNumb_join=join('+',@JWNumb);
				my @ECKNumb=();
				while ($gene_synonym=~/ECK(\d{4})/g)
				{
					push (@ECKNumb,$&);
				}
				$ECKNumb_join=join('+',@ECKNumb);
				print OUT "$b[3]\t$b[4]\t$Biotype{$parent}\t$strand\t$gene\t$product\t$BNumb_join\t$JWNumb_join\t$ECKNumb_join\n";
				$nb++;
			}
		}
	}
}

close REF;
close OUT;
