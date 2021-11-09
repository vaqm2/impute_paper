#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use Bio::DB::HTS::Tabix;
use Getopt::Long;
use Statistics::Descriptive;

#####################################################################
##################### PURPOSE #######################################
# Takes a VCF of true genotypes, beagle formatted VCF of impute 
# dosages, the reference VCF used, all gzipped files to compute
# SNP and sample rsquared, rsquared by BIN and an aggregate of
# true genotypes vs imputed dosage for plotting bias
#################### CHANGES MADE ###################################
# Date: April 1 2020
# Takes a single extension string input to name out files
# Generates cumulative counts of dosage vs impute insted of verbose
# Generates SNP and Sample based R2 metrics
# Counts minor alleles based on allele count in reference
# ------------------------------------------------------------------
# Date June 8 2020
# Changed from allele counts to allele frequencies to index r2
# ------------------------------------------------------------------
# Date August 27 2020
# Made changes to calculate per SNP and INDIVIDUAL deviations from 
# truth
# ------------------------------------------------------------------
# Date Feb 24 2021
# Changed tabix lib path to point to conda installation
# Moved to Genome DK
###################### AUTHOR #######################################
# Vivek Appadurai | PhD Student | IBP | vivek.appadurai@regionh.dk
##################################################################### 

my ($gzVcfTruth, $gzVcfImpute, $gzVcfRef, $ext);

GetOptions("truth=s"  => \$gzVcfTruth,
	   "impute=s" => \$gzVcfImpute,
	   "ref=s"    => \$gzVcfRef,
	   "out=s"    => \$ext);

my $usage = "\n"."USAGE: perl $0"." ";
$usage   .= "--truth truth.vcf.gz"." "; 
$usage   .= "--impute imputed.vcf.gz"." ";
$usage   .= "--ref reference.vcf.gz"." ";
$usage   .= "--out OutputPrefix"."\n\n";

if(not defined $gzVcfTruth)
{
    print $usage;
    exit;
}

if(not defined $gzVcfImpute)
{
    print $usage;
    exit;
}

if(not defined $gzVcfRef)
{
    print $usage;
    exit;
}

if(not defined $ext)
{
    print $usage;
    exit;
}

my @truth_samples;
my @impute_samples;
my $impute_sample_index = {};
my $bins                = {};
my $aggregate           = {};
my $samples             = {};

open(IN, "zcat $gzVcfTruth |");

while(<IN>)
{
    chomp($_);
    if($_ =~ /^\#CHROM/)
    {
	@truth_samples = split(/\s+/, $_);
	last;
    }
}

close(IN);

open(IN, "zcat $gzVcfImpute |");

while(<IN>)
{
    chomp($_);
    if($_ =~ /^\#CHROM/)
    {
	@impute_samples = split(/\s+/, $_);

	for my $index(0..$#impute_samples)
	{
	    my $impute_sample = $impute_samples[$index];
	    $impute_sample_index->{$impute_sample} = $index;
	}
	last;
    }
}

for my $index(9..$#truth_samples)
{
    if(exists $impute_samples[$index])
    {
	my $sample = $truth_samples[$index];
	
	$samples->{$sample}->{n}       = 0;
	$samples->{$sample}->{sum_x}   = 0;
	$samples->{$sample}->{sum_y}   = 0;
	$samples->{$sample}->{sum_xy}  = 0;
	$samples->{$sample}->{sum_xsq} = 0;
	$samples->{$sample}->{sum_ysq} = 0;
	$samples->{$sample}->{bias_0}  = "";
	$samples->{$sample}->{bias_1}  = "";
	$samples->{$sample}->{bias_2}  = "";
    }
    else
    {
	next;
    }
}

close(IN);

my $vcf     = Bio::DB::HTS::Tabix->new(filename => $gzVcfImpute);
my $ref_vcf = Bio::DB::HTS::Tabix->new(filename => $gzVcfRef);
my $snp_out = IO::File->new("> $ext.mcorr.txt") || die "ERROR: Cannot create output file $ext.mcorr.txt"."\n";

print $snp_out "CHR"." ";
print $snp_out "POS"." ";
print $snp_out "RsID"." ";
print $snp_out "REF"." ";
print $snp_out "ALT"." ";
print $snp_out "MAC"." ";
print $snp_out "N"." ";
print $snp_out "IMPUTE_Rsq"." ";
print $snp_out "EMPIRICAL_Rsq"." ";
print $snp_out "MEAN_BIAS_0"." ";
print $snp_out "MEDIAN_BIAS_0"." ";
print $snp_out "SD_BIAS_0"." ";
print $snp_out "MEAN_BIAS_1"." ";
print $snp_out "MEDIAN_BIAS_1"." ";
print $snp_out "SD_BIAS_1"." ";
print $snp_out "MEAN_BIAS_2"." ";
print $snp_out "MEDIAN_BIAS_2"." ";
print $snp_out "SD_BIAS_2"."\n";

open(IN, "zcat $gzVcfTruth |");

while(<IN>)
{
    chomp($_);
    if($_ =~ /^\#/)
    {
	next;
    }
    else
    {
	my @lineContents = split(/\s+/, $_);
	my $chrom        = $lineContents[0];
	my $position     = $lineContents[1];
#	my $start        = $position - 1;
	my $rsId         = $lineContents[2];	
	my $ref          = $lineContents[3];
	my $alt          = $lineContents[4];	
	my $imputeQuery  = $vcf->query("$chrom:$position-$position");
	my $refQuery     = $ref_vcf->query("$chrom:$position-$position");
	my $minorAllele  = 1;
	my $maf;
	my $impute_rsq;

	# SNP correlation

	my $snp_sum_x   = 0;
	my $snp_sum_y   = 0;
	my $snp_sum_xy  = 0;
	my $snp_sum_xsq = 0;
	my $snp_sum_ysq = 0;
	my $snp_n       = 0;
	my $bias_0      = "";
	my $bias_1      = "";
	my $bias_2      = "";

	#Get Minor Allele Count from Reference
	
	while(my $refMatch = $refQuery->next)
	{
	    my @refContents  = split(/\s+/, $refMatch);
	    my $refAllele    = $refContents[3];
	    my $altAllele    = $refContents[4];

	    if($ref ne $refAllele || $alt ne $altAllele)
	    {
		next;
	    }
	    else
	    {
		my @infoContents    = split(/\;/, $refContents[7]);
		my $alleleFrequency = $infoContents[2];
		$alleleFrequency    =~ s/^AF\=//;

		if($alleleFrequency > 0.5)
		{
		    $maf = 1 - $alleleFrequency;
		    $minorAllele = 0;
		}
		else
		{
		    $maf = $alleleFrequency;
		}
		last;
	    }
	}						
	    
	while(my $match = $imputeQuery->next)
	{
	    my @matchContents = split(/\s+/, $match);
	    my $match_ref     = $matchContents[3];
	    my $match_alt     = $matchContents[4];
	    $impute_rsq       = $matchContents[7];
	    $impute_rsq       =~ s/\;.*$//;
	    $impute_rsq       =~ s/^DR2\=//;
	    
	    if($ref ne $match_ref || $alt ne $match_alt)
	    {
		next;
	    }
	    else
	    {
		my $bin;

		if($maf == 0)
		{
		    next;
		}
		elsif($maf > 0 && $maf <= 0.0001)
		{
		    $bin = 0.0001;
		}
		elsif($maf > 0.0001 && $maf <= 0.0005)
		{
		    $bin = 0.0005;
		}
		elsif($maf > 0.0005 && $maf <= 0.001)
		{
		    $bin = 0.001;
		}
		elsif($maf > 0.001 && $maf <= 0.005)
		{
		    $bin = 0.005;
		}
		elsif($maf > 0.005 && $maf <= 0.01)
		{
		    $bin = 0.01;
		}
		elsif($maf > 0.01 && $maf <= 0.05)
		{
		    $bin = 0.05;
		}
		elsif($maf > 0.05 && $maf <= 0.1)
		{
		    $bin = 0.1;
		}
		elsif($maf > 0.1 && $maf <= 0.2)
		{
		    $bin = 0.2
		}
		elsif($maf > 0.2 && $maf <= 0.5)
		{
		    $bin = 0.5;
		}
		else
		{
		    print "ERROR: Erroneous MAF: $maf at $chrom:$position:$ref:$alt"."\n";
		    exit;
		}
	    
		for my $index(9..$#lineContents)
		{
		    my $truth_sample = $truth_samples[$index];
		    my $truth        = $lineContents[$index];
		    $truth           =~ s/\:.*$//;
		    my $impute;

		    if($truth eq "./.")
		    {
			next;
		    }
		    else
		    {
			$truth = alleleNum($truth);
		    }

		    if(exists $impute_sample_index->{$truth_sample})
		    {
			my $impute_index     = $impute_sample_index->{$truth_sample};
			my @imputeGtContents = split(/\:/, $matchContents[$impute_index]);
			$impute              = $imputeGtContents[1];
		    }
		    else
		    {
			next;
		    }

		    if($minorAllele == 0) #Change GT and Dosage to minor allele if Alt in VCF is Major Allele
		    {
			$truth  = 2 - $truth;
			$impute = 2 - $impute;
			$impute = sprintf("%0.2f", $impute);
		    }

		    # Update BIN aggregate counts
		    
		    my $aggregate_key = $bin."_".$truth."_".$impute;

		    if(exists $aggregate->{$aggregate_key})
		    {
			$aggregate->{$aggregate_key}->{count}++;			
		    }
		    else
		    {
			$aggregate->{$aggregate_key}->{count} = 1;
		    }

		    # Load BIN counts 
		    
		    if(exists $bins->{$bin})
		    {
			$bins->{$bin}->{n}++;
			$bins->{$bin}->{sum_x}   += $truth;
			$bins->{$bin}->{sum_y}   += $impute;
			$bins->{$bin}->{sum_xy}  += $truth * $impute;
			$bins->{$bin}->{sum_xsq} += $truth * $truth;
			$bins->{$bin}->{sum_ysq} += $impute * $impute;
		    }
		    else
		    {
			$bins->{$bin}->{n}       = 1;
			$bins->{$bin}->{sum_x}   = $truth;
			$bins->{$bin}->{sum_y}   = $impute;
			$bins->{$bin}->{sum_xy}  = $truth * $impute;
			$bins->{$bin}->{sum_xsq} = $truth * $truth;
			$bins->{$bin}->{sum_ysq} = $impute * $impute;
		    }

		    # Load SAMPLE counts

		    $samples->{$truth_sample}->{n}       += 1;
		    $samples->{$truth_sample}->{sum_x}   += $truth;
		    $samples->{$truth_sample}->{sum_y}   += $impute;
		    $samples->{$truth_sample}->{sum_xy}  += $truth * $impute;
		    $samples->{$truth_sample}->{sum_xsq} += $truth * $truth;
		    $samples->{$truth_sample}->{sum_ysq} += $impute * $impute;

		    my $bias = $truth - $impute;
		    $bias    = sprintf("%0.2f", $bias);

		    if($truth == 0)
		    {
			$samples->{$truth_sample}->{bias_0} .= $bias.",";
			$bias_0 .= $bias.",";
		    }
		    elsif($truth == 1)
		    {
			$samples->{$truth_sample}->{bias_1} .= $bias.",";
			$bias_1 .= $bias.",";
		    }
		    elsif($truth == 2)
		    {
			$samples->{$truth_sample}->{bias_2} .= $bias.",";
			$bias_2 .= $bias.",";
		    }

		    # Load SNP counts

		    $snp_n       += 1;
		    $snp_sum_x   += $truth;
		    $snp_sum_y   += $impute;
		    $snp_sum_xy  += $truth * $impute;
		    $snp_sum_xsq += $truth * $truth;
		    $snp_sum_ysq += $impute * $impute;
		}

		# Calculate SNP r^2
		
		my $snp_rsq = rSquared($snp_n, $snp_sum_x, $snp_sum_y, $snp_sum_xy, $snp_sum_xsq, $snp_sum_ysq);

		# Calculate SNP bias

		my $mean_bias_0   = calcMean($bias_0);
		my $median_bias_0 = calcMedian($bias_0);
		my $sd_bias_0     = calcSD($bias_0);
		my $mean_bias_1   = calcMean($bias_1);
		my $median_bias_1 = calcMedian($bias_1);
		my $sd_bias_1     = calcSD($bias_1);
		my $mean_bias_2   = calcMean($bias_2);
		my $median_bias_2 = calcMedian($bias_2);
		my $sd_bias_2     = calcSD($bias_2);		
		
		print $snp_out $chrom." ";
		print $snp_out $position." ";
		print $snp_out $rsId." ";
		print $snp_out $ref." ";
		print $snp_out $alt." ";
		print $snp_out $maf." ";
		print $snp_out $snp_n." ";
		print $snp_out $impute_rsq." ";
		print $snp_out $snp_rsq." ";
		print $snp_out $mean_bias_0." ";
		print $snp_out $median_bias_0." ";
		print $snp_out $sd_bias_0." ";
		print $snp_out $mean_bias_1." ";
		print $snp_out $median_bias_1." ";
		print $snp_out $sd_bias_1." ";
		print $snp_out $mean_bias_2." ";
		print $snp_out $median_bias_2." ";
		print $snp_out $sd_bias_2."\n";
	    }
	}
    }
}

close(IN);
$vcf->close;
$ref_vcf->close;
$snp_out->close;

my $out1_fh = IO::File->new("> $ext.BinAggregate.txt") || die "ERROR: Cannot create file: $ext.Aggregate.txt"."\n";

print $out1_fh "BIN"." ";
print $out1_fh "TRUTH"." ";
print $out1_fh "IMPUTE"." ";
print $out1_fh "COUNT"."\n";

for my $index(keys %$aggregate)
{
    my ($bin, $truth, $impute) = split(/\_/, $index);

    print $out1_fh $bin." ";
    print $out1_fh $truth." ";
    print $out1_fh $impute." ";
    print $out1_fh $aggregate->{$index}->{count}."\n";
}

$out1_fh->close;

my $out2_fh = IO::File->new("> $ext.BinCorr.txt") || die "ERROR: Cannot create file: $ext.BinCorr.txt"."\n";

print $out2_fh "MaxMAF"." ";
print $out2_fh "N"." ";
print $out2_fh "Rsq"."\n";

for my $index(sort {$a <=> $b} keys %$bins)
{
    my $n       = $bins->{$index}->{n};
    my $sum_x   = $bins->{$index}->{sum_x};
    my $sum_y   = $bins->{$index}->{sum_y};
    my $sum_xy  = $bins->{$index}->{sum_xy};
    my $sum_xsq = $bins->{$index}->{sum_xsq};
    my $sum_ysq = $bins->{$index}->{sum_ysq};
    my $rsq     = rSquared($n, $sum_x, $sum_y, $sum_xy, $sum_xsq, $sum_ysq);
    
    print $out2_fh $index." ";
    print $out2_fh $n." ";
    print $out2_fh $rsq."\n";
}

$out2_fh->close;

my $out3_fh = IO::File->new("> $ext.icorr.txt") || die "ERROR: Cannot create file: $ext.icorr.txt"."\n";

print $out3_fh "SAMPLE"." ";
print $out3_fh "N"." ";
print $out3_fh "Rsq"." ";
print $out3_fh "MEAN_BIAS_0"." ";
print $out3_fh "MEDIAN_BIAS_0"." ";
print $out3_fh "SD_BIAS_0"." ";
print $out3_fh "MEAN_BIAS_1"." ";
print $out3_fh "MEDIAN_BIAS_1"." ";
print $out3_fh "SD_BIAS_1"." ";
print $out3_fh "MEAN_BIAS_2"." ";
print $out3_fh "MEDIAN_BIAS_2"." ";
print $out3_fh "SD_BIAS_2"."\n";

for my $index(sort keys %$samples)
{
    my $n             = $samples->{$index}->{n};
    my $sum_x         = $samples->{$index}->{sum_x};
    my $sum_y         = $samples->{$index}->{sum_y};
    my $sum_xy        = $samples->{$index}->{sum_xy};
    my $sum_xsq       = $samples->{$index}->{sum_xsq};
    my $sum_ysq       = $samples->{$index}->{sum_ysq};
    my $rsq           = rSquared($n, $sum_x, $sum_y, $sum_xy, $sum_xsq, $sum_ysq);
    my $mean_bias_0   = calcMean($samples->{$index}->{bias_0});
    my $median_bias_0 = calcMedian($samples->{$index}->{bias_0});
    my $sd_bias_0     = calcSD($samples->{$index}->{bias_0});
    my $mean_bias_1   = calcMean($samples->{$index}->{bias_1});
    my $median_bias_1 = calcMedian($samples->{$index}->{bias_1});
    my $sd_bias_1     = calcSD($samples->{$index}->{bias_1});
    my $mean_bias_2   = calcMean($samples->{$index}->{bias_2});
    my $median_bias_2 = calcMedian($samples->{$index}->{bias_2});
    my $sd_bias_2     = calcSD($samples->{$index}->{bias_2});

    print $out3_fh $index." ";
    print $out3_fh $n." ";
    print $out3_fh $rsq." ";
    print $out3_fh $mean_bias_0." ";
    print $out3_fh $median_bias_0." ";
    print $out3_fh $sd_bias_0." ";
    print $out3_fh $mean_bias_1." ";
    print $out3_fh $median_bias_1." ";
    print $out3_fh $sd_bias_1." ";
    print $out3_fh $mean_bias_2." ";
    print $out3_fh $median_bias_2." ";
    print $out3_fh $sd_bias_2."\n";
}

$out3_fh->close;

sub alleleNum
{
    my ($ref, $alt) = split(/\//, $_[0]);
    my $numAlleles   = $ref + $alt;
    return $numAlleles;
}

sub rSquared
{
    my $n       = $_[0];
    my $sum_x   = $_[1];
    my $sum_y   = $_[2];
    my $sum_xy  = $_[3];
    my $sum_xsq = $_[4];
    my $sum_ysq = $_[5];
    my $a       = $n * $sum_xy;
    my $b       = $sum_x * $sum_y;
    my $c       = (($n * $sum_xsq) - ($sum_x * $sum_x));
    my $d       = (($n * $sum_ysq) - ($sum_y * $sum_y));

    if($c == 0 || $d == 0)
    {
	return "NA";
    }
    else
    {    
	my $r       = ($a - $b) / sqrt($c * $d);
	my $rsq     = $r * $r;
	return $rsq;
    }
}

sub calcMean
{
    my $string = shift;
    $string    =~ s/\,$//;
    my @array  = split(/\,/, $string);

    if(!@array)
    {
	return "NA";
    }
    else
    {
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $mean = $stat->mean();
	$mean    = sprintf("%0.2f", $mean);
	return $mean;
    }
}

sub calcMedian
{
     my $string = shift;
    $string    =~ s/\,$//;
    my @array  = split(/\,/, $string);

    if(!@array)
    {
	return "NA";
    }
     else
     {
	 my $stat = Statistics::Descriptive::Full->new();
	 $stat->add_data(@array);
	 my $median = $stat->median();
	 $median    = sprintf("%0.2f", $median);
	 return $median;
     }
}

sub calcSD
{
    my $string = shift;
    $string    =~ s/\,$//;
    my @array  = split(/\,/, $string);

    if(!@array)
    {
	return "NA";
    }
    else
    {
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $sd = $stat->standard_deviation();
	$sd    = sprintf("%0.4f", $sd);
	return $sd;
    }
}
