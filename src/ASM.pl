#!/usr/bin/perl -w

# cgmaptools - ASM.pl
#
# Copyright (C) Ping Zhu
# Contact: Ping Zhu <pingzhu.work@gmail.com>
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;

=pod

=head1 DESCRIPTION

    Allele specific methylated region/site calling
    * Fisher exact test for site calling.
    * Students' t-test for region calling.

=head1 USAGE 

    cgmaptools asm [options] -r <ref.fa> -b <input.bam> -l <snp.vcf>
    (aka ASM)

    Options:
    -r    Samtools indexed reference genome seqeunce, fasta format. eg. hg19.fa
          - use samtools to index reference first: samtools faidx hg19.fa
    -b    Samtools indexed Bam format file.
          - use samtools to index bam file first: samtools index <input.bam>
    -l    SNPs in vcf file format.
    -s    Path to samtools eg. /home/user/bin/samtools
          - by defualt, we try to search samtools in your system PATH,
    -o    Output results to file. [default: STDOUT]
    -t    C context. [default: CG]
          - available context: C, CG, CH, CW, CC, CA, CT, CHG, CHH
    -m    Specify calling mode. [default: asr]
          - alternative: ass
          - asr: allele specific methylated region
          - ass: allele specific methylated site
    -d    Minimum number of read for each allele linked site to call ass. [default: 3]
          - ass specific.
    -n    Minimum number of C site each allele linked to call asr. [default: 2]
          - asr specific.
    -D    Minimum read depth for C site to call methylation level when calling asr. [default: 1]
          - asr specific.
    -L    Low methylation level threshold. [default: 0.2]
          - allele linked region [or site] with low methylation level should be no greater than this threshold.
    -H    High methylation level threshold. [default: 0.8]
          - allele linked region[or site] with high methylation level should be no less than this threshold.
    -q    Adjusted p value using Benjamini & Hochberg (1995) ("BH" or its alias "fdr"). [default: 0.05]
    -h    Help message.

=head1 AUTHOR

    Contact:     Zhu, Ping; pingzhu.work@gmail.com
    Last update: 2016-12-07

=cut

## Parsing arguments from command line
my ( $ref, $BAM, $vcf, $samtools, $out, $Ctype, $mode, $ass_Cdepth, $asm_Cnum, $asm_Cdepth, $lowMet, $highMet, $fdr, $help );

GetOptions(
    'r=s' => \$ref,
    'b=s' => \$BAM,
    'l=s' => \$vcf,
    's:s' => \$samtools,
    'o:s' => \$out,
    't:s' => \$Ctype,
    'm:s' => \$mode,
    'd:i' => \$ass_Cdepth, 
    'n:i' => \$asm_Cnum,
    'D:i' => \$asm_Cdepth,
    'L:f' => \$lowMet,
    'H:f' => \$highMet,
    'q:f' => \$fdr,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 2, -output => \*STDERR } ) if ( $help );
( $ref and $BAM and $vcf ) or pod2usage();


## Set default 
$samtools ||= `which samtools`;
chomp $samtools;
$samtools or pod2usage();
$out ||= "";
$Ctype ||= "CG";
$mode ||= "asr";
$ass_Cdepth ||= 3;
$asm_Cnum ||= 2;
$asm_Cdepth ||= 1;
$lowMet ||= 0.2;
$highMet ||= 0.8;
$fdr ||= 0.05;

## check parameters
my %Ctype_support = (
    "C" => 1, "CG" => 1, "CH" => 1, "CW" => 1, "CC" => 1,
    "CA" => 1, "CT" => 1, "CHG" => 1, "CHH" => 1
);
( exists $Ctype_support{$Ctype} ) or pod2usage("\nError: C context $Ctype is not supported!\n");

( $mode eq "ass" or $mode eq "asr" ) or pod2usage("\nError: $mode mode is not supported!\n");

( $ass_Cdepth >= 1 ) or pod2usage("\nError: Minimum read depth for each allele linked site to call ass should be no less than 1.\n");

( $asm_Cnum >= 2 ) or pod2usage("\nError: Minimum number of C site each allele linked to call asr should be no less than 2.\n");

( $asm_Cdepth >= 1 ) or pod2usage("\nError: Minimum read depth for C site to call methylation level when calling asr should be no less than 1.\n");

( $lowMet < 1 and $lowMet >= 0 ) or pod2usage("\nError: Low methylation level should range in [0,1).\n");

( $highMet <= 1 and $highMet > 0 ) or pod2usage("\nError: High methylation level should range in (0,1].\n");

( $lowMet < $highMet ) or pod2usage("\nError: High methylation level threshould should be greater than low methylation level threshould.\n");

( $fdr <= 1 and $fdr >= 0 ) or pod2usage("\nError: Adjusted p value should be range in [0,1].\n");

## Methylation data preparation for calling
my $SNP_effective_linked;

## SNPs
my $format_column = 8;
open IN,"$vcf" or die ("No such file $vcf.\n");
while(<IN>){
    chomp;
    next if (/^##/); # skip meta description

    #if (/^#CHR/) { # header line
    #    my @line = split /\s+/,$_;
    #    foreach my $i (0..$#line){
    #        if ( $line[$i] eq "FORMAT" ){
    #            $format_column = $i;
    #            last;
    #        }
    #    }
    #    next;
    #}

    ## data lines
    my @a = split /\s+/,$_;
    my ( $chr, $pos, $ref, $alt, $format, $gt ) = ( $a[0], $a[1], $a[3], $a[4], $a[$format_column], $a[$format_column+1] );

    #print "$pos\n";

    # alleles
    next if (length($ref) > 1); # only retain SNP site
    my @alleles = ($ref); # 0 allele is reference 
    my $indel_flag = 0; # indel flag
    foreach my $nt ( split /,/,$alt){
        if ( length($nt) != 1 or $nt !~ /[ATCG]/){ # scan indel or unresolved sites
            $indel_flag = 1;
            last;
        }
        push ( @alleles, $nt );
    }
    next if ( $indel_flag == 1 ); # only retian SNP site

    # genotype
    my @gt_index = (split /[:|\/|=]/,$gt)[0..1];
    my ($allele1, $allele2) = ( uc $alleles[ $gt_index[0] ], uc $alleles[ $gt_index[1] ] );
    next if ( $allele1 eq $allele2 ); # only retain heterozygous SNP site (htSNP)

    # determine region of SNP linked reads
    my $snp_reference;
    my ( $start, $end ) = find_region_for_snp("$chr:$pos-$pos");
    my $region = "$chr:$start-$end";
    #print "$region\n";

    # find C sites on reference genome sequence
    my %Csites;
    findC_on_reference($region, \%Csites, $pos);

    my @csites = sort {$a<=>$b} keys %Csites;
    #print "@csites\n";

    # methylation of linked sites
    my %c_met;
    my $genotype = join "", sort($allele1, $allele2);
    my ( $min_cover, $max_cover ) = cmet_in_region($pos, $genotype, $region, \%c_met, \%Csites);

    # format linked C site
    my $format_C = format_C_site(\%c_met, \%Csites, $min_cover, $max_cover, $pos);

    #print $format_C."\n";

    ## Allele specific methylation
    # asr: Allele specific methylated region
    if ( $mode eq "asr" ){
        my $this_SNP_prepare = asr_call_prepare($format_C, $pos, $allele1, $allele2);
        if ( $this_SNP_prepare ne "NA" ){
            $SNP_effective_linked .= "$chr,$pos,$ref,$allele1,$allele2,$this_SNP_prepare\n";
        }
    } 
    # ass: Allele specific methylated site
    elsif ( $mode eq "ass" ){
        my $this_SNP_prepare = ass_call_prepare($format_C, $pos, $allele1, $allele2);
        if ( $this_SNP_prepare ne "NA" ){
            my @C_sites = split /;/,$this_SNP_prepare;
            foreach (@C_sites){
                $SNP_effective_linked .= "$chr,$pos,$ref,$allele1,$allele2,$_\n";
            }
        }
    }
}
close IN;

# effective linked C sites
if (defined($SNP_effective_linked)){
    if ( $mode eq "asr" ){
        asr_test( $SNP_effective_linked );
    } elsif ( $mode eq "ass" ){
        ass_test( $SNP_effective_linked );
    }
}
# no effective linked C sites
else {
    pod2usage("\nError: No effective C sites linked to SNPs!\n");
}


## subroutines

sub format_C_site {
    my ( $cmet_tmp, $Csites_tmp, $min_cover_tmp, $max_cover_tmp, $pos_tmp) = @_;
    my $format_C_tmp = "Strand,ReadID";
    my @sites = sort {$a<=>$b} keys %$Csites_tmp;

    # for CG context, symmetrical on +/- strands
    # record them as a CG dinucleotide 
    if ( $Ctype eq "CG" ){
        map { if ( $$Csites_tmp{$_} eq "G" and $_ != $pos_tmp ) { delete $$Csites_tmp{$_} } } @sites; # delete G site at + strand
        @sites = sort {$a<=>$b} keys %$Csites_tmp;
    }

    # remove C sites outside of [ $min_cover_tmp, $max_cover_tmp ]
    map {if ( ($_ >= $min_cover_tmp and $_ <= $max_cover_tmp ) or $_ == $pos_tmp ) { $format_C_tmp .= ",$_" } } @sites;
    $format_C_tmp .= ";";

    foreach my $strand_tmp ( "+", "-" ){
        foreach my $id_tmp ( keys %{$$cmet_tmp{$strand_tmp}} ){
            $format_C_tmp .= "$strand_tmp,$id_tmp";
            foreach (@sites){
                if ( ($_ >= $min_cover_tmp and $_ <= $max_cover_tmp ) or $_ == $pos_tmp ){
                    my $met = exists $$cmet_tmp{$strand_tmp}{$id_tmp}{$_} ? $$cmet_tmp{$strand_tmp}{$id_tmp}{$_} : "NA";
                    $format_C_tmp .= ",$met";
                }
            }
            $format_C_tmp .= ";";
        }
    }

    # return format C site
    return($format_C_tmp);
}


sub asr_call_prepare {
    my ( $format_C_tmp, $pos_tmp, $allele1_tmp, $allele2_tmp ) = @_;
    my @linked_reads = split /;/,$format_C_tmp;

    # find SNP column
    my $snp_column;
    my @header = split /,/,$linked_reads[0];
    foreach my $site ( 2..$#header ){
        if ( $header[$site] == $pos_tmp ){
            $snp_column = $site;
            last;
        }
    }

    # Allele index
    my $acgt = "ACGT";
    my $allele1_pos = index($acgt, $allele1_tmp) + 2;
    my $allele2_pos = index($acgt, $allele2_tmp) + 2;

    # C site
    my %C_stat;
    foreach my $read ( @linked_reads[ 1..$#linked_reads ] ){
        my @C_sites = split /,/,$read;

        # ignore nt not in alleles
        next if ( $C_sites[ $snp_column ] != $allele1_pos and $C_sites[ $snp_column ] != $allele2_pos );

        foreach my $site ( 2..$#C_sites ){
            next if ( $C_sites[ $site ]  eq "NA" or $site == $snp_column ); # ignore NA and SNP column
            $C_stat{$C_sites[ $snp_column ]}{$site} .= "$C_sites[$site],";
        }
    }

    ## check threshold
    # 1. Minimum number of reads to call methylation at C site
    foreach my $allele_index_tmp ( keys %C_stat ){ # foreach allele in ACGT ([ 2,3,4,5 ] positions)
        foreach my $site_tmp ( keys %{$C_stat{$allele_index_tmp}} ){
            $C_stat{$allele_index_tmp}{$site_tmp} =~ s/,$//;
            my @support_met = split /,/,$C_stat{$allele_index_tmp}{$site_tmp};
            if ( @support_met < $asm_Cdepth ){
                delete $C_stat{$allele_index_tmp}{$site_tmp};
            } else {
                my ( $met_sum, $met_this_site );
                map { $met_sum += $_ } @support_met;
                $met_this_site = int( $met_sum / @support_met * 100 ) / 100;
                $C_stat{$allele_index_tmp}{$site_tmp} = $met_this_site;
            }
        }
    }
    # 2. Minimum number of C sites in this region
    my $met_groups;
    foreach my $allele_tmp ( $allele1_pos, $allele2_pos ){
        if ( (not exists $C_stat{$allele_tmp}) or (keys %{$C_stat{$allele_tmp}}) < $asm_Cnum ){
            return("NA");
        } else {
            foreach my $site_tmp ( sort {$a<=>$b} keys %{$C_stat{$allele_tmp}} ){
                $met_groups .= $C_stat{$allele_tmp}{$site_tmp}."-";
            }
        }
        $met_groups =~ s/-$//;
        $met_groups .= ","; # separate each allele with ','
    }
    $met_groups =~ s/,$//;
    return( $met_groups );
}


sub ass_call_prepare {
    my ( $format_C_tmp, $pos_tmp, $allele1_tmp, $allele2_tmp ) = @_;
    my @linked_reads = split /;/,$format_C_tmp;

    # find SNP column
    my $snp_column;
    my @header = split /,/,$linked_reads[0];
    foreach my $site ( 2..$#header ){
        if ( $header[$site] == $pos_tmp ){
            $snp_column = $site;
            last;
        }
    }

    # Allele index
    my $acgt = "ACGT";
    my $allele1_pos = index($acgt, $allele1_tmp) + 2;
    my $allele2_pos = index($acgt, $allele2_tmp) + 2;

    # C site
    my %C_stat;
    foreach my $read ( @linked_reads[ 1..$#linked_reads ] ){
        my @C_sites = split /,/,$read;

        # ignore nt not in alleles
        next if ( $C_sites[ $snp_column ] != $allele1_pos and $C_sites[ $snp_column ] != $allele2_pos );

        foreach my $site ( 2..$#C_sites ){
            next if ( $C_sites[ $site ]  eq "NA" or $site == $snp_column ); # ignore NA and SNP column
            $C_stat{$header[$site]}{$C_sites[ $snp_column ]}{"Met"} ++ if ( $C_sites[$site] == 1 );
            $C_stat{$header[$site]}{$C_sites[ $snp_column ]}{"unMet"} ++ if ( $C_sites[$site] == 0 );
        }
    }

    ## check threshold
    # Minimum number of reads covered at C site
    foreach my $site_tmp ( keys %C_stat ){ # foreach allele in ACGT ([ 2,3,4,5 ] positions)
        foreach my $allele_index_tmp ( keys %{$C_stat{$site_tmp}} ){
            if ( not exists $C_stat{$site_tmp}{$allele_index_tmp}{"Met"} ){
                $C_stat{$site_tmp}{$allele_index_tmp}{"Met"} = 0;
            } 
            if ( not exists $C_stat{$site_tmp}{$allele_index_tmp}{"unMet"} ){
                $C_stat{$site_tmp}{$allele_index_tmp}{"unMet"} = 0;
            }

            # number of reads provide methylation status
            my $support_met = $C_stat{$site_tmp}{$allele_index_tmp}{"Met"} + $C_stat{$site_tmp}{$allele_index_tmp}{"unMet"};

            if ( $support_met < $ass_Cdepth ){
                delete $C_stat{$site_tmp};
                last;
            } else {
                $C_stat{$site_tmp}{$allele_index_tmp} = "$C_stat{$site_tmp}{$allele_index_tmp}{'Met'}-$C_stat{$site_tmp}{$allele_index_tmp}{'unMet'}";
            }
        }
    }
    # prepare return
    my $met_groups;
    foreach my $site_tmp ( sort {$a<=>$b} keys %C_stat){
        my $this_met;
        foreach my $allele_tmp ( $allele1_pos, $allele2_pos ){
            if ( not exists $C_stat{$site_tmp}{$allele_tmp} ){
                return("NA");
            }
            $this_met .= $C_stat{$site_tmp}{$allele_tmp}.","; # separate each allele with ","
        }
        $this_met =~ s/,$//;
        $met_groups .= "$site_tmp,$this_met;";
    }
    
    if ( defined( $met_groups ) ){
        $met_groups =~ s/;$//;
        return( $met_groups );
    } else {
        return("NA");
    }
}

sub asr_test {
    my ( $linked_C_sites ) = @_;
    my $Rcode = <<EOF;
    # avoid warning mesages glabally!
    options(warn=-1)
    data <- read.csv(text="$linked_C_sites" , head=F, sep=",", stringsAsFactors=F, colClasses=c("character", "numeric", "character", "character","character", "character","character"))

    names(data) <- c("Chr", "Pos", "Ref", "Allele1", "Allele2", "Allele1_linked_C", "Allele2_linked_C")
    columnNum <- ncol(data) 
    data\$Allele1_linked_C_met <- sapply( data\$Allele1_linked_C, function(x){
        format(mean(as.numeric(unlist(strsplit( x, split="-")))), digits=2)
    })
    data\$Allele2_linked_C_met <- sapply( data\$Allele2_linked_C, function(x){
        format(mean(as.numeric(unlist(strsplit( x, split="-")))), digits=2)
    })

    data\$pvalue <- apply(data[, c(columnNum-1, columnNum) ], 1, function(x){
        allele1_met <- as.numeric(unlist(strsplit(x[1], split="-")))
        allele2_met <- as.numeric(unlist(strsplit(x[2], split="-")))

        # in case of error "data are essentially constant"
        if ( sd( allele1_met ) == 0 & sd( allele2_met ) == 0 & mean( allele1_met ) == mean( allele2_met ) ){
            return(1)
        } else {
            while (sum(abs(allele1_met - mean(allele1_met))) == 0){
                allele1_met <- as.numeric(abs(jitter(allele1_met, amount = 0.001)))
            }
            return(sprintf("%.2e", t.test(allele1_met, allele2_met)\$p.value))
        }
    })

    # thresholds for ASM
    lowMet <- as.numeric("$lowMet")
    highMet <- as.numeric("$highMet")
    fdr <- as.numeric("$fdr")

    data\$fdr <- sprintf("%.2e", p.adjust(data\$pvalue, method="fdr"))
    data\$ASM <- "FALSE"
    data\$ASM[ data\$Allele1_linked_C_met >= highMet & data\$Allele2_linked_C_met <= lowMet & as.numeric(data\$fdr) <= fdr ] <- "TRUE"
    data\$ASM[ data\$Allele2_linked_C_met >= highMet & data\$Allele1_linked_C_met <= lowMet & as.numeric(data\$fdr) <= fdr ] <- "TRUE"

    file.out <- "$out"
    if (file.out == ""){
        write.table(data, stdout(), quote=F, sep="\t", col.names=T, row.names=F)
    } else {
        write.table(data, file.out, quote=F, sep="\t", col.names=T, row.names=F)
    }


EOF

open R,"|R --vanilla --slave" or die $!;
print R $Rcode;
close R;

}

sub ass_test {
    my ( $linked_C_sites ) = @_;
    my $Rcode = <<EOF;
    # avoid warning mesages glabally!
    options(warn=-1)
    data <- read.csv(text="$linked_C_sites" , head=F, sep=",", stringsAsFactors=F, colClasses=c("character", "numeric", "character", "character","character", "numeric", "character","character"))

    names(data) <- c("Chr", "SNP_Pos", "Ref", "Allele1", "Allele2", "C_Pos", "Allele1_linked_C", "Allele2_linked_C")
    columnNum <- ncol(data) 
    data\$Allele1_linked_C_met <- sapply( data\$Allele1_linked_C, function(x){
        met <- as.numeric( unlist( strsplit( x, split="-") ) )
        sprintf("%.2f", met[1] / ( met[1] + met[2] ))
    })
    data\$Allele2_linked_C_met <- sapply( data\$Allele2_linked_C, function(x){
        met <- as.numeric( unlist( strsplit( x, split="-") ) )
        sprintf("%.2f", met[1] / ( met[1] + met[2]))
    })

    data\$pvalue <- apply(data[, c(columnNum-1, columnNum) ] , 1, function(x){
        met1 <- as.numeric( unlist( strsplit(x[1], split="-") ) )
        met2 <- as.numeric( unlist( strsplit(x[2], split="-") ) )
        sprintf("%.2e", fisher.test(matrix(c(met1, met2), nrow=2), conf.int=F)\$p.value)
    })

    # thresholds for ASM
    lowMet <- as.numeric("$lowMet")
    highMet <- as.numeric("$highMet")
    fdr <- as.numeric("$fdr")

    data\$fdr <- sprintf("%.2e", p.adjust(as.numeric(data\$pvalue), method="fdr"))

    data\$ASM <- "FALSE"
    data\$ASM[ data\$Allele1_linked_C_met >= highMet & data\$Allele2_linked_C_met <= lowMet & as.numeric(data\$fdr) <= fdr ] <- "TRUE"
    data\$ASM[ data\$Allele2_linked_C_met >= highMet & data\$Allele1_linked_C_met <= lowMet & as.numeric(data\$fdr) <= fdr ] <- "TRUE"


    file.out <- "$out"
    if (file.out == ""){
        write.table(data, stdout(), quote=F, sep="\t", col.names=T, row.names=F)
    } else {
        write.table(data, file.out, quote=F, sep="\t", col.names=T, row.names=F)
    }

EOF

open R,"|R --vanilla --slave" or die $!;
print R $Rcode;
close R;
    
}


sub findC_on_reference {
    my ( $region_tmp, $Csites_tmp, $pos_tmp ) = @_;
    # loading sequences in $region_tmp
    my ( $chr_tmp, $start_tmp, $end_tmp ) = split /:|-/,$region_tmp;
    open IN_tmp,"$samtools faidx $ref $region_tmp | " or die $!;
    # example data format 
    # `samtools faidx hg19.fa chr7:10022-10044`
    # >chr7:10022-10044
    # accctaaccctaaccctaaccct
    <IN_tmp>;
    my $ref_seq = uc join "", <IN_tmp>;
    $ref_seq =~ s/\n//g; # remove line break in fasta file
    close IN_tmp;

    my $length_of_Ccontext = length( $Ctype );

    # find reference nt on SNP site
    my $ref_nt_pos = $pos_tmp - $start_tmp;
    my $ref_nt = substr( $ref_seq, $ref_nt_pos, 1 );
    ${$Csites_tmp}{$pos_tmp} = $ref_nt;

    # find C sites
    while( $ref_seq =~ /C|G/g ){
        my $c_pos = pos( $ref_seq ) - 1; 
        my $c_base = substr( $ref_seq, $c_pos, 1 );

        my $c_context; # c context in the reference genome
        if ( $c_base eq "C" ){ # C on reference sequence
            # eg. AAAA'C'GAAAA, here we matched C on the top strand
            #     TTTT'G'CTTTT 
            $c_context = substr( $ref_seq, $c_pos, $length_of_Ccontext );
        } else {
            # eg. AAAAC'G'AAAA, here we matched C on the bottom strand
            #     TTTTG'C'TTTT 
            my $match_start = ( $c_pos - $length_of_Ccontext + 1 ) < 0 ? 0 : ( $c_pos - $length_of_Ccontext + 1 );
            my $substr_length = ( $c_pos - $length_of_Ccontext + 1 ) < 0 ? $c_pos : $length_of_Ccontext;
            $c_context = substr( $ref_seq, $match_start, $substr_length );
            $c_context = reverse( $c_context ); # reverse
            $c_context =~ tr/ACGT/TGCA/; # complement
        }

        my %c_context_sets;
		
		# require Switch.pm installed. discarded
		#switch ( $c_context ) {
        #    case /^C/ { $c_context_sets{"C"} = 1; }
        #    case /CA/ { $c_context_sets{"CA"} = 1; }
        #    case /CC/ { $c_context_sets{"CC"} = 1; }
        #    case /CG/ { $c_context_sets{"CG"} = 1; }
        #    case /CT/ { $c_context_sets{"CT"} = 1; }
        #    case /C[AT]/ { $c_context_sets{"CW"} = 1; }
        #    case /C[ATC]/ { $c_context_sets{"CH"} = 1; }
        #    case /C[ATC]G/ { $c_context_sets{"CHG"} = 1; }
        #    case /C[ATC][ATC]/ { $c_context_sets{"CHH"} = 1; }
        #}
		
        if ( $c_context =~ /^C/ ) { $c_context_sets{"C"} = 1; }
        if ( $c_context =~ /CA/ ) { $c_context_sets{"CA"} = 1; }
        if ( $c_context =~ /CC/ ) { $c_context_sets{"CC"} = 1; }
        if ( $c_context =~ /CG/ ) { $c_context_sets{"CG"} = 1; }
        if ( $c_context =~ /CT/ ) { $c_context_sets{"CT"} = 1; }
        if ( $c_context =~ /C[AT]/ ) { $c_context_sets{"CW"} = 1; }
        if ( $c_context =~ /C[ATC]/ ) { $c_context_sets{"CH"} = 1; }
        if ( $c_context =~ /C[ATC]G/ ) { $c_context_sets{"CHG"} = 1; }
        if ( $c_context =~ /C[ATC][ATC]/ ) { $c_context_sets{"CHH"} = 1; }

        # record this position if the c context fits user defined type
        # default: CG
        if ( exists $c_context_sets{$Ctype} ){
            my $ref_c_pos = $start_tmp + $c_pos;
            ${$Csites_tmp}{$ref_c_pos} = $c_base;
        }
    }
}

sub cmet_in_region {
    my ( $pos_tmp, $genotype_tmp, $region_tmp, $cmet_tmp, $Csites_tmp ) = @_;
    my ( $chr_tmp, $start_tmp, $end_tmp ) = split /:|-/,$region_tmp;

    # region redefine
    # ignore sites outside $region_tmp
    my ( $min_cover_tmp, $max_cover_tmp ) = (0) x 2;

    # Extract reads in specified region_tmp
    open IN_tmp,"$samtools view $BAM $region_tmp | sort -nk 4 |" or die $!; # sort mapped reads by coordinate
    while(<IN_tmp>){
        chomp;
        my @map = split /\s+/,$_;
        my ($id, $flag, $pos, $cigar, $seq) = ($map[0], $map[1], $map[3], $map[5], $map[9]);
        # 140M
        # read map to minus/plus strand
        # determine CT or AG to count for methylation level
        my $strand = $flag & 0x10 ? "-" : "+";
        my ($met_base, $unmet_base) = $strand eq "+" ? ("C", "T") : ("G", "A");

        my $seq_pos = 0;
        my $outer_flag = 0; # position within region_tmp flag 

        # parse mapping with Insertion/Deletion
        while($cigar =~ /(\d+)(\D)/g){
            my $match_number = $1;
            my $cigar_flag = $2;
            if ($cigar_flag eq "M"){ # matched region_tmp_tmp
                foreach my $tmp_pos (0..($match_number - 1)){ 
                    my $ref_pos_tmp = $pos + $tmp_pos;
                    #my $ref_pos_tmp_potentialG = $ref_pos_tmp - 1;
                    next if ($ref_pos_tmp < $start_tmp);
                    if ($ref_pos_tmp > $end_tmp){
                        $outer_flag = 1;
                        last;
                    }

                    next if (not exists $$Csites_tmp{$ref_pos_tmp}); # only focused on Csites_tmp, which includes SNP site

                    # SNP site
                    if ( $ref_pos_tmp == $pos_tmp ) { # for allele specific display
                        my $acgt = "ACGT";
                        my $read_base = substr($seq, ($seq_pos+$tmp_pos), 1);
                        my $base_pos = index($acgt, $read_base) + 2;
						
                        # ignore indistinguishable reads caused by Bisulfite coversion
                        # + strand, T, may derived from C or T
                        # - strand, A, mey derived from G or A
                        if ( ( $strand eq "+" and $genotype_tmp =~ /CT/ ) or ( $strand eq "-" and $genotype_tmp =~ /AG/ )){
							$base_pos += 0.5; # vague reads are labled by base_pos + 0.5
                        } elsif ( $strand eq "+" and $genotype_tmp =~ /C/ and $read_base eq "T" ){
                            $read_base = "C";
                            $base_pos = index($acgt, $read_base) + 2;
                        } elsif ( $strand eq "-" and $genotype_tmp =~ /G/ and $read_base eq "A" ) {
                            $read_base = "G";
                            $base_pos = index($acgt, $read_base) + 2;
                        }

                        $$cmet_tmp{$strand}{$id}{$ref_pos_tmp} = $base_pos; 
                        next; # do not call methylation for SNP position
                    }

                    # C site
                    next if ($$Csites_tmp{$ref_pos_tmp} eq "C" and $strand eq "-"); # for C site, only search on "+" strand
                    next if ($$Csites_tmp{$ref_pos_tmp} eq "G" and $strand eq "+"); # for G site, only search on "-" strand

                    my $read_base = substr($seq, ($seq_pos+$tmp_pos), 1);
                    my $digit_met;
                    if ($read_base eq $met_base){ # met
                        $digit_met = 1;
                    } elsif ($read_base eq $unmet_base){ # unmet
                        $digit_met = 0;
                    } else { # neither equal Met or unMet base
                    }

                    if ( defined($digit_met) ){
                        # for CG context, symmetrical on +/- strands
                        # record it as a CG dinucleotide mode
                        if ( $Ctype eq "CG" and $$Csites_tmp{$ref_pos_tmp} eq "G" ){
                            $ref_pos_tmp --;
							next if ( $ref_pos_tmp == $pos_tmp );
                        }
                        # for - strand, CG context, avoid reassign for SNP position
                        next if ( exists $$cmet_tmp{$strand}{$id}{$ref_pos_tmp} );
                        $$cmet_tmp{$strand}{$id}{$ref_pos_tmp} = $digit_met;

                        # min covered position in region_tmp
                        if ($min_cover_tmp == 0 or $ref_pos_tmp < $min_cover_tmp){
                            $min_cover_tmp = $ref_pos_tmp;
                        }
                        # max covered position in region_tmp
                        if ($max_cover_tmp == 0 or $ref_pos_tmp > $max_cover_tmp ){
                            $max_cover_tmp = $ref_pos_tmp;
                        }
                    }
                }
                last if ($outer_flag == 1); # out of region_tmp
                $pos += $match_number; 
                $seq_pos += $match_number;
            } elsif ($cigar_flag eq "I"){ # insertion
                $seq_pos += $match_number;
            } elsif ($cigar_flag eq "D"){ # deletion
                $pos += $match_number;
            } else { # other cigar character
            }
        }
        # delete reads without coverage at SNP site
        delete $$cmet_tmp{$strand}{$id} if ( not exists $$cmet_tmp{$strand}{$id}{$pos_tmp} );
    }
    close IN_tmp;
    
    return( $min_cover_tmp, $max_cover_tmp);
}

sub find_region_for_snp {
    my ( $region_tmp ) = @_;
    open IN_tmp,"$samtools view $BAM $region_tmp | sort -nk 4 | sed -n -e '1p' -e '\$p' | " or die $!; # sort mapped reads by coordinate
    my ($start_tmp, $end_tmp) = (0) x 2;
    while(<IN_tmp>){
        chomp;
        my @a = split /\s+/,$_;
        my $pos = $a[3];
        my $length = length($a[9]);
        my $end_this = $pos + $length - 1;
        if ($start_tmp == 0) {
            $start_tmp = $pos;
        }
        if ($end_tmp < $end_this){
            $end_tmp = $end_this;
        }
    }
    close IN_tmp;
    return($start_tmp, $end_tmp);
}
