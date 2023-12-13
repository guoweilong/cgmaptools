#!/usr/bin/perl -w

# cgmaptools - mCTanghulu.pl
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
#use Switch 'fallthrough';

=pod

=head1 DESCRIPTION

    Circle plot representing DNA methylation of each C [defualt CpG] site
    on each mapped reads.

=head1 USAGE

    cgmaptools tanghulu [options] -r <ref> -b <bam> -l chr1:133-144
    or: cgmaptools tanghulu [options] -r <ref> -b <bam> -l chr1:133
    (aka mCTanghulu)

    Options:
    -r    Samtools indexed reference genome seqeunce, fasta format. eg. hg19.fa
          - use samtools to index reference: samtools faidx <hg19.fa>
    -b    Samtools indexed Bam file to view.
          - use samtools to index bam file: samtools index <input.bam>
    -l    Region in which to display DNA methylation.
          - or specify a single position (eg. heterozygous SNP site), we will show allele specific methylation.
    -s    Path to samtools eg. /home/user/bin/samtools
          - by defualt, we try to search samtools in your system PATH.
    -o    Output results to file [default: CirclePlot.Ctype.region.Date.pdf].
    -t    C context. [default: CG]
          - available context: C, CG, CH, CW, CC, CA, CT, CHG, CHH
    -d    Ouput device. [default: pdf]
          - alternative: png
    -c    Seperate reads by chain. [default: OFF]
          - specify this option to turn ON.
    -v    Show vague allele linked reads. [ default: OFF]
    -g    Genotype of heterozygous SNP site.
          - This option provides two alleles of htSNP site. eg. AT
          - The genotype information can be used to reduce vague alleles.
          - This option is specific to display methylation in allele specific mode.
    -D    Minimum number of reads (depth) covered in this region or allele linked. [default: 0|OFF]
    -C    Minimum number of C (specified type) covered in this region or allele linked. [default: 0|OFF]
    -W    Width of graphics reigon in inches. [default: 4]
    -H    Height of graphics reigon in inches. [default: 4]
    -R    Resolution in ppi. [default: 300]
          - only available for png device.
    -h    Help message.

=head1 AUTHOR

    Contact:     Zhu, Ping; pingzhu.work@gmail.com
    Last update: 2016-12-07

=cut

## Parsing arguments from command line
my ( $ref, $BAM, $region, $samtools, $out, $Ctype, $device, $byChain, $vague, $gt, $depth, $Cnum, $figure_width, $figure_height, $resolution, $help );

GetOptions(
    'r=s' => \$ref,
    'b=s' => \$BAM,
    'l=s' => \$region,
    's:s' => \$samtools,
    'o:s' => \$out,
    't:s' => \$Ctype,
    'd:s' => \$device,
    'c'   => \$byChain,
	'v'   => \$vague,
    'g:s' => \$gt,
	'D:i' => \$depth,
	'C:i' => \$Cnum,
    'W:i' => \$figure_width,
    'H:i' => \$figure_height,
    'R:i' => \$resolution,
    'h|help' => \$help
);

## Print usage  
pod2usage( { -verbose => 2, -output => \*STDERR } ) if ( $help );
( $ref and $BAM and $region ) or pod2usage();


## Set default 
$samtools ||= `which samtools`;
chomp $samtools;
$samtools or pod2usage();
$out ||= "";
$Ctype ||= "CG";
$device ||= "pdf";
$byChain ||= 0;
$vague ||= 0;
$depth ||= 0;
$Cnum ||= 0;
# default: automatic width, height
$figure_width ||= 0;
$figure_height ||= 0;
$resolution ||= 300;

# check parameters
my %Ctype_support = (
    "C" => 1, "CG" => 1, "CH" => 1, "CW" => 1, "CC" => 1,
    "CA" => 1, "CT" => 1, "CHG" => 1, "CHH" => 1
);
( exists $Ctype_support{$Ctype} ) or pod2usage("\nError: C context $Ctype is not supported!\n");
( $device =~ /pdf|png/ ) or pod2usage("\nError: Device $device is not supported!\n");

if ( $out ne "" ){ # trim suffix
	$out =~ s/\.pdf|\.png$//;
}

my ( $chr, $start, $end ) = split /:|-/,$region;

# redefine region if SNP position specified
my $snp_flag = 0;
my $snp_reference = "";
if ( ! defined($end) ){ # show allele specific methylation
	$end = $start;
	$region = "$chr:$start-$end";
    my ( $start_region, $end_region ) = find_region_for_snp();
    $snp_flag = 1;
    $region = "$chr:$start_region-$end_region";
}

# find C sites on reference genome sequence
my %Csites;
findC_on_reference();

# extract reads in region or allele linked reads
my %c_met;
my $min_cover = 0;
my $max_cover = 0;
c_met_in_region();

# check number of reads
if ( $depth > 0 ){
	my $depth_plus = exists ( $c_met{"+"} ) ? keys %{$c_met{"+"}} : 0;
	my $depth_minus = exists ( $c_met{"-"} ) ? keys %{$c_met{"-"}} : 0;
	my $depth_cover = $depth_plus + $depth_minus;
	( $depth_cover >= $depth ) or die("Number of reads covered $depth_cover less than $depth.\n");
}

# output 
my @covered_c = keys %c_met;


my $output;
if (@covered_c > 0){
    #print "ReadID";
    if ( $snp_flag ) { # add allele position
        $Csites{$start} = "C";
    }
    $output .= "Strand,ReadID";
    my @sites = sort {$a<=>$b} keys %Csites;
    #print "@sites\n";
    # for CG context, symmetrical on +/- strands
    # record it as a CG dinucleotide mode
    if ( $Ctype eq "CG" ){
        foreach ( @sites ){
            delete $Csites{$_} if ( $Csites{$_} eq "G" );
        }
        @sites = sort {$a<=>$b} keys %Csites;
    }
    
    # genomic track
	my $countC = 0;
    map {if ( ($_ >= $min_cover and $_ <= $max_cover) or $_ == $start ) {$output .= ",$_"; $countC ++}} @sites;
    $output .= ";";

	# check number of C covered
	if ( $Cnum > 0 ){
		( $countC >= $Cnum ) or die("Number of $Ctype $countC covered less than $Cnum.\n");
	}

	# at least 1 C covered reads
	( $countC >= 1 ) or die("No reads covered $Ctype at $region\n");
	if ( $snp_flag ) {
		( $countC >= 2 ) or die("No reads covered $Ctype at $chr:$start\n");
	}

    foreach my $strand_tmp ( "+", "-" ){
        foreach my $id_tmp (keys %{$c_met{$strand_tmp}}){
            $output .= "$strand_tmp,$id_tmp";
            foreach (@sites){
                if ( ($_ >= $min_cover and $_ <= $max_cover) or $_ == $start ){
                    my $met = exists $c_met{$strand_tmp}{$id_tmp}{$_} ? $c_met{$strand_tmp}{$id_tmp}{$_} : "NA";
                    $output .= ",$met";
                }
            }
            $output .= ";";
        }
    }

    #print $output."\n";
    
    # plot circles
my $Rcode = <<EOF;
	# avoid warning messages globally!
	options(warn=-1)
    data <- "$output"
    data.split <- unlist(strsplit(data, split = "[;]"))
    data.split <- sapply(data.split, strsplit, split=",")
    data.split<- as.data.frame(t(as.data.frame(data.split)), stringsAsFactors=F)
    names(data.split) <- data.split[1,]
    data.split <- data.split[-1,] 
    #row.names(data.split) <- data.split[,2]
    #data.split[,1] <- NULL
    data.split[, 1] <- factor( data.split[, 1], levels=c("+", "-"), ordered=T)

	# ignore strand and readID
	# when calculate methylation and position
	ignore.column <- c(1, 2)

	##### sort by allele, methylation, position
	#
    # check snp specification 
    snp <- as.numeric("$snp_flag")
    start <- as.numeric("$start")

	# alleles
	if ( snp ){	
		allele <- as.numeric(data.split[, names(data.split) == start])
		allele.column <- (1:ncol(data.split))[names(data.split) == start]
		# delete allels column
		# data.split <- data.split[ , names(data.split) != start]
		# ignore.column <- c(1, 2, allele.column)
	}
	
	# methylation level of reads
    met <- apply(data.split[, - ignore.column], 1, function(x){
		#mean(as.numeric(x), na.rm=T)
		x <- as.numeric(x)
		total.c <- sum( x == 0 | x == 1, na.rm=T)
		met.c <- sum( x == 1, na.rm=T)
		if ( total.c == 0 ) {
			NA
		} else {
			met.c / total.c 
		}
    })

	# positions on reference of reads
    pos <- unlist(apply(data.split[, - ignore.column], 1, function(x){
          min((1:ncol(data.split))[!is.na(as.numeric(x))])
    }))

	# sort
    if ( snp ) { # by allele
        data.split <- data.split[ order( allele, - met, pos), ]
    } else { # no allele
        data.split <- data.split[ order( - met, pos), ]
    }

    # sort by chain 
	byChain <- as.numeric("$byChain")
    # data.split <- data.split[ order( data.split[,1], met, - pos), ]

    # ctype and region
    ctype <- "$Ctype"
    chr <- "$chr"
    start <- "$start"
    end <- "$end"
    
    # output device
    device.out <- "$device"
    file.out <- "$out"
    if ( file.out == ""){
        file.out <- paste("CirclePlot.", ctype, ".", chr, "_", start, "_", end, ".", Sys.Date(), sep="")
		if ( snp ) {
        	file.out <- paste("CirclePlot.", ctype, ".", chr, "_", start, ".", Sys.Date(), sep="")
		}
    }

    # number of circles 
    xpos <- 1:(ncol(data.split)-2)
    # number of reads
    #read.num <- max(c(sum(data.split[,1] == "+"), sum(data.split[,1] == "-")))
    read.num <- nrow(data.split)

    # compare #reads and #circles to adjust pointsize
    max.num <- max(c(xpos, read.num))

    # change point size
    # or figure width, height
    # according to points and read number 

    # width and height 
    figure.width <- as.numeric("$figure_width")
    figure.height <- as.numeric("$figure_height")

    if ( figure.width == 0 & figure.height == 0 ){ # use automatic default width, height using pointsize as 12
        # default point size 12
        pointsize <- 12 
        
        # defalut figure width, height depend on points and reads
        figure.width <- pointsize * length(xpos) / 72
        figure.height <- pointsize * read.num / 72

        if ( figure.width < 4 ) { # at least 4 inch
            figure.width <- 4
        }
        if ( figure.width > 15 ) { # at most 15 inch
            figure.width <- 15
        }

        if ( figure.height < 4 ) { # at least 4 inch
            figure.height <- 4
        }
        if ( figure.height > 15 ) { # at most 15 inch
            figure.height <- 15
        }

        figure.max <- max( figure.width, figure.height )
        if ( figure.max == 15 ){
            pointsize <- figure.max * 72 / max.num
        }
    } else {
        if ( figure.width != 0 ){ # defined width only, automatic define height
            pointsize <- figure.width * 72 / length(xpos)
            figure.height <- pointsize * read.num / 72
            if ( figure.height > 15 ) { # automatic maximum height
                figure.height <- 15
            }
        } else if ( figure.height != 0 ){ # defined height only
            pointsize <- figure.height * 72 / length(xpos)
            figure.width <- pointsize * length(xpos) / 72
            if ( figure.width > 15 ) { # automatic maximum width
                figure.width <- 15
            }
        } else { # defined both width and height
        }
    }

    # new figure
    if ( device.out == "pdf" ){
        pdf(paste(file.out, ".pdf", sep=""), width = figure.width, height = figure.height, pointsize = pointsize)
    } else if ( device.out == "png" ){
        resolution <- as.numeric("$resolution")
        png(paste(file.out, ".png", sep=""), width = figure.width, height = figure.height, pointsize = pointsize, units="in", res=resolution)
    }

	main.title <- paste(chr, ":", start, "-", end, sep="")
	if ( snp ) {
		main.title <- paste(chr, ":", start, sep="")
	}
    plot(x=0, xlim=c(0, max(xpos)+1), ylim=c(0, read.num + 1), type="n", bty="n", xaxt="n", yaxt="n", xlab = "", main = main.title, ylab = "")

    # genomic track
	allele.cex <- 0.6
    y.genomic <- 0
    col.genomic <- rep("grey", length(xpos))
    #####
    # allele specific only
    if ( snp ){
        ref.base <- "$snp_reference"
        col.genomic[ allele.column - 2 ] <- NA
    }
    #####
    segments(x0=1, y0= y.genomic, x1= max(xpos), y1=y.genomic, col="grey") 
    points(x=xpos, y=rep(y.genomic, length(xpos)), pch=19, col="white")
    points(x=xpos, y=rep(y.genomic, length(xpos)), pch=19, col=col.genomic)

    #####
    # allele specific only
    if ( snp ){
        text(x=xpos[ allele.column - 2 ], y=y.genomic, ref.base, xpd=T, col="grey", cex=allele.cex)
    }
    #####

    text(x=1/max(xpos)*0.4, y=y.genomic, adj=c(1, NA), "Genomic track", cex=0.8, col="black", xpd=T)

    # reads in region
    #for( chain in c("+", "-")){
        #data.split.tmp <- data.split[ data.split[,1] == chain, -c(1,2)]
        data.split.tmp <- data.split[ , -c(1,2)]
        data.split.chain <- data.split[ , 1]
        for(i in 1:nrow(data.split.tmp)){
            #ypos.tmp <- i
            ypos.tmp <- (nrow(data.split.tmp):1)[i]
            #if ( chain == "-" ){
            #    ypos.tmp <- -(nrow(data.split.tmp):1)[i]
            #}

            value.tmp <- as.numeric(data.split.tmp[i, ])
            pch.tmp <- rep(NA, length(value.tmp))
            pch.tmp[ value.tmp == 1] <- 19
            pch.tmp[ value.tmp == 0] <- 21
            col.tmp <- rep("white", length(value.tmp))
            col.tmp[ is.na(value.tmp)] <- NA

            #####
            # separate chain by color
            segments.col <- "black"
            if ( byChain & data.split.chain[i] == "-"){
                segments.col <- "grey"
            }
            #####

            segments(x0= min(xpos[!is.na(value.tmp)]), y0= ypos.tmp, x1= max(xpos[!is.na(value.tmp)]), y1=ypos.tmp, col=segments.col) 
            points(x=xpos, y=rep(ypos.tmp, length(xpos)), pch=19, xpd=T, col=col.tmp)
            points(x=xpos, y=rep(ypos.tmp, length(xpos)), pch=pch.tmp, xpd=T, col="black")

            # alleles
            text(x=xpos[ value.tmp == 2], y=ypos.tmp, "A", xpd=T, cex=allele.cex, col="darkgreen")
            text(x=xpos[ value.tmp == 3], y=ypos.tmp, "C", xpd=T, cex=allele.cex, col="blue")
            text(x=xpos[ value.tmp == 4], y=ypos.tmp, "G", xpd=T, cex=allele.cex, col="orange")
            text(x=xpos[ value.tmp == 5], y=ypos.tmp, "T", xpd=T, cex=allele.cex, col="red")
			# vague allele
			#rect(xleft=xpos[ value.tmp == 2.5] - 0.5, ybottom=rep(ypos.tmp-0.5, length(xpos[ value.tmp == 2.5])), xright=xpos[ value.tmp == 2.5] +0.5, ytop= rep(ypos.tmp+0.5, length(xpos[ value.tmp == 2.5])), xpd=T, cex=allele.cex, col="grey", border="white")
			col.vague <- "grey"
			pch.vague <- 15
			if ( sum(value.tmp == 2.5, na.rm=T) > 0 ){
				xpos.tmp <- as.numeric(na.omit(xpos[value.tmp == 2.5]))
				points(x=xpos.tmp, y=ypos.tmp, pch = pch.vague, xpd=T, col=col.vague)
				text(x=xpos[ value.tmp == 2.5], y=ypos.tmp, "A", xpd=T, cex=allele.cex, col="darkgreen")
			}

			#rect(xleft=xpos[ value.tmp == 3.5] - 0.5, ybottom=rep(ypos.tmp-0.5, length(xpos[ value.tmp == 3.5])), xright=xpos[ value.tmp == 3.5] +0.5, ytop= rep(ypos.tmp+0.5, length(xpos[ value.tmp == 3.5])), xpd=T, cex=allele.cex, col="grey", border="white")
			if ( sum(value.tmp == 3.5, na.rm=T) > 0 ){
				xpos.tmp <- as.numeric(na.omit(xpos[value.tmp == 3.5]))
				points(x=xpos.tmp, y=ypos.tmp, pch = pch.vague, xpd=T, col=col.vague)
				text(x=xpos[ value.tmp == 3.5], y=ypos.tmp, "C", xpd=T, cex=allele.cex, col="blue")
			}

			#rect(xleft=xpos[ value.tmp == 4.5] - 0.5, ybottom=rep(ypos.tmp-0.5, length(xpos[ value.tmp == 4.5])), xright=xpos[ value.tmp == 4.5] +0.5, ytop= rep(ypos.tmp+0.5, length(xpos[ value.tmp == 4.5])), xpd=T, cex=allele.cex, col="grey", border="white")
			if ( sum(value.tmp == 4.5, na.rm=T) > 0 ){
				xpos.tmp <- as.numeric(na.omit(xpos[value.tmp == 4.5]))
				points(x=xpos.tmp, y=ypos.tmp, pch = pch.vague, xpd=T, col=col.vague)
				text(x=xpos[ value.tmp == 4.5], y=ypos.tmp, "G", xpd=T, cex=allele.cex, col="orange")
			}

			#rect(xleft=xpos[ value.tmp == 5.5] - 0.5, ybottom=rep(ypos.tmp-0.5, length(xpos[ value.tmp == 5.5])), xright=xpos[ value.tmp == 5.5] +0.5, ytop= rep(ypos.tmp+0.5, length(xpos[ value.tmp == 5.5])), xpd=T, cex=allele.cex, col="grey", border="white")
			if ( sum(value.tmp == 5.5, na.rm=T) > 0 ){
				xpos.tmp <- as.numeric(na.omit(xpos[value.tmp == 5.5]))
				points(x=xpos.tmp, y=ypos.tmp, pch = pch.vague, xpd=T, col=col.vague)
				text(x=xpos[ value.tmp == 5.5], y=ypos.tmp, "T", xpd=T, cex=allele.cex, col="red")
			}
        }
    #}
	# invisible null device message
    invisible(dev.off())

EOF

open R,"|R --vanilla --slave" or die $!;
print R $Rcode;
close R;

} else {
    pod2usage("\nWarning: No reads covered $Ctype in $region\n\tYou can try another region!\n")
}



sub findC_on_reference {
    # loading sequences in $region
    my ( $chr_tmp, $start_tmp, $end_tmp ) = split /:|-/,$region;
    open IN,"$samtools faidx $ref $region | " or die $!;
    # example data format 
    # `samtools faidx hg19.fa chr7:10022-10044`
    # >chr7:10022-10044
    # accctaaccctaaccctaaccct
    <IN>;
    my $ref_seq = uc join "", <IN>;
    $ref_seq =~ s/\n//g; # remove line break in fasta file
    close IN;

    ########
    # allele specific 
    if ( $snp_flag ){
        $snp_reference = substr( $ref_seq, ( $start - $start_tmp ), 1);
    }
    #########
    
    #print "$ref_seq\n";

    my $length_of_Ccontext = length( $Ctype );

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
            $Csites{$ref_c_pos} = $c_base;
        }
    }
}

sub c_met_in_region {
    my ( $chr_tmp, $start_tmp, $end_tmp ) = split /:|-/,$region;

    # Extract reads in specified region
    open IN,"$samtools view $BAM $region | sort -nk 4 |" or die $!; # sort mapped reads by coordinate
    while(<IN>){
        chomp;
        my @map = split /\s+/,$_;
        my ($id, $flag, $pos, $cigar, $seq) = ($map[0], $map[1], $map[3], $map[5], $map[9]);
        # 140M
        # read map to minus/plus strand
        # determine CT or AG to count for methylation level
        my $strand = $flag & 0x10 ? "-" : "+";
        my ($met_base, $unmet_base) = $strand eq "+" ? ("C", "T") : ("G", "A");
        #print "$strand\t$_\n";

        #my ($position, $met); # covered C position, methylation
        my $seq_pos = 0;
        my $outer_flag = 0; # position within region_tmp flag 

        # parse mapping with Insertion/Deletion
        while($cigar =~ /(\d+)(\D)/g){
            my $match_number = $1;
            my $cigar_flag = $2;
            if ($cigar_flag eq "M"){ # matched region_tmp
                foreach my $tmp_pos (0..($match_number - 1)){ 
                    my $ref_pos_tmp = $pos + $tmp_pos;
                    #my $ref_pos_tmp_potentialG = $ref_pos_tmp - 1;
                    next if ($ref_pos_tmp < $start_tmp);
                    if ($ref_pos_tmp > $end_tmp){
                        $outer_flag = 1;
                        last;
                    }

                    if ( $snp_flag and $ref_pos_tmp == $start ) { # for allele specific display
                        my $acgt = "ACGT";
                        my $read_base = substr($seq, ($seq_pos+$tmp_pos), 1);
                        my $base_pos = index($acgt, $read_base) + 2;
						
                        # ignore indistinguishable reads caused by Bisulfite coversion
                        # + strand, T
                        # - strand, A
                        if ( defined($gt) ){ # defined genotype helps to reduce vague alleles
                            if ( ( $strand eq "+" and $gt =~ /CT|TC/ ) or ( $strand eq "-" and $gt =~ /AG|GA/ )){
                                $base_pos += 0.5; # vague reads are labled by base_pos + 0.5
                                if (! $vague){ # do not show vague reads
                                    next;
                                }
                            } elsif ( $strand eq "+" and $gt =~ /C/ and $read_base eq "T" ){
                                $read_base = "C";
                                $base_pos = index($acgt, $read_base) + 2;
                            } elsif ( $strand eq "-" and $gt =~ /G/ and $read_base eq "A" ) {
                                $read_base = "G";
                                $base_pos = index($acgt, $read_base) + 2;
                            }
                        } else { # undefined genotypes
                            if ( ($strand eq "+" and $read_base eq "T") or ($strand eq "-" and $read_base eq "A")){
                                $base_pos += 0.5; # show vague reads
                                if (! $vague){ # do not show vague reads
                                    next;
                                }
                            }
                        }

                        $c_met{$strand}{$id}{$ref_pos_tmp} = $base_pos; 
                        #$c_met{$strand}{$id}{$ref_pos_tmp} = $read_base; 
                        next;
                    }

                    next if (not exists $Csites{$ref_pos_tmp}); # only focused on Csites
                    next if ($Csites{$ref_pos_tmp} eq "C" and $strand eq "-"); # for C site, only search on "+" strand
                    next if ($Csites{$ref_pos_tmp} eq "G" and $strand eq "+"); # for G site, only search on "-" strand

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
                        if ( $Ctype eq "CG" and $Csites{$ref_pos_tmp} eq "G" ){
                            $ref_pos_tmp --;
							next if ( $snp_flag and $ref_pos_tmp == $start );
                        }
                        # for - strand, CG context, avoid reassign for SNP position
                        next if ( exists $c_met{$strand}{$id}{$ref_pos_tmp} );
                        $c_met{$strand}{$id}{$ref_pos_tmp} = $digit_met;

                        # min covered position in region_tmp
                        if ($min_cover == 0 or $ref_pos_tmp < $min_cover){
                            $min_cover = $ref_pos_tmp;
                        }
                        # max covered position in region_tmp
                        if ($max_cover == 0 or $ref_pos_tmp > $max_cover){
                            $max_cover = $ref_pos_tmp;
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
        if ( $snp_flag ) { # delete reads indistinguishable 
            delete $c_met{$strand}{$id} if ( not exists $c_met{$strand}{$id}{$start} );
        }
    }
    close IN;
}

sub find_region_for_snp {
    open IN,"$samtools view $BAM $region | sort -nk 4 | sed -n -e '1p' -e '\$p' | " or die $!; # sort mapped reads by coordinate
    my ($start_tmp, $end_tmp) = (0) x 2;
    while(<IN>){
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
    return($start_tmp, $end_tmp);
}
