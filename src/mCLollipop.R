#!/usr/bin/env Rscript



# cgmaptools - mCLollipop
#
# Copyright (C) Weilong Guo
# Contact: Weilong Guo <guoweilong@126.com>
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

# Guo, Weilong; guoweilong@126.com; 2015-05-07

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("optparse")){
  warning("Detect package \"optparse\" is not installed in your R enviroment.")
  warning("Trying to install the \"optparse\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("optparse")
}

# Argument
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file, use STDIN if ommited, multiple-chr is not suggested"),
  make_option(c("-a", "--annotation"), dest = "annofile", default = "",
              help="[opt] annotation file name, refFlat format"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file"),
  make_option(c("-f", "--format"), dest = "format", default = "pdf",
              help = "[opt] the format for output figure: pdf (default), png, eps"),
  make_option(c("-l", "--left"), dest = "left", default = "",
              help = "[opt] Left-most position, use the 1st position if omitted"),
  make_option(c("-r", "--right"), dest = "right", default = "",
              help = "[opt] Right-most position, use the last position of input if omitted"),
  make_option(c("-c", "--chr"), dest = "chr", default = "",
              help = "[opt] chromosome name, use the chr in 1st line of input file if omitted"),
  make_option(c("-s", "--site"), dest = "sitefile", default = "", 
              help = "[opt] file of site to be marked"),
  make_option(c("-b", "--bed"), dest = "bedfile", default = "", 
              help = "[opt] BED file for region to be markered"),
  make_option(c("-t", "--title"), dest = "title", default = "",
              help = "[opt] text shown on title"),
  make_option(c("-w", "--width"), dest = "fig_width", default = 8,
              help = "[opt] width (in inch). Default: 8."),
  make_option(c("--height"), dest = "fig_height", default = 8,
              help = "[opt] height (in inch). Default: 8.")
)
#
parser <- OptionParser(usage = "cgmaptools lollipop [options] file",
     option_list=option_list, description = "      (aka mCLollipop) \
Description: Plot local mC level for multiple samples \
Contact:     Guo, Weilong; guoweilong@126.com\
Last Update: 2018-04-10 \
Example: \
    mCLollipop [-i input] -o gene.png \
-Input Format (-i)\
    Can be output by \"cgmaptools mergelist tomatrix\". Use STDIN if omitted.\
    The 1st line (header line) is required.\
    Example: \
       chr     pos     tag1    tag2    tag3\
       Chr1    111403  0.30    nan     0.80\
       Chr1    111406  0.66    0.40    0.60\
-Site File (-s)\
    >= 3 columns, the 1st line (header line) is required, using R color name or \"NaN\". \
    To show specific sites (such as DMS, SNV) at the bottom as triangles. \
    Example: \
        chr   pos       A_vs_B  B_vs_C  A_vs_C\
        chr1  13116801  NaN     NaN     darkgreen\
        chr1  13116899  NaN     red     NaN\
-Region File (-b)\
    the first 4 columns are required.\
    To show specific region (such as DMR, Repeats) at the bottom as blocks. \
    Example: \
        chr1  213941196  213942363  hyper-DMR \
        chr1  213942363  213943530  hypo-DMR \
    #   chr   left       right      region-description \
-annotation file (-a), refFlat Format:\
    To show the structure of genes/transcripts. One-line in annotation, one-track in figure. \
    Example: \
        GeneA   TransA  chr2  +	     1000      2000       1100    1950     3     1100,1500,1700,  1200,1580,1950,\
    #   GeneID  TrandID ChrID Strand TransLeft TransRight CDSLeft CDSRight nExon ExonLefts        ExonRights\
    "
)
#
arguments <- parse_args(parser)
opt <- arguments$options
infile <- arguments$infile
#
if(infile == "") { # default, STDIN
  infile = file("stdin")
  #print_help(parser)
  #stop(sprintf("input file is not specified"))
}  else {
  if(file.access(infile) == -1) {
    stop(sprintf("Specified file ( %s ) does not exist", infile))
  }
}
#
outfile = arguments$outfile
annofile = arguments$annofile
figure.format = arguments$format
#
if( outfile == "") {
  outfile = paste( infile, figure.format, sep=".")
}
#
# Read the input file
T = read.table(infile, header=TRUE)
N_Sample = ncol(T)-2
#
chr=T[1,1]
if(is.na(chr)){
  chr=arguments$chr
}
mC = T[T[,1]==chr,c(-1,-2)]
pos = t(T[,2])
#
if( arguments$left=="" || arguments$right=="" ) {
  border = (max(pos)-min(pos))/8
  LeftMost = min(pos)#-border
  RightMost = max(pos)#+border
} else {
  LeftMost = as.integer(arguments$left)
  RightMost = as.integer(arguments$right)
}
#
if(LeftMost==RightMost){
  LeftMost = LeftMost-5000
  RightMost = RightMost+5000
}
#
yscale=1.1
Xscale = ceiling(sqrt(100000.0/(RightMost-LeftMost)))/10.0
Xscale
# ------------------------
# Draw the gene annotation using refFlat
if (annofile != ""){
  R = read.table(annofile)
  N_anno = nrow(R)
} else {
  N_anno = 0
}
# 
anno_height = 0.5
if (arguments$right=="") {
  figure.width = 8
} else {
  figure.width = arguments$fig_width
}
#
# ------------------------
# Read the site information
sitefile = arguments$sitefile
if (sitefile != "") {
  if ( file.info(sitefile)$size==0 ) {
    N_site = 0
    Ncol_site = 0
  } else {
    S = read.table(sitefile, header = TRUE)
    N_site = nrow(S)
    Ncol_site = ncol(S)-2
  }
} else {
  N_site = 0
  Ncol_site = 0
}
site_height = 0.5
#
# ------------------------
# Read the BED information
bedfile = arguments$bedfile
if (bedfile != "") {
  if ( file.info(bedfile)$size==0 ) {
    N_bed = 0
    Ncol_bed = 0
  } else {
    BED = read.table(bedfile, header = FALSE)
    N_bed = nrow(BED)
    Ncol_bed = ncol(BED)
  }
} else {
  N_bed = 0
  Ncol_bed = 0
}
bed_height = 0.25
# ========================
# Draw the figures
figure.width = as.integer(arguments$fig_width)
figure.height = as.integer(arguments$fig_height)
if (figure.format == "pdf") {
  pdf(outfile, width = figure.width, height = figure.height )
} else if (figure.format == "eps") {
  postscript(outfile, width = figure.width, height = figure.height) 
} else { # png
  png(outfile, width = figure.width*72, height = figure.height*72 )
}
# -------------------------
par( mai = c(1, 2, 1, 0.5) )
# Set the xlim and ylim
plot(NULL, NULL, 
     xlim = c(LeftMost, RightMost), 
     ylim = c(1 - Ncol_site*site_height - N_anno*anno_height - N_bed*bed_height, 
              N_Sample+2 )*yscale, 
     xlab = "", ylab = "",
     yaxt = "n", frame.plot = FALSE,
     main = arguments$title,
     cex.axis = 1.5, cex = 1.5, cex.main = 1.8 )
# -------------------------
#ColPanel = c( "#2E2EB2", 
#              "#0033CC", "#009999", "#00CC00", "#99FF00", "#FFFF00", 
#              "#FFCC00", "#FFCC00", "#FF9900", "#FF6600", "#FF0000")
ColPanel = rev(rainbow(13, end=0.7))[c(-6,-8)]
#ColPanel = topo.colors(11)
# -------------------------
# Draw the DNA strands
abline( h = c(1:N_Sample)*yscale, lwd = 2 )
# library("gplots")
for ( k in c(1:N_Sample) ) {
  x = pos[!is.na(mC[,k])]
  # Draw points for sites with enough coverages
  #points( x, rep(k, length(x))*yscale, pch = 20, col = rainbow(11, end = 0.7)[11], cex = 1.5 )
  # Draw the lines
  for( i in c(1:length(pos))[!is.na(mC[,k])] ){
    # lines( c(pos[i], pos[i]), c(0, mC[i,k])+k*yscale , lwd=5)
    lines( c(pos[i], pos[i]), c(0, mC[i,k])+k*yscale , lwd = 3*Xscale, 
           col = ColPanel[as.integer(mC[i,k]*10)+1] )
    points( pos[i], mC[i,k]+k*yscale , cex = 2*Xscale, pch = 20,
           col = ColPanel[as.integer(mC[i,k]*10)+1] )
  }
}
# -------------------------
# Draw the color legends
text = rep("",11)
text[seq(1,11,2)] = rev(as.character(c(0:10)/10)[seq(1,11,2)])
legend("topright", fill = rev(ColPanel), text, cex = 0.8,
       y.intersp=0.4, border = rev(ColPanel), bty = 'n')
# -------------------------
# Annotate samples
SampleName = colnames(T)[c(-1, -2)]
SampleNamePos=LeftMost-(RightMost-LeftMost)/20
text( rep(SampleNamePos, N_Sample), (c(1:N_Sample)+0.4)*yscale, SampleName, 
      col="black", xpd=TRUE, adj = c(1, 0), cex = 1.2 )
# -------------------------
# Draw the ruler bar
for ( k in c(1:N_Sample) ) {
  RulerWidth = (RightMost-LeftMost)/100
  RulerPos = LeftMost-RulerWidth*4
  lines(rep(RulerPos,2), rep(k*yscale, 2)+c(0, 1), lwd=2)
  lines(rep(RulerPos,2)+c(0,RulerWidth), rep(k*yscale+1, 2), lwd = 2)
  text(rep(RulerPos+RulerWidth*2, 2), rep(k*yscale, 2)+c(0.05, 0.85), c("0", "1"), cex = 0.8,
      adj=c(0.5,0))
}
# ==========================
# Text information

# Text ex:  "chr2: 90,936-91,653"
text( mean( c(LeftMost, RightMost) ), (N_Sample+1.75)*yscale, cex = 1.5,
      paste(chr, " : ", 
            format(LeftMost, big.mark = ","), " - ", 
            format(RightMost, big.mark = ","), sep = "") )
#
# ==========================
# Draw the gene annotation
y_anno_start = 1.3- site_height*Ncol_site
if(N_anno > 0){
  for (i in c(1:N_anno) ) {
    Left = as.integer(strsplit(as.character(R[i,10]), split = ",")[[1]])
    Right = as.integer(strsplit(as.character(R[i,11]), split = ",")[[1]])
    N_exon = length(Left)
    y_center = y_anno_start-anno_height*i
    y_width = 0.04
    GeneCol = "black"
    rect( Left, rep(y_center-y_width, N_exon), 
          Right, rep(y_center+y_width, N_exon), col = GeneCol )
    # Draw the CDS Regions
    CDSLeft = as.integer(R[i,7])
    CDSRight = as.integer(R[i,8])
    for (j in c(1:N_exon) ) {
      # 3' end partially overlap with CDS
      if ( CDSLeft>Left[j] && CDSLeft<Right[j] && CDSRight>Right[j] ) {
        rect( CDSLeft, y_center-y_width*2, Right[j], y_center+y_width*2, col = GeneCol )
      }
      # Whole exon in CDS
      if ( Left[j]>CDSLeft && Right[j]<CDSRight ) {
        rect( Left[j], y_center-y_width*2, Right[j], y_center+y_width*2, col = GeneCol ) 
      }
      # 5' end partially overlap with CDS
      if ( CDSRight>Left[j] && CDSRight<Right[j] && CDSLeft<Left[j] ) {
        rect( Left[j], y_center-y_width*2, CDSRight, y_center+y_width*2, col = GeneCol )
      }
      # CDS all inside the exon
      if ( Left[j]<CDSLeft && CDSRight<Right[j] ) {
        rect( CDSLeft, y_center-y_width*2, Right, y_center+y_width*2, col = GeneCol )         
      }
    }
    rect( Left, rep(y_center-y_width, N_exon), 
          Right, rep(y_center+y_width, N_exon), col = GeneCol )
    #
    lines( c(min(Left), max(Right)), c(y_center, y_center), col = GeneCol  )
    strand = R[i,4]
    StrandSign=(c(1,-1)[c("+","-")==strand])
    lines( rep(c(min(Left), max(Right))[c("+","-")==strand],2), 
           c(y_center, y_center+y_width*3), col = GeneCol  )
    x0 = c(min(Left), max(Right))[c("+","-")==strand]
    x1 = x0 + (RightMost-LeftMost)/15*StrandSign
    arrows( x0,   y_center+y_width*3, x1,  y_center+y_width*3, length = 0.08, col = GeneCol  )
    # Gene Name
    #text( x0-(RightMost-LeftMost)/10*StrandSign, y_center,  R[i,2], col = GeneCol)
    x0[x0<LeftMost] = LeftMost
    x0[x0>RightMost] = RightMost
    text( x0-(RightMost-LeftMost)/50*StrandSign, y_center,  R[i,2], col = GeneCol,
          xpd = TRUE, adj = c(0.5+0.5*StrandSign, 0), cex = 0.8)
  }
}
#
# ==========================
# Draw the marker sites
if(N_site > 0) {
  Ncol_site = ncol(S)-2
  if (Ncol_site>0) {
    colname_site = colnames(S)[c(-1,-2)]
    type_name_pos = SampleNamePos
    for (i in c(1:Ncol_site)) {
      y_pos = 1.3 - site_height*i
      x_list = (S[,2+i]!="NaN")
      points(S[x_list,2], rep(y_pos, sum(x_list)), pch = 17, col = as.character(S[x_list,2+i]), cex = 1.5)
      text(type_name_pos, y_pos, colname_site[i], col = "grey", cex = 0.8)
    }
#    points(S[,2], rep(0.8, N_site), pch = 17, col = as.character(S[,3]), cex = 1.5)
  }    
}
#
# ==========================
# Draw the marker region
y_bed_start = 1.3 - site_height*Ncol_site - anno_height*N_anno
if(N_bed > 0) {
  if (Ncol_bed >= 4) {
    for (i in c(1:N_bed)) {
      #if (chr == BED[i,1]) {
        y_pos = y_bed_start-bed_height*(i+0.5)
        x_list = BED[i,c(2,3)]
        bedCol="grey"
        if(Ncol_bed >=5) {
        	bedCol = as.character(BED[i,5])
        }
        rect( x_list[1], y_pos+bed_height*0.4, 
              x_list[2], y_pos-bed_height*0.4, 
              border = bedCol, col = bedCol)
        text( x_list[1]-(RightMost-LeftMost)/10, y_pos,
              BED[i,4], col = "black", cex = 0.8)
      #}
    }
    #    points(S[,2], rep(0.8, N_site), pch = 17, col = as.character(S[,3]), cex = 1.5)
  }    
}
#
invisible(dev.off())
#
