#!/usr/bin/env Rscript

# cgmaptools - mCFragRegView.R
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


is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("RColorBrewer")){
  warning("Detect package \"RColorBrewer\" is not installed in your R enviroment.")
  warning("Trying to install the \"RColorBrewer\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("RColorBrewer")
}
if(!is.installed("optparse")){
  warning("Detect package \"optparse\" is not installed in your R enviroment.")
  warning("Trying to install the \"optparse\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("optparse")
}

## libraries
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-r", "--ratio"), dest = "range.ratio", default = 5,
              help="[opt] range ratio between target region and flanking region in plot. [default: 5]"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file name. [default: FragRegView.SysDate.pdf"),
  make_option(c("-W", "--width"), dest = "figure.width", default = 7,
              help = "[opt] width of figure (inch). [default: 7]"),
  make_option(c("-H", "--height"), dest = "figure.height", default = 7,
              help = "[opt] height of figure (inch). [default: 7]"),
  make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
              help = "[opt] format of output figure. Alternative: png. [default: pdf]"),
  make_option(c("-R","--resolution"), dest = "figure.resolution", default = 300,
              help = "[opt] Resolution in ppi. Only available for png format. [default: 300]")
)

parser <- OptionParser(usage = "cgmaptools fragreg [options]",
                       option_list=option_list, description = "      (aka mCFragRegView) \
Description: Plot methylation dynamics of target and flanking region for multiple samples \
Contact:     Zhu, Ping; pingzhu.work@gmail.com\
Last update: 2018-02-12\
Example: \
  FragRegView.R -i input -r 5 -o genebody.pdf \
-Input File Format: \
  1st line is the header.\
  Each row contains methylation measurements of a sample. \
  The user may need to use shell script to generate following format \
 based on the results of \"cgmaptools mfg\". \
Example: \
  Sample  Up1  Up2  ...  Region1  Region2 ...  Down1  Down2  ...\ 
  Sample1 0.1  0.1  ...  0.2      0.2     ...  0.3    0.3    ...\
  Sample2 0.1  0.1  ...  0.2      0.2     ...  0.3    0.3    ...\
"
)
## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
if(infile == ""){ # default, STDIN
  infile <- file("stdin")
} else { # user specified
  if( file.access(infile) == -1){ # file not exists
    print_help(parser)
  }
}
range.ratio <- arguments$range.ratio # range ratio (Region/Flanking)
outfile <- arguments$outfile
if(outfile == ""){ # default, "FragRegView.Date"
  outfile <- paste("FragRegView", Sys.Date(), sep=".")
} else { # user specified
  outfile <- gsub(".pdf$|.png$", "", outfile, perl=T)
}
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
figure.format <- arguments$figure.format
if(! figure.format %in% c("pdf", "png")){ # format not support
  print_help(parser)
} else {
  outfile <- paste(outfile, figure.format, sep = ".")
}
figure.resolution <- arguments$figure.resolution

## Input format
## keyword: Up
## keyword: Region
## keyword: Down
## Sample   Up1  Up2  Up3  ...  Region1  Region2  Region3  ...  Down1  Down2  Down3  ...
## Sample1  0.1  0.1  0.1  ...  0.2      0.2      0.2      ...  0.5    0.5    0.5    ...
## Sample2  0.1  0.1  0.1  ...  0.2      0.2      0.2      ...  0.5    0.5    0.5    ...

# ## test data
# up <- seq(from=0.2, to=0.1, length.out = 10)
# down <- seq(from=0.1, to=0.2, length.out = 10)
# region <- seq(from=0.05, to=0.3, length.out = 10)
# sample1 <- c(up, region, down)
# data <- as.data.frame(matrix(c("Sample1", sample1, 
#                                "Sampel2", sample1 * 2, 
#                                "Sample3", sample1 * 3), 
#                              nrow=3, byrow = T), stringsAsFactors = F)
# names(data) <- c("Sample", paste("Up", 1:10, sep=""), 
#                  paste("Region", 1:10, sep=""), 
#                  paste("Down", 1:10, sep = ""))
# row.names(data) <- data$Sample
# data$Sample <- NULL
# write.table(x=data, file="test.txt", sep="\t", quote = F, row.names = T, col.names = T)

data <- read.table(file = infile, header = T, stringsAsFactors = F)
row.names(data) <- data[,1]
data[,1] <- NULL

## search for "Region" windows
## by names of the data 
region.columns <- grep("R", names(data))
if (is.na(region.columns[1])){ ## ERROR: can not determine windows in Region
}

## range ratio (Up : Region : Down): 1:5:1
#range.ratio <- 5
# Up
x.pos.up <- seq(from=1, length.out = region.columns[1]-1)
# Region
x.pos.region <- seq(from=region.columns[1], by=range.ratio, length.out = length(region.columns))
# Down
x.pos.down <- seq(from=x.pos.region[length(x.pos.region)]+1, length.out = region.columns[1]-1)
x.pos <- c(x.pos.up, x.pos.region, x.pos.down)

## figure
if (figure.format == "png"){
  png(outfile, height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
  pdf(outfile, height = figure.height, width = figure.width)
}

## figure margin
par(oma=c(2,2,2,2))

## empty figure panel
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="DNA methylation level (%)", 
     xlim=c(0, max(x.pos)+1), ylim=c(0, 1))
## x axis
# Up
rect(xleft = 1, ybottom = -0.035, xright = x.pos[region.columns[1]], ytop = -0.015, col="grey", border = "grey", xpd=T)
# Down
rect(xleft = x.pos[region.columns[length(region.columns)]], ybottom = -0.035, xright = max(x.pos), ytop = -0.015, col="grey", border = "grey", xpd=T)
# Region
rect(xleft = x.pos[region.columns[1]], ybottom = -0.04, xright = x.pos[region.columns[length(region.columns)]], ytop = -0.01, col="blue", border = "blue", xpd=T)
# ticks
x0.region <- seq(from=x.pos.region[1], to=x.pos.region[length(x.pos.region)], length.out = 6)
segments(x0=x0.region, x1 = x0.region, y0=-0.04, y1=-0.045, xpd=T, lwd=2)
text(x=x0.region, y=-0.06, paste(seq(0,100, length.out = 6), "%", sep=""), adj=c(0.5, 1), xpd=T, lwd=0.8)
## y axis
axis(side=2, at=seq(0,1, length.out = 6), labels = seq(0, 100, length.out = 6), las=1)
## colors
col.panel <- colorRampPalette(brewer.pal(n=8, name = "Dark2"),bias=0.1,space="rgb")
if ( nrow(data) > 8 ){
    col.sample <- col.panel(nrow(data))
} else {
    col.sample <- brewer.pal(n=nrow(data), name = "Dark2")
}

## Methylation ~ samples
for( i in 1:nrow(data)){
  lines(x=x.pos, y=data[i,], col=col.sample[i], lwd=2)
}

## legends
legend(x=0, y = 1, fill = c("blue", "grey"), 
       border = c("blue", "grey"),
       legend = c("Target region", "Flanking region"), 
       bty="n", cex=0.8)
legend(x=0, y=1.1, col=col.sample,
       legend = row.names(data),
       lwd=2, bty="n", cex=0.6, ncol = 3, xpd=T)

invisible(dev.off())



