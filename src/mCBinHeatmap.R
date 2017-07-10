#!/usr/bin/env Rscript

# cgmaptools - mCBinHeatmap.R
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
if(!is.installed("gplots")){
  install.packages("gplots")
}
## libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  #make_option(c("-r", "--methylation_right"), dest = "methylation.right", default = TRUE,
  #            help="[opt] show average methylation level on the right panel."),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file name. [default: mCBinHeatmap.SysDate.pdf]"),
  make_option(c("-c", "--cluster"), dest = "cluster", default = FALSE, action="store_true",
              help = "[opt] cluster samples by methylation in regions. [default: FALSE]"),
  make_option(c("-l","--colorLow"), dest = "low.color", default = "cyan3",
              help = "[opt] color used for the lowest methylation value. [default: cyan3]"),
  make_option(c("-m","--colorMid"), dest = "mid.color", default = "",
              help = "[opt] color used for the middle methylation value. [default: null]"),
  make_option(c("-b","--colorHigh"), dest = "high.color", default = "coral2",
              help = "[opt] color used for the highest methylation value. [default: coral2]"),
  make_option(c("-n","--colorNumber"), dest = "num.color", default = 10,
              help = "[opt] desired number of color elements in the panel. [default: 10]"),
  make_option(c("-W","--width"), dest = "figure.width", default = 7,
              help = "[opt] width of figure (inch). [default: 7]"),
  make_option(c("-H","--height"), dest = "figure.height", default = 5,
              help = "[opt] height of figure (inch). [default: 7]"),
  make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
              help = "[opt] format of output figure. Alternative: png. [default: pdf]"),
  make_option(c("-R","--resolution"), dest = "figure.resolution", default = 300,
              help = "[opt] Resolution in ppi. Only available for png format. [default: 300]")
)

parser <- OptionParser(usage = "cgmaptools heatmap [options]",
                       option_list=option_list, description = "      (aka mCBinHeatmap)\
Description: Plot methylation dynamics of target region for multiple samples [heatmap]\
Contact:     Zhu, Ping; pingzhu.work@gmail.com\
Last update: 2016-12-07\
Example: \
  mCBinHeatmap.R -i input -m white -o chr1.xxx-xxx.pdf \
  -Input File Format: \
  1st line is the header.\
  Each column contains methylation measurements of a sample. \
  Example: \
  Region  Sample1  Sample2 ...  \
  Region1 0.1      0.1     ...  \
  Region2 0.1      0.1     ...  \
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
outfile <- arguments$outfile
if(outfile == ""){ # default, "FragRegView.Date"
  outfile <- paste("mCBinHeatmap", Sys.Date(), sep=".")
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
low.color <- arguments$low.color
mid.color <- arguments$mid.color
if(mid.color == ""){
    mid.color <- NA
}
high.color <- arguments$high.color
num.color <- arguments$num.color

data <- read.table(file = infile, header = T, stringsAsFactors = F, check.names = FALSE)
row.names(data) <- data[,1]
data[,1] <- NULL

## figure device
if (figure.format == "png"){
  png(outfile, height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
  pdf(outfile, height = figure.height, width = figure.width)
}

## cluster specified
cluster <- arguments$cluster
if (cluster){
    ## figure margin 
    par(mar=c(8, 0, 0, 0), oma=c(0,4,2,3))
    layout(matrix(1:2, nrow=1), widths=c(1,5))
    correlation <- cor(data, method="spearman", use="pairwise.complete.obs")
    hc <- hclust(d = as.dist(1-correlation), method = "ward.D2")
    dendrogram <- as.dendrogram(hc)
    plot(dendrogram, horiz=T, leaflab="none", yaxt = "none", dLeaf = 0, yaxs="i")

    # reorder sample index
    data <- data[, order.dendrogram(dendrogram)]
} else {
    ## figure margin 
    par(mar=c(5, 6, 4.1, 2.1))
}

## empty figure panel
plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
     xlab="", ylab="", 
     xlim=c(1, nrow(data) * 1.3), ylim=c(1, ncol(data)+1), xaxs="i", yaxs="i")

## methylation to colors
if ( is.na(mid.color) ){
    color.rects <- colorpanel(low=low.color, high=high.color, n=num.color)
} else {
    color.rects <- colorpanel(low=low.color, mid=mid.color, high=high.color, n=num.color)
}
methylation.interval <- seq(0, 1, length.out=num.color+1)
row.height <- 0.8
rightPanel.extension <- 1.15
for(i in 1:ncol(data)){
  #y.i <- (ncol(data):1)[i]
  y.i <- i
  x.color <- color.rects[ findInterval(data[,i], methylation.interval, all.inside=T)]

  # heatmap body
  rect(xleft = 1:nrow(data), 
       ybottom = rep(y.i, nrow(data)) + 0.1, 
       xright = 1:nrow(data) + 1, 
       ytop = rep(y.i, nrow(data)) + 0.1 + row.height,
       col = x.color,
       border = x.color
       )
    
  # methylation level on the right panel, barplot
  met.mean <- mean(as.numeric(data[,i]), na.rm = T)
  rect(xleft = nrow(data)+1,
       ybottom = y.i + 0.1,
       xright = met.mean * nrow(data)*(rightPanel.extension - 1) + nrow(data) + 1,
       ytop = y.i + 0.1 + row.height,
       col = "#66CD00",
       border = "#66CD00"
       )
  text(x=mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1,
       y=y.i+0.1+row.height/2,
       paste(format(met.mean*100,digits = 3), "%", sep = ""),
       cex = 0.7, adj = c(0, 0.5)
       )
}

## mean methylation label on right panel
segments(x0 = nrow(data)+1, # horiz segments
         y0 = 1 + 0.1, 
         x1 = nrow(data)*rightPanel.extension + 1,
         y1 = 1 + 0.1
         )
segments(x0 = c(nrow(data)+1, nrow(data)*rightPanel.extension + 1), # vertical ticks
         y0 = 1 + 0.1,
         x1 = c(nrow(data)+1, nrow(data)*rightPanel.extension + 1),
         y1 = 1 + 0.1 - ncol(data) * 0.008
         )
text(x = c(nrow(data)+1, nrow(data)*rightPanel.extension + 1), # value label
     y = 1 + 0.1 - ncol(data) * 0.01,
     paste(c(0,100), "%", sep = ""),
     cex = 0.8, xpd=T, adj = c(0.5, 1)
     )
text(x=nrow(data)+ 1 + nrow(data)*(rightPanel.extension-1)/2,
     y = 1 + 0.1 - ncol(data) * 0.08,
     "Average\nDNA methylation",
     cex=0.8, xpd=T, adj = c(0.5, 1)
     )


## sample label
sampleLabel.col <- "black"
if (cluster){
    text(x=1 + 0.1, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(0,0.5), col=sampleLabel.col)
} else {
    text(x=1, y=1:ncol(data)+row.height/2 + 0.1, names(data), cex=0.8, xpd=T, adj=c(1,0.5))
}


## methylation value legend
colorLegend.x.min <- nrow(data)*0.1
colorLegend.x.max <- nrow(data)*0.3
colorLegend.x <- seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = num.color+1)
colorLegend.y <- 0.8
rect(xleft = colorLegend.x[-length(colorLegend.x)],
     xright = colorLegend.x[-1],
     ybottom = colorLegend.y, 
     ytop = colorLegend.y - ncol(data) * 0.05,
     col = color.rects,
     border = color.rects,
     xpd = T
)
# ticks
segments(x0=seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = 3),
         x1 = seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = 3),
         y0 = rep(colorLegend.y - ncol(data) * 0.05, 3),
         y1 = rep(colorLegend.y - ncol(data) * 0.058, 3),
         xpd = T
)
# value label
text(x=seq(from=colorLegend.x.min, to=colorLegend.x.max, length.out = 3),
     y=rep(colorLegend.y - ncol(data) * 0.057, 3),
     paste(seq(0,100,length.out = 3), "%", sep = ""),
     adj = c(0.5, 1),
     cex=0.8,
     xpd = T
)
# text label
#text(x = colorLegend.x.min + (colorLegend.x.max - colorLegend.x.min)/ 2,
#     y = colorLegend.y + ncol(data) * 0.015,
#     "DNA methylation",
#     adj = c(0.5, 0),
#     cex = 0.8,
#     xpd = T
#)
                
invisible(dev.off())



