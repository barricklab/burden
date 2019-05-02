#!/usr/bin/env Rscript

require(dplyr)
require(ggplot2)
require(tidyr)
require(gridExtra)
require(cowplot)
require(readr)

require(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input file. You may separate values by commas to analyze multiple input files together", metavar="input.csv")
)

usage_string = paste(
  "summary-graph.R -i input1,input2,input3\n\n",
  sep = ""
) 

opt_parser = OptionParser(usage=usage_string, option_list=option_list)
opt = parse_args(opt_parser)

## Allow loading of multiple different results
input.file.names = strsplit(opt$input,",")

## Alternative that pastes in expected file name
##input.file.names = paste0(strsplit(opt$input, ","), ".rates.summary.tsv")


#debug
input.file.names = c("exp003.rates.summary.tsv","exp005.rates.summary.tsv")

###
# Special sample names for metadata
#  "blank" well will be used to create an average blank value
#  "ignore" well will be removed before analysis

readings = c("growth", "GFP") 

all.data = data.frame()
for (this.file.name in input.file.names) {
  this.data <- read_tsv(this.file.name)
  if("graph" %in% colnames(this.data)) {
    this.data = this.data %>% filter(graph==1)
  }
  all.data = rbind(all.data, this.data)
  
  #remove nongraphable ones
  
}

all.data$replicate=as.factor(all.data$replicate)

## Need to programmatically define the error bars from variable names...
ggplot(all.data, aes_(x=as.name(paste0(readings[1], ".rate")), y=as.name(paste0(readings[2], ".rate")), color=as.name("strain"), shape=as.name("replicate")))  +
  geom_errorbarh(aes(xmin=growth.rate-growth.rate.sd, xmax=growth.rate+growth.rate.sd)) +
  geom_errorbar(aes(ymin=GFP.rate-GFP.rate.sd, ymax=GFP.rate+GFP.rate.sd)) + 
  geom_point(size=5)   +
  scale_x_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) + 
  scale_y_continuous(limits = c(0, max(all.data$GFP.rate+all.data$GFP.rate.sd))) + 
  geom_abline(intercept=0, slope = 600/1.4)

ggsave(sub(".tsv", ".pdf", input.file.names[1]))

## fit line and show entire range

fit.data = all.data
if("fit" %in% colnames(fit.data)) {
  fit.data = this.data %>% filter(fit==1)
}

fit_fixed_zero = lm(GFP.rate~growth.rate + 0, fit.data)
slope_fixed_zero = coef(fit_fixed_zero)

ggplot(all.data, aes_(x=as.name(paste0(readings[1], ".rate")), y=as.name(paste0(readings[2], ".rate")), color=as.name("strain"), shape=as.name("replicate")))  +
  geom_errorbarh(aes(xmin=growth.rate-growth.rate.sd, xmax=growth.rate+growth.rate.sd)) +
  geom_errorbar(aes(ymin=GFP.rate-GFP.rate.sd, ymax=GFP.rate+GFP.rate.sd)) + 
  geom_point(size=5)  +
  scale_x_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) + 
  scale_y_continuous(limits = c(0, max(all.data$GFP.rate+all.data$GFP.rate.sd))) + 
  geom_abline(intercept=0, slope = slope_fixed_zero)

ggsave(sub(".tsv", ".fixed_zero.pdf", input.file.names[1]))


all.data$replicate=as.factor(all.data$replicate)
ggplot(all.data, aes_(x=as.name("strain"), y=as.name(paste0(readings[1], ".rate")), fill=as.name("replicate")))  +  geom_bar(size=3, stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=growth.rate-growth.rate.sd, ymax=growth.rate+growth.rate.sd), position=position_dodge()) + 
  scale_y_continuous(limits = c(0, max(all.data$growth.rate+all.data$growth.rate.sd))) + 

ggsave(sub(".tsv", ".growth.rates.pdf", input.file.names[1]))

