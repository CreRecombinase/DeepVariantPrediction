#!/usr/bin/env Rscript
library(optparse)
library(biomaRt)

option_list = list(
  make_option("--input", type="character", default=NULL, 
              help="input SNP list where rsID header is marker and chromosome information is in chr, referenece_allele and other_allele are also specified.", metavar="character"),
  make_option("--output", type="character", default=NULL, 
              help="Name of output file", metavar="character"),
  make_option('--group_size', type='double', default=1000)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data <- read.table(opt$input, sep = '\t', header = T)
rs_list <- as.character(data$marker)

snp_mart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
snp_attributes <- c("refsnp_id", "chr_name", "chrom_start", "allele_1")
snp_size <- length(rs_list)
ngroups <- ceiling(snp_size / group_size)
snp_locations <- c()
for(i in 1 : ngroups){
  start <- (i - 1) * group_size + 1
  end <- i * group_size
  end <- min(end, snp_size)
  snp_locations_i <- getBM(attributes=snp_attributes, filters="snp_filter", 
                      values=rs_list[start:end], mart=snp_mart)
  snp_locations <- rbind(snp_locations, snp_locations_i)
}
