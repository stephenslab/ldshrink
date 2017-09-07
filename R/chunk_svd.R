#' match_SNP
#' Given a dataframe of ranges, match another dataframe of SNPs to that range
#' @
#' @param break_df A dataframe having the following columns:
#' chr (character or integer)
#' start (integer)
#' stop (integer)
#' range_id (integer or character)
#' @param snp_df A dataframe having the following columns
#' chr (character or integer)
#' pos (integer)
#' snp_id (character or integer)
#' @param match_at_start indicate whether you want to boundary SNPs to match at start or stop
#' @return 
#' A dataframe joining the two dataframes when start<=pos<=stop
match_SNP <- function(break_df,snp_df,match_at_start=T){
  library(dplyr)
  library(purrr)
  
  # Whether or not break_df and snp_df use character or integer for chr, we'll use character
  break_cols <- colnames(break_df)
  snp_cols <- colnames(snp_df)
  stopifnot(all(c("chr","start","stop","range_id") %in% break_cols))
  stopifnot(all(c("chr","pos","snp_id") %in% snp_cols))
  
  snp_chr <- is.character(snp_df$chr)
  break_chr <- is.character(break_df$chr)
  stopifnot((!snp_chr & !break_chr) | (snp_chr & break_chr))
  
  stopifnot(nrow(anti_join(snp_df, break_df, by = "chr")) == 0)
  if (match_at_start) {
    break_df <- mutate(break_df,tstop = stop - 1)
    break_df <- mutate(break_df,tstart = start)
  }else{
    break_df <- mutate(break_df,tstart = start + 1)
    break_df <- mutate(break_df,tstop = stop)
  }
  #Start by just subsetting break_df by chr in snp_df
  ch_break_df <- semi_join(break_df, snp_df, by = "chr")
  ch_breakl <- split(ch_break_df,ch_break_df$chr)
  snp_breakl <- split(snp_df,snp_df$chr)
  
  stopifnot(length(ch_breakl) == length(snp_breakl))
  stopifnot(all(names(ch_breakl) == names(snp_breakl)))
  
  match_df <- bind_rows(map2(ch_breakl,snp_breakl,function(ch_break_df,snp_df){
    ld_ranges <- IRanges::IRanges(ch_break_df$tstart,ch_break_df$tstop,names = ch_break_df$range_id)
    pos_ranges <- IRanges::IRanges(snp_df$pos,snp_df$pos,names = snp_df$snp_id)
    ol <-   IRanges::findOverlapPairs(pos_ranges,ld_ranges)
    match_df <- data_frame(snp_id = ol@first@NAMES,range_id = ol@second@NAMES) %>% 
      inner_join(snp_df,by = "snp_id") %>% 
      inner_join(ch_break_df,by = c("range_id","chr"))
    return(match_df)
  }))
  n_snp_matches <- group_by(match_df,chr) %>% summarise(nmatches = n())
  n_snp_ct <- group_by(snp_df,chr) %>% summarise(nsnp=n())
  snp_match_sum <- inner_join(n_snp_matches,n_snp_ct,by="chr")
  n_unmapped <- filter(snp_match_sum,nmatches<nsnp) %>% summarise(n_unmapped=sum(nmatches)-sum(nsnp)) %>% pull(n_unmapped)
  n_dups <- filter(snp_match_sum,nmatches>nsnp) %>% summarise(n_dups=sum(nmatches)-sum(nsnp)) %>% pull(n_dups)
  if (n_unmapped > 0){
    warning(paste0(n_unmapped, " SNPs did not map to any region"))
  }
  if (n_dups > 0) {
    warning(paste0(n_dups, " SNPs mapped to multiple regions"))
  }
  match_df <- group_by(match_df,chr) %>% summarise(p=n_distinct(snp_id)) %>% inner_join(match_df)
  
  return(match_df %>% select(-tstart,-tstop))
}


#' read_SNPinfo
#' reads SNPinfo from an HDF5 file
#' @param snpfile
#' @param chr_to_char boolean specifying whether to convert chrom to character type (by prefixing "chr"),default is TRUE
#' @param extra_cols character vector specifying other columns you'd like returned 
#' @return dataframe with (at least) `chr`,`pos` and optionally `extra_cols`
read_SNPinfo <- function(snpfile,chr_to_char=T, extra_cols = NULL, id_col=NULL){
  library(RcppEigenH5)
  library(dplyr)
  # snpfile_dsets <- RcppEigenH5::h5ls(snpfile)
  
  pos <- read_ivec(snpfile,"/","pos")
  chr <- read_ivec(snpfile,"/","chr")
  snp_df <- data_frame(pos = pos, chr = chr) %>% mutate(snp_id = paste0("SNP: ",1:n()))
  if (chr_to_char) {
    snp_df <- mutate(snp_df,chr = paste0("chr",chr))
  }
  return(snp_df)
}



#' chunk_R_matfile
#' Given a path to a HDF5 file with a sparse LD matrix and 
#' a path to a bed file giving breakpoints, chunk the LD matrix (keeping track of where the chunks are)
#' @param sparse_Rfile path to HDF5 file. Must have group "R",data "pos"(integer), and data "chr"(integer)
#' @param break_df path to bed file. must have column chr, 
chunk_R_evd_matfile <- function(match_df){
  library(tidyr)
  library(purrr)
  sparse_Rfile <- unique(match_df$sparse_Rfile)
  p <- unique(match_df$p)
  
  stopifnot(file.exists(sparse_Rfile),length(sparse_Rfile)==1,
            length(p)==1,length(unique(match_df$chr))==1)
  R <- read_ccs_h5(sparse_Rfile,groupname = "R",dims=c(p,p))
  match_l <- split(match_df,match_df$range_id)
  match_R <- map(match_l,function(x,sparse_Rfile,R){
    pos_id <- as.integer(gsub("[^0-9]+","",x$snp_id))
    stopifnot(all(!is.na(pos_id)))
    nx <- nest(x,-range_id,-start,-stop,-chr,.key="snpinfo")
    return(mutate(nx,R=list(R[pos_id,pos_id])))
  },R = R) %>% bind_rows()
  match_R <- rowwise(match_R) %>% mutate(evdR = list(eigen(R))) %>% ungroup()
  Ql <- transpose(match_R$evdR)
  match_R <- select(match_R,-evdR) %>% mutate(Q=Ql$vectors,d=Ql$values)
  return(match_R)
}



brute_write_SVD <- function(svd_chunk,outfile){
  library(dplyr)
  library(RcppEigenH5)
  library(purrr)
    svd_chunk <- mutate(svd_chunk,
                        chr=as.integer(gsub("chr","",chr)),
                        range_id=as.integer(gsub("range: ","",range_id)))
    chr <- unique(svd_chunk$chr)
    stopifnot(length(chr)==1)
    write_ivec_h5(h5file = outfile,groupname = "LDchunks",dataname="start",data = svd_chunk$start,deflate_level = 4)
    write_ivec_h5(h5file = outfile,groupname = "LDchunks",dataname="stop",data = svd_chunk$stop,deflate_level = 4)
    write_ivec_h5(h5file = outfile,groupname = "LDchunks",dataname="range_id",data = svd_chunk$range_id,deflate_level = 4)
    write_data_string_attr_h5(h5filename = outfile,h5_groupname = "LDchunks",h5_dataname = "range_id",h5_attr_name = "prefix",h5_attr_value = "range: ")
    write_group_int_attr_h5(h5filename = outfile,h5_groupname = "/",h5_attr_name = "chr",h5_attr_value = chr)
    iwalk(svd_chunk$snpinfo,function(x,y,outfile,range_id){
      snp_id <- as.integer(gsub("SNP: ","",x$snp_id))
      pos <- x$pos
      write_ivec_h5(h5file = outfile, groupname = range_id[y],dataname = "snp_id",data = snp_id,deflate_level = 4)
      write_ivec_h5(h5file = outfile, groupname = range_id[y],dataname = "pos",data = pos,deflate_level = 4)
    },outfile=outfile,range_id=svd_chunk$range_id)
    iwalk(svd_chunk$Q,function(x,y,outfile,range_id){
      write_mat_h5(h5file = outfile, groupname = range_id[y],dataname = "Q",data = x,deflate_level = 4)
    },outfile=outfile,range_id=svd_chunk$range_id)
    iwalk(svd_chunk$d,function(x,y,outfile,range_id){
      write_dvec_h5(h5file = outfile, groupname = range_id[y],dataname = "d",data = x,deflate_level = 4)
    },outfile=outfile,range_id=svd_chunk$range_id)
    return(T)
}

chunk_all_LD <- function(sparse_Rfiles,breakfiles,outfiles=NULL){
  library(dplyr)
  library(purrr)
  if(is.null(outfiles)){
    sparse_dir <- dirname(sparse_Rfiles)
    sparse_files <- basename(sparse_Rfiles)
    outfiles <- file.path(sparse_dir,paste0("EVD_",sparse_files))
  }
  stopifnot(length(outfiles) == length(sparse_Rfiles))
  break_df <- map_df(breakfiles,readr::read_delim,delim = "\t", trim_ws = T) %>% mutate(range_id=paste0("range: ", 1:n()))
  snp_df <- map_df(sparse_Rfiles,function(x){read_SNPinfo(x) %>% mutate(sparse_Rfile=x)})
  match_df <- match_SNP(break_df,snp_df)
  
  match_df_l <- split(match_df,match_df$chr)
  names(outfiles) <- gsub(".+_1kg_([0-9]+).mat","chr\\1",outfiles)
  iwalk(match_df_l,function(x,y,outfiles){
    cat(y," ",outfiles[y],"\n")
    brute_write_SVD(chunk_R_evd_matfile(x),outfiles[y])
  },outfiles=outfiles)
}





