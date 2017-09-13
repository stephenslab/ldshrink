read_SNPinfo_gds <- function(gds,alleles=F,MAF=F){

    tdf <- tibble::data_frame(SNP=seqGetData(gds,var.name="annotation/id"),
                              snp_id=seqGetData(gds,var.name="variant.id"),
                              chr=seqGetData(gds,var.name="chromosome"),
                              pos=seqGetData(gds,var.name="position"))
    if(alleles){
      tdf <- mutate(tdf,allele=seqGetData(gds,var.name="allele"))
    }
    if(MAF){
      tdf <- mutate(tdf,MAF=seqAlleleFreq(gds))
    }
    return(tdf)
}



#read_

read_si_ql <- function(chunk_evdf){
  library(RcppEigenH5)
  library(tidyr)
  library(dplyr)
  retl <- list()
  retl[["LDinfo"]] <- read_df_h5(chunk_evdf,"LDinfo") %>% nest(-region_id)
  region_id <-  retl[["LDinfo"]][["region_id"]]
  retl[["Dl"]] <- lapply(region_id,read_dvec,h5file=chunk_evdf,dataname="D")
  names(retl[["Dl"]]) <- region_id
  retl[["Ql"]] <- lapply(region_id,read_2d_mat_h5,h5file=chunk_evdf,dataname="Q")
  names(retl[["Ql"]]) <- region_id
  return(retl)
}

read_SNPinfo_ldsc <- function(gds){
  return(tibble::data_frame(CHR=seqGetData(gds,var.name="chromosome"),
                            SNP=seqGetData(gds,var.name="annotation/id"),
                            BP=seqGetData(gds,var.name="position"),
                            CM=seqGetData(gds,var.name="annotation/info/map"),
                            MAF=seqAlleleFreq(gds)))
}






is_SNV <- function(allelevec){
  stopifnot(is.character(allelevec))
  return(sapply(strsplit(allelevec,split=",",fixed=T,useBytes = T),function(x){max(sapply(x,nchar))})==1)
}


import_panel_data <- function(temp_gds=NULL,
                              map_df=NULL,
                              output_file=NULL,
                              overwrite=F,
                              parallel=1,
                              ld_break_file=NULL){
  library(SeqArray)
  library(dplyr)
  library(purrr)
  library(readr)
  stopifnot(!is.null(output_file),
            !is.null(map_df),
            all(file.exists(temp_gds)),
            !is.null(ld_break_file))

    stopifnot(file.exists(ld_break_file))

  if(file.exists(output_file)){
    stopifnot(overwrite) 
    file.remove(output_file)
  }
  seqMerge(temp_gds,output_file,storage.option="LZ4_RA.fast")
  
  gds <- seqOpen(output_file,readonly = F)
  leg_df <- read_SNPinfo_gds(gds) %>% left_join(map_df)
  stopifnot(nrow(leg_df)==length(seqGetData(gds,"variant.id")))


  gdsfmt::add.gdsn(index.gdsn(gds,"annotation/info"),"map",leg_df$map,replace=T,compress="LZ4_RA.fast")
  cat("Adding Chunk Delimiters\n")
  seqClose(gds)
  add_chunk_gds(output_file,ld_break_file)
}


add_chunk_gds <- function(gds_file,region_bed_file){
  library(readr)
  library(dplyr)
  library(SeqArray)
  stopifnot(file.exists(region_bed_file),file.exists(gds_file))
  gds <- seqOpen(gds_file,readonly = F)
  region_bed <- read_delim(region_bed_file,delim="\t",trim_ws = T) %>% mutate(range_id=as.character(1:n()),chr=gsub("chr","",chr))
#  region_range <- GenomicRanges::GRanges(seqnames=region_bed$chr,ranges=IRanges::IRanges(start=region_bed$start,end = region_bed$stop))
  snp_df <- read_SNPinfo_gds(gds) %>% mutate(snp_id=as.character(snp_id))
  match_df <- match_SNP(region_bed,snp_df) %>% arrange(as.integer(snp_id))
  stopifnot(all(snp_df$snp_id==match_df$snp_id))
  gdsfmt::add.gdsn(index.gdsn(gds,"annotation/info"),"LD_chunk",as.integer(match_df$range_id),replace=T,compress="LZ4_RA.fast")
  seqClose(gds)
  
}

calc_LD_gds <- function(gds,m=85,Ne=11490.672741,cutoff=1e-3){
  map_dat <- seqGetData(gds,var.name="annotation/info/map")
  H <- 2.0-seqGetData(gds,var.name="$dosage")
  stopifnot(max(H)==1)
  return(calcLD(hmata = H,mapa = map_dat,m = m,Ne = Ne,cutoff = cutoff))
}

is_haplo <- function(gds,checkSNPs=F){
  samples <- seqGetData(gds,"sample.id")
  sample_trim <- substr(samples,1,nchar(samples)-2)
  if(all(names(table(table(sample_trim)))=="2")){
    return(T)
  }
  return(F)
}

    

chunkwise_LDshrink <- function(gds_file,region_id=1,outfile=NULL,m=85,Ne=11490.672741,cutoff=1e-3,evd=T){
  library(SeqArray)
  library(RcppEigenH5)
  library(dplyr)
  stopifnot(file.exists(gds_file),!is.null(outfile))
  gds <- seqOpen(gds_file,readonly = T)
  LD_chunks <- seqGetData(gds,var.name="annotation/info/LD_chunk")
  stopifnot(sum(LD_chunks%in%region_id)>0)
  seqSetFilter(gds,variant.sel=LD_chunks %in% region_id)

  si <- read_SNPinfo_gds(gds) %>% mutate(region_id=region_id)
  R <- calc_LD_gds(gds,m=m,Ne=Ne,cutoff=cutoff)
  # map_dat <- seqGetData(gds,var.name="annotation/info/map")
  # stopifnot(max(H)==1)
  # 
  # R <- calcLD(hmata = t(H),mapa = map_dat,m = m,Ne = Ne,cutoff = cutoff)
  
  RcppEigenH5::write_df_h5(df=si,groupname = "LDinfo",outfile = outfile,deflate_level = 4)
  tls <- h5ls(outfile)
  stopifnot(any("/LDinfo/SNP" %in% tls))
  RcppEigenH5::write_mat_h5(outfile,groupname = "LD",dataname = "R",data = R,deflate_level = 4,doTranspose = F)
  if(evd){
    eigenR <- eigen(R)
    stopifnot(min(eigenR$values)>0)
    write_mat_h5(outfile,groupname="EVD",dataname = "Q",data = eigenR$vectors,deflate_level = 0L)
    write_dvec_h5(outfile,groupname="EVD",dataname = "D",data = eigenR$values,deflate_level = 2L)
  }
  return(dim(R))
}


chunkwise_LDshrink_ldsc <- function(gds_file,chrom=19,out_dir=NULL,m=85,Ne=11490.672741,cutoff=1e-3){
  library(SeqArray)
  library(RcppEigenH5)
  library(dplyr)


  stopifnot(!is.null(gds_file),
            file.exists(gds_file),
            !is.null(out_dir),dir.exists(out_dir))
  outfile <- file.path(out_dir,paste0(chrom,".l2.ldscore.gz"))
  soutfile <- file.path(out_dir,paste0(chrom,".l2.M_5_50"))
  gds <- seqOpen(gds_file,readonly = T)
  # chroms <- seqGetData(gds,"chromosome")
  seqSetFilterChrom(gds,include=as.character(chrom))
  LD_chunks <- seqGetData(gds,var.name="annotation/info/LD_chunk")
  region_ids <- unique(LD_chunks)
  mfilt <- seqGetFilter(gds)
  resl <- list()
  library(progress)
  pb <- progress::progress_bar$new(total=length(region_ids))
  for(i in 1:length(region_ids)){
    region_id <- region_ids[i]
    tr <- LD_chunks %in% region_id
    seqSetFilter(gds,variant.sel=tr,action="push+intersect")
    si <- read_SNPinfo_ldsc(gds) 
    R <- calc_LD_gds(gds,m = m,Ne = Ne,cutoff = cutoff)
    si <- mutate(si,L2=colSums(R^2)-1)
    resl[[i]] <- si
    seqSetFilter(gds,action="pop")
    tfilt <- seqGetFilter(gds)
    stopifnot(all.equal(mfilt,tfilt))
    pb$tick()
  }
  resdf <- bind_rows(resl)
  write_delim(resdf,path=outfile,delim = "\t")
  nc <- filter(resdf,MAF>0.05) %>% nrow
  write(x = nc,file = soutfile)
  return(dim(R))
}


# 
# haplo_2_geno<-function(x,snps_in_rows=T){
#   if(snps_in_rows){
#     return(apply(array(c(x),dim=c(2,ncol(x)/2,nrow(x))),c(2,3),sum))
#   }
#   else{
#     return(apply(array(c(t(x)),dim=c(2,nrow(x)/2,ncol(x))),c(2,3),sum))
#   }
# }
  


