#' Obtain info about variants in a gds object
#' 
#' 
#' @param gds an `SeqArray` `gds` object (obtained by using the `seqOpen` function in the package `SeqArray`)
#' @param alleles A boolean indicating whether alleles should be returned 
#' @param MAF A boolean indicating whether MAF should be computed and returned of not
#' @param region_id A boolean indicating whether to return the LD block of each SNP
read_SNPinfo_gds <- function(gds,alleles=F,MAF=F,region_id=F){

    tdf <- tibble::data_frame(SNP=SeqArray::seqGetData(gds,var.name="annotation/id"),
                              snp_id=SeqArray::seqGetData(gds,var.name="variant.id"),
                              chr=SeqArray::seqGetData(gds,var.name="chromosome"),
                              pos=SeqArray::seqGetData(gds,var.name="position"))
    if(alleles){
      tdf <- dplyr::mutate(tdf,allele=SeqArray::seqGetData(gds,var.name="allele"))
    }
    if(MAF){
      tdf <- dplyr::mutate(tdf,MAF=SeqArray::seqAlleleFreq(gds))
    }
    if(region_id){
      tdf <- dplyr::mutate(tdf,region_id=SeqArray::seqGetData(gds,"annotation/info/LD_chunk"))
    }
    return(tdf)
}


subset_gds <- function(gds,info_df,region_id=F){
  si_df <- dplyr::inner_join(info_df, read_SNPinfo_gds(gds,region_id=region_id)) %>% 
    dplyr::distinct(snp_id, .keep_all = T) %>% dplyr::arrange(snp_id)
  SeqArray::seqSetFilter(gds,variant.sel = si_df$snp_id)
  return(si_df)
}

mult_eigen_chunks <- function(chunk_evdf,uhmat){
  retl <- list()
  retl[["LDinfo"]] <- RcppEigenH5::read_df_h5(chunk_evdf,"LDinfo")
  region_id <-  unique(retl[["LDinfo"]][["region_id"]])

  retl[["D"]] <-  unlist(lapply(region_id,RcppEigenH5::read_dvec,h5file=chunk_evdf,dataname="D"))
  stopifnot(length(retl$LDinfo$SNP)==length(retl$D))
  stopifnot(!is.unsorted(retl$LDinfo$snp_id))
  retl$LDinfo$snp_id <- retl$LDinfo$snp_id-min(retl$LDinfo$snp_id)+1
  return(sparse_matmul(chunk_evdf,as.character(region_id),rep("Q",length(region_id)),uhmat,do_transpose = T))
}


read_SNPinfo_ldsc <- function(gds){
  return(tibble::data_frame(CHR=SeqArray::seqGetData(gds,var.name="chromosome"),
                            SNP=SeqArray::seqGetData(gds,var.name="annotation/id"),
                            BP=SeqArray::seqGetData(gds,var.name="position"),
                            CM=SeqArray::seqGetData(gds,var.name="annotation/info/map"),
                            MAF=SeqArray::seqAlleleFreq(gds)))
}






is_SNV <- function(allelevec){
  stopifnot(is.character(allelevec))
  return(sapply(strsplit(allelevec,split=",",fixed=T,useBytes = T),function(x){max(sapply(x,nchar))})==1)
}



seqIMPUTE2GDS <- function(hap.fn,leg.fn,sample.fn,map.fn,out.gdsfn,chrom,compress.geno="LZ4_RA.fast",compress.annotation="LZ4_RA.fast"){
  
  map_dat <- readr::read_delim(map.fn,delim=" ",col_names = c("ID","pos","map")) %>% dplyr::mutate(chrom=chrom)
  sample.id <- scan(sample.fn,what=character(),sep="\n")
  sample.id <- c(sapply(sample.id,function(x)c(paste0(x,"-1"),paste0(x,"-2"))))
  leg_df <- readr::read_delim(leg.fn,delim=" ") %>% dplyr::mutate(allele=paste0(allele0,",",allele1),chrom=chrom)
  geno_mat <- readr::read_delim(hap.fn,delim=" ",col_names = sample.id) %>% dplyr::mutate(ID=leg_df$ID) 
  leg_df <- dplyr::inner_join(leg_df,map_dat)
  geno_mat <- dplyr::semi_join(geno_mat,leg_df)
  geno_mat <- data.matrix(dplyr::select(geno_mat,-ID))
  p <- nrow(geno_mat)
  tfile <- tempfile()
  SNPRelate::snpgdsCreateGeno(gds.fn = tfile,genmat = geno_mat,
                   sample.id = sample.id,
                   snp.id = 1:p,
                   snp.rs.id = leg_df$ID,
                   snp.chromosome = leg_df$chrom,
                   snp.position = leg_df$pos,
                   snp.allele = leg_df$allele,
                   compress.annotation = compress.annotation,
                   compress.geno = compress.geno)

  SeqArray::seqSNP2GDS(gds.fn =tfile ,out.fn =out.gdsfn ,storage.option ="LZ4_RA.fast")
  gds <- SeqArray::seqOpen(out.gdsfn,readonly = F)
  gdsfmt::add.gdsn(gdsfmt::index.gdsn(gds,"annotation/info"),"map",leg_df$map,replace=T,compress="LZ4_RA.fast")  
  SeqArray::seqClose(gds)

}


import_panel_data <- function(temp_gds=NULL,
                              map_df=NULL,
                              output_file=NULL,
                              overwrite=F,
                              parallel=1,
                              ld_break_file=NULL){

  stopifnot(!is.null(output_file),
            !is.null(map_df),
            all(file.exists(temp_gds)),
            !is.null(ld_break_file))

    stopifnot(file.exists(ld_break_file))

  if(file.exists(output_file)){
    stopifnot(overwrite) 
    file.remove(output_file)
  }
 SeqArray::seqMerge(temp_gds,output_file,storage.option="LZ4_RA.fast")
  
  gds <-SeqArray::seqOpen(output_file,readonly = F)
  leg_df <- read_SNPinfo_gds(gds) %>% dplyr::left_join(map_df)
  stopifnot(nrow(leg_df)==length(SeqArray::seqGetData(gds,"variant.id")))


  gdsfmt::add.gdsn(gdsfmt::index.gdsn(gds,"annotation/info"),"map",leg_df$map,replace=T,compress="LZ4_RA.fast")
  cat("Adding Chunk Delimiters\n")
 SeqArray::seqClose(gds)
  add_chunk_gds(output_file,ld_break_file)
}


add_chunk_gds <- function(gds_file,region_bed_file){

  stopifnot(file.exists(region_bed_file),file.exists(gds_file))
  gds <- SeqArray::seqOpen(gds_file,readonly = F)
  region_bed <- readr::read_delim(region_bed_file,delim="\t",trim_ws = T) %>% dplyr::mutate(range_id=as.character(1:n()),chr=gsub("chr","",chr))
#  region_range <- GenomicRanges::GRanges(seqnames=region_bed$chr,ranges=IRanges::IRanges(start=region_bed$start,end = region_bed$stop))
  snp_df <- read_SNPinfo_gds(gds) %>% dplyr::mutate(snp_id=as.character(snp_id))
  match_df <- match_SNP(region_bed,snp_df) %>% dplyr::arrange(as.integer(snp_id))
  stopifnot(all(snp_df$snp_id==match_df$snp_id))
  gdsfmt::add.gdsn(gdsfmt::index.gdsn(gds,"annotation/info"),"LD_chunk",as.integer(match_df$range_id),replace=T,compress="LZ4_RA.fast")
  SeqArray::seqClose(gds)
}

calc_LD_gds <- function(gds,m=85,Ne=11490.672741,cutoff=1e-3){
  map_dat <- SeqArray::seqGetData(gds,var.name="annotation/info/map")
  H <- 2.0-SeqArray::seqGetData(gds,var.name="$dosage")
  stopifnot(max(H)==1)
  return(calcLD(hmata = H,mapa = map_dat,m = m,Ne = Ne,cutoff = cutoff))
}

is_haplo <- function(gds,checkSNPs=F){
  samples <- SeqArray::seqGetData(gds,"sample.id")
  sample_trim <- substr(samples,1,nchar(samples)-2)
  if(all(names(table(table(sample_trim)))=="2")){
    return(T)
  }
  return(F)
}


# allchunk_LDshrink(gds_file,region_ids=NULL,outfiles,m=85,Ne=11490.672741,cutoff=1e-3,evd=T)

    

chunkwise_LDshrink <- function(gds_file,region_id=1,outfile=NULL,m=85,Ne=11490.672741,cutoff=1e-3,evd=T){
  
  stopifnot(file.exists(gds_file),!is.null(outfile))
  gds <- SeqArray::seqOpen(gds_file,readonly = T)
  LD_chunks <- SeqArray::seqGetData(gds,var.name="annotation/info/LD_chunk")
  stopifnot(sum(LD_chunks%in%region_id)>0)
  SeqArray::seqSetFilter(gds,variant.sel=LD_chunks %in% region_id)

  si <- read_SNPinfo_gds(gds) %>% dplyr::mutate(region_id=region_id)
  R <- calc_LD_gds(gds,m=m,Ne=Ne,cutoff=cutoff)
  
  
  RcppEigenH5::write_df_h5(df=si,groupname = "LDinfo",outfile = outfile,deflate_level = 4)
  tls <- RcppEigenH5::h5ls(outfile)
  stopifnot(any("/LDinfo/SNP" %in% tls))
  RcppEigenH5::write_mat_h5(outfile,groupname = "LD",dataname = "R",data = R,deflate_level = 4,doTranspose = F)
  if(evd){
    eigenR <- eigen(R)
    stopifnot(min(eigenR$values)>0)
    RcppEigenH5::write_mat_h5(outfile,groupname="EVD",dataname = "Q",data = eigenR$vectors,deflate_level = 0L)
    RcppEigenH5::write_dvec_h5(outfile,groupname="EVD",dataname = "D",data = eigenR$values,deflate_level = 2L)
  }
  return(dim(R))
}


chunkwise_LDshrink_ldsc <- function(gds_file,chrom=19,out_dir=NULL,m=85,Ne=11490.672741,cutoff=1e-3){

  stopifnot(!is.null(gds_file),
            file.exists(gds_file),
            !is.null(out_dir),dir.exists(out_dir))
  outfile <- file.path(out_dir,paste0(chrom,".l2.ldscore.gz"))
  soutfile <- file.path(out_dir,paste0(chrom,".l2.M_5_50"))
  gds <-SeqArray::seqOpen(gds_file,readonly = T)
  # chroms <-SeqArray::seqGetData(gds,"chromosome")
 SeqArray::seqSetFilterChrom(gds,include=as.character(chrom))
  LD_chunks <-SeqArray::seqGetData(gds,var.name="annotation/info/LD_chunk")
  region_ids <- unique(LD_chunks)
  mfilt <-SeqArray::seqGetFilter(gds)
  resl <- list()
  # pb <- progress::progress_bar$new(total=length(region_ids))
  for(i in 1:length(region_ids)){
    region_id <- region_ids[i]
    tr <- LD_chunks %in% region_id
   SeqArray::seqSetFilter(gds,variant.sel=tr,action="push+intersect")
    si <- read_SNPinfo_ldsc(gds) 
    R <- calc_LD_gds(gds,m = m,Ne = Ne,cutoff = cutoff)
    si <- dplyr::mutate(si,L2=colSums(R^2)-1)
    resl[[i]] <- si
   SeqArray::seqSetFilter(gds,action="pop")
    tfilt <-SeqArray::seqGetFilter(gds)
    stopifnot(all.equal(mfilt,tfilt))
    # pb$tick()
  }
  resdf <- bind_rows(resl)
  readr::write_delim(resdf,path=outfile,delim = "\t")
  nc <- filter(resdf,MAF>0.05) %>% nrow
  write(x = nc,file = soutfile)
  return(dim(R))
}




