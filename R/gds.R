read_SNPinfo_gds <- function(gds,alleles=F){
  if(!alleles){
    return(tibble::data_frame(SNP=seqGetData(gds,var.name="annotation/id"),
                              snp_id=seqGetData(gds,var.name="variant.id"),
                              chr=seqGetData(gds,var.name="chromosome"),
                              pos=seqGetData(gds,var.name="position")))
  }else{
    return(tibble::data_frame(SNP=seqGetData(gds,var.name="annotation/id"),
                              snp_id=seqGetData(gds,var.name="variant.id"),
                              chr=seqGetData(gds,var.name="chromosome"),
                              pos=seqGetData(gds,var.name="position"),
                              allele=seqGetData(gds,var.name="allele")))
  }
  
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
                              panel_ind_file=NULL,
                              ld_break_file=NULL,
                              subset_MAF=0,remove_multi_allelic=T,snv_only=T){
  library(SeqArray)
  library(dplyr)
  library(purrr)
  library(readr)
  stopifnot(!is.null(output_file),
            !is.null(panel_ind_file),
            file.exists(panel_ind_file),
            file.exists(temp_gds),
            !is.null(map_df),!is.null(ld_break_file))

    stopifnot(file.exists(ld_break_file))

  if(file.exists(output_file)){
    stopifnot(overwrite)
    file.remove(output_file)
  }
  sub_ind=scan(panel_ind_file,what=character())
  gc()
#  temp_gds <- tempfile()

  gds <- seqOpen(temp_gds,readonly = F)
  seqSetFilter(gds,sample.id=sub_ind)
  seqSetFilterCond(gds,maf=subset_MAF,parallel = parallel,.progress = T)

  gds_snpinfo <- read_SNPinfo_gds(gds,alleles = T) %>% mutate(Num_Alleles=seqNumAllele(gds),is_SNV=is_SNV(allele))
  
  if(remove_multi_allelic){
    gds_snpinfo <- filter(gds_snpinfo,Num_Alleles==2)
  }
  if(snv_only){
    gds_snpinfo <- filter(gds_snpinfo,is_SNV)
  }
  b_gds_snpinfo <- inner_join(gds_snpinfo,map_df)
  
  seqSetFilter(gds,variant.id=b_gds_snpinfo$snp_id,action="intersect")
  seqExport(gdsfile = gds,out.fn = output_file)
  seqClose(gds)
  cat("Adding mapAnnotation\n")
  gds <- seqOpen(output_file,readonly = F)
  gdsfmt::add.gdsn(index.gdsn(gds,"annotation/info"),"map",b_gds_snpinfo$map,replace=T,compress="LZ4_RA.fast")
  seqClose(gds)
  cat("Adding Chunk Delimiters\n")
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
chunkwise_LDshrink <- function(gds_file,region_id=1,outfile=NULL,m=85,Ne=11490.672741,cutoff=1e-3,evd=T){
  library(SeqArray)
  library(RcppEigenH5)
  library(dplyr)
  stopifnot(file.exists(gds_file),!is.null(outfile))
  gds <- seqOpen(gds_file,readonly = T)
  LD_chunks <- seqGetData(gds,var.name="annotation/info/LD_chunk")
  stopifnot(sum(LD_chunks%in%region_id)>0)
  seqSetFilter(gds,variant.sel=LD_chunks %in% region_id)
  map_dat <- seqGetData(gds,var.name="annotation/info/map")
  si <- read_SNPinfo_gds(gds) %>% mutate(region_id=region_id)

  H <- matrix(as.numeric(c(seqGetData(gds,var.name="genotype"))),length(map_dat),byrow = T)
  stopifnot(max(H)==1)

  R <- calcLD(hmata = t(H),mapa = map_dat,m = m,Ne = Ne,cutoff = cutoff)
  RcppEigenH5::write_df_h5(df=si,groupname = "LDinfo",outfile = outfile,deflate_level = 4)
  RcppEigenH5::write_mat_h5(outfile,groupname = "LD",dataname = "R",data = R,deflate_level = 4,doTranspose = F)
  if(evd){
    eigenR <- eigen(R)
    stopifnot(min(eigenR$values)>0)
    write_mat_h5(outfile,groupname="EVD",dataname = "Q",data = eigenR$vectors,deflate_level = 0L)
    write_dvec_h5(outfile,groupname="EVD",dataname = "D",data = eigenR$values,deflate_level = 2L)
  }
  return(dim(R))
}



