
#' Collapse small LD blocks into adjacent blocks
#'
#' @param snp_df a dataframe with three (integer) columns, `chr`,`pos`, and `region_id`,
#'  specifying the position of each SNP as well as it's current LD block assignment
#' @param min_block_size the minimum number of SNPs per block.  
#' Blocks with fewer than this number will be merged with adjacent blocks (blocks are not merged across chromosome boundaries)
#' @return a copy of the dataframe, where the `region_id` vector has been modifed to reflect the merged LD blocks
#' @export
#'
#' @examples
collapse_ld_region <- function(snp_df,break_df=NULL,min_block_size=1){
  if(is.null(break_df)){
    data("break_df")
  }
  break_df <- spread_ld_region(break_df)
  tot_chrom <- unique(break_df$chr)
  tbreak_df <- tibble::data_frame(chr=tot_chrom)
  # chrom_ct <- 

  snp_df <- dplyr::group_by(snp_df,chr) %>% 
    dplyr::summarise(chrom_min=n()) %>% 
    dplyr::inner_join(snp_df) %>% dplyr::mutate(min_size=ifelse(chrom_min<min_block_size,chrom_min,min_block_size)) %>% dplyr::select(-chrom_min)
  # if(any(chrom_ct$n_snp<min_block_size)){
  #   stop("not all chrom have enough SNPs for min_block_size")
  # }
  nbreak_df <- dplyr::mutate(break_df,next_start=dplyr::lead(start,default=.Machine$integer.max),prev_stop=dplyr::lag(stop,default = 0))
  snp_ct <- dplyr::group_by(snp_df,chr,region_id,min_size)%>% dplyr::summarise(snp_ct=n()) %>% dplyr::ungroup() %>% dplyr::filter(snp_ct<min_size)
  while(nrow(snp_ct)>0){
    nbreak_df <- dplyr::anti_join(nbreak_df,snp_ct) %>% 
      dplyr::full_join(tbreak_df) %>%
      dplyr::mutate(start=ifelse(is.na(start),0,start),
                    stop=ifelse(is.na(stop),.Machine$integer.max,stop)
                    ) %>% 
      dplyr::arrange(chr,start,stop) %>% 
      dplyr::mutate(region_id=1:n()) %>%
      dplyr::group_by(chr) %>% 
      dplyr::mutate(next_start=dplyr::lead(start,default = .Machine$integer.max),
                    prev_stop=dplyr::lag(stop,default = 0))  %>% 
      dplyr::ungroup()
    stopifnot(all(nbreak_df$next_start>=nbreak_df$start),all(nbreak_df$prev_stop<=nbreak_df$start))
    
    stot_chrom <- length(unique(nbreak_df$chr))
    stopifnot(stot_chrom==length(tot_chrom))

    nbreak_df <- nbreak_df %>%
      dplyr::mutate(stop=ifelse(next_start==stop,stop,ifelse(next_start==.Machine$integer.max,.Machine$integer.max,stop+(next_start-stop)/2))
                    ,start=ifelse(prev_stop==start,start,ifelse(prev_stop==0,0,start-(start-prev_stop)/2))) %>%
      dplyr::arrange(chr,start,stop,region_id) %>% dplyr::select(-next_start,-prev_stop)
    snp_df <- assign_snp_block(snp_df,nbreak_df,assign_all = T)
    snp_ct <- dplyr::group_by(snp_df,chr,region_id,min_size)%>% dplyr::summarise(snp_ct=n()) %>% dplyr::ungroup() %>% dplyr::filter(snp_ct<min_size)
  }
  return(dplyr::select(snp_df,-min_size))
}

#' Expand LD region boundaries to fill gaps
#'
#' @param break_df 
#'
#' @return
#' @export
#'
#' @examples
spread_ld_region <- function(nbreak_df){
  break_dfl <- split(nbreak_df,nbreak_df$chr)
  nbreak_df <- purrr::map_df(break_dfl,function(x){
    dplyr::mutate(x,start=ifelse(start==min(start),0,start),stop=ifelse(stop==max(stop),.Machine$integer.max,stop))
  })
  
  return(nbreak_df)
}





#' #' Thin genome by removing adjacent LD blocks
#' #'
#' #' @param nbreak_df 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' thin_ld_region <- function(nbreak_df){
#'   
#' }

#' Assign SNPs to LD blocks
#'
#' @param snp_df 
#' @param break_df dataframe of LD blocks, if `NULL`, use precomputed (EUR) ldshrink LD blocks
#' @param assign_all whether to throw an error if a SNP cannot be assigned to a block, or to assign it to block `NA`
#'
#' @return
#' @export
#'
#' @examples
assign_snp_block <- function(snp_df,break_df=NULL,assign_all=T){
  if(is.null(break_df)){
    data("break_df")
  }
  return(dplyr::mutate(snp_df,region_id=set_ld_region(break_df,snp_df,assign_all = assign_all)))
}


assign_map <- function(snp_df,map_df){
  u_chr <- dplyr::distinct(snp_df,chr)
  snp_dfl <- split(snp_df,snp_df$chr)
  map_dfl <- dplyr::semi_join(map_df,u_chr) %>% split(.$chr)
  stopifnot(all(names(map_dfl)==names(snp_dfl)))
  return(purrr::map2_df(map_dfl,snp_dfl,~dplyr::mutate(.y,map=interpolate_map(.x$map,.x$pos,.y$pos))))
}




filter_pvv <- function(D_df,pvv=1){
  dplyr::group_by(D_df,region_id) %>%
    dplyr::mutate(rel_D=D/sum(D),cumsum_D=cumsum(rel_D),chunk_id=1:n(),chunksize=n()) %>%
    dplyr::filter(cumsum_D<=max(c(pvv,min(cumsum_D)))) %>%
    dplyr::ungroup() %>%return()
    # dplyr::summarise(chunk_size=first(which(cumsum_D>=pvv)),total_size=n()) %>%
}




chunkwise_LDshrink_h5 <- function(input_file,
                                  output_file,
                                  snp_df,
                                  m=85,
                                  Ne=11490.672741,
                                  cutoff=1e-3,
                                  evd=T,
                                  svd=T,
                                  df=F,
                                  r2cutoff=0.01){
  
  stopifnot(file.exists(input_file),
            !file.exists(output_file),
            !is.null(snp_df[["region_id"]]),
            !is.null(snp_df[["map"]]))
  if(is.null(snp_df[["snp_id"]])){
    snp_df <- dplyr::mutate(snp_df,snp_id=1:n())
  }
  p <- nrow(snp_df)
  
  dosage_dims <-EigenH5::get_dims_h5(input_file,"/","dosage")
  pm <- dosage_dims-p
  stopifnot(!all(pm<0))
  SNPfirst <-  which.min(abs(pm))==1
  mapf <- tempfile()
  EigenH5::write_vector_h5(filename = mapf,groupname = "SNPinfo",dataname = "map",data = snp_df$map)
  # EigenH5::write_vector_h5(filename = mapf,groupname = "SNPinfo",dataname = "SNP",data = snp_df$SNP)
  
  snp_df <- dplyr::mutate(snp_df,ld_snp_id=1:n())
  
  input_dff <- EigenH5::split_chunk_df(snp_df,pos_id=snp_id,group_id=region_id) %>% dplyr::mutate(chunk_group=region_id) %>% dplyr::filter(!is.na(chunk_group))
  map_dff <- EigenH5::split_chunk_df(snp_df,pos_id=ld_snp_id,group_id=region_id) %>% dplyr::mutate(chunk_group=region_id) %>% dplyr::filter(!is.na(chunk_group))
  map_dff <- dplyr::mutate(map_dff,filenames=mapf,groupnames="SNPinfo",datanames="map") %>% dplyr::mutate(col_offsets=0L,col_chunksizes=1L)
  data_dff <- dplyr::mutate(input_dff,filenames=input_file,groupnames="/",datanames="dosage")
  SNP_dff <-  dplyr::mutate(input_dff,filenames=input_file,groupnames="SNPinfo",datanames="SNP",col_offsets=0L,col_chunksizes=1L)
  
  EigenH5::write_df_h5(df=snp_df,groupname = "LDinfo",outfile = output_file)
  
  if(SNPfirst){
    data_dff <- data_dff%>% dplyr::select(-dplyr::starts_with("col_"))
  }else{
    data_dff <- data_dff%>% dplyr::select(-dplyr::starts_with("row_"))
  }
  
  in_dff <- dplyr::bind_rows(map_dff,data_dff)
  in_dff <- dplyr::mutate(in_dff,row_offsets=0L,row_chunksizes=-1L)
  in_dff <- dplyr::bind_rows(in_dff,SNP_dff)
  
  output_dff <-dplyr::mutate(input_dff,create_dynamic=F,row_offsets=0L,col_offsets=0L,row_c_chunksizes=row_chunksizes,col_c_chunksizes=col_chunksizes,datatypes="numeric")
  R_dff <- dplyr::mutate(output_dff,filenames=output_file,groupnames=paste0("LD/",region_id),datanames="R")
  L2_dff <- dplyr::mutate(output_dff,filenames=output_file,groupnames=paste0("L2/",region_id),datanames="L2",col_chunksizes=1L,col_c_chunksizes=1L)
  Q_dff <- dplyr::filter(output_dff,evd) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("EVD/",region_id),datanames="Q")
  D_dff <- dplyr::filter(output_dff,evd) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("EVD/",region_id),datanames="D",col_chunksizes=1L,col_c_chunksizes=1L)
  

  rowsnp_dff <- dplyr::filter(output_dff,df) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("LD_DF/",region_id,"/LD"),
                                                               datanames="rowSNP",
                                                               col_chunksizes=1L,
                                                               col_c_chunksizes=1L,
                                                               create_dynamic=T)
  colsnp_dff <- dplyr::filter(output_dff,df) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("LD_DF/",region_id,"/LD"),
                                                               datanames="colSNP",
                                                               col_chunksizes=1L,
                                                               col_c_chunksizes=1L,
                                                               create_dynamic=T)
  r2_dff <- dplyr::filter(output_dff,df) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("LD_DF/",region_id,"/LD"),
                                                           datanames="r2",
                                                           col_chunksizes=1L,
                                                           col_c_chunksizes=1L,
                                                           create_dynamic=T)
  
  out_data_dff <- dplyr::bind_rows(R_dff,L2_dff,Q_dff,D_dff,rowsnp_dff,colsnp_dff,r2_dff)

  
  dplyr::filter(out_data_dff,!create_dynamic) %>% EigenH5::create_mat_l()
  
  
  
  calc_LD_chunk_h5(input_dff = in_dff,output_dff = out_data_dff,m=m,Ne=Ne,cutoff=cutoff,SNPfirst=SNPfirst,evd=evd,df=df,r2cutoff=r2cutoff)
  if(svd){
    run_svd <-function(filenames,groupnames,datanames,chunk_group,row_offsets=0,col_offsets=0,col_chunksizes=NULL,row_chunksizes=NULL,...){
      tX <- EigenH5::read_mat_h5(filenames,
                                 groupnames,
                                 datanames,
                                 offset_rows = row_offsets,
                                 offset_col=col_offsets,
                                 chunksize_rows=row_chunksizes,
                                 chunksize_cols=col_chunksizes)
      if(SNPfirst){
        tX <- t(tX)
      }
      tX <- scale(tX,center=T,scale=T)
      tp <- ncol(tX)
      svdR <- svd(tX)
      # if(tp<=3){
      #   svdR <- svd(tX)
      # }else{
      #   svdR <-RSpectra::svds(tX,k=tp,nu=0,nv=tp)
      # }
      EigenH5::write_vector_h5(filename = output_file,groupname = paste0("SVD/",chunk_group),dataname="d",svdR$d)
      EigenH5::write_matrix_h5(filename = output_file,groupname = paste0("SVD/",chunk_group),dataname="V",svdR$v)
      return(T)
    }
    trun_svd <- purrr::lift_dl(run_svd)
    dplyr::rowwise(data_dff) %>% dplyr::do(tibble::data_frame(res=trun_svd(.)))
  }
}






