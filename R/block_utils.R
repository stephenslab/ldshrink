
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
    dplyr::inner_join(snp_df) %>%
    dplyr::mutate(min_size=ifelse(chrom_min<min_block_size,chrom_min,min_block_size)) %>%
    dplyr::select(-chrom_min)

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


assign_map <- function(snp_df,map_df,parallel=F){
  u_chr <- dplyr::distinct(snp_df,chr)
  snp_dfl <- split(snp_df,snp_df$chr)
  map_dfl <- dplyr::semi_join(map_df,u_chr,by="chr") %>% split(.$chr)
  stopifnot(all(names(map_dfl)==names(snp_dfl)))
  if(parallel){
    retl <- purrr::map2(map_dfl,snp_dfl,function(x,y){
      return(future::future(interpolate_map(x$map,x$pos,y$pos),packages="LDshrink"))
    })
    tot_c <- length(retl)
    cat("Futures assigned\n")
    tot_f <- sum(purrr::map_lgl(retl,resolved))
    pb <- progress::progress_bar$new(total = tot_c)
    while(tot_f<tot_c){
      tot_f <- sum(purrr::map_lgl(retl,future::resolved))
      pb$update(tot_f/tot_c)
      Sys.sleep(0.5)
    }
    retdf <- map2_df(snp_dfl,retl,~dplyr::mutate(.x,map=future::value(.y)))
  }else{
    retdf <- purrr::map2_df(map_dfl,snp_dfl,~dplyr::mutate(.y,map=interpolate_map(.x$map,.x$pos,.y$pos)))
  }
  return(retdf)
}











LDshrink_evd <- function(panel,map=NULL,m=85,
                            Ne=11490.672741,
                            cutoff=1e-3,
                         useLDshrink=T,na.rm=T){
#    stopifnot(ncol(panel)==len
    isGeno <- max(panel,na.rm = na.rm)>1

    # S <- cov(scale(panel,center=T,scale=F))
    if(useLDshrink){
        S <- LDshrink(haplo_panel = panel,
                      map_data = map,
                      m = m,
                      Ne = Ne,
                      cutoff = cutoff,
                      cov_2_cor = T,
                      na.rm = T)
    }else{
      S <- cor(panel,use = "complete.obs")
    }
    L2 <- colSums(S^2)-1
    evdR <- eigen(S)
    return(list(R=S,L2=L2,D=evdR$values,Q=evdR$vectors))
}
    


chunkwise_LD_h5 <- function(input_file,
                            output_file,
                            snp_df,
                            m=85,
                            Ne=11490.672741,
                            cutoff=1e-3,
                            useLDshrink=T){


}
        
    ## if(SNPfirst){
    ##     dosage_l <- EigenH5::read_mat_row_futures(h5filename = input_file,filepath="dosage",subset_rows = split(snp_df$snp_id,snp_df$region_id),doTranspose = T)
    ## }else{
    ##     dosage_l <- EigenH5::read_mat_col_futures(h5filename = input_file,filepath="dosage",subset_cols = split(snp_df$snp_id,snp_df$region_id),doTranspose = F)
    ## }
    ## retl <- purrr::map2(dosage_l,map(snp_dfl,"map"),~future::future({chunkwise_LD(future::value(.x),.y,m=m,Ne=Ne,cutoff=cutoff,LDshrink=LDshrink)}))

    

  ##       cov(scale(future::value(.x),center=T,scale=F))))
  ## if(LDshrink){
  ##   map_l <- split(snp_df$map,snp_df$region_id)
  ##   R_l <- purrr::map2(cov_l,map_l,~future::future({

  ## }
 
  # iwalk(R_l,function(x){
  #   EigenH5::write_matrix_h5(filename = output_file,groupname = "")
  # })
  
  
  # 
  # input_dff <- EigenH5::split_chunk_df(snp_df,pos_id=snp_id,group_id=region_id,rowsel = SNPfirst,colsel = !SNPfirst) %>% dplyr::mutate(chunk_group=region_id) %>% dplyr::filter(!is.na(chunk_group))
  # map_dff <- EigenH5::split_chunk_df(snp_df,pos_id=ld_snp_id,group_id=region_id) %>% dplyr::mutate(chunk_group=region_id) %>% dplyr::filter(!is.na(chunk_group))
  # map_dff <- dplyr::mutate(map_dff,filenames=mapf,groupnames="SNPinfo",datanames="map") %>% dplyr::mutate(col_offsets=0L,col_chunksizes=1L)
  # data_dff <- dplyr::mutate(input_dff,filenames=input_file,groupnames="/",datanames="dosage")
  # SNP_dff <-  dplyr::mutate(input_dff,filenames=input_file,groupnames="SNPinfo",datanames="SNP",col_offsets=0L,col_chunksizes=1L)
  # 
  # EigenH5::write_df_h5(df=snp_df,groupname = "LDinfo",outfile = output_file)
  # 
  # if(SNPfirst){
  #   data_dff <- data_dff%>% dplyr::mutate(col_offsets=0L,col_chunksizes=N)
  # }else{
  #   data_dff <- data_dff%>% dplyr::mutate(row_offsets=0L,row_chunksizes=N)
  # }
  # 
  # in_dff <- dplyr::bind_rows(map_dff,data_dff)
  # 
  # in_dff <- dplyr::bind_rows(in_dff,SNP_dff)
  # chunksize_max<-1024
  # if(SNPfirst){
  #   output_dff <-dplyr::mutate(input_dff,create_dynamic=F,row_offsets=0L,col_offsets=0L,datatypes="numeric",col_chunksizes=row_chunksizes,
  #                              col_c_chunksizes=pmin(chunksize_max,col_chunksizes),
  #                              row_c_chunksizes=col_c_chunksizes)
  # }else{
  #   output_dff <-dplyr::mutate(input_dff,create_dynamic=F,row_offsets=0L,col_offsets=0L,datatypes="numeric",row_chunksizes=col_chunksizes,
  #                              row_c_chunksizes=pmin(chunksize_max,row_chunksizes),
  #                              col_c_chunksizes=row_c_chunksizes)
  # }
  ## R_dff <- dplyr::mutate(output_dff,filenames=output_file,groupnames=paste0("LD/",region_id),datanames="R")
  ## L2_dff <- dplyr::mutate(output_dff,filenames=output_file,groupnames=paste0("L2/",region_id),datanames="L2",col_chunksizes=1L,col_c_chunksizes=1L)
  ## Q_dff <- dplyr::filter(output_dff,evd) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("EVD/",region_id),datanames="Q")
  ## D_dff <- dplyr::filter(output_dff,evd) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("EVD/",region_id),datanames="D",col_chunksizes=1L,col_c_chunksizes=1L)
  ## svd_d_dff <- dplyr::filter(output_dff,svd) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("SVD/",region_id),
  ##                                                              datanames="d",
  ##                                                              col_chunksizes=1L,
  ##                                                              col_c_chunksizes=1L,
  ##                                                              create_dynamic=T)
 
  ## svd_v_dff <- dplyr::filter(output_dff,svd) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("SVD/",region_id),datanames="V")
  

  ## rowsnp_dff <- dplyr::filter(output_dff,df) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("LD_DF/",region_id,"/LD"),
  ##                                                              datanames="rowSNP",
  ##                                                              col_chunksizes=1L,
  ##                                                              col_c_chunksizes=1L,
  ##                                                              create_dynamic=T)
  ## colsnp_dff <- dplyr::filter(output_dff,df) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("LD_DF/",region_id,"/LD"),
  ##                                                              datanames="colSNP",
  ##                                                              col_chunksizes=1L,
  ##                                                              col_c_chunksizes=1L,
  ##                                                              create_dynamic=T)
  ## r2_dff <- dplyr::filter(output_dff,df) %>% dplyr::mutate(filenames=output_file,groupnames=paste0("LD_DF/",region_id,"/LD"),
  ##                                                          datanames="r2",
  ##                                                          col_chunksizes=1L,
  ##                                                          col_c_chunksizes=1L,
  ##                                                          create_dynamic=T)
  
  ## out_data_dff <- dplyr::bind_rows(R_dff,
  ##                                  L2_dff,
  ##                                  Q_dff,
  ##                                  D_dff,
  ##                                  rowsnp_dff,
  ##                                  colsnp_dff,
  ##                                  r2_dff,
  ##                                  svd_v_dff,
  ##                                  svd_d_dff)

  
  ## dplyr::filter(out_data_dff,!create_dynamic) %>% EigenH5::create_mat_l()
  
  ## calc_LD_chunk_h5(input_dff = in_dff,output_dff = out_data_dff,m=m,Ne=Ne,cutoff=cutoff,SNPfirst=SNPfirst,evd=evd,svd=svd,df=df,r2cutoff=r2cutoff)
#}


#' #' Title
#' #'
#' #' @param snp_df_1 
#' #' @param snp_df_2 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' merge_snpsets <- function(snp_df_1,snp_df_2){
#'   
#'   intersect_snp <- dplyr::inner_join(snp_df_1,snp_df_2,by=c("SNP","chr","pos"))
#'   
#'   
#' }





