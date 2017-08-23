#Wrappers for ldetect


part_chrom <- function(map_file_gz,out_file,panel_size){
  script_file <- system.file("ldetect_files/partition_chromosome.py", package = "LDshrink") 
  sys_res <- system(command = paste("python3",script_file,map_file_gz,panel_size,out_file,sep = " "))
}

# calc_cov_ld <- function(ref_panel,map_file_gz,cov_gz,
#                         ind_list,
#                         outfile,
#                         Ne=11418,
#                         cutoff=1e-7){
#   
# }

ldetect_partition <- function(map_dat,panel_size=379){
  tmap_file <- tempfile()
  map_file <- gzfile(tmap_file,open = "w")
  readr::write_delim(map_dat,path = map_file,delim = " ",col_names = F)
  close(map_file)
  out_file <- tempfile()
  part_chrom(map_file_gz = tmap_file, out_file = out_file, panel_size = panel_size)
  splits <- readr::read_delim(out_file,delim = " ", col_names = c("startpos","endpos"))
  file.remove(c(out_file,tmap_file))
  return(tibble::as_data_frame(splits))
}  