haplo_reader <- function(filename,region_id,...){
  # gfile <- createfn.gds("test.gds")
  gfile <- gdsfmt::openfn.gds(filename = filename,readonly = T,allow.duplicate = T)
  snp_df <- tibble::data_frame(
                               SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "annotation/id")), 
                               snp_id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "variant.id")), 
                               chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "chromosome")), 
                               pos = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gfile, "position")))
}
