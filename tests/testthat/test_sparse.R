context("sparse")
# test_that("Sparse and Dense implementations are equivalent for non-ldshrink", {
#   n <- 500
#   p <- 1100
#   
#   haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
#   tmap <- runif(p)
#   
#   mapdat <- cumsum(tmap)
#   
#   tLD <- cor(haplomat)
#   sLD <- ldshrink::sparse_ldshrink(data = haplomat,
#                                    mapd = mapdat,
#                                    indices=0:(p-1),
#                                    m=85,
#                                    Ne=11490.672741,
#                                    cutoff = 0.001,
#                                    total_size=p,
#                                    useldshrink = F)
#   dsLD <- as.matrix(sLD)
#   
# 
#   expect_equivalent(tLD, dsLD)
# })
# 
# 
# 
# test_that("Sparse and Dense implementations are equivalent", {
#   n <- 5000
#   p <- 1100
#   
#   haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
#   nhaplomat <- haplomat+0
#   tmap <- runif(p)
#   
#   mapdat <- cumsum(tmap)
#   
#   tLD <- ldshrink::ldshrink(haplomat, mapdat, na.rm = F, m=85, Ne=11490.672741, cutoff = 0.001)
#   sLD <- ldshrink::sparse_ldshrink(data = nhaplomat, indices=0:(p-1),total_size=p,
#                          mapd = mapdat, m=85, Ne=11490.672741, cutoff = 0.001, useldshrink=T, progress=F)
#   dsLD <- as.matrix(sLD)
#   
#   summary(c(tLD-dsLD))
#   expect_equivalent(tLD, dsLD, 1e-4)
# })
# 
# 
# 
# 
# 
# 
# test_that("Sparse implementations are equivalent even when chunked", {
#   n <- 501
#   p <- 1100
#   haplomat <- matrix(sample(0:1, n*p, replace = T), n, p)
#   tmap <- runif(p)
#   mapdat <- cumsum(tmap)
#   rs <- as.character(1:p)
#   shap <- scale(haplomat, center=T, scale=F)
#   sLD <- ld2df(data = haplomat,
#                mapd = mapdat,
#                 rsid = as.character(1:p),
#                m=85, Ne=11490.672741, cutoff = 0.001, r2cutoff = 0, useldshrink=T) %>%
#     dplyr::arrange(rowsnp, colsnp, r)
#   hap_a <- haplomat[, 1:500]
#   hap_b <- haplomat[, -(1:500)]
#   map_a <- mapdat[1:500]
#   map_b <- mapdat[-(1:500)]
#   rs_a <- rs[1:500]
#   rs_b <- rs[-(1:500)]
#   
#   df_ab <- ld2df_p(data_a = hap_a,
#                    data_b = hap_b,
#                    mapd_a = map_a,
#                    mapd_b = map_b,
#                    rsid_a = rs_a,
#                    rsid_b = rs_b,
#                    m=85,
#                    Ne=11490.672741,
#                    cutoff = 0.001, r2cutoff = 0,
#                    useldshrink = T)
#   df_aa <- ld2df(data = hap_a,
#                  mapd = map_a,
#                  rsid = rs_a,
#                  m=85,
#                  Ne=11490.672741,
#                  cutoff = 0.001, r2cutoff = 0,
#                  useldshrink=T)
#   df_bb <- ld2df(data = hap_b,
#                  mapd = map_b,
#                  rsid = rs_b,
#                  m=85,
#                  Ne=11490.672741,
#                  cutoff = 0.001,r2cutoff = 0,
#                  useldshrink=T)
#   df_p <- dplyr::bind_rows(df_aa,df_ab,df_bb) %>% dplyr::arrange(rowsnp,colsnp,r)
#   expect_equal(df_p,sLD)
# })
