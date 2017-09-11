
library(dplyr)
library(readr)
library(RSSp)

inf <- snakemake@input[["logf"]]
tpf <- snakemake@input[["tparamf"]]
fgeneid <- snakemake@params[["fgeneid"]]
outf <- snakemake@output[["logf"]]

#save.image()
res_df <- map2_dfr(inf,fgeneid,function(filen,genen){
    parse_ldsc_h2log(filen) %>%
        mutate(fgeneid = as.character(genen)) %>% filter(Variable!="Ratio")%>% select(-SD)
}
) %>% spread(Variable,Est)


res_df <- map2_dfr(tpf,fgeneid,function(filen,genen){
    read_delim(filen, delim = "\t") %>%
        mutate(fgeneid = as.character(genen))
}
)%>%inner_join(res_df)
write_delim(res_df, path = outf, delim = "\t")
