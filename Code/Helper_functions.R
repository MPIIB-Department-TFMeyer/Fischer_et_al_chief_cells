#########################################################
# assumes that dge_results is a data.table or data.frame object
dedup_DGE_results <- function(dge_results, by="pvalue", logfc_col = "logFC", pval_col = "P.Value", adjp_col = "adj.P.Val", id_col = "GeneSymbol", probeid_col="ProbeName") {
  probes = dge_results[[probeid_col]]
  symbol = dge_results[[id_col]]
  logfc = dge_results[[logfc_col]]
  pval = dge_results[[pval_col]]
  adjp = dge_results[[adjp_col]]
  
  probes_by_symbol = tapply(probes, symbol, function(x) x)
  if (by=="logfc") {
    max_fc = function(x) {which(abs(x) == max(abs(x),na.rm=T))}
    max_fc_ind_per_symbol = lapply(split(logfc,symbol), max_fc)
    sel_probes = mapply( function(x,y) x[y], probes_by_symbol, max_fc_ind_per_symbol)
  } else if (by == "pvalue") {
    minp = function(x) {which(x == min(x,na.rm=T))}
    minp_ind_per_symbol = lapply(split(pval,symbol), minp)
    sel_probes = mapply( function(x,y) x[y], probes_by_symbol, minp_ind_per_symbol)    
  }
  
  return( dge_results[dge_results[[probeid_col]] %in% unlist(sel_probes), ] )
}

suppressMessages(require(RColorBrewer))
RedWhtBlue_orig = brewer.pal(n = 7, name = "RdBu")
RedGyBlue = RedWhtBlue_orig
RedGyBlue[4] = "grey90"
RedGreyBluePal = colorRampPalette(rev(RedGyBlue))(100)

RedGyBlue2 = brewer.pal(n = 7, name = "RdYlBu")
RedGyBlue2[4] = "grey90"
RedGreyBluePal2 = colorRampPalette(rev(RedGyBlue2))(100)

RedGreyBluePal3 = colorRampPalette(c(brewer.pal(n = 7, name = "RdYlBu")[7],  "grey90" ,brewer.pal(n = 7, name = "RdYlBu")[1]))(100)

RedWhtBlue_orig = brewer.pal(n = 7, name = "RdBu")
RedGyBlue_light = RedWhtBlue_orig
RedGyBlue_light[4] = "grey90"
RedGyBlue_light[1] = "red"
RedGyBlue_light[7] = "blue"
RedGreyBlueLightPal = colorRampPalette(rev(RedGyBlue_light))(100)

darker_rgrb = apply(round(t(col2rgb(RedGyBlue_light)) / 1.5), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
RedGreyBlueDarkPal = colorRampPalette(rev(darker_rgrb))(100)


suppressMessages(library(org.Hs.eg.db))
entrez_to_symbol = function(entrezids) {
  tmp = as.character(entrezids)
  tmp = tmp[!is.na(tmp)]
  tmp_df= dcast(select(org.Hs.eg.db, keys=tmp, keytype="ENTREZID", columns="SYMBOL"), ENTREZID ~ ., value.var="SYMBOL", fun.aggregate = function(x) paste(x[!is.na(x)], collapse=",") )
  colnames(tmp_df) = c("EntrezID", "GeneSymbol")  
  return(tmp_df)
}

symbol_to_entrez = function(symbols) {
  tmp = as.character(symbols)
  tmp = tmp[!is.na(tmp)]
  tmp_df= dcast(select(org.Hs.eg.db, keys=tmp, keytype="SYMBOL", columns="ENTREZID"), SYMBOL ~ ., value.var="ENTREZID", fun.aggregate = function(x) paste(x[!is.na(x)], collapse=",") )
  colnames(tmp_df) = c("GeneSymbol", "EntrezID")  
  return(tmp_df)
}

probename_to_symbol_agilent = function(probenames) {
  tmp = normalized$genes[normalized$genes$ProbeName %in% probenames, c("ProbeName", "GeneSymbol")]
  rownames(tmp) = tmp$ProbeName
  return(tmp[probenames,])
}

probename_to_entrez_agilent = function(probenames) {
  tmp = normalized$genes[normalized$genes$ProbeName %in% probenames, c("ProbeName", "EntrezID")]
  rownames(tmp) = tmp$ProbeName
  return(tmp[probenames,])
}



