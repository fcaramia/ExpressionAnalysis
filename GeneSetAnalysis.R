library(GSVA)
library(Biobase)
library(GSEABase)
library(GSA)
library(EGSEA)
library(biomaRt)


GetHumanGeneSymbolsfromENSEMBL <- function(ensembl.symbols = NULL)
{
  mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")
  res = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
              filters = "ensembl_gene_id", values = ensembl.symbols, mart = mart.hs)
  
  return(res)
}

GetMouseOrtho <- function(human.symbols = NULL)
{
  mart.hs <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensembl.ids = getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
              filters = "hgnc_symbol", values = human.symbols, mart = mart.hs)
  ensembl.ids = getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), 
                      filters = "ensembl_gene_id", values = ensembl.ids$ensembl_gene_id, mart = mart.hs)
  
  
  mart.hs <- useMart("ensembl", "mmusculus_gene_ensembl")
  res = getBM(attributes = c("mgi_symbol", "ensembl_gene_id"), 
              filters = "ensembl_gene_id", values = ensembl.ids$mmusculus_homolog_ensembl_gene, mart = mart.hs)
  
  return(res)
}

GetHumanOrtho <- function(mouse.symbols = NULL)
{
  mart.hs <- useMart("ensembl", "mmusculus_gene_ensembl")
  ensembl.ids = getBM(attributes = c("mgi_symbol", "ensembl_gene_id"), 
              filters = "mgi_symbol", values = mouse.symbols, mart = mart.hs)
  
  ensembl.ids = getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), 
                      filters = "ensembl_gene_id", values = ensembl.ids$ensembl_gene_id, mart = mart.hs)
  
  mart.hs <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  res = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                      filters = "ensembl_gene_id", values = ensembl.ids$hsapiens_homolog_ensembl_gene, mart = mart.hs)
  
  
  
  
  return(res)
}

GetEntrezIDs <- function(symbols = NULL)
{
  mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")
  res = getBM(attributes = c("hgnc_symbol", "entrezgene"), 
              filters = "hgnc_symbol", values = symbols, mart = mart.hs)
  
  return(res)
}

GetEntrezIDsMouse <- function(symbols = NULL)
{
  mart.hs <- useMart("ENSEMBL_MART_MOUSE", "mc57bl6nj_gene_ensembl")
  res = getBM(attributes = c("mgi_symbol", "entrezgene"), 
              filters = "mgi_symbol", values = symbols, mart = mart.hs)
  
  return(res)
}

GetGeneSymbolsHuman <- function(entrezgenes = NULL)
{
  mart.hs <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl",verbose = T,host = 'asia.ensembl.org')
  res = getBM(attributes = c("hgnc_symbol", "entrezgene"), 
              filters = "entrezgene", values = entrezgenes, mart = mart.hs)
  
  
  return(res)
}


RunGsva <- function(mat, gmt,method="gsva",rnaseq=F,ssgsea.norm=T,abs.ranking=F){
  # rnaseq=T # ONLY FOR RAW COUNT INPUT
  test = gsva(mat,gmt[[1]],method=method,rnaseq=rnaseq,ssgsea.norm=ssgsea.norm,abs.ranking=abs.ranking)
  #save_gsva(test,out,method)
}

SaveGsva <- function(gsva_out, file.name,method){
  if(method=="gsva"){
    es.obs = as.data.frame(gsva_out$es.obs)
  }
  else if(method=="ssgsea"){
    es.obs = as.data.frame(gsva_out)
  }
  es.obs$GeneSet = rownames(es.obs)
  es.obs = es.obs[,c(ncol(es.obs),1:(ncol(es.obs)-1))]
  write.table(es.obs,file=file.name,sep=",",quote=F,col.names=T,row.names=F)
}


GetDefaultGMT <- function(gmt)
{
  if(gmt == 'C7')res = GSA.read.gmt('~/Documents/Work/analysis/GeneSets/c7.all.v6.0.symbols.gmt')
  if(gmt == 'C5')res = GSA.read.gmt('~/Documents/Work/analysis/GeneSets/c5.all.v6.0.symbols.gmt')
  if(gmt == 'H')res = GSA.read.gmt('~/Documents/Work/analysis/GeneSets/h.all.v6.0.symbols.gmt')
  names(res[[1]]) = res[[2]]
  return(res)
}

BuildIndexes <- function(genes,spec='human',msigdb.gsets='all')
{
  
  
  gs.annots = buildIdx(entrezIDs = genes, species = spec, 
                       msigdb.gsets = msigdb.gsets, go.part = TRUE)
  
  return(gs.annots)
}

