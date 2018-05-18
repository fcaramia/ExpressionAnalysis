library(gplots)
library(RColorBrewer)
library(calibrate)

PrintDiffExpTable <- function(contrast, fit, dir, contrast.name, gene.annot){
  top.t = topTable(fit,coef = contrast.name, sort.by = 'p', number = Inf, p.value = 1,adjust.method = 'BH')
  top.t$ENTREZ.ID = as.integer(rownames(top.t))
  top.t$GENE.SYMBOL = gene.annot[match(rownames(top.t),gene.annot[,'Entrez.Gene.ID']),'Approved.Symbol']
  top.t$CHROMOSOME = gene.annot[match(rownames(top.t),gene.annot[,'Entrez.Gene.ID']),'Chrm']
  
  write.csv(top.t, paste(dir,contrast.name,'.csv',sep = ''), row.names = F)
  return(top.t)
}

PrintDiffExpTable2 <- function(contrast, fit, dir, contrast.name, gene.annot){
  top.t = topTable(fit,coef = contrast.name, sort.by = 'p', number = Inf, p.value = 1,adjust.method = 'BH')
  top.t$SYMBOL = rownames(top.t)
  top.t$SYMBOL.EXCEL.VIEW = paste('`',rownames(top.t),sep = '')
  top.t$CHROMOSOME = gene.annot[match(rownames(top.t),gene.annot[,'Approved.Symbol']),'Chrm']
  
  write.csv(top.t, paste(dir,contrast.name,'.csv',sep = ''), row.names = F)
  return(top.t)
}


mydistfun = function(d) {
  return (as.dist((1-cor(t(d))/2)))  # correlation, good for expression data
}

myhclustfun = function(d) {
  return (hclust(d,"average"))
}

createHeatMap = function(dat.ext,plot.lbl, labels, legend = F, batch=NULL, hclustering = T, font.size = .5,labels.as.cols =F) {
  if(!is.null(batch)){
    dat.ext = removeBatchEffect(dat.ext,batch)
  }
  if(nrow(dat.ext)<2)
  {
    return("")
  } 
  
  means=apply(dat.ext,1,mean)
  
  g.clust = as.dendrogram(myhclustfun(mydistfun(as.matrix(dat.ext))))
  samp.clust = as.dendrogram(myhclustfun(mydistfun(as.matrix(t(dat.ext)))))    

  if(!hclustering){
    samp.clust = NULL
  }
  
  x=as.matrix(sweep(dat.ext,1,means))
  if(labels.as.cols==T){
    colnames(x) = labels
  }
  breaks.up=quantile(x,0.98)
  breaks.down=quantile(x,0.02)
  breaks=seq(breaks.down,breaks.up,(breaks.up-breaks.down)/75)
  
  known.labels.factor=factor(labels)
  #lw = c(0.1,4)
  #### modified a bit to get colours for legends
  labels.cols.all = sample(colours(),length(known.labels.factor))
  label.cols = labels.cols.all[known.labels.factor]
  col1 = colorRampPalette(c("yellow",'black','blue'))
  #lmat.layout = rbind(4:3,2:1)
  
  heatmap.2(x,dendrogram="column",
            Rowv=g.clust,
            Colv=samp.clust,
            breaks=breaks,ColSideColors=label.cols,
            distfun=NULL,hclustfun=NULL,scale="none",col=col1,srtCol = 15,
            trace="none",cexCol=1,density.info="none",#margins=c(15,12),lwid= lw,
            main=plot.lbl, key = T,  cexRow = font.size)
  if (legend == T)legend(x=1,y=0,legend=levels(known.labels.factor),fill=labels.cols.all,cex=0.5,xjust=1,xpd=NA)
  
}

PrintHeatMap <- function(mat, contrast, title, annot, top.table, dir, cont.col, font.size = 0.5, legend = F, hclust = T, labels.as.cols = F){
  #######Set data
  top.genes = unique(rownames(top.table))
  top.genes = top.genes[which(top.genes%in%rownames(mat))]
  cont.groups = unlist(strsplit(contrast,'vs'))
  if('All'%in%cont.groups){
    samples = as.vector( annot[,'Sample.ID'])
    labels = as.vector(annot[,cont.col]  )
  }else{
    samples = as.vector(annot[which(annot[,cont.col]%in%cont.groups),'Sample.ID']) 
    labels = annot[which(annot[,cont.col]%in%cont.groups),cont.col]
  }
  top.data = mat[top.genes,samples]
  
  createHeatMap(top.data,title,labels, font.size = font.size, legend = legend, hclustering = hclust, labels.as.cols = labels.as.cols)
}




