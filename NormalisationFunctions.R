library(limma)
library(NanoStringNorm)
library(edgeR)
library(rgl)
library(ggplot2)
library(ggrepel)
SaveNormResult <- function(data, file.name, dir){
  #Clean Sample names
  write.csv(data,paste(dir,file.name,'.csv',sep = ''))
  
}

DoPlots <-function(data,annot, color.cols,spec=NULL, samp.prefix, dir.prefix, dir, raw.data,
                   mds.genes = 100, batch=NULL, legend = T, plot.nobatch=T, plot.3d = F) {
  
  pdf(paste(dir,dir.prefix,'.pdf',sep = ''))
  ####Before and After norms####
  par(mfrow= c(2,2))
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  plot(density(log2(raw.data+0.5)),main = paste(dir.prefix,'Raw',sep = '-'))
  plot(density(data),main = paste(dir.prefix,'Norm',sep = '-'))
  ####Boxplots#####
  boxplot(log2(raw.data+0.5), xlab="", ylab="Log2 Intensities",las=2,main="Unnormalised logIntensities")
  ## Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(log2(raw.data+0.5)),col="blue")
  col.box = 'gold'
  if('BATCH'%in%colnames(annot)){
    col.box = rainbow(max(annot$BATCH))[annot$BATCH]
  }
  
  boxplot(data, xlab="", ylab="Log2 Intensities",las=2,main="Normalised logIntensities", col = col.box)
  ## Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(data),col="blue")
  ###Norm Plots####
  meanSdPlot(data.matrix(data))
  
  #par(mfrow=c(1,1))
  for(i in color.cols){
    t1 = unlist(i)[1]
    t2 = unlist(i)[2]
    if(is.na(t2))
    {
      t2 = t1
    }
    annot$color = '#FFFFFFFF'
    colors =  data.frame(color=rainbow(length(levels(as.factor(annot[,t1]))),start = 0,end = 1),val=levels(as.factor(annot[,t1])))
    #colors =  data.frame(color=colorRampPalette(brewer.pal(12,'Paired'))(length(levels(as.factor(annot[,t1])))),val=levels(as.factor(annot[,t1])))
    # colorP =  colorRampPalette(brewer.pal(12,'Paired'))(200)
    # 
    # colors =  data.frame(color=sample(colorP,length(levels(as.factor(annot[,t1])))),val=levels(as.factor(annot[,t1])))
    # 
    #colorRampPalette(c('white','black'))(5)
    for (j in levels(as.factor(annot[,t1]))){
      
      annot[which(annot[,i]==j),'color'] = as.vector(colors[which(colors$val==j),'color'])
      
    }
    
    ###MDS Plots###
    m = matrix(c(1,1,2,2), 2, 2, byrow = TRUE)
    layout(mat = m, heights = c(0.6,0.4))
    plotMDS(x=data, top=mds.genes,main=paste(dir.prefix,i,sep = '-'),plot = T,
            pch = 16,
            col=annot[match(colnames(data),annot$Sample.ID),'color']
            ,gene.selection="pairwise");
    if(legend == T){
      par(mar=c(0,0,0,0))
      plot(1, type='n',axes=F,xlab = '', ylab='')
      legend('top',inset = 0,legend = as.vector(colors[,'val']),fill = as.vector(colors[,'color']),title = t1,ncol = 4 )
      
    }
    
    if(plot.3d ==T){
      
      dist =  plotMDS(data,top = mds.genes, ndim = 3)$cmdscale.out
      r3dDefaults$windowRect = c(0,50,1000,1000)
      plot3d(x=dist[,1],y = dist[,2], z= dist[,3],col = annot[match(colnames(data),annot$Sample.ID),'color'],
             type = 's', size = 3, lwd = 4, box =T, colvar=NA, xlab = 'x',ylab = 'y', zlab = 'z',radius = 0.02,aspect = T)
      #text3d(x=dist[,1],y = dist[,2], z= dist[,3]+0.04, annot[match(colnames(data),annot$Sample.ID),t2], 
      #       fontweight = 'bold',cex=1, col = annot[match(colnames(data),annot$Sample.ID),'color'])
      #s= spin3d(axis = c(0,0,1),rpm = 5)
      #play3d(s,duration = 10)
      #movie3d(s, duration = 10, dir = dir, clean = F, convert = T)
      legend3d('bottomright',legend = as.vector(colors[,'val']),fill = as.vector(colors[,'color']),
               title = t1, pch = 8 ,cex = 0.7)
      rgl.postscript(paste(dir,dir.prefix,'3D.pdf',sep = ''),fmt = 'pdf')
    }
    
    if(plot.nobatch==T){
      if(!is.null(batch)){
        no.batch = removeBatchEffect(data,annot[,batch])
        
        plotMDS(no.batch, top=mds.genes,main=paste(dir.prefix,i,'No Batch',sep = '-'),
                labels = annot[match(colnames(data),annot$Sample.ID),t2],
                col=annot[match(colnames(data),annot$Sample.ID),'color']
                ,gene.selection="pairwise", legend = T);
        if(legend == T){
          legend('bottomright',levels(as.factor(annot[,t1])),fill = levels(as.factor(annot[,'color'])),title = t1 )
          
        }
      }
    }
  }
  bckp.a = annot
  bckp.d = data
  for(s in spec)
  {
    for(j in levels(as.factor(annot[,s])))
    {
      annot = bckp.a
      data = bckp.d
      annot = as.data.frame(as.matrix(annot[which(annot[,s]==j),]))
      data = data[,as.vector(annot$Sample.ID)]
      if(dim(annot)[1]>2)
      {
        
        for(i in color.cols)
        {
          t1 = unlist(i)[1]
          t2 = unlist(i)[2]
          annot$color = 0
          col = 1
          for (j in levels(as.factor(annot[,t1]))){
            annot[which(annot[,i]==j),'color'] = col
            col = col+1
          }
          ###MDS Plots###
          plotMDS(data, top=mds.genes,main=paste(dir.prefix,i,'Spec',s,j,sep = '-'),
                  labels = annot[,t2],
                  col=annot[match(colnames(data),annot$Sample.ID),'color']
                  ,gene.selection="pairwise");
          if(legend == T){
            legend('bottomright',levels(as.factor(annot[,t1])),fill = levels(as.factor(annot[,'color'])),title = t1 )
          
          }
          if(!is.null(batch)&&s!=batch&&plot.nobatch==T)
          {
          
            no.batch = removeBatchEffect(data,annot[,batch])
            
            plotMDS(no.batch, top=mds.genes,main=paste(dir.prefix,i,'Spec',s,'No Batch',sep = '-'),
                    labels = annot[,t2],
                    col=annot[match(colnames(data),annot$Sample.ID),'color']
                    ,gene.selection="pairwise", legend = T);
            if(legend == T){
              legend('bottomright',levels(as.factor(annot[,t1])),fill = levels(as.factor(annot[,'color'])) ,title = t1)
            }
          }
          
        }
      }
    }
  }
  
  
  dev.off()
  
}
DoNSNPlots <- function(data,annot, samp.prefix, dir.prefix, dir){
  #Clean Sample names

  pdf(paste(dir,dir.prefix,'.pdf',sep = ''))
  Plot.NanoStringNorm(data,plot.type = 'all', label.best.guess = T, label.as.legend = T)
  dev.off()
  
  
}
DoCPMNorm <- function(raw.counts){
  
  ###Filter low read genes####
  L = min(colSums(raw.counts)) # min library size
  P = round(ncol(raw.counts)*0.20) #20% population
  dge = DGEList(counts = raw.counts)
  keep <- rowSums(cpm(dge) > 5/L*1e6) > P
  dge = DGEList(dge[keep,,keep.lib.sizes=F]) 
  ###########################
  ###Normalise####
  dge = calcNormFactors(dge)
  logCPM = cpm(dge,prior.count = 3, log = T)
  ###############
  
  return(logCPM)
  
}

NormaliseVoom <- function(dge.data = dge.data, design=NULL, quanseq =F, pop.pctg = 0.60) {
  if (is.null(dge.data)) {
    return()
  }
  ###Filter low read genes####
  L = min(colSums(dge.data$counts)) # min library size
  P = round(ncol(dge.data$counts)*pop.pctg) #% population
  print(paste('Keep gene if low in less than',P,'samples'))
  keep <- rowSums(cpm(dge.data) > 5/L*1e6) > P
  
  if (quanseq == T){
    keep <- rowSums(cpm(dge.data) > 5) > P
  }
  dge = DGEList(dge.data[keep,,keep.lib.sizes=F]) 
  ##Normalise
  dge = calcNormFactors(dge)
  print('Normalising')
  v <- voom(dge, design = design, plot = T)
  return(v)
  
}


NormaliseNSN <- function(raw.data, norm.method.string = NULL, plot = F, plot.dir = NULL, endo.gene.matrix = T) {
  if (is.null(norm.method.string)) {
      #####Need to identify best method
      print("No Norm Method selected")
      return()
  }
  split.string = strsplit(norm.method.string,'_')
  res <- NanoStringNorm(raw.data, CodeCount = split.string[[1]][1] , Background = split.string[[1]][2], 
                 SampleContent = split.string[[1]][3] , OtherNorm = split.string[[1]][4] , 
                 return.matrix.of.endogenous.probes = endo.gene.matrix, take.log = T )

  return(res)
}

DoNormComp <- function(raw.data, reps) {
  norm.comp.results = norm.comp(
    raw.data, 
    NULL, 
    replicates = reps, 
    CodeCount.methods = c('none', 'sum', 'geo.mean'),
    Background.methods = c('none','mean', 'mean.2sd','max'), 
    SampleContent.methods = c('none','housekeeping.sum', 'housekeeping.geo.mean', 
                              'total.sum','top.mean', 'top.geo.mean', 'low.cv.geo.mean'),
    OtherNorm.methods = c('none','quantile','zscore', 'rank.normal', 'vsn'),
    histogram = FALSE, 
    verbose = TRUE, 
    icc.method = "mixed")
  
  return(norm.comp.results)
  
}

PrintVolcano <- function(dat,x.col,y.col,name.col,plot.title ='Volcano ',labels = NULL,
             x.title = 'log Fold Change',y.title = 'FDR',y.cut = 0.05,x.cut = .58, font.size=.5) {
  
  dat[,'col'] = 'black'
  dat[which(abs(dat[,x.col])>=x.cut),'col'] = 'blue'
  dat[which(dat[,y.col]<y.cut),'col'] = 'red'
  dat[which(dat[,y.col]<y.cut&abs(dat[,x.col])>=x.cut),'col'] = 'green'
  dat[which(dat[,y.col]>=y.cut|abs(dat[,x.col])<x.cut),name.col] = ""
  
  
  p <- ggplot(dat,aes(dat[,x.col],-log10(dat[,y.col])))
  print(
    p + geom_point(aes(dat[,x.col],-log10(dat[,y.col]),color=dat[,'col']),show.legend = F) +  
      coord_cartesian(xlim=c(-1,1)) + labs(x=x.title,y=y.title) +
      scale_color_manual(values = c('black','orange','green','red')) +
      #  theme(legend.title = element_text(size=18), legend.text = element_text(size=18)) #+
       geom_text_repel(aes(dat[,x.col],-log10(dat[,y.col]),label=dat[,name.col] ),
                         size =3, fontface='bold',show.legend = F,
                         segment.colour = "grey50")
  )
  
  # with(dat,plot(eval(x.col),-log10(eval(y.col)),pch=20,main=plot.title,
  #               xlab= x.title, ylab=y.title))
  # #Add colored dots
  # with(subset(dat,eval(y.col)<y.cut),points(eval(x.col),-log10(eval(y.col)),pch = 20,col='red'))
  # with(subset(dat,abs(eval(x.col))>=x.cut),points(eval(x.col),-log10(eval(y.col)),pch = 20,col='blue'))
  # with(subset(dat,abs(eval(x.col))>=x.cut&eval(y.col)<y.cut),
  #      points(eval(x.col),-log10(eval(y.col)),pch = 20,col='green'))
  # 
  # library(calibrate)
  # if(is.null(labels)){
  #   with(subset(dat,abs(eval(x.col))>=x.cut&eval(y.col)<y.cut),
  #        textxy(eval(x.col),-log10(eval(y.col)),labs = eval(name.col),offset = 0,cex=font.size))
  #   
  # }
  # else{
  #   with(subset(dat,eval(name.col)%in%labels),
  #        textxy(eval(x.col),-log10(eval(y.col)),labs = eval(name.col),offset = 0,cex=font.size))
  # }
  # 
}

