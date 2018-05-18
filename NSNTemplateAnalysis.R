source("~/Documents/Work/analysis/NanoStringPipeline/Normalisation_Functions.R")
source("~/Documents/Work/analysis/NanoStringPipeline/DifferentialExpression.R")
###Variables#### <-----------Action
input.dir = "~/Documents/Work/Projects/Heloise/LentigoMaligna/Combined.Analysis/data/"
#prefix = 'X20160824_160824_Heloise_Pan_Cancer_[[:alnum:]]{2}_'
contrast.file = 'contrast.csv'
output.dir = "~/Documents/Work/Projects/Heloise/LentigoMaligna/Combined.Analysis/ResNormalVSTumorsWVolcanos/"
dir.create(output.dir,showWarnings = F)
###Read Data####

raw.data=read.markup.RCC(
  rcc.path = input.dir,
  rcc.pattern = "*.RCC|*.rcc",
  exclude = NULL,
  include = NULL)

##Read contrast file####
contrast.annot = read.csv(paste(input.dir,contrast.file, sep = ''), as.is = T)
#colnames(raw.data$x)=gsub(prefix,'',colnames(raw.data$x))
#contrast.annot$SAMPLE=gsub(prefix,'',contrast.annot$SAMPLE)
#colnames(raw.data$header)=gsub(prefix,'',colnames(raw.data$header))

###Create dge counts file#####
dge.data.raw = data.matrix(raw.data$x[,-c(1:3)])
rownames(dge.data.raw) = raw.data$x$Name
dge.data <- DGEList(counts=dge.data.raw)


##Get replicates##
reps = NULL
if ("REPLICATE" %in% colnames(contrast.annot)){
  reps = contrast.annot[match(rownames(dge.data$samples),contrast.annot$SAMPLE),'REPLICATE']
}

####Normalise Voom######
####Create Design####### <-------------- Action
#aux <- factor(paste(contrast.annot$PATIENT,contrast.annot$COND,sep="."))
#design <- model.matrix(~0+contrast.annot$SAMPLE_TYPE+contrast.annot$CLEREANCE)
#colnames(design) <- c("Normal","Tumor","Cle")
#rownames(design) <- contrast.annot$SAMPLE
#print (design)
#voom.norm = NormaliseVoom(dge.data = dge.data, design)
####Normalise NSN####### <-------------- Action 
#comp.result = DoNormComp(raw.data, reps)
nsn.norm.method = 'sum_none_housekeeping.geo.mean_vsn'
#nsn.norm.method = 'sum_max_top.geo.mean_none'

nsn.norm = NormaliseNSN(raw.data = raw.data, nsn.norm.method)
#nsn.norm.full = NormaliseNSN(raw.data, nsn.norm.method, endo.gene.matrix = F )
####Plot MDS, Densities, libraries##### 
#DoPlots(voom.norm$E, contrast.annot, color.cols = c('TREATMENT','PATIENT','BATCH', 'INFLAMMATION','CLEREANCE','SAMPLE_TYPE'), 
#        samp.prefix = prefix, dir.prefix = "VOOM", output.dir, dge.data$counts)

#DoPlots(data.matrix(nsn.norm), contrast.annot, color.cols = c('TREATMENT','PATIENT','BATCH', 'INFLAMMATION','CLEREANCE','SAMPLE_TYPE'), 
#        samp.prefix = prefix, dir.prefix = paste("NSN",nsn.norm.method,sep = '.'), output.dir, dge.data$counts)

#DoNSNPlots(nsn.norm.full,contrast.annot,samp.prefix = prefix, dir.prefix = 'NSN.PLOTS', dir = output.dir)
#####Save Normalised results#########
#SaveNormResult(voom.norm$E, samp.prefix = prefix, dir.prefix = "Norm.VOOM",  dir = output.dir)
SaveNormResult(data.matrix(nsn.norm), samp.prefix = prefix, dir.prefix = "Norm.NSN",  dir = output.dir)

#####Choose data###################### <-------------Action
dge.data = DGEList(nsn.norm)

#######Modify Samples to use if necessary############## <-----------Action
samples.to.use = contrast.annot[which(contrast.annot$USE == 'yes'),'SAMPLE']
dge.data = dge.data[,samples.to.use]
contrast.annot = contrast.annot[which(contrast.annot$SAMPLE%in%colnames(dge.data)),]


# #####If technical replicates present. Average###### <-----------Action
# if(!is.null(reps)){
#   dge.data = avearrays(dge.data,ID = reps)
#   #Sample names post average
#   names = c("RP1_TM","RP1_IT","RP1_NT","RP2_TM","RP2_IT","RP2_NT","RP3_TM","RP3_IT","RP3_NT")
#   colnames(dge.data) = names
#   contrast.annot$SAMPLE=gsub('(.*)_(.*)_.*','\\1_\\2',contrast.annot$SAMPLE)
#   contrast.annot = contrast.annot[which(contrast.annot$SAMPLE%in%colnames(dge.data)),]
# }



######Define the design #############################  <------------Action 
contrast.annot$TYPE.CLE = paste(contrast.annot$SAMPLE_TYPE,contrast.annot$CLEREANCE, sep = '.')

DoPlots(raw.data = dge.data.raw ,annot =  contrast.annot, color.cols = c('TYPE.CLE'), 
        samp.prefix = prefix, dir.prefix = "FILTER", dir = output.dir,data =dge.data$counts)



design <- model.matrix(~0+contrast.annot$TYPE.CLE)
colnames(design) = levels(factor(contrast.annot$TYPE.CLE))
rownames(design) = contrast.annot$SAMPLE
#####When samples from same patient do correlation step ####### <---------Action
corfit <- duplicateCorrelation(dge.data$counts,design = design, block = contrast.annot$PATIENT)
fit <- lmFit(dge.data$counts,design = design, block = contrast.annot$PATIENT, correlation = corfit$consensus.correlation)

#####Make Contrasts############## <-----------Action
cm <- makeContrasts(
        NormalvsTumorCle = ((Normal.cle+Normal.no_cle)/2 - Tumor.cle),
        NormalvsTumorNoCle = ((Normal.cle+Normal.no_cle)/2 - Tumor.no_cle),
        NormalClevsTumorCle = (Normal.cle - Tumor.cle),
        NormalNoClevsTumorNoCle = (Normal.no_cle - Tumor.no_cle),
        levels = design
)
cms = c('NormalvsTumorCle',"NormalvsTumorNoCle","NormalClevsTumorCle","NormalNoClevsTumorNoCle")
#cms=c('NormalvsTumorCle')
fit2 = contrasts.fit(fit,cm)
fit2 = eBayes(fit2)

contrast.column = 'TYPE.CLE'

for (i in cms){
  ####Print Table for each contrast####
  top.t = PrintDiffExpTable(contrast = i, fit = fit2, dir = output.dir)
  top.t$Symbol = rownames(top.t)
  
  ####Print Heatmaps######
  top.genes=top.t[which(top.t$P.Value<0.05),]
  pdf(paste(output.dir,i,'.pdf',sep = ''))
  if(nrow(top.genes)>10){
    #Print Vocano pval
    PrintVolcano(top.t,expression(logFC),expression(P.Value),expression(Symbol),'Volcano pval'
                 ,x.title = 'log Fold Change',y.title = 'p-value',pval.val = .05,fc.val = .18)
    top10.genes=top.genes[1:10,]
    PrintHeatMap(mat = dge.data$counts, design = design,contrast.values =  cm[,i],contrast.name = i,
                 contrast.column = contrast.column,title = paste('All p.value<0.05\n',i),annot =contrast.annot ,
                 top.table = top.genes,dir = output.dir)
    PrintHeatMap(mat = dge.data$counts, design = design,contrast.values =  cm[,i],contrast.name = i,
                 contrast.column = contrast.column,title = paste('Top 10 p.value<0.05\n',i),annot =contrast.annot ,
                 top.table = top10.genes,dir = output.dir)
    if (nrow(top.genes)>=50){
      top50.genes=top.genes[1:50,]
      PrintHeatMap(mat = dge.data$counts, design = design,contrast.values =  cm[,i],contrast.name = i,
                   contrast.column = contrast.column,title = paste('Top 50 p.value<0.05\n',i),annot =contrast.annot ,
                   top.table = top50.genes,dir = output.dir)

    }

  }
  top.genes=top.t[which(top.t$adj.P.Val<0.05),]
  if(nrow(top.genes)>10){
    #Print Volcano adj.pval
    PrintVolcano(top.t,expression(logFC),expression(adj.P.Val),expression(Symbol),'Volcano FDR'
                 ,x.title = 'log Fold Change',y.title = 'FDR',pval.val = .05,fc.val = .18)
    top10.genes=top.genes[1:10,]
    PrintHeatMap(mat = dge.data$counts, design = design,contrast.values =  cm[,i],contrast.name = i,
                 contrast.column = contrast.column,title = paste('All adj.p.value<0.05\n',i),annot =contrast.annot ,
                 top.table = top.genes,dir = output.dir)
    PrintHeatMap(mat = dge.data$counts, design = design,contrast.values =  cm[,i],contrast.name = i,
                 contrast.column = contrast.column,title = paste('Top 10 adj.p.value<0.05\n',i),annot =contrast.annot ,
                 top.table = top10.genes,dir = output.dir)
    if (nrow(top.genes)>=50){
      top50.genes=top.genes[1:50,]
      PrintHeatMap(mat = dge.data$counts, design = design,contrast.values =  cm[,i],contrast.name = i,
                   contrast.column = contrast.column,title = paste('Top 50 adj.p.value<0.05\n',i),annot =contrast.annot ,
                   top.table = top50.genes,dir = output.dir)


    }

  }
  dev.off()
}
