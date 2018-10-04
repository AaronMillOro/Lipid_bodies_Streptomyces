library("ade4")
library("gdata")
library("lme4")
library("nlme")
library("car")
library("gplots")
library("gdata")
library("made4")
library("clValid")
library("lattice")
library("ggplot2")
library("reshape2")


setwd("~/Documents/2016_02_17 Corps lipidiques MJ/R_analysis")

#******	 Spectral Counting analysis of proteins from lipid bodies  ****
# import data obtained from X!TandemPipeline

	#emPAI<- read.csv2(file="emPAI_lipid_bodies.csv", header = TRUE, sep="\t", dec=",", stringsAsFactors=FALSE) #443 prots
	dataSC<- read.csv2(file="SC_lipid_bodies_norm.csv", header = TRUE, sep="\t", dec=",", stringsAsFactors=FALSE) #443 prots
	voies<- read.csv2(file="paths.tsv", header = TRUE, sep=",", stringsAsFactors=FALSE)

# generate metadata	
	r<- (colnames(dataSC[6:17]))
	a=strsplit(unique(r), "_", fixed=TRUE)
	esp=NULL
	rep=NULL
	temps=NULL
	for(i in 1:length(a)){
	  esp=c(esp, a[[i]][1])
	  temps=c(temps, a[[i]][2])
	  rep=c(rep, a[[i]][3])
	  print(i)
	}
	metadata=cbind.data.frame(msrunfile=unique(r), esp=esp, temps=temps, rep=rep)
	metadata$esp.temps=as.factor(paste(metadata$esp,metadata$temps,sep='-'))
	
# Generate a dataframe containing all the informations
	test<- stack(dataSC[,6:17])
	names(test)[1] <- "spectra"
	names(test)[2] <- "msrunfile"
	test=merge(test, metadata, "msrunfile")
	
	tab.sc<-cbind.data.frame(test, protein=rep(dataSC$Top.Protein.ID, 12))
	tab.sc<-cbind.data.frame(tab.sc, desc=rep(dataSC$Top.Protein.Description, 12))
	head(tab.sc)
	
# Filter proteins showing low ration between conditions
	
	drop.low.ratio=data.frame(dataSC$Top.Protein.ID)
	names(drop.low.ratio) <- sub("dataSC.Top.Protein.ID", "prot", names(drop.low.ratio))
	
	proteines = levels(tab.sc$protein) # 2355 prots
	min.ratio = 2
	for (i in 1:length(proteines)){ 
	  low.ratio=tab.sc[tab.sc$protein==proteines[i],] 
	  low.ratio=drop.levels(low.ratio) 
	  tab.ratios = aggregate(low.ratio$spectra, list(low.ratio$protein, low.ratio$esp.temps), FUN = mean)
	  maxvalue = max(tab.ratios$x)
	  minvalue = min(tab.ratios$x) 
	  ratio = maxvalue/minvalue
	  if (ratio == Inf)
	    ratio = maxvalue/(minvalue+1)
	  if (ratio >= min.ratio)
	    drop.low.ratio$ratio[i] <-2
	  else
	    drop.low.ratio$ratio[i] <-0
	  print(i) 
	}	  
	
	good_spectra=drop.low.ratio$prot[drop.low.ratio$ratio>1]
	good_spectra=drop.levels(good_spectra)
	SC=tab.sc[tab.sc$protein %in% good_spectra,]
	SC = drop.levels(SC) 
	str(SC)
	
	###################
	###  384 prots ###
	###################
	
# GLM model and multiple ANOVA tests
	proteines = levels(SC$protein)
	resultglm = NULL 
	for (i in 1:length(proteines))
	  {
	  sub=SC[SC$protein==proteines[i],] 
	  sub=drop.levels(sub) 
	  model=glm(spectra~esp+temps+rep, family="quasipoisson", data=sub) 
	  test=anova(model, test="Chisq") 
	  resultglm=rbind.data.frame(resultglm, cbind.data.frame(prot=proteines[i], pesp=test[[5]][2], ptemps=test[[5]][3], prep=test[[5]][4])) 
	  print(i) 
	}	  
	 
	resultglm$fdr.esp=p.adjust(resultglm$pesp,method="fdr")
	resultglm$fdr.temps=p.adjust(resultglm$ptemps,method="fdr")

	signif.esp=resultglm$prot[resultglm$fdr.esp<0.01]
	signif.esp = drop.levels(signif.esp)
	signif.temps=resultglm$prot[resultglm$temps<0.01]
	signif.temps = drop.levels(signif.temps)
	liste_prot_signif = union(signif.esp, signif.temps)
	length(liste_prot_signif)
	
	spectral.count.glm.signif = SC[which(SC$protein %in% liste_prot_signif),]
	spectral.count.glm.signif = drop.levels(spectral.count.glm.signif)
	length(unique(spectral.count.glm.signif$protein)) ## 200 prots
  spectral.count.glm.signif_INTACT <- spectral.count.glm.signif
  

  #levels(spectral.count.glm.signif$esp)[levels(spectral.count.glm.signif$esp)=="Coel"] <- "M145"
  #levels(spectral.count.glm.signif$esp)[levels(spectral.count.glm.signif$esp)=="Livi"] <- "TK24"
  #spectral.count.glm.signif$esp.temps <- as.factor (paste(spectral.count.glm.signif$esp, spectral.count.glm.signif$temps, sep="-" ))
	
	################# Export data Spectral Count

	spec_signif=tapply(spectral.count.glm.signif$spectra,list(spectral.count.glm.signif$protein,spectral.count.glm.signif$esp.temps),FUN=mean)
	spec_signif = as.data.frame(spec_signif)
	# by.y =0 parce que il n'a pas le nom des proteines, donc utilise le nom des lignes
	test1 =merge(resultglm, spec_signif, by.x="prot", by.y=0)
  test1 = merge (test1, voies [,-2], by=c("prot"), all.x=TRUE)
	write.table(test1,"prots_signif.tsv",sep="\t",row.names=F,col.names=T)
	
	colnames(spectral.count.glm.signif)[7]="prot"
	spectral.count.glm.signif=merge(spectral.count.glm.signif,voies [,-c(2)], by=c("prot"), all.x=TRUE)
	spectral.count.glm.signif <- drop.levels(spectral.count.glm.signif)
	spectral.count.glm.signif$Names <- as.factor(spectral.count.glm.signif$Names)
	spectral.count.glm.signif$Sub_class <- as.factor(spectral.count.glm.signif$Sub_class)

# boxplots	
	formule1=formula("spectra ~ esp")
	formule3=formula("spectra~esp.temps")
	
	pdf(file="lipid_bodies_boxplots_signif.pdf", width=10,height=6)
	for (i in 1:length(unique(spectral.count.glm.signif$prot))) {
	  subSC= spectral.count.glm.signif[spectral.count.glm.signif$prot==levels(spectral.count.glm.signif$prot)[i],]
	  par(mfrow=c(1,2))
	  boxplot(formule1,subSC,las=2,col=c("blue","red", 'darkgreen'),main=unique(subSC$Names),ylab="Spectral Count")
	  boxplot(formule3,subSC,las=2,col=c("blue","blue","red","red","darkgreen","darkgreen"), main=unique(subSC$Sub_class),ylab="Spectral Count")
	}
	dev.off()
	
	tab.acpSC = tapply (spectral.count.glm.signif$spectra, list(spectral.count.glm.signif$prot,spectral.count.glm.signif$msrunfile), FUN=mean)

	####### Principal components analysis  + heatmap
	
	quanti_data_acp = na.omit(tab.acpSC)
	quanti_data_acp = t(quanti_data_acp)
	z <- dudi.pca(quanti_data_acp,center = T, scale = T, scannf = F, nf = 4)
	sm = sum(z$ei)
	pound = round((z$e/sm*100),digits = 1)
	acp = z$li
	acp$msrunfile=row.names(acp)
	acp = merge(acp,metadata,by= c("msrunfile"),all.x=TRUE,all.y=FALSE)
	
	pdf(file="ACP_signif.pdf", width = 13, height = 7)
	par(mfrow=c(1,2))
	plot(acp$Axis1, acp$Axis2,type="n", xlab=paste("Axe1(",pound[1],"%)",sep=" "),ylab=paste("Axe2(",pound[2],"%)",sep=" "))
	text(acp$Axis1, acp$Axis2, acp$msrunfile, col=c(acp$esp), cex = 0.9)
	abline(h=0, v=0)
  plot(acp$Axis1,acp$Axis3,type="n", xlab=paste("Axe1(",pound[1],"%)",sep=" "),ylab=paste("Axe3(",pound[3],"%)",sep=" "))
  text(acp$Axis1, acp$Axis3, acp$msrunfile, col=c(acp$esp), cex=0.9)
  abline(h=0, v=0)
	dev.off()
  
	# Class vec contains the values ib Spectral counting plus the metabolic pathways
	heatmap_SC<- read.csv2("class_vec.csv", header = TRUE, sep="\t", dec=",")
	pdf(file="Heatplot_signif.pdf", width = 7, height = 23)
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$CARBON, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Carbon")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$CELL_DIVISION_WALL, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Cell division/wall")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$NITROGEN, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Nitrogen")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$ENERGY, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Respiratory chain")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$METABOLITES, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Secondary metabolites")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$TRANSPORT, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Transport")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$TRANSLA_TRANSCRIP, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Translation/transcription")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$SIGNAL, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Signaling")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$OTHERS, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Other")
  heatplot(heatmap_SC[,2:7], margins=c(5,20) ,distfun="euclidean", dend="row", cexRow= 0.6, cex=0.7, classvec=heatmap_SC$UNKNOWN, classvecCol=c("white","darkgreen"),labRow=heatmap_SC$Names, main="Unknown")
  dev.off()
