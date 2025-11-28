library(vcfR)
library(adegenet)
library(StAMPP)

#genlight opject batis
x_files <- list.files(pattern="batis_nothin.vcf")
for(a in 1:length(x_files)) {
	# read vcf
	a_rep <- read.vcfR(x_files[a])
	
	# convert to genlight
	a_rep <- vcfR2genlight(a_rep)}

saveRDS(a_rep, file="a_rep.rds")
readRDS("a_rep.rds")


#genlight opject cossypha
x_files <- list.files(pattern="cossypha_nothin.vcf")
for(a in 1:length(x_files)) {
	# read vcf
	a_rep <- read.vcfR(x_files[a])
	
	# convert to genlight
	a_rep <- vcfR2genlight(a_rep)}

saveRDS(a_rep, file="cossypha_rep.rds")
readRDS("cossypha_rep.rds")


#genlight opject chamaetylas
x_files <- list.files(pattern="chamaetylas_nothin.vcf")
for(a in 1:length(x_files)) {
	# read vcf
	a_rep <- read.vcfR(x_files[a])
	
	# convert to genlight
	a_rep <- vcfR2genlight(a_rep)}

saveRDS(a_rep, file="chamaetylas_rep.rds")
readRDS("chamaetylas_rep.rds")

#genlight opject phylloscopus
x_files <- list.files(pattern="phylloscopus_nothin.vcf")
for(a in 1:length(x_files)) {
	# read vcf
	a_rep <- read.vcfR(x_files[a])
	
	# convert to genlight
	a_rep <- vcfR2genlight(a_rep)}

saveRDS(a_rep, file="phylloscopus_rep.rds")
readRDS("phylloscopus_rep.rds")

#Batis
source("reich_fst.r")
gl = readRDS("a_rep.rds")
batis_fst <- reich.fst(gl = gl)
batis_fst$fsts
write.csv(batis_fst$fsts,"/common/cooperlab/phuss58/albertine/12_fst/batis_fst.csv",quote=F)

#Cossypha 
source("reich_fst.r")
gl = readRDS("cossypha_rep.rds")
cossypha_fst <- reich.fst(gl = gl)
cossypha_fst$fsts
write.csv(cossypha_fst$fsts,"/common/cooperlab/phuss58/albertine/12_fst/cossypha_fst.csv",quote=F)

#Chamaetylas
source("reich_fst.r")
gl = readRDS("chamaetylas_rep.rds")
chamaetylas_fst <- reich.fst(gl = gl)
chamaetylas_fst$fsts
write.csv(chamaetylas_fst$fsts,"/common/cooperlab/phuss58/albertine/12_fst/chamaetylas_fst.csv",quote=F)

#Phylloscopus
source("reich_fst.r")
gl = readRDS("phylloscopus_rep.rds")
phylloscopus_fst <- reich.fst(gl = gl)
phylloscopus_fst$fsts
write.csv(phylloscopus_fst$fsts,"/common/cooperlab/phuss58/albertine/12_fst/phylloscopus_fst.csv",quote=F)

## reich fst estimator
## vectorized version
## input=genlight object
## FST will be calculated between pops in genlight object
## specify number of bootstraps using "bootstrap=100"

reich.fst <- function(gl, bootstrap=FALSE, plot=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    install.packages("dplyr")
    library(dplyr, character.only=T)
  }
  if (!require("dartR",character.only=T, quietly=T)) {
    install.packages("dartR")
    library(dartR, character.only=T)
  }
  if (!require("combinat",character.only=T, quietly=T)) {
    install.packages("combinat")
  }
  
  nloc <- gl@n.loc
  npop <- length(gl@ind.names)

  fsts <- matrix(nrow=npop,
                 ncol=npop,
                 dimnames=list(gl@ind.names,gl@ind.names))
  
  if (bootstrap != FALSE){
    n.bs <- bootstrap
    bs <- data.frame(matrix(nrow=nrow(combinat::combn2(gl@ind.names)),
                            ncol=n.bs+5))
  }
  
  k <- 0
  
  for (p1 in gl@ind.names){
    for (p2 in gl@ind.names){
      if (which (gl@ind.names == p1) < (which gl@ind.names == p2)) {
        k <- 1+k
        
        pop1 <- gl.keep.ind(gl, p1, mono.rm=FALSE, v=0)
        pop2 <- gl.keep.ind(gl, p2, mono.rm=FALSE, v=0)
        
        a1 <- colSums2(as.matrix(pop1),na.rm=T)
        a2 <- colSums2(as.matrix(pop2),na.rm=T)
        n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
        n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
        
        h1 <- (a1*(n1-a1))/(n1*(n1-1))
        h2 <- (a2*(n2-a2))/(n2*(n2-1))
        
        N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
        D <- N + h1 + h2
        
        F <- sum(N, na.rm=T)/sum(D, na.rm=T)
        fsts[p2,p1] <- F
        if (verbose == TRUE) {
        	print(paste("Pop1: ",p1,", Pop2: ",p2,", Reich FST: ",F,sep=""))
        }
        
        if (bootstrap != FALSE) {
          if (verbose == TRUE) {
            print("beginning bootstrapping")
          }
          
          bs[k,1:3] <- c(p2,p1,as.numeric(F))
          
          for (i in 1:n.bs){
            loci <- sample((1:nloc), nloc, replace=TRUE)
          
            pop1.bs <- matrix(as.matrix(pop1)[,loci],
                              ncol=length(loci))
            pop2.bs <- matrix(as.matrix(pop2)[,loci],
                              ncol=length(loci))
          
            a1 <- colSums2(as.matrix(pop1.bs),na.rm=T)
            a2 <- colSums2(as.matrix(pop2.bs),na.rm=T)
            n1 <- apply(as.matrix(pop1.bs),2,function(x) 2*sum(!is.na(x)))
            n2 <- apply(as.matrix(pop2.bs),2,function(x) 2*sum(!is.na(x)))
          
            h1 <- (a1*(n1-a1))/(n1*(n1-1))
            h2 <- (a2*(n2-a2))/(n2*(n2-1))
          
            N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
            D <- N + h1 + h2
          
            F.bs <- sum(N, na.rm=T)/sum(D, na.rm=T)
            bs[k,i+5] <- F.bs
          }
          if (verbose == TRUE){
            print(paste("bootstrapping 95% CI: ",
                        quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),"-",
                        quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T)))
          }
          
          bs[k,4:5] <- c(quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),
                         quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T))
        }
        
      }
    }
  }
  
  fsts[fsts < 0] <- 0
  
  if (bootstrap != FALSE){
    colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","min_CI","max_CI")
    fst.list <- list(fsts,bs)
    names(fst.list) <- c("fsts","bootstraps")
    
    if (plot == TRUE){
      print("drawing plot with bootstraps")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- bs[,1:5]
      plot.data$fst_estimate <- as.numeric(plot.data$fst_estimate)
      plot.data$min_CI <- as.numeric(plot.data$min_CI)
      plot.data$max_CI <- as.numeric(plot.data$max_CI)
      plot.data$pop_pair <- paste(plot.data$pop1,plot.data$pop2,sep="_")
      plot.data$signif <- case_when(plot.data$min_CI > 0 ~ TRUE,
                                    TRUE ~ FALSE)

      
      bs.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate,col=signif)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=min_CI,ymax=max_CI),width=0.1,linewidth=1) + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(bs.plot)
    }
  } else {
    fst.list <- list(fsts)
    names(fst.list) <- "fsts"
    
    if (plot == TRUE){
      print("drawing plot without bootstraps")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- data.frame(combinat::combn2(row.names(fsts)),
                              fst_estimate=fsts[lower.tri(fsts)])
      plot.data$pop_pair <- paste(plot.data$X1,plot.data$X2,sep="_")

      fst.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(fst.plot)
    }
  }
  
  return(fst.list)
  beepr::beep()
}