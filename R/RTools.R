InstallPackages=function(x){
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      # check the package is a bioconductor package and install it if required
      if (length(BiocManager::available(paste0("^",i,"$")))>0){
        BiocManager::install(i)
      }
      else{
        #  If package was not able to be loaded then re-install
        install.packages( i , dependencies = TRUE )
      }
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}
Packagelist=c("optimbase",
              "rmarkdown",
              "flexdashboard",
              "ggnewscale",
              "kableExtra",
              "pheatmap",
              "Glimma",
              "dendextend",
              "data.table", 
              "ggplot2",
              "MatrixEQTL",
              "stringr", 
              "grid", 
              "gridExtra", 
              "limma",
              "sva",
              "plotly", 
              "RColorBrewer",
              "tidyverse",
              "ggpubr",
              "svglite",
              "topGO",
              "R.utils",
              "hash",
              "doParallel",
              "ggbiplot",
              "devtools",
              "clusterProfiler",
              "AnnotationForge",
              "ViSEAGO",
              "ggupset",
              "wesanderson",
              "enrichplot",
              "pcr",
              "ggvenn",
              "rtracklayer",
              "Hmisc",
              "STRINGdb",
              "org.Hs.eg.db",
              "M3C",
              "networkD3")
# InstallPackages(Packagelist)
remotes::install_github("vqv/ggbiplot")

# library(optimbase)
library(data.table)
# library(ggplot2)
library(data.table)
library(stringi)
library(stringr)

# Simple tools
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# Theme of ggplot single graph

mynamestheme <- ggplot2::theme(plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5 ), 
                      legend.title = ggplot2::element_text( family = "Helvetica"), 
                      legend.text = ggplot2::element_text(family = "Helvetica"), 
                      axis.title = ggplot2::element_text(family = "Helvetica", size = (15)),
                      axis.text = ggplot2::element_text(family = "Helvetica", size = (15), color = "black"),
                      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                      panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
mynamesthemeSmall <- ggplot2::theme(plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", size = (12),hjust = 0.5 ), 
                      legend.title = ggplot2::element_text( family = "Helvetica"), 
                      legend.text = ggplot2::element_text(family = "Helvetica"), 
                      axis.title = ggplot2::element_text(family = "Helvetica", size = (10)),
                      axis.text = ggplot2::element_text(family = "Helvetica", size = (10), color = "black"),
                      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                      panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))

ThemeZim <- ggplot2::theme(plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", size = (18),hjust = 0.5 ), 
                      legend.title = ggplot2::element_text( family = "Helvetica"), 
                      legend.text = ggplot2::element_text(family = "Helvetica"), 
                      axis.title = ggplot2::element_text(family = "Helvetica", size = (18)),
                      axis.text = ggplot2::element_text(family = "Helvetica", size = (18), color = "black"),
                      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                      panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))

# Theme of ggplot graph pdf
PDFtheme <- ggplot2::theme(plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", size = (10),hjust = 0.5 ), 
                      legend.title = ggplot2::element_text( family = "Helvetica"), 
                      legend.text = ggplot2::element_text(family = "Helvetica"), 
                      axis.title = ggplot2::element_text(family = "Helvetica", size = (10)),
                      axis.text = ggplot2::element_text(family = "Helvetica", size = (10), color = "black"),
                      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                      panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))

ThemeGO <- ggplot2::theme(plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", size = (15),hjust = 0.5 ), 
                      legend.title = ggplot2::element_text( family = "Helvetica"), 
                      legend.text = ggplot2::element_text(family = "Helvetica"), 
                      axis.title = ggplot2::element_text(family = "Helvetica", size = (15)),
                      axis.text = ggplot2::element_text(family = "Helvetica", size = (15), color = "black"),
                      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                      panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))


# Moving average 
mav <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)}

Radar=function(Mat, Normalisation=T, Annotation=F){
  Nv=ncol(Mat)
  Ni=nrow(Mat)
  # NormMat=t(apply(Mat, 1, rescale))
  
  if (isTRUE(Normalisation)){
    NormMat=t(apply(Mat, 1, function(x)((x-min(x))/sum(x-min(x)))))
    Mat=NormMat
  }
  else{
    NormMat=t(apply(Mat, 1, function(x)((x-min(x)))))
    Mat=NormMat
  }
  
  D=dist(t(Mat))
  IndexVar=order.hclust(hclust(D))
  Mat=Mat[,IndexVar]
  # matrix of value for axis 
  AxisVal1=zeros(1,Nv)
  AxisVal2=ones(1, Nv)*max(Mat)
  AxisVal=rbind(AxisVal1, AxisVal2)
  
  
  # Get angle step for axis
  StepAngles=2*pi/(Nv)
  #Get the range of angles
  Angles=1:Nv*StepAngles
  Mat_Angle=t(matrix(rep(Angles, Ni),nrow = Nv))
  colnames(Mat_Angle)=colnames(Mat)
  Mat_AngleAxis=Mat_Angle[1:2,]
  # create the x and y coordinates
  Xaxis=as.data.table(cos(Mat_AngleAxis)*AxisVal)
  Yaxis=as.data.table(sin(Mat_AngleAxis)*AxisVal)
  
  XMat=as.data.table(cos(Mat_Angle)*Mat)
  YMat=as.data.table(sin(Mat_Angle)*Mat)
  
  # DataTBMat=cbind(melt(XMat),melt(YMat))
  # DataTBMat=DataTBMat[,c(1,2,4)]
  
  # get the centroid of dots
  
  # DataTBMatmeanX=rowMeans(XMat)
  # DataTBMatmeanY=rowMeans(YMat)
  
  DataTBMatmeanX=apply(XMat, 1, mean)
  DataTBMatmeanY=apply(YMat, 1, mean)
  #  factor to reduce axis coordinate
  Factor=max(abs(DataTBMatmeanX))/max(Xaxis)
  Factor2=max(abs(DataTBMatmeanY))/max(Yaxis)
  Xaxis=Xaxis*Factor*1.3
  Yaxis=Yaxis*Factor2*1.3
  # setnames(DataTBMat,c("Var","X", "Y"))
  
  # convert to DataTable
  DataTBAxis=cbind(melt(Xaxis),melt(Yaxis))
  DataTBAxis=DataTBAxis[,c(1,2,4)]
  setnames(DataTBAxis,c("Var","X", "Y"))
  
  if (isFALSE(Annotation)){
    P=ggplot(DataTBAxis, aes(x = X, y=Y, color=Var))+
      geom_point()+
      geom_smooth(method=lm)
    P=P+geom_point(data =as.data.frame(cbind(DataTBMatmeanX,DataTBMatmeanY)), aes(x=DataTBMatmeanX, y=DataTBMatmeanY) , inherit.aes = F)
    print(P)
  }
  else{
    MatAnnot=as.data.frame(cbind(cbind(DataTBMatmeanX,DataTBMatmeanY), Annotation))
    MatAnnot$Annotation=as.factor(MatAnnot$Annotation)
    P=ggplot(DataTBAxis, aes(x = X, y=Y, color=Var))+
      geom_point()+
      geom_smooth(method=lm)
    P=P+geom_point(data=MatAnnot, aes(x = DataTBMatmeanX, y=DataTBMatmeanY, color=Annotation), inherit.aes = F)+
      scale_color_brewer(palette="Dark2")
    # P=ggplot(data =MatAnnot, aes(x=DataTBMatmeanX, y=DataTBMatmeanY, color=Annotation))+geom_point()+ scale_color_brewer(palette="Dark2")
    print(P)
  }
 
}

bigcor <- function(
  x, 
  y = NULL,
  fun = c("cor", "cov"), 
  size = 2000, 
  verbose = TRUE, 
  ...)
{
  fun <- match.arg(fun)
  if (fun == "cor") FUN <- cor else FUN <- cov
  if (fun == "cor") STR <- "Correlation" else STR <- "Covariance" 
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
  
  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)
  
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST  
  NBLOCKS <- NCOL %/% size
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  if (is.null(y)) resMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))  
  else resMAT <- ff(vmode = "double", dim = c(NCOL, YCOL))
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)  
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  
  ## initiate time counter
  timeINIT <- proc.time() 
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]    
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]    
    
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) cat(sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ", i, STR,  COMB[1],
                               COMB[2], length(G1),  length(G2)))      
      RES <- FUN(x[, G1], x[, G2], ...)
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES) 
    } else ## if y = smaller matrix or vector  
    {
      if (verbose) cat(sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ", i, STR,  COMB[1],
                               length(G1),  YCOL))    
      RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES             
    }
    
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      cat(timeNOW[3], "s\n")
    }
    
    gc()
  } 
  
  return(resMAT)
}
