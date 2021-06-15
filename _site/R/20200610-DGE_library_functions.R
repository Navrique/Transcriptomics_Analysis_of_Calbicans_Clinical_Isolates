# The aim of this file is to define a library of functions that can be used to perform a differential gene expression analysis
library(limma)
library(edgeR)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(sva)
library(pheatmap)
library(hash)
library(tidyr)
library(grDevices)
library(Glimma)
# library(foreach)
# library(doParallel)
library(dendextend)
library(gridExtra)
library(ggbiplot)
library(Hmisc)
library(enrichplot)
# library(M3C)
# PathPrefix="/home/eduranda/Documents/SwitchDrive/CHUV/Sanglard_1/Notebook"
PathPrefix="/media/edurandau/DATA/Eric_Durandau/Sanglard/Labbook"
# source(file.path(PathPrefix,"Workbench/Bioinformatics/R_Tools/RTools.R"))
# source("/media/edurandau/DATA/Eric_Durandau/Sanglard/Labbook/Workbench/Bioinformatics/R_Tools/RTools.R")
# source("/home/eduranda/Documents/SwitchDrive/CHUV/Sanglard_1/Notebook/Workbench/Bioinformatics/R_Tools/RTools.R")
# function to rename the column of a read count table both "ListOfStringToReplace" and "ListOfNameUsed2Replace" should be of the same size and in the right order
SampleID2RealName=function(Table, ListOfStringToReplace, ListOfNameUsed2Replace){
  ColNames=colnames(Table)
  for (i in 1:length(ListOfNameUsed2Replace)){
    bol=str_detect(ColNames, paste0("_",ListOfStringToReplace[i],".b"))
    # ColNames[bol]=str_replace(ColNames[bol],paste0("_",ListOfStringToReplace[i],".") ,paste0("_",ListOfNameUsed2Replace[i],"."))
    ColNames[bol]=ListOfNameUsed2Replace[i]
    # print(paste0("_",ListOfStringToReplace[i],"."))
    # print(paste0("_",ListOfNameUsed2Replace[i],"."))
  }
  Output=Table
  colnames(Output)=ColNames
  # print(colnames(Table))
  # print(colnames(Output))
  return(Output)
}

LoadRawReadCountsCg=function(Path2ReadCount,
                             MetaDataFromSampleName='(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)',
                             ExcludeSamples='(?<Isolate>NI)_(?<replicate>\\d+)_',
                             skip=1,
                             Path2Ref="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/C_glabrata", GeneID=""){
  # Load the count table and extract the metadata from the sample names using "LoadRawReadCounts()". The output of the function is a "hash()" object that will contains gene information, read counts and Metadata. All of these table are produce from the count table itself and the table containing the EntrezID.
  
  ## The function has three optional arguments used to extract the metadata from the sample name "MetaDataFromSampleName" and the other to remove sample experiment from the data set: "ExcludeSamples". Both need to be RegEx patterns. The third argument is to determine whether a some line have to be skipped while loading the read count table.
  
  ## !!! the function has been created using this sample dataset. !!! To be fully compatible, samples names in the count table has to fit a specific REGEX pattern : "(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)". If this is not the case, the user has the possibility to add a first optionals argument to the function that can be a pattern that fit the Sample name pattern. A second optional argument can be specified to remove some sample from the the final tables. Once again this is a REGEX pattern.
  ## 
  
  
  
  # define the location of the Count Data Table file. 
  # define a regex expression that will enable to extract sample information from the file name.
  RegexPattern=MetaDataFromSampleName
  # !!!!!!The regex pattern enables to reject directly the non infected (NI) sample
  RegexPattern2=ExcludeSamples
  
  # load the Read count table and extract the colnames
  SeqData= read.delim(Path2ReadCount, sep = "\t",skip = skip)
  # Remove specified condition using the "MetaDataFromSampleName"
  SeqData=SeqData[,!str_detect(colnames(SeqData), RegexPattern2 )]
  
  # Set the rowname as the Gene_ID
  rownames(SeqData)=SeqData$Geneid
  
  # HTML decode the Note and the Gene column. the note column has not been properly imported some unwanted character need to be converted
  # define the list of character and corresponding replacement
  NoteCorr=SeqData$Note
  Char=c("%28", "%29", "%20", "%3B", "%3B","%3A", "%2C", "%2D", "%2F", "%27")
  CoChar=c("(", ")", " ", ";", ":",",", ",", "-", "/", "'")
  for(i in 1:length(Char)) {
    NoteCorr=str_replace_all(NoteCorr, Char[i],CoChar[i])
  }
  if (length(NoteCorr)!=0){
    SeqData$Note=NoteCorr  
  }
  
  # remove genes that are not expressed. This will include genes with the chromosome B
  
  
  
  ## Using NCBI Table
  ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  AKFeatureTable=read.delim(file.path(Path2Ref, "NCBI_GeneList.txt"))
  # remove some unused character from the "alias" column
  AKFeatureTable$Aliases=str_remove_all(AKFeatureTable$Aliases,"CAALFM_")
  # recover the ORF19 assembly 22 mapping
  ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "Assembly_mapping.tab"))
  ## get Assembly ID and remove the "_" character
  Assembly22=str_remove_all(ORF19_Asssembly22ID$locus_tag.CBS138,"_")
  ## get the ORF19IDs
  orf19IDS=ORF19_Asssembly22ID$locus_tag.DSY562
  Coord22= match(Assembly22, AKFeatureTable$Aliases)
  Coord19=match(orf19IDS, AKFeatureTable$Aliases)
  ORF19_Asssembly22ID$EntrezID=NA
  ORF19_Asssembly22ID$EntrezID[!is.na(Coord22)]=AKFeatureTable$GeneID[Coord22[!is.na(Coord22)]]
  ORF19_Asssembly22ID$EntrezID[!is.na(Coord19)]=AKFeatureTable$GeneID[Coord19[!is.na(Coord19)]]
  rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$locus_tag.DSY562
  
  ## ADDING ENTREZ_ID TO SEQ DATA TABLE
  ## Using AK Table
  ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  # AKFeatureTable=read.delim(file.path(Path2Ref,"setA", "C_albicans_genes_allids_desc.tsv"))
  # rownames(AKFeatureTable)=AKFeatureTable$ASSEMBLY22_ID
  # recover the ORF19 assembly 22 mapping
  # ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "ORF19_Assembly22_mapping.tab"))
  # rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$ASSEMBLY22_ID
  # # Add the EntrezID to the ORF19 assembly22 map table
  # ORF19_Asssembly22ID=merge(ORF19_Asssembly22ID, AKFeatureTable["EntrezID"], by="row.names", all=T)
  # row.names(ORF19_Asssembly22ID)=ORF19_Asssembly22ID$Row.names
  # ORF19_Asssembly22ID=ORF19_Asssembly22ID[,colnames(ORF19_Asssembly22ID)!= "Row.names"]
  
  
  # merge the readcount table with the ORF19 Assembly table
  if (GeneID=="CAG"){
    SeqData=merge(SeqData,ORF19_Asssembly22ID,by.x="row.names",by.y="locus_tag.CBS138",  all = TRUE)
  }else{
    SeqData=merge(SeqData,ORF19_Asssembly22ID,by="row.names", all = TRUE)  
  }
  
  
  
  # use the "stringr" library to find and extract metadata from sample names
  # extract only the content of the colname containing the sample name and clean it up.
  bol=str_detect(colnames(SeqData), RegexPattern)
  # remove undesired column based on the second regex pattern
  bol2=str_detect(colnames(SeqData), RegexPattern2 )
  bol=bol&!bol2
  ExtractSampleName=str_extract(colnames(SeqData), RegexPattern)
  ExtractSampleName[bol2]=NA
  colnames(SeqData)[bol]=ExtractSampleName[bol]
  
  # colnames(SeqData)=str_remove_all(colnames(SeqData), "(Alignement_No_rRNA_)|(.bam)")
  # get colnames
  Samples=colnames(SeqData)[bol]
  Metadata=str_match(Samples, RegexPattern )
  # define the corresponding variable name for the metadata
  VarMetadata=c("Sample", str_remove_all(str_extract_all(RegexPattern, "<(\\w+)>")[[1]], "<|>"))
  # remove empty lines
  Metadata=Metadata[rowSums(is.na(Metadata)) !=ncol(Metadata),]
  # rename column and use the Sample column as index, reorder the "SeqData" Sample column such like match the row name of the Metadata tbl.
  colnames(Metadata)=VarMetadata
  rownames(Metadata)=Metadata[,"Sample"]
  #convert to dataframe
  Metadata=data.frame(Metadata, stringsAsFactors = T)
  #build the condition column
  if (length(colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")])>1){
    Metadata=unite(Metadata, "Conditions", colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")], remove = F)
  }
  else{
    Metadata$Conditions=Metadata[,colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")]]
  }
  Metadata$Conditions=factor(Metadata$Conditions)
  # Add a color colum. A unique colors for each condition
  Metadata$Color=Metadata$Conditions
  levels(Metadata$Color)=rainbow(length(levels(Metadata$Conditions)))
  
  # REMOVE NOT EXPRESSED GENES
  SeqData= SeqData[!is.na(SeqData$Geneid),]
  SeqData=SeqData[rowSums(SeqData[,rownames(Metadata)])!=0,]
  
  #  Set the ORF19 column as rownames
  rownames(SeqData)=SeqData$Geneid
  # Fill Empty Gene values with orf19 ID
  levels(SeqData$Gene)=c(levels(SeqData$Gene), levels(SeqData$Geneid[is.na(SeqData$Gene)]))
  SeqData$Gene[is.na(SeqData$Gene)]=SeqData$Geneid[is.na(SeqData$Gene)]
  
  # Correct the Gene columns. remove Html links characteres
  Genes=SeqData$Gene
  for(i in 1:length(Char)) {
    Genes=str_replace_all(Genes, Char[i],CoChar[i])
  }
  SeqData$Gene=Genes
  # use the "SeqData$Gene" column as rownames
  rownames(SeqData)=Genes
  
  #  Split the SEQ table in two: one with infos one with countTable
  RawCountTabl=SeqData[rownames(Metadata)]
  SeqData=SeqData[is.na(match(colnames(SeqData), rownames(Metadata)))]
  # convert MetaData tbl to a dataframe
  Metadata=data.frame(Metadata)
  # Confirm RawCountTabl Colnames are in the same order as the rownames of the Metadata tbl
  RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]=RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]
  
  # Setup Output Variable in a hash() collection
  Output=hash()
  Output["MetaData"]=data.frame(Metadata, stringsAsFactors = T)
  Output["ReadCounts"]=RawCountTabl
  Output["Genes"]=SeqData
  return(Output)
}

LoadRawReadCountsHs=function(Path2ReadCount,
                             MetaDataFromSampleName='(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)',
                             ExcludeSamples='(?<Isolate>NI)_(?<replicate>\\d+)_',
                             skip=1,Path2Ref="" ){
  # Load the count table and extract the metadata from the sample names using "LoadRawReadCounts()". The output of the function is a "hash()" object that will contains gene information, read counts and Metadata. All of these table are produce from the count table itself and the table containing the EntrezID.
  
  ## The function has three optional arguments used to extract the metadata from the sample name "MetaDataFromSampleName" and the other to remove sample experiment from the data set: "ExcludeSamples". Both need to be RegEx patterns. The third argument is to determine whether a some line have to be skipped while loading the read count table.
  
  ## !!! the function has been created using this sample dataset. !!! To be fully compatible, samples names in the count table has to fit a specific REGEX pattern : "(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)". If this is not the case, the user has the possibility to add a first optionals argument to the function that can be a pattern that fit the Sample name pattern. A second optional argument can be specified to remove some sample from the the final tables. Once again this is a REGEX pattern.
  ## 
  
  
  
  # define the location of the Count Data Table file. 
  # define a regex expression that will enable to extract sample information from the file name.
  RegexPattern=MetaDataFromSampleName
  # !!!!!!The regex pattern enables to reject directly the non infected (NI) sample
  RegexPattern2=ExcludeSamples
  
  # load the Read count table and extract the colnames
  SeqData= read.delim(Path2ReadCount, sep = "\t",skip = skip)
  # Remove specified condition using the "MetaDataFromSampleName"
  SeqData=SeqData[,!str_detect(colnames(SeqData), RegexPattern2 )]
  
  # Set the rowname as the Gene_ID
  rownames(SeqData)=SeqData$Geneid
  
  # HTML decode the Note and the Gene column. the note column has not been properly imported some unwanted character need to be converted
  # define the list of character and corresponding replacement
  NoteCorr=SeqData$Note
  Char=c("%28", "%29", "%20", "%3B", "%3B","%3A", "%2C", "%2D", "%2F", "%27")
  CoChar=c("(", ")", " ", ";", ":",",", ",", "-", "/", "'")
  for(i in 1:length(Char)) {
    NoteCorr=str_replace_all(NoteCorr, Char[i],CoChar[i])
  }
  SeqData$Note=NoteCorr
  # remove genes that are not expressed. This will include genes with the chromosome B
  
  
  
  # ## Using NCBI Table
  # ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  # AKFeatureTable=read.delim(file.path(Path2Ref, "NCBI_GeneList.txt"))
  # # remove some unused character from the "alias" column
  # AKFeatureTable$Aliases=str_remove_all(AKFeatureTable$Aliases,"CAALFM_")
  # # recover the ORF19 assembly 22 mapping
  # ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "Assembly_mapping.tab"))
  # ## get Assembly ID and remove the "_" character
  # Assembly22=str_remove_all(ORF19_Asssembly22ID$locus_tag.CBS138,"_")
  # ## get the ORF19IDs
  # orf19IDS=ORF19_Asssembly22ID$locus_tag.DSY562
  # Coord22= match(Assembly22, AKFeatureTable$Aliases)
  # Coord19=match(orf19IDS, AKFeatureTable$Aliases)
  # ORF19_Asssembly22ID$EntrezID=NA
  # ORF19_Asssembly22ID$EntrezID[!is.na(Coord22)]=AKFeatureTable$GeneID[Coord22[!is.na(Coord22)]]
  # ORF19_Asssembly22ID$EntrezID[!is.na(Coord19)]=AKFeatureTable$GeneID[Coord19[!is.na(Coord19)]]
  # rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$locus_tag.DSY562
  # 
  ## Extracting the EntrezID from the dbxref
  Entrez=str_extract(SeqData$Dbxref, "GeneID:\\d+")
  Entrez=str_remove(Entrez, "GeneID:")
  Entrez=str_remove(Entrez, ",")
  SeqData$EntrezID=as.numeric(Entrez)
  
  ## ADDING ENTREZ_ID TO SEQ DATA TABLE
  ## Using AK Table
  ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  # AKFeatureTable=read.delim(file.path(Path2Ref,"setA", "C_albicans_genes_allids_desc.tsv"))
  # rownames(AKFeatureTable)=AKFeatureTable$ASSEMBLY22_ID
  # recover the ORF19 assembly 22 mapping
  # ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "ORF19_Assembly22_mapping.tab"))
  # rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$ASSEMBLY22_ID
  # # Add the EntrezID to the ORF19 assembly22 map table
  # ORF19_Asssembly22ID=merge(ORF19_Asssembly22ID, AKFeatureTable["EntrezID"], by="row.names", all=T)
  # row.names(ORF19_Asssembly22ID)=ORF19_Asssembly22ID$Row.names
  # ORF19_Asssembly22ID=ORF19_Asssembly22ID[,colnames(ORF19_Asssembly22ID)!= "Row.names"]
  
  
  # merge the readcount table with the ORF19 Assembly table
  SeqData$ORF19_Asssembly22ID=SeqData$EntrezID
  
  
  # use the "stringr" library to find and extract metadata from sample names
  # extract only the content of the colname containing the sample name and clean it up.
  bol=str_detect(colnames(SeqData), RegexPattern)
  # remove undesired column based on the second regex pattern
  bol2=str_detect(colnames(SeqData), RegexPattern2 )
  bol=bol&!bol2
  ExtractSampleName=str_extract(colnames(SeqData), RegexPattern)
  ExtractSampleName[bol2]=NA
  colnames(SeqData)[bol]=ExtractSampleName[bol]
  
  # colnames(SeqData)=str_remove_all(colnames(SeqData), "(Alignement_No_rRNA_)|(.bam)")
  # get colnames
  Samples=colnames(SeqData)[bol]
  Metadata=str_match(Samples, RegexPattern )
  # define the corresponding variable name for the metadata
  VarMetadata=c("Sample", str_remove_all(str_extract_all(RegexPattern, "<(\\w+)>")[[1]], "<|>"))
  # remove empty lines
  Metadata=Metadata[rowSums(is.na(Metadata)) !=ncol(Metadata),]
  # rename column and use the Sample column as index, reorder the "SeqData" Sample column such like match the row name of the Metadata tbl.
  colnames(Metadata)=VarMetadata
  rownames(Metadata)=Metadata[,"Sample"]
  #convert to dataframe
  Metadata=data.frame(Metadata, stringsAsFactors = T)
  #build the condition column
  if (length(colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")])>1){
    Metadata=unite(Metadata, "Conditions", colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")], remove = F)
  }
  else{
    Metadata$Conditions=Metadata[,colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")]]
  }
  Metadata$Conditions=factor(Metadata$Conditions)
  # Add a color colum. A unique colors for each condition
  Metadata$Color=Metadata$Conditions
  levels(Metadata$Color)=rainbow(length(levels(Metadata$Conditions)))
  
  # REMOVE NOT EXPRESSED GENES
  SeqData= SeqData[!is.na(SeqData$Geneid),]
  SeqData=SeqData[rowSums(SeqData[,rownames(Metadata)])!=0,]
  
  #  Set the ORF19 column as rownames
  rownames(SeqData)=SeqData$Geneid
  # Fill Empty Gene values with orf19 ID
  levels(SeqData$Gene)=c(levels(SeqData$Gene), levels(SeqData$Geneid[is.na(SeqData$Gene)]))
  SeqData$Gene[is.na(SeqData$Gene)]=SeqData$Geneid[is.na(SeqData$Gene)]
  
  # Correct the Gene columns. remove Html links characteres
  Genes=SeqData$Gene
  for(i in 1:length(Char)) {
    Genes=str_replace_all(Genes, Char[i],CoChar[i])
  }
  SeqData$Gene=Genes
  # use the "SeqData$Gene" column as rownames
  rownames(SeqData)=Genes
  
  #  Split the SEQ table in two: one with infos one with countTable
  RawCountTabl=SeqData[rownames(Metadata)]
  SeqData=SeqData[is.na(match(colnames(SeqData), rownames(Metadata)))]
  # convert MetaData tbl to a dataframe
  Metadata=data.frame(Metadata)
  # Confirm RawCountTabl Colnames are in the same order as the rownames of the Metadata tbl
  RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]=RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]
  
  # Setup Output Variable in a hash() collection
  Output=hash()
  Output["MetaData"]=data.frame(Metadata, stringsAsFactors = T)
  Output["ReadCounts"]=RawCountTabl
  Output["Genes"]=SeqData
  return(Output)
}

LoadRawReadCounts=function(Path2ReadCount,
                           MetaDataFromSampleName='(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)',
                           ExcludeSamples='(?<Isolate>NI)_(?<replicate>\\d+)_',
                           skip=1,
                           Path2Ref=file.path("..","Data")){
  # Load the count table and extract the metadata from the sample names using "LoadRawReadCounts()". The output of the function is a "hash()" object that will contains gene information, read counts and Metadata. All of these table are produce from the count table itself and the table containing the EntrezID.
  
  ## The function has three optional arguments used to extract the metadata from the sample name "MetaDataFromSampleName" and the other to remove sample experiment from the data set: "ExcludeSamples". Both need to be RegEx patterns. The third argument is to determine whether a some line have to be skipped while loading the read count table.
  
  ## !!! the function has been created using this sample dataset. !!! To be fully compatible, samples names in the count table has to fit a specific REGEX pattern : "(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)". If this is not the case, the user has the possibility to add a first optionals argument to the function that can be a pattern that fit the Sample name pattern. A second optional argument can be specified to remove some sample from the the final tables. Once again this is a REGEX pattern.
  ## 
  
    
  
  # define the location of the Count Data Table file. 
  # define a regex expression that will enable to extract sample information from the file name.
  RegexPattern=MetaDataFromSampleName
  # !!!!!!The regex pattern enables to reject directly the non infected (NI) sample
  RegexPattern2=ExcludeSamples
  
  # load the Read count table and extract the colnames
  SeqData= read.delim(Path2ReadCount, sep = "\t",skip = skip)
  # Remove specified condition using the "MetaDataFromSampleName"
  SeqData=SeqData[,!str_detect(colnames(SeqData), RegexPattern2 )]
  
  # Set the rowname as the Gene_ID
  rownames(SeqData)=SeqData$Geneid
  
  # HTML decode the Note and the Gene column. the note column has not been properly imported some unwanted character need to be converted
  # define the list of character and corresponding replacement
  NoteCorr=SeqData$Note
  Char=c("%28", "%29", "%20", "%3B", "%3B","%3A", "%2C", "%2D", "%2F", "%27")
  CoChar=c("(", ")", " ", ";", ":",",", ",", "-", "/", "'")
  for(i in 1:length(Char)) {
    NoteCorr=str_replace_all(NoteCorr, Char[i],CoChar[i])
  }
  SeqData$Note=NoteCorr
  # remove genes that are not expressed. This will include genes with the chromosome B
  
 
  
    ## Using NCBI Table
    ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  AKFeatureTable=read.delim(file.path(Path2Ref, "20201107-NCBI_GeneList.tab"))
  # remove some unused character from the "alias" column
  AKFeatureTable$Aliases=str_remove_all(AKFeatureTable$Aliases,"CAALFM_")
  # recover the ORF19 assembly 22 mapping
  ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "20200519-Assembly_mapping.tab"))
  ## get Assembly22 ID and remove the "_" character
  Assembly22=str_remove_all(ORF19_Asssembly22ID$ASSEMBLY22_ID,"_")
  ## get the ORF19IDs
  orf19IDS=ORF19_Asssembly22ID$ORF19_ID
  Coord22= match(Assembly22, AKFeatureTable$Aliases)
  Coord19=match(orf19IDS, AKFeatureTable$Aliases)
  ORF19_Asssembly22ID$EntrezID=NA
  ORF19_Asssembly22ID$EntrezID[!is.na(Coord22)]=AKFeatureTable$GeneID[Coord22[!is.na(Coord22)]]
  ORF19_Asssembly22ID$EntrezID[!is.na(Coord19)]=AKFeatureTable$GeneID[Coord19[!is.na(Coord19)]]
  rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$ASSEMBLY22_ID
  
  ## ADDING ENTREZ_ID TO SEQ DATA TABLE
  ## Using AK Table
  ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  # AKFeatureTable=read.delim(file.path(Path2Ref,"setA", "C_albicans_genes_allids_desc.tsv"))
  # rownames(AKFeatureTable)=AKFeatureTable$ASSEMBLY22_ID
  # recover the ORF19 assembly 22 mapping
  # ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "ORF19_Assembly22_mapping.tab"))
  # rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$ASSEMBLY22_ID
  # # Add the EntrezID to the ORF19 assembly22 map table
  # ORF19_Asssembly22ID=merge(ORF19_Asssembly22ID, AKFeatureTable["EntrezID"], by="row.names", all=T)
  # row.names(ORF19_Asssembly22ID)=ORF19_Asssembly22ID$Row.names
  # ORF19_Asssembly22ID=ORF19_Asssembly22ID[,colnames(ORF19_Asssembly22ID)!= "Row.names"]
  
  
  # merge the readcount table with the ORF19 Assembly table
  SeqData=merge(SeqData,ORF19_Asssembly22ID,by="row.names", all = TRUE)
  
  # use the "stringr" library to find and extract metadata from sample names
  # extract only the content of the colname containing the sample name and clean it up.
  bol=str_detect(colnames(SeqData), RegexPattern)
  ExtractSampleName=str_extract(colnames(SeqData), RegexPattern)
  colnames(SeqData)[bol]=ExtractSampleName[bol]
  # colnames(SeqData)=str_remove_all(colnames(SeqData), "(Alignement_No_rRNA_)|(.bam)")
  # get colnames
  Samples=colnames(SeqData)
  Metadata=str_match(Samples, RegexPattern )
  # define the corresponding variable name for the metadata
  VarMetadata=c("Sample", str_remove_all(str_extract_all(RegexPattern, "<(\\w+)>")[[1]], "<|>"))
  # remove empty lines
  Metadata=Metadata[rowSums(is.na(Metadata)) !=ncol(Metadata),]
  # rename column and use the Sample column as index, reorder the "SeqData" Sample column such like match the row name of the Metadata tbl.
  colnames(Metadata)=VarMetadata
  rownames(Metadata)=Metadata[,"Sample"]
  #convert to dataframe
  Metadata=data.frame(Metadata, stringsAsFactors = T)
  #build the condition colum
  
  Metadata=unite(Metadata, "Conditions", colnames(Metadata)[!str_detect(colnames(Metadata), "ample|eplica|batch")], remove = F)
  Metadata$Conditions=factor(Metadata$Conditions)
  # Add a color colum. A unique colors for each condition
  Metadata$Color=Metadata$Conditions
  levels(Metadata$Color)=rainbow(length(levels(Metadata$Conditions)))
  
  # REMOVE NOT EXPRESSED GENES
  SeqData= SeqData[!is.na(SeqData$Geneid),]
  SeqData=SeqData[rowSums(SeqData[,rownames(Metadata)])!=0,]
  
  #  Set the ORF19 column as rownames
  rownames(SeqData)=SeqData$ORF19_ID
  # Fill Empty Gene values with orf19 ID
  levels(SeqData$Gene)=c(levels(SeqData$Gene), levels(SeqData$ORF19_ID[is.na(SeqData$Gene)]))
  SeqData$Gene[is.na(SeqData$Gene)]=SeqData$ORF19_ID[is.na(SeqData$Gene)]

  # Correct the Gene columns. remove Html links characteres
  Genes=SeqData$Gene
  for(i in 1:length(Char)) {
    Genes=str_replace_all(Genes, Char[i],CoChar[i])
  }
  SeqData$Gene=Genes
  # use the "SeqData$Gene" column as rownames
  rownames(SeqData)=Genes
  
  #  Split the SEQ table in two: one with infos one with countTable
  RawCountTabl=SeqData[rownames(Metadata)]
  SeqData=SeqData[is.na(match(colnames(SeqData), rownames(Metadata)))]
  # convert MetaData tbl to a dataframe
  Metadata=data.frame(Metadata)
  # Confirm RawCountTabl Colnames are in the same order as the rownames of the Metadata tbl
  RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]=RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]

  # Setup Output Variable in a hash() collection
  Output=hash()
  Output["MetaData"]=data.frame(Metadata, stringsAsFactors = T)
  Output["ReadCounts"]=RawCountTabl
  Output["Genes"]=SeqData
  return(Output)
}

LoadRawReadCountsSimp=function(Path2ReadCount,
                             MetaDataFromSampleName='(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)',
                             ExcludeSamples='(?<Isolate>NI)_(?<replicate>\\d+)_',
                             skip=1){
  # Load the count table and extract the metadata from the sample names using "LoadRawReadCounts()". The output of the function is a "hash()" object that will contains gene information, read counts and Metadata. All of these table are produce from the count table itself and the table containing the EntrezID.
  
  ## The function has three optional arguments used to extract the metadata from the sample name "MetaDataFromSampleName" and the other to remove sample experiment from the data set: "ExcludeSamples". Both need to be RegEx patterns. The third argument is to determine whether a some line have to be skipped while loading the read count table.
  
  ## !!! the function has been created using this sample dataset. !!! To be fully compatible, samples names in the count table has to fit a specific REGEX pattern : "(?<Isolate>\\w+|\\d+)_(?<replicate>\\d+)_(?<time>\\d+)". If this is not the case, the user has the possibility to add a first optionals argument to the function that can be a pattern that fit the Sample name pattern. A second optional argument can be specified to remove some sample from the the final tables. Once again this is a REGEX pattern.
  ## 
  
  
  
  # define the location of the Count Data Table file. 
  # define a regex expression that will enable to extract sample information from the file name.
  RegexPattern=MetaDataFromSampleName
  # !!!!!!The regex pattern enables to reject directly the non infected (NI) sample
  RegexPattern2=ExcludeSamples
  
  # load the Read count table and extract the colnames
  SeqData= read.delim(Path2ReadCount, sep = "\t",skip = skip)
  # Remove specified condition using the "MetaDataFromSampleName"
  SeqData=SeqData[,!str_detect(colnames(SeqData), RegexPattern2 )]
  
  # Set the rowname as the Gene_ID
  rownames(SeqData)=SeqData$Geneid
  
  # HTML decode the Note and the Gene column. the note column has not been properly imported some unwanted character need to be converted
  # define the list of character and corresponding replacement
  NoteCorr=SeqData$Note
  Char=c("%28", "%29", "%20", "%3B", "%3B","%3A", "%2C", "%2D", "%2F", "%27")
  CoChar=c("(", ")", " ", ";", ":",",", ",", "-", "/", "'")
  for(i in 1:length(Char)) {
    NoteCorr=str_replace_all(NoteCorr, Char[i],CoChar[i])
  }
  if (length(NoteCorr)!=0){
    SeqData$Note=NoteCorr  
  }
  
  # remove genes that are not expressed. This will include genes with the chromosome B
  
  # ## Using NCBI Table
  # ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  # AKFeatureTable=read.delim(file.path(Path2Ref, "NCBI_GeneList.txt"))
  # # remove some unused character from the "alias" column
  # AKFeatureTable$Aliases=str_remove_all(AKFeatureTable$Aliases,"CAALFM_")
  # # recover the ORF19 assembly 22 mapping
  # ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "Assembly_mapping.tab"))
  # ## get Assembly ID and remove the "_" character
  # Assembly22=str_remove_all(ORF19_Asssembly22ID$locus_tag.CBS138,"_")
  # ## get the ORF19IDs
  # orf19IDS=ORF19_Asssembly22ID$locus_tag.DSY562
  # Coord22= match(Assembly22, AKFeatureTable$Aliases)
  # Coord19=match(orf19IDS, AKFeatureTable$Aliases)
  # ORF19_Asssembly22ID$EntrezID=NA
  # ORF19_Asssembly22ID$EntrezID[!is.na(Coord22)]=AKFeatureTable$GeneID[Coord22[!is.na(Coord22)]]
  # ORF19_Asssembly22ID$EntrezID[!is.na(Coord19)]=AKFeatureTable$GeneID[Coord19[!is.na(Coord19)]]
  # rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$locus_tag.DSY562
  
  ## ADDING ENTREZ_ID TO SEQ DATA TABLE
  ## Using AK Table
  ## load the gene allele ID table build by AK. Contains all IDs for each genes.
  # AKFeatureTable=read.delim(file.path(Path2Ref,"setA", "C_albicans_genes_allids_desc.tsv"))
  # rownames(AKFeatureTable)=AKFeatureTable$ASSEMBLY22_ID
  # recover the ORF19 assembly 22 mapping
  # ORF19_Asssembly22ID=read.delim(file.path(Path2Ref, "ORF19_Assembly22_mapping.tab"))
  # rownames(ORF19_Asssembly22ID) =ORF19_Asssembly22ID$ASSEMBLY22_ID
  # # Add the EntrezID to the ORF19 assembly22 map table
  # ORF19_Asssembly22ID=merge(ORF19_Asssembly22ID, AKFeatureTable["EntrezID"], by="row.names", all=T)
  # row.names(ORF19_Asssembly22ID)=ORF19_Asssembly22ID$Row.names
  # ORF19_Asssembly22ID=ORF19_Asssembly22ID[,colnames(ORF19_Asssembly22ID)!= "Row.names"]
  
  
  # # merge the readcount table with the ORF19 Assembly table
  # if (GeneID=="CAG"){
  #   SeqData=merge(SeqData,ORF19_Asssembly22ID,by.x="row.names",by.y="locus_tag.CBS138",  all = TRUE)
  # }else{
  #   SeqData=merge(SeqData,ORF19_Asssembly22ID,by="row.names", all = TRUE)  
  # }
  # 
  
  
  # use the "stringr" library to find and extract metadata from sample names
  # extract only the content of the colname containing the sample name and clean it up.
  bol=str_detect(colnames(SeqData), RegexPattern)
  # remove undesired column based on the second regex pattern
  bol2=str_detect(colnames(SeqData), RegexPattern2 )
  bol=bol&!bol2
  ExtractSampleName=str_extract(colnames(SeqData), RegexPattern)
  ExtractSampleName[bol2]=NA
  colnames(SeqData)[bol]=ExtractSampleName[bol]
  
  # colnames(SeqData)=str_remove_all(colnames(SeqData), "(Alignement_No_rRNA_)|(.bam)")
  # get colnames
  Samples=colnames(SeqData)[bol]
  Metadata=str_match(Samples, RegexPattern )
  # define the corresponding variable name for the metadata
  VarMetadata=c("Sample", str_remove_all(str_extract_all(RegexPattern, "<(\\w+)>")[[1]], "<|>"))
  # remove empty lines
  Metadata=Metadata[rowSums(is.na(Metadata)) !=ncol(Metadata),]
  # rename column and use the Sample column as index, reorder the "SeqData" Sample column such like match the row name of the Metadata tbl.
  colnames(Metadata)=VarMetadata
  rownames(Metadata)=Metadata[,"Sample"]
  #convert to dataframe
  Metadata=data.frame(Metadata, stringsAsFactors = T)
  #build the condition column
  if (length(colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")])>1){
    Metadata=unite(Metadata, "Conditions", colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")], remove = F)
  }
  else{
    Metadata$Conditions=Metadata[,colnames(Metadata)[ !str_detect(colnames(Metadata), "ample|eplica|batch")]]
  }
  Metadata$Conditions=factor(Metadata$Conditions)
  # Add a color colum. A unique colors for each condition
  Metadata$Color=Metadata$Conditions
  levels(Metadata$Color)=rainbow(length(levels(Metadata$Conditions)))
  
  # REMOVE NOT EXPRESSED GENES
  SeqData= SeqData[!is.na(SeqData$Geneid),]
  SeqData=SeqData[rowSums(SeqData[,rownames(Metadata)])!=0,]
  
  #  Set the ORF19 column as rownames
  rownames(SeqData)=SeqData$Geneid
  # Fill Empty Gene values with orf19 ID
  levels(SeqData$Gene)=c(levels(SeqData$Gene), levels(SeqData$Geneid[is.na(SeqData$Gene)]))
  SeqData$Gene[is.na(SeqData$Gene)]=SeqData$Geneid[is.na(SeqData$Gene)]
  
  # Correct the Gene columns. remove Html links characteres
  Genes=SeqData$Gene
  for(i in 1:length(Char)) {
    Genes=str_replace_all(Genes, Char[i],CoChar[i])
  }
  SeqData$Gene=Genes
  # use the "SeqData$Gene" column as rownames
  rownames(SeqData)=Genes
  
  #  Split the SEQ table in two: one with infos one with countTable
  RawCountTabl=SeqData[rownames(Metadata)]
  SeqData=SeqData[is.na(match(colnames(SeqData), rownames(Metadata)))]
  # convert MetaData tbl to a dataframe
  Metadata=data.frame(Metadata)
  # Confirm RawCountTabl Colnames are in the same order as the rownames of the Metadata tbl
  RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]=RawCountTabl[,match(rownames(Metadata), colnames(RawCountTabl))]
  
  # Setup Output Variable in a hash() collection
  Output=hash()
  Output["MetaData"]=data.frame(Metadata, stringsAsFactors = T)
  Output["ReadCounts"]=RawCountTabl
  Output["Genes"]=SeqData
  return(Output)
}


# The output of the previous function is a collection of variable that can be used as input for data Normalisation
ReadCountNormalisation=function(ReadCounts, Metadata, Factor,  Genes,AdditionalFactor=c(),Plots=T){
  # create output varible
  Output=hash()
  # Convert the datatable to a DGElist object
  # DGEList enables to create an object of type S4 it. running the "DGEList("Data")" will create the object from a matrix of count. 
  # Rows need to be genes and Column correspond to the sample. We can view the names of Col and Row using dimname.DGEList(Object).
  # # The object is organised like a structure. I does can contains different information enabling to groups data. The subsection "sample" is generated from the initial table. The table contains a table with the number of total reads, the normalisation factor and groups. note that, again dim() colnames(), rownames() functions can be used to extract the relevant information by giving the Objec$subsection ars argument
  if (sum(colnames(Metadata)=="condition")==0){
    Metadata$condition = "FP"
  }
  DGE_RawCount=DGEList(ReadCounts)
# Add gene annotation
  # DGE_RawCount$genes = merge(DGE_RawCount$genes, Genes, by="row.names", all=T)
  DGE_RawCount$genes=Genes
  # add our own groups to the object structure into the "sample" subsection.
  DGE_RawCount$samples$group_Cond=Factor
  DGE_RawCount$samples$Color=Factor
  for (i in AdditionalFactor){
    DGE_RawCount$samples[i]=Metadata[i]
  }
  # levels(DGE_RawCount$samples$Color)=brewer.pal(length(levels(Factor)), "Paired")
  levels(DGE_RawCount$samples$Color)=rainbow(length(levels(Factor)))
  DGE_RawCount$samples$Color=as.character(DGE_RawCount$samples$Color)
  # Add Annotation (Gene information) to the DGElist object
    # add the gene info to the DGE list. Add the content of the "SeqData" variable to the "gene" subvariable of the DGE list.
  # DGE_RawCount$gene=SeqData
  
                          # Normalisation of Data : remove lowly expressed genes
    # transform the raw reads count in CPM. Note that de function add an extra offset of 2/meanLibSize.
  CPM=cpm(DGE_RawCount)
    # same operation but returning log2 value. 2 reads are added to each observation to avoid "log0" if the some feature do note have reads.
  Log2CPM= cpm(DGE_RawCount, log = T, prior.count = 2)
    # Calculating the average and median sequencing depth per million among all samples. This will be used to calculate a cutoff while comparing the expression pattern of the different sample
  Mean_SeqDepth <- mean(DGE_RawCount$samples$lib.size) * 1e-6 # for threshold of the filter GeneExpress
  Med_SeqDepth <- median(DGE_RawCount$samples$lib.size) * 1e-6 # for offset added while calculating the cpm.
  
  # Remove lowly expressed genes. Using the "filterByExpr" function from the "egdeR" package, one can remove lowly expressed genes automatically while keeping genes with worthwhile counts.
    # First view the number of non expressed genes (no reads at all)
  NbNoneExpGenes=table(sum(rowSums(DGE_RawCount$counts==0)==ncol(DGE_RawCount$counts)))
  # print(str_c(as.character(NbNoneExpGenes), " genes above ", as.character(nrow(DGE_RawCount$counts)), " has been filtered out because they do not have any reads in any conditions"))
    # calculate the cutoff used by the filterByExpr() function()
  Log_FilterCutoff=log2(10/Med_SeqDepth +2/Mean_SeqDepth)
    # finding lowly expressed genes using filterByExpr from the edgeR library. This will check count value among a given groups of sample. Here we choose the "group_condition" as factor. this return a one column binary table with gene name as index. The filter a consider a threshold at 10 reads relative to the median library size  (10 reads/MedianLibSize per million) + the CPM offset (2)
  
  Expressed_Genes <- filterByExpr(DGE_RawCount, group=DGE_RawCount$samples$group_Cond)
    # Update the DGE object by removing lowly expressed genes from the table. Correct also the lib.size variable
  DGE_Filtered_RawCount=DGE_RawCount[Expressed_Genes,,keep.lib.sizes=F]
  # Add to the Output
  Output["Filter_Genes"]=DGE_Filtered_RawCount
  
  
  #                   Data normalisation : TMM
  # Remove the compositions bias using the TMM (Trimed mean  of M-value (Robinson and Oshlack 2010))
  ## Normalizing gene expression distributions, remove Composition Bias using the TMM method.
  ## Artificially add noise to a copy of the Filtered count table for comparison.
  DGE3_cp=DGE_Filtered_RawCount
  DGE3_cp$counts[,1]=ceiling(DGE3_cp$counts[,1]*0.05)
  DGE3_cp$counts[,2]=DGE3_cp$counts[,2]*5
  lcpm_DGE3_cp=cpm(DGE3_cp, log = T, prior.count = 2)
  # normalise the DGE3_cp
  DGE3_cp=calcNormFactors(DGE3_cp,method="TMM")
  # Correct the real Data
  DGE3_TMM_FiltCount=calcNormFactors(DGE_Filtered_RawCount, method = "TMM") 
  # Calculate the logCPM values or the corrected data
  LCP_DGE3_TMM_FilterCount=cpm(DGE3_TMM_FiltCount, log = T, prior.count = 2)
  Output["TMM_correction"]=DGE3_TMM_FiltCount
  
                          # Voom transform of counts  Removing heteroscedascity from count data
  
  # The basic principle is to remove a bias related to the relationship between gene expression variance and gene expression mean. The more condition/replicate we have the better it is. voom aim at prepare the data for a further linear fit see chapter 6.2 at https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
  
  #Setting up design matrix for the condition. The input of the model matrix must be bolean matrix. row correspond to each samples and col correspond to group, here group_cond. this is a manner to tell voom which samples must be treated together while performing the transformation. 
  ModelCond=model.matrix(~0+Factor)
  # perform the voom transform. We end up with another object
  Voom_DGE3 <- limma::voom(DGE3_TMM_FiltCount, ModelCond, plot=F, normalize="quantile")
  # Voom_DGE3 <- voom(DGE3_TMM_FiltCount, ModelCond, plot=F)
  # perform the linear fit and plot the resulting correction
  vfit <- limma::lmFit(Voom_DGE3, ModelCond)
  # vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  Output["Voom_Raw"]=Voom_DGE3
  
                                # Adjusting for the Batch effects by Unsupervised clustering of samples. 
# This step might not be essential since all samples are from the same batch. This can be interesting if we do have sequencing files generated from different runs or lanes. a GOOD example would be 
  # We will use the SVA package to infer the batch effects. Very quickly sva estimate the number of covariate (surrogate variable) and then adjust value to remove the noise associated to such covariate.  
  # Then adjust for the batch effects using Surrogate variables ### Estimate the surrogate variable and run the SVA
  # build the design matrix with the condition you care about
  mod_SurrogateVar=model.matrix(~0 +Factor, data=DGE3_TMM_FiltCount$samples)
  mod0=model.matrix(~1, data=DGE3_TMM_FiltCount$samples)
  # #Using SVA package to infer the batch effects
  #The number of factors can also be estimated using the num.sv function
  # calculate the number of confounding factor ()
  # This step estimate the number of confounding factors then used to find the surrogate variables. this is the unsupervised part
  Nb_LatentFactors = num.sv(Voom_DGE3$E,mod_SurrogateVar)# we apply the sva function to estimate the number of latent variables using
  # print(paste0("number of latent variable = ", as.character(Nb_LatentFactors)))
  # Surrogate_Var=sva(Voom_DGE3$E,mod_SurrogateVar,mod0, n.sv = Nb_LatentFactors)
  Surrogate_Var=sva(Voom_DGE3$E,mod_SurrogateVar,mod0)# if not specified n.sv can be added manually 
  # print(c("number of Surrogate variable = ", as.character(Surrogate_Var)))
  # Cleanup the data
  VoomTrans_ExpressVal = cleaningP(Voom_DGE3$E, mod_SurrogateVar, Surrogate_Var)
  #Replacing the counts from original voom object with batch corrected expression values
  Voom_DGE3_cp=Voom_DGE3
  Voom_DGE3_cp$E=VoomTrans_ExpressVal
  Output["Norm_Data"]=Voom_DGE3_cp
  Output["Expressed_Genes"]=DGE_RawCount$counts>Log_FilterCutoff
  
#                                             !!!!!!  PLOT SECTION
  if (is.character(Plots) | isTRUE(Plots)){
    # Get the number of Sample from the DGEobject
    Nbsamples <- ncol(DGE_Filtered_RawCount)
    
   
    # set a new figure
    
                  # !!!!!!! Filtered genes based on a low expression!!!!
    # Plot the Log-CPM value to appreciate the efficiency of the normalization using the logCPM. This enables to see how comparable samples are. The par() function is used to combine multiple plot on the same canvas.
    # PLot the Raw data as a density plot
    # Plot the first column
    DT_LOG2CPM=as.data.table(Log2CPM, keep.rownames = T)
    Melt_DT_LPG2CPM=melt(DT_LOG2CPM)
    Melt_DT_LPG2CPM$Corr=rep("Raw Data", nrow(Melt_DT_LPG2CPM))
    DT_FILT_LOG2CPM=as.data.table(cpm(DGE_Filtered_RawCount, log = T, prior.count = 2), keep.rownames=T)
    Melt_DT_FILT_LPG2CPM=melt(DT_FILT_LOG2CPM)
    Melt_DT_FILT_LPG2CPM$Corr=rep("Filtered Data", nrow(Melt_DT_FILT_LPG2CPM))
    Binded=rbind(Melt_DT_LPG2CPM,Melt_DT_FILT_LPG2CPM)
    Binded$variable_Simp=""
    
    for (i in 1:nrow(Metadata)){
      Binded$variable_Simp[Binded$variable==Metadata$Sample[i]]=as.character(Metadata$Conditions[i])
    }

    
    
    P11=ggplot(Binded, aes(x=value, Fill=variable, color=variable_Simp))+
      geom_density()+
      geom_vline(aes(xintercept=Log_FilterCutoff),
                 color="blue", linetype="dashed", size=1)+
      facet_wrap(~Corr)+
      labs(title = "Low expression filtering", y="Gene density", x="# reads [log2(cpm)]", col="Conditions")+mynamesthemeSmall

    
    if (!is.logical(Plots)){
      PlotExportPath=file.path(Plots,"Normalisation","Plots")
      dir.create(file.path(Plots,"Normalisation"),showWarnings = F)
      dir.create(PlotExportPath,showWarnings = F)
      ggsave(filename = file.path(PlotExportPath,"Density_filtering_Low_Exp_Genes.svg"),width = 10, height = 5)
    }
    
    # Plot the library size
    # Plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.
    # The names argument tells the barplot to use the sample names on the x-axis
    # The las argument rotates the axis names
    # if (!is.logical(Plots)){
    #   svg(filename = file.path(PlotExportPath,"Lirary_size.svg"),width = 10, height = 10, pointsize = 10)
    # }
    # par(mfrow=c(1,1))
    Data=DGE_Filtered_RawCount$samples
    Data$Samples=rownames(Data)
    Data$Replicate=Metadata$replicate
    
    P2=ggplot(Data,aes(x=group_Cond, y=lib.size, fill=Replicate))+
      geom_bar(stat="identity", position=position_dodge())+
      ggtitle("Library size")+ mynamestheme+
      theme(axis.text.x = element_text(angle = 90))+labs(title = "Library size", y="# reads [-]", x="Sample")
    
    if (!is.logical(Plots)){
      ggsave(plot = P2,filename = file.path(PlotExportPath,"Library_size.svg"),width = 8, height = 8)
    }
    # barplot(DGE_Filtered_RawCount$samples$lib.size,names=str_replace_all(DGE_Filtered_RawCount$samples$group_Cond, "_", " "),las=2, cex.axis = 0.8,cex.names = 1)
    # # Add a title to the plot
    # title("Barplot of library sizes")
    # if (!is.logical(Plots)){
    #   dev.off()
    # }
    
    #                           # !!!!!!!!! TMM based data correction
    # if (!is.logical(Plots)){
    #   svg(filename = file.path(PlotExportPath,"TMM_Normalisation.svg"),width = 20, height = 10, pointsize = 10)
    # }
    # 
    #Create a two plots canvas
    # par(mfrow=c(1,2))
    # order the variable based on the condition name and recover the index
    SortedNamesIdx=sort.int(DGE3_cp$samples$group_Cond, index.return = T)
    Raw=as.data.table(lcpm_DGE3_cp[,SortedNamesIdx$ix])
    Raw$Corr="Raw data"
    TMM=as.data.table(cpm(DGE3_cp,log = T, prior.count = 2)[,SortedNamesIdx$ix])
    TMM$Corr="TMM"
    Binded=rbind(Raw,TMM)
    Melt_Binded=melt(as.data.table(Binded))
    Melt_Binded$variable_Simp=""
    for (i in 1:nrow(Metadata)){
      Melt_Binded$variable_Simp[Melt_Binded$variable==Metadata$Sample[i]]=as.character(Metadata$Conditions[i])
    }
    P3=ggplot(Melt_Binded, aes(x=variable,y =value, fill=variable_Simp))+
      geom_boxplot()+ facet_wrap(~Corr)+
      theme(axis.text.x = element_text(angle = 90))+mynamesthemeSmall+
      labs(title = "TMM correction", y="Gene reads count [log(cpm)]", x="Sample")
      
    if (!is.logical(Plots)){
      ggsave(plot = P3,filename = file.path(PlotExportPath,"TMM_Normalisation.svg"),width = 20, height = 10)
    }
    # Plot the Boxplot for the lcpm_DGE3_cp unormalized
    # boxplot(lcpm_DGE3_cp[,SortedNamesIdx$ix], las = 2, col= DGE3_TMM_FiltCount$samples$Color[SortedNamesIdx$ix], main="",names=str_replace_all(DGE_Filtered_RawCount$samples$group_Cond[SortedNamesIdx$ix], "_", " "))
    # title(main ="Un-normalised data", ylab = "Log-cpm")
    # 
    # # Plot the normalized data, see whether the artificial noise has been corrected
    # boxplot(cpm(DGE3_cp,log = T, prior.count = 2)[,SortedNamesIdx$ix], las=2,col=DGE3_TMM_FiltCount$samples$Color[SortedNamesIdx$ix], main="", names=str_replace_all(DGE_Filtered_RawCount$samples$group_Cond[SortedNamesIdx$ix], "_", " "))
    # title(main="Corrected TMM normalisation")
    # if (!is.logical(Plots)){
    #   dev.off()
    # }
   
    #                       !!!!!!!!!!!! Voom correction
    if (!is.logical(Plots)){
      svg(filename = file.path(PlotExportPath,"Voom_correction.svg"),width = 20, height = 10, pointsize = 10)
    }
    # par(mfrow=c(1,2))
    P41=voom(DGE3_TMM_FiltCount, ModelCond, save.plot = T)
    
    # P42=recordPlot(plotSA(efit, main="Final model: Mean-variance trend",))
    VoomModelPlot=as.data.table(efit$genes)
    VoomModelPlot[, `:=`(x=P41$voom.xy$x, y=P41$voom.xy$y, ablineX=P41$voom.line$x,ablineY=P41$voom.line$y, Correction="Not Corrected")]
    VoomCorrplot=data.table(efit$genes)
    VoomCorrplot[,`:=`(x=efit$Amean , y = sqrt(efit$sigma), ablineX=efit$Amean, ablineY=sqrt(sqrt(efit$s2.prior)), Correction="Corrected")]
    
    P4=ggplot(rbind(VoomModelPlot,VoomCorrplot), aes(x=x, y=y), size=0.7)+
      facet_wrap(~Correction,scales = "free")+
      geom_point()+mynamestheme+
      geom_line(inherit.aes = F,aes(x=ablineX, y=ablineY),color="red", size=1.5)+labs(x="Average expression [log2(cpm)]", y="Square root variance/sigma [-]")
    
    if (!is.logical(Plots)){
      dev.off()
    }
    
    #                       !!!!!!!!!!!!Batch effect correction
    
    
    Home=plotMDS(Voom_DGE3$E,plot = F)
    Home2=plotMDS(VoomTrans_ExpressVal,plot = F)
    
    Val=lapply(list(Home, Home2),function(x){
      Tbl=cbind(data.frame(Dim1=x$x, Dim2=x$y),
                Metadata)
    })
    PCA=bind_rows(Val)
    PCA$Corrected="With batch effects"
    PCA$Corrected[duplicated(PCA$Sample)]="Batch effect corrected"
    
    P5=ggplot(PCA,aes(x=Dim1, y=Dim2, color=Conditions))+
      facet_wrap(~Corrected)+
      geom_point(aes(Replicate=replicate, Condition=condition))+mynamestheme+
      labs(x="Leading LogFC Dim 1",y="Leading LogFC Dim 1")
    if (!is.logical(Plots)){
      ggsave(plot = P5,filename = file.path(PlotExportPath,"BatchEffectRemoval_PCA.svg"),width = 20, height = 10)
    }
  
    # observing the density of batch corrected and batch not corrected values
    VoomDGE3_melt_Not_Corr=melt(as.data.table(Voom_DGE3$E, keep.rownames = T), id.vars = "rn",variable.name = "sample")
    VoomDGE3_melt_Not_Corr=VoomDGE3_melt_Not_Corr[Metadata, on=.(sample=Sample)]
    VoomDGE3_melt_Not_Corr[,Batch:="Not_Corrected"]
    VoomDGE3_melt_Corr=melt(as.data.table(VoomTrans_ExpressVal, keep.rownames = T), id.vars = "rn",variable.name = "sample")
    VoomDGE3_melt_Corr=VoomDGE3_melt_Corr[Metadata, on=.(sample=Sample)]
    VoomDGE3_melt_Corr[,Batch:="Corrected"]
    Voom_DGE3_melt_bind=rbind(VoomDGE3_melt_Not_Corr, VoomDGE3_melt_Corr)
    
    P6=ggplot(Voom_DGE3_melt_bind, aes(x=value, color=Conditions))+
      facet_wrap(~Batch)+
      geom_density(aes(Replicate=replicate, Condition=condition))+mynamestheme+
      labs(x="Voom-lcpm [log2(cpm)]",y="Density [-]")
    
    if (!is.logical(Plots)){
      # dev.off()
      ggsave(plot = P6,filename = file.path(PlotExportPath,"BatchEffectRemoval_Density.svg"),width = 20, height = 10)
    }
    #             !!!!! Visualise the efficiency of the batch correction using a heatmap
    # Hierarchical heatmap by condition
    # Transform the normalized counts from the dds_smoc2 object using the vst() function with the blind argument and save to vsd_smoc2
    # Extract the matrix of transformed normalized counts from the vsd_smoc2 object using the assay() function and save as vsd_mat_smoc2.
    # Calculate the correlation values between samples and save to vsd_cor_smoc2.
    # Create a heatmap of the correlation values using pheatmap() with an annotation bar designating condition from the smoc2_metadata data frame.
    #compute the correlation value between samples
    if (!is.logical(Plots)){
      svg(filename = file.path(PlotExportPath,"BatchEffectRmv_NotRm_HeatMap.svg"),width = 20, height = 20, pointsize = 10)
    }
    CorCondWiseRaw= cor(Voom_DGE3$E)
    CorCondWiseCorr= cor(VoomTrans_ExpressVal)
    # Plot the heatmap
    Annotation= data.frame(Voom_DGE3$targets$group_Cond,Voom_DGE3$targets[AdditionalFactor], row.names = rownames(Voom_DGE3$targets))
    rownames(Annotation)=rownames(CorCondWiseRaw)
    colnames(Annotation)=c("Sample_groups", AdditionalFactor)
    # par(mfrow=c(1,2))
    P7=pheatmap(CorCondWiseRaw, annotation =Annotation ,treeheight_row=20, treeheight_col = 20,  main= "Not corrected")[[4]]
    if (!is.logical(Plots)){
      dev.off()
    }
    if (!is.logical(Plots)){
      svg(filename = file.path(PlotExportPath,"BatchEffectRmv_HeatMap.svg"),width = 20, height = 20, pointsize = 10)
    }
    
    P8=pheatmap(CorCondWiseCorr, annotation =Annotation,treeheight_row=20, treeheight_col = 20,main= "Corrected" )[[4]]
    if (!is.logical(Plots)){
      dev.off()
    }
    # P9=recordPlot(grid.arrange(arrangeGrob(grobs= list(P7,P8),ncol=2)))
    Output["GGplot"]=list("Low_Expression_Filtering"=P11,
                          "Library_size"=P2,
                          "TMM_Correction"=P3,
                          "Voom_transform"=P4,
                          "Batch_effect_PCA"=P5,
                          "Batch_effect_density"=P6,
                          "Batch_effect_heatmap"=list(P7,P8))
  }
  
  
  
  return(Output)
}

#Function to get a clean data after removing the batch effects. 
cleaningP = function(x, mod_f1_f2, sva1_f1_f2,  P=ncol(mod_f1_f2)) {
  X=cbind(mod_f1_f2,sva1_f1_f2$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(x))
  cleany=x-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

Genes_Diff_Expressed =function(NormData,Factor4Design,Group4InteractivePlot=Factor4Design, RefCond=levels(Factor4Design) , CompareCond=levels(Factor4Design), AdditionalFactor=c() , Method="", Plots=T, lfc_cuttoff = 0.5, pvalue_cuttoff = 0.05){
  # NormData = product of the normalisation step. This is an object output of the ReadCountNormalisation function
    # Factor4Design = a factor describing how data should be associated. The list should tell to which group a sample is associated. Usually, the group information can be found in the "NormData" object.
  # RefCond= list of condition from the "Factor4Design" that should be used as reference. By default every levels of the "Factor4Design" are used as reference for comparison with all other condition this is basically testing all possible combination
  # CompareCond= Same as "RefCond" but for the list of condition that have to be compared to the reference.
  # Method = Tells how to build deal with the list of reference and compared conditions. 
      # The "" is the default and will lead to compare each condition from the "ComparedCond" to each condition in the "Reference"
      # "Paired" will assemble pairewise conditions from "ComparedCond" and "Reference". Both "ComparedCond" and "Reference" need to have the same size.
      # The Defaut method is ""
  # Plot = whether to plot or not the result and eventually save the figure
  # setup the Output variable
  Output=hash()
  #                       Differential Gene expression
  #Setting up design matrix. This will define which samples are replicates.
  Mat_Model_Cond=model.matrix(~0+Factor4Design)
  #Modify the Column name in the "Mat_Model_Cond" to make it more easy to read.
  colnames(Mat_Model_Cond) <- gsub("Factor4Design", "X", colnames(Mat_Model_Cond))
  # Create a contrast map. This matrix define which samples are compared. This matrix define who is the reference and who is the compared condition.
  # Define the list of condition to compare and the list of references.
  References=RefCond
  ToCompare=CompareCond
  # Build the contrast map
  ContrastMap=BuildContrastMap(References,ToCompare,Method = Method)
  
  # Add a X to the rowname of the contrast map if needed and check that the resulting rowname fit the colnames of the Mat Model
  # get the colnames
  RowNameContrast=rownames(ContrastMap)
  Annotation=data.frame(rownames(ContrastMap))
  extraAnnotDF=data.frame(NormData$targets[!duplicated(NormData$targets$group_Cond),c("group_Cond",AdditionalFactor)])
  colnames(extraAnnotDF)=c("group_Cond",AdditionalFactor)
  rownames(extraAnnotDF)=extraAnnotDF$group_Cond
  Annotation=merge(Annotation, extraAnnotDF, by.x = colnames(Annotation), by.y = "group_Cond")
  colnames(Annotation)=colnames(extraAnnotDF)
  rownames(Annotation)=Annotation$group_Cond
  RowNameAnnot=rownames(Annotation)
  
  RowNameContrast[!str_detect(RowNameContrast, "^X")]=str_replace_all(RowNameContrast[!str_detect(RowNameContrast, "^X")], "^","X")
  RowNameAnnot[!str_detect(RowNameAnnot, "^X")]=str_replace_all(RowNameAnnot[!str_detect(RowNameAnnot, "^X")], "^","X")
  rownames(ContrastMap)=RowNameContrast
  rownames(Annotation)=RowNameAnnot
  # colnames()=str_detect(rownames(ContrastMap))
  
  #Fit the linear model
  LMFit_Data <- lmFit(NormData, Mat_Model_Cond) #produce the linear model
  LMFit_Data=LMFit_Data[,rownames(ContrastMap)]
  ContrastMap=ContrastMap[colnames(LMFit_Data$coefficients),]
  Contrast_LMFit_Data <- contrasts.fit(LMFit_Data, contrasts=ContrastMap )
  
  efit_Data <- eBayes(Contrast_LMFit_Data) #extract the relevant statistics
  # efit_Data$genes=efit_Data$genes[,c("Genename", "GENE_NAME","ORF19_ID","ASSEMBLY22_ID", "EntrezID", "Note")]
  tfit=treat(Contrast_LMFit_Data)
  # get the differentially expressed genes
  Find_DiffReg_genes <- decideTests(efit_Data,lfc = lfc_cuttoff, p.value = pvalue_cuttoff, adjust.method = "BH")
 
  #summary(decideTests(efit_nosv_isolate_f1_f2))
  Summary_DiffReg_genes=summary(Find_DiffReg_genes)
  Find_DiffReg_genes=as.data.frame(Find_DiffReg_genes)
  
  Coefs=colnames(efit_Data$coefficients)
  DEgenes=hash()
  
  Ref_Meta=Annotation[match(RefCond,Annotation$group_Cond),]
  Test_Meta=Annotation[match(CompareCond,Annotation$group_Cond),]
  if (nrow(Ref_Meta)!=nrow(Test_Meta)){
    Annotation2=as.data.frame(Test_Meta$group_Cond)
    colnames(Annotation2)="Tested_Condition"
  }else{
    BolAnnotCol=colSums(Ref_Meta==Test_Meta)==nrow(Ref_Meta)
    Annotation2=Test_Meta[,AdditionalFactor]
  }
  
  if(is.data.frame(Annotation2)){
  }
  else{
    Annotation2=as.data.frame(Ref_Meta[,BolAnnotCol])
    colnames(Annotation2)=names(BolAnnotCol)[BolAnnotCol]
  }
  rownames(Annotation2)=Coefs
  ind=1
  for (i in Coefs){
  # cl=parallel::makeForkCluster(10)
  # doParallel::registerDoParallel(cl)
  # foreach (j =1:length(Coefs), .combine=c) %dopar% {
  # i=Coefs[j]
   # print(i)
    DEgenes[i]=topTable(efit_Data, coef=i, number=nrow(NormData$E), sort.by = "logFC",adjust.method = "BH",p.value = pvalue_cuttoff,lfc = lfc_cuttoff)
    basal.vs.ml <- topTreat(tfit, coef=i, n=Inf,adjust.method = "BH", sort.by = "logFC",p.value = pvalue_cuttoff,lfc = lfc_cuttoff)
    ind=ind+1
    # plot the VolcanoPlot
    
    if (is.character(Plots)){
      if (!dir.exists(Plots)){
        dir.create(Plots,showWarnings = F)
      }
      dir.create(Plots,showWarnings = F)
      PlotExportPath=Plots
        if (!dir.exists(file.path(PlotExportPath,i))){
          dir.create(PlotExportPath,showWarnings = F)
          glXYPlot(x=efit_Data$coefficients[,i], y=efit_Data$lods[,i], xlab="logFC", ylab="-log10(Pval)",status=Find_DiffReg_genes[,i],counts=NormData$E, groups=Group4InteractivePlot[,1], side.main="ID",sample.cols =Group4InteractivePlot[,2],main=i, folder = i, path = PlotExportPath, launch = F)
          glMDPlot(efit_Data, coef=i, status=Find_DiffReg_genes[,i], main=i,
                   side.main="Gene", counts=NormData$E,sample.cols =Group4InteractivePlot[,2] , groups=Group4InteractivePlot[,1],ylim=c(--10,13), launch= F, path = PlotExportPath,folder = i, anno =  efit_Data$genes)
         }
      write.table(DEgenes[[i]],file.path(Plots,i,str_c("TopTable", ".tsv")), sep = "\t", row.names = FALSE)
      write.table(basal.vs.ml,file.path(Plots,i,str_c("Treat_table", "_t.tsv")), sep = "\t", row.names = FALSE)
      }
    
   }
  doParallel::stopImplicitCluster()
  Output=DEgenes
 
                             # PLOT SECTION
  if (is.character(Plots) | is.logical(Plots)){
    if (!isFALSE(Plots)) {
      if (!is.logical(Plots)){
        svg(filename = file.path(PlotExportPath,"Correlation_HeatMap.svg"),width = 20, height = 60, pointsize = 10)
      }
      # Plot the correlation between conditions using the LMFitData. this is to observe how much the global expression pattern differ from one condition to another
      #compute the correlation value between samples
      CorCondWise= cor(LMFit_Data$coefficients)
      # Annotation= data.frame(colnames(LMFit_Data$coefficients), row.names = colnames(LMFit_Data$coefficients))
      # rownames(Annotation)=rownames(CorCondWise)
      # colnames(Annotation)="Conditions"
      pheatmap(CorCondWise, annotation =Annotation )
      # rownames(ContrastMap)=rownames(ContrastMap)[match(rownames(ContrastMap), colnames(LMFit_Data))] # reorder the ContrastMap rowname such like it is in the same order as the column of the LMFit_Data  
      if (!is.logical(Plots)){
        dev.off()
      }
    }
    
  }
 if (is.character(Plots) | is.logical(Plots)){
   
   if (!isFALSE(Plots)) {
     if (!is.logical(Plots)){
       svg(filename = file.path(PlotExportPath,"GeneCluster.svg"),width = 20, height = 60, pointsize = 10)
     }
     # Plot the hierarchical clustering of genes. This is to see if some genes cluster with each other in the selected conditions
     SignificantGene <- decideTests(efit_Data,lfc = lfc_cuttoff ,p.value = pvalue_cuttoff, adjust.method = "BH")
     # get the coordinate of the genes differencially expressed
     SignificantGene= as.data.frame(SignificantGene)
     BolSignificantGene= ListTopGenesAmongConditions(SignificantGene)
     Title="Genes Differentially Expressed \nin All Conditions"
     if (length(BolSignificantGene)==100){
       Title="Top 100 Genes Differentially Expressed \namong Conditions"
     }
     lfc_sign=as.data.frame(efit_Data$coefficients)
     
     lfc_sign=lfc_sign[BolSignificantGene,]
     # row_dend =hclust(dist(lfc_sign, method = "minkowski")) # row clustering
     # col_dend = hclust(dist(t(lfc_sign), method = "minkowski"),) # column clustering
     row_dend =dist(lfc_sign, method = "minkowski") # row clustering
     col_dend =dist(t(lfc_sign), method = "minkowski") # column clustering
     # Ref_Meta=Annotation[match(RefCond,Annotation$group_Cond),]
     # Test_Meta=Annotation[match(CompareCond,Annotation$group_Cond),]
     # BolAnnotCol=colSums(Ref_Meta==Test_Meta)==nrow(Ref_Meta)
     # Annotation2=Ref_Meta[,BolAnnotCol]
     # rownames(Annotation2)=colnames(lfc_sign)
     pheatmap(lfc_sign,
              cluster_rows = T,
              annotation_col = Annotation2,
              clustering_distance_rows = row_dend,
              clustering_distance_cols = col_dend,
              clustering_method = "average")
     # pheatmap(lfc_sign, cluster_rows = T,annotation_col = Annotation2, main = Title)
     if (!is.logical(Plots)){
       dev.off()
     }
   }
   
    
   
 }
  # Build PCA of genes. This is to see if some group of genes may clusters together
  if (is.character(Plots) | is.logical(Plots)){
    # if (!is.logical(Plots)){
    #   svg(filename = file.path(PlotExportPath,"PCAs.svg"),width = 20, height = 20, pointsize = 10)
    # }
    if (!isFALSE(Plots)) {
      # if (!is.logical(Plots)){
      #   svg(filename = file.path(PlotExportPath,"PCAs.svg"),width = 20, height = 20, pointsize = 10)
      # }
      # MeanlogFCini=rowMeans(efit_Data$coefficients)
      MeanlogFCini=rowMeans(efit_Data$coefficients)
      Percents=quantile(MeanlogFCini, probs = c(0.005, .9995))
      MeanlogFC=MeanlogFCini
      # MeanlogFC=efit_Data$Amean
      MeanlogFC[MeanlogFCini<Percents[1]]=rep(Percents[1], sum(MeanlogFCini<Percents[1]))
      MeanlogFC[MeanlogFCini>Percents[2]]=rep(Percents[2], sum(MeanlogFCini>Percents[2]))
      
      SignficantivelyDiff=decideTests(efit_Data, p.value = 0.05, lfc = 0.5, adjust.method = "BH")
      SignficantivelyDiff_string=decideTests(tfit, p.value = 0.05, adjust.method = "BH")
      
      # select genes
      BolOnceSignficantivelyDiff=rowSums(abs(SignficantivelyDiff))!=0
      BolAllSignficantivelyDiff=ListTopGenesAmongConditions(SignficantivelyDiff)
      BolLFC_Once_above_1=rowSums(abs(SignficantivelyDiff_string))!=0
      BolLFC_All_above_1=ListTopGenesAmongConditions(SignficantivelyDiff_string)
      
      Title1 = "p-value < 0.05 \n all among conditions"
      Title2= "p-value < 0.05, logFC>1.2  \n all among conditions"
      
      if (length(BolAllSignficantivelyDiff)==100){
        Title1=str_replace(Title1,"all", "Top 100 genes")
      }
      if (length(BolLFC_All_above_1)==100){
        Title2=str_replace(Title2,"all", "Top 100 genes")
      }
      
      PC1=prcomp(efit_Data$coefficients[BolOnceSignficantivelyDiff,])
      PC2=prcomp(efit_Data$coefficients[BolLFC_Once_above_1,])
      PC3=prcomp(efit_Data$coefficients[BolAllSignficantivelyDiff,])
      PC4=prcomp(efit_Data$coefficients[BolLFC_All_above_1,])
      
      
      # P1=ggbiplot(PC1,choices=2:3,groups = MeanlogFC[BolOnceSignficantivelyDiff]) +
      # P1=ggbiplot(PC1,choices=1:2,groups = MeanlogFC[BolOnceSignficantivelyDiff]) +
      PCA_Data=as.data.frame(cbind(PC1$x, MeanlogFC[BolOnceSignficantivelyDiff]))
      colnames(PCA_Data)[ncol(PCA_Data)]="V3"
      P1=ggplot(PCA_Data, aes(x=PC1, y = PC2, color=V3)) +
        geom_point() +
        ThemeZim + 
        scale_color_gradient( name="Average logFC",
                              low="blue", high="red") + ggtitle("p-value < 0.05 \n once among conditions")
      
      # P2=ggbiplot(PC2,choices=2:3,groups = MeanlogFC[BolLFC_Once_above_1]) +
      
      # P2=ggbiplot(PC2,choices=1:2,groups = MeanlogFC[BolLFC_Once_above_1]) +
      PCA_Data=as.data.frame(cbind(PC2$x, MeanlogFC[BolLFC_Once_above_1]))
      colnames(PCA_Data)[ncol(PCA_Data)]="V3"
      P2=ggplot(PCA_Data, aes(x=PC1, y = PC2, color=V3)) +
        geom_point() +
        ThemeZim + 
        scale_color_gradient(  name="Average logFC",
                               low="blue", high="red")+
        ggtitle("p-value < 0.05, logFC>1.2  \n once among conditions")
      
      # P3=ggbiplot(PC3,choices=2:3,groups = MeanlogFC[BolAllSignficantivelyDiff]) +
      # P3=ggbiplot(PC3,choices=1:2,groups = MeanlogFC[BolAllSignficantivelyDiff]) +
        PCA_Data=as.data.frame(cbind(PC3$x, MeanlogFC[BolAllSignficantivelyDiff]))
        colnames(PCA_Data)[ncol(PCA_Data)]="V3"
      P3=ggplot(PCA_Data, aes(x=PC1, y = PC2, color=V3)) +
        geom_point() +
        ThemeZim + 
        scale_color_gradient(   name="Average logFC",
                                low="blue", high="red") + 
        ggtitle(Title1)
      
      # P4=ggbiplot(PC4,choices=2:3,groups = MeanlogFC[BolLFC_All_above_1]) +
      PCA_Data=as.data.frame(cbind(PC4$x, MeanlogFC[BolLFC_All_above_1]))
      colnames(PCA_Data)[ncol(PCA_Data)]="V3"
      P4=ggplot(PCA_Data, aes(x=PC1, y = PC2, color=V3)) +
        geom_point() +
        ThemeZim + 
        scale_color_gradient( name="Average logFC",
                              low="blue", high="red")+
        ggtitle(Title2)
      # ggbiplot(PC4,choices=1:2,groups = MeanlogFC[BolLFC_All_above_1])%>%ggsave(file = file.path(PlotExportPath,"PCAs.svg"),  width=15, height=15)
      
      # grid.arrange(arrangeGrob(P1, P2, ncol=2),arrangeGrob(P3,P4 , ncol=2), nrow = 2)
      if (!is.logical(Plots)){
        grid.arrange(P1, P2,P3,P4 , ncol=2, nrow = 2)%>%ggsave(file = file.path(PlotExportPath,"PCAs.svg"),  width=15, height=15)
      }
      else{
        grid.arrange(arrangeGrob(P1, P2, ncol=2),arrangeGrob(P3,P4 , ncol=2), nrow = 2)
      }
    }
    
    
  }
  
  Output$MetaData=Annotation2
  Output$MetaLM=Annotation
  Output[["eBayes_Stat"]]=efit_Data
  Output[["treat_Stat"]]=tfit
  Output[["eBayes_Summary"]]=Summary_DiffReg_genes
  Output[["treat_summary"]]=summary(decideTests(tfit,p.value = 0.05, adjust.method = "BH"))
  Output[["LM_Fit"]]=LMFit_Data
  return(Output)
}

BuildContrastMap = function(Sample,Ref,  Method = "c" ) {
  # used to build a contrast map further used for deferentially expressed genes. The basic principle consist in is to build a list of string that will describe which condition will be compared.
  # for a detailed description of argument see the "Genes_Diff_Expressed" function
  Ref=make.names(Ref)
  Sample=make.names(Sample)
  if (str_detect(  Method, "aired")) {
    ComparisonStr= str_c(Ref, Sample, sep = "-")
  }
  else{
    ComparisonStr=vector()
    for (i in Ref) {
      Temp= str_c(rep(i,length(Sample)),Sample, sep = "-")
      ComparisonStr=c(ComparisonStr, Temp)
    }
  }
  UniqueCompa=as.character(unique(ComparisonStr))
  # print(UniqueCompa)
  Lev=unique(as.character(c(Ref,Sample)))
  # print(Lev)
  Contrast=makeContrasts(contrasts = UniqueCompa, levels = Lev)
  
  Contrast=Contrast[,colSums(Contrast==0)!=nrow(Contrast)]
  if (is.null(nrow(Contrast))){Contrast=makeContrasts(contrasts = UniqueCompa, levels = Lev)}
  return(Contrast)
}

MergeDEResults=function(HashObj,VarName="Gene", Method="c"){
  Output=data.frame()
  K=keys(HashObj)
  K=K[K!="MetaData"]
  if (Method=="c"){
    GeneInfo=HashObj[[K[1]]][, 1:(ncol(HashObj[[K[1]]])-6)]
    Output=GeneInfo
    Cond=data.frame(rep(NA,nrow(Output)))
    rownames(Cond)=Output[,VarName]
    rownames(Output)=Output[,VarName]
    for (i in K){
      colnames(Cond)=i
      tbl=HashObj[[i]][, (ncol(HashObj[[i]])-5):ncol(HashObj[[i]])]
      Match=match(Output[,VarName],HashObj[[i]][,VarName])
      rownames(tbl)=Output[,VarName]
      Output=cbind(Output,Cond,tbl[Match,])
      # Output=cbind(Output,Cond,tbl)
    }
  }
  else{
    for (i in K){
      tbl=HashObj[[i]]
      Cond=data.frame(rep(i,nrow(tbl)))
      Output=rbind(Output,cbind(Cond,tbl))
      # Output=cbind(Output,Cond,tbl)
    }
  }
  Output=data.frame(Output, stringsAsFactors = T)
  colnames(Output)[1]="Condition"
  return(Output)
}

MergeDEResSpecifiComp=function(DGE_Res,Condition,RefConditions, CommonVar="Gene", OutputFolderPath=""){
  Output=hash()
  Renamed_selected_cond=hash()
  Metadata=DGE_Res$MetaData
  Metadata=sapply(Metadata, as.character)
  Metadata=as.data.frame(Metadata)
  ind=1
  for (i in Condition){
    bolRefCond=Metadata[[i]]==RefConditions[ind]
    # BolOtherCondition=colnames(Metadata)!=i
    for (j in levels(Metadata[[i]])[levels(Metadata[[i]])!=RefConditions[ind]]){
      bolTestCond=Metadata[[i]]==j
      Match=match(paste0(Metadata[bolTestCond,i],collapse = ""),
                  paste0(Metadata[bolRefCond,i],collapse = ""))
      # print(i)
      # print(paste(Metadata$group_Cond[bolTestCond],
      #              Metadata$group_Cond[bolRefCond],sep = "-"))
      SelectedComp=paste(Metadata$group_Cond[bolTestCond],
                         Metadata$group_Cond[bolRefCond],sep = "-")
      SelectedCompName=paste(Metadata$group_Cond[bolRefCond],
                             Metadata$group_Cond[bolTestCond],sep = "_Vs_")
      
      SelectDGE=DGE_Res[SelectedComp]
      for (k in 1:length(SelectedComp)){
        Renamed_selected_cond[SelectedCompName[k]]=SelectDGE[[SelectedComp[k]]]
      }
      Output[[i]]=MergeDEResults(Renamed_selected_cond,CommonVar)
    }
    ind=ind+1
  }
return(Output)
  # droplevels
}

BuildDEComparaison=function(Meta, UniqueVar, ReferencesVal){
  Output=data.frame()
  Renamed_selected_cond=hash()
  
  Metadata=sapply(Meta, as.character)
  Metadata=as.data.frame(Metadata)
  CondVar=colnames(Meta)[!str_detect(colnames(Meta),UniqueVar)]
  ind=1
  for (i in CondVar){
    bolRefCond=Metadata[,i]==ReferencesVal[ind]
    # BolOtherCondition=colnames(Metadata)!=i
    for (j in unique(Metadata[,i])[unique(Metadata[,i])!=ReferencesVal[ind]]){
      bolTestCond=Metadata[,i]==j
      
      # associate the other conditions in a single string and find row with similar pattern among ref and Test
      RefString=unite(as.data.frame(Metadata[,CondVar[CondVar!=i]]),"Concat",sep = "")
      Match=match(RefString$Concat[bolTestCond],
                  RefString$Concat[bolRefCond])
      Match=Match[!is.na(Match)]
      # print(i)
      # print(paste(Metadata[bolTestCond,UniqueVar],
                   # Metadata[bolRefCond,UniqueVar][Match],sep = "-"))
      
      Output=rbind(Output,
                   data.frame(unite(data.frame(Metadata[bolTestCond,UniqueVar],
                                               Metadata[bolRefCond,UniqueVar][Match]), "Name",sep = "_-_")
                     ,Metadata[bolTestCond,UniqueVar],
                                Metadata[bolRefCond,UniqueVar][Match],
                         rep(i, length(Match))))
    }
    ind=ind+1
  }
  colnames(Output)[2:ncol(Output)]=c("Test","Ref","Test_Cond")
  return(Output)
}


ListTopGenesAmongConditions=function(decideTestTbl, n=100,force.Top.Ranking=F){
  Bol=rowSums(abs(decideTestTbl))==ncol(decideTestTbl)
  if (sum(Bol)<5 | force.Top.Ranking==T){
    SumDiff=rowSums(abs(decideTestTbl))
    SortedSumDiff=sort(SumDiff,decreasing = T)
    Bol=names(SortedSumDiff)[1:n]
    print(paste("Top", as.character(n),  "Genes Differentially Expressed \nin",
                as.character(round(mean(SortedSumDiff[Bol]))) ,
                "Conditions"))
    Output=Bol
  }
  else{
    Output=rownames(decideTestTbl)[Bol]
  }
  return(Output)
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

BuildCorrelationMatrix=function(table, method="pearson", plot=T){
  # CorMat=cor(table, method =method)
  CorMat=rcorr(as.matrix(table) ,type = method)
  FlattenCorrMatrix=flattenCorrMatrix(CorMat$r, CorMat$P)
  Output=list(
    CorrMat=CorMat,
    Summary=FlattenCorrMatrix
  )
  if (plot){
    P=pheatmap(CorMat$r,cluster_cols = T,
               clustering_distance_rows="canberra", clustering_distance_cols ="canberra" ,clustering_method = "average", main = "Correlation Matrix",cutree_rows = 20)
    Output[["HeatMap"]]=P
  }
  return(Output)
}


# 
# Meta=RawData$MetaData
# Meta2=Meta[!duplicated(Meta$Conditions),]
# Cond=Meta2[,c("Conditions","Isolate","condition","time")]
# Ref=c("CEC3621", "FP", 0)
# Sample="Conditions"
# 
# Map=BuildDEComparaison(Cond,Sample, Ref)
# 
# DETest=Genes_Diff_Expressed(NormData$Norm_Data,NormData$Norm_Data$targets$group_Cond,Group4InteractivePlot = data.frame(NormData$Norm_Data$targets$time, NormData$Norm_Data$targets$Isolate), Plots = file.path(OutputFolderPath,"CondIsolate"), AdditionalFactor=c("condition", "time", "Isolate"), RefCond = Map$Ref[Map$Test_Cond=="Isolate"],CompareCond = Map$Test[Map$Test_Cond=="Isolate"],Method = "Paired" )

# row_dend = hclust(dist(df)) # row clustering
# col_dend = hclust(dist(t(df))) # column clustering
# Heatmap(df, name = "mtcars", 
#         row_names_gp = gpar(fontsize = 6.5),
#         cluster_rows = color_branches(row_dend, k = 4),
#         cluster_columns = color_branches(col_dend, k = 2))

