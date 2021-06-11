#  this set of function enable to perform GO analysis using TopGO library. The manual for this package can be found at https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

# The topGO package is designed to facilitate semi-automated enrichment analysis for Gene Ontology (GO) terms.   The  process  consists  of  input  of  normalised  gene  expression  measurements,  gene-wise  correlationor differential expression analysis, enrichment analysis of GO terms, interpretation and visualisation of the results. One of the main advantages of topGO is the unified gene set testing framework it offers.  Besides providing an  easy  to  use  set  of  functions  for  performing  GO  enrichment  analysis,  it  also  enables  the  user  to  easily implement  new  statistical  tests  or  new  algorithms  that  deal  with  the  GO  graph  structure.   This  unified framework also facilitates the comparison between different GO enrichment methodologies. There are a number of test statistics and algorithms dealing with the GO graph structured ready to use in topGO. The elim and weight algorithms were introduced in Alexa et al. (2006).  The default algorithm used by the topGO package is a mixture between the elim and theweightalgorithms and it will be referred asweight01.TheparentChildalgorithm was introduced by Grossmann et al. (2007).We assume the user has a good understanding of GO, see Consortium (2001), and is familiar with gene setenrichment tests.  Also this document requires basic knowledge of R language.The  next  section  presents  a  quick  tour  intotopGOand  is  thought  to  be  independent  of  the  rest  of  thismanuscript.  The remaining sections provide details on the functions used in the sample section as well as showing more advance functionality implemented in the topGO package. 

library(topGO)
library(R.utils)
# library(mgsa)
# library(ViSEAGO)
library(STRINGdb)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
# Choose the destination Folder for download and db build 

Dest="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/candida_ref"

BuildCandidaGOdb=function(SaveDir=Dest){
  # Gene universe file http://www.candidagenome.org/download/go/gene_association.cgd.gz
  Path2GOannotationFile="http://www.candidagenome.org/download/go/gene_association.cgd.gz"
  
  # Define a file name for the gzip file
  gzipFileName="gene_association.cgd.gz"
  FileName="gene_association.cgd"
  
  # Download gzip file from the Candida db site
  download.file(Path2GOannotationFile,destfile = file.path(SaveDir, gzipFileName))
  
  #Unzip the file containing associations
  gunzip(file.path(SaveDir, gzipFileName), overwrite=T)
  
  #  import annotation table
  UnzipFile=filePath(SaveDir,FileName)
  ImportTable=read.delim(UnzipFile, skip = 21, sep = "\t",header = F)
  # remove the last 2 column
  ImportTable=ImportTable[,1:(ncol(ImportTable)-2)]
  
  # Columns are:					Contents:
  #   
  #   1) DB						- database contributing the file (always "CGD" for this file)
  # 2) DB_Object_ID				- CGDID
  # 3) DB_Object_Symbol				- see below
  # 4) Qualifier 			(optional)	- 'NOT', 'contributes_to', or 'colocalizes_with' qualifier 
  # for a GO annotation, when needed
  # 5) GO ID					- unique numeric identifier for the GO term
  # 6) DB:Reference(|DB:Reference)			- the reference associated with the GO annotation
  # 7) Evidence					- the evidence code for the GO annotation
  # 8) With (or) From 		(optional)	- any With or From qualifier for the GO annotation
  # 9) Aspect					- which ontology the GO term belongs in (see note below)
  # 10) DB_Object_Name(|Name) 	(optional)	- a name for the gene product in words, e.g. 'acid phosphatase'
  # 11) DB_Object_Synonym(|Synonym) (optional)	- see below
  # 12) DB_Object_Type				- type of object annotated, e.g. gene, protein, etc.
  # 13) taxon(|taxon)				- taxonomic identifier of species encoding gene product
  # 14) Date					- date GO annotation was made
  # 15) Assigned_by					- source of the annotation
  ColNames=c("DB","Ca_db_id","gene_id","Qualifier","GOID","DB_Reference","evidence","With","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxid","Date","Assigned_by")
  colnames(ImportTable)=ColNames
  
  # Choose a taxon based on the Taxon ID 
  # ! Organism 1: Candida albicans SC5314, genome version: A22-s07-m01-r124, taxon: 237561
  # ! Organism 2: Candida glabrata CBS138, genome version: s03-m01-r09, taxon: 284593
  # ! Organism 3: Candida auris B8441, genome version: s01-m01-r12, taxon: 498019
  # ! Organism 4: Candida dubliniensis CD36, genome version: s01-m02-r29, taxon: 573826
  # ! Organism 5: Candida parapsilosis CDC317, genome version: s01-m03-r46, taxon: 578454
  # ! Organism 6: Candida orthopsilosis Co 90-125, taxon: 1136231
  # ! Organism 7: Candida albicans WO-1, taxon: 294748
  # ! Organism 8: Debaryomyces hansenii CBS767, taxon: 284592
  # ! Organism 9: Lodderomyces elongisporus NRLL YB-4239, taxon: 379508
  # ! Organism 10: Candida tropicalis MYA-3404, taxon: 294747
  # ! Organism 11: Candida lusitaniae CBS6936, taxon: 36911
  # ! Organism 12: Candida lusitaniae ATCC 42720, taxon: 306902
  # ! Organism 13: Candida auris B11221, taxon: 498019
  # ! Organism 14: Candida guilliermondii ATCC 6260, taxon: 294746
  # ! Date created: Wed Nov 18 01:21:21 2020
  
# reformate the table
  ImportTable$taxid=str_remove(ImportTable$taxid,"taxon:")
  ImportTable$taxid=as.numeric(as.character(ImportTable$taxid))
  # Bol=str_detect(ImportTable$taxid,as.character("Taxon"))
  # keep line containing GO ID that are known in the GO.db
  BolGO=ImportTable$GOID %in% keys(GO.db)
  BolDup=duplicated(ImportTable[BolGO,c("taxid","gene_id","GOID","evidence")])
  ImportTable=ImportTable[BolGO, c("taxid","gene_id","GOID","evidence")]
  
  # write the table to a new file
  write.table(ImportTable[!BolDup,],file.path(SaveDir,"Ca_GO_db.tsv"),sep = "\t", row.names = F)
  return(file.path(SaveDir,"Ca_GO_db.tsv"))
}

# BuildCandidaGOdb()

# BuildDBUsingAnnotationForge

# Making use of makeOrgPackage()

# Sometimes you may not find what you need at NCBI, most commonly this is because they may just not have enough data about the organism you are interested in. But often other resources will have annotation data that you want to make into an organism package. When this happens you can use the much more general makeOrgPackage() function. This function takes more arguments, but it does not rely on NCBI in order to run. We do however still ask for you to provice a tax ID for the metadata (even though we are not using it to look up data from NCBI). Many of the other arguments are also the same as the makeOrgPackageFromNCBI() function. But a key difference is that the 1st argument for this function is (…). For that argument, we want you to provide named arguments corresponding to data.frames of data. Each named argument will become a table name in the resulting database, and each field name (for the data.frames) will become the field names of the database as well as the names looked up by the columns() and keytypes() methods. With the exception of any table that is named by the goTable argument (more on this below), there are not too many restrictions on what kind of data you can put into the data.frame. But one rule you must follow is that the 1st collumn of each data.frame has to correspond to a central gene ID and be labeled “GID”.
# 
# Finally, the goTable method is also new. That argument indicates when one of the data.frames contains GO information. If you choose to use this argument, makeOrgPackage() will post-process your GO data to 1) remove IDs that are too new and 2) create a second table to also represent the GOALL, EVIDENCEALL and ONTOLOGYALL fields for the select method etc. However to use the goTable argument, you have to follow a strict convention with the data. Such a data.frame must have three columns only and these must correspond to the gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE” Below is an example that parses an example file into three data.frame and that makes use of the goTable argument.
BuildCaGOdb_CP=function(ChromosomalFeatureFile="Data/20200528-C_albicans_SC5314_A22_current_chromosomal_feature.tab",
                     GeneAssociationFile="Data/20201123-CGD_gene_association.cgd",
                     NCBIGeneList="Data/20201107-NCBI_GeneList.txt",
                     OuputDir="Data", TaxonID="237561"){
  # Path2Files="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/candida_ref"
  FChromosomalfeatures=file.path(ChromosomalFeatureFile)
  FGeneAssociation=GeneAssociationFile
  
  TblNCBI=read.delim(NCBIGeneList, header = T, sep = "\t")
  
  # Feature name (mandatory); this is the primary systematic name, if available
  # 2.  Gene name (locus name)
  # 3.  Aliases (multiples separated by |)
  # 4.  Feature type
  # 5.  Chromosome
  # 6.  Start Coordinate
  # 7.  Stop Coordinate
  # 8.  Strand 
  # 9.  Primary CGDID
  # 10. Secondary CGDID (if any)
  # 11. Description
  # 12. Date Created
  # 13. Sequence Coordinate Version Date (if any)
  # 14. Blank
  # 15. Blank
  # 16. Date of gene name reservation (if any).
  # 17. Has the reserved gene name become the standard name? (Y/N)
  # 18. Name of S. cerevisiae ortholog(s) (multiples separated by |)
  Ca_ChromosomalFeatures=read.delim(FChromosomalfeatures, skip = 8, header = F, col.names = c("FeatureName","Gene_name", "Alias", "Feature_Type","Chromosome","Start","End","Strand","Primary_ID", "Secondary_ID", "Description", "Date_Created","Seq.Coor.Version_Date", "Blank1", "Blank2", "Date_Name_Reservation","ReserverdName","Sc_ortholog" ), stringsAsFactors = F )
  
  GeneAssociation=read.delim(FGeneAssociation, skip = 21, sep = "\t",header = F , stringsAsFactors = F)
  # remove the last 2 column
  GeneAssociation=GeneAssociation[,1:(ncol(GeneAssociation)-2)]
  # Add colnames
  colnames(GeneAssociation)=c("DB","Ca_db_id","gene_id","Qualifier","GOID","DB_Reference","evidence","With","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxid","Date","Assigned_by")
  # include Associations for the desired taxon
  GeneAssociation=GeneAssociation[str_detect(GeneAssociation$taxid,TaxonID),]
  # Include the Association corresponding to the taxon wanted
  # GID, Symbole, GeneName
  GeneList=unique(GeneAssociation$Ca_db_id)
  
  Ca_ChromosomalFeatures=Ca_ChromosomalFeatures[Ca_ChromosomalFeatures$Primary_ID%in%GeneList,]
  GeneAssociation=GeneAssociation[GeneAssociation$Ca_db_id%in%GeneList & str_detect(GeneAssociation$taxid, TaxonID),]
  
  Ca_ChromosomalFeatures$NCBI_Key=str_remove_all(Ca_ChromosomalFeatures$FeatureName, "_")
  Ca_ChromosomalFeatures$NCBI_Key=paste0("CAALFM_", Ca_ChromosomalFeatures$NCBI_Key)
  Ca_ChromosomalFeatures=merge(Ca_ChromosomalFeatures, TblNCBI[, c("Aliases", "GeneID")], by.x="NCBI_Key",by.y="Aliases")
  # tbl1=Ca_ChromosomalFeatures[, c("Primary_ID", "FeatureName", "Gene_name", "Description")]
  tbl1=Ca_ChromosomalFeatures[, c("GeneID", "FeatureName", "Gene_name", "Description")]
  tbl1$Gene_name[tbl1$Gene_name==""]=as.character(tbl1$Gene_name[tbl1$Gene_name==""])
  tbl1$Gene_name[tbl1$Gene_name==""]=as.character(tbl1$FeatureName[tbl1$Gene_name==""])
  tbl1=tbl1[,c(1,3,4)]
  colnames(tbl1)=c("GID","SYMBOL","GENENAME")
  # colnames(tbl1)=c("ENTREZID","SYMBOL","GENENAME")
  
  tbl2=Ca_ChromosomalFeatures[, c("GeneID", "Chromosome")]
  colnames(tbl2)=c("GID","CHROMOSOME")
  # colnames(tbl2)=c("ENTREZID","CHROMOSOME")
  
  GeneAssociation=merge(GeneAssociation, Ca_ChromosomalFeatures[,c("Primary_ID", "GeneID")], by.x="Ca_db_id", by.y="Primary_ID")
  
  tbl3=GeneAssociation[!duplicated(GeneAssociation[, c("GeneID", "GOID", "evidence")]), c("GeneID", "GOID", "evidence")]
  colnames(tbl3)=c("GID","GO","EVIDENCE")
  # colnames(tbl3)=c("ENTREZID","GO","EVIDENCE")
  
  Output=AnnotationForge::makeOrgPackage(gene_info=tbl1, chromosome=tbl2, go=tbl3,
                                  version="0.1",
                                  maintainer="Eric Durandau <eric.durandau@chuv.ch>",
                                  author="Eric Durandau <eric.durandau@chuv.ch>",
                                  outputDir = OuputDir,
                                  tax_id="237561",
                                  genus="Candida",
                                  species="albicans",
                                  goTable="go" )
  return(Output)
}
# makeOrgPackageFromNCBI(version = "0.1",
#                        author = "Eric Durandau <eric.durandau@chuv.ch>",
#                        maintainer = "Eric Durandau <eric.durandau@chuv.ch>",
#                        outputDir = "/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/Databases/AnnotationForge",
#                        tax_id = "5476",
#                        genus = "Candida",
#                        species = "albicans")

# BuildCaGOdb_CP()
# install.packages("/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/Databases/AnnotationForge/org.Calbicans.eg.db/",repos = NULL)

# ##### Section relative to stringdb
# define the taxon
Taxon=5476
# Load the stringdb object
string_db_Ca <- STRINGdb$new( version="11", species=Taxon,score_threshold=200, input_directory="../Data/STRINGdb")

string_db_Cg <- STRINGdb$new( version="11", species=5478,score_threshold=200, input_directory="../Data/STRINGdb")

string_db_Hs <- STRINGdb$new( version="11", species=9606,score_threshold=200, input_directory="../Data/STRINGdb")

PlotSTRINGdbNetwork=function(Table,Cluster="all", VarCluster="HC_Clusters", StringIDVar="STRING_id", AliasVar="Gene",MaxConnection=150, string_db=string_db_Ca, Plot=T){
  RemoveRedundantConnexion=function(Table){
    ConcatInt=paste0(Table[,1], Table[,2])
    Table=Table[!duplicated(ConcatInt),]
    return(Table)
  }
  
  if(Cluster=="all"){
    GenesSetStringID=Table[,StringIDVar]
    GenesSetAlias=Table[,AliasVar]
  }
  else{
    Bol=Table[,VarCluster]==Cluster
    GenesSetStringID=Table[Bol,StringIDVar]
    GenesSetAlias=Table[Bol,AliasVar]
  }
  
  # Get interaction for the selected gene set
  Interaction=string_db$get_interactions(GenesSetStringID)
  # remove redundant connexion
  Interaction=RemoveRedundantConnexion(Interaction)
  # Conserve connexion including the geneSet
  Bol= Interaction[,1]%in%GenesSetStringID & Interaction[,2]%in%GenesSetStringID
  Interaction=Interaction[Bol,]
  
  # Add neighbor if not enough connexion
  if (nrow(Interaction)<MaxConnection){
    # calculate the missing number of connection
    MissingConnexion=MaxConnection-nrow(Interaction)
    # get neighborhood
    UniqueID=string_db$get_neighbors(GenesSetStringID)
    Concatint=string_db$get_interactions(c(UniqueID,GenesSetStringID))
    Concatint=RemoveRedundantConnexion(Concatint)
    CurrentConnection=nrow(Concatint)
    # get connection between neighborhood
    IntraConnect=string_db$get_interactions(UniqueID)
    IntraConnect=RemoveRedundantConnexion(IntraConnect)
    IntraConnect=IntraConnect[order(IntraConnect$combined_score,decreasing = T),]
    PrevConnec=CurrentConnection+1
    Ind=1
    while(CurrentConnection>MaxConnection && Ind<200 ){
      PrevConnec=CurrentConnection
      # Get coordinate value above 5 percent
      IntraConnect=IntraConnect[1:round(nrow(IntraConnect)*0.95),]
      # Get unique ID of extra Genes
      UniqueID=unique(c(IntraConnect[, 1], IntraConnect[, 2]))
      # Rebuild connection with the initial dataset
      Concatint=string_db$get_interactions(c(UniqueID,GenesSetStringID))
      Concatint=RemoveRedundantConnexion(Concatint)
      CurrentConnection=nrow(Concatint)
      # print(CurrentConnection)
      # print(Ind)
      # print(length(UniqueID))
      Ind=Ind+1
    }
    # CurrentConnection
    Title=Title=paste(as.character(MaxConnection), "STRINGdb gene connection with extension")
  }
  else {
    Concatint=Interaction
    # restrict the connection to the maximum connection number
    Concatint=Concatint[order(Concatint$combined_score,decreasing = T),]
    Concatint=Concatint[1:MaxConnection,]
    Title=Title=paste("Best", as.character(MaxConnection), "STRINGdb gene connection")
  }
  
  
  Interaction=Concatint
  
  Coord= match(Interaction$from, Table[,StringIDVar])
  Coord2= match(Interaction$to, Table[,StringIDVar])
  
  Interaction$from=Table[Coord,AliasVar]
  Interaction$to=Table[Coord2,AliasVar]
  
  Interaction=Interaction[!is.na(Interaction$from) & !is.na(Interaction$to),]
  
  colnames(Interaction)=c("source","target","importance")
  Interaction$importance=Interaction$importance/100
  

  Counts=as.data.table(table(c(Interaction$source, Interaction$target)))
  Counts$Colors="gray"
  Counts$Colors[Counts$V1%in%GenesSetAlias]="chartreuse"
  setnames(Counts,c("source", "Count","Colors"))
  Counts$FontColors=Counts$Colors
  Counts$FontColors[Counts$Colors=="gray"]="black"
  Counts$FontColors[!Counts$Colors=="gray"]="darkgreen"

  
  network <- igraph::graph.data.frame(d=Interaction, vertices = Counts, directed=F)
  
  Size=((Counts$Count-min(Counts$Count)+1)/(max(Counts$Count)-min(Counts$Count)+1))*7
  Importance=edge.width=igraph::E(network)$importance
  
  l <-igraph::layout_with_fr(network)
  
  par(xpd = T, mar = par()$mar + c(0,0,0,7))
  igraph::V(network)$label.cex=0.8
  igraph::V(network)$label.font[igraph::V(network)$Colors=="chartreuse"]=2
  
  SeqN=7
  if (!(length(unique(igraph::E(network)$importance)) > SeqN)){
    SeqN=length(unique(igraph::E(network)$importance))
  }
  ImportanceSeq=seq(1, length(igraph::E(network)$importance), length.out = SeqN)
  
  if (Plot==T){
    plot(network,
         Layout=l, 
         edge.width=((Importance-min(Importance))/max(Importance-min(Importance)))*5+0.5,
         vertex.color=igraph::V(network)$Colors,
         vertex.label.font=igraph::V(network)$label.font,
         vertex.label.degree=-pi/2,
         vertex.label.dist=0.9,
         vertex.label.color=igraph::V(network)$FontColors,
         vertex.size=Size, 
         edge.curved=.2,
         main= Title)
    # ggplot(ggnetwork::fortify(network) )
    # Add a legend
    legend(1.2, 1.2, title =  "Connections",legend=as.character(levels(as.factor(igraph::V(network)$Count))) , bty = "n", pch=20 , pt.cex = sort(unique(Size))/1.7, cex = 0.7, horiz = FALSE)
    
    legend(1.2, 0.1, title =  "Type",legend=c("Gene from cluster", "Extended connection") , bty = "n", pch=20 , cex = 0.7, horiz = FALSE, col=c("chartreuse", "gray"),pt.cex =2)
    
    legend(1.2, 0.6, title =  "Link Confidence",legend=as.character(levels(as.factor(igraph::E(network)$importance*100))[ImportanceSeq])  , lwd = sort(unique(((Importance-min(Importance))/max(Importance-min(Importance)))*7+0.5))[ImportanceSeq], cex = 0.7, horiz = FALSE, bty = "n")
    
    par(mar=c(4, 3, 1, 1 + 0.1))
  }else{
    return(Interaction) 
  }
}

BuildCandidaGOdb=function(SaveDir=Dest){
  # Gene universe file http://www.candidagenome.org/download/go/gene_association.cgd.gz
  Path2GOannotationFile="http://www.candidagenome.org/download/go/gene_association.cgd.gz"
  
  # Define a file name for the gzip file
  gzipFileName="gene_association.cgd.gz"
  FileName="gene_association.cgd"
  
  # Download gzip file from the Candida db site
  download.file(Path2GOannotationFile,destfile = file.path(SaveDir, gzipFileName))
  
  #Unzip the file containing associations
  gunzip(file.path(SaveDir, gzipFileName), overwrite=T)
  
  #  import annotation table
  UnzipFile=filePath(SaveDir,FileName)
  ImportTable=read.delim(UnzipFile, skip = 21, sep = "\t",header = F)
  # remove the last 2 column
  ImportTable=ImportTable[,1:(ncol(ImportTable)-2)]
  
  # Columns are:					Contents:
  #   
  #   1) DB						- database contributing the file (always "CGD" for this file)
  # 2) DB_Object_ID				- CGDID
  # 3) DB_Object_Symbol				- see below
  # 4) Qualifier 			(optional)	- 'NOT', 'contributes_to', or 'colocalizes_with' qualifier 
  # for a GO annotation, when needed
  # 5) GO ID					- unique numeric identifier for the GO term
  # 6) DB:Reference(|DB:Reference)			- the reference associated with the GO annotation
  # 7) Evidence					- the evidence code for the GO annotation
  # 8) With (or) From 		(optional)	- any With or From qualifier for the GO annotation
  # 9) Aspect					- which ontology the GO term belongs in (see note below)
  # 10) DB_Object_Name(|Name) 	(optional)	- a name for the gene product in words, e.g. 'acid phosphatase'
  # 11) DB_Object_Synonym(|Synonym) (optional)	- see below
  # 12) DB_Object_Type				- type of object annotated, e.g. gene, protein, etc.
  # 13) taxon(|taxon)				- taxonomic identifier of species encoding gene product
  # 14) Date					- date GO annotation was made
  # 15) Assigned_by					- source of the annotation
  ColNames=c("DB","Ca_db_id","gene_id","Qualifier","GOID","DB_Reference","evidence","With","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxid","Date","Assigned_by")
  colnames(ImportTable)=ColNames
  
  # Choose a taxon based on the Taxon ID 
  # ! Organism 1: Candida albicans SC5314, genome version: A22-s07-m01-r124, taxon: 237561
  # ! Organism 2: Candida glabrata CBS138, genome version: s03-m01-r09, taxon: 284593
  # ! Organism 3: Candida auris B8441, genome version: s01-m01-r12, taxon: 498019
  # ! Organism 4: Candida dubliniensis CD36, genome version: s01-m02-r29, taxon: 573826
  # ! Organism 5: Candida parapsilosis CDC317, genome version: s01-m03-r46, taxon: 578454
  # ! Organism 6: Candida orthopsilosis Co 90-125, taxon: 1136231
  # ! Organism 7: Candida albicans WO-1, taxon: 294748
  # ! Organism 8: Debaryomyces hansenii CBS767, taxon: 284592
  # ! Organism 9: Lodderomyces elongisporus NRLL YB-4239, taxon: 379508
  # ! Organism 10: Candida tropicalis MYA-3404, taxon: 294747
  # ! Organism 11: Candida lusitaniae CBS6936, taxon: 36911
  # ! Organism 12: Candida lusitaniae ATCC 42720, taxon: 306902
  # ! Organism 13: Candida auris B11221, taxon: 498019
  # ! Organism 14: Candida guilliermondii ATCC 6260, taxon: 294746
  # ! Date created: Wed Nov 18 01:21:21 2020
  
  # reformate the table
  ImportTable$taxid=str_remove(ImportTable$taxid,"taxon:")
  ImportTable$taxid=as.numeric(as.character(ImportTable$taxid))
  # Bol=str_detect(ImportTable$taxid,as.character("Taxon"))
  # keep line containing GO ID that are known in the GO.db
  BolGO=ImportTable$GOID %in% keys(GO.db)
  BolDup=duplicated(ImportTable[BolGO,c("taxid","gene_id","GOID","evidence")])
  ImportTable=ImportTable[BolGO, c("taxid","gene_id","GOID","evidence")]
  
  # write the table to a new file
  write.table(ImportTable[!BolDup,],file.path(SaveDir,"Ca_GO_db.tsv"),sep = "\t", row.names = F)
  return(file.path(SaveDir,"Ca_GO_db.tsv"))
}

# BuildCandidaGOdb()

# BuildDBUsingAnnotationForge

# Making use of makeOrgPackage()

# Sometimes you may not find what you need at NCBI, most commonly this is because they may just not have enough data about the organism you are interested in. But often other resources will have annotation data that you want to make into an organism package. When this happens you can use the much more general makeOrgPackage() function. This function takes more arguments, but it does not rely on NCBI in order to run. We do however still ask for you to provice a tax ID for the metadata (even though we are not using it to look up data from NCBI). Many of the other arguments are also the same as the makeOrgPackageFromNCBI() function. But a key difference is that the 1st argument for this function is (…). For that argument, we want you to provide named arguments corresponding to data.frames of data. Each named argument will become a table name in the resulting database, and each field name (for the data.frames) will become the field names of the database as well as the names looked up by the columns() and keytypes() methods. With the exception of any table that is named by the goTable argument (more on this below), there are not too many restrictions on what kind of data you can put into the data.frame. But one rule you must follow is that the 1st collumn of each data.frame has to correspond to a central gene ID and be labeled “GID”.
# 
# Finally, the goTable method is also new. That argument indicates when one of the data.frames contains GO information. If you choose to use this argument, makeOrgPackage() will post-process your GO data to 1) remove IDs that are too new and 2) create a second table to also represent the GOALL, EVIDENCEALL and ONTOLOGYALL fields for the select method etc. However to use the goTable argument, you have to follow a strict convention with the data. Such a data.frame must have three columns only and these must correspond to the gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE” Below is an example that parses an example file into three data.frame and that makes use of the goTable argument.
BuildCgGOdb_CP=function(ChromosomalFeatureFile="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/C_glabrata/Assembly_mapping.tab",
                        GeneAssociationFile="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/CGD_gene_association.cgd",
                        OuputDir="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/Databases/AnnotationForge", TaxonID="284593"){
  # Path2Files="/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/candida_ref"
  FChromosomalfeatures=file.path(ChromosomalFeatureFile)
  FGeneAssociation=GeneAssociationFile
  
  # Feature name (mandatory); this is the primary systematic name, if available
  # 2.  Gene name (locus name)
  # 3.  Aliases (multiples separated by |)
  # 4.  Feature type
  # 5.  Chromosome
  # 6.  Start Coordinate
  # 7.  Stop Coordinate
  # 8.  Strand 
  # 9.  Primary CGDID
  # 10. Secondary CGDID (if any)
  # 11. Description
  # 12. Date Created
  # 13. Sequence Coordinate Version Date (if any)
  # 14. Blank
  # 15. Blank
  # 16. Date of gene name reservation (if any).
  # 17. Has the reserved gene name become the standard name? (Y/N)
  # 18. Name of S. cerevisiae ortholog(s) (multiples separated by |)
  
  Ca_ChromosomalFeatures=read.delim(FChromosomalfeatures, header = T, stringsAsFactors = F )
  
  # Recover the ID of 
  
  GeneAssociation=read.delim(FGeneAssociation, skip = 21, sep = "\t",header = F , stringsAsFactors = F)
  # remove the last 2 column
  GeneAssociation=GeneAssociation[,1:(ncol(GeneAssociation)-2)]
  # Add colnames
  colnames(GeneAssociation)=c("DB","Ca_db_id","gene_id","Qualifier","GOID","DB_Reference","evidence","With","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxid","Date","Assigned_by")
  # include Associations for the desired taxon
  GeneAssociation=GeneAssociation[str_detect(GeneAssociation$taxid,TaxonID),]
  # Include the Association corresponding to the taxon wanted
  # GID, Symbole, GeneName
  GeneAssociation$CBS138_ID=str_extract(GeneAssociation$DB_Object_Synonym, "(CAGL\\w+|Cagl\\w+)")
  GeneList=unique(GeneAssociation$CBS138_ID)
  Ca_ChromosomalFeatures=Ca_ChromosomalFeatures[Ca_ChromosomalFeatures$locus_tag.CBS138%in%GeneList,]
  Ca_ChromosomalFeatures$Name[str_detect(Ca_ChromosomalFeatures$Name, "^CDS$")]=Ca_ChromosomalFeatures$locus_tag.CBS138[str_detect(Ca_ChromosomalFeatures$Name, "^CDS$")]
  GeneAssociation=GeneAssociation[GeneAssociation$CBS138_ID%in%GeneList & str_detect(GeneAssociation$taxid, TaxonID),]
  
  tbl1=Ca_ChromosomalFeatures[, c("locus_tag.CBS138", "locus_tag.DSY562", "Name", "note.CBS138")]
  # tbl1$Gene_name[tbl1$Gene_name==""]=as.character(tbl1$Gene_name[tbl1$Gene_name==""])
  # tbl1$Gene_name[tbl1$Gene_name==""]=as.character(tbl1$FeatureName[tbl1$Gene_name==""])
  tbl1=tbl1[,c(1,3,3)]
  colnames(tbl1)=c("GID","SYMBOL","GENENAME")
  # colnames(tbl1)=c("ENTREZ","SYMBOL","GENENAME")
  
  tbl2=Ca_ChromosomalFeatures[, c("locus_tag.CBS138", "Chromosome.contig.name.DSY562")]
  colnames(tbl2)=c("GID","CHROMOSOME")
  # colnames(tbl2)=c("ENTREZ","CHROMOSOME")
  
  tbl3=GeneAssociation[!duplicated(GeneAssociation[, c("CBS138_ID", "GOID", "evidence")]), c("CBS138_ID", "GOID", "evidence")]
  colnames(tbl3)=c("GID","GO","EVIDENCE")
  # colnames(tbl3)=c("ENTREZ","GO","EVIDENCE")
  
  Output=AnnotationForge::makeOrgPackage(gene_info=tbl1, chromosome=tbl2, go=tbl3,
                                         version="0.1",
                                         maintainer="Eric Durandau <eric.durandau@chuv.ch>",
                                         author="Eric Durandau <eric.durandau@chuv.ch>",
                                         outputDir = "/media/edurandau/DATA/Eric_Durandau/Sanglard/Data/Sequencing/Databases/AnnotationForge",
                                         tax_id=TaxonID,
                                         genus="Candida",
                                         species="glabrata",
                                         goTable="go")
  return(Output)
}


GOE=function(GeneList, Universe, OrgDb, pAdjustMethod = "BH", pvalueCutoffE  = 0.05,qvalueCutoffE  = 0.2, keyType = "GID", ONT=c("CC", "BP", "MF"), MinGOCount=3, MaxGOCount=50){
  # GOData=list(
  #   "BP"= GOSemSim::godata(OrgDb ,ont="BP", keytype = keyType, computeIC = F),
  #   "CC"= GOSemSim::godata(OrgDb ,ont="CC", keytype = keyType, computeIC = F),
  #   "MF"=GOSemSim::godata(OrgDb ,ont="MF", keytype = keyType, computeIC = F)
  # )
  
  GOSum=data.frame()
  PlotList=list()
  ind=1
  GenesID=names(GeneList)
  if (length(GeneList)>5){
      
      # Over-representation test(Boyle et al. 2004) were implemented in clusterProfiler. For calculation details and explanation of paramters, please refer to the vignette of DOSE.
      
    for (i in ONT){
      ego <- clusterProfiler::enrichGO(gene= GenesID,
                                       OrgDb         = OrgDb,
                                       ont           = i,
                                       keyType = keyType,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = pvalueCutoffE,
                                       qvalueCutoff  = qvalueCutoffE,
                                       readable = T,
                                       universe = Universe,
                                       minGSSize = MinGOCount,
                                       maxGSSize = MaxGOCount
      )
        
        ## enrichGO test the whole GO corpus and enriched result may contains very general terms. With dropGO function, user can remove specific GO terms or GO level from results obtained from both enrichGO and compareCluster.
        # gego=dropGO(ego)
        
        ## reduce redundancy of enriched GO terms
        ## GO is organized in parent-child structure, thus a parent term can be overlap with a large proportion with all its child terms. This can result in redundant findings. To solve this issue, clusterProfiler implement simplify method to reduce redundant GO terms from the outputs of enrichGO and gseGO. The function internally called GOSemSim (Yu et al. 2010) to calculate semantic similarity among GO terms and remove those highly similar terms by keeping one representative term. An example can be found in the blog post.
        
      d1 <- GOSemSim::godata(OrgDb ,ont=i, keytype = keyType, computeIC = F)
      
        # d2 <- godata('org.Calbicans.eg.db', ont="MF", keytype = "GID")
        # d3 <- godata('org.Calbicans.eg.db', ont="CC", keytype = "GID")
        # ego2 <- pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 40)
        # ego2 <- pairwise_termsim(ego2, method = "Wang", semData = d2)
        # ego2 <- pairwise_termsim(ego2,method = "Wang", semData = d3)
        
      
      # sego=ego
      # if (nrow(ego@result)<3){
      #   sego=ego
      # }
      x=tryCatch(pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 25), error = function(e){e})
      IsError=class(x)[[1]]!="enrichResult"
      
      if (IsError | length(ego$ID)<3){
        
      }else {
        sego <- pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 25)
        # sego@result$Description=str_wrap(sego@result$Description,width = 40)
        
        
        
        P4=emapplot(sego,showCategory = 10,
                    layout="kk", 
                    label_format = 5,
                    cex_label_category = 0.5,
                    cex_category = 0.5,
                    cex_line = 0.5,
                    group_category = T,
                    node_label="group")
        
        sego@result$Description=str_wrap(sego@result$Description,width = 30)
        sego2=simplify(sego)
        
        Title=i
        # PlotList[[paste0("Clust_", as.character(ClusterNb))]][[paste0("Plot_P4_",i)]]= P4
      
      
      # barplot(ego2, showCategory=40,
      # color = "p.adjust" )
        P1=enrichplot::dotplot(sego2,showCategory=15, title= i, font.size=7)
        # P1=barplot(sego2, showCategory=30,
        #            color = "p.adjust", font.size = 7, title = Title)
        
        # dotplot(ego2,color = "p.adjust", showCategory=40,title = paste("Cluster", ClusterNb))
        P2=goplot(sego,color = "p.adjust", showCategory=5 ,title = paste("Cluster", ClusterNb), title=i)
        # cnetplot(ego2, categorySize="p.adjust", showCategory = 12)
        P3=cnetplot(
          sego2,
          showCategory = 10,
          foldChange = GeneList,
          layout = "kk",
          colorEdge = T,
          circular = FALSE,
          node_label = "gene",
          cex_category = 0.5,
          cex_gene =0.5,
          node_label_size = 0.5,
          cex_label_category = 0.5,
          categorySize="pvalue",
          cex_label_gene = 0.5, title=i ) + scale_color_gradient(low = "blue", high = "red")+labs(color="Highest\nlog2(FC)")+ ggtitle(i)
          
          
          # emapplot_cluster(ego2,showCategory = 40
          # , layout="kk")
        Sum=sego2@result
        Sum$GO_Class=i
        GOSum=rbind(GOSum,Sum)
        
        PlotList[[i]]=list(
          Barplot=P1,
          GoPlot=P2,
          CNetPlot=P3,
          EnriMapPlot=P4,
          ego=sego2
          
          
        )
        
      }
    }
    PlotList[["Summary"]]=GOSum
  }
  
  return(PlotList)
}

GO2GeneTableFromPlotList=function(List,ONT=c("MF","CC", "BP"), Path=F){
  
  
}


GSEA_hm=function(GeneList, OrgDb, pAdjustMethod = "BH",pValCutOff=0.05,  keyType = "GID", ONT=c("CC", "BP", "MF"), MinGOCount=2, MaxGOCount=500){
  # GOData=list(
  #   "BP"= GOSemSim::godata(OrgDb ,ont="BP", keytype = keyType, computeIC = F),
  #   "CC"= GOSemSim::godata(OrgDb ,ont="CC", keytype = keyType, computeIC = F),
  #   "MF"=GOSemSim::godata(OrgDb ,ont="MF", keytype = keyType, computeIC = F)
  # )
  
  GOSum=data.frame()
  PlotList=list()
  ind=1
  GenesID=names(GeneList)
  if (length(GeneList)>5){
    
    # Over-representation test(Boyle et al. 2004) were implemented in clusterProfiler. For calculation details and explanation of paramters, please refer to the vignette of DOSE.
    
    for (i in ONT){
      ego <- clusterProfiler::gseGO(gene= GeneList,
                                       OrgDb         = OrgDb,
                                       ont           = i,
                                       keyType = keyType,
                                       pAdjustMethod = "BH",
                                       minGSSize = MinGOCount,
                                       maxGSSize = MaxGOCount,
                                    pvalueCutoff = pValCutOff
      )
      
      ## enrichGO test the whole GO corpus and enriched result may contains very general terms. With dropGO function, user can remove specific GO terms or GO level from results obtained from both enrichGO and compareCluster.
      # gego=dropGO(ego)
      
      ## reduce redundancy of enriched GO terms
      ## GO is organized in parent-child structure, thus a parent term can be overlap with a large proportion with all its child terms. This can result in redundant findings. To solve this issue, clusterProfiler implement simplify method to reduce redundant GO terms from the outputs of enrichGO and gseGO. The function internally called GOSemSim (Yu et al. 2010) to calculate semantic similarity among GO terms and remove those highly similar terms by keeping one representative term. An example can be found in the blog post.
      
      d1 <- GOSemSim::godata(OrgDb ,ont=i, keytype = keyType, computeIC = F)
      
      # d2 <- godata('org.Calbicans.eg.db', ont="MF", keytype = "GID")
      # d3 <- godata('org.Calbicans.eg.db', ont="CC", keytype = "GID")
      # ego2 <- pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 40)
      # ego2 <- pairwise_termsim(ego2, method = "Wang", semData = d2)
      # ego2 <- pairwise_termsim(ego2,method = "Wang", semData = d3)
      
      
      # sego=ego
      # if (nrow(ego@result)<3){
      #   sego=ego
      # }
      x=tryCatch(pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 25), error = function(e){e})
      IsError=class(x)[[1]]!="gseaResult"
      
      if (IsError | length(ego$ID)<3){
        
      }else {
        sego <- pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 25)
        sego@result$Description=str_wrap(sego@result$Description,width = 40)
        sego2=simplify(sego)
        sego2=setReadable(sego2,OrgDb )
        
        
        P4=emapplot_cluster(sego,showCategory = 25
                            , layout="kk", label_format = 5, title=i)
        Title=i
        # PlotList[[paste0("Clust_", as.character(ClusterNb))]][[paste0("Plot_P4_",i)]]= P4
        
        
        # barplot(ego2, showCategory=40,
        # color = "p.adjust" )
        P1=enrichplot::dotplot(sego2,showCategory=20, title= i)
        # P1=barplot(sego2, showCategory=30,
        #            color = "p.adjust", font.size = 7, title = Title)
        
        # dotplot(ego2,color = "p.adjust", showCategory=40,title = paste("Cluster", ClusterNb))
        # P2=goplot(sego,color = "p.adjust", showCategory=5 ,title = paste("Cluster", ClusterNb), title=i)
        
        # cnetplot(ego2, categorySize="p.adjust", showCategory = 12)
        P3=cnetplot(
          sego2,
          showCategory = 8,
          foldChange = GeneList,
          layout = "kk",
          colorEdge = T,
          circular = FALSE,
          node_label = "gene",
          cex_category = 0.8,
          cex_gene = 1,
          node_label_size = NULL,
          cex_label_category = 1,
          cex_label_gene = 0.5, title=i ) + scale_color_gradient(low = "blue", high = "red")+labs(color="log2(FC)")
        
        
        # emapplot_cluster(ego2,showCategory = 40
        # , layout="kk")
        Sum=sego2@result
        Sum$GO_Class=i
        GOSum=rbind(GOSum,Sum)
        
        PlotList[[i]]=list(
          Barplot=P1,
          CNetPlot=P3,
          EnriMapPlot=P4,
          ego=sego2
          
          
        )
        
      }
    }
    PlotList[["Summary"]]=GOSum
  }
  
  return(PlotList)
}

GSEA_KEGG=function(GeneList, Organism='hsa',pValCutOff=0.05,  keyType = "ncbi-geneid", MinGOCount=2, MaxGOCount=500){
  # GOData=list(
  #   "BP"= GOSemSim::godata(OrgDb ,ont="BP", keytype = keyType, computeIC = F),
  #   "CC"= GOSemSim::godata(OrgDb ,ont="CC", keytype = keyType, computeIC = F),
  #   "MF"=GOSemSim::godata(OrgDb ,ont="MF", keytype = keyType, computeIC = F)
  # )
  GeneList=sort(GeneList, decreasing = T)
  GOSum=data.frame()
  PlotList=list()
  ind=1
  GenesID=names(GeneList)
  if (length(GeneList)>5){
    
    # Over-representation test(Boyle et al. 2004) were implemented in clusterProfiler. For calculation details and explanation of paramters, please refer to the vignette of DOSE.
    

    ego <- clusterProfiler::gseKEGG(gene= GeneList,
                                  organism = Organism,
                                  keyType = keyType,
                                  pAdjustMethod = "BH",
                                  minGSSize = MinGOCount,
                                  maxGSSize = MaxGOCount,
                                  pvalueCutoff = pValCutOff,
    )
    
    ## enrichGO test the whole GO corpus and enriched result may contains very general terms. With dropGO function, user can remove specific GO terms or GO level from results obtained from both enrichGO and compareCluster.
    # gego=dropGO(ego)
    
    ## reduce redundancy of enriched GO terms
    ## GO is organized in parent-child structure, thus a parent term can be overlap with a large proportion with all its child terms. This can result in redundant findings. To solve this issue, clusterProfiler implement simplify method to reduce redundant GO terms from the outputs of enrichGO and gseGO. The function internally called GOSemSim (Yu et al. 2010) to calculate semantic similarity among GO terms and remove those highly similar terms by keeping one representative term. An example can be found in the blog post.
    
    # d1 <- GOSemSim::godata(OrgDb ,ont=i, keytype = keyType, computeIC = F)
    
    # d2 <- godata('org.Calbicans.eg.db', ont="MF", keytype = "GID")
    # d3 <- godata('org.Calbicans.eg.db', ont="CC", keytype = "GID")
    # ego2 <- pairwise_termsim(ego, method = "Wang", semData = d1, showCategory = 40)
    # ego2 <- pairwise_termsim(ego2, method = "Wang", semData = d2)
    # ego2 <- pairwise_termsim(ego2,method = "Wang", semData = d3)
    
    
    # sego=ego
    # if (nrow(ego@result)<3){
    #   sego=ego
    # }
    
    if ( length(ego$ID)<3){
      
    }else {
     
      ego@result$Description=str_wrap(ego@result$Description,width = 40)
     
      sego2=ego      
      
      # Title=i
      # PlotList[[paste0("Clust_", as.character(ClusterNb))]][[paste0("Plot_P4_",i)]]= P4
      
      
      # barplot(ego2, showCategory=40,
      # color = "p.adjust" )
      P1=enrichplot::dotplot(sego2,showCategory=20, title= "Over-represented KEGG pathway")
      # P1=barplot(sego2, showCategory=30,
      #            color = "p.adjust", font.size = 7, title = Title)
      
      # dotplot(ego2,color = "p.adjust", showCategory=40,title = paste("Cluster", ClusterNb))
      # P2=goplot(sego,color = "p.adjust", showCategory=5 ,title = paste("Cluster", ClusterNb), title="Major En")
      
      # cnetplot(ego2, categorySize="p.adjust", showCategory = 12)
      P3=cnetplot(
        sego2,
        showCategory = 20,
        foldChange = GeneList,
        layout = "kk",
        colorEdge = T,
        circular = FALSE,
        node_label = "gene",
        cex_category = 0.8,
        cex_gene = 1,
        node_label_size = NULL,
        cex_label_category = 1,
        cex_label_gene = 0.5, title= "Major KEGG gene association" ) + scale_color_gradient(low = "blue", high = "red")+labs(color="log2(FC)")
      
      
      # emapplot_cluster(ego2,showCategory = 40
      # , layout="kk")
      Sum=sego2@result
      GOSum=rbind(GOSum,Sum)
      
      PlotList=list(
        Barplot=P1,
        CNetPlot=P3,
        ego=sego2
        )
    }
    PlotList[["Summary"]]=GOSum
  }
  
  return(PlotList)
}