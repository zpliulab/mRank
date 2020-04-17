NOA <- function(upfile,species="Yeast",idinput,reffile=NULL,evidences="ALL"){
 # This function is a R implement for the noa annotation analysis, which contains more species and ids, and also the result show.
 
 # upfile: upload file to analysis which contains a gene list or gene network
 # reffile: reference gene list or reference gene network file
 ####### threshold: cut off values for P-value
 # species: the sepcies
 # idinput: the input id, which is necessary for the id-mapping process
 # evidences: the GO ontology evidences, default, all used
 
 # step 1: whether to do id mapping
 speinfor <- species.list(species) # get the species information
 
 # step 2: gene list analysis or gene network analysis
 anaset <- read.file(upfile,reffile) # read the upload file to the analysis set

  glmresult <- NULL
  noaresult <- NULL
 if(anaset$method==1){
  glmresult <- GLM(anaset,speinfor,idinput,evidences)
 }
 if(anaset$method==2){
  noaresult <- noa.test(anaset,speinfor,idinput,upfile,reffile,evidences)
  anasetglm <- noa2glm(anaset,reffile)
  glmresult <- GLM(anasetglm,speinfor,idinput,evidences)
 }
    
 list(noaresult=noaresult,glmresult=glmresult)
 
}

noa2glm <- function(anaset,reffile){
 anasetglm <- list()
 anasetglm$genepair <- rbind(anaset$genepair[,1],anaset$genepair[,2])
 anasetglm$method <- 1
 anasetglm$geneset <- anaset$geneset
 if(!is.null(reffile)){
   anasetglm$refgenepair <- rbind(anaset$refgenepair[,1],anaset$refgenepair[,2])
   anasetglm$refgeneset <- anaset$refgeneset
 }
 
 anasetglm
}

noa.test <- function(anaset,speinfor,idinput,upfile,reffile,evidences){
 
 # id mapping and call term table
 if(!(idinput %in% speinfor$idmaps)){stop("The input id for the corresponding species are not avaiable!");}
   id1 <- idinput
   id2 <- speinfor$annotationid
   species <- speinfor$name
   # id mapping for geneset
   anaset$mapgeneset <- id.mapping(anaset$geneset,id1,id2,species)
   unmapset1 <- anaset$mapgeneset[anaset$mapgeneset[,2]=="gap",1]
   # call term table 
   termtable <- gene.term(anaset$mapgeneset[,2],species,evidences)
   termtable <- map.assign(termtable,anaset$mapgeneset,method=1)
   testnet <- read.net(upfile)
   testnet$matrix <- testnet$matrix + t(testnet$matrix)
   testnet$matrix[testnet$matrix>1] <- 1
   diag(testnet$matrix) <- 0
   
   # id mapping and call term table for reference genes
   if(is.null(anaset$refgeneset)){
	anaset$maprefgeneset <- anaset$mapgeneset
	anaset$refgeneset <- anaset$geneset
	reftermtable <- termtable
	refnet <- testnet
	refnet$matrix <- matrix(1,refnet$size,refnet$size)
	diag(refnet$matrix) <- 0  
	dimnames(refnet$matrix) <- dimnames(testnet$matrix)

	anaset$unmappedgene <- unmapset1
   }else{
    anaset$maprefgeneset <- id.mapping(anaset$refgeneset,id1,id2,species)
	reftermtable <- gene.term(anaset$maprefgeneset[,2],species,evidences)
	reftermtable <- map.assign(reftermtable,anaset$maprefgeneset,method=1)
	refnet <- read.net(reffile)
	refnet$matrix <- refnet$matrix + t(refnet$matrix)
    refnet$matrix[refnet$matrix>1] <- 1
    diag(refnet$matrix) <- 0
	
	unmapset2 <- anaset$maprefgeneset[anaset$maprefgeneset[,2]=="gap",1]
    anaset$unmappedgene <- union(unmapset1,unmapset2)
   }
     
 # hypergeometric test  
 #termn <- length(unique(termtable[,2])) # all test term number
 term <- matrix(0,0,8,dimnames=list(c(),c("name","Ontology","pvalue","O","T","G","R","annotated edges")))
  
 ontology <- c("BP","CC","MF")
 for(i in ontology){
  termattr <- edge.test(termtable,reftermtable,testnet,refnet,i)
  term <- rbind(term,termattr)
 }
 
 anaset$term <- term
 
 anaset
}

edge.test <- function(termtable,reftermtable,testnet,refnet,Ontology){
 termtable <- termtable[termtable[,4]==Ontology,]
 reftermtable <- reftermtable[reftermtable[,4]==Ontology,]
 
 tname <- as.vector(unique(termtable[,2])) # term name
 termn <- length(tname)
 kterm <- 0 # iter for real useful terms
 
 tmatrix <- matrix(0,termn,8,dimnames=list(c(),c("name","Ontology","pvalue","x","M","X","N","annotated edges")))
 N <- sum(refnet$matrix)/2
 M <- sum(testnet$matrix)/2
 for(i in 1:termn){
  genetest <- unique(termtable[termtable[,2]==tname[i],1])
  x <- sum(testnet$matrix[genetest,genetest])/2
  generef <- unique(reftermtable[reftermtable[,2]==tname[i],1])
  X <- sum(refnet$matrix[generef,generef])/2
  
  if(x!=0){
   kterm <- kterm + 1
   genematrix <- testnet$matrix[genetest,genetest]
   p <- phyper(x-1, X, N-X, M, lower.tail = FALSE, log.p = FALSE)
   annogene <- ""
   ind <- which(genematrix==1,arr.ind=TRUE) 
   for(j in 1:(2*x)){
    if(ind[j,1]>ind[j,2]){
    geneedge <- paste(genetest[ind[j,1]],genetest[ind[j,2]],sep=" ")
    annogene <- paste(annogene,geneedge,";",sep="")
    }}
   tmatrix[kterm,] <- c(tname[i],Ontology,p,x,M,X,N,annogene)
  }
  
  }
  
  if(kterm>0){
	tmatrix <- tmatrix[1:kterm,]
	#print(kterm)
	pvalues <- as.numeric(tmatrix[,3])
	sorder <- sort(pvalues,index.return=TRUE)
	pvalues <- sorder$x
	# corrected for P values
	pvalues <- pvalues * kterm
	pvalues[pvalues>1] <- 1
	tmatrix <- tmatrix[sorder$ix,]
	tmatrix[,3] <- pvalues
  }else{tmatrix <- NULL;}
  
  tmatrix
}

GLM <- function(anaset,speinfor,idinput,evidences){

 # id mapping
 if(!(idinput %in% speinfor$idmaps)){stop("The input id for the corresponding species are not avaiable!");}
   id1 <- idinput
   id2 <- speinfor$annotationid
   species <- speinfor$name
   # id mapping for geneset
   anaset$mapgeneset <- id.mapping(anaset$geneset,id1,id2,species)
   # id mapping for reference genes
   flag <- 0
   if(is.null(anaset$refgeneset)){
    refmap <- toupper(as.matrix(idmap.table(id1,id2,species)))
	anaset$maprefgeneset <- refmap[,c(id1,id2)]
	anaset$refgeneset <- unique(refmap[,id1])
	flag <- 1
   }else{
    anaset$maprefgeneset <- id.mapping(anaset$refgeneset,id1,id2,species)
   }
 # call term table
 termtable <- gene.term(anaset$mapgeneset[,2],species,evidences)
 reftermtable <- gene.term(anaset$maprefgeneset[,2],species,evidences)
   # map assign for term table 
   termtable <- map.assign(termtable,anaset$mapgeneset,method=1)
   reftermtable <- map.assign(reftermtable,anaset$maprefgeneset,method=1)
 # hypergeometric test 
 unmapset1 <- anaset$mapgeneset[anaset$mapgeneset[,2]=="gap",1]
 if(flag==0){
  unmapset2 <- anaset$maprefgeneset[anaset$maprefgeneset[,2]=="gap",1]
  anaset$unmappedgene <- union(unmapset1,unmapset2)
  }else{anaset$unmappedgene <- unmapset1}
  
 #termn <- length(unique(termtable[,2])) # all test term number
 term <- matrix(0,0,8,dimnames=list(c(),c("name","Ontology","pvalue","O","T","G","R","annotated gene")))
  
 ontology <- c("BP","CC","MF")
 for(i in ontology){
  termattr <- term.attributes(termtable,reftermtable,i)
  term <- rbind(term,termattr)
 }
 
 anaset$term <- term
 
 anaset

}

term.attributes <- function(termtable,reftermtable,Ontology){  
 termtable <- termtable[termtable[,4]==Ontology,]
 reftermtable <- reftermtable[reftermtable[,4]==Ontology,]
 
 tname <- as.vector(unique(termtable[,2])) # term name
 termn <- length(tname)
 tmatrix <- matrix(0,termn,8,dimnames=list(c(),c("name","Ontology","pvalue","x","M","X","N","annotated gene")))
 N <- length(unique(reftermtable[,1]))
 M <- length(unique(termtable[,1]))
 #print(termn)
 for(i in 1:termn){
  X <- length(unique(reftermtable[reftermtable[,2]==tname[i],1]))
  annogeneset <- unique(termtable[termtable[,2]==tname[i],1])
  x <- length(annogeneset)
  p <- phyper(x-1, X, N-X, M, lower.tail = FALSE, log.p = FALSE)
  annogene <- c()
  for(j in 1:x){
   annogene <- paste(annogene,annogeneset[j],sep=" ")
  }
  tmatrix[i,] <- c(tname[i],Ontology,p,x,M,X,N,annogene)
  }
  
  if(termn>0){
	pvalues <- as.numeric(tmatrix[,3])
	sorder <- sort(pvalues,index.return=TRUE)
	pvalues <- sorder$x
	# corrected for P values
	pvalues <- pvalues * termn
	pvalues[pvalues>1] <- 1
	tmatrix <- tmatrix[sorder$ix,]
	tmatrix[,3] <- pvalues
  }else{tmatrix <- NULL;}
  
  tmatrix
}

gene.term <- function(geneset,species,evidences="ALL"){
 ## evidences== ALL
 ## c("EXP","IDA","IPI","IMP","IGI","IEP","ISS","ISO","ISA","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","IC","ND","","IEA","NR")
 termtable <- c()
 if(species=="Yeast"){library(org.Sc.sgd.db);x <- org.Sc.sgdGO2ALLORFS;}
 if(species=="Human"){library(org.Hs.eg.db);x <- org.Hs.egGO2ALLEGS;}
 if(species=="Fly"){library(org.Dm.eg.db);x <- org.Dm.egGO2ALLEGS;}
 if(species=="Worm"){library(org.Ce.eg.db);x <- org.Ce.egGO2ALLEGS;}
 if(species=="Mouse"){library(org.Mm.eg.db);x <- org.Mm.egGO2ALLEGS;}
 if(species=="Arabidopsis"){library(org.At.tair.db);x <- org.At.tairGO2ALLTAIRS;}
 if(species=="Rat"){library(org.Rn.eg.db);x <- org.Rn.egGO2ALLEGS;}
 if(species=="Zebrafish"){library(org.Dr.eg.db);x <- org.Dr.egGO2ALLEGS;}
 if(species=="Bovine"){library(org.Bt.eg.db);x <- org.Bt.egGO2ALLEGS;}
 if(species=="Canine"){library(org.Cf.eg.db);x <- org.Cf.egGO2ALLEGS;}
 if(species=="Anopheles"){library(org.Ag.eg.db);x <- org.Ag.egGO2ALLEGS;}
 if(species=="E coli strain Sakai"){library(org.EcSakai.eg.db);x <- org.EcSakai.egGO2ALLEGS;}
 if(species=="Chicken"){library(org.Gg.eg.db);x <- org.Gg.egGO2ALLEGS;}
 if(species=="Chimp"){library(org.Pt.eg.db);x <- org.Pt.egGO2ALLEGS;}
 if(species=="Malaria"){library(org.Pf.plasmo.db);x <- org.Pf.plasmoGO2ALLORFS;}
 if(species=="Rhesus"){library(org.Mmu.eg.db);x <- org.Mmu.egGO2ALLEGS;}
 if(species=="Streptomyces coelicolor"){library(org.Sco.eg.db);x <- org.Sco.egGO2ALLEGS;}
 if(species=="Pig"){library(org.Ss.eg.db);x <- org.Ss.egGO2ALLEGS;}
 if(species=="Xenopus"){library(org.Xl.eg.db);x <- org.Xl.egGO2ALLEGS;}
 if(species=="E coli strain K12"){library(org.EcK12.eg.db);x <- org.EcK12.egGO2ALLEGS;}
 
  xx <- as.matrix(toTable(x))
  xx[,1] <- toupper(xx[,1])
  if(!is.null(geneset)){
   sub0 <- xx[,1] %in% geneset
   termtable <- xx[sub0,]
  }else{
   termtable <- xx
  }
 
 if(evidences!="ALL"){
  sub1 <- termtable[,3] %in% evidences
  termtable <- termtable[sub1,]
 }
 
 termtable
}

map.assign <- function(tertable,map,method=1){
 
 tertable <- as.matrix(tertable)
 map <- as.matrix(map)
 if(method==1){
  sub1 <- match(tertable[,1], map[,2])
  tertable[,1] <- map[sub1,1]
 }

 tertable
}

read.file <- function(upfile,reffile){

 # gene list: we delete the repeat genes as one
 # gene edge: we do not delete the repeat edges: note!!!!!
 anaset <- list()
 
 anaset$genepair <- toupper(as.matrix(read.table(upfile)))
 anaset$method <- dim(anaset$genepair)[2]
 
 if(anaset$method==1){
  anaset$geneset <- unique(anaset$genepair)
 }else{
  anaset$geneset <- unique(union(anaset$genepair[,1],anaset$genepair[,2]))
 }
 
 if(is.null(reffile)) { anaset$refgeneset <- NULL; anaset$refgenepair <- NULL }
 else{
  anaset$refgenepair <- toupper(as.matrix(read.table(reffile)))
  if(anaset$method==1){
   anaset$refgeneset <- unique(anaset$refgenepair)
  }else{
   anaset$refgeneset <- unique(union(anaset$refgenepair[,1],anaset$refgenepair[,2]))
  }
 }
 
 anaset
 
}

id.mapping <- function(genelist,id1,id2,species){
 # genelist: gene list to do id mapping
 # id mapping for a gene list from id1 to id2
 # species: the species
 
 # we simplely deal the one to many or many to one to one-to-one
 if(is.null(genelist)){ # the robust represent
  mapgene <- NULL
 }else if(id1==id2){
  mapgene <- cbind(genelist,genelist)
  colnames(mapgene) <- c(id1,id2)
 }else{
  mapgene <- matrix("gap",length(genelist),2)
  X <- toupper(as.matrix(idmap.table(id1,id2,species)))
  sub <- genelist %in% X[,id1]  ## many to one or one to many ?????
  sub0 <- !sub
  mappedtable <- X[X[,id1] %in% genelist,]
  mappedtable <- unique.table(mappedtable)
  if(any(sub0)) mapgene <- rbind(cbind(mappedtable[,id1],mappedtable[,id2]),cbind(genelist[sub0],"gap"))
  else mapgene <- cbind(mappedtable[,id1],mappedtable[,id2])
  colnames(mapgene) <- c(id1,id2)
 }
 
 mapgene
}

unique.table <- function(mappedtable){
 # make the unique mapping for the id mapping results
 
	# the duplicated elements sub
	subdup0 <- !duplicated(mappedtable[,1])
	subdup1 <- !duplicated(mappedtable[,2])
	# the unique matching sub
	subdup <- subdup0 & subdup1
	# the unique mathcing pairs
	unipair <- mappedtable[subdup,]

	unipair

}

idmap.table <- function(id1,id2,species){
 # 
 X <- c()
 if(species=="Yeast"){
  library(org.Sc.sgd.db)
  if(id1=="ALIAS" & id2=="ORF" ) {X <- toTable(org.Sc.sgdALIAS);colnames(X) <- c(id2,id1);}
  if(id1=="COMMON" & id2=="ORF" ) {X <- toTable(org.Sc.sgdCOMMON2ORF);colnames(X) <- c(id1,id2);}
  if(id1=="GENENAME" & id2=="ORF" ){X <- toTable(org.Sc.sgdGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SGD" & id2=="ORF" ){X <- toTable(org.Sc.sgdSGD);colnames(X) <- c(id2,id1);}
  if(id1=="SMART" & id2=="ORF" ){X <- toTable(org.Sc.sgdSMART);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="ORF" ){X <- toTable(org.Sc.sgdUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="ORF" ){X <- toTable(org.Sc.sgdREFSEQ);colnames(X) <- c(id2,id1);}
  if(id1=="PUBMED" & id2=="ORF" ){X <- toTable(org.Sc.sgdPMID);colnames(X) <- c(id2,id1);}
  if(id1=="Entrez" & id2=="ORF" ){X <- toTable(org.Sc.sgdENTREZID);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="ORF" ){X <- toTable(org.Sc.sgdENSEMBL);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Malaria"){
  library(org.Pf.plasmo.db)
  if(id1=="ALIAS" & id2=="ORF" ){X <- toTable(org.Pf.plasmoALIAS2ORF);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="ORF" ){X <- toTable(org.Pf.plasmoGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="ORF" ){X <- toTable(org.Pf.plasmoSYMBOL);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Arabidopsis"){
  library(org.At.tair.db)
  if(id1=="Entrez" & id2=="TAIR" ){X <- toTable(org.At.tairENTREZID);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="TAIR" ){X <- toTable(org.At.tairGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="TAIR" ){X <- toTable(org.At.tairSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="TAIR" ){X <- toTable(org.At.tairREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Anopheles"){
  library(org.Ag.eg.db)
  #c("ACCNUM","ENSEMBL","ENSEMBLPROT","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez");
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Ag.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Ag.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Ag.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Ag.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Ag.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Ag.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Ag.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="E coli strain Sakai"){
  #c("ACCNUM","ALIAS","GENENAME","SYMBOL","REFSEQ","Entrez");
  library(org.EcSakai.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.EcSakai.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.EcSakai.egALIAS2EG);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.EcSakai.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.EcSakai.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.EcSakai.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Chimp"){
  #c("ACCNUM","ENSEMBL","ENSEMBLPROT","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez")
  library(org.Pt.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Pt.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Pt.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Pt.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Pt.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Pt.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Pt.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Pt.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Rhesus"){
  #c("ACCNUM","ENSEMBL","ENSEMBLPROT","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez");
  library(org.Mmu.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Mmu.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Mmu.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Mmu.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Mmu.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Mmu.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Mmu.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Mmu.egREFSEQ);colnames(X) <- c(id2,id1);}
 }

 if(species=="Streptomyces coelicolor"){
  #c("ACCNUM","ALIAS","GENENAME","SYMBOL","REFSEQ","Entrez");
  library(org.Sco.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Sco.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Sco.egALIAS2EG);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Sco.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Sco.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Sco.egREFSEQ);colnames(X) <- c(id2,id1);}
  } 

 if(species=="Pig"){
  #c("ACCNUM","ALIAS","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez");
  library(org.Ss.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Ss.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Ss.egALIAS2EG);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Ss.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Ss.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Ss.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Ss.egREFSEQ);colnames(X) <- c(id2,id1);}
  }
 
 if(species=="Xenopus"){
  #c("ACCNUM","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez");
  library(org.Xl.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Xl.egACCNUM);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Xl.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Xl.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Xl.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Xl.egREFSEQ);colnames(X) <- c(id2,id1);}
  }  
 
 if(species=="E coli strain K12"){
  #c("ACCNUM","ALIAS","GENENAME","SYMBOL","REFSEQ","Entrez");
  library(org.EcK12.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.EcK12.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.EcK12.egALIAS2EG);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.EcK12.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.EcK12.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.EcK12.egREFSEQ);colnames(X) <- c(id2,id1);}
  }  

 #"Human","Fly","Worm","Mouse","Rat","Zebrafish","Bovine","Canine","Chicken"
 if(species=="Human"){
 #c("ACCNUM","ALIAS","ENSEMBL","ENSEMBLPROT","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez");
 library(org.Hs.eg.db)
 if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Hs.egACCNUM);colnames(X) <- c(id2,id1);}
 if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Hs.egALIAS2EG);colnames(X) <- c(id2,id1);}
 if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Hs.egENSEMBL);colnames(X) <- c(id2,id1);}
 if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Hs.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
 #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Hs.egGENENAME);colnames(X) <- c(id2,id1);}
 if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Hs.egSYMBOL);colnames(X) <- c(id2,id1);}
 if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Hs.egUNIPROT);colnames(X) <- c(id2,id1);}
 if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Hs.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Fly"){
  library(org.Dm.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Dm.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Dm.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Dm.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Dm.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Dm.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Dm.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Dm.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Dm.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Worm"){
  library(org.Ce.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Ce.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Ce.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Ce.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Ce.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Ce.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Ce.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Ce.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Ce.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Mouse"){
  library(org.Mm.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Mm.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Mm.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Mm.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Mm.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Mm.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Mm.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Mm.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Mm.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Rat"){
  library(org.Rn.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Rn.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Rn.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Rn.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Rn.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Rn.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Rn.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Rn.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Rn.egREFSEQ);colnames(X) <- c(id2,id1);}
 }

 if(species=="Zebrafish"){
  library(org.Dr.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Dr.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Dr.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Dr.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Dr.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Dr.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Dr.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Dr.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Dr.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Bovine"){
  library(org.Bt.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Bt.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Bt.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Bt.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Bt.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Bt.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Bt.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Bt.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Bt.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Canine"){
  library(org.Cf.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Cf.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Cf.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Cf.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Cf.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Cf.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Cf.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Cf.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Cf.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 if(species=="Chicken"){
  library(org.Gg.eg.db)
  if(id1=="ACCNUM" & id2=="Entrez" ){X <- toTable(org.Gg.egACCNUM);colnames(X) <- c(id2,id1);}
  if(id1=="ALIAS" & id2=="Entrez" ){X <- toTable(org.Gg.egALIAS2EG);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBL" & id2=="Entrez" ){X <- toTable(org.Gg.egENSEMBL);colnames(X) <- c(id2,id1);}
  if(id1=="ENSEMBLPROT" & id2=="Entrez" ){X <- toTable(org.Gg.egENSEMBLPROT);colnames(X) <- c(id2,id1);}
  #if(id1=="GENENAME" & id2=="Entrez" ){X <- toTable(org.Gg.egGENENAME);colnames(X) <- c(id2,id1);}
  if(id1=="SYMBOL" & id2=="Entrez" ){X <- toTable(org.Gg.egSYMBOL);colnames(X) <- c(id2,id1);}
  if(id1=="UNIPROT" & id2=="Entrez" ){X <- toTable(org.Gg.egUNIPROT);colnames(X) <- c(id2,id1);}
  if(id1=="REFSEQ" & id2=="Entrez" ){X <- toTable(org.Gg.egREFSEQ);colnames(X) <- c(id2,id1);}
 }
 
 X
}

species.list <- function(species){
 # species: the species
 # this function gives the specific information about a species
 
 speinfor <- list()
 speinfor$name <- species
 speinfor$annotationid <- "Entrez";
 
 idmaps0 <- c("ACCNUM","ALIAS","ENSEMBL","ENSEMBLPROT","SYMBOL","UNIPROT","REFSEQ","Entrez");
 specie0 <- c("Human","Fly","Worm","Mouse","Rat","Zebrafish","Bovine","Canine","Chicken");
 if(species%in%specie0) speinfor$idmaps <-  idmaps0
 else{
 
 if(species=="Yeast"){speinfor$annotationid <- "ORF";
  speinfor$idmaps <-  c("ALIAS","ENSEMBL","COMMON","GENENAME","SGD","SMART","UNIPROT","REFSEQ","PUBMED","Entrez","ORF");}

 if(species=="Malaria"){speinfor$annotationid <- "ORF";
  speinfor$idmaps <-  c("ALIAS","SYMBOL","ORF");}  
  
 if(species=="Arabidopsis"){speinfor$annotationid <- "TAIR";
  speinfor$idmaps <-  c("Entrez","SYMBOL","REFSEQ","TAIR");} 

 if(species=="Anopheles"){
  speinfor$idmaps <-  c("ACCNUM","ENSEMBL","ENSEMBLPROT","GENENAME","SYMBOL","UNIPROT","REFSEQ","Entrez");} 

 if(species=="E coli strain Sakai"){
  speinfor$idmaps <-  c("ACCNUM","ALIAS","SYMBOL","REFSEQ","Entrez");} 

 if(species=="Chimp"){
  speinfor$idmaps <-  c("ACCNUM","ENSEMBL","ENSEMBLPROT","SYMBOL","UNIPROT","REFSEQ","Entrez");} 

 if(species=="Rhesus"){
  speinfor$idmaps <-  c("ACCNUM","ENSEMBL","ENSEMBLPROT","SYMBOL","UNIPROT","REFSEQ","Entrez");} 

 if(species=="Streptomyces coelicolor"){
  speinfor$idmaps <-  c("ACCNUM","ALIAS","SYMBOL","REFSEQ","Entrez");} 

 if(species=="Pig"){
  speinfor$idmaps <-  c("ACCNUM","ALIAS","SYMBOL","UNIPROT","REFSEQ","Entrez");}
 
 if(species=="Xenopus"){
  speinfor$idmaps <-  c("ACCNUM","SYMBOL","UNIPROT","REFSEQ","Entrez");}  
 
 if(species=="E coli strain K12"){
  speinfor$idmaps <-  c("ACCNUM","ALIAS","SYMBOL","REFSEQ","Entrez");}  

 }
 
 speinfor

}

read.net <- function(file){
	net.text <- as.matrix(read.table(file, fill=T, as.is=T, col.names=1:max(count.fields(file))))
	net.text <- toupper(net.text)
	net.node <- unique(as.character(net.text))
	net.node <- net.node[net.node != ""]
	net.size <- length(net.node)
	net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,-1]))
	net.edge <- net.edge[net.edge[,2] != "", ]
	net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
	net.matrix[net.edge] <- 1
	list(size=net.size, node=net.node, matrix=net.matrix)
}

Readme <- function(){
 cat("This is an R package of an implement NOA (Network Ontology Analysis) method. \n In this package, we add new species and gene input ids to NOA, \n  and also give an indenpendent id-mapping module for each species. \n What's more, the result we ouput includes more detailed information \n which are not included in the previous verion of NOA and its webserver. \n  More details can be found by type 'species.information()', 'idmapping.information()', 'result.information()'.", "\n");
 #return;
}

functionrelated <- function(){
# a interior function for the relationship among all the functions used in this files
# NOA (read.file and species.list)
# noa.test					noa2glm					GLM  
# edge.test	read.net								term.attributes
# id.mapping				map.assign		gene.term =>(noa.test and GLM)
# unique.table idmap.table
}

species.information <- function(){
 print("ALL the species information used in NOA package!")	
 data("species")
 print(species)
 print("You can install the corresponding species' annotation packages using the follwing commands.")
 data("dependences")
 print(dependences)
}

idmapping.information <- function(){
 print("ALL the id maps used in NOA package!")	
 data("idmaps")
 print(idmaps)
}

result.information <- function(){
 cat("noaresult: all the information of NOA method, if the input is a gene list, it is NULL. \n
 glmresult: all the information of GLM method. \n
 genepair: the input gene nework edges or gene list of upfile. \n     
 method: 1 means that GLM(gene list method); 2 means that NOA. \n
 geneset: the unique gene set contains in the upfile.
 refgenepair: the reference gene network edges or gene list of reffile. \n 
 If reffile is NULL, for GLM, the refgenepair is all the genes in the species, for NOA, the reference network is the clique network. \n
 refgeneset: the unique gene set contains in the refgenepair.\n
 mapgeneset: the mapped gene set for geneset. In enrichment analysis, we should map input gene id into the annotated gene id. \n
 maprefgeneset: the mapped gene set for refgeneset. \n
 unmappedgene: the gene set that can not map to the annotated id. \n
 term: the enrichment analysis result which is a matrix. Each row is corresponding to one term.\n");
}
if(F){
install.annotation.packages <- function(){
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("org.Sc.sgd.db")
biocLite("org.Dm.eg.db")
biocLite("org.Ce.eg.db")
biocLite("org.Mm.eg.db")
biocLite("org.At.tair.db")
biocLite("org.Rn.eg.db")
biocLite("org.Dr.eg.db")
biocLite("org.Bt.eg.db")
biocLite("org.Cf.eg.db")
biocLite("org.Ag.eg.db")
biocLite("org.EcSakai.eg.db")
biocLite("org.Gg.eg.db")
biocLite("org.Pt.eg.db")
biocLite("org.Pf.plasmo.db")
biocLite("org.Mmu.eg.db")
biocLite("org.Sco.eg.db")
biocLite("org.Ss.eg.db")
biocLite("org.Xl.eg.db")
biocLite("org.EcK12.eg.db")
  }
}