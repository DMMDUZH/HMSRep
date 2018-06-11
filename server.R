rm(list=ls()) 
library(ggplot2)
library(shiny)
#library(plyr)
library(dplyr)
library(Hmisc)
#library(XLConnect)
library(calibrate)
library(RSQLite)
library(VennDiagram)
library(data.table)
library(plotly)
library(RColorBrewer)
library(ComplexHeatmap)
library(gplots)
library(DescTools)
library(heatmaply)
library(clusterProfiler)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
library(DOSE)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(lettercase)

shinyServer(function(input, output,session) {
	options(shiny.maxRequestSize=10000*1024^2)
	outersect <- function(x, y, ...) {
  		big.vec <- c(x, y, ...)
  		duplicates <- big.vec[duplicated(big.vec)]
  		setdiff(big.vec, unique(duplicates))
	}
	#This function computes GO enrichment
	goEnrich <- function(globalProt, localProt, organism, goterm){
		#function(globalProt, localProt, geneList, organism){
		egy <- 0
		egx <- 0
		orgDB <- 0
		y.prot <- globalProt
		x.prot <- localProt
		if(organism == "Mouse"){
		egy <- bitr(y.prot, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
		egx <- bitr(x.prot, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
		orgDB <- "org.Mm.eg.db"
		}
		else if(organism == "Human"){
		print("y.prot")
		egy <- bitr(y.prot, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		print("x.prot")
                egx <- bitr(x.prot, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		orgDB <- "org.Hs.eg.db"
		}
		else{
			print("No organim is selected")
		}
		#names(geneList) <- plyr::mapvalues(names(geneList), from = egx$SYMBOL, to = egx$ENTREZID)
		gene <- egx$ENTREZID
		uniList <- egy$ENTREZID
		goEnDataGT <- enrichGoGT(gene, uniList, orgDB, goterm)
                return(goEnDataGT)
	}
	#This function computes GO enrichment
	enrichGoGT <- function(gene, uniList, orgDB, goterm){
                ego.gt <- enrichGO(gene    = gene,
                universe      = uniList,
                OrgDb         = orgDB,
                ont           = goterm,
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                minGSSize = 3,
                readable      = TRUE)

                ego.gt.r <- ego.gt[order(ego.gt$ID),]
                ego.gt.r <- ego.gt.r[!duplicated(ego.gt.r$ID),]
                ego.gt.ro <- ego.gt.r[order(ego.gt.r$p.adjust),]
                goEnDataGT <- ego.gt.ro
                return(goEnDataGT)
	}
	#Stack bar preparation
	freqDistrCalc <- function(stackBart, expType){
		stackBar1 <- data.frame(matrix(ncol = 3, nrow = length(stackBart$groups)))
                colnames(stackBar1) <- c("Freq1", "pep_var_mod_res", "groups")
		stackBar1$Freq1 <- stackBart$Freq
		stackBar1$pep_var_mod_res <- stackBart$pep_var_mod_res
		stackBar1$groups <- stackBart$groups
		reshaped <- reshape(stackBar1, varying = 1, sep = "", direction = 'long')
		reshaped$cat <- ''
		reshaped[reshaped$time ==1,]$cat <- expType
		myBlues <- c("dodgerblue7", "blue")
		#reshaped$cat <- factor(reshaped$cat, levels = reshaped$cat[order(reshaped$time)])
		return(reshaped)
	}
	#Function that is not in use at the moment
	pepVarModConf <- function(locConfTab, pep_var_mod_res, groups){
		pepVarModConf <- rename(count(locConfTab, pep_var_mod_res, groups), Freq = n)
		pepVarModConf <- data.frame(pepVarModConf)
		pepVarModConf$groups[pepVarModConf$groups == 1] <- "< 60%"
		pepVarModConf$groups[pepVarModConf$groups == 2] <- "> 60% and < 95%"
		pepVarModConf$groups[pepVarModConf$groups == 3] <- ">= 95%"
		return(pepVarModConf)
	}
#Initializing some of the dataframes needed
msmsFiles1 <- list()
msmsFiles2 <- list()
stackBar <- list()
msmsFiles3 <- list()
msmsProts1 <- list()
msmsPep1 <- list()
msmsPep11 <- list()
msmsPep12 <- list()
msmsHeat1 <- list()
msmsHeat2 <- list()
msmsLocConf1 <- list()
msmsLocConf2 <- list()
msmsPepScore1 <- list()
msmsPepScore2 <- list()
msmsProts2 <- list()
msmsProts3 <- list()
msmsPep2 <- list()
msmsPep3 <- list()
msmsPep22 <- list()
msmsPep32 <- list()
msmsPsm1 <- list()
msmsPsm2 <- list()
msmsPsm3 <- list()
msmsFiles4 <- list()
msmsFiles5 <- list()
msmsFiles6 <- list()
msmsFiles7 <- list()
msmsFiles8 <- list()
msmsFiles51 <- list()
msmsFiles81 <- list()
msmsFiles41 <- list()
	
	#Event reactive that reads the uploaded files and store the number of files loaded in a list
	#This also stores the number of the files
	csvdata1 <- eventReactive(input$choice,{
    		infile1 <- input$csvfile1
		req(infile1)
    		if (is.null(infile1)) {
      			return(NULL)
    		}
		else {
		numfiles <<- nrow(infile1)
        		for (i in 1:numfiles){
				#Reading csv files to the data objects in the list
				msmsFiles1[[i]] <<- read.csv(input$csvfile1[[i, 'datapath']], header = input$header, sep = input$sep, quote = input$quote, fill = TRUE)
				#Below is a testing code that has no effect on the HMSRep script
				msmsName <- paste("msms1.", i, sep = "")
				assign(msmsName, read.csv(input$csvfile1[[i, 'datapath']], header = input$header, sep = input$sep, quote = input$quote, fill = TRUE))
			}
		}
	})
	#Reactive that is not neessary for the function of the script
	create_tabs1 <- reactive({
		dataTab1 <- csvdata1()
                for (i in 1:numfiles){
                        msmsName2 <- paste("msms2.", i, sep = "")
			msmsFiles2[[i]] <<-  msmsFiles1[[i]]
		}
        })
	#Reactive that does most of the work in this script
	create_tabs3 <- reactive({
		organism <- input$organism
		if(organism == "Mouse"){
		globalProt <<- read.table("mm_universal_list.txt", header = TRUE, sep="\t")
		}
		else if(organism == "Human"){
		globalProt <<- read.table("hs_universal_list.txt", header = TRUE, sep="\t")
		}
		msmsFiles <- csvdata1()
                searchVar <- input$searchtype
                searchOrderVar <- input$searchorder
                for(i in 1:numfiles){
                msms4.2 <- msmsFiles1[[i]]
		#Initial filtering to accept the records having mascot score more than or equal to 20 and that there is a protein name allocated for the record
                msms4.2.1 <- msms4.2[msms4.2$pep_score >= 20 & msms4.2$prot_acc !="",]
		#Remove quotes from modification column if there is any
		msms4.2.1$pep_var_mod <- gsub("\"","",msms4.2.1$pep_var_mod)
		#Filter only the search type user selected
		msms4.2.2 <- msms4.2.1[msms4.2.1$pep_var_mod == searchVar,]
		
		msms4.2.3 <- msms4.2.2[FALSE,]
		#Remove corrupted records
		if("X.query_number." %in% colnames(msms4.2.2)){
			msms4.2.3 <- msms4.2.2[msms4.2.2$pep_expect <= 0.05 & msms4.2.2$X.query_number. != "emPAI",]
		}
		else{
			msms4.2.3 <- msms4.2.2[msms4.2.2$pep_expect <= 0.05 & msms4.2.2$query_number != "emPAI",]
		}	
		#If modification confidence cell is empty write 99.99% in the cell
		msms4.2.3$pep_var_mod_conf <- ifelse(msms4.2.3$pep_var_mod_conf == "", "99.99%", as.character(msms4.2.3$pep_var_mod_conf)) 
		#Create a numeric column out of the localization confidence column
		msms4.2.3$pep_var_mod_conf.1 <- as.numeric(as.character(sapply(strsplit(as.character(msms4.2.3$pep_var_mod_conf), "%"), `[`, 1)))

		msms4.2.4 <- msms4.2.3[FALSE,]
		msms4.2.5 <- msms4.2.3[FALSE,]

		#Sort and remove duplicate querries depending on the csv species
		if("X.query_number." %in% colnames(msms4.2.3)){
			msms4.2.4 <- msms4.2.3[order(msms4.2.3$X.query_number., -msms4.2.3$pep_var_mod_conf.1),]
			msms4.2.5 <- msms4.2.4[!duplicated(msms4.2.4$X.query_number.),]
		}
		else{
			msms4.2.4 <- msms4.2.3[order(msms4.2.3$query_number, -msms4.2.3$pep_var_mod_conf.1),]
                        msms4.2.5 <- msms4.2.4[!duplicated(msms4.2.4$query_number),]
		}

		msmsFiles3[[i]] <<- msms4.2.5

		#Extract gene symbol and store in a separate column
		ifelse(str_count(as.character(msmsFiles3[[i]]$prot_acc), "\\|")>1, protName <- sapply(strsplit(as.character(msmsFiles3[[i]]$prot_acc), '\\|'), '[', 3), protName <- sapply(strsplit(as.character(msmsFiles3[[i]]$prot_acc), '\\|'), '[', 2))
		msmsFiles3[[i]]$protName <<- sapply(strsplit(as.character(protName), '_'), '[', 1)
		#remove contaminats
		msmsFiles3[[i]] <<- msmsFiles3[[i]][msmsFiles3[[i]]$protName != "ZZ" & msmsFiles3[[i]]$protName != '"',]

		#Further process gene symbol column just created
		msmsFiles3[[i]]$protName1 <<- gsub(".*GN=","",msmsFiles3[[i]]$prot_desc)
		msmsFiles3[[i]]$protName1 <<- gsub("\\ PE=.*","",msmsFiles3[[i]]$protName1)
		msmsFiles3[[i]]$protName1 <<- gsub("\"", "", msmsFiles3[[i]]$protName1, fixed=TRUE)

		#Find exact position of the modification residue within each peptide		
		pos <- (unlist(gregexpr(pattern =searchOrderVar,msmsFiles3[[i]]$pep_var_mod_pos.1))-2)

		n <- length(msmsFiles3[[i]]$pep_seq.1)

		#Next we try to collapse duplicated peptides
		msmsFiles3[[i]]$pep_seq.2 <<- gsub("\\(|\\)", "",msmsFiles3[[i]]$pep_seq.1)
		msmsFiles3[[i]]$pep_seq.by <<- gsub("\\(|\\)", "",msmsFiles3[[i]]$pep_seq.1)
		msmsFiles3[[i]]$pep_seq.ds <<- gsub("\\(|\\)", "",msmsFiles3[[i]]$pep_seq.1)
		msmsFiles3[[i]]$pep_seq.us <<- gsub("\\(|\\)", "",msmsFiles3[[i]]$pep_seq.1)

		msmsFiles3[[i]] <<- msmsFiles3[[i]][order(nchar(msmsFiles3[[i]]$pep_seq.2), -msmsFiles3[[i]]$pep_var_mod_conf.1),]

		msmsFiles4[[i]] <<- msmsFiles3[[i]]

for(k in 1:length(msmsFiles3[[i]]$pep_seq.2)){
        pepSub1 <- paste0(".*(", msmsFiles3[[i]]$pep_seq.2[[k]], ").*")
        pepSub2 <- paste0(msmsFiles3[[i]]$pep_seq.2[[k]], ".*")
        pepSub3 <- paste0(".*", msmsFiles3[[i]]$pep_seq.2[[k]])
        for(j in 1:length(msmsFiles3[[i]]$pep_seq.2)){
                msmsFiles4[[i]]$pep_seq.by[[j]] <<- gsub(pepSub1,"\\1",as.character(msmsFiles4[[i]]$pep_seq.by[[j]]))
                msmsFiles4[[i]]$pep_seq.us[[j]] <<- gsub(pepSub2,"",as.character(msmsFiles4[[i]]$pep_seq.us[[j]]))
                msmsFiles4[[i]]$pep_seq.ds[[j]] <<- gsub(pepSub3,"",as.character(msmsFiles4[[i]]$pep_seq.ds[[j]]))
        }
}

	msmsFiles4[[i]]$pep_seq.2.len <<- str_length(msmsFiles4[[i]]$pep_seq.2)
	msmsFiles4[[i]]$pep_seq.by.len <<- str_length(msmsFiles4[[i]]$pep_seq.by)
	msmsFiles4[[i]]$pep_seq.us.len <<- str_length(msmsFiles4[[i]]$pep_seq.us)
	msmsFiles4[[i]]$pep_seq.ds.len <<- str_length(msmsFiles4[[i]]$pep_seq.ds)

	msmsFiles4[[i]]$pep_var_mod_pos_num <<- (unlist(gregexpr(pattern =searchOrderVar,msmsFiles4[[i]]$pep_var_mod_pos.1))-2)
	msmsFiles4[[i]]$pep_var_mod_res <<- substring(msmsFiles4[[i]]$pep_seq.2, msmsFiles4[[i]]$pep_var_mod_pos_num, msmsFiles4[[i]]$pep_var_mod_pos_num)
	msmsFiles4[[i]]$pep_var_mod_res_num <<- paste0( msmsFiles4[[i]]$pep_var_mod_res, msmsFiles4[[i]]$pep_var_mod_pos_num)
	#msmsFiles4[[i]]$prot_var_mod_res_num <<- (unlist(regexpr(pattern=msmsFiles4[[i]]$pep_seq.2, msmsFiles4[[i]]$prot_seq, perl=TRUE)) + (msmsFiles4[[i]]$pep_var_mod_pos_num -1))

	msmsFiles4[[i]]$prot_var_mod_pos_num <<- mapply(function(x, y, z) (unlist(regexpr(pattern=x, y, perl=TRUE)) + (z - 1)), msmsFiles4[[i]]$pep_seq.2, msmsFiles4[[i]]$prot_seq, msmsFiles4[[i]]$pep_var_mod_pos_num)

	msmsFiles4[[i]]$prot_var_mod_res_num <<- paste0( msmsFiles4[[i]]$pep_var_mod_res, msmsFiles4[[i]]$prot_var_mod_pos_num)

	#msmsFiles4[[i]]$prot_var_mod_res_num <<- StrPos(msmsFiles4[[i]]$prot_seq, pattern=msmsFiles4[[i]]$pep_seq.2, pos = 1)

	msmsFiles4[[i]]$pep_var_mod_res_site <<- msmsFiles4[[i]]$pep_var_mod_pos_num

for(k in 1:length(msmsFiles4[[i]]$pep_seq.us.len)){
	if(msmsFiles4[[i]]$pep_seq.us.len[[k]] > 0 & msmsFiles4[[i]]$pep_seq.us.len[[k]] >= msmsFiles4[[i]]$pep_var_mod_pos_num[[k]]){
        	msmsFiles4[[i]]$pep_var_mod_res_site[[k]] <<- msmsFiles4[[i]]$pep_var_mod_pos_num[[k]]
#		msmsFiles4[[i]]$prot_var_mod_res_num[[k]] <<- (unlist(gregexpr(pattern=msmsFiles4[[i]]$pep_seq.2[[k]], msmsFiles4[[i]]$prot_seq[[k]])) + (msmsFiles4[[i]]$pep_var_mod_pos_num[[k]] -1))
        	}
	else if(msmsFiles4[[i]]$pep_seq.us.len[[k]] > 0 & msmsFiles4[[i]]$pep_seq.us.len[[k]] < msmsFiles4[[i]]$pep_var_mod_pos_num[[k]]){
        	msmsFiles4[[i]]$pep_var_mod_res_site[[k]] <<- msmsFiles4[[i]]$pep_var_mod_pos_num[[k]] - msmsFiles4[[i]]$pep_seq.us.len[[k]]
        	}
	else {
        	msmsFiles4$pep_var_mod_res_site[[k]] <<- msmsFiles4$pep_var_mod_pos_num[[k]]
        }

}

	cx1 <- grep("pep_seq.by", colnames(msmsFiles4[[i]]))

	cx2 <- grep("pep_var_mod_res_site", colnames(msmsFiles4[[i]]))

	cx3 <- grep("pep_seq.2", colnames(msmsFiles4[[i]]))

	msmsFiles5[[i]] <<- msmsFiles4[[i]][!duplicated(msmsFiles4[[i]][c(cx1,cx2)]),]

	#msmsFiles41[[i]] <<- with(msmsFiles4[[i]], ave(as.numeric(S, Cnty, FUN=function(x) length(unique(x)))) 	

	msmsProts1[[i]] <<- msmsFiles5[[i]][!duplicated(msmsFiles5[[i]]$protName1),]

	msmsFiles6.t <- data.frame(protName=msmsFiles5[[i]]$protName1, pepSeq=msmsFiles5[[i]]$pep_seq.2, pepSeq.upstream=msmsFiles5[[i]]$pep_seq.us, pepSeq.stem=msmsFiles5[[i]]$pep_seq.by,pepSeq.downstream=msmsFiles5[[i]]$pep_seq.ds,modRes=msmsFiles5[[i]]$pep_var_mod_res_num, protPos=msmsFiles5[[i]]$prot_var_mod_res_num, protSeq=msmsFiles5[[i]]$prot_seq) 

	msmsFiles6[[i]] <<- msmsFiles6.t

	msmsFiles7[[i]] <<- msmsFiles4[[i]][!duplicated(msmsFiles4[[i]][c(cx3,cx2)]),]

	msmsFiles8.t <- data.frame(protName=msmsFiles7[[i]]$protName1, pepSeq=msmsFiles7[[i]]$pep_seq.2, pepSeq.upstream=msmsFiles7[[i]]$pep_seq.us, pepSeq.stem=msmsFiles7[[i]]$pep_seq.by,pepSeq.downstream=msmsFiles7[[i]]$pep_seq.ds,modRes=msmsFiles7[[i]]$pep_var_mod_res_num, protPos=msmsFiles7[[i]]$prot_var_mod_res_num, protSeq=msmsFiles7[[i]]$prot_seq)

	msmsFiles8[[i]] <<- msmsFiles8.t

                msmsPep11[[i]] <<- msmsFiles5[[i]]
		msmsPep12[[i]] <<- msmsFiles7[[i]]

		msmsHeat1[[i]] <<- data.frame(protName=character(n), sample=as.numeric(n))
                msmsHeat1[[i]]$protName  <<- as.character(msmsFiles3[[i]]$protName1)
		msmsHeat1[[i]]$sample <<- rep(i,nrow(msmsHeat1[[i]]))
		msmsHeat2[[i]] <<- as.data.frame(table(msmsHeat1[[i]][,1]))
		colnames(msmsHeat2[[i]])[1] <<- "protName"
		colnames(msmsHeat2[[i]])[2] <<- as.character(msmsHeat1[[i]][1,2])
		locConfTab <- data.frame(pep_seq=character(n),pep_score=as.numeric(n), pep_var_mod_conf=as.numeric(n), pep_var_mod_pos=as.numeric(n), pep_var_mod_pos_num=as.numeric(n), pep_var_length=as.numeric(n),pep_var_mod_res=character(n),stringsAsFactors=FALSE)
		pep.seq <- as.character(msmsFiles3[[i]]$pep_seq)
		x <- as.character(msmsFiles3[[i]]$pep_var_mod_conf)
		pep.var.mod.conf <- as.numeric(sapply(strsplit(as.character(x), "%"), `[`, 1))
		#msmsFiles3[[i]]$pep_var_mod_conf.n <- as.numeric(sapply(strsplit(as.character(msmsFiles3[[i]]$pep_var_mod_conf), "%"), `[`, 1))

		msmsLocConf1[[i]] <<- data.frame(pep_seq=character(n),pep_score=as.numeric(n), pep_var_mod_conf=as.numeric(n))
		msmsLocConf1[[i]]$pep_seq <<- msmsFiles3[[i]]$pep_seq
		msmsLocConf1[[i]]$pep_score <<- msmsFiles3[[i]]$pep_score
		msmsLocConf1[[i]]$pep_var_mod_conf <<- msmsFiles3[[i]]$pep_var_mod_conf.1
		msmsLocConf1[[i]]$group <<- cut(msmsLocConf1[[i]]$pep_var_mod_conf,seq(from = 0, to = 100, by = 10))
		
		msmsLocConf11 <- as.data.frame(table(msmsLocConf1[[i]][,4]))
		colnames(msmsLocConf11)[1] <- "locConf"
                colnames(msmsLocConf11)[2] <- "psmCounts"
		msmsLocConf2[[i]] <<- msmsLocConf11[msmsLocConf11$psmCounts >=1,]

		n1 <- length(msmsFiles5[[i]]$pep_seq.2)
		
		msmsProts2[[i]] <<- data.frame(protName=character(n), pep_seq=character(n),pep_score=as.numeric(n),pep_var_mod_conf=as.numeric(n))
                msmsProts2[[i]]$pep_seq <<- msmsFiles3[[i]]$pep_seq
		msmsProts2[[i]]$protName <<- msmsFiles3[[i]]$protName1
                msmsProts2[[i]]$pep_score <<- msmsFiles3[[i]]$pep_score
                msmsProts2[[i]]$pep_var_mod_conf <<- msmsFiles3[[i]]$pep_var_mod_conf.1

                msmsProts21 <- as.data.frame(table(msmsProts2[[i]][,1]))
                colnames(msmsProts21)[1] <- "protName"
                colnames(msmsProts21)[2] <- "pepCount"
                msmsProts3[[i]] <<- msmsProts21

		n5 <- length(msmsFiles5[[i]]$pep_seq.1)

		msmsPep2[[i]] <<- data.frame(pep_seq=character(n5),protName=character(n5), pep_score=as.numeric(n5),pep_var_mod_conf=as.numeric(n5))
                msmsPep2[[i]]$pep_seq <<- msmsFiles5[[i]]$pep_seq
                msmsPep2[[i]]$protName <<- msmsFiles5[[i]]$protName1
                msmsPep2[[i]]$pep_score <<- msmsFiles5[[i]]$pep_score
                msmsPep2[[i]]$pep_var_mod_conf <<- msmsFiles5[[i]]$pep_var_mod_conf.1

                msmsPep21 <- as.data.frame(table(droplevels(msmsPep2[[i]][,1])))
                colnames(msmsPep21)[1] <- "pepSeq"
                colnames(msmsPep21)[2] <- "pepCount"
                msmsPep3[[i]] <<- msmsPep21

		n7 <- length(msmsFiles7[[i]]$pep_seq.1)

		msmsPep22[[i]] <<- data.frame(pep_seq=character(n7),protName=character(n7), pep_score=as.numeric(n7),pep_var_mod_conf=as.numeric(n7))
                msmsPep22[[i]]$pep_seq <<- msmsFiles7[[i]]$pep_seq
                msmsPep22[[i]]$protName <<- msmsFiles7[[i]]$protName1
                msmsPep22[[i]]$pep_score <<- msmsFiles7[[i]]$pep_score
                msmsPep22[[i]]$pep_var_mod_conf <<- msmsFiles7[[i]]$pep_var_mod_conf.1

                msmsPep212 <- as.data.frame(table(droplevels(msmsPep22[[i]][,1])))
                colnames(msmsPep212)[1] <- "pepSeq"
                colnames(msmsPep212)[2] <- "pepCount"
                msmsPep32[[i]] <<- msmsPep212

		msmsPepScore1[[i]] <<- data.frame(pep_seq=character(n),pep_score=as.numeric(n), pep_var_mod_conf=as.numeric(n))
		msmsPepScore1[[i]]$pep_seq <<- msmsFiles3[[i]]$pep_seq
                msmsPepScore1[[i]]$pep_score <<- msmsFiles3[[i]]$pep_score
                msmsPepScore1[[i]]$pep_var_mod_conf <<- msmsFiles3[[i]]$pep_var_mod_conf.1
		msmsPepScore1[[i]]$group <<- cut(msmsPepScore1[[i]]$pep_score, seq(from = 20, to = max(msmsPepScore1[[i]]$pep_score), by = 10))

                msmsPepScore11 <- as.data.frame(table(msmsPepScore1[[i]][,4]))
                colnames(msmsPepScore11)[1] <- "pepScore"
               	colnames(msmsPepScore11)[2] <- "psmCounts"
                msmsPepScore2[[i]] <<- msmsPepScore11[msmsPepScore11$psmCounts >=1,]

		

		locConfTab$pep_seq <- msmsFiles3[[i]]$pep_seq.2
		locConfTab$pep_score <- msmsFiles3[[i]]$pep_score
		locConfTab$pep_var_mod_conf <- msmsFiles3[[i]]$pep_var_mod_conf.1
		locConfTab$pep_var_mod_pos <- msmsFiles3[[i]]$pep_var_mod_pos.1
		locConfTab$pep_var_mod_pos_num <- (unlist(gregexpr(pattern =searchOrderVar,as.character(msmsFiles3[[i]]$pep_var_mod_pos.1)))-2)
		locConfTab$pep_var_length <- lapply(pep.seq, function(x) max(nchar(x)))
		locConfTab$pep_var_mod_res <- substring(locConfTab$pep_seq, locConfTab$pep_var_mod_pos_num, locConfTab$pep_var_mod_pos_num)
		locConfTab$groups <- as.numeric(cut2(locConfTab$pep_var_mod_conf, g=3))
		locConfTab$groups[locConfTab$pep_var_mod_conf < 60 ] <- 1
		locConfTab$groups[locConfTab$pep_var_mod_conf > 60 & locConfTab$pep_var_mod_conf  < 95] <- 2
		locConfTab$groups[locConfTab$pep_var_mod_conf > 95 ] <- 3
		locConfTab1 <- locConfTab
		pepVarModConf <- dplyr::rename(count(locConfTab, pep_var_mod_res, groups), Freq = n)
		pepVarModConf <- data.frame(pepVarModConf)
		pepVarModConf$groups[pepVarModConf$groups == 1] <- "< 60%"
		pepVarModConf$groups[pepVarModConf$groups == 2] <- "> 60% and < 95%"
		pepVarModConf$groups[pepVarModConf$groups == 3] <- "> 95%"
		stackBar[[i]] <<- data.frame(matrix(ncol = 3, nrow = length(pepVarModConf$groups)))
		colnames(stackBar[[i]]) <<- c("groups", "Freq", "pep_var_mod_res")
		stackBar[[i]]$groups <<- pepVarModConf$groups
		stackBar[[i]]$Freq <<- pepVarModConf$Freq
		stackBar[[i]]$pep_var_mod_res <<- pepVarModConf$pep_var_mod_res
		}
        })
	
	observe({
		 updateCheckboxGroupInput(session, "select1",
      			label = "Columns in msms_file1 to show",
      			choices = names(create_tabs1()),
      			selected = names(create_tabs1())
    			)
		})

	output$mytable1 <- renderDataTable({
		testre<-create_tabs1()
		columns <- names(msmsFiles2[[1]])
      		if (!is.null(input$select1)) {
        		columns <- input$select1
      		}
   		msmsFiles2[[1]][,columns,drop=FALSE] 
  	},caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)
	        output$mytable2 <- renderDataTable({
		testre<-create_tabs1()
		if(numfiles>=2){
                columns <- names(msmsFiles2[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles2[[2]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)
		output$mytable3 <- renderDataTable({
                testre<-create_tabs1()
                if(numfiles>=3){
		columns <- names(msmsFiles2[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles2[[3]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)
		output$mytable4 <- renderDataTable({
                testre<-create_tabs1()
                if(numfiles>=4){
		columns <- names(msmsFiles2[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles2[[4]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$mytable11 <- renderDataTable({
		testre<-create_tabs3()
		if(numfiles>=1){
                columns <- names(msmsProts3[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[1]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$mytable21 <- renderDataTable({
		testre<-create_tabs3()
                if(numfiles>=2){
		columns <- names(msmsProts3[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[2]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable31 <- renderDataTable({
		testre<-create_tabs3()
                if(numfiles>=3){
		columns <- names(msmsProts3[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[3]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable41 <- renderDataTable({
		testre<-create_tabs3()
		if(numfiles>=4){
                columns <- names(msmsProts3[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[4]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)
	output$mytable12 <- renderDataTable({
		testre<-create_tabs3()
                columns <- names(stackBar[[1]])
                stackBar[[1]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$mytable13 <- renderDataTable({
                testre<-create_tabs3()
                columns <- names(stackBar[[1]])
                stackBar[[1]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable22 <- renderDataTable({
		testre<-create_tabs3()
                columns <- names(stackBar[[2]])
                stackBar[[2]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable32 <- renderDataTable({
		testre<-create_tabs3()
                columns <- names(stackBar[[3]])
                stackBar[[3]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable42 <- renderDataTable({
		testre<-create_tabs3()
                columns <- names(stackBar[[4]])
                stackBar[[4]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable1 <- DT::renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles6[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                DT::datatable(msmsFiles6[[1]][,columns,drop=FALSE],options = list(searchHighlight = TRUE, search = list(search = 'YQKSTELLIR')))
		#formatStyle('Sepal.Width',color = styleInterval(c(3.4, 3.8), c('white', 'blue', 'red')),backgroundColor = styleInterval(3.4, c('gray', 'yellow')))
		#,options = list(searchHighlight = TRUE, search = list(search = msmsFiles6[[1]]$pep_seq.2)))
		#[,columns,drop=FALSE]
		#,options = list(searchHighlight = TRUE, search = list(search = msmsFiles6[[1]]$pep_seq.2)))
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable11 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles8[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles8[[1]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable2 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles6[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles6[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable21 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles8[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles8[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable3 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles6[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles6[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable31 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles8[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles8[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable4 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles6[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles6[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable41 <- renderDataTable({
                testre<-create_tabs3()
                if(numfiles>=1){
                columns <- names(msmsFiles8[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles8[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$freqdistr1 <- renderPlotly({
		expType1 <- input$exptype1
		testre<-create_tabs3()
		stackBart <- stackBar[[1]]
		reshaped <- freqDistrCalc(stackBart, expType1)
		p <- ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
  		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
  		xlab("ADPr assigned residues") +
  		ylab("Number of ADPr assignments") +
  		labs(title = "mascot localization confidence") +
  		scale_fill_manual(values=c("grey68", "grey48", "black")) +
  		theme(panel.background = element_rect(fill = "white"),
  		panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
  		panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
		axis.text.x = element_text(colour="white"))
		p1 <- ggplotly(p)
		p1$elementId <- NULL
		print(p1)		
	})

		output$psmLocConf1 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsLocConf2[[1]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))

	})

		output$psmPepScore1 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsPepScore2[[1]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmPepScore2 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsPepScore2[[2]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmPepScore3 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsPepScore2[[3]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmPepScore4 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsPepScore2[[4]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
			
		output$psmLocConf2 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsLocConf2[[2]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
		geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmLocConf3 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsLocConf2[[3]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
		geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmLocConf4 <- renderPlot({
                testre<-create_tabs3()
                ggplot(data=msmsLocConf2[[4]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
		geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })

		output$freqdistr12 <- renderPlot({
		expType1 <- input$exptype1
                testre<-create_tabs3()
                stackBart <- stackBar[[1]]
                reshaped <- freqDistrCalc(stackBart, expType1)
                ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
                xlab("ADPr assigned residues") +
                ylab("Number of ADPr assignments") +
                labs(title = "mascot localization confidence") +
                scale_fill_manual(values=c("grey68", "grey48", "black")) +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
output$freqdistr2 <- renderPlotly({
		expType2 <- input$exptype2
		testre<-create_tabs3()
		stackBart <- stackBar[[2]]
		reshaped <- freqDistrCalc(stackBart, expType2)
		p <- ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
  		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
  		xlab("ADPr assigned residues") +
  		ylab("Number of ADPr assignments") +
  		labs(title = "mascot localization confidence") +
  		scale_fill_manual(values=c("grey68", "grey48", "black")) +
  		theme(panel.background = element_rect(fill = "white"),
  		panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
  		panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
		axis.text.x = element_text(colour="white"))
		p1 <- ggplotly(p)
		p1$elementId <- NULL
		print(p1)		
	})

		output$freqdistr22 <- renderPlot({
		expType2 <- input$exptype2
                testre<-create_tabs3()
                stackBart <- stackBar[[2]]
		reshaped <- freqDistrCalc(stackBart, expType2)
                ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
                xlab("ADPr assigned residues") +
                ylab("Number of ADPr assignments") +
                labs(title = "mascot localization confidence") +
                scale_fill_manual(values=c("grey68", "grey48", "black")) +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
output$freqdistr3 <- renderPlotly({
		expType3 <- input$exptype3
		testre<-create_tabs3()
		stackBart <- stackBar[[3]]
		reshaped <- freqDistrCalc(stackBart, expType3)
		p <- ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
  		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
  		xlab("ADPr assigned residues") +
  		ylab("Number of ADPr assignments") +
  		labs(title = "mascot localization confidence") +
  		scale_fill_manual(values=c("grey68", "grey48", "black")) +
  		theme(panel.background = element_rect(fill = "white"),
  		panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
  		panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
		axis.text.x = element_text(colour="white"))
		p1 <- ggplotly(p)
		p1$elementId <- NULL
		print(p1)		
	})

		output$freqdistr32 <- renderPlot({
		expType3 <- input$exptype3
                testre<-create_tabs3()
                stackBart <- stackBar[[3]]
		reshaped <- freqDistrCalc(stackBart, expType3)
                ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
                xlab("ADPr assigned residues") +
                ylab("Number of ADPr assignments") +
                labs(title = "mascot localization confidence") +
                scale_fill_manual(values=c("grey68", "grey48", "black")) +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
	output$freqdistr4 <- renderPlotly({
		expType4 <- input$exptype4
		testre<-create_tabs3()
		stackBart <- stackBar[[4]]
		reshaped <- freqDistrCalc(stackBart, expType4)
		p <- ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
  		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
  		xlab("ADPr assigned residues") +
  		ylab("Number of ADPr assignments") +
  		labs(title = "mascot localization confidence") +
  		scale_fill_manual(values=c("grey68", "grey48", "black")) +
  		theme(panel.background = element_rect(fill = "white"),
  		panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
  		panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
		axis.text.x = element_text(colour="white"))
		p1 <- ggplotly(p)
		p1$elementId <- NULL
		print(p1)		
	})

		output$freqdistr42 <- renderPlot({
		expType4 <- input$exptype4
                testre<-create_tabs3()
                stackBart <- stackBar[[4]]
		reshaped <- freqDistrCalc(stackBart, expType4)
                ggplot(reshaped, aes(x = cat, y = Freq, fill = groups)) +
		geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ pep_var_mod_res) +
                xlab("ADPr assigned residues") +
                ylab("Number of ADPr assignments") +
                labs(title = "mascot localization confidence") +
                scale_fill_manual(values=c("grey68", "grey48", "black")) +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "black"),
                panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
	output$vennDiaProt <- renderPlot({
	expType1 <- input$exptype1
	expType2 <- input$exptype2
	expType3 <- input$exptype3
	expType4 <- input$exptype4
	testre<-create_tabs3()
	if(numfiles==4){
          	msmsProts11<-msmsProts1[[1]]$protName
		msmsProts21<-msmsProts1[[2]]$protName
		msmsProts31<-msmsProts1[[3]]$protName 
          	msmsProts41<-msmsProts1[[4]]$protName                
		
		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
          	InterSec13<-Reduce(intersect, list(msmsProts11,msmsProts31))
          	InterSec14<-Reduce(intersect, list(msmsProts11,msmsProts41))
         	InterSec23<-Reduce(intersect, list(msmsProts21,msmsProts31))
          	InterSec24<-Reduce(intersect, list(msmsProts21,msmsProts41))
          	InterSec34<-Reduce(intersect, list(msmsProts31,msmsProts41))
          	InterSec123<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31))
          	InterSec124<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts41))
          	InterSec134<-Reduce(intersect, list(msmsProts11,msmsProts31,msmsProts41))
          	InterSec234<-Reduce(intersect, list(msmsProts21,msmsProts31,msmsProts41))
          	InterSec1234<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31,msmsProts41))
		venn.plot <- draw.quad.venn(
          	area1 = length(msmsProts1[[1]]$protName),
		area2 = length(msmsProts1[[2]]$protName),
          	area3 = length(msmsProts1[[3]]$protName),
		area4 = length(msmsProts1[[4]]$protName),
                        n12 = length(InterSec12),
                        n13 = length(InterSec13),
                        n14 = length(InterSec14),
                        n23 = length(InterSec23),
                        n24 = length(InterSec24),
                        n34 = length(InterSec34),
                        n123 = length(InterSec123),
                        n124 = length(InterSec124),
                        n134 = length(InterSec134),
                        n234 = length(InterSec234),
                        n1234 = length(InterSec1234),
			category = c("A", "B", "C", "D"),
                        fill = c("orange", "green", "darkorchid1", "cornflowerblue"),
                        lty = "blank",
                        col = c("orange", "cornflowerblue", "darkorchid1", "green"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen", "darkorchid1", "cornflowerblue")
                )
                mainT = "Overalps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3," \nD = ",expType4,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	if(numfiles==3){
          msmsProts11<-msmsProts1[[1]]$protName
		msmsProts21<-msmsProts1[[2]]$protName
		msmsProts31<-msmsProts1[[3]]$protName 
		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
          InterSec13<-Reduce(intersect, list(msmsProts11,msmsProts31))
         	InterSec23<-Reduce(intersect, list(msmsProts21,msmsProts31))
          InterSec123<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31))
		venn.plot <- draw.triple.venn(
          area1 = length(msmsProts1[[1]]$protName),
		area2 = length(msmsProts1[[2]]$protName),
          area3 = length(msmsProts1[[3]]$protName),
                        n12 = length(InterSec12),
                        n13 = length(InterSec13),
                        n23 = length(InterSec23),
                        n123 = length(InterSec123),
			category = c("A", "B", "C"),
                        fill = c("orange", "green", "darkorchid1"),
                        lty = "blank",
                        col = c("orange", "cornflowerblue", "darkorchid1"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen", "darkorchid1")
                )
                mainT = "Overalps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3, sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	if(numfiles==2){
          	msmsProts11<-msmsProts1[[1]]$protName
		msmsProts21<-msmsProts1[[2]]$protName
		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
		venn.plot <- draw.pairwise.venn(
          area1 = length(msmsProts1[[1]]$protName),
		area2 = length(msmsProts1[[2]]$protName),
                        cross.area = length(InterSec12),
			category = c("A", "B"),
                        fill = c("orange", "green"),
                        lty = "blank",
                        col = c("orange", "cornflowerblue"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen")
                )
                mainT = "Overalps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	if(numfiles==1){
                msmsProts11<-msmsProts1[[1]]$protName
                msmsProts21<-msmsProts1[[2]]$protName
                InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
                venn.plot <- draw.pairwise.venn(
          	area1 = length(msmsProts1[[1]]$protName),
                area2 = 0,
                        cross.area = 0,
                        category = c("A", "B"),
                        fill = c("orange", "green"),
                        lty = "blank",
                        col = c("orange", "cornflowerblue"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen")
                )
                mainT = "Overalps of ADP-ribosylated proteins"
                bottomT = paste("A = ",expType1,sep="")
                require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
        })
	output$vennDiaPep <- renderPlot({
	expType1 <- input$exptype1
        expType2 <- input$exptype2
        expType3 <- input$exptype3
        expType4 <- input$exptype4
	testre<-create_tabs3()
		if(numfiles==4){
                msmsProts11<-msmsPep11[[1]]$pepSeq
		msmsProts21<-msmsPep11[[2]]$pepSeq
		msmsProts31<-msmsPep11[[3]]$pepSeq
		msmsProts41<-msmsPep11[[4]]$pepSeq 
                #msmsProts41<-msmsPep1[[4]]$pep_seq                
		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
                InterSec13<-Reduce(intersect, list(msmsProts11,msmsProts31))
                InterSec14<-Reduce(intersect, list(msmsProts11,msmsProts41))
                InterSec23<-Reduce(intersect, list(msmsProts21,msmsProts31))
                InterSec24<-Reduce(intersect, list(msmsProts21,msmsProts41))
                InterSec34<-Reduce(intersect, list(msmsProts31,msmsProts41))
                InterSec123<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31))
                InterSec124<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts41))
                InterSec134<-Reduce(intersect, list(msmsProts11,msmsProts31,msmsProts41))
                InterSec234<-Reduce(intersect, list(msmsProts21,msmsProts31,msmsProts41))
                InterSec1234<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31,msmsProts41))

		venn.plot <- draw.quad.venn(
                        area1 = length(msmsPep11[[1]]$pepSeq),
			area2 = length(msmsPep11[[2]]$pepSeq),
                        area3 = length(msmsPep11[[3]]$pepSeq),
			area4 = length(msmsPep11[[4]]$pepSeq),
                        n12 = length(InterSec12),
                        n13 = length(InterSec13),
                        n14 = length(InterSec14),
                        n23 = length(InterSec23),
                        n24 = length(InterSec24),
                        n34 = length(InterSec34),
                        n123 = length(InterSec123),
                        n124 = length(InterSec124),
                        n134 = length(InterSec134),
                        n234 = length(InterSec234),
                        n1234 = length(InterSec1234),
			category = c("A", "B", "C", "D"),
                        fill = c("orange", "green", "darkorchid1", "cornflowerblue"),
                        lty = "blank",
                        col = c("orange", "cornflowerblue", "darkorchid1", "green"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen", "darkorchid1", "cornflowerblue")
                )
                mainT = "Overlaps of ADP-ribosylated peptides"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3," \nD = ",expType4,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	if(numfiles==3){
                msmsProts11<-msmsPep11[[1]]$pep_seq
		msmsProts21<-msmsPep11[[2]]$pep_seq
		msmsProts31<-msmsPep11[[3]]$pep_seq 
		

		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
                InterSec13<-Reduce(intersect, list(msmsProts11,msmsProts31))
                InterSec23<-Reduce(intersect, list(msmsProts21,msmsProts31))
                InterSec123<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31))

		venn.plot <- draw.triple.venn(
                        area1 = length(msmsPep11[[1]]$pep_seq),
			area2 = length(msmsPep11[[2]]$pep_seq),
                        area3 = length(msmsPep11[[3]]$pep_seq),
                        n12 = length(InterSec12),
                        n13 = length(InterSec13),
                        n23 = length(InterSec23),
                        n123 = length(InterSec123),
			category = c("A", "B", "C"),
                        fill = c("orange", "green", "darkorchid1"),
                        lty = "blank",
                        col = c("orange", "darkorchid1", "green"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen", "darkorchid1")
                )
                mainT = "Overlaps of ADP-ribosylated peptides"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	if(numfiles==2){
                msmsProts11<-msmsPep11[[1]]$pep_seq
		msmsProts21<-msmsPep11[[2]]$pep_seq
		

		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))

		venn.plot <- draw.pairwise.venn(
                        area1 = length(msmsPep11[[1]]$pep_seq),
			area2 = length(msmsPep11[[2]]$pep_seq),
                        cross.area = length(InterSec12),
			category = c("A", "B"),
                        fill = c("orange", "green"),
                        lty = "blank",
                        col = c("orange", "green"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen")
                )
                mainT = "Overlaps of ADP-ribosylated peptides"
		bottomT = paste("A = ",expType1," \nB = ",expType2,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	if(numfiles==1){
                msmsProts11<-msmsPep11[[1]]$pep_seq
                msmsProts21<-msmsPep11[[2]]$pep_seq


                InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))

                venn.plot <- draw.pairwise.venn(
                        area1 = length(msmsPep11[[1]]$pep_seq),
                        area2 = 0,
                        cross.area = 0,
                        category = c("A", "B"),
                        fill = c("orange", "green"),
                        lty = "blank",
                        col = c("orange", "green"),
                        cex = 1.5,
                        cat.cex = 2,
                        main = "Genes having increase in abundance",
                        main.cex = 2,
                        main.col = "black",
                        ind=FALSE,
                cat.col = c("darkorange", "darkgreen")
                )
                mainT = "Overlaps of ADP-ribosylated peptides"
                bottomT = paste("A = ",expType1,sep="")
                require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
        }
        })
output$heatmapProt <- renderPlot({
		expType1 <- input$exptype1
        	expType2 <- input$exptype2
        	expType3 <- input$exptype3
        	expType4 <- input$exptype4
		testre<-create_tabs3()
                if(numfiles==4){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab1 <- merge(x=voltab2, y=msmsHeat2[[4]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab1[,-1]
                row.names(voltab) <- voltab1[,1]

                colnames(voltab)<-c(expType1, expType2, expType3, expType4)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
		if(numfiles==3){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab2[,-1]
                row.names(voltab) <- voltab2[,1]

                colnames(voltab)<-c(expType1, expType2, expType3)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
		if(numfiles==2){
		voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab3[,-1]
                row.names(voltab) <- voltab3[,1]

                colnames(voltab)<-c(expType1, expType2)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
                my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                col_breaks = c(seq(0,0.5,length=100),
                seq(0.6,0.9,length=100),
                seq(1,2,length=100))
        	heatmap.2(voltab.mat, cellnote = voltab.mat, main = "No of psms detected per protein", notecol="black", density.info="none", trace="none", margins =c(12,9), col=my_palette, breaks=col_breaks, dendrogram="row", Colv="NA")
  })
	output$heatmapProtr <- renderPlotly({
		expType1 <- input$exptype1
        	expType2 <- input$exptype2
        	expType3 <- input$exptype3
        	expType4 <- input$exptype4
		testre<-create_tabs3()
                if(numfiles==4){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab1 <- merge(x=voltab2, y=msmsHeat2[[4]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab1[,-1]
                rownames(voltab) <- as.character(voltab1[,1])

                colnames(voltab)<-c(expType1, expType2, expType3, expType4)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
		if(numfiles==3){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab2[,-1]
                rownames(voltab) <- as.character(voltab2[,1])

                colnames(voltab)<-c(expType1, expType2, expType3)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
		if(numfiles==2){
		voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab3[,-1]
                rownames(voltab) <- as.character(voltab3[,1])

                colnames(voltab)<-c(expType1, expType2)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
                my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                col_breaks = c(seq(0,0.5,length=100),
                seq(0.6,0.9,length=100),
                seq(1,2,length=100))
        	#heatmap.2(voltab.mat, cellnote = voltab.mat, main = "No of psms detected per protein", notecol="black", density.info="none", trace="none", margins =c(12,9), col=my_palette, breaks=col_breaks, dendrogram="row", Colv="NA")
		#rwb <- colorRampPalette(colors = c("red", "yellow", "green"))
		heatmaply(voltab.mat, xlab = "Number of peptides", ylab = "Proteins", labRow= rownames(voltab.mat), labCol=colnames(voltab.mat),
                main = "Difference in ADPr peptides detected among samples",
                margins = c(60,100,40,20), scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "red", high = "green", midpoint = 2, limits = c(0, 30)))
		#colors = my_palette, limits = c(0, 10))        
  })

	output$psmBar1 <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs3()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfiles==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		psmCounts = c(length(msmsFiles5[[1]]$pep_seq), length(msmsFiles5[[2]]$pep_seq),length(msmsFiles5[[3]]$pep_seq),length(msmsFiles5[[4]]$pep_seq))
		)				
		ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
    		geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
    		scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		}
		else if(numfiles==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq), length(msmsFiles5[[2]]$pep_seq),length(msmsFiles5[[3]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq), length(msmsFiles5[[2]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
                scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
	})

	 output$psmBar2 <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs3()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfiles==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq), length(msmsFiles3[[2]]$pep_seq),length(msmsFiles3[[3]]$pep_seq),length(msmsFiles3[[4]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
                scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
                else if(numfiles==3){
		dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq), length(msmsFiles3[[2]]$pep_seq),length(msmsFiles3[[3]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
                else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq), length(msmsFiles3[[2]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
		scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
                else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
        })


	output$pepBar1 <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs3()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfiles==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		pepCounts = c(length(msmsPep11[[1]]$pepSeq), length(msmsPep11[[2]]$pepSeq),length(msmsPep11[[3]]$pepSeq),length(msmsPep11[[4]]$pepSeq))
		)				
		ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
    		geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths ignored)") +
    		scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		}
		else if(numfiles==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		pepCounts = c(length(msmsPep11[[1]]$pep_seq), length(msmsPep11[[2]]$pep_seq),length(msmsPep11[[3]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths ignored)") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		pepCounts = c(length(msmsPep11[[1]]$pep_seq), length(msmsPep11[[2]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths ignored)") +
                scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		pepCounts = c(length(msmsPep11[[1]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths ignored)") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
	})

	output$pepBar2 <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs3()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfiles==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq), length(msmsPep12[[2]]$pepSeq),length(msmsPep12[[3]]$pepSeq),length(msmsPep12[[4]]$pepSeq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths considered)") +
                scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
                else if(numfiles==3){
                dat1 <- data.frame(
		 sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                pepCounts = c(length(msmsPep12[[1]]$pep_seq), length(msmsPep12[[2]]$pep_seq),length(msmsPep12[[3]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths considered)") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
                else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                pepCounts = c(length(msmsPep12[[1]]$pep_seq), length(msmsPep12[[2]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths considered)") +
                scale_fill_manual(values=c("white", "grey68"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
                else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                pepCounts = c(length(msmsPep12[[1]]$pep_seq))
                )
                ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(lengths considered)") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
        })

	output$protBar <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4	
		testre<-create_tabs3()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfiles==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		protCounts = c(length(msmsProts1[[1]]$protName), length(msmsProts1[[2]]$protName),length(msmsProts1[[3]]$protName), length(msmsProts1[[4]]$protName))
		)				
		ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
    		geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
    		scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		}
		else if(numfiles==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		protCounts = c(length(msmsProts1[[1]]$protName), length(msmsProts1[[2]]$protName),length(msmsProts1[[3]]$protName))
                )
                ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		protCounts = c(length(msmsProts1[[1]]$protName), length(msmsProts1[[2]]$protName))
                )
                ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
                scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		protCounts = c(length(msmsProts1[[1]]$protName))
                )
                ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
	})
	output$psmBar1r <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs3()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfiles==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		psmCounts = c(length(msmsFiles5[[1]]$pep_seq), length(msmsFiles5[[2]]$pep_seq),length(msmsFiles5[[3]]$pep_seq),length(msmsFiles5[[4]]$pep_seq))
		)				
		p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
    		geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
    		scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
		}
		else if(numfiles==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq), length(msmsFiles5[[2]]$pep_seq),length(msmsFiles5[[3]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
		else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq), length(msmsFiles5[[2]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
                scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
		else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique PSMs") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
	})

	output$psmBar2r <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs3()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfiles==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq), length(msmsFiles3[[2]]$pep_seq),length(msmsFiles3[[3]]$pep_seq),length(msmsFiles3[[4]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
                scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
                else if(numfiles==3){
		dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq), length(msmsFiles3[[2]]$pep_seq),length(msmsFiles3[[3]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
                else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq), length(msmsFiles3[[2]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
		scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
                else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=psmCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Total PSMs") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
                }
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        })
	output$pepBar1r <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs3()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfiles==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		pepCounts = c(length(msmsPep11[[1]]$pepSeq), length(msmsPep11[[2]]$pepSeq),length(msmsPep11[[3]]$pepSeq),length(msmsPep11[[4]]$pepSeq))
		)				
		p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
    		geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length ignored)") +
    		scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
		}
		else if(numfiles==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		pepCounts = c(length(msmsPep11[[1]]$pep_seq), length(msmsPep11[[2]]$pep_seq),length(msmsPep11[[3]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length ignored)") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
		else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		pepCounts = c(length(msmsPep11[[1]]$pep_seq), length(msmsPep11[[2]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length ignored)") +
                scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
		else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		pepCounts = c(length(msmsPep11[[1]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length ignored)") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
	})

	output$pepBar2r <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs3()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfiles==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq), length(msmsPep12[[2]]$pepSeq),length(msmsPep12[[3]]$pepSeq),length(msmsPep12[[4]]$pepSeq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length considered)") +
                scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
                else if(numfiles==3){
                dat1 <- data.frame(
		 sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                pepCounts = c(length(msmsPep12[[1]]$pep_seq), length(msmsPep12[[2]]$pep_seq),length(msmsPep12[[3]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length considered)") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
                else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                pepCounts = c(length(msmsPep12[[1]]$pep_seq), length(msmsPep12[[2]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length considered)") +
                scale_fill_manual(values=c("white", "grey68"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
                else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                pepCounts = c(length(msmsPep12[[1]]$pep_seq))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=pepCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique peptides(Length considered)") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
        })

	output$protBarr <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs3()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfiles==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		protCounts = c(length(msmsProts3[[1]]$protName), length(msmsProts3[[2]]$protName),length(msmsProts3[[3]]$protName), length(msmsProts3[[4]]$protName))
		)				
		p <- ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
    		geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
    		scale_fill_manual(values=c("white", "grey68", "grey48", "black"))+
		theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
		}
		else if(numfiles==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		protCounts = c(length(msmsProts1[[1]]$protName), length(msmsProts1[[2]]$protName),length(msmsProts1[[3]]$protName))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
                scale_fill_manual(values=c("white", "grey68", "grey48"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
		else if(numfiles==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		protCounts = c(length(msmsProts1[[1]]$protName), length(msmsProts1[[2]]$protName))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
                scale_fill_manual(values=c("white", "grey68"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
		else if(numfiles==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		protCounts = c(length(msmsProts1[[1]]$protName))
                )
                p <- ggplot(data=dat1, aes(x=sample, y=protCounts, fill=sam)) +
                geom_bar(stat="identity", position=position_dodge(), colour=c("black")) +
		labs(title = "Unique proteins") +
                scale_fill_manual(values=c("white"))+
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
                }
	})

	output$goPlot1BP <- renderPlotly({
		expType1 <- input$exptype1
                testre<-create_tabs3()
		globalProt.1 <- globalProt[,1]
		organism <- input$organism
		if(numfiles>=1){
		localProt <- msmsHeat2[[1]]$protName[-1]
		goterm = "BP"
		goEnDataBP <- goEnrich(globalProt.1, localProt, organism, goterm)
		p <- ggplot(data = goEnDataBP, aes(x = Description,y=Count,
              	fill=p.adjust)) +
       		geom_bar(stat = 'identity') +
        	coord_flip() +
        	scale_x_discrete(limits = goEnDataBP$Description)+
		labs(title = paste(expType1," - BP enrichemnt", sep="")) +
        	ylab('Number of genes') +
        	xlab('')+
        	theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
        	scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
        	theme(panel.background = element_blank())
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
	}
	})
	output$goPlot1MF <- renderPlotly({
                expType1 <- input$exptype1
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles>=1){
                localProt <- msmsHeat2[[1]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[1]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[1]])
		goterm = "MF"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataMF <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataMF, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataMF$Description)+
		labs(title = paste(expType1," - MF enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
		p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	output$goPlot1CC <- renderPlotly({
                expType1 <- input$exptype1
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles>=1){
                localProt <- msmsHeat2[[1]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[1]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[1]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataCC, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataCC$Description)+
                labs(title = paste(expType1," - CC enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	output$goPlot2BP <- renderPlotly({
                expType2 <- input$exptype2
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles>=2){
                localProt <- msmsHeat2[[2]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataBP, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataBP$Description)+
                labs(title = paste(expType2," - BP enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	output$goPlot2MF <- renderPlotly({
                expType2 <- input$exptype2
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles>=2){
                localProt <- msmsHeat2[[2]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataMF, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataMF$Description)+
                labs(title = paste(expType2," - MF enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })

	output$goPlot3BP <- renderPlotly({
                expType3 <- input$exptype3
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles>=3){
                localProt <- msmsHeat2[[3]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataBP, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataBP$Description)+
                labs(title = paste(expType3," - BP enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	output$goPlot3MF <- renderPlotly({
                expType3 <- input$exptype3
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles>=1){
                localProt <- msmsHeat2[[3]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataMF, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataMF$Description)+
                labs(title = paste(expType3," - MF enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	output$goPlot4BP <- renderPlotly({
                expType4 <- input$exptype4
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles==4){
                localProt <- msmsHeat2[[4]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataBP, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataBP$Description)+
                labs(title = paste(expType4," - BP enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	output$goPlot4MF <- renderPlotly({
                expType4 <- input$exptype4
                testre<-create_tabs3()
                globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfiles==4){
                localProt <- msmsHeat2[[4]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(globalProt.1, localProt, organism, goterm)
                p <- ggplot(data = goEnDataMF, aes(x = Description,y=Count,
                fill=p.adjust)) +
                geom_bar(stat = 'identity') +
                coord_flip() +
                scale_x_discrete(limits = goEnDataMF$Description)+
                labs(title = paste(expType4," - MF enrichemnt", sep="")) +
                ylab('Number of genes') +
                xlab('')+
                theme(text = element_text(size=6), axis.text=element_text(size=6,face="bold"))+
                scale_fill_gradient2(low='grey8', mid='grey38', high='grey87', space='Lab', name = "p-value")+
                theme(panel.background = element_blank())
                p1 <- ggplotly(p)
                p1$elementId <- NULL
                print(p1)
        }
        })
	#output$text1 <- renderText({
	#	testre<-create_tabs2()
	#	table(names(modres))
     #   })

})




