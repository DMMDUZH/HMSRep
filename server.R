rm(list=ls()) 
if (!require("heatmaply")){
library(devtools)
install_github('talgalili/heatmaply')
}
if (!require("clusterProfiler")){
library(devtools)
install_github('GuangchuangYu/DOSE')
install_github('GuangchuangYu/enrichplot')
install_github('GuangchuangYu/clusterProfiler')
}
if (!require("org.Mm.eg.db")){
source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
biocLite("org.Hs.eg.db")
}

if (!require("lettercase")){
source("https://bioconductor.org/biocLite.R")
biocLite("lettercase")
}

if (!require("jsonlite")){
source("https://bioconductor.org/biocLite.R")
biocLite("jsonlite")
}

if (!require("stringr")){
source("https://bioconductor.org/biocLite.R")
biocLite("stringr")
}

if (!require("shinyBS")){
source("https://bioconductor.org/biocLite.R")
biocLite("shinyBS")
}

library(shinyBS)
library(stringr)
library(jsonlite)
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
	goEnrich <- function(localProt, organism, goterm, pVal){
		#function(localProt, geneList, organism){
		organism <- input$organism
                if(organism == "Mouse"){
                globalProt <<- read.table("mm_universal_list.txt", header = TRUE, sep="\t")
                }
                else if(organism == "Human"){
                globalProt <<- read.table("hs_universal_list.txt", header = TRUE, sep="\t")
                }
		globalProt.1 <- globalProt[,1]
		egy <- 0
		egx <- 0
		orgDB <- 0
		y.prot <- globalProt.1
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
		goEnDataGT <- enrichGoGT(gene, uniList, orgDB, goterm, pVal)
                return(goEnDataGT)
	}
	#This function computes GO enrichment
	enrichGoGT <- function(gene, uniList, orgDB, goterm, pVal){
                ego.gt <- enrichGO(gene    = gene,
                universe      = uniList,
                OrgDb         = orgDB,
                ont           = goterm,
                pAdjustMethod = "fdr",
                pvalueCutoff  = pVal,
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
EThcDFiles1 <- list()
HCDFiles1 <- list()
EThcDFiles2 <- list()
HCDFiles2 <- list()
EThcDFiles3 <- list()
HCDFiles3 <- list()
EThcDFiles4 <- list()
HCDFiles4 <- list()
EThcDFiles5 <- list()
HCDFiles5 <- list()
EThcDFiles6 <- list()
HCDFiles6 <- list()
EThcDFiles7 <- list()
HCDFiles7 <- list()
EThcDFiles8 <- list()
HCDFiles8 <- list()
EThcDFiles41 <- list()
HCDFiles9 <- list()
EThcDFiles51 <- list()
HCDFiles10 <- list()
EThcDFiles81 <- list()
HCDFiles11 <- list()
HCDFiles12 <- list()

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
msmsFiles9 <- list()
msmsFiles10 <- list()
msmsFiles11 <- list()
msmsFiles12 <- list()

stackBarH <- list()
HCDProts1 <- list()
HCDPep1 <- list()
HCDPep11 <- list()
HCDPep12 <- list()
HCDHeat1 <- list()
HCDHeat2 <- list()
HCDLocConf1 <- list()
HCDLocConf2 <- list()
HCDPepScore1 <- list()
HCDPepScore2 <- list()
HCDProts2 <- list()
HCDProts3 <- list()
HCDPep2 <- list()
HCDPep3 <- list()
HCDPep22 <- list()
HCDPep32 <- list()
HCDPsm1 <- list()
HCDPsm2 <- list()
HCDPsm3 <- list()

HCDEThcD1 <- list()
HCDEThcD2 <- list()
HCDEThcD3 <- list()
HCDEThcD4 <- list()
#globalProt <- vector()
	
	#Event reactive that reads the uploaded files and store the number of files loaded in a list
	#This also stores the number of the files
	csvdata1 <- eventReactive(input$choice,{
    		infile1 <- input$csvfile1
		req(infile1)
    		if (is.null(infile1)) {
      			return(NULL)
    		}
		else {
		numfilesE <<- nrow(infile1)
        		for (i in 1:numfilesE){
				#Reading csv files to the data objects in the list
				EThcDFiles1[[i]] <<- read.csv(input$csvfile1[[i, 'datapath']], header = input$header, sep = input$sep, quote = input$quote, fill = TRUE)
				#Below is a testing code that has no effect on the HMSRep script
				msmsName <- paste("msms1.", i, sep = "")
				assign(msmsName, read.csv(input$csvfile1[[i, 'datapath']], header = input$header, sep = input$sep, quote = input$quote, fill = TRUE))
			}
		}
	})
	csvdata2 <- eventReactive(input$choice,{
                infile2 <- input$csvfile2
                req(infile2)
                if (is.null(infile2)) {
                        return(NULL)
                }
                else {
                numfilesH <<- nrow(infile2)
                        for (i in 1:numfilesH){
                                #Reading csv files to the data objects in the list
                                HCDFiles1[[i]] <<- read.csv(input$csvfile2[[i, 'datapath']], header = input$header, sep = input$sep, quote = input$quote, fill = TRUE)
                                #Below is a testing code that has no effect on the HMSRep script
                                msmsName <- paste("msms2.", i, sep = "")
                                assign(msmsName, read.csv(input$csvfile2[[i, 'datapath']], header = input$header, sep = input$sep, quote = input$quote, fill = TRUE))
                        }
                }
        })
	#Reactive that is not neessary for the function of the script
	create_tabsE <- reactive({
		dataTab1 <- csvdata1()
                for (i in 1:numfilesE){
                        msmsName2 <- paste("msms2.", i, sep = "")
			if("X.query_number." %in% colnames(EThcDFiles1[[i]])){
				EThcDFiles2.t <- data.frame(protAcc=EThcDFiles1[[i]]$prot_acc, protDesc=EThcDFiles1[[i]]$prot_desc, pepSeq1=EThcDFiles1[[i]]$pep_seq,pepScore1=EThcDFiles1[[i]]$pep_score,pepExpect1=EThcDFiles1[[i]]$pep_expect,pepVarModConf=EThcDFiles1[[i]]$pep_var_mod_conf,pepVarMod=EThcDFiles1[[i]]$pep_var_mod,protAcc=EThcDFiles1[[i]]$prot_acc,pepSeq2=EThcDFiles1[[i]]$pep_seq.1,pepStart2=EThcDFiles1[[i]]$pep_res_before,pepEnd2=EThcDFiles1[[i]]$pep_res_after,protSeq=EThcDFiles1[[i]]$prot_seq,protScore=EThcDFiles1[[i]]$prot_score,protCover=EThcDFiles1[[i]]$prot_cover,pepScanTitle=EThcDFiles1[[i]]$pep_scan_title,queryNumber=EThcDFiles1[[i]]$X.query_number.,scanNumberRange=EThcDFiles1[[i]]$X.Scan.number.range.,pepMiss=EThcDFiles1[[i]]$pep_miss,pepVarModPos=EThcDFiles1[[i]]$pep_var_mod_pos, newScanTitle=paste(strsplit(paste(strsplit(as.character(EThcDFiles1[[i]]$X.StringTitle.),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\"))
			}
			else{
				EThcDFiles2.t <- data.frame(protAcc=EThcDFiles1[[i]]$prot_acc, protDesc=EThcDFiles1[[i]]$prot_desc, pepSeq1=EThcDFiles1[[i]]$pep_seq,pepScore1=EThcDFiles1[[i]]$pep_score,pepExpect1=EThcDFiles1[[i]]$pep_expect,pepVarModConf=EThcDFiles1[[i]]$pep_var_mod_conf,pepVarMod=EThcDFiles1[[i]]$pep_var_mod,protAcc=EThcDFiles1[[i]]$prot_acc,pepSeq2=EThcDFiles1[[i]]$pep_seq.1,pepStart2=EThcDFiles1[[i]]$pep_res_before,pepEnd2=EThcDFiles1[[i]]$pep_res_after,protSeq=EThcDFiles1[[i]]$prot_seq,protScore=EThcDFiles1[[i]]$prot_score,protCover=EThcDFiles1[[i]]$prot_cover,pepScanTitle=EThcDFiles1[[i]]$pep_scan_title,queryNumber=EThcDFiles1[[i]]$query_number,scanNumberRange=EThcDFiles1[[i]]$Scan.number.range,pepMiss=EThcDFiles1[[i]]$pep_miss,pepVarModPos=EThcDFiles1[[i]]$pep_var_mod_pos,newScanTitle=paste(strsplit(paste(strsplit(as.character(EThcDFiles1[[i]]$StringTitle),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\"))
			}
			EThcDFiles2[[i]] <<-  EThcDFiles2.t
		}
        })
	create_tabsH <- reactive({
                dataTab2 <- csvdata2()
                for (i in 1:numfilesH){
                        msmsName2 <- paste("msms2.", i, sep = "")
			if("X.query_number." %in% colnames(HCDFiles1[[i]])){
				HCDFiles2.t <- data.frame(protAcc=HCDFiles1[[i]]$prot_acc, protDesc=HCDFiles1[[i]]$prot_desc, pepSeq1=HCDFiles1[[i]]$pep_seq,pepScore1=HCDFiles1[[i]]$pep_score,pepExpect1=HCDFiles1[[i]]$pep_expect,pepVarModConf=HCDFiles1[[i]]$pep_var_mod_conf,pepVarMod=HCDFiles1[[i]]$pep_var_mod,protAcc=HCDFiles1[[i]]$prot_acc,pepSeq2=HCDFiles1[[i]]$pep_seq.1,pepStart2=HCDFiles1[[i]]$pep_res_before,pepEnd2=HCDFiles1[[i]]$pep_res_after,protSeq=HCDFiles1[[i]]$prot_seq,protScore=HCDFiles1[[i]]$prot_score,protCover=HCDFiles1[[i]]$prot_cover,pepScanTitle=HCDFiles1[[i]]$pep_scan_title,queryNumber=HCDFiles1[[i]]$X.query_number.,scanNumberRange=HCDFiles1[[i]]$X.Scan.number.range.,pepMiss=HCDFiles1[[i]]$pep_miss,pepVarModPos=HCDFiles1[[i]]$pep_var_mod_pos,newScanTitle=paste(strsplit(paste(strsplit(as.character(HCDFiles1[[i]]$X.StringTitle.),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\"))
			}
			else{
				HCDFiles2.t <- data.frame(protAcc=HCDFiles1[[i]]$prot_acc, protDesc=HCDFiles1[[i]]$prot_desc, pepSeq1=HCDFiles1[[i]]$pep_seq,pepScore1=HCDFiles1[[i]]$pep_score,pepExpect1=HCDFiles1[[i]]$pep_expect,pepVarModConf=HCDFiles1[[i]]$pep_var_mod_conf,pepVarMod=HCDFiles1[[i]]$pep_var_mod,protAcc=HCDFiles1[[i]]$prot_acc,pepSeq2=HCDFiles1[[i]]$pep_seq.1,pepStart2=HCDFiles1[[i]]$pep_res_before,pepEnd2=HCDFiles1[[i]]$pep_res_after,protSeq=HCDFiles1[[i]]$prot_seq,protScore=HCDFiles1[[i]]$prot_score,protCover=HCDFiles1[[i]]$prot_cover,pepScanTitle=HCDFiles1[[i]]$pep_scan_title,queryNumber=HCDFiles1[[i]]$query_number,scanNumberRange=HCDFiles1[[i]]$Scan.number.range,pepMiss=HCDFiles1[[i]]$pep_miss,pepVarModPos=HCDFiles1[[i]]$pep_var_mod_pos, newScanTitle=paste(strsplit(paste(strsplit(as.character(HCDFiles1[[i]]$StringTitle),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\"))
			}
                        HCDFiles2[[i]] <<-  HCDFiles2.t
                }
        })
	#Reactive that does most of the work in this script
	create_tabs3 <- reactive({
		organism <- input$organism
		#if(organism == "Mouse" && length(globalProt) == 0){
		#globalProt <<- read.table("mm_universal_list.txt", header = TRUE, sep="\t")
		#}
		#else if(organism == "Human" && length(globalProt) == 0){
		#globalProt <<- read.table("hs_universal_list.txt", header = TRUE, sep="\t")
		#}
		msmsFiles <- csvdata1()
                #searchVar1 <- input$searchtype1
		if(numfilesE >= 1){
		searchType <- read.table("ethcdSearch.txt", header=TRUE, sep="\t") 
                #searchOrderVar1 <- input$searchorder1
                for(i in 1:numfilesE){
                msms4.2 <- EThcDFiles1[[i]]
		#Initial filtering to accept the records having mascot score more than or equal to 20 and that there is a protein name allocated for the record
                msms4.2.1 <- msms4.2[msms4.2$pep_score >= 20 & msms4.2$prot_acc !="",]
		#Remove quotes from modification column if there is any
		msms4.2.1$pep_var_mod <- gsub("\"","",msms4.2.1$pep_var_mod)
		#Filter only the search type user selected
		for(v in 1:length(msms4.2.1$pep_var_mod)){
			for(j in  1:length(searchType[,1])){
				if(msms4.2.1$pep_var_mod[v] == searchType[j,1]){
					searchVar1 <- searchType[j,1]
					msms4.2.2 <- msms4.2.1[msms4.2.1$pep_var_mod == searchVar1,]
					v<-0
					break
				}
			}
			if(v==0){
				break
			}
		}
		for(v in 1:length(msms4.2.2$pep_var_mod_pos.1)){
			for(k in 1:5){
				if(grepl(k, msms4.2.2$pep_var_mod_pos.1[v])){
					searchOrderVar1 <- k
					v<-0
					break
				}
			}
			if(v==0){
                                break
                        }
		}
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
		ifelse(str_count(as.character(msmsFiles3[[i]]$prot_acc), "\\|")>1, protName <- sapply(strsplit(as.character(msmsFiles3[[i]]$prot_acc), '\\|'), '[', 2), protName <- sapply(strsplit(as.character(msmsFiles3[[i]]$prot_acc), '\\|'), '[', 1))
		msmsFiles3[[i]]$protName <<- sapply(strsplit(as.character(protName), '-'), '[', 1)
		#msmsFiles3[[i]]$protName <<- sapply(msmsFiles3[[i]]$protName, function(x) toString(tags$a(href=paste0("'","https://www.uniprot.org/uniprot/", x,"'"," target='_blank'"), x)))
		msmsFiles3[[i]]$protName <<- paste0("<a href='",  "https://www.uniprot.org/uniprot/", msmsFiles3[[i]]$protName, "' target='_blank'>",msmsFiles3[[i]]$protName,"</a>")
		#remove contaminats
		msmsFiles3[[i]] <<- msmsFiles3[[i]][msmsFiles3[[i]]$protName != "ZZ" & msmsFiles3[[i]]$protName != '"',]

		#Further process gene symbol column just created
		msmsFiles3[[i]]$protName1 <<- gsub(".*GN=","",msmsFiles3[[i]]$prot_desc)
		msmsFiles3[[i]]$protName1 <<- gsub("\\ PE=.*","",msmsFiles3[[i]]$protName1)
		msmsFiles3[[i]]$protName1 <<- gsub("\"", "", msmsFiles3[[i]]$protName1, fixed=TRUE)

		#Find exact position of the modification residue within each peptide		
		pos <- (unlist(gregexpr(pattern =searchOrderVar1,msmsFiles3[[i]]$pep_var_mod_pos.1))-2)

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

	msmsFiles4[[i]]$pep_var_mod_pos_num <<- (unlist(gregexpr(pattern =searchOrderVar1,msmsFiles4[[i]]$pep_var_mod_pos.1))-2)
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

	msmsFiles4[[i]]$fragMeth <<- "EThcD"

	msmsFiles5[[i]] <<- msmsFiles4[[i]][!duplicated(msmsFiles4[[i]][c(cx1,cx2)]),]

	msmsFiles7[[i]] <<- msmsFiles4[[i]][!duplicated(msmsFiles4[[i]][c(cx3,cx2)]),]

	msmsFiles9[[i]] <<- plyr::count(msmsFiles4[[i]], c("pep_seq.by", "pep_var_mod_res_site"))

	msmsFiles10[[i]] <<- plyr::count(msmsFiles4[[i]], c("pep_seq.2", "pep_var_mod_res_site"))

	#msmsFiles41[[i]] <<- with(msmsFiles4[[i]], ave(as.numeric(S, Cnty, FUN=function(x) length(unique(x)))) 	

	#msmsProts1[[i]] <<- msmsFiles5[[i]][!duplicated(msmsFiles5[[i]]$protName1),]

	msmsFiles6.t <- data.frame(protID=msmsFiles5[[i]]$protName, protName=msmsFiles5[[i]]$protName1, pepSeq=msmsFiles5[[i]]$pep_seq.2, pepSeq.upstream=msmsFiles5[[i]]$pep_seq.us, pepSeq.stem=msmsFiles5[[i]]$pep_seq.by,pepSeq.downstream=msmsFiles5[[i]]$pep_seq.ds, modRes=msmsFiles5[[i]]$pep_var_mod_res_num, protPos=msmsFiles5[[i]]$prot_var_mod_res_num, protSeq=msmsFiles5[[i]]$prot_seq, fragMeth=msmsFiles5[[i]]$fragMeth) 

	msmsFiles11.5 <- merge(msmsFiles5[[i]], msmsFiles9[[i]], by=c("pep_seq.by","pep_var_mod_res_site"))

	msmsFiles12.7 <- merge(msmsFiles7[[i]], msmsFiles10[[i]], by=c("pep_seq.2","pep_var_mod_res_site"))	

	msmsFiles11.9 <- data.frame(protID=msmsFiles11.5$protName, protName=msmsFiles11.5$protName1, pepSeq=msmsFiles11.5$pep_seq.2, pepSeq.upstream=msmsFiles11.5$pep_seq.us, pepSeq.stem=msmsFiles11.5$pep_seq.by, pepSeq.downstream=msmsFiles11.5$pep_seq.ds, modRes=msmsFiles11.5$pep_var_mod_res_num, Freq=msmsFiles11.5$freq, pepScore=msmsFiles11.5$pep_score, pepExpect=msmsFiles11.5$pep_expect, pepVarModConf=msmsFiles11.5$pep_var_mod_conf, protPos=msmsFiles11.5$prot_var_mod_res_num, protSeq=msmsFiles11.5$prot_seq, fragMeth=msmsFiles11.5$fragMeth, protSeq1=msmsFiles11.5$prot_seq, stringsAsFactors=F)

	msmsFiles11.9$protName <- paste0("<a href='https://string-db.org/api/image/network?identifiers=", msmsFiles11.9$protName, "&required_score=700' target='_blank'>", msmsFiles11.9$protName,"</a><br /> <a href='https://www.proteinatlas.org/search/",msmsFiles11.9$protName,"' target='_blank'>",msmsFiles11.9$protName,"</a>")

	#msmsFiles11.9$pepSeq.stem <- paste0("<a href='", "https://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=", msmsFiles11.9$protSeq1,"&eCutOff=10.0&matrix=BLOSUM62&outputFormat=HTML' target='_blank'>", msmsFiles11.9$pepSeq.stem,"</a>")

	msmsFiles12.10 <- data.frame(protID=msmsFiles12.7$protName, protName=msmsFiles12.7$protName1, pepSeq=msmsFiles12.7$pep_seq.2, pepSeq.upstream=msmsFiles12.7$pep_seq.us, pepSeq.stem=msmsFiles12.7$pep_seq.by, pepSeq.downstream=msmsFiles12.7$pep_seq.ds, modRes=msmsFiles12.7$pep_var_mod_res_num, Freq=msmsFiles12.7$freq, pepScore=msmsFiles12.7$pep_score, pepExpect=msmsFiles12.7$pep_expect, pepVarModConf=msmsFiles12.7$pep_var_mod_conf, protPos=msmsFiles12.7$prot_var_mod_res_num, protSeq=msmsFiles12.7$prot_seq, fragMeth=msmsFiles12.7$fragMeth, stringsAsFactors=F)

	msmsFiles12.10$protName <- paste0("<a href='https://string-db.org/api/image/network?identifiers=", msmsFiles12.10$protName, "&required_score=700' target='_blank'>", msmsFiles12.10$protName,"</a><br /> <a href='https://www.proteinatlas.org/search/",msmsFiles12.10$protName,"' target='_blank'>",msmsFiles12.10$protName,"</a>")

	msmsFiles11.9$protSeq <- paste0("<mark style=\"background-color: #339fff\">",msmsFiles11.9$protSeq,"</mark>")

	msmsFiles12.10$protSeq <- paste0("<mark style=\"background-color: #339fff\">",msmsFiles12.10$protSeq,"</mark>")

	msmsFiles11.9$protSeq <- str_replace(msmsFiles11.9$protSeq, msmsFiles11.9$pepSeq, paste0("</mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",msmsFiles11.9$pepSeq,"</span></mark><mark style=\"background-color: #339fff\">")) 

	msmsFiles12.10$protSeq <- str_replace(msmsFiles12.10$protSeq, msmsFiles12.10$pepSeq, paste0("</mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",msmsFiles12.10$pepSeq,"</span></mark><mark style=\"background-color: #339fff\">"))

	msmsFiles11.9$modPos <- as.numeric(str_extract(msmsFiles11.9$protPos, "[0-9]+"))

	msmsFiles11.9$modPosPep <- as.numeric(str_extract(msmsFiles11.9$modRes, "[0-9]+"))

	msmsFiles12.10$modPos <- as.numeric(str_extract(msmsFiles12.10$protPos, "[0-9]+"))

        msmsFiles12.10$modPosPep <- as.numeric(str_extract(msmsFiles12.10$modRes, "[0-9]+"))

	msmsFiles11.9$protSeq <- paste0(substring(msmsFiles11.9$protSeq, 1, msmsFiles11.9$modPos+114), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(msmsFiles11.9$protSeq, msmsFiles11.9$modPos+115, msmsFiles11.9$modPos+115),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(msmsFiles11.9$protSeq, msmsFiles11.9$modPos+116),sep="")

	msmsFiles11.9$pepSeq <-paste0("<mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",msmsFiles11.9$pepSeq,"</span></mark>")

	msmsFiles11.9$pepSeq <- paste0(substring(msmsFiles11.9$pepSeq, 1, msmsFiles11.9$modPosPep+67), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(msmsFiles11.9$pepSeq, msmsFiles11.9$modPosPep+68, msmsFiles11.9$modPosPep+68),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(msmsFiles11.9$pepSeq, msmsFiles11.9$modPosPep+69),sep="")

	msmsFiles12.10$protSeq <- paste0(substring(msmsFiles12.10$protSeq, 1, msmsFiles12.10$modPos+114), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(msmsFiles12.10$protSeq, msmsFiles12.10$modPos+115, msmsFiles12.10$modPos+115),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(msmsFiles12.10$protSeq, msmsFiles12.10$modPos+116),sep="")

	msmsFiles12.10$pepSeq <-paste0("<mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",msmsFiles12.10$pepSeq,"</span></mark>")

	msmsFiles12.10$pepSeq <- paste0(substring(msmsFiles12.10$pepSeq, 1, msmsFiles12.10$modPosPep+67), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(msmsFiles12.10$pepSeq, msmsFiles12.10$modPosPep+68, msmsFiles12.10$modPosPep+68),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(msmsFiles12.10$pepSeq, msmsFiles12.10$modPosPep+69),sep="")

	msmsFiles6[[i]] <<- msmsFiles6.t

	msmsFiles11[[i]] <<- msmsFiles11.9

	msmsFiles12[[i]] <<- msmsFiles12.10

	msmsProts1[[i]] <<- msmsFiles6.t[!duplicated(msmsFiles6.t$protName),]

	msmsFiles8.t <- data.frame(protID=msmsFiles7[[i]]$protName, protName=msmsFiles7[[i]]$protName1, pepSeq=msmsFiles7[[i]]$pep_seq.2, pepSeq.upstream=msmsFiles7[[i]]$pep_seq.us, pepSeq.stem=msmsFiles7[[i]]$pep_seq.by,pepSeq.downstream=msmsFiles7[[i]]$pep_seq.ds,modRes=msmsFiles7[[i]]$pep_var_mod_res_num, protPos=msmsFiles7[[i]]$prot_var_mod_res_num, protSeq=msmsFiles7[[i]]$prot_seq, fragMeth=msmsFiles7[[i]]$fragMeth)

	msmsFiles8[[i]] <<- msmsFiles8.t

                msmsPep11[[i]] <<- msmsFiles6[[i]]
		msmsPep12[[i]] <<- msmsFiles8[[i]]


		ngo <- length(msmsFiles5[[i]]$protName1)

		msmsHeat1[[i]] <<- data.frame(protName=character(ngo), sample=as.numeric(ngo))
                msmsHeat1[[i]]$protName  <<- as.character(msmsFiles5[[i]]$protName1)
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
		locConfTab$pep_var_mod_pos_num <- (unlist(gregexpr(pattern =searchOrderVar1,as.character(msmsFiles3[[i]]$pep_var_mod_pos.1)))-2)
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
	}
	})

	create_tabs4 <- reactive({
               	organism <- input$organism
                #if(organism == "Mouse" && length(globalProt) == 0){
                #globalProt <<- read.table("mm_universal_list.txt", header = TRUE, sep="\t")
                #}
                #else if(organism == "Human" && length(globalProt) == 0){
                #globalProt <<- read.table("hs_universal_list.txt", header = TRUE, sep="\t")
                #}
 		HCDFiles <- csvdata2()
                #searchVar2 <- input$searchtype2
		if(numfilesH >= 1){
                #searchOrderVar2 <- input$searchorder2
		searchType <- read.table("hcdSearch.txt", header=TRUE, sep="\t")
                for(i in 1:numfilesH){
                HCD4.2 <- HCDFiles1[[i]]
		#Initial filtering to accept the records having mascot score more than or equal to 20 and that there is a protein name allocated for the record
                HCD4.2.1 <- HCD4.2[HCD4.2$pep_score >= 20 & HCD4.2$prot_acc !="",]
		#Remove quotes from modification column if there is any
		HCD4.2.1$pep_var_mod <- gsub("\"","", HCD4.2.1$pep_var_mod)
		#Filter only the search type user selected
		#HCD4.2.2 <- HCD4.2.1[HCD4.2.1$pep_var_mod == searchVar2,]
		for(v in 1:length(HCD4.2.1$pep_var_mod)){
                        for(j in  1:length(searchType[,1])){
                                if(HCD4.2.1$pep_var_mod[v] == searchType[j,1]){
                                        searchVar2 <- searchType[j,1]
                                        HCD4.2.2 <- HCD4.2.1[HCD4.2.1$pep_var_mod == searchVar2,]
                                        v<-0
                                        break
                                }
                        }
                        if(v==0){
                                break
                        }
                }
                for(v in 1:length(HCD4.2.2$pep_var_mod_pos.1)){
                        for(k in 1:5){
                                if(grepl(k, HCD4.2.2$pep_var_mod_pos.1[v])){
                                        searchOrderVar2 <- k
                                        v<-0
                                        break
                                }
                        }
                        if(v==0){
                                break
                        }
                }
		HCD4.2.3 <- HCD4.2.2[FALSE,]
		#Remove corrupted records
		if("X.query_number." %in% colnames(HCD4.2.2)){
			HCD4.2.3 <- HCD4.2.2[HCD4.2.2$pep_expect <= 0.05 & HCD4.2.2$X.query_number. != "emPAI",]
		}
		else{
			HCD4.2.3 <- HCD4.2.2[HCD4.2.2$pep_expect <= 0.05 & HCD4.2.2$query_number != "emPAI",]
		}	
		#If modification confidence cell is empty write 99.99% in the cell
		HCD4.2.3$pep_var_mod_conf <- ifelse(HCD4.2.3$pep_var_mod_conf == "", "99.99%", as.character(HCD4.2.3$pep_var_mod_conf)) 
		#Create a numeric column out of the localization confidence column
		HCD4.2.3$pep_var_mod_conf.1 <- as.numeric(as.character(sapply(strsplit(as.character(HCD4.2.3$pep_var_mod_conf), "%"), `[`, 1)))

		HCD4.2.4 <- HCD4.2.3[FALSE,]
		HCD4.2.5 <- HCD4.2.3[FALSE,]

		#Sort and remove duplicate querries depending on the csv species
		if("X.query_number." %in% colnames(HCD4.2.3)){
			HCD4.2.4 <- HCD4.2.3[order(HCD4.2.3$X.query_number., -HCD4.2.3$pep_var_mod_conf.1),]
			HCD4.2.5 <- HCD4.2.4[!duplicated(HCD4.2.4$X.query_number.),]
		}
		else{
			HCD4.2.4 <- HCD4.2.3[order(HCD4.2.3$query_number, -HCD4.2.3$pep_var_mod_conf.1),]
                        HCD4.2.5 <- HCD4.2.4[!duplicated(HCD4.2.4$query_number),]
		}

		HCDFiles3[[i]] <<- HCD4.2.5

		#Extract gene symbol and store in a separate column
		ifelse(str_count(as.character(HCDFiles3[[i]]$prot_acc), "\\|")>1, protName <- sapply(strsplit(as.character(HCDFiles3[[i]]$prot_acc), '\\|'), '[', 2), protName <- sapply(strsplit(as.character(HCDFiles3[[i]]$prot_acc), '\\|'), '[', 1))
		HCDFiles3[[i]]$protName <<- sapply(strsplit(as.character(protName), '_'), '[', 1)
		HCDFiles3[[i]]$protName <<- paste0("<a href='",  "https://www.uniprot.org/uniprot/", HCDFiles3[[i]]$protName, "' target='_blank'>",HCDFiles3[[i]]$protName,"</a>")
		#remove contaminats
		HCDFiles3[[i]] <<- HCDFiles3[[i]][HCDFiles3[[i]]$protName != "ZZ" & HCDFiles3[[i]]$protName != '"',]

		#Further process gene symbol column just created
		HCDFiles3[[i]]$protName1 <<- gsub(".*GN=","",HCDFiles3[[i]]$prot_desc)
		HCDFiles3[[i]]$protName1 <<- gsub("\\ PE=.*","",HCDFiles3[[i]]$protName1)
		HCDFiles3[[i]]$protName1 <<- gsub("\"", "", HCDFiles3[[i]]$protName1, fixed=TRUE)

		#Find exact position of the modification residue within each peptide		
		pos <- (unlist(gregexpr(pattern =searchOrderVar2,HCDFiles3[[i]]$pep_var_mod_pos.1))-2)

		n <- length(HCDFiles3[[i]]$pep_seq.1)

		#Next we try to collapse duplicated peptides
		HCDFiles3[[i]]$pep_seq.2 <<- gsub("\\(|\\)", "",HCDFiles3[[i]]$pep_seq.1)
		HCDFiles3[[i]]$pep_seq.by <<- gsub("\\(|\\)", "",HCDFiles3[[i]]$pep_seq.1)
		HCDFiles3[[i]]$pep_seq.ds <<- gsub("\\(|\\)", "",HCDFiles3[[i]]$pep_seq.1)
		HCDFiles3[[i]]$pep_seq.us <<- gsub("\\(|\\)", "",HCDFiles3[[i]]$pep_seq.1)

		HCDFiles3[[i]] <<- HCDFiles3[[i]][order(nchar(HCDFiles3[[i]]$pep_seq.2), -HCDFiles3[[i]]$pep_var_mod_conf.1),]

		HCDFiles4[[i]] <<- HCDFiles3[[i]]

	for(k in 1:length(HCDFiles3[[i]]$pep_seq.2)){
        	pepSub1 <- paste0(".*(", HCDFiles3[[i]]$pep_seq.2[[k]], ").*")
        	pepSub2 <- paste0(HCDFiles3[[i]]$pep_seq.2[[k]], ".*")
        	pepSub3 <- paste0(".*", HCDFiles3[[i]]$pep_seq.2[[k]])
        	for(j in 1:length(HCDFiles3[[i]]$pep_seq.2)){
                	HCDFiles4[[i]]$pep_seq.by[[j]] <<- gsub(pepSub1,"\\1",as.character(HCDFiles4[[i]]$pep_seq.by[[j]]))
                	HCDFiles4[[i]]$pep_seq.us[[j]] <<- gsub(pepSub2,"",as.character(HCDFiles4[[i]]$pep_seq.us[[j]]))
                	HCDFiles4[[i]]$pep_seq.ds[[j]] <<- gsub(pepSub3,"",as.character(HCDFiles4[[i]]$pep_seq.ds[[j]]))
        	}
	}

	HCDFiles4[[i]]$pep_seq.2.len <<- str_length(HCDFiles4[[i]]$pep_seq.2)
	HCDFiles4[[i]]$pep_seq.by.len <<- str_length(HCDFiles4[[i]]$pep_seq.by)
	HCDFiles4[[i]]$pep_seq.us.len <<- str_length(HCDFiles4[[i]]$pep_seq.us)
	HCDFiles4[[i]]$pep_seq.ds.len <<- str_length(HCDFiles4[[i]]$pep_seq.ds)

	HCDFiles4[[i]]$pep_var_mod_pos_num <<- (unlist(gregexpr(pattern =searchOrderVar2,HCDFiles4[[i]]$pep_var_mod_pos.1))-2)
	HCDFiles4[[i]]$pep_var_mod_res <<- substring(HCDFiles4[[i]]$pep_seq.2, HCDFiles4[[i]]$pep_var_mod_pos_num, HCDFiles4[[i]]$pep_var_mod_pos_num)
	HCDFiles4[[i]]$pep_var_mod_res_num <<- paste0( HCDFiles4[[i]]$pep_var_mod_res, HCDFiles4[[i]]$pep_var_mod_pos_num)
	#HCDFiles4[[i]]$prot_var_mod_res_num <<- (unlist(regexpr(pattern=HCDFiles4[[i]]$pep_seq.2, HCDFiles4[[i]]$prot_seq, perl=TRUE)) + (HCDFiles4[[i]]$pep_var_mod_pos_num -1))

	HCDFiles4[[i]]$prot_var_mod_pos_num <<- mapply(function(x, y, z) (unlist(regexpr(pattern=x, y, perl=TRUE)) + (z - 1)), HCDFiles4[[i]]$pep_seq.2, HCDFiles4[[i]]$prot_seq, HCDFiles4[[i]]$pep_var_mod_pos_num)

	HCDFiles4[[i]]$prot_var_mod_res_num <<- paste0( HCDFiles4[[i]]$pep_var_mod_res, HCDFiles4[[i]]$prot_var_mod_pos_num)

	#HCDFiles4[[i]]$prot_var_mod_res_num <<- StrPos(HCDFiles4[[i]]$prot_seq, pattern=HCDFiles4[[i]]$pep_seq.2, pos = 1)

	HCDFiles4[[i]]$pep_var_mod_res_site <<- HCDFiles4[[i]]$pep_var_mod_pos_num

	for(k in 1:length(HCDFiles4[[i]]$pep_seq.us.len)){
		if(HCDFiles4[[i]]$pep_seq.us.len[[k]] > 0 & HCDFiles4[[i]]$pep_seq.us.len[[k]] >= HCDFiles4[[i]]$pep_var_mod_pos_num[[k]]){
        		HCDFiles4[[i]]$pep_var_mod_res_site[[k]] <<- HCDFiles4[[i]]$pep_var_mod_pos_num[[k]]
        	}
		else if(HCDFiles4[[i]]$pep_seq.us.len[[k]] > 0 & HCDFiles4[[i]]$pep_seq.us.len[[k]] < HCDFiles4[[i]]$pep_var_mod_pos_num[[k]]){
        		HCDFiles4[[i]]$pep_var_mod_res_site[[k]] <<- HCDFiles4[[i]]$pep_var_mod_pos_num[[k]] - HCDFiles4[[i]]$pep_seq.us.len[[k]]
        	}
		else {
        		HCDFiles4$pep_var_mod_res_site[[k]] <<- HCDFiles4$pep_var_mod_pos_num[[k]]
        	}
	}

	cx1 <- grep("pep_seq.by", colnames(HCDFiles4[[i]]))

	cx2 <- grep("pep_var_mod_res_site", colnames(HCDFiles4[[i]]))

	cx3 <- grep("pep_seq.2", colnames(HCDFiles4[[i]]))

	HCDFiles4[[i]]$fragMeth <<- "HCD" 

	HCDFiles5[[i]] <<- HCDFiles4[[i]][!duplicated(HCDFiles4[[i]][c(cx1,cx2)]),]

	HCDFiles7[[i]] <<- HCDFiles4[[i]][!duplicated(HCDFiles4[[i]][c(cx3,cx2)]),]	

	HCDFiles9[[i]] <<- plyr::count(HCDFiles4[[i]], c("pep_seq.by", "pep_var_mod_res_site"))

	HCDFiles10[[i]] <<- plyr::count(HCDFiles4[[i]], c("pep_seq.2", "pep_var_mod_res_site"))

	#HCDFiles41[[i]] <<- with(HCDFiles4[[i]], ave(as.numeric(S, Cnty, FUN=function(x) length(unique(x)))) 	

	#HCDProts1[[i]] <<- HCDFiles5[[i]][!duplicated(HCDFiles5[[i]]$protName1),]

	HCDFiles6.t <- data.frame(protID=HCDFiles5[[i]]$protName, protName=HCDFiles5[[i]]$protName1, pepSeq=HCDFiles5[[i]]$pep_seq.2, pepSeq.upstream=HCDFiles5[[i]]$pep_seq.us, pepSeq.stem=HCDFiles5[[i]]$pep_seq.by,pepSeq.downstream=HCDFiles5[[i]]$pep_seq.ds,modRes=HCDFiles5[[i]]$pep_var_mod_res_num, protPos=HCDFiles5[[i]]$prot_var_mod_res_num, protSeq=HCDFiles5[[i]]$prot_seq, fragMeth=HCDFiles5[[i]]$fragMeth) 

	HCDFiles11.5 <- merge(HCDFiles5[[i]], HCDFiles9[[i]], by=c("pep_seq.by","pep_var_mod_res_site"))

	HCDFiles12.7 <- merge(HCDFiles7[[i]], HCDFiles10[[i]], by=c("pep_seq.2","pep_var_mod_res_site"))	

	HCDFiles11.9 <- data.frame(protID=HCDFiles11.5$protName,protName=HCDFiles11.5$protName1, pepSeq=HCDFiles11.5$pep_seq.2, pepSeq.upstream=HCDFiles11.5$pep_seq.us, pepSeq.stem=HCDFiles11.5$pep_seq.by, pepSeq.downstream=HCDFiles11.5$pep_seq.ds, modRes=HCDFiles11.5$pep_var_mod_res_num, Freq=HCDFiles11.5$freq, pepScore=HCDFiles11.5$pep_score, pepExpect=HCDFiles11.5$pep_expect, pepVarModConf=HCDFiles11.5$pep_var_mod_conf, protPos=HCDFiles11.5$prot_var_mod_res_num, protSeq=HCDFiles11.5$prot_seq, fragMeth=HCDFiles11.5$fragMeth, protSeq1=HCDFiles11.5$prot_seq, stringsAsFactors=F)

	HCDFiles11.9$protName <- paste0("<a href='https://string-db.org/api/image/network?identifiers=", HCDFiles11.9$protName, "&required_score=700' target='_blank'>", HCDFiles11.9$protName,"</a><br /> <a href='https://www.proteinatlas.org/search/", HCDFiles11.9$protName,"' target='_blank'>", HCDFiles11.9$protName,"</a>")

	HCDFiles12.10 <- data.frame(protID=HCDFiles12.7$protName,protName=HCDFiles12.7$protName1, pepSeq=HCDFiles12.7$pep_seq.2, pepSeq.upstream=HCDFiles12.7$pep_seq.us, pepSeq.stem=HCDFiles12.7$pep_seq.by, pepSeq.downstream=HCDFiles12.7$pep_seq.ds, modRes=HCDFiles12.7$pep_var_mod_res_num, Freq=HCDFiles12.7$freq, pepScore=HCDFiles12.7$pep_score, pepExpect=HCDFiles12.7$pep_expect, pepVarModConf=HCDFiles12.7$pep_var_mod_conf, protPos=HCDFiles12.7$prot_var_mod_res_num, protSeq=HCDFiles12.7$prot_seq, fragMeth=HCDFiles12.7$fragMeth, stringsAsFactors=F)

	HCDFiles12.10$protName <- paste0("<a href='https://string-db.org/api/image/network?identifiers=", HCDFiles12.10$protName, "&required_score=700' target='_blank'>", HCDFiles12.10$protName,"</a><br /> <a href='https://www.proteinatlas.org/search/", HCDFiles12.10$protName,"' target='_blank'>", HCDFiles12.10$protName,"</a>")	

	HCDFiles11.9$protSeq <- paste0("<mark style=\"background-color: #339fff\">", HCDFiles11.9$protSeq,"</mark>")

	HCDFiles12.10$protSeq <- paste0("<mark style=\"background-color: #339fff\">", HCDFiles12.10$protSeq,"</mark>")

	HCDFiles11.9$protSeq <- str_replace(HCDFiles11.9$protSeq, HCDFiles11.9$pepSeq, paste0("</mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">", HCDFiles11.9$pepSeq,"</span></mark><mark style=\"background-color: #339fff\">")) 

	HCDFiles12.10$protSeq <- str_replace(HCDFiles12.10$protSeq, HCDFiles12.10$pepSeq, paste0("</mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",HCDFiles12.10$pepSeq,"</span></mark><mark style=\"background-color: #339fff\">"))

	HCDFiles11.9$modPos <- as.numeric(str_extract(HCDFiles11.9$protPos, "[0-9]+"))

	HCDFiles11.9$modPosPep <- as.numeric(str_extract(HCDFiles11.9$modRes, "[0-9]+"))

	HCDFiles12.10$modPos <- as.numeric(str_extract(HCDFiles12.10$protPos, "[0-9]+"))

        HCDFiles12.10$modPosPep <- as.numeric(str_extract(HCDFiles12.10$modRes, "[0-9]+"))

	HCDFiles11.9$protSeq <- paste0(substring(HCDFiles11.9$protSeq, 1, HCDFiles11.9$modPos+114), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(HCDFiles11.9$protSeq, HCDFiles11.9$modPos+115, HCDFiles11.9$modPos+115),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(HCDFiles11.9$protSeq, HCDFiles11.9$modPos+116),sep="")

	HCDFiles11.9$pepSeq <-paste0("<mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">", HCDFiles11.9$pepSeq,"</span></mark>")

	HCDFiles11.9$pepSeq <- paste(substring(HCDFiles11.9$pepSeq, 1, HCDFiles11.9$modPosPep+67), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(HCDFiles11.9$pepSeq, HCDFiles11.9$modPosPep+68, HCDFiles11.9$modPosPep+68),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(HCDFiles11.9$pepSeq, HCDFiles11.9$modPosPep+69),sep="")

	HCDFiles12.10$protSeq <- paste0(substring(HCDFiles12.10$protSeq, 1, HCDFiles12.10$modPos+114), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(HCDFiles12.10$protSeq, HCDFiles12.10$modPos+115, HCDFiles12.10$modPos+115),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(HCDFiles12.10$protSeq, HCDFiles12.10$modPos+116),sep="")

	HCDFiles12.10$pepSeq <-paste0("<mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">", HCDFiles12.10$pepSeq,"</span></mark>")

	HCDFiles12.10$pepSeq <- paste(substring(HCDFiles12.10$pepSeq, 1, HCDFiles12.10$modPosPep+67), "</span></mark><mark style=\"background-color: #fb0603\"><span style=\"color:#6a0403\">",substring(HCDFiles12.10$pepSeq, HCDFiles12.10$modPosPep+68, HCDFiles12.10$modPosPep+68),"</span></mark><mark style=\"background-color: #33ffb8\"><span style=\"color:#01583a\">",substring(HCDFiles12.10$pepSeq, HCDFiles12.10$modPosPep+69),sep="")

	HCDFiles6[[i]] <<- HCDFiles6.t

	HCDFiles11[[i]] <<- HCDFiles11.9

	HCDFiles12[[i]] <<- HCDFiles12.10

	HCDProts1[[i]] <<- HCDFiles6.t[!duplicated(HCDFiles6.t$protName),]

	HCDFiles8.t <- data.frame(protID=HCDFiles7[[i]]$protName,protName=HCDFiles7[[i]]$protName1, pepSeq=HCDFiles7[[i]]$pep_seq.2, pepSeq.upstream=HCDFiles7[[i]]$pep_seq.us, pepSeq.stem=HCDFiles7[[i]]$pep_seq.by,pepSeq.downstream=HCDFiles7[[i]]$pep_seq.ds,modRes=HCDFiles7[[i]]$pep_var_mod_res_num, protPos=HCDFiles7[[i]]$prot_var_mod_res_num, protSeq=HCDFiles7[[i]]$prot_seq, fragMeth=HCDFiles7[[i]]$fragMeth)

	HCDFiles8[[i]] <<- HCDFiles8.t

                HCDPep11[[i]] <<- HCDFiles6[[i]]
		HCDPep12[[i]] <<- HCDFiles8[[i]]

		ngo <- length(HCDFiles5[[i]]$protName1)

		HCDHeat1[[i]] <<- data.frame(protName=character(ngo), sample=as.numeric(ngo))
                HCDHeat1[[i]]$protName  <<- as.character(HCDFiles5[[i]]$protName1)
		HCDHeat1[[i]]$sample <<- rep(i,nrow(HCDHeat1[[i]]))
		HCDHeat2[[i]] <<- as.data.frame(table(HCDHeat1[[i]][,1]))
		colnames(HCDHeat2[[i]])[1] <<- "protName"
		colnames(HCDHeat2[[i]])[2] <<- as.character(HCDHeat1[[i]][1,2])
		locConfTabH <- data.frame(pep_seq=character(n),pep_score=as.numeric(n), pep_var_mod_conf=as.numeric(n), pep_var_mod_pos=as.numeric(n), pep_var_mod_pos_num=as.numeric(n), pep_var_length=as.numeric(n),pep_var_mod_res=character(n),stringsAsFactors=FALSE)
		pep.seq <- as.character(HCDFiles3[[i]]$pep_seq)
		x <- as.character(HCDFiles3[[i]]$pep_var_mod_conf)
		pep.var.mod.conf <- as.numeric(sapply(strsplit(as.character(x), "%"), `[`, 1))
		#HCDFiles3[[i]]$pep_var_mod_conf.n <- as.numeric(sapply(strsplit(as.character(HCDFiles3[[i]]$pep_var_mod_conf), "%"), `[`, 1))

		HCDLocConf1[[i]] <<- data.frame(pep_seq=character(n),pep_score=as.numeric(n), pep_var_mod_conf=as.numeric(n))
		HCDLocConf1[[i]]$pep_seq <<- HCDFiles3[[i]]$pep_seq
		HCDLocConf1[[i]]$pep_score <<- HCDFiles3[[i]]$pep_score
		HCDLocConf1[[i]]$pep_var_mod_conf <<- HCDFiles3[[i]]$pep_var_mod_conf.1
		HCDLocConf1[[i]]$group <<- cut(HCDLocConf1[[i]]$pep_var_mod_conf,seq(from = 0, to = 100, by = 10))
		
		HCDLocConf11 <- as.data.frame(table(HCDLocConf1[[i]][,4]))
		colnames(HCDLocConf11)[1] <- "locConf"
                colnames(HCDLocConf11)[2] <- "psmCounts"
		HCDLocConf2[[i]] <<- HCDLocConf11[HCDLocConf11$psmCounts >=1,]

		n1 <- length(HCDFiles5[[i]]$pep_seq.2)
		
		HCDProts2[[i]] <<- data.frame(protName=character(n), pep_seq=character(n),pep_score=as.numeric(n),pep_var_mod_conf=as.numeric(n))
               	HCDProts2[[i]]$pep_seq <<- HCDFiles3[[i]]$pep_seq
		HCDProts2[[i]]$protName <<- HCDFiles3[[i]]$protName1
                HCDProts2[[i]]$pep_score <<- HCDFiles3[[i]]$pep_score
                HCDProts2[[i]]$pep_var_mod_conf <<- HCDFiles3[[i]]$pep_var_mod_conf.1

                HCDProts21 <- as.data.frame(table(HCDProts2[[i]][,1]))
                colnames(HCDProts21)[1] <- "protName"
                colnames(HCDProts21)[2] <- "pepCount"
                HCDProts3[[i]] <<- HCDProts21

		n5 <- length(HCDFiles5[[i]]$pep_seq.1)

		HCDPep2[[i]] <<- data.frame(pep_seq=character(n5),protName=character(n5), pep_score=as.numeric(n5),pep_var_mod_conf=as.numeric(n5))
                HCDPep2[[i]]$pep_seq <<- HCDFiles5[[i]]$pep_seq
                HCDPep2[[i]]$protName <<- HCDFiles5[[i]]$protName1
                HCDPep2[[i]]$pep_score <<- HCDFiles5[[i]]$pep_score
                HCDPep2[[i]]$pep_var_mod_conf <<- HCDFiles5[[i]]$pep_var_mod_conf.1

                HCDPep21 <- as.data.frame(table(droplevels(HCDPep2[[i]][,1])))
                colnames(HCDPep21)[1] <- "pepSeq"
                colnames(HCDPep21)[2] <- "pepCount"
                HCDPep3[[i]] <<- HCDPep21

		n7 <- length(HCDFiles7[[i]]$pep_seq.1)

		HCDPep22[[i]] <<- data.frame(pep_seq=character(n7),protName=character(n7), pep_score=as.numeric(n7),pep_var_mod_conf=as.numeric(n7))
                HCDPep22[[i]]$pep_seq <<- HCDFiles7[[i]]$pep_seq
                HCDPep22[[i]]$protName <<- HCDFiles7[[i]]$protName1
                HCDPep22[[i]]$pep_score <<- HCDFiles7[[i]]$pep_score
                HCDPep22[[i]]$pep_var_mod_conf <<- HCDFiles7[[i]]$pep_var_mod_conf.1

                HCDPep212 <- as.data.frame(table(droplevels(HCDPep22[[i]][,1])))
                colnames(HCDPep212)[1] <- "pepSeq"
                colnames(HCDPep212)[2] <- "pepCount"
                HCDPep32[[i]] <<- HCDPep212

		HCDPepScore1[[i]] <<- data.frame(pep_seq=character(n),pep_score=as.numeric(n), pep_var_mod_conf=as.numeric(n))
		HCDPepScore1[[i]]$pep_seq <<- HCDFiles3[[i]]$pep_seq
                HCDPepScore1[[i]]$pep_score <<- HCDFiles3[[i]]$pep_score
                HCDPepScore1[[i]]$pep_var_mod_conf <<- HCDFiles3[[i]]$pep_var_mod_conf.1
		HCDPepScore1[[i]]$group <<- cut(HCDPepScore1[[i]]$pep_score, seq(from = 20, to = max(HCDPepScore1[[i]]$pep_score), by = 10))

                HCDPepScore11 <- as.data.frame(table(HCDPepScore1[[i]][,4]))
                colnames(HCDPepScore11)[1] <- "pepScore"
               	colnames(HCDPepScore11)[2] <- "psmCounts"
                HCDPepScore2[[i]] <<- HCDPepScore11[HCDPepScore11$psmCounts >=1,]

		

		locConfTabH$pep_seq <- HCDFiles3[[i]]$pep_seq.2
		locConfTabH$pep_score <- HCDFiles3[[i]]$pep_score
		locConfTabH$pep_var_mod_conf <- HCDFiles3[[i]]$pep_var_mod_conf.1
		locConfTabH$pep_var_mod_pos <- HCDFiles3[[i]]$pep_var_mod_pos.1
		locConfTabH$pep_var_mod_pos_num <- (unlist(gregexpr(pattern =searchOrderVar2,as.character(HCDFiles3[[i]]$pep_var_mod_pos.1)))-2)
		locConfTabH$pep_var_length <- lapply(pep.seq, function(x) max(nchar(x)))
		locConfTabH$pep_var_mod_res <- substring(locConfTabH$pep_seq, locConfTabH$pep_var_mod_pos_num, locConfTabH$pep_var_mod_pos_num)
		locConfTabH$groups <- as.numeric(cut2(locConfTabH$pep_var_mod_conf, g=3))
		locConfTabH$groups[locConfTabH$pep_var_mod_conf < 60 ] <- 1
		locConfTabH$groups[locConfTabH$pep_var_mod_conf > 60 & locConfTabH$pep_var_mod_conf  < 95] <- 2
		locConfTabH$groups[locConfTabH$pep_var_mod_conf > 95 ] <- 3
		locConfTabH1 <- locConfTabH
		pepVarModConfH <- dplyr::rename(count(locConfTabH, pep_var_mod_res, groups), Freq = n)
		pepVarModConfH <- data.frame(pepVarModConfH)
		pepVarModConfH$groups[pepVarModConfH$groups == 1] <- "< 60%"
		pepVarModConfH$groups[pepVarModConfH$groups == 2] <- "> 60% and < 95%"
		pepVarModConfH$groups[pepVarModConfH$groups == 3] <- "> 95%"
		stackBarH[[i]] <<- data.frame(matrix(ncol = 3, nrow = length(pepVarModConfH$groups)))
		colnames(stackBarH[[i]]) <<- c("groups", "Freq", "pep_var_mod_res")
		stackBarH[[i]]$groups <<- pepVarModConfH$groups
		stackBarH[[i]]$Freq <<- pepVarModConfH$Freq
		stackBarH[[i]]$pep_var_mod_res <<- pepVarModConfH$pep_var_mod_res
		}
	}	
        })

	create_tabs5 <- reactive({
		testreE<-create_tabs3()
		testreH<-create_tabs4()
                #searchVar1 <- input$searchtype1
		#searchVar2 <- input$searchtype2
		if(numfilesE >= 1 && numfilesH >= 1){
			for(i in 1:numfilesE){
				for(j in 1:numfilesH){	
					if("X.StringTitle." %in% colnames(msmsFiles3[[i]]) & "X.StringTitle." %in% colnames(HCDFiles3[[j]])){
						EThcD.ti <- paste(strsplit(paste(strsplit(as.character(msmsFiles3[[i]]$X.StringTitle.[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
						HCD.ti <- paste(strsplit(paste(strsplit(as.character(HCDFiles3[[j]]$X.StringTitle.[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
						if(EThcD.ti == HCD.ti){
                        				#HCDEThcD1[[i]] <- merge(x=msmsFiles6[[i]], y=HCDFiles6[[i]], by="X.StringTitle.", all.x=TRUE, all.y=TRUE)
							HCDEThcD1[[i]] <<- rbind(msmsFiles11[[i]], HCDFiles11[[j]])
							#colnames(HCDEThcD1[[i]]) 
						}
					}
					else if(!("X.StringTitle." %in% colnames(msmsFiles3[[i]])) & "X.StringTitle." %in% colnames(HCDFiles3[[j]])){
						EThcD.ti <- paste(strsplit(paste(strsplit(as.character(msmsFiles3[[i]]$StringTitle[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
						HCD.ti <- paste(strsplit(paste(strsplit(as.character(HCDFiles3[[j]]$X.StringTitle.[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
                                                if(EThcD.ti == HCD.ti){
							#HCDEThcD1[[i]] <- merge(x=msmsFiles6[[i]], y=HCDFiles6[[i]], by="StringTitle", all.x=TRUE, all.y=TRUE)
							HCDEThcD1[[i]] <<- rbind(msmsFiles11[[i]], HCDFiles11[[j]])
						}
					}
					else if("X.StringTitle." %in% colnames(msmsFiles3[[i]]) & !("X.StringTitle." %in% colnames(HCDFiles3[[j]]))){
						EThcD.ti <- paste(strsplit(paste(strsplit(as.character(msmsFiles3[[i]]$X.StringTitle.[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
						HCD.ti <- paste(strsplit(paste(strsplit(as.character(HCDFiles3[[j]]$StringTitle[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
                                                if(EThcD.ti == HCD.ti){
							#HCDEThcD1[[i]] <- merge(x=msmsFiles6[[i]], y=HCDFiles6[[i]], by="StringTitle", all.x=TRUE, all.y=TRUE)
							HCDEThcD1[[i]] <<- rbind(msmsFiles11[[i]], HCDFiles11[[j]])
						}
					}
                			else{
						EThcD.ti <- paste(strsplit(paste(strsplit(as.character(msmsFiles3[[i]]$StringTitle[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
						HCD.ti <- paste(strsplit(paste(strsplit(as.character(HCDFiles3[[j]]$StringTitle[1]),"\\\\")[[1]][6:7],collapse = "\\\\"), "~")[[1]][1], collapse = "\\\\")
                                                if(EThcD.ti == HCD.ti){
							#HCDEThcD1[[i]] <- merge(x=msmsFiles6[[i]], y=HCDFiles6[[i]], by="StringTitle", all.x=TRUE, all.y=TRUE)
							HCDEThcD1[[i]] <<- rbind(msmsFiles11[[i]], HCDFiles11[[j]])
						}
					}
                		}
			}
		}
	})
	
	observe({
		 updateCheckboxGroupInput(session, "select1",
      			label = "Columns in msms_file1 to show",
      			choices = names(create_tabsE()),
      			selected = names(create_tabsE())
    			)
		})

	output$mytable1 <- renderDataTable({
		testre<-create_tabsE()
		columns <- names(EThcDFiles2[[1]])
      		if (!is.null(input$select1)) {
        		columns <- input$select1
      		}
   		EThcDFiles2[[1]][,columns,drop=FALSE] 
  	},caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)
	output$mytable2 <- renderDataTable({
		testre<-create_tabsE()
		if(numfilesE>=2){
                columns <- names(EThcDFiles2[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                EThcDFiles2[[2]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)
	output$mytable3 <- renderDataTable({
                testre<-create_tabsE()
                if(numfilesE>=3){
		columns <- names(EThcDFiles2[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                EThcDFiles2[[3]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)
	output$mytable4 <- renderDataTable({
                testre<-create_tabsE()
                if(numfilesE>=4){
		columns <- names(EThcDFiles2[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                EThcDFiles2[[4]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(pageLength = 10,
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)

	output$mytable11 <- renderDataTable({
		testre<-create_tabs3()
		if(numfilesE>=1){
                columns <- names(msmsProts3[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[1]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsProts3[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$mytable21 <- renderDataTable({
		testre<-create_tabs3()
                if(numfilesE>=2){
		columns <- names(msmsProts3[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[2]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsProts3[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable31 <- renderDataTable({
		testre<-create_tabs3()
                if(numfilesE>=3){
		columns <- names(msmsProts3[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[3]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsProts3[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytable41 <- renderDataTable({
		testre<-create_tabs3()
		if(numfilesE>=4){
                columns <- names(msmsProts3[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsProts3[[4]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsProts3[[4]]),
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
        options = list(pageLength = nrow(stackBar[[1]]),
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
        options = list(pageLength = nrow(stackBar[[1]]),
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
        options = list(pageLength = nrow(stackBar[[2]]),
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
        options = list(pageLength = nrow(stackBar[[3]]),
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
        options = list(pageLength = nrow(stackBar[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

#output$uniquePepTable1 <- DT::renderDataTable({
output$uniquePepTable1 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                #columns <- names(msmsFiles11[[1]])
		colnames(msmsFiles11[[1]]) <- str_wrap(colnames(msmsFiles11[[1]]),width = 10)
		columns <- names(msmsFiles11[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
		msmsFiles11[[1]][,columns,drop=FALSE] 
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
        escape = FALSE,
	extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles11[[1]]),
	#	rowCallback = JS(
        #              "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
        #              "var full_text = aData[2]",
        #              "$('td:eq(2)', nRow).attr('title', full_text);",
        #              "$('td:eq(2)', nRow).tooltip({","'delay': 0,","'track': true,","'fade': 250,","});",
	#	      "}"),
	#	escape = FALSE,
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)


output$HCDEThcDTab1 <- renderDataTable({
                testre<-create_tabs5()
                if(numfilesE>=1){
                columns <- names(HCDEThcD1[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDEThcD1[[1]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDEThcD1[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$HCDEThcDTab21 <- renderDataTable({
                testre<-create_tabs5()
                if(numfilesE>=2){
                columns <- names(HCDEThcD1[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDEThcD1[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDEThcD1[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$HCDEThcDTab31 <- renderDataTable({
                testre<-create_tabs5()
                if(numfilesE>=3){
                columns <- names(HCDEThcD1[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDEThcD1[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDEThcD1[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$HCDEThcDTab41 <- renderDataTable({
                testre<-create_tabs5()
                if(numfilesE==4){
                columns <- names(HCDEThcD1[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDEThcD1[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDEThcD1[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable11 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles12[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles12[[1]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
        escape = FALSE,
	extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles12[[1]]),
		#rowCallback=JS(
                #"function(row,data) {
                #data$protSeq=data$protSeq.replace('LAKYNQLLR','<span style=\"color:red\">LAKYNQLLR</span>');
                #$('td:eq(0)', row).html(data$protSeq);
                #}"
                #),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable2 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles11[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles11[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles11[[2]]),
		#rowCallback=JS(
                #"function(row,data) {
                #data$protSeq=data$protSeq.replace('LAKYNQLLR','<span style=\"color:red\">LAKYNQLLR</span>');
                #$('td:eq(0)', row).html(data$protSeq);
                #}"
                #),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable21 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles12[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles12[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
        escape = FALSE,
	extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles12[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable3 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles11[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles11[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles11[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable31 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles12[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles12[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles12[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable4 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles12[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles12[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles12[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTable41 <- renderDataTable({
                testre<-create_tabs3()
                if(numfilesE>=1){
                columns <- names(msmsFiles12[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                msmsFiles12[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(msmsFiles12[[4]]),
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
	if(numfilesE==4){
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
                mainT = "Overlaps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3," \nD = ",expType4,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	else if(numfilesE==3){
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
                mainT = "Overlaps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3, sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	else if(numfilesE==2){
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
                mainT = "Overlaps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	else if(numfilesE==1){
                msmsProts11<-msmsProts1[[1]]$protName
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
                mainT = "Overlaps of ADP-ribosylated proteins"
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
		if(numfilesE==4){
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
	else if(numfilesE==3){
                msmsProts11<-msmsPep11[[1]]$pepSeq
		msmsProts21<-msmsPep11[[2]]$pepSeq
		msmsProts31<-msmsPep11[[3]]$pepSeq 
		

		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))
                InterSec13<-Reduce(intersect, list(msmsProts11,msmsProts31))
                InterSec23<-Reduce(intersect, list(msmsProts21,msmsProts31))
                InterSec123<-Reduce(intersect, list(msmsProts11,msmsProts21,msmsProts31))

		venn.plot <- draw.triple.venn(
                        area1 = length(msmsPep11[[1]]$pepSeq),
			area2 = length(msmsPep11[[2]]$pepSeq),
                        area3 = length(msmsPep11[[3]]$pepSeq),
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
	else if(numfilesE==2){
                msmsProts11<-msmsPep11[[1]]$pepSeq
		msmsProts21<-msmsPep11[[2]]$pepSeq
		

		InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))

		venn.plot <- draw.pairwise.venn(
                        area1 = length(msmsPep11[[1]]$pepSeq),
			area2 = length(msmsPep11[[2]]$pepSeq),
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
	else if(numfilesE==1){
                msmsProts11<-msmsPep11[[1]]$pepSeq
                InterSec12<-Reduce(intersect, list(msmsProts11,msmsProts21))

                venn.plot <- draw.pairwise.venn(
                        area1 = length(msmsPep11[[1]]$pepSeq),
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
                if(numfilesE==4){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab1 <- merge(x=voltab2, y=msmsHeat2[[4]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab1[,-1]
                row.names(voltab) <- voltab1[,1]

                colnames(voltab)<-c(expType1, expType2, expType3, expType4)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
		if(numfilesE==3){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab2[,-1]
                row.names(voltab) <- voltab2[,1]

                colnames(voltab)<-c(expType1, expType2, expType3)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
		if(numfilesE==2){
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
                if(numfilesE==4){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab1 <- merge(x=voltab2, y=msmsHeat2[[4]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab1[,-1]
                rownames(voltab) <- as.character(voltab1[,1])

                colnames(voltab)<-c(expType1, expType2, expType3, expType4)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
		if(numfilesE==3){
                voltab3 <- merge(x=msmsHeat2[[1]], y=msmsHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=msmsHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab2[,-1]
                rownames(voltab) <- as.character(voltab2[,1])

                colnames(voltab)<-c(expType1, expType2, expType3)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
		if(numfilesE==2){
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
		if(numfilesE==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2), length(msmsFiles5[[2]]$pep_seq.2),length(msmsFiles5[[3]]$pep_seq.2),length(msmsFiles5[[4]]$pep_seq.2))
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
		else if(numfilesE==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2), length(msmsFiles5[[2]]$pep_seq.2),length(msmsFiles5[[3]]$pep_seq.2))
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
		else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2), length(msmsFiles5[[2]]$pep_seq.2))
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
		else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2))
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
                if(numfilesE==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2), length(msmsFiles3[[2]]$pep_seq.2),length(msmsFiles3[[3]]$pep_seq.2),length(msmsFiles3[[4]]$pep_seq.2))
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
                else if(numfilesE==3){
		dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2), length(msmsFiles3[[2]]$pep_seq.2),length(msmsFiles3[[3]]$pep_seq.2))
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
                else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2), length(msmsFiles3[[2]]$pep_seq.2))
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
                else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2))
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
		if(numfilesE==4){
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
		else if(numfilesE==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		pepCounts = c(length(msmsPep11[[1]]$pepSeq), length(msmsPep11[[2]]$pepSeq),length(msmsPep11[[3]]$pepSeq))
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
		else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		pepCounts = c(length(msmsPep11[[1]]$pepSeq), length(msmsPep11[[2]]$pepSeq))
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
		else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		pepCounts = c(length(msmsPep11[[1]]$pepSeq))
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
                if(numfilesE==4){
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
                else if(numfilesE==3){
                dat1 <- data.frame(
		 sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq), length(msmsPep12[[2]]$pepSeq),length(msmsPep12[[3]]$pepSeq))
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
                else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq), length(msmsPep12[[2]]$pepSeq))
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
                else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq))
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
		if(numfilesE==4){
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
		else if(numfilesE==3){
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
		else if(numfilesE==2){
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
		else if(numfilesE==1){
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
		if(numfilesE==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2), length(msmsFiles5[[2]]$pep_seq.2),length(msmsFiles5[[3]]$pep_seq.2),length(msmsFiles5[[4]]$pep_seq.2))
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
		else if(numfilesE==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2), length(msmsFiles5[[2]]$pep_seq.2),length(msmsFiles5[[3]]$pep_seq.2))
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
		else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2), length(msmsFiles5[[2]]$pep_seq.2))
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
		else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles5[[1]]$pep_seq.2))
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
                if(numfilesE==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2), length(msmsFiles3[[2]]$pep_seq.2),length(msmsFiles3[[3]]$pep_seq.2),length(msmsFiles3[[4]]$pep_seq.2))
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
                else if(numfilesE==3){
		dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2), length(msmsFiles3[[2]]$pep_seq.2),length(msmsFiles3[[3]]$pep_seq.2))
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
                else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2), length(msmsFiles3[[2]]$pep_seq.2))
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
                else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(msmsFiles3[[1]]$pep_seq.2))
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
		if(numfilesE==4){
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
		else if(numfilesE==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		pepCounts = c(length(msmsPep11[[1]]$pepSeq), length(msmsPep11[[2]]$pepSeq),length(msmsPep11[[3]]$pepSeq))
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
		else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		pepCounts = c(length(msmsPep11[[1]]$pepSeq), length(msmsPep11[[2]]$pepSeq))
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
		else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		pepCounts = c(length(msmsPep11[[1]]$pepSeq))
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
                if(numfilesE==4){
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
                else if(numfilesE==3){
                dat1 <- data.frame(
		 sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq), length(msmsPep12[[2]]$pepSeq),length(msmsPep12[[3]]$pepSeq))
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
                else if(numfilesE==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq), length(msmsPep12[[2]]$pepSeq))
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
                else if(numfilesE==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                pepCounts = c(length(msmsPep12[[1]]$pepSeq))
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
		if(numfilesE==4){
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
		else if(numfilesE==3){
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
		else if(numfilesE==2){
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
		else if(numfilesE==1){
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
                pValBP <- input$pval1
		testre<-create_tabs3()
		#globalProt.1 <- globalProt[,1]
		organism <- input$organism
		if(numfilesE>=1){
		localProt <- msmsHeat2[[1]]$protName[-1]
		goterm = "BP"
		#goEnDataBP <- goEnrich(globalProt.1, localProt, organism, goterm, pValBP)
		goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=1){
                localProt <- msmsHeat2[[1]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[1]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[1]])
		goterm = "MF"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=1){
                localProt <- msmsHeat2[[1]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[1]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[1]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
		pValBP <- input$pval1
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=2){
                localProt <- msmsHeat2[[2]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=2){
                localProt <- msmsHeat2[[2]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot2CC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=2){
                localProt <- msmsHeat2[[2]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[2]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[2]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
	output$goPlot3BP <- renderPlotly({
                expType3 <- input$exptype3
		pValBP <- input$pval1
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=3){
                localProt <- msmsHeat2[[3]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=3){
                localProt <- msmsHeat2[[3]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot3CC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE>=3){
                localProt <- msmsHeat2[[3]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[3]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[3]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
	output$goPlot4BP <- renderPlotly({
                expType4 <- input$exptype4
		pValBP <- input$pval1
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE==4){
                localProt <- msmsHeat2[[4]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE==4){
                localProt <- msmsHeat2[[4]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot4CC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesE==4){
                localProt <- msmsHeat2[[4]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[4]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[4]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
	#output$text1 <- renderText({
	#	testre<-create_tabs2()
	#	table(names(modres))
     #   })
	output$mytableH1 <- renderDataTable({
		testre<-create_tabsH()
		columns <- names(HCDFiles2[[1]])
      		if (!is.null(input$select1)) {
        		columns <- input$select1
      		}
   		HCDFiles2[[1]][,columns,drop=FALSE] 
  	},caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)
	output$mytableH2 <- renderDataTable({
		testre<-create_tabsH()
		if(numfilesH>=2){
                columns <- names(HCDFiles2[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles2[[2]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)
	output$mytableH3 <- renderDataTable({
                testre<-create_tabsH()
                if(numfilesH>=3){
		columns <- names(HCDFiles2[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles2[[3]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)
	output$mytableH4 <- renderDataTable({
                testre<-create_tabsH()
                if(numfilesH>=4){
		columns <- names(HCDFiles2[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles2[[4]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having site localization information"),
        extensions = 'Buttons',
        options = list(
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
	)

	output$mytableH11 <- renderDataTable({
		testre<-create_tabs4()
		if(numfilesH>=1){
                columns <- names(HCDProts3[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDProts3[[1]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDProts3[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$mytableH21 <- renderDataTable({
		testre<-create_tabs4()
                if(numfilesH>=2){
		columns <- names(HCDProts3[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDProts3[[2]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDProts3[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytableH31 <- renderDataTable({
		testre<-create_tabs4()
                if(numfilesH>=3){
		columns <- names(HCDProts3[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDProts3[[3]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDProts3[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytableH41 <- renderDataTable({
		testre<-create_tabs4()
		if(numfilesH>=4){
                columns <- names(HCDProts3[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDProts3[[4]][,columns,drop=FALSE]
		}
        },caption=paste("MSMS table having peptide info of single ADPr modification"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDProts3[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)
	output$mytableH12 <- renderDataTable({
		testre<-create_tabs4()
                columns <- names(stackBarH[[1]])
                stackBarH[[1]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(stackBarH[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$mytableH13 <- renderDataTable({
                testre<-create_tabs4()
                columns <- names(stackBarH[[1]])
                stackBarH[[1]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(stackBarH[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytableH22 <- renderDataTable({
		testre<-create_tabs4()
                columns <- names(stackBarH[[2]])
                stackBarH[[2]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(stackBarH[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytableH32 <- renderDataTable({
		testre<-create_tabs4()
                columns <- names(stackBarH[[3]])
                stackBarH[[3]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(stackBarH[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$mytableH42 <- renderDataTable({
		testre<-create_tabs4()
                columns <- names(stackBarH[[4]])
                stackBarH[[4]][,columns,drop=FALSE]
        },caption=paste("Distribution of site localization based on the confidence score"),
        extensions = 'Buttons',
        options = list(pageLength = nrow(stackBarH[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH1 <- DT::renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles11[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
		HCDFiles11[[1]][,columns,drop=FALSE]
                #DT::datatable(HCDFiles6[[1]][,columns,drop=FALSE],options = list(searchHighlight = TRUE, search = list(search = 'YQKSTELLIR')))
		#formatStyle('Sepal.Width',color = styleInterval(c(3.4, 3.8), c('white', 'blue', 'red')),backgroundColor = styleInterval(3.4, c('gray', 'yellow')))
		#,options = list(searchHighlight = TRUE, search = list(search = HCDFiles6[[1]]$pep_seq.2)))
		#[,columns,drop=FALSE]
		#,options = list(searchHighlight = TRUE, search = list(search = HCDFiles6[[1]]$pep_seq.2)))
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
        escape = FALSE,
	extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles11[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH11 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles12[[1]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles12[[1]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles12[[1]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH2 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles11[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles11[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles11[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH21 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles12[[2]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles12[[2]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles12[[2]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH3 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles11[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles11[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles11[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH31 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles12[[3]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles12[[3]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles12[[3]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH4 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles11[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles11[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths ignored)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles11[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

output$uniquePepTableH41 <- renderDataTable({
                testre<-create_tabs4()
                if(numfilesH>=1){
                columns <- names(HCDFiles12[[4]])
                if (!is.null(input$select1)) {
                        columns <- input$select1
                }
                HCDFiles12[[4]][,columns,drop=FALSE]
                }
        },caption=paste("Unique peptides with asingle ADPr modification (lengths considered)"),
	escape = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = nrow(HCDFiles12[[4]]),
                "dom" = 'T<"clear">lBfrtip',
                buttons = list('copy', list(extend='csv',filename="diff_reg_genes"), list(extend='excel',filename="diff_reg_genes"), list(extend='pdf',filename="diff_reg_genes"), 'print')
                )
)

	output$freqdistrH1 <- renderPlotly({
		expType1 <- input$exptype1
		testre<-create_tabs4()
		stackBart <- stackBarH[[1]]
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

		output$psmLocConfH1 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDLocConf2[[1]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))

	})

		output$psmPepScoreH1 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDPepScore2[[1]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmPepScoreH2 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDPepScore2[[2]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmPepScoreH3 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDPepScore2[[3]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmPepScoreH4 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDPepScore2[[4]], aes(x=pepScore, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
                geom_smooth()+
		labs(title = "Distribution of PSMS across mascot score") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
			
		output$psmLocConfH2 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDLocConf2[[2]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
		geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmLocConfH3 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDLocConf2[[3]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
		geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })
		output$psmLocConfH4 <- renderPlot({
                testre<-create_tabs4()
                ggplot(data=HCDLocConf2[[4]], aes(x=locConf, y=psmCounts, group=1)) +
                geom_line()+
                geom_point()+
		geom_smooth()+
		labs(title = "Distribution of PSMS across localization confidence") +
                theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = 'dotted',
                                colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1))
        })

		output$freqdistrH12 <- renderPlot({
		expType1 <- input$exptype1
                testre<-create_tabs4()
                stackBart <- stackBarH[[1]]
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
output$freqdistrH2 <- renderPlotly({
		expType2 <- input$exptype2
		testre<-create_tabs4()
		stackBart <- stackBarH[[2]]
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

		output$freqdistrH22 <- renderPlot({
		expType2 <- input$exptype2
                testre<-create_tabs4()
                stackBart <- stackBarH[[2]]
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
output$freqdistrH3 <- renderPlotly({
		expType3 <- input$exptype3
		testre<-create_tabs4()
		stackBart <- stackBarH[[3]]
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

		output$freqdistrH32 <- renderPlot({
		expType3 <- input$exptype3
                testre<-create_tabs4()
                stackBart <- stackBarH[[3]]
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
	output$freqdistrH4 <- renderPlotly({
		expType4 <- input$exptype4
		testre<-create_tabs4()
		stackBart <- stackBarH[[4]]
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

		output$freqdistrH42 <- renderPlot({
		expType4 <- input$exptype4
                testre<-create_tabs4()
                stackBart <- stackBarH[[4]]
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
	output$vennDiaProtH <- renderPlot({
	expType1 <- input$exptype1
	expType2 <- input$exptype2
	expType3 <- input$exptype3
	expType4 <- input$exptype4
	testre<-create_tabs4()
	if(numfilesH==4){
          	HCDProts11<-HCDProts1[[1]]$protName
		HCDProts21<-HCDProts1[[2]]$protName
		HCDProts31<-HCDProts1[[3]]$protName 
          	HCDProts41<-HCDProts1[[4]]$protName                
		
		InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
          	InterSec13<-Reduce(intersect, list(HCDProts11,HCDProts31))
          	InterSec14<-Reduce(intersect, list(HCDProts11,HCDProts41))
         	InterSec23<-Reduce(intersect, list(HCDProts21,HCDProts31))
          	InterSec24<-Reduce(intersect, list(HCDProts21,HCDProts41))
          	InterSec34<-Reduce(intersect, list(HCDProts31,HCDProts41))
          	InterSec123<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts31))
          	InterSec124<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts41))
          	InterSec134<-Reduce(intersect, list(HCDProts11,HCDProts31,HCDProts41))
          	InterSec234<-Reduce(intersect, list(HCDProts21,HCDProts31,HCDProts41))
          	InterSec1234<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts31,HCDProts41))
		venn.plot <- draw.quad.venn(
          	area1 = length(HCDProts1[[1]]$protName),
		area2 = length(HCDProts1[[2]]$protName),
          	area3 = length(HCDProts1[[3]]$protName),
		area4 = length(HCDProts1[[4]]$protName),
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
                mainT = "Overlaps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3," \nD = ",expType4,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	else if(numfilesH==3){
          	HCDProts11<-HCDProts1[[1]]$protName
		HCDProts21<-HCDProts1[[2]]$protName
		HCDProts31<-HCDProts1[[3]]$protName 
		InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
          	InterSec13<-Reduce(intersect, list(HCDProts11,HCDProts31))
         	InterSec23<-Reduce(intersect, list(HCDProts21,HCDProts31))
          	InterSec123<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts31))
		venn.plot <- draw.triple.venn(
          	area1 = length(HCDProts1[[1]]$protName),
		area2 = length(HCDProts1[[2]]$protName),
          	area3 = length(HCDProts1[[3]]$protName),
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
                mainT = "Overlaps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2," \nC = ",expType3, sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	else if(numfilesH==2){
          	HCDProts11<-as.character(HCDProts1[[1]]$protName)
		HCDProts21<-as.character(HCDProts1[[2]]$protName)
		InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
		venn.plot <- draw.pairwise.venn(
          	area1 = length(HCDProts1[[1]]$protName),
		area2 = length(HCDProts1[[2]]$protName),
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
                mainT = "Overlaps of ADP-ribosylated proteins"
		bottomT = paste("A = ",expType1," \nB = ",expType2,sep="")
		require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
	else if(numfilesH==1){
                HCDProts11<-HCDProts1[[1]]$protName
                InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
                venn.plot <- draw.pairwise.venn(
          	area1 = length(HCDProts1[[1]]$protName),
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
                mainT = "Overlaps of ADP-ribosylated proteins"
                bottomT = paste("A = ",expType1,sep="")
                require(gridExtra)
                grid.arrange(gTree(children=venn.plot), top=mainT, bottom=textGrob(bottomT))
	}
        })
	output$vennDiaPepH <- renderPlot({
	expType1 <- input$exptype1
        expType2 <- input$exptype2
        expType3 <- input$exptype3
        expType4 <- input$exptype4
	testre<-create_tabs4()
		if(numfilesH==4){
                HCDProts11<-HCDPep11[[1]]$pepSeq
		HCDProts21<-HCDPep11[[2]]$pepSeq
		HCDProts31<-HCDPep11[[3]]$pepSeq
		HCDProts41<-HCDPep11[[4]]$pepSeq 
                #HCDProts41<-HCDPep1[[4]]$pep_seq                
		InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
                InterSec13<-Reduce(intersect, list(HCDProts11,HCDProts31))
                InterSec14<-Reduce(intersect, list(HCDProts11,HCDProts41))
                InterSec23<-Reduce(intersect, list(HCDProts21,HCDProts31))
                InterSec24<-Reduce(intersect, list(HCDProts21,HCDProts41))
                InterSec34<-Reduce(intersect, list(HCDProts31,HCDProts41))
                InterSec123<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts31))
                InterSec124<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts41))
                InterSec134<-Reduce(intersect, list(HCDProts11,HCDProts31,HCDProts41))
                InterSec234<-Reduce(intersect, list(HCDProts21,HCDProts31,HCDProts41))
                InterSec1234<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts31,HCDProts41))

		venn.plot <- draw.quad.venn(
                        area1 = length(HCDPep11[[1]]$pepSeq),
			area2 = length(HCDPep11[[2]]$pepSeq),
                        area3 = length(HCDPep11[[3]]$pepSeq),
			area4 = length(HCDPep11[[4]]$pepSeq),
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
	else if(numfilesH==3){
                HCDProts11<-HCDPep11[[1]]$pepSeq
		HCDProts21<-HCDPep11[[2]]$pepSeq
		HCDProts31<-HCDPep11[[3]]$pepSeq 
		

		InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
                InterSec13<-Reduce(intersect, list(HCDProts11,HCDProts31))
                InterSec23<-Reduce(intersect, list(HCDProts21,HCDProts31))
                InterSec123<-Reduce(intersect, list(HCDProts11,HCDProts21,HCDProts31))

		venn.plot <- draw.triple.venn(
                        area1 = length(HCDPep11[[1]]$pepSeq),
			area2 = length(HCDPep11[[2]]$pepSeq),
                        area3 = length(HCDPep11[[3]]$pepSeq),
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
	else if(numfilesH==2){
                HCDProts11<-HCDPep11[[1]]$pepSeq
		HCDProts21<-HCDPep11[[2]]$pepSeq
		InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))
		venn.plot <- draw.pairwise.venn(
                        area1 = length(HCDPep11[[1]]$pepSeq),
			area2 = length(HCDPep11[[2]]$pepSeq),
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
	else if(numfilesH==1){
                HCDProts11<-HCDPep11[[1]]$pepSeq
                InterSec12<-Reduce(intersect, list(HCDProts11,HCDProts21))

                venn.plot <- draw.pairwise.venn(
                        area1 = length(HCDPep11[[1]]$pepSeq),
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
output$heatmapProtH <- renderPlot({
		expType1 <- input$exptype1
        	expType2 <- input$exptype2
        	expType3 <- input$exptype3
        	expType4 <- input$exptype4
		testre<-create_tabs4()
                if(numfilesH==4){
                voltab3 <- merge(x=HCDHeat2[[1]], y=HCDHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=HCDHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab1 <- merge(x=voltab2, y=HCDHeat2[[4]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab1[,-1]
                row.names(voltab) <- voltab1[,1]

                colnames(voltab)<-c(expType1, expType2, expType3, expType4)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
		if(numfilesH==3){
                voltab3 <- merge(x=HCDHeat2[[1]], y=HCDHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=HCDHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab2[,-1]
                row.names(voltab) <- voltab2[,1]

                colnames(voltab)<-c(expType1, expType2, expType3)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
		if(numfilesH==2){
		voltab3 <- merge(x=HCDHeat2[[1]], y=HCDHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
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
	output$heatmapProtHr <- renderPlotly({
		expType1 <- input$exptype1
        	expType2 <- input$exptype2
        	expType3 <- input$exptype3
        	expType4 <- input$exptype4
		testre<-create_tabs4()
                if(numfilesH==4){
                voltab3 <- merge(x=HCDHeat2[[1]], y=HCDHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=HCDHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab1 <- merge(x=voltab2, y=HCDHeat2[[4]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab1[,-1]
                rownames(voltab) <- as.character(voltab1[,1])

                colnames(voltab)<-c(expType1, expType2, expType3, expType4)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
		}
		if(numfilesH==3){
                voltab3 <- merge(x=HCDHeat2[[1]], y=HCDHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab2 <- merge(x=voltab3, y=HCDHeat2[[3]], by="protName", all.x=TRUE, all.y=TRUE)
                voltab <- voltab2[,-1]
                rownames(voltab) <- as.character(voltab2[,1])

                colnames(voltab)<-c(expType1, expType2, expType3)
                voltab.mat<-as.matrix(voltab)
                voltab.mat[is.na(voltab.mat)] <- 0
                }
		if(numfilesH==2){
		voltab3 <- merge(x=HCDHeat2[[1]], y=HCDHeat2[[2]], by="protName", all.x=TRUE, all.y=TRUE)
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

	output$psmBar1H <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs4()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfilesH==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2), length(HCDFiles5[[2]]$pep_seq.2),length(HCDFiles5[[3]]$pep_seq.2),length(HCDFiles5[[4]]$pep_seq.2))
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
		else if(numfilesH==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2), length(HCDFiles5[[2]]$pep_seq.2),length(HCDFiles5[[3]]$pep_seq.2))
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
		else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2), length(HCDFiles5[[2]]$pep_seq.2))
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
		else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2))
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

	 output$psmBar2H <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs4()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfilesH==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2), length(HCDFiles3[[2]]$pep_seq.2),length(HCDFiles3[[3]]$pep_seq.2),length(HCDFiles3[[4]]$pep_seq.2))
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
                else if(numfilesH==3){
		dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2), length(HCDFiles3[[2]]$pep_seq.2),length(HCDFiles3[[3]]$pep_seq.2))
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
                else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2), length(HCDFiles3[[2]]$pep_seq.2))
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
                else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2))
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


	output$pepBar1H <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs4()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfilesH==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		pepCounts = c(length(HCDPep11[[1]]$pepSeq), length(HCDPep11[[2]]$pepSeq),length(HCDPep11[[3]]$pepSeq),length(HCDPep11[[4]]$pepSeq))
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
		else if(numfilesH==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		pepCounts = c(length(HCDPep11[[1]]$pepSeq), length(HCDPep11[[2]]$pepSeq),length(HCDPep11[[3]]$pepSeq))
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
		else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		pepCounts = c(length(HCDPep11[[1]]$pepSeq), length(HCDPep11[[2]]$pepSeq))
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
		else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		pepCounts = c(length(HCDPep11[[1]]$pepSeq))
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

	output$pepBar2H <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs4()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfilesH==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq), length(HCDPep12[[2]]$pepSeq),length(HCDPep12[[3]]$pepSeq),length(HCDPep12[[4]]$pepSeq))
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
                else if(numfilesH==3){
                dat1 <- data.frame(
		 sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq), length(HCDPep12[[2]]$pepSeq),length(HCDPep12[[3]]$pepSeq))
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
                else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq), length(HCDPep12[[2]]$pepSeq))
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
                else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq))
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

	output$protBarH <- renderPlot({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4	
		testre<-create_tabs4()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfilesH==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		protCounts = c(length(HCDProts1[[1]]$protName), length(HCDProts1[[2]]$protName),length(HCDProts1[[3]]$protName), length(HCDProts1[[4]]$protName))
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
		else if(numfilesH==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		protCounts = c(length(HCDProts1[[1]]$protName), length(HCDProts1[[2]]$protName),length(HCDProts1[[3]]$protName))
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
		else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		protCounts = c(length(HCDProts1[[1]]$protName), length(HCDProts1[[2]]$protName))
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
		else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		protCounts = c(length(HCDProts1[[1]]$protName))
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
	output$psmBar1Hr <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs4()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfilesH==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2), length(HCDFiles5[[2]]$pep_seq.2),length(HCDFiles5[[3]]$pep_seq.2),length(HCDFiles5[[4]]$pep_seq.2))
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
		else if(numfilesH==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2), length(HCDFiles5[[2]]$pep_seq.2),length(HCDFiles5[[3]]$pep_seq.2))
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
		else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2), length(HCDFiles5[[2]]$pep_seq.2))
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
		else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(HCDFiles5[[1]]$pep_seq.2))
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

	output$psmBar2Hr <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs4()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfilesH==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2), length(HCDFiles3[[2]]$pep_seq.2),length(HCDFiles3[[3]]$pep_seq.2),length(HCDFiles3[[4]]$pep_seq.2))
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
                else if(numfilesH==3){
		dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2), length(HCDFiles3[[2]]$pep_seq.2),length(HCDFiles3[[3]]$pep_seq.2))
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
                else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2), length(HCDFiles3[[2]]$pep_seq.2))
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
                else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                psmCounts = c(length(HCDFiles3[[1]]$pep_seq.2))
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
	output$pepBar1Hr <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs4()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfilesH==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		pepCounts = c(length(HCDPep11[[1]]$pepSeq), length(HCDPep11[[2]]$pepSeq),length(HCDPep11[[3]]$pepSeq),length(HCDPep11[[4]]$pepSeq))
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
		else if(numfilesH==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		pepCounts = c(length(HCDPep11[[1]]$pepSeq), length(HCDPep11[[2]]$pepSeq),length(HCDPep11[[3]]$pepSeq))
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
		else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		pepCounts = c(length(HCDPep11[[1]]$pepSeq), length(HCDPep11[[2]]$pepSeq))
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
		else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		pepCounts = c(length(HCDPep11[[1]]$pepSeq))
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

	output$pepBar2Hr <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
                testre<-create_tabs4()
                ggplotColours <- function(n = 6, h = c(0, 360) + 15){
                if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
                hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
                }
                if(numfilesH==4){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3, expType4)),
                sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq), length(HCDPep12[[2]]$pepSeq),length(HCDPep12[[3]]$pepSeq),length(HCDPep12[[4]]$pepSeq))
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
                else if(numfilesH==3){
                dat1 <- data.frame(
		 sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq), length(HCDPep12[[2]]$pepSeq),length(HCDPep12[[3]]$pepSeq))
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
                else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq), length(HCDPep12[[2]]$pepSeq))
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
                else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
                pepCounts = c(length(HCDPep12[[1]]$pepSeq))
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

	output$protBarHr <- renderPlotly({
		expType1 <- input$exptype1
                expType2 <- input$exptype2
                expType3 <- input$exptype3
                expType4 <- input$exptype4
		testre<-create_tabs4()
		ggplotColours <- function(n = 6, h = c(0, 360) + 15){
 		if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  		hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
		}	
		if(numfilesH==4){
		dat1 <- data.frame(
    		sam = factor(c(expType1, expType2, expType3, expType4)),
    		sample = factor(c(expType1, expType2, expType3, expType4), levels=c(expType1, expType2, expType3, expType4)),
    		protCounts = c(length(HCDProts3[[1]]$protName), length(HCDProts3[[2]]$protName),length(HCDProts3[[3]]$protName), length(HCDProts3[[4]]$protName))
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
		else if(numfilesH==3){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2, expType3)),
                sample = factor(c(expType1, expType2, expType3), levels=c(expType1, expType2, expType3)),
		protCounts = c(length(HCDProts1[[1]]$protName), length(HCDProts1[[2]]$protName),length(HCDProts1[[3]]$protName))
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
		else if(numfilesH==2){
                dat1 <- data.frame(
                sam = factor(c(expType1, expType2)),
                sample = factor(c(expType1, expType2), levels=c(expType1, expType2)),
		protCounts = c(length(HCDProts1[[1]]$protName), length(HCDProts1[[2]]$protName))
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
		else if(numfilesH==1){
                dat1 <- data.frame(
                sam = factor(c(expType1)),
                sample = factor(c(expType1), levels=c(expType1)),
		protCounts = c(length(HCDProts1[[1]]$protName))
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
	output$goPlot1HBP <- renderPlotly({
		expType1 <- input$exptype1
                pValBP <- input$pval1
		testre<-create_tabs3()
		#globalProt.1 <- globalProt[,1]
		organism <- input$organism
		if(numfilesH>=1){
		localProt <- HCDHeat2[[1]]$protName[-1]
		goterm = "BP"
		#goEnDataBP <- goEnrich(globalProt.1, localProt, organism, goterm, pValBP)
		goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
	output$goPlot1HMF <- renderPlotly({
                expType1 <- input$exptype1
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=1){
                localProt <- HCDHeat2[[1]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[1]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[1]])
		goterm = "MF"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot1HCC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=1){
                localProt <- HCDHeat2[[1]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[1]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[1]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
	output$goPlot2HBP <- renderPlotly({
                expType2 <- input$exptype2
		pValBP <- input$pval1
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=2){
                localProt <- HCDHeat2[[2]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
	output$goPlot2HMF <- renderPlotly({
                expType2 <- input$exptype2
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=2){
                localProt <- HCDHeat2[[2]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot2HCC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=2){
                localProt <- HCDHeat2[[2]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[2]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[2]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
	output$goPlot3HBP <- renderPlotly({
                expType3 <- input$exptype3
		pValBP <- input$pval1
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=3){
                localProt <- HCDHeat2[[3]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
	output$goPlot3HMF <- renderPlotly({
                expType3 <- input$exptype3
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=3){
                localProt <- HCDHeat2[[3]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot3HCC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH>=3){
                localProt <- HCDHeat2[[3]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[3]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[3]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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
	output$goPlot4HBP <- renderPlotly({
                expType4 <- input$exptype4
		pValBP <- input$pval1
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH==4){
                localProt <- HCDHeat2[[4]]$protName[-1]
                goterm = "BP"
                goEnDataBP <- goEnrich(localProt, organism, goterm, pValBP)
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
	output$goPlot4HMF <- renderPlotly({
                expType4 <- input$exptype4
		pValMF <- input$pval2
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH==4){
                localProt <- HCDHeat2[[4]]$protName[-1]
                goterm = "MF"
                goEnDataMF <- goEnrich(localProt, organism, goterm, pValMF)
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
	output$goPlot4HCC <- renderPlotly({
                expType1 <- input$exptype1
		pValCC <- input$pval3
                testre<-create_tabs3()
                #globalProt.1 <- globalProt[,1]
                organism <- input$organism
                if(numfilesH==4){
                localProt <- HCDHeat2[[4]]$protName[-1]
                #localProt.1 <- str_decapitalize(msmsHeat2[[4]]$protName[-1])
                #localProt <- str_ucfirst(localProt.1)
                #geneList <- t(msmsHeat2[[4]])
                goterm = "CC"
                #goEnData3 <- goEnrich(globalProt.1, localProt, geneList, organism)
                goEnDataCC <- goEnrich(localProt, organism, goterm, pValCC)
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

})




