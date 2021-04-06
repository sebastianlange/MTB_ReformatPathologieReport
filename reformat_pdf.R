library(tabulizer)
library(dplyr)
library(pdftools)
library(readr)
encoding="latin1"
Sys.setlocale("LC_ALL", "de_DE.UTF-8")

files.ls=list.files(pattern=".pdf$")
files.ls=files.ls[grep("Histopathologischer Befundbericht",files.ls)]
files.ls=rev(files.ls)
#files.ls=files.ls[105:length(files.ls)]

for (filename in files.ls) 
{
  tryCatch({
	rm("Metadata")
	rm("Variants")
	rm("text")
	rm("TMB")
	rm("MSI")
	rm("CNV")
	rm("Translokationen")

	Metadata=data.frame("CCC"=NA,Vorname=NA, Nachname=NA, Geburtsdatum=NA, Eingangsnummer=NA, "Entität"=NA,"Blockmaterial"=NA,"Tumorzellgehalt"=NA,"Qualität"=NA,"Kit"=NA)
	text=pdf_text(filename) %>% readr::read_lines()
	Metadata[,"CCC"]=strsplit(strsplit(filename,"-")[[1]][2],"_")[[1]][1]
	Metadata[,"Vorname"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Vorname",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Nachname"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Name",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Geburtsdatum"]=gsub("\\s+"," ",trimws(strsplit(text[grep("geb. am",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Eingangsnummer"]=gsub("\\s+"," ",trimws(strsplit(text[grep("E-Nr",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Entität"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Klinische Angaben:",text)],":")[[1]][2],which="left"))
	Metadata[,"Blockmaterial"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Blockmaterial",text)],":")[[1]][2],which="left"))
	Metadata[,"Tumorzellgehalt"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Tumorzellgehalt des untersuchten Materials",text)],":")[[1]][2],which="left"))
	Metadata[,"Qualität"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Qualität der extrahierten Nukleinsäuren",text)],":")[[1]][2],which="left"))
	Metadata[,"Kit"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Durchführung der Analyse mittel",text)],":")[[1]][2],which="left"))
 	Metadata=rbind(colnames(Metadata),Metadata[1,])
	txt <- pdf_text(filename)

	Variants <- grep("COSMIC", txt, ignore.case = TRUE)[1]
	TMB <- grep("Anzahl Alterationen zur TMB-Bestimmung",txt, ignore.case = TRUE)[1]
	MSI <- grep("Anzahl analysierbarer Mikrosatellitenstellen ",txt, ignore.case = TRUE)[1]
	CNV <- grep("Amplifikation/Deletion", txt, ignore.case = TRUE)[1]
	Translokationen <- grep("Fusionstranskript", txt, ignore.case = TRUE)[1]

	table <- extract_tables(filename, method="stream",pages=TMB)
	TMB=table[[grep("Anzahl Alterationen zur",table)]]
	TMB=data.frame(TMB)
	TMB=t(TMB)
	rownames(TMB)=c()

	table <- extract_tables(filename,method = "stream", pages=MSI)
	MSI=table[[grep("Anzahl analysierbarer Mikrosatellitenstellen",table)]]
	MSI=data.frame(MSI)
	MSI=t(MSI)

	if (is.na(Variants)) {
		Variants="In allen untersuchten Genbereichen (s.u.) konnten keine Mutationen nachgewiesen werden."
	} else {
		table <- extract_tables(filename,method = "stream", pages=Variants)
		Variants=table[[grep("COSMIC",table)]]
		if(ncol(Variants)==4){
			Variants[1,3]=gsub("\r","",Variants[1,3])
			if(Variants[2,1]==""){Variants=Variants[-2,]}
			Variants[1,3]="COSMIC-Datenbank v90"
			Variants=data.frame(Variants)
			colnames(Variants)=Variants[1,]
			Variants=Variants[2:nrow(Variants),]
			Variants[,"Nukleotid"]=data.frame(sapply(lapply(strsplit(Variants[,"Ergebnis"],","),`[`, 1),c))
			Variants[,"Protein"]=data.frame(sapply(lapply(strsplit(Variants[,"Ergebnis"],","),`[`, 2),c))
			Variants[,"Protein"]=data.frame(sapply(lapply(strsplit(Variants[,"Ergebnis"]," "),`[`, 2),c))
			Variants[,"Allelfrequenz"]=data.frame(sapply(lapply(strsplit(Variants[,"Ergebnis"]," "),`[`, 3),c)) 
			Variants[,"Allelfrequenz"]=gsub("\\(","",Variants[,"Allelfrequenz"])
			Variants[,"Allelfrequenz"]=gsub("\\)","",Variants[,"Allelfrequenz"])
			Variants[,"Allelfrequenz"]=gsub("%","",Variants[,"Allelfrequenz"])
			Variants=Variants[,c("Gen","Nukleotid","Protein","Allelfrequenz",Variants[1,3],"Kommentar")]
		} else if(ncol(Variants)>4){
			if(Variants[2,2]=="nummer"||Variants[2,1]=="nummer"||Variants[2,3]=="Datenbank v90"||Variants[2,3]=="v90"){
				Variants[1,1]=paste0(Variants[1,1],Variants[2,1])
				Variants[1,2]=paste0(Variants[1,2],Variants[2,2])
				Variants[1,6]=paste0(Variants[1,6],Variants[2,6])
				Variants[1,7]=paste0(Variants[1,7],Variants[2,7])
				Variants=Variants[-c(1,2),]
				colnames(Variants)=c("Gen","Referenznummer","Exon","cDNA","Protein","Allelfrequenz","COSMIC-Datenbank v90","Klasse")
			}
		}
	}

	if (is.na(CNV)) {
		CNV="In allen untersuchten Genbereichen (s.u.) konnten keine Amplifikationen oder Deletionen nachgewiesen werden."
	} else {
		table <- extract_tables(filename,method = "stream", pages=CNV)
		CNV=table[[grep("Amplifikation/Deletion",table)]]
		if(CNV[2,2]=="nummer"){
		CNV[2,1]=paste0(CNV[1,2],CNV[2,2])
		CNV=CNV[-2,]
		}
	}

	if (is.na(Translokationen)) {
		Translokationen="In allen untersuchten Genbereichen (s.u.) konnten keine Translokationen nachgewiesen werden."
	} else {
		table <- extract_tables(filename, pages=Translokationen)
		Translokationen=table[[grep("Fusion",table)]]
		if(Translokationen[2,2]=="nummer"){
			Translokationen[1,2]=paste0(Translokationen[1,2],Translokationen[2,2])
			Translokationen[1,4]=paste0(Translokationen[1,4],Translokationen[2,4])
			Translokationen[1,5]=paste0(Translokationen[1,5],Translokationen[2,5])
			Translokationen=Translokationen[-2,]
			Translokationen=Translokationen[,1:6]
			#colnames(Translokationen)=Translokationen[1,]
			#Translokationen=Translokationen[-1,]
			Translokationen=Translokationen[-3,]
		}
	}

	write.table(cbind(data.frame(CCC=replicate(NROW(t(Metadata)),Metadata[2,1])),t(Metadata)),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(cbind(data.frame(CCC=replicate(NROW(t(TMB)),Metadata[2,1])),t(TMB)),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	write.table(cbind(data.frame(CCC=replicate(NROW(t(MSI)),Metadata[2,1])),t(MSI)),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	write.table(cbind(data.frame(CCC=replicate(NROW(Variants),Metadata[2,1])),Variants),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	if(CNV=="In allen untersuchten Genbereichen (s.u.) konnten keine Amplifikationen oder Deletionen nachgewiesen werden."){
		write.table(cbind(data.frame(CCC=replicate(NROW(CNV),Metadata[2,1])),CNV),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	} else {
		write.table(cbind(data.frame(CCC=replicate(NROW(CNV)-1,Metadata[2,1])),CNV[2:NROW(CNV)]),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	}
	write.table(cbind(data.frame(CCC=replicate(NROW(t(Translokationen)),Metadata[2,1])),t(Translokationen)),paste0("CCC-",Metadata[2,1],".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	
	write.table(Metadata[,1:4],paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(Metadata[,5:9],paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(TMB,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(MSI,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE,append=TRUE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	if (nrow(t(Variants)) == 1 )
	{
		write.table(Variants,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
	} else {
		write.table(Variants,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=TRUE)
	}
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(CNV,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(Translokationen,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)

	write.table(Metadata[,1:4],paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=TRUE)
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
	write.table(Metadata[,5:9],paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=TRUE,append=TRUE)
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
	write.table(TMB,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=FALSE,append=TRUE)
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
	write.table(MSI,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=FALSE,append=TRUE)
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
		if (nrow(t(Variants)) == 1 )
	{
		write.table(Variants,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=FALSE, append=TRUE)
	} else {
		write.table(Variants,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=TRUE, append=TRUE)
	}
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
	write.table(CNV,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=FALSE, append=TRUE)
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
	write.table(Translokationen,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=FALSE, append=TRUE)
	}, error=function(e){
		cat("ERROR :",conditionMessage(e), "\n")
		file.create(paste0(filename,".txt"))
		file.create(paste0(filename,"_SAP.txt"))
		file.create(paste0("CCC-",Metadata[2,1],".txt"))
	})
}