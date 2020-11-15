library(tabulizer)
library(dplyr)
library(pdftools)
library(readr)
encoding="latin1"
Sys.setlocale("LC_ALL", "de_DE.UTF-8")

files.ls=list.files(pattern=".pdf$")
files.ls=files.ls[grep("Histopathologischer Befundbericht",files.ls)]
files.ls=files.ls[135:140]

for (filename in files.ls) {
	rm("Metadata")
	rm("Variants")
	rm("text")
	rm("TMB")
	rm("MSI")
	rm("CNV")
	rm("Translokationen")

	Metadata=data.frame(Vorname=NA, Nachname=NA, Geburtsdatum=NA, Eingangsnummer=NA, "Entität"=NA,"Blockmaterial"=NA,"Tumorzellgehalt"=NA,"Qualität"=NA,"Kit"=NA)
	text=pdf_text(filename) %>% readr::read_lines()
	Metadata[,"Vorname"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Vorname",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Nachname"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Name",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Geburtsdatum"]=gsub("\\s+"," ",trimws(strsplit(text[grep("geb. am",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Eingangsnummer"]=gsub("\\s+"," ",trimws(strsplit(text[grep("E-Nr",text)][2],":")[[1]][2],which="left"))
	Metadata[,"Entität"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Klinische Angaben:",text)],":")[[1]][2],which="left"))
	Metadata[,"Blockmaterial"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Blockmaterial",text)],":")[[1]][2],which="left"))
	Metadata[,"Tumorzellgehalt"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Tumorzellgehalt des untersuchten Materials",text)],":")[[1]][2],which="left"))
	Metadata[,"Qualität"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Qualität der extrahierten Nukleinsäuren",text)],":")[[1]][2],which="left"))
	Metadata[,"Kit"]=gsub("\\s+"," ",trimws(strsplit(text[grep("Durchführung der Analyse mittel",text)],":")[[1]][2],which="left"))
 
	txt <- pdf_text(filename)

	Variants <- grep("COSMIC", txt, ignore.case = TRUE)[1]
	TMB <- grep("Anzahl Alterationen zur TMB-Bestimmung",txt, ignore.case = TRUE)[1]
	MSI <- grep("Anzahl analysierbarer Mikrosatellitenstellen ",txt, ignore.case = TRUE)[1]
	CNV <- grep("Amplifikation/Deletion", txt, ignore.case = TRUE)[1]
	Translokationen <- grep("Gen 1", txt, ignore.case = TRUE)[1]

	table <- extract_tables(filename,method = "lattice", pages=TMB)
	TMB=table[[grep("Anzahl Alterationen zur",table)]]
	TMB=data.frame(TMB)
	colnames(TMB)=TMB[,1]
	TMB[1,1]=TMB[1,2]
	TMB[1,2]=TMB[2,2]
	TMB=TMB[1,]

	table <- extract_tables(filename,method = "stream", pages=MSI)
	MSI=table[[grep("Anzahl analysierbarer Mikrosatellitenstellen",table)]]
	MSI=data.frame(MSI)
	MSI=t(MSI)

	if (is.na(Variants)) {
		Variants="In allen untersuchten Genbereichen (s.u.) konnten keine Mutationen nachgewiesen werden."
	} else {
		table <- extract_tables(filename,method = "stream", pages=Variants)
		Variants=table[[grep("COSMIC",table)]]
		Variants[1,3]=gsub("\r","",Variants[1,3])
		Variants=Variants[-2,]
		Variants[,3]="COSMIC-Datenbank v90"
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
		Variants=Variants[,c("Gen","Nukleotid","Protein","Allelfrequenz",colnames(Variants)[3],"Kommentar")]
	}

	if (is.na(CNV)) {
		CNV="In allen untersuchten Genbereichen (s.u.) konnten keine Amplifikationen oder Deletionen nachgewiesen werden."
	} else {
		table <- extract_tables(filename,method = "stream", pages=CNV)
		CNV=table[[grep("Amplifikation/Deletion",table)]]
	}

	if (is.na(Translokationen)) {
		Translokationen="In allen untersuchten Genbereichen (s.u.) konnten keine Translokationen nachgewiesen werden."
	} else {
		table <- extract_tables(filename,method = "stream", pages=Translokationen)
		Translokationen=table[[grep("Fusion",table)]]
	}

	write.table(Metadata[,1:4],paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(Metadata[,5:9],paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE,append=TRUE)
	cat("\n", file = paste0(filename,".txt"), append = TRUE)
	write.table(TMB,paste0(filename,".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE,append=T)
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
	write.table(TMB,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=TRUE,append=T)
	cat("\n", file = paste0(filename,"_SAP.txt"), append = TRUE)
	write.table(MSI,paste0(filename,"_SAP.txt"), quote=FALSE, sep=" \\\\ ", row.names=FALSE, col.names=F,append=TRUE)
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
}