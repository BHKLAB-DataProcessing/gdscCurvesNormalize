#!/usr/bin/env Rscript

library(PharmacoGx)
library(PharmacoGxPrivate)
library(data.table)
# library(readxl)

args = commandArgs(trailingOnly=TRUE)

version <- args[[1]]

switch(version, v1 = {
	myInFile <- "GDSC1_public_raw_data_17Jul19.csv"
	myOutPrefix <- "/pfs/out/gdscv1"
	control.column <- "NC-0"
	}, v2 = {
	myInFile <- "GDSC2_public_raw_data_17Jul19.csv"
	myOutPrefix <- "/pfs/out/gdscv2"
	control.column <- "NC-1"
	})


gdsc_sens <- fread(file.path("/pfs/in", myInFile))

## BARCODE is the plate/batch ID. Everything with same barcode is a techinical replicate
## NB: v1 has only NC-0 controls, v2 has NC-0 and NC-1
## NB: One cell line per plate is plated. 

# gdscv1_sens[,unique(BARCODE)]

gdsc_sens[,TAG.1 := sapply(strsplit(TAG, split="-"), function(x) return(x[[1]]))]


## Taking mean if more than 5 measurements, otherwise median

gdsc_sens[,Viability := {
	control <- .SD[TAG==..control.column, ifelse(length(INTENSITY)>5,mean(INTENSITY),median(INTENSITY))]
	background <- .SD[TAG.1=="B",ifelse(length(INTENSITY)>5,mean(INTENSITY),median(INTENSITY))]
	Viability <- (INTENSITY - background)/(control-background)*100
	list(Viability = Viability)}, .(MASTER_CELL_ID, BARCODE) ]


# gdscv2_sens[,Viability := {
# 	control <- .SD[TAG=="NC-1", ifelse(length(INTENSITY)>5,mean(INTENSITY),median(INTENSITY))]
# 	background <- .SD[TAG.1=="B",ifelse(length(INTENSITY)>5,mean(INTENSITY),median(INTENSITY))]
# 	Viability <- (INTENSITY - background)/(control-background)
# 	list(Viability = Viability)}, .(MASTER_CELL_ID, BARCODE) ]

load("drugInfo.RData")
load("cellInfo.RData")

gdsc_sens[,drugid := ..drug.info[match(DRUG_ID, ..drug.info$DRUG_ID), "unique.drugid"]]
gdsc_sens[,cellid := ..cell.info[match(COSMIC_ID, ..cell.info$COSMIC.identifier), "unique.cellid"]]


# gdscv1_sens[grepl("L|R", x=TAG.1),]

# gdscv2_sens[grepl("L|R", x=TAG.1),]

## Unique ID formulation will be: cellid_drugid_plate_seeding_dens_assay_duration_tag.1
## Will add tag later because I want to average the different technical replicates with the same concentration 
## range but keep the ones with different ranges separate
## Changed my mind, using different tags differently seems unfair

## Tag.1 is the "curve id"

gdsc_exps <- gdsc_sens[grepl("L|R", x=TAG.1),]
rm(gdsc_sens); gc()

gdsc_exps[,exp_id := paste(cellid, drugid, BARCODE, "SeedDens", SEEDING_DENSITY, "assay", ASSAY, "dur", DURATION, sep="_")]

gdsc_raw <- gdsc_exps[,.(CONC,Viability, exp_id)]

gdsc_raw_list <- split(gdsc_raw, by="exp_id")

for(i in seq_along(gdsc_raw_list)){
	gdsc_raw_list[[i]][["exp_id"]] <- NULL
}

for(i in seq_along(gdsc_raw_list)){
	gdsc_raw_list[[i]] <- gdsc_raw_list[[i]][,lapply(.SD,median),CONC]
}


## Currently treating curves that have same plate but different ranges as a single experiment

## Calculate sizes needed for the raw array

max.num.conc <- max(sapply(gdsc_raw_list, NROW))

gdsc_raw <- array(data=NA_real_, dim=c(length(gdsc_raw_list),max.num.conc,2))

for(i in seq_along(gdsc_raw_list)){
	setorder(gdsc_raw_list[[i]], CONC)
	gdsc_raw[i,seq_len(NROW(gdsc_raw_list[[i]])),1] <- gdsc_raw_list[[i]][,CONC]
	gdsc_raw[i,seq_len(NROW(gdsc_raw_list[[i]])),2] <- gdsc_raw_list[[i]][,Viability]
}

dimnames(gdsc_raw)[[2]] <- paste("dose", 1:NCOL(gdsc_raw), sep=".")
dimnames(gdsc_raw)[[3]] <- c("Dose", "Viability") 


## Now lets make the info tables complete. 


gdsc_exps[,SITE := ..drug.info$SCREENING_SITE[match(gdsc_exps$DRUG_ID, ..drug.info$DRUG_ID)]]

gdsc_exps[, c("MAX.CONC", "MIN.CONC") :=  .(max(CONC), min(CONC)), exp_id]

gdsc_info <- gdsc_exps[,.(cellid, drugid, MAX.CONC, MIN.CONC, ASSAY, DURATION, BARCODE, SITE, SCAN_ID, DATE_CREATED, SCAN_DATE, CELL_ID, MASTER_CELL_ID, COSMIC_ID, CELL_LINE_NAME, SEEDING_DENSITY, DRUGSET_ID, DRUG_ID, "TAG.1s" = paste(unique(TAG.1), collapse="///")), exp_id]

gdsc_info <- unique(gdsc_info)

gdsc_info <- as.data.frame(gdsc_info)


saveRDS(gdsc_info, file=paste0(myOutPrefix, "_sens_info.rds"))
saveRDS(gdsc_raw, file=paste0(myOutPrefix, "_sens_raw.rds"))


sens.info <- gdsc_info
sens.raw <- gdsc_raw


rownames(sens.info) <- sens.info$exp_id
rownames(sens.raw) <- sens.info$exp_id

sens.raw.x <- parallel::splitIndices(nrow(sens.raw), floor(nrow(sens.raw)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(sens.raw.x)){

  slce <- sens.raw[sens.raw.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/gdsc_raw_sens_", i, ".rds"))

}


