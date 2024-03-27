####################### Lung Cancer ########################
## obtain significant SNVs found in lung cancer
snplist = fread("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-single.txt")
snplist$Allele1 = toupper(snplist$Allele1)
snplist$Allele2 = toupper(snplist$Allele2)
########## WGS #############
ukb = fread("/home/biostat/DiskB/jldai/BIB_revise/ukb_biochemistry.txt")
ukb = select(ukb,id_WGS,LC)

wgs_150k = fread("/home/biostat/DiskB/jldai/Part1-相关性分析/WGS_150k.txt",header = F)
wgs_150k = unite(wgs_150k,id_WGS,c("V1","V2"),sep = "_")
sample = wgs_150k$id_WGS
sample = sample[-1]

ukb_wgs = ukb[ukb$id_WGS %in% sample,]
ukb_wgs = separate(ukb_wgs,id_WGS,into = c("#FID","IID"),sep = "_")

case = ukb_wgs[ukb_wgs$LC==1,c(1,2)]
case$SEX = NA
control = ukb_wgs[ukb_wgs$LC==0,c(1,2)]
control$SEX = NA
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case.txt",col.names = T,row.names = F,quote = F)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control.txt",col.names = T,row.names = F,quote = F)

snplist$MarkerName = paste0("chr",snplist$MarkerName,":SG")
extract = snplist[,MarkerName]
write.table(snplist,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/LC-single.txt",col.names = T,row.names = F,quote = F)
write.table(extract,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/LC-extract.txt",col.names = F,row.names = F,quote = F)

#### plink-WGS ####
## case ##
for chr in 5 14 15 16; do
plink2 --pfile /home/biostat/DiskB/WGS150k/WGS_150k.chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/LC-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr${chr}
done

## control ##
for chr in 5 14 15 16; do
plink2 --pfile /home/biostat/DiskB/WGS150k/WGS_150k.chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/LC-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr${chr}
done
########## Onco #################
pheno = fread("/home/sshen/Disk_m1/DiskB/onco/phenotype_29097_20220413.txt")
ID = pheno[pheno$ctype!="Control",ID]
case = data.frame(v1=0,ID)
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/case.txt",col.names = F,row.names = F,quote = F)
ID = pheno[pheno$ctype=="Control",ID]
control = data.frame(v1=0,ID)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/control.txt",col.names = F,row.names = F,quote = F)

snp_list = snplist
snp_list$MarkerName = paste0("chr",snplist$MarkerName)
snp_list = unite(snp_list,ID,c("MarkerName","Allele1","Allele2"),sep = ":")
write.table(snp_list,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/LC-extract.txt",col.names = F,row.names = F,quote = F)
#### plink-Onco ####
## case ##
for chr in 5 14 15 16; do
plink2 --bfile /home/sshen/Disk_m1/DiskB/onco/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/LC-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/case/chr${chr}
done

## control ##
for chr in 5 14 15 16; do
plink2 --bfile /home/sshen/Disk_m1/DiskB/onco/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/LC-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/control/chr${chr}
done

########## TRICL ########
pheno = fread("/home/sshen/Disk_m1/DiskB/TRICL/covar_TRICL_20220420.txt")
ID = pheno[pheno$ctype!="Control",ID]
case = data.frame(v1=0,ID)
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/case.txt",col.names = F,row.names = F,quote = F)
ID = pheno[pheno$ctype=="Control",ID]
control = data.frame(v1=0,ID)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/control.txt",col.names = F,row.names = F,quote = F)

snp_list = snplist
snp_list$MarkerName = paste0("chr",snplist$MarkerName)
snp_list = unite(snp_list,ID,c("MarkerName","Allele1","Allele2"),sep = ":")
write.table(snp_list,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/LC-extract.txt",col.names = F,row.names = F,quote = F)

#### plink-TRICL ####
## case ##
for chr in 5 14 15 16; do
plink2 --bfile /home/sshen/Disk_m1/DiskB/TRICL/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/LC-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/case/chr${chr}
done

## control ##
for chr in 5 14 15 16; do
plink2 --bfile /home/sshen/Disk_m1/DiskB/TRICL/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/LC-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/control/chr${chr}
done

########## PLCO #############
pheno = fread("/home/sshen/Disk_m1/DiskB/PLCO/plco_euro_20230908.txt")
ID = pheno[pheno$ctype!="Control",ID]
case = data.frame(v1=0,ID)
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/case.txt",col.names = F,row.names = F,quote = F)
ID = pheno[pheno$ctype=="Control",ID]
control = data.frame(v1=0,ID)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/control.txt",col.names = F,row.names = F,quote = F)

snp_list = snplist
snp_list$MarkerName = paste0("chr",snplist$MarkerName)
snp_list = unite(snp_list,ID,c("MarkerName","Allele1","Allele2"),sep = ":")
write.table(snp_list,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/LC-extract.txt",col.names = F,row.names = F,quote = F)
#### plink-PLCO ####
## case ##
for chr in 5 14 15 16; do
plink2 --bfile /home/sshen/Disk_m1/DiskB/PLCO/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/LC-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/case/chr${chr}
done

## control ##
for chr in 5 14 15 16; do
plink2 --bfile /home/sshen/Disk_m1/DiskB/PLCO/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/LC-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/control/chr${chr}
done

########## test for cases #########
#### WGS vs Onco ####
chr = c(14,15,16)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/case/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[match(other_id,wgs_id)]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/case/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id

    wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-WGS_Onco_case.txt",col.names = T,row.names = F,quote = F)

#### WGS vs TRICL ####
chr = c(14,15,16)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/case/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[match(other_id,wgs_id)]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/case/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id

    wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-WGS_TRICL_case.txt",col.names = T,row.names = F,quote = F)

#### WGS vs PLCO ####
chr = c(14,15,16)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/case/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[match(other_id,wgs_id)]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/case/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/case/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id
  
  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-WGS_PLCO_case.txt",col.names = T,row.names = F,quote = F)

########## test for control #######
#### WGS vs Onco ####
chr = c(14,15,16)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/control/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[match(other_id,wgs_id)]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/onco/control/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id

  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-WGS_Onco_control.txt",col.names = T,row.names = F,quote = F)


#### WGS vs TRICL ####
chr = c(14,15,16)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/control/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[match(other_id,wgs_id)]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/TRICL/control/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id
  
  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-WGS_TRICL_control.txt",col.names = T,row.names = F,quote = F)

#### WGS vs PLCO ####
chr = c(14,15,16)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/control/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[match(other_id,wgs_id)]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs/control/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/PLCO/control/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id
  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/LC-WGS_PLCO_control.txt",col.names = T,row.names = F,quote = F)

####################### Ovarian Cancer ########################
## obtain significant SNVs found in Ovarian Cancer
snplist = fread("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-single.txt")
snplist$Allele1 = toupper(snplist$Allele1)
snplist$Allele2 = toupper(snplist$Allele2)

########## WGS #############
ukb = fread("/home/biostat/DiskB/jldai/BIB_revise/ukb_biochemistry.txt")
ukb = select(ukb,id_WGS,OV)

wgs_150k = fread("/home/biostat/DiskB/jldai/Part1-相关性分析/WGS_150k.txt",header = F)
wgs_150k = unite(wgs_150k,id_WGS,c("V1","V2"),sep = "_")
sample = wgs_150k$id_WGS
sample = sample[-1]

ukb_wgs = ukb[ukb$id_WGS %in% sample,]
ukb_wgs = separate(ukb_wgs,id_WGS,into = c("#FID","IID"),sep = "_")

case = subset(ukb_wgs,OV==1)[,1:2]
case$SEX = NA
control = subset(ukb_wgs,OV==0)[,1:2]
control$SEX = NA
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case.txt",col.names = T,row.names = F,quote = F)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control.txt",col.names = T,row.names = F,quote = F)

snplist$MarkerName = paste0("chr",snplist$MarkerName,":SG")
extract = snplist[,MarkerName]
write.table(snplist,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/OV-single.txt",col.names = T,row.names = F,quote = F)
write.table(extract,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/OV-extract.txt",col.names = F,row.names = F,quote = F)

#### plink-WGS ####
## case ##
for chr in 6 9 10 17; do
plink2 --pfile /home/biostat/DiskB/WGS150k/WGS_150k.chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/OV-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr${chr}
done

## control ##
for chr in 6 9 10 17; do
plink2 --pfile /home/biostat/DiskB/WGS150k/WGS_150k.chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/OV-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr${chr}
done
########## CIMBA #################
pheno = fread("/home/sshen/Disk_m1/Disk_m3/CIMBA/covar_CIMBA_OV_euro.txt")
ID = pheno[pheno$OV==1,SUBJECT_ID]
case = data.frame(v1=0,ID)
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/case.txt",col.names = F,row.names = F,quote = F)
ID = pheno[pheno$OV==0,SUBJECT_ID]
control = data.frame(v1=0,ID)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/control.txt",col.names = F,row.names = F,quote = F)

snp_list = snplist
snp_list$MarkerName = paste0("chr",snplist$MarkerName)
snp_list = unite(snp_list,ID,c("MarkerName","Allele1","Allele2"),sep = ":")
write.table(snp_list,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/OV-extract.txt",col.names = F,row.names = F,quote = F)

#### plink-CIMBA ####
## case ##
for chr in 6 9 10 17; do
plink2 --bfile /home/sshen/Disk_m1/Disk_m3/CIMBA/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/OV-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/case/chr${chr}
done

## control ##
for chr in 6 9 10 17; do
plink2 --bfile /home/sshen/Disk_m1/Disk_m3/CIMBA/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/OV-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/control/chr${chr}
done
########## FIOC-ONCO #################
pheno = fread("/home/sshen/Disk_m1/Disk_m3/OV/pheno_onco_euro.txt")
ID = pheno[pheno$OV==1,ID]
case = data.frame(v1=0,ID)
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/case.txt",col.names = F,row.names = F,quote = F)
ID = pheno[pheno$OV==0,ID]
control = data.frame(v1=0,ID)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/control.txt",col.names = F,row.names = F,quote = F)

snp_list = snplist
snp_list$MarkerName = paste0("chr",snplist$MarkerName)
snp_list = unite(snp_list,ID,c("MarkerName","Allele1","Allele2"),sep = ":")
write.table(snp_list,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/OV-extract.txt",col.names = F,row.names = F,quote = F)

#### plink-FIOC-ONCO ####
## case ##
for chr in 6 9 10 17; do
plink2 --bfile /home/sshen/Disk_m1/Disk_m3/OV/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/OV-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/case/chr${chr}
done

## control ##
for chr in 6 9 10 17; do
plink2 --bfile /home/sshen/Disk_m1/Disk_m3/OV/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/OV-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/control/chr${chr}
done


########## FIOC-EXOME #################
pheno = fread("/home/sshen/Disk_m1/Disk_m3/OVe/pheno_exome_euro.txt")
ID = pheno[pheno$OV==1,ID]
case = data.frame(v1=0,ID)
write.table(case,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/case.txt",col.names = F,row.names = F,quote = F)
ID = pheno[pheno$OV==0,ID]
control = data.frame(v1=0,ID)
write.table(control,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/control.txt",col.names = F,row.names = F,quote = F)

snp_list = snplist
snp_list$MarkerName = paste0("chr",snplist$MarkerName)
snp_list = unite(snp_list,ID,c("MarkerName","Allele1","Allele2"),sep = ":")
write.table(snp_list,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/OV-extract.txt",col.names = F,row.names = F,quote = F)

#### plink-FIOC-EXOME ####
## case ##
for chr in 6 9 10 17; do
plink2 --bfile /home/sshen/Disk_m1/Disk_m3/OVe/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/case.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/OV-extract.txt \
--make-bed \
-out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/case/chr${chr}
done

## control ##
for chr in 6 9 10 17; do
plink2 --bfile /home/sshen/Disk_m1/Disk_m3/OVe/plink/chr${chr} \
--keep /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/control.txt \
--extract /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/OV-extract.txt \
--make-bed \
--out /home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/control/chr${chr}
done



########## test for cases #########
#### WGS vs CIMBA ####
chr = c(6,9,17)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/case/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[na.omit(match(other_id,wgs_id))]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/case/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id

  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-WGS_CIMBA_case.txt",col.names = T,row.names = F,quote = F)

#### WGS vs FIOC-ONCO ####
chr = c(6,9,17)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/case/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[na.omit(match(other_id,wgs_id))]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/case/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id
  
  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-WGS_FOCI-onco_case.txt",col.names = T,row.names = F,quote = F)

#### WGS vs FIOC-EXOME ####
chr = c(6,9,17)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/case/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[na.omit(match(other_id,wgs_id))]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/case/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/case/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id

  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-WGS_FOCI-exome_case.txt",col.names = T,row.names = F,quote = F)

########## test for controls #########
#### WGS vs CIMBA ####
chr = c(6,9,17)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/control/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[na.omit(match(other_id,wgs_id))]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/CIMBA/control/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id

  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-WGS_CIMBA_control.txt",col.names = T,row.names = F,quote = F)

#### WGS vs FIOC-ONCO ####
chr = c(6,9,17)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/control/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[na.omit(match(other_id,wgs_id))]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-onco/control/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id
  
  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-WGS_FOCI-onco_control.txt",col.names = T,row.names = F,quote = F)

#### WGS vs FIOC-EXOME ####
chr = c(6,9,17)
result_all = NULL
for ( i in chr ) {
  wgs_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr",i,".bim"))
  wgs_bim$V2 = gsub(":SG","",wgs_bim$V2)
  wgs_bim = unite(wgs_bim,ID,c("V2","V6","V5"),sep=":")
  wgs_id = wgs_bim$ID
  
  other_bim = fread(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/control/chr",i,".bim"))
  other_id = other_bim$V2
  id = wgs_id[na.omit(match(other_id,wgs_id))]
  
  wgs_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/wgs_ov/control/chr",i),simple_names=T))
  other_bed = as.matrix(BEDMatrix(paste0("/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/FOCI-exome/control/chr",i),simple_names=T))
  colnames(wgs_bed) = wgs_id
  
  wgs_bed = wgs_bed[drop=FALSE,,match(id,wgs_id)]
  
  for (j in 1:length(id)) {
    wgs = wgs_bed[,j]
    other = other_bed[,j] 
    a = length(wgs[wgs==0])
    b = length(other[other==0])
    c = length(wgs[wgs==1])
    d = length(other[other==1])
    e = length(wgs[wgs==2])
    f = length(other[other==2])
    
    compare = matrix(c(a,b,c,d,e,f),nr=2,
                     dimnames = list(c("wgs","other"),c("type0","type1","type2")))
    result = fisher.test(compare,simulate.p.value = TRUE)
    snp_id = id[j]
    p = result$p.value
    temp = data.frame(snp_id,p)
    result_all = rbind(result_all,temp)
  }
}
write.table(result_all,file = "/home/biostat/DiskB/jldai/BIB_revise/sensitive_analysis/OV-WGS_FOCI-exome_control.txt",col.names = T,row.names = F,quote = F)




