#!/usr/bin/env Rscript

# ========================
# === Sample Filtering === 
# ========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
basketFile = args[1] # basketFile="data/basket/20240307_4017567/data/ukb678162.RData"
kinshipFile = args[2] # kinshipFile="data/genetics/2024/meta/ukb42032_rel_s487950.dat"
panFile= args[3] # panFile="data/basket/20210327_2008261/data/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
bridgeFile = args[4] # bridgeFile="data/basket/20210327_2008261/data/ukb42032bridge31063.txt"
targetDir = args[5] # targetDir="results/sample/"

logInfo = paste0('\n--- Sample Filtering: Settings ---',
               '\nbasketFile: ', basketFile,
               '\nkinshipFile: ', kinshipFile,
               '\npanFile: ', panFile,
               '\nbridgeFile: ', bridgeFile,
               '\ntargetDir: ', targetDir,'\n')
message(logInfo)

# load required packages
for (pkg in c('dplyr','igraph','data.table')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# -------------------------------------------
# --- Filter data released until Feb 2024 ---
# -------------------------------------------
message('\nSelecting individuals.')

# load basket data and CAT12 statistics of preprocessed T1 files
message(sprintf(' - loading basket file %s', basketFile))
load(basketFile)
message(sprintf(' - %d columns and %d rows.', dim(bd)[1],dim(bd)[2]))

# apply filter criteria
# - reported and genetic sex mismatch
# - sex aneuploidy
# - outliers in heterozygosity and missingness rates
# - kinship information available
bd2024 = bd
message(' - applying filter criteria.')
sex = as.numeric(bd2024$f.31.0.0)==as.numeric(bd2024$f.22001.0.0)
sex[is.na(sex)] = FALSE
noAneuploidy = is.na(as.numeric(bd2024$f.22019.0.0))
noHetMiss = is.na(as.numeric(bd2024$f.22027.0.0))
kinship = rep(TRUE, nrow(bd2024)); kinship[is.na(bd2024$f.22021)] = FALSE; kinship[as.numeric(bd2024$f.22021)==1] = FALSE
bd2024 = bd2024[sex & noAneuploidy & noHetMiss & kinship,]
message(sprintf(' - %d cases remaining.', dim(bd2024)[1]))

# define function to select unrelated individuals
select.unrelated <- function(kinshipTable, iid, seed, keep = NULL) {

  # set random number seed
  set.seed(seed)

  # shrink kinshipTable to pairs of individuals included in sample of interest
  kinshipTable <- kinshipTable[(kinshipTable$ID1 %in% iid) & (kinshipTable$ID2 %in% iid),]

  # identify dyadic pairs
  message(' - identify and remove dyadic pair kinships.')
  cat <- c(kinshipTable$ID1,kinshipTable$ID2)
  nodeSubjects <- cat[duplicated(cat)]
  dyadic = kinshipTable[!(kinshipTable$ID1 %in% nodeSubjects) & !(kinshipTable$ID2 %in% nodeSubjects), ]

  # remove one out of two subjects in dyadic pair
  exclude = sample(1:2, nrow(dyadic), replace=T)
  dyadicExclude = c()
  for (i in 1:nrow(dyadic)) {
    if (!is.null(keep) & sum(dyadic[i,c('ID1','ID2')] %in% keep) > 0) { dyadicExclude = c(dyadicExclude, dyadic[i,which(!(dyadic[i,c('ID1','ID2')] %in% keep))]) }
    else { dyadicExclude = c(dyadicExclude, dyadic[i,exclude[i]]) }
  }

  # identify subjects with multiple relationships
  multi = kinshipTable[kinshipTable$ID1 %in% nodeSubjects | kinshipTable$ID2 %in% nodeSubjects, ]
  multi$ID1 = factor(multi$ID1)
  multi$ID2 = factor(multi$ID2)
  adj = get.adjacency(graph.edgelist(as.matrix(multi[,c("ID1","ID2")]), directed = FALSE))

  # Identify trios
  message(' - identify and remove trio kinships.')
  trios = c(NULL) # identify trios
  degrees = c(NULL)
  pb = txtProgressBar(min = 1, max = nrow(adj), style = 3)
  for (i in 1:nrow(adj)) {
  setTxtProgressBar(pb, i)
    trio = c(row.names(adj)[i])
    degree = c(sum(adj[,i]))
    
    for (j in which(adj[,i]==1)) {
      if (!(row.names(adj)[j] %in% trio) && !(row.names(adj)[j] %in% trios)) {
        trio = c(trio, row.names(adj)[j])
        degree = c(degree, sum(adj[,j]))
        if (length(trio) > 3) { break }
        
        for (k in which(adj[,j]==1)) {
          if (!(row.names(adj)[k] %in% trio) && !(row.names(adj)[k] %in% trios)) {
            trio = c(trio, row.names(adj)[k])
            degree = c(degree, sum(adj[,k]))
            if (length(trio) > 3) { break }

            for (l in which(adj[,k]==1)) {
              if (!(row.names(adj)[l] %in% trio) && !(row.names(adj)[l] %in% trios)) {
                trio = c(trio, row.names(adj)[l])
                degree = c(degree, sum(adj[,l]))
                if (length(trio) > 3) { break }
              }
            }
          }
        }
      }
    }
    if (length(trio) == 3) {
      trios = rbind(trios, trio)
      degrees = rbind(degrees, degree)
    }
  }

  # remove "hub" subject from type 1 trios (A related with B, B related with C) 
  triosType1 = matrix(trios[rowSums(degrees) == 4,], ncol = 3)
  degreesType1 = matrix(degrees[rowSums(degrees) == 4,], ncol = 3)
  triosType1Exclude = c()
  for (i in 1:nrow(triosType1)) {
    if (!is.null(keep) & sum(triosType1[i,] %in% keep) > 0 & identical(as.numeric(degreesType1[i,which(triosType1[i,] %in% keep)]),2)) { 
        triosType1Exclude = c(triosType1Exclude, triosType1[i,which(!(triosType1[i,] %in% keep))]) }
    else {
        triosType1Exclude = c(triosType1Exclude, triosType1[i, degreesType1[i,]==2]) }
  }

  # randomly remove 2 subjects from type 2 trios (a, b and c related to each other)
  triosType2 = matrix(trios[rowSums(degrees) == 6,], ncol = 3)
  select = sample(1:3, nrow(triosType2), replace=T)
  triosType2Exclude = c()
  for (i in 1:nrow(triosType2)) {
    if (!is.null(keep) & sum(triosType2[i,] %in% keep) > 0) { triosType2Exclude = c(triosType2Exclude, triosType2[i,which(!(triosType2[i] %in% keep))]) }
    else { triosType2Exclude = c(triosType2Exclude, triosType2[i, -select[i]]) }
  }

  # get list of remaining subjects
  kinshipExclude = c(dyadicExclude, triosType1Exclude, triosType2Exclude)
  iid = iid[!(iid %in% kinshipExclude)]

  # identify remaining kinships
  message('\n - identify and remove remaining kinships.')
  kinshipTable <- kinshipTable[(kinshipTable$ID1 %in% iid) & (kinshipTable$ID2 %in% iid),]
  cat = c(kinshipTable$ID1,kinshipTable$ID2)
  multi = kinshipTable[kinshipTable$ID1 %in% cat | kinshipTable$ID2 %in% cat, ]
  multi$ID1 = factor(multi$ID1)
  multi$ID2 = factor(multi$ID2)
  adj = get.adjacency(graph.edgelist(as.matrix(multi[,c("ID1","ID2")]), directed = FALSE))

  # iteratively remove subject with highest degree until no kinship remains
    # a) remove subjects related to 'keep' subjects
    if (!is.null(keep) & sum(row.names(adj) %in% keep) > 0) {
      idx = c()
      for (i in which(row.names(adj) %in% keep)) {
        idx = c(idx, as.numeric(which(adj[i,] == 1)))
      }
      idx = unique(idx)
      kinshipExclude = c(kinshipExclude,row.names(adj)[idx])
      adj = adj[-idx,-idx]
    }
    # b) remove remaining relationships
    edges = sum(adj)
    pb = txtProgressBar(min = 1, max = edges, style = 3)
    while (sum(adj) > 0) {
      setTxtProgressBar(pb, edges-sum(adj))
      idx = sample.int(nrow(adj),nrow(adj)) # randomly shuffle participants in case of ties
      adj = adj[idx,idx]
      subExclude = row.names(adj)[Matrix::rowSums(adj) %>% order(decreasing = T)][1]
      kinshipExclude = c(kinshipExclude,subExclude)
      idx = which(row.names(adj) %in% subExclude)
      adj = adj[-idx,-idx]
    }

  # return list of individuals to keep
  message('\n')
  iid = iid[!(iid %in% kinshipExclude)]
  return(iid)
}

# apply function and get list of unrelated individuals
message(sprintf(' - loading kinship file %s', kinshipFile))
kinshipTable = read.table(kinshipFile, head=TRUE)
unrelated = select.unrelated(kinshipTable,bd2024$f.eid,86609)
bd = bd[bd$f.eid %in% unrelated,]

# Sanity check: Test for two subjects in a row in relatedness table
if (sum(kinshipTable$ID1 %in% bd$f.eid & kinshipTable$ID2 %in% bd$f.eid) == 0) {
  message(sprintf(' - relatedness successfully removed\n - %d cases remaining.', dim(bd)[1]))
} else {
  message(' - relatedness has not been removed, exiting script.')
  quit()
}

# =============================
# === add pan ancestry data ===
# =============================
message('\nAdding pan ancestry data.')

# load dataset
message(sprintf(' - loading %s',panFile))
pan = read.delim(panFile, sep = '\t', header = TRUE)
message(sprintf(' - loading %s',bridgeFile))
bridge = read.delim(bridgeFile, sep = ' ', header = FALSE)

# prepare pan-ancestry data
# remove duplicates with related flag 'true'
# add pan to dataset
dupids = pan$s[duplicated(pan$s)]
pan = pan[!(pan$s %in% dupids & pan$related == 'true'),]
names(bridge) = c('f.eid', 's')
pan = inner_join(pan, bridge, 's')
bd = left_join(bd, pan, 'f.eid')

# ==============================================================
# === Create variable file and ancestry-stratified iid files ===
# ==============================================================
message('\nCreating variable file and ancestry-stratified iid files.')

# get caucasian ancestry and pan-ancestry population data
message(' - preparing output variables.')
caucasian = as.numeric(bd$f.22006)
caucasian[is.na(caucasian)] = 0
pan = bd$pop

# discovery or replication?
bd$discovery = NA
bd$discovery[caucasian==1] = 1
bd$discovery[caucasian==0 & !is.na(pan)] = 0
bd = bd[!is.na(bd$discovery),]
message(sprintf(' - discovery n = %d | replication n = %d',sum(bd$discovery==1),sum(bd$discovery==0)))

# calculate exact age
birth_year_month = paste(as.numeric(bd$f.34.0.0),as.numeric(bd$f.52.0.0),"01", sep="-")
age = as.numeric((bd$f.53.0.0 - as.Date(birth_year_month,"%Y-%m-%d"))/365.24219052)

# prepare assessment center
ac.f = factor(bd$f.54.0.0)
ac.dummies = model.matrix(~ac.f)
ac = matrix(NA, nrow = length(ac.f), ncol = ncol(ac.dummies)-1)
ac[!is.na(ac.f),] = ac.dummies[,c(2:ncol(ac.dummies))]
ac = data.frame(ac)
names(ac) = paste0('ac',1:ncol(ac))
ac.dummy = ac

# genotyping array dichotom
array = bd$f.22000.0.0
array[array>0] = 1
array[array<0] = 0

# Covs data.frame
df = data.frame(
          FID = bd$f.eid,
          IID = bd$f.eid,
          wba = bd$discovery,
          pan = bd$pop,
          sex = as.numeric(bd$f.31.0.0),
          age = age,
          age2 = age^2,
          ageXsex = age*as.numeric(bd$f.31.0.0),
          age2Xsex = age^2*as.numeric(bd$f.31.0.0),
          headScale = bd$f.25000.2.0,
          x = bd$f.25756.2.0,
          y = bd$f.25757.2.0, 
          z = bd$f.25758.2.0,
          rfMRI = bd$f.25741.2.0, 
          tfMRI = bd$f.25742.2.0,
          array = array,
          ac.dummy) %>%
      cbind(data.frame(
          PC1 = bd$f.22009.0.1, PC2 = bd$f.22009.0.2, PC3 = bd$f.22009.0.3, PC4 = bd$f.22009.0.4, PC5 = bd$f.22009.0.5,
          PC6 = bd$f.22009.0.6, PC7 = bd$f.22009.0.7, PC8 = bd$f.22009.0.8, PC9 = bd$f.22009.0.9, PC10 = bd$f.22009.0.10,
          PC11 = bd$f.22009.0.11, PC12 = bd$f.22009.0.12, PC13 = bd$f.22009.0.13, PC14 = bd$f.22009.0.14, PC15 = bd$f.22009.0.15,
          PC16 = bd$f.22009.0.16, PC17 = bd$f.22009.0.17, PC18 = bd$f.22009.0.18, PC19 = bd$f.22009.0.19, PC20 = bd$f.22009.0.20,
          PC21 = bd$f.22009.0.21, PC22 = bd$f.22009.0.22, PC23 = bd$f.22009.0.23, PC24 = bd$f.22009.0.24, PC25 = bd$f.22009.0.25,
          PC26 = bd$f.22009.0.26, PC27 = bd$f.22009.0.27, PC28 = bd$f.22009.0.28, PC29 = bd$f.22009.0.29, PC30 = bd$f.22009.0.30,
          PC31 = bd$f.22009.0.31, PC32 = bd$f.22009.0.32, PC33 = bd$f.22009.0.33, PC34 = bd$f.22009.0.34, PC35 = bd$f.22009.0.35,
          PC36 = bd$f.22009.0.36, PC37 = bd$f.22009.0.37, PC38 = bd$f.22009.0.38, PC39 = bd$f.22009.0.39, PC40 = bd$f.22009.0.40,
          PanC1 = bd$PC1, PanC2 = bd$PC2, PanC3 = bd$PC3, PanC4 = bd$PC4, PanC5 = bd$PC5, 
          PanC6 = bd$PC6, PanC7 = bd$PC7, PanC8 = bd$PC8, PanC9 = bd$PC9, PanC10 = bd$PC10,
          PanC11 = bd$PC11, PanC12 = bd$PC12, PanC13 = bd$PC13, PanC14 = bd$PC14, PanC15 = bd$PC15,
          PanC16 = bd$PC16, PanC17 = bd$PC17, PanC18 = bd$PC18, PanC19 = bd$PC19, PanC20 = bd$PC20))

# write subject IIDs and vars
message(' - writing output files.')
system(sprintf('mkdir -p %s',targetDir))
data.table::fwrite(df, file = sprintf('%s/r2024.vars.gz',targetDir), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, compress = 'gzip')
data.table::fwrite(df[df$wba == 1,] , file = sprintf('%s/r2024.vars.discovery.gz',targetDir), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, compress = 'gzip')
data.table::fwrite(df[df$wba == 0,] , file = sprintf('%s/r2024.vars.replication.gz',targetDir), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, compress = 'gzip')
write.table(data.frame(FID = df$FID[df$wba==1], IID = df$IID[df$wba==1]), file = sprintf('%s/iid.discovery.txt',targetDir), quote = F, sep = '\t', row.names = F)
write.table(data.frame(FID = df$FID[df$wba==0], IID = df$IID[df$wba==0]), file = sprintf('%s/iid.replication.txt',targetDir), quote = F, sep = '\t', row.names = F)

# write subject IIDs for each ancestry (replication sample only)
repl = df[df$wba==0,]
for (anc in names(table(repl$pan))) {
  tmp = data.frame(FID = repl$FID[repl$pan==anc & !is.na(repl$pan)], IID = repl$IID[repl$pan==anc & !is.na(repl$pan)])
  write.table(tmp, file = sprintf('%s/iid.replication.%s.txt',targetDir,anc), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}
system(sprintf('chmod 750 %s/iid*',targetDir))
system(sprintf('chmod 750 %s/r2024*',targetDir))
message('\n--- Completed: Sample Filtering ---')

