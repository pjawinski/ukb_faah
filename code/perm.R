#!/usr/bin/env Rscript

# ==================================================
# === Run permutation test for rank observations ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop(paste0('expected 8 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inFile = args[1] # inFile="results/combined/rg.info.txt"
outFile = args[2] # outFile="results/combined/ranktest"
vals1 = args[3] # vals1="subj_sd_dashti_2019_p,subj_chron_jones_2019_p,subj_ins_jansen_2019_p,subj_ds_wang_2019_p,subj_dn_dashti_2021_p"
vals2 = args[4] # vals2="obj_duration_jones_2019b_p,obj_timing_jones_2019b_p,obj_quality_jones_2019b_p,obj_efficiency_jones_2019b_p,obj_inactivity_jones_2019b_p"
vals1label = args[5] # vals1label='subjective'
vals2label = args[6] # vals2label='objective'
nperms = as.numeric(args[7]) # nperms=10000
absolute = args[8] # absolute=TRUE
  
logInfo = paste0('\n--- Run permutation test for rank observations ---',
               '\ninFile: ', inFile,
               '\noutFile: ', outFile,
               '\nvals1: ', vals1,
               '\nvals2: ', vals2,
               '\nvals1label: ', vals1label,
               '\nvals2label: ', vals2label,
               '\nnperms: ', nperms,'\n')
message(logInfo)

# attach required packages
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
vals1 = str_split(vals1, ',')[[1]]
vals2 = str_split(vals2, ',')[[1]]

# load data
df = read.delim(inFile, header = T, sep = '\t')

# merge values of vector 1
for (i in 1:length(vals1)) {
  tmp = data.frame(label = vals1label, vals = df[,vals1[i]])
  if (i == 1) { v1 = tmp } else { v1 = rbind(v1,tmp) }
}

# merge values of vector 2
for (i in 1:length(vals2)) {
  tmp = data.frame(label = vals2label, vals = df[,vals2[i]])
  if (i == 1) { v2 = tmp } else { v2 = rbind(v2,tmp) }
}

# combine lists of vectors and rank them
obs = rbind(v1,v2)
if (absolute=="TRUE") { 
  message(' - considering absolute values for calculations.')
  obs$vals = abs(obs$vals)
}
obs = obs[order(obs$vals),]
obs$rank = rank(obs$vals,ties.method = 'average')

# calculate observed mean ranks
obs.mean.v1 = mean(obs$rank[obs$label == vals1label]) 
obs.mean.v2 = mean(obs$rank[obs$label == vals2label]) 
obs.diff = obs.mean.v2 - obs.mean.v1
logInfo = paste0(logInfo,
                 '\n--- calculations ---',
                 '\nobserved meanRank v1: ', obs.mean.v1,
                 '\nobserved meanRank v2: ', obs.mean.v2,
                 '\nobserved difference: ', obs.diff)

# now randomly shuffle ranks
set.seed(1234)
message(' - run permutations.')
exp = data.frame(matrix(nrow = nperms,ncol = 3))
names(exp) = c('mean.v1','mean.v2','diff')
pb = txtProgressBar(min = 1, max = nperms, style = 3); k = 0
for (i in 1:nperms) {
  setTxtProgressBar(pb, i)
  tmp = obs
  tmp$rank = sample(x = tmp$rank, size = nrow(tmp), replace = F)
  exp[i,1] = mean(tmp$rank[obs$label == vals1label])
  exp[i,2] = mean(tmp$rank[obs$label == vals2label])
  exp[i,3] = exp[i,2] - exp[i,1] 
}

# calculate p-values
logInfo = paste0(logInfo,
                 '\nnperms: ', nperms,
                 '\np-value for observed difference (two-tailed): ', sum(abs(exp$diff) > abs(obs.diff))/nperms,
                 '\np-value for v1<v2: ', sum(exp$diff > obs.diff)/nperms,
                 '\np-value for v1>v2: ', sum(exp$diff < obs.diff)/nperms,
                 '\np-value for v1 (higher rank than expected): ', sum(exp$mean.v1 > obs.mean.v1)/nperms,
                 '\np-value for v1 (lower rank than expected): ', sum(exp$mean.v1 < obs.mean.v1)/nperms,
                 '\np-value for v2 (higher rank than expected): ', sum(exp$mean.v2 > obs.mean.v2)/nperms,
                 '\np-value for v2 (lower rank than expected): ', sum(exp$mean.v2 < obs.mean.v2)/nperms)

# save results
message(sprintf('\n - writing %s.log',outFile))
sink(sprintf('%s.log',outFile)) 
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s.log', outFile))
message(sprintf(' - writing %s.exp.txt',outFile))
write.table(exp, file = sprintf('%s.exp.txt',outFile), sep = '\t', col.names = T, row.names = F, quote = F)
message(paste0('--- Completed: Run permutation test for rank observations ---\n'))
