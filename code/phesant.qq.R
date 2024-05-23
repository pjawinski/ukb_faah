#!/usr/bin/env Rscript

# =======================================
# === draw qq plots of phewas results ===
# =======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=13) {
  stop(paste0('expected 13 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inputFile = args[1] # inputFile="results/combined/phesant.summary.txt"
outFile = args[2] # outFile="results/combined/phewas.qq.png"
categoryCol = args[3] # categoryCol="custom_category"
categories = args[4] # categories="Blood and urine assays,Cancer register,Cognitive function,Diet by 24-hour recall,Digestive health,Family history and early life factors,Hospital Inpatient - Administration,Imaging,Lifestyle and environment,Maternity and sex-specific factors,Medical history and conditions,Medications,Mental health,Operations,Physical activity,Physical measures,Sociodemographics,Work environment"
categoriesRename = args[5] # categoriesRename="Blood and urine,Cancer register,Cognitive function,Diet,Digestive health,Family history,Hospital Inpatient,Imaging,Lifestyle,Maternity,Medical history,Medications,Mental health,Operations,Physical activity,Physical measures,Sociodemographics,Work environment"
pCol = args[6] # pCol="pvalue"
ylim = as.numeric(args[7]) # upper y axis limit; ylim = 6
ysteps = as.numeric(args[8]) # y axis breaks; ysteps = 2
xlim = as.numeric(args[9]) # upper y axis limit; xlim = 4.5
xsteps = as.numeric(args[10]) # y axis breaks; xsteps = 2
patchCols = as.numeric(args[11]) # patchwork columns; patchCols = 6
width = as.numeric(args[12]) # plot width; width = 7
height = as.numeric(args[13]) # plot height, height = 4.5

message(paste0('\n--- Create qq plots of phewas results ---',
               '\ninputFile: ', inputFile,
               '\noutFile: ', outFile,
               '\ncategoryCol: ', categoryCol,
               '\ncategories: ', categories,
               '\ncategoriesRename: ', categoriesRename,
               '\npCol: ', pCol,
               '\nylim: ', ylim,
               '\nysteps ', ysteps,
               '\nxlim: ', xlim,
               '\nxsteps ', xsteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr','ggplot2','patchwork','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
categories = str_split(categories, ',')[[1]]
categoriesRename = str_split(categoriesRename, ',')[[1]]

# load data
df = read.delim(inputFile, header = T, sep = '\t')

# make plot
df$p_log10 = -log10(df[,pCol])

# calculate FDR < 0.05
df$fdr = p.adjust(df[,pCol], method = 'BH')
df$fdrLogical = df$fdr < 0.05

# get categories
df = df[df[[categoryCol]] %in% categories,]
for (i in 1:length(categoriesRename)) {
  df[[categoryCol]][df[[categoryCol]]==categories[i]] = categoriesRename[i]
}
df[[categoryCol]] = factor(df[[categoryCol]],
                            levels = categoriesRename,
                            labels = categoriesRename)

# --------------------
# --- draw qq plot ---
# --------------------

# create new data frame for plot
dfplot = df
dfplot = dfplot[order(dfplot[[categoryCol]]),]
ncategories = length(unique(dfplot$custom_category))

# define qq plot function
qq.plot = function(pvalues, qq_title, xlim, xsteps, ylim, ysteps, type) {
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  
  # calculated expected values, transform to -log10 scale
  n = length(pvalues)+1
  exp.x = -log10((rank(pvalues, ties.method="first")-.5)/n) # exp.x2 = -log10((rank(pvalues2, ties.method="first")-.5)/n)
  pvalues = -log10(pvalues) # pvalues2 = -log10(pvalues2)
  
  #this is a helper function to draw the confidence interval
  conf.points=1000
  conf.alpha = 0.05
  conf.points = min(conf.points, n-1);
  mpts<-matrix(nrow=conf.points*2, ncol=2)
  for(i in seq(from=1, to=conf.points)) {
    mpts[i,1]<- -log10((i-.5)/n)
    mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
    mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
    mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
  }
  
  # define dataframes for ggplot
  qqdata = data.frame(pvalues, exp.x) # qqdata = data.frame(pvalues, pvalues2, exp.x, exp.x2)
  confdata = data.frame(x = mpts[,1], y = mpts[,2])
  
  # make ggplot
  if (type == "xaxis") {
    qqtheme = theme(plot.margin=margin(2,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = -2, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_blank(), # axis.title.y = element_text(size=8,margin = margin(t = 0, r = 1, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(linewidth = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  } else if
  (type == "yaxis") {
    qqtheme = theme(plot.margin=margin(2,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = -2, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_blank(), # axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_text(size=8,margin = margin(t = 0, r = 5, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(linewidth= 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  } else if
  (type == "both") {
    qqtheme = theme(plot.margin=margin(2,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = -2, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_text(size=8,margin = margin(t = 0, r = 5, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(linewidth = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  } else if
  (type == "none") {
    qqtheme = theme(plot.margin=margin(2,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = -2, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_blank(), # axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_blank(), # axis.title.y = element_text(size=8,margin = margin(t = 0, r = 1, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(linewidth = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  }
  
  ggplot(data=qqdata, aes(x=exp.x, y=pvalues)) +
    ggtitle(qq_title) +
    geom_polygon(data = confdata, aes(x = x, y = y), cex = 1, fill = "grey90", show.legend = FALSE, alpha = 1) +
    geom_line(aes(x = exp.x, y = exp.x)) + 
    geom_point(cex = 0.5, col=c("royalblue3"), show.legend = FALSE, alpha = 1) +
    #geom_text(data = data.frame("NA"), label = "Mental health", x=1.5, y=11, hjust = 0.5, size = 4, fontface='bold') +
    #geom_point(aes(x=exp.x2, y=pvalues2), cex = 0.5, col=c("red"), show.legend = FALSE, alpha = 1) +
    #stat_smooth(aes(x = exp.x, y = exp.x), method="lm",fullrange=TRUE) +
    scale_x_continuous(limits = c(0, xlim), breaks = seq(0,xlim,xsteps), name = expression("expected -log"[10]*"("*italic(p)*")")) + # gsub("\\.", " ",
    scale_y_continuous(limits = c(0, ylim), breaks = seq(0,ylim,ysteps), name = expression("observed -log"[10]*"("*italic(p)*")")) + # gsub("\\.", " ",
    geom_segment(aes(x=0,xend=xlim,y=-Inf,yend=-Inf), colour = "black", linewidth = 0.25) +
    geom_segment(aes(y=0,yend=ylim,x=-Inf,xend=-Inf), colour = "black", linewidth = 0.25) +
    theme_bw() +
    qqtheme
}
  
# create qq plot
nboth=length(categoriesRename)-length(categoriesRename)%%patchCols+1
if (nboth>length(categoriesRename)) { nboth = nboth-patchCols }
for (i in 1:length(categoriesRename)) {
  if (i==nboth) {
    tmp = qq.plot(pvalues = dfplot[[pCol]][dfplot[[categoryCol]]==categoriesRename[i]], qq_title = categoriesRename[i], xlim, xsteps, ylim, ysteps, type = 'both')
  } else if ((i+patchCols-1)%%patchCols==0) {
    tmp = qq.plot(pvalues = dfplot[[pCol]][dfplot[[categoryCol]]==categoriesRename[i]], qq_title = categoriesRename[i], xlim, xsteps, ylim, ysteps, type = 'yaxis')
  } else if (i>nboth) {
    tmp = qq.plot(pvalues = dfplot[[pCol]][dfplot[[categoryCol]]==categoriesRename[i]], qq_title = categoriesRename[i], xlim, xsteps, ylim, ysteps, type = 'xaxis')
  } else {
    tmp = qq.plot(pvalues = dfplot[[pCol]][dfplot[[categoryCol]]==categoriesRename[i]], qq_title = categoriesRename[i], xlim, xsteps, ylim, ysteps, type = 'none')
  }
  assign(sprintf('qq.%d',i),tmp)
}

# merge plots
for (i in 1:length(categoriesRename)) {
  if (i==1) {
    qqplots = qq.1
  } else {
    qqplots = qqplots + get(sprintf('qq.%d',i))
  }
}
if (length(categoriesRename)>patchCols) {
  qqplots = qqplots + plot_layout(ncol = patchCols)
} else {
  qqplots = qqplots + plot_layout(ncol = length(categoriesRename))
}


# save plot  
message(sprintf('Saving %s', outFile))
ggsave(outFile, qqplots, width = width, height = height, units = "in", dpi = 300)
system(sprintf('chmod 770 %s', outFile))
message('-- Completed: Create qq plots of phewas results ---')


