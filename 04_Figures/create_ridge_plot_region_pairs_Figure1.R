# Author: Gang Chen
# modified: Jonas Steinhäuser

# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(data.table)
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(scales)
library(brms)

# Select R-Workspace file from output folder of running MBA analysis
# alternatively, select R-Workspace from ./03_MBA_analysis/exemplary_output folder 
workspace_path <- file.choose()
load(workspace_path)


# this is the actual plotting function that will be used to create the ridge plot
stackPlot <- function(dat, xlim, labx) {
  data <- data.frame(dat)
  data$X <- NULL
  nobj=dim(data)[1]
  
  # rename columns with ROI list
  # comment out lines 27-32 to create the plot with the region names used in the dataframe
  dimnames(dat)[[2]] <- c("amygdala-AI", "dACC-AI", "dACC-amygdala", "dmPFC-AI", "dmPFC-amygdala", "dmPFC-dACC",
                          "PCC-AI", "PCC-amygdala", "PCC-dACC", "PCC-dmPFC", "PMI-AI", "PMI-amygdala",
                          "PMI-dACC", "PMI-dmPFC", "PMI-PCC", "TP-AI",  "TP-amygdala",  "TP-dACC", "TP-dmPFC",
                          "TP-PCC", "TP-PMI", "vmPFC-AI", "vmPFC-amygdala", "vmPFC-dACC", "vmPFC-dmPFC",
                          "vmPFC-PCC", "vmPFC-PMI", "vmPFC-TP", "dlPFC-AI", "dlPFC-amygdala", "dlPFC-dACC",
                          "dlPFC-dmPFC", "dlPFC-PCC", "dlPFC-PMI", "dlPFC-TP", "dlPFC-vmPFC")
  rois <- dimnames(dat)[[2]]
  colnames(data) <- rois
  data_stats <- data.frame(1:length(rois))
  
  # create ROI column instead of numerics to match threat table above
  data_stats$ROI <- rois
  data_stats$mean <- colMeans(data)  # median: quantile(x, probs=.5)
  data_stats$P <- colSums(data > 0)/nobj
  data_stats$Pn <- ifelse(data_stats$P < .5, 1-data_stats$P, data_stats$P)
  # this will order the distributions correctly
  data_stats <- data_stats[order(data_stats$mean),]
  
  data_trans <- as.data.frame(t(as.matrix(data)))
  # add two more columns
  data_trans <- tibble::rownames_to_column(data_trans, "ROI")
  data_trans$X <- 1:nrow(data_trans)
  
  # merge values & stats into one table by ROI
  data_merge <- merge(data_stats, data_trans, by = "ROI")
  data_merge <- data_merge[order(data_merge$X),]
  #browser()
  # Transform data into long form: Melt dataframe by ROI
  data_long <- reshape2::melt(data_trans, id=c("ROI","X"))
  data_long <- data_long[order(data_long$X),]
  
  #clunky, but for now stats by ensuring orders are all the same and repeating each value nobj times. no success for alternatives. 
  data_long$mean <- rep(data_merge$mean, each = nobj)
  data_long$P <- rep(data_merge$P, each =nobj)
  data_long$Pn <- rep(data_merge$Pn, each =nobj)
  data_long$gray.vals <- rep(data_merge$gray.vals, each =nobj)
  
  ########################  G R A P H I N G  #########################################################
  #######################     V A R I A B L E S    ######################################################
  # set your labels here so you don't have to change within the plot below: 
  y.axis.labs <- data_stats$ROI                              # y axis labels
  sec.y.axis.labs <- round(data_stats$P,2)                   # second y axis labels (probabilities) - Rounded to 2 decimals
  
  ################# X AXIS LABELS ###########################################################
  # X AXIS LABELS NEED TO CHANGE TO CORRESPOND TO DATA SET! UNCOMMENT WHICHEVER MATCHES
  x.axis.labs <- NULL                                # x axis labels  INTERACTION, not sure what to put.
  x.labs.pos <- NULL                                 # a axis position INTERACTION, change when labels decided
  
  ######## T I T L E S #############################################################
  #data.name <- tl
  graph.title <- "Interaction (% signal change)"                                   # graph title 
  legend.title <- "P+"                              # legend title
  y.axis.title <- NULL                                       # for now ...
  x.axis.title <- NULL                                       # for now...
  
  ########################## D A T A  ##############################################################
  # GRAPH DATA 
  dataset <- data_long
  x.values <- data_long$value                               # x values
  y.values <- data_long$ROI                                 # y values
  y.values.RO <- data_long$value                            # values to reorder Y by
  distrib.fill <- data_long$P                       # fill graph with probabilities
  group <- data_long$ROI
  
  ######################### S A V E  ################################################
  # SAVE SETTINGS -- Currently low res and jpeg to make it easier to share
  # adjusting height + width will throw font sizes out of wack: need change (see other aspects below) 
  
  dpi <- 300
  units <- "in"                                           # "in", "cm", or "mm"
  height <- 5
  width <- 9
  file.type <- ".pdf"                   # can be ".jpeg",".pdf",".png",".bmp",".tiff",etc
  
  ############################### O T H E R  #################################################
  #gradient.colors<-c("#41245C","yellow","gray","gray","blue","#C9182B") # change gradient colors  
  gradient.colors <- c("blue","cyan","gray","gray","yellow","#C9182B")  # change gradient colors here 
  ROI.label.size <- 13                 # adjust ROI and probability y-axis font size
  P.label.size <- 13
  title.size <- 20                         # adjust graph title size 
  x.axis.size <- 15                                        # adjust x-axis label sizes
  
  ##################  G R A P H  ########################################
  # change information about the graph and add other characteristics using ggplot and ggridges
  ggplot(dataset, aes(x = x.values, 
                      y = as.numeric(reorder(y.values, y.values.RO)), 
                      fill = distrib.fill, 
                      group = group))   +
    guides(fill = guide_colorbar(barwidth = 1,             #legend characteristics
                                 barheight = 20,
                                 nbin = 100, # can change # bins to play around with color gradient
                                 frame.colour = "black",
                                 frame.linewidth = 1.5,
                                 ticks.colour = "black",
                                 title.position = "top",
                                 title.hjust = 0.5)) +
    #geom_density_ridges() +                            # scale = spacing, alpha = transparency
    stat_density_ridges(quantile_lines = TRUE,         # divide into two quantiles (show mean)
                        quantiles = 2,
                        size = .6,
                        alpha = .8,
                        scale = 2,
                        color = "black") +
    geom_vline(xintercept = 0,                        #create line at X = 0
               linetype="solid",
               alpha = 1,
               size = 1,
               color = "green") +
    scale_fill_gradientn(
      colors = gradient.colors,                       # set gradient
      limits = c(0,1),                                # scale size
      name = legend.title,
      breaks = c(0,0.05,0.1,0.9,0.95,1),
      expand = expand_scale(0),
      labels = c("0","0.05","0.1","0.9", "0.95","1")
    ) +
    scale_y_continuous(breaks = 1:length(rois), # A VERY HACK-Y WAY TO HAVE TWO Y AXES W DISCRETE DATA
                       labels = y.axis.labs,   # Trick ggplot into thinking data is continuous...
                       sec.axis = sec_axis(~.,  # Second axis to show probabilities
                                           breaks = 1:length(rois),
                                           labels = sec.y.axis.labs)) +
    theme_ridges(font_size = ROI.label.size, grid = TRUE, center_axis_labels = TRUE) +  # theme info
    #ggtitle(graph.title)+                                                   # graph title
    annotate("text", x=0.06, y=1.5, label=labx, size=5.5)+
    theme(
      plot.title = element_text(vjust = -0.5, size = title.size),   # plot title size and position
      #axis.title.x=element_text(vjust=0, hjust=0.5),
      axis.text.y.left = element_text(size=ROI.label.size),        # y-axis text size
      axis.text.y.right = element_text(size = P.label.size),       # y-axis info for right axis
      axis.text.x = element_text(size = x.axis.size),   # x-axis text size/face
      #axis.text.x = element_text(size = x.axis.size, face = "bold"),   # x-axis text size/face
      legend.title.align = 5,
      #legend.text = element_text(face = "bold"),
      legend.title = element_text(size = 15))+
    #legend.title = element_text(face = "bold", size = 15))+
    labs(
      #x = 'interaction (% signal change)',                 # Add or not add X and Y labels
      x = NULL,
      y = NULL) +
    scale_x_continuous(labels=waiver(), limits = xlim)
  #scale_x_continuous(breaks = x.labs.pos, labels = c(x.axis.labs),limits = xlim)
  # ggsave(file = paste(period, "_", data.name, file.type, sep=""), width=width, height=height, units = units, dpi = dpi)
}


# set up data structure to create ridge plot 
ns <- lop$iterations*lop$chains/2
aa <- fixef(fm, summary = FALSE)
bb <- ranef(fm, summary = FALSE)
rL <- dimnames(bb$mmROI1ROI2)[[2]]
cnt <- function(x, ns) return(sum(x>0)/ns)

ww <- function(aa, bb, tm, nR) {
  ps0 <- array(apply(bb[['mmROI1ROI2']][,,tm], 2, "+", bb[['mmROI1ROI2']][,,tm]), c(ns, nR, nR))
  ps <- apply(ps0, c(2,3), '+', aa[,tm])
  dimnames(ps) <- list(1:ns, rL, rL)
  tmp <- ps
  sel1 <- match(dimnames(bb$`ROI1:ROI2`)[[2]], outer(dimnames(ps)[[2]],dimnames(ps)[[3]], function(x,y) paste(x,y,sep="_")))
  sel2 <- match(dimnames(bb$`ROI1:ROI2`)[[2]], outer(dimnames(ps)[[2]],dimnames(ps)[[3]], function(x,y) paste(y,x,sep="_")))
  ad <- function(tt,bb,s1,s2) {tt[s1] <- tt[s1] + bb; tt[s2] <- tt[s2] + bb; return(tt)}
  for(ii in 1:ns) tmp[ii,,] <- ad(tmp[ii,,], bb$`ROI1:ROI2`[ii,,tm], sel1, sel2)
  ps <- tmp
  return(ps)
}

vv <- function(ps, ns, nR) {
  mm <- apply(ps, c(2,3), mean)
  for(ii in 1:nR) for(jj in 1:nR) ps[,ii,jj] <- sqrt(2)*(ps[,ii,jj] - mm[ii,jj]) + mm[ii,jj]
  return(ps)
}

lvl <- levels(lop$dataTable[[lop$EOIc]])  # levels
nl <- nlevels(lop$dataTable[[lop$EOIc]])  # number of levels: last level is the reference in deviation coding
ps <- array(0, dim=c(nl, ns, nR, nR)) # posterior samples
for(jj in 1:(nl-1)) ps[jj,,,] <- ww(aa, bb, paste0(lop$EOIc[ii],jj), nR)
ps[nl,,,] <- ww(aa, bb, 'Intercept', nR)
psa <- array(0, dim=c(nl, ns, nR, nR)) # posterior samples adjusted
for(jj in 1:(nl-1)) {
  psa[jj,,,] <- ps[nl,,,] + ps[jj,,,]
  psa[nl,,,] <- psa[nl,,,] + ps[jj,,,]
}
psa[nl,,,] <- ps[nl,,,] - psa[nl,,,]  # reference level
dimnames(psa)[[3]] <- rL
dimnames(psa)[[4]] <- rL

for(jj in 1:(nl-1)) for(kk in (jj+1):nl)
  xx <- vv(psa[jj,,,] - psa[kk,,,], ns, nR)
dimnames(xx)[[1]] <- 1:ns
nE <- nR*(nR-1)/2
ef <- NULL
ll <- 0
for(ii in 2:dim(xx)[2]) for(jj in 1:(ii-1)) {
  ll <- ll + 1
  ef <- cbind(ef, xx[,ii,jj])
  dimnames(ef)[[2]][ll] <- paste0(dimnames(xx)[[2]][ii], '_', dimnames(xx)[[3]][jj])
}

dev.new(width=10, height=20) # set the plot window size
stackPlot(ef, c(-0.1, 0.26), labx='') # set your x-axis range: c(-0.5, 0.5)

# save as raw pdf figure
ggsave("Figure1_ridge_plot_region_pairs_raw.pdf", device = "pdf", dpi = "print")