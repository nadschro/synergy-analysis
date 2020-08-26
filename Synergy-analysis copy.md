---
title: "Analysis framework for evaluating synergy driving gene expression"
output: html_document
author: Nadine Schrode
date: 02/28/2020
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#### 1. Create directories
Create a “results” folder in your current working directory. All plots (PDF) and tables (CSV) will be saved in this folder.
Create a “genesets” folder, containing gene set groups of interest in the gmt format.


## Preparing R environment
#### 2.	Name your experiment. 
This will be used to save result files and plot titles. In this example data set, we study the perturbation of 4 genes, in a hiPSC model of schizophrenia.
```{r}
experiment.title="SCZ"
```

#### 3.	Load R packages
```{r}
pacman::p_load(limma, edgeR, pheatmap, RColorBrewer, ggplot2, ggpubr, 
               qvalue,  plyr, wesanderson, GSEABase, grid, scales, WebGestaltR, stringr)
```

#### 4.	Run custom functions

The mds() function is based on plotMDS() in the limma package. When given a DGE object and a column in the meta data table, containing groups of interest, it produces an multidimensional scaling plot, colored by the provided groups.
```{r}
mds <- function(normDGE, metacol, title){
  mcol <- as.factor(metacol)
  col <- rainbow(length(levels(mcol)), 1, 0.8, alpha = 0.5)[mcol]
  plotMDS(normDGE, col = col, pch = 16, cex = 2)
  legend("center", 
         fill = rainbow(length(levels(mcol)), 1, 0.8), 
         legend = levels(mcol), 
         horiz = F, 
         bty = "o", 
         box.col="grey", 
         xpd=TRUE)
  title(main=title)
}
```

The cameraplusplots() function is based on camera() in the limma package. When given a contrast in form of a named vector (such as a column of the contrast matrix), a list of gene set groups, a voom object, a design matrix and a color palette for the gene set groups, it produces a scatter plot of all tested gene sets and their adjusted P-values as well as a bar graph of the 10 most significant gene sets, colored by gene set group/category.

```{r}
cameraplusplots <- function(contrast, genesetlist, vobject, design, catcolors, title){ 
  tmp.list <- list()
  cam <- data.frame(matrix(ncol = 5, nrow = 0))
  for (i in 1:length(genesetlist)){ 
    cam.s <- camera(vobject, genesetlist[[i]], design, contrast = contrast, inter.gene.cor = 0.01) 
    tmp.list[[i]] <- cam.s
    names(tmp.list)[i] <- names(genesetlist)[i] 
    tmp.list[[i]]$category <- names(tmp.list[i]) 
    colnames(cam) <- names(tmp.list[[1]]) 
    cam <- rbind.data.frame(cam, tmp.list[[i]]) 
    print(paste0("Gene set categories run: ", i))
  }
  cam$neglogFDR <- -log10(cam$FDR) 
  ## for plotting purposes only: 
  cam$dirNeglogFDR <- cam$neglogFDR
  cam[(cam$Direction == "Down"), "dirNeglogFDR"] <- -cam[(cam$Direction == "Down"), "neglogFDR"]
  grob <- grobTree(textGrob(c("UP","DOWN"), x = c(0.94, 0.89), y = c(0.95, 0.05), hjust = 0, gp = gpar(fontsize = 13))) 
  q <- ggplot(aes(x = cam$category, y = dirNeglogFDR, color = category), data = cam) +
    scale_color_manual(values = catcolors) +
    geom_jitter(aes(size = NGenes, alpha = neglogFDR), pch = 19, show.legend = F) +
    scale_size_continuous(range = c(4,16)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    geom_hline(yintercept = c(-1.3, 1.3), color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-10, 10), oob = squish, labels = abs) +
    labs(x = "Gene set categories", y = "-log10(FDR)", title = title) +
    theme_bw(14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    annotation_custom(grob) 
  print(q)
  cam$geneSet <- row.names(cam)
  cam10 <- as.data.frame(cam)
  cam10 <- cam10[order(cam10$FDR),]
  cam10 <- cam10[1:10,]
  grob <- grobTree(textGrob(c("DOWN","UP"), x = c(0.03, 0.9),  y=c(0.025), hjust = 0, gp = gpar(fontsize = 9, col = "grey60")))
  g <- ggplot(aes(x = geneSet, y = dirNeglogFDR, fill = category), data = cam10) +
    geom_col()+
    aes(reorder(stringr::str_wrap(geneSet, 60),-FDR), dirNeglogFDR) +
    xlab(NULL) +
    geom_hline(yintercept = c(-1.3, 1.3), color = "red", alpha = 0.3) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-10, 10), oob = squish, labels = abs) +
    labs(y = "-log10(FDR)", title = title) +
    scale_fill_manual(values = catcolors) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    annotation_custom(grob) 
  print(g) 
  return(cam)
}
```

The oraplot() function takes the data frame resulting from over-representation analysis using WebGestaltR and a color palette as input and returns a bar graph of the 10 most significant gene sets.
```{r}
oraplot <- function(orares, catcolors, name){
  orares.n <- orares[order(orares$FDR), ]
  orares.n <- orares.n[1:10, ]
  orares.n$neglogFDR <- -log10(orares.n$FDR)
  orares.n <- orares.n[orares.n$neglogFDR>0, ]
  orares.n$geneSet <- gsub("_", " ", orares.n$geneSet)
  g <- ggplot(aes(x=reorder(str_wrap(geneSet, 60), neglogFDR), y = neglogFDR, fill = category), data = orares.n) +
    geom_col() + 
    geom_hline(yintercept = 1.3, color = "red", alpha = 0.5) +
    labs(y = "-log10(FDR)", x = "", title = paste0(name)) +
    scale_fill_manual(values = catcolors) +
    coord_flip() +
    theme_bw(11)
  return(g)
}
```

The power.compare.logFC() function is given the variances of the combinatorial perturbation and the additive model comparisons (sig1, sig2; variance in the additive model is usually higher, in proportion to the number of individual perturbations). Further, the number of samples used to determine the variances (N), a vector of sample numbers of interest (N_other), a significance cutoff (alpha) and the number of tests performed (n_tests, usually the number of transcripts). N_other/N is the relative sample size and the variance of logFC1 - logFC2 is sig12 + sig22. Since the standard error of the mean is inversely proportional to sqrt(N), multiplying the sample size by F decreases the SE by sqrt(F). On the variance scale, this corresponds to dividing by n_scale.
```{r}
power.compare.logFC <- function( sig1, sig2, N, N_other = c(2,4,6,8,10), alpha = 0.05, n_tests = 20000){
  d <- seq(0, 3, length.out=1000)
  alpha_multiple <- alpha / n_tests
  df <- lapply( N_other/N, function(n_scale){
    sigSq <- (sig1^2 + sig2^2) / n_scale
    cutoff <- qnorm( alpha_multiple/2, 0, sd = sqrt(sigSq), lower.tail = FALSE)
    p1 <- pnorm(-1*cutoff, d, sqrt(sigSq))
    p2 <- 1-pnorm(cutoff, d, sqrt(sigSq))
    data.frame(n_scale, d, power=p1+p2)
  })
  df <- do.call("rbind", df)
  ggplot(df, aes(d, power, color = as.factor(n_scale*N))) + 
    geom_line() + 
    theme_bw(14) + 
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + 
    ylim(0, 1) + 
    scale_color_discrete("Samples") + 
    xlab(bquote(abs(logFC[observed] - logFC[expected]))) + 
    ggtitle("Power versus difference in logFC")
}
```

The categorize.synergy() function is given the combined matrix of log2FC values of the additive model, the combinatorial perturbation and the synergistic effect differential expression results. It creates a new column in the resulting data frame, assigning synergy categories to each gene.
```{r}
categorize.synergy <- function(logFCmatrix, meanSE){
  m <- logFCmatrix
  m$magnitude.syn <- NA
  for (i in 1:length(m$Gene_name)){
    if (m$Synergistic.logFC[i] > meanSE){
      if (m$Additive.logFC[i] < -meanSE){
        if (m$Combinatorial.logFC[i] > meanSE){
          m$magnitude.syn[i] = "more.up"
        } else m$magnitude.syn[i] = "less.down"
      } else m$magnitude.syn[i] = "more.up"
    }
    else if (m$Synergistic.logFC[i] < -meanSE){
      if (m$Additive.logFC[i] > meanSE){
        if (m$Combinatorial.logFC[i] < -meanSE){
          m$magnitude.syn[i] = "more.down"
        } else m$magnitude.syn[i] = "less.up"
      } else m$magnitude.syn[i] = "more.down"
    } else m$magnitude.syn[i] = "same"
  }
  m$magnitude.syn <- as.factor(m$magnitude.syn)
  return(m)
}
```

The stratify.by.syn.cat() function is given a subset of interest of the table created by the categorize.synergy() function. It creates a list object containing vectors of gene symbols by synergy category, which is used as input for over-representation analysis with WebGestaltR.
```{r}
stratify.by.syn.cat <- function(log2FC.matrix.sub){
  synergy.cat.list <- list("less.down" = as.character(log2FC.matrix.sub[
    log2FC.matrix.sub$magnitude.syn == "less.down", "Gene_name"]),
    "less.up" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "less.up", "Gene_name"]),
    "more.down" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "more.down", "Gene_name"]),
    "more.up" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "more.up", "Gene_name"]),
    "same" = as.character(log2FC.matrix.sub[
      log2FC.matrix.sub$magnitude.syn == "same", "Gene_name"]))
  return(synergy.cat.list)
}
```

## Loading and preprocessing data
#### 5.	Load data 
Read in metadata (here "Meta_CombinatorialPerturbationProject.csv") and expression data (raw counts, here "Counts_CombinatorialPerturbationProject.csv") and match the order (this step is very important to avoid later confusion).
```{r}
meta <- read.csv("Meta_CombinatorialPerturbationProject.csv", row.names = 1)
counts <- read.csv("Counts_CombinatorialPerturbationProject.csv", row.names=1)
meta <- meta[match(colnames(counts), row.names(meta)), ]
```

#### 6.	Filter lowly expressed genes
Plot counts over cpm(counts) and visually inspect the graph. Here, we aim to keep roughly 10 counts in 4 or more samples (number dependent on total number of samples/replicates) (Fig. 5A).
```{r}
#pdf(paste0("results/", experiment.title, "-1_cpm-counts.pdf"))
plot(cpm(counts)[, 1], counts[, 1], ylim = c(0, 50), xlim = c(0, 3))
abline(h = 10, col = "red")
abline(v = 0.25, col = "red")
#dev.off()
keep <- rowSums(cpm(counts[]) > 0.25) >= 4
gExpr <- counts[keep, ]
dim(gExpr)
```

#### 7.	Create DGEList object including gene annotation
```{r}
y <- DGEList(gExpr)
y <- calcNormFactors(y)
anno <- read.csv("anno.csv")
row.names(anno) <- anno$ensembl
anno <- anno[match(rownames(y), rownames(anno)), ]
y$genes <- anno
```

#### 8.	Create diagnostic plot:
Plot multidimensional scaling plot to assess variables that should be added as covariates (Fig. 5B).
```{r}
#pdf(paste0("results/", experiment.title, "-2_mds.pdf"))
for (i in 1:length(colnames(meta))){
  mds(y, meta[ ,i], colnames(meta)[i])
}
plotMDS(y)
#dev.off()
```

## Fitting linear model
#### 9.	Design model
Define the linear model to be used in the differential expression analysis. Add the variable of interest as well as all variables (metadata columns) that resulted from the visual inspection of the MDS plot in Step 8. In this example “mod.gene” represents our variable of interest, while the MDS plot showed “line” and “donor” as additional covariates. 
```{r}
design <- model.matrix(~ 0 + mod.gene + line + donor, meta)
```
```{r}
#### Optional: Subsequently remove the column name from the created design matrix to ease its use.
colnames(design)
colnames(design) <- gsub("mod.gene", "", colnames(design))
colnames(design)
```

#### 10.	Voom transform 
The voom() function prepares the data for linear modeling by converting counts to logCPM and computing weights for heteroscedasticity adjustment. It also creates a diagnostic plot for the mean-variance trend (Fig. 5C). Visually inspect it for fit. Increase the cutoff for lowly expressed genes to improve fit if necessary.
```{r}
v <- voom(y, design, plot = TRUE, save.plot = TRUE)
```
```{r, fig.show="hide"}
## Save the plot (optional)
pdf(paste0("results/", experiment.title, "-3_voom.pdf"))
plot(v$voom.xy, type = "p", pch=20, cex=0.16, 
     main = "voom: Mean-variance trend", 
     xlab = "log2( count size + 0.5 )", 
     ylab = "Sqrt( standard deviation )")
lines(v$voom.line, col="red")
dev.off()
```

#### 11.	Fit model
```{r}
fit <- lmFit(v, design)
```

#### 12.	Define group comparisons (contrasts) 
Define your comparisons of interest. Each comparison is assigned a name and its function using elementary algebra. In addition to the standard comparisons of the individual and combined gene perturbations with their control, we also add equations for the additive and the synergistic model.
```{r}
cont.matrix <- makeContrasts(
  SNAP91a = sanp91 - ctrl,
  TSNARE1a = tsnare1 - ctrl,
  CLCN3a = clcn3 - ctrl,
  FURINi = furin - ctrl,
  Additive = sanp91 + tsnare1 + clcn3 + furin - 4*ctrl,
  Combinatorial = all - all.ctrl,
  Synergy = all - sanp91 - tsnare1 - clcn3 - furin - all.ctrl + 4*ctrl,
  levels = design)
```

#### 13.	Optional: Visualize the contrasts in a heatmap.
```{r}
cont.p <- t(cont.matrix)
h <- pheatmap(cont.p,
              display_numbers = T, number_format = "%.0f",
              breaks = seq(-3, 1, by = 0.5),
              color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(12), 
              cluster_cols = F, cluster_rows = F)
print(h)
```

#### 14.	Calculate coefficients and standard errors for each contrast.
```{r}
fit.cont <- contrasts.fit(fit, cont.matrix)
```

## Assessing differential expression
#### 15.	Empirical Bayes moderation.
```{r}
fit.cont <- eBayes(fit.cont)
plotSA(fit.cont, main = "Final model: Mean-variance trend", ylab = "Sqrt( standard deviation )")
summa.fit = decideTests(fit.cont, adjust.method = "fdr")
```

#### 16.	Save DEG result tables
Create list of all results. This will be used later to run analyses for all comparisons.
```{r}
res.list <- list()
for (i in 1:length(colnames(fit.cont$contrasts))){
  x <- topTable(fit.cont, coef = i, sort.by = "p", n = Inf, confint = T)
  res.list[[i]] <- x
  names(res.list)[i] <- colnames(fit.cont$contrasts)[i]
  write.csv(x, paste0("results/", experiment.title, "_DEGs_", 
                      colnames(fit.cont$contrasts)[i], ".csv"))
}
#### Optional: Create separate variable for each list entry.
list2env(res.list, globalenv())
```

#### 17.	Save DEG result plots
Create volcano and mean difference (MA) plots for each contrast, with all significant DE genes/ the top 10 genes highlighted, respectively (Fig. 5D).
```{r}
#pdf(paste0("results/", experiment.title, "-4_volcano-md-plots.pdf"))
par(mfrow = c(1, 2))
for (i in 1:length(colnames(fit.cont$contrasts))){
  plotMD(fit.cont, coef = i, status = summa.fit[, i], values = c(-1, 1))
  volcanoplot(fit.cont, coef = i, highlight = 10, 
              names = fit.cont$genes$Gene_name)
}
#dev.off()
par(mfrow = c(1, 1))
```

Optional: Plot expression for top 3 DE genes in each contrast in all samples.
Change the meta$mod.gene variable (red) to your variable of interest
"scale_x_dicrete" and "scale_fill_manual" (blue) are commented out and optional additions for aesthetics. *If used, their variable names have to be changed to reflect the data at hand.*
```{r optional DEG plots, fig.show="hide", results="hide"}
#pdf(paste0("results/", experiment.title, "-5_top3-expression-plots.pdf"))
par(mfrow = c(1, 1))
for (i in 1:length(colnames(fit.cont$contrasts))){
  x <- topTable(fit.cont, coef = i, sort.by = "p", n = Inf)
  cat("  \n\n### Plotting",  colnames(fit.cont$contrasts)[i], "  \n\n")
  for (j in 1:3){
    deg <- as.character(x[j,"ensembl"])
    p <- qplot(meta$mod.gene, v$E[deg, ], 
               geom = "boxplot", fill = meta$mod.gene, 
               ylab = "Normalized expression", xlab = "group", 
               main = paste0(j, ". DEG: ", as.character(x[j, "Gene_name"]))) +
      geom_jitter() +
      #scale_x_discrete(limits = c("ctrl", "sanp91", "tsnare1", "clcn3", "furin", "all.ctrl", "all")) +
      #scale_fill_manual(values = (c("orchid4", "grey", "steelblue", "grey", "firebrick", "blue", "darkblue"))) +
      rotate_x_text(angle = 45) +
      theme_bw(14)+
      theme(legend.position = "none", 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
            print(p)
  }
}
dev.off()
```

## Determining power to detect synergistic effects
#### 18.	Calculate power. Power calculations should be performed on a similar dataset in advance.
Calculate mean standard error for all measured comparisons.
```{r}
SE <- sqrt(fit.cont$s2.post) * fit.cont$stdev.unscaled
#### Calculate power. Choose the SE matrix column names that represent the additive and the combinatorial perturbation to calculate the median standard error. Then use them to run the power.compare.logFC() function, which creates a power plot (Fig. 6A).
colnames(SE)
sig1 <- median(SE[,"Additive"])
sig2 <- median(SE[,"Combinatorial"])
g <- power.compare.logFC(sig1, sig2, N = 4, N_other = c(4, 6, 8, 10, 14), 
                         alpha = 0.05, n_tests = 20000)
#pdf(paste0("results/", experiment.title, "-6_synergy-power.pdf"))
print(g)
#dev.off()
```

## Determining the extent of synergy
#### 19.	Determine extent of synergy.
Calculate synergy coefficient and percentage of synergistic DE genes (FDR < 10%). Plot histogram of all synergistic P-values to visualize the distribution (Fig. 6B).
```{r}
synergy.pvalues <- res.list$Synergy$P.Value
pi1 <- 1 - qvalue(synergy.pvalues)$pi0
print(pi1)
#pdf(paste0("results/", experiment.title, "-7_synergy-coefficient.pdf"))
plot.new()
text(0.4,0.75,labels=paste0("\n",round(pi1 * 100, 2), 
                            "% non-null \np-values and \n", 
                            round(sum(res.list$Synergy$adj.P.Val < 0.1) *
                                    100/length(res.list$Synergy$ensembl), 2), 
                            " % of genes with \nsynergy FDR < 0.1"))
hist(synergy.pvalues)
#dev.off()
```

## Identification and categorization of synergistic genes
#### 20.	Define synergistic effect categories.
Determine an expression cutoff range. Here +/- the mean standard error for all empirically measured comparisons is used.
```{r}
meanSE = mean(SE[,c(1,2,3,4,6)])
```

Create a table combining log2 fold change columns of additive, combinatorial and synergy contrasts. Make certain the vectors (marked red) refer to the column indices of ensembl ID, gene name, logFC and adjusted P-value in the DEG result tables.
```{r}
colnames(res.list$Additive)
log2FC.matrix <- Reduce(function(x,y) merge(x,y,by=c("ensembl", "Gene_name"), all = TRUE),
                        list(res.list$Additive[,c(1,2,4,9)],
                             res.list$Combinatorial[,c(1,2,4,9)],
                             res.list$Synergy[,c(1,2,4,9)]))
colnames(log2FC.matrix)
colnames(log2FC.matrix) <- c("Ensembl", "Gene_name", 
                             "Additive.logFC", "Additive.FDR", 
                             "Combinatorial.logFC", "Combinatorial.FDR", 
                             "Synergistic.logFC", "Synergistic.FDR")
rownames(log2FC.matrix) <- log2FC.matrix$Ensembl
```

Add a column assigning synergistic categories to each gene using the categorize.synergy() function, which takes the previously created matrix and the expression cutoff into account.
```{r}
log2FC.matrix <- categorize.synergy(log2FC.matrix, meanSE)
write.csv(log2FC.matrix, paste0("results/", experiment.title, "_LogFC-FDR-synergy_matrix.csv"))
genes.per.category <- count(log2FC.matrix, "magnitude.syn")
print(genes.per.category)
```

#### 21.	Visualize the categories in a pie chart.
Create a table containing sums of synergistic category genes.
```{r}
genes.per.category$category <- factor(genes.per.category$magnitude.syn,
                                      levels=c("same", "less.up", "less.down", "more.up", "more.down"))
genes.per.category$percent <- paste0(round(genes.per.category$freq*100/sum(genes.per.category$freq), 0), " %")
write.csv(genes.per.category, paste0("results/", experiment.title, "_gene-count_synergy-categories.csv"))
```

Plot a pie chart (Fig. 6C).
```{r}
zissou <- wes_palette("Zissou1", 6, type = "continuous")
q <- ggplot(genes.per.category, aes(x = "", y = freq, fill = category)) +
  geom_col() +
  coord_polar("y", start=0) + 
  scale_fill_manual(values=zissou) +
  theme_void()
#pdf(paste0("results/", experiment.title, "-8_synergy-categories_pie-chart.pdf"))
print(q)
#dev.off()
```

#### 22.	Visualize the categories in heatmaps
Create heatmaps of the log2FC in the additive and the combinatorial comparisons for each synergy category (Fig. 6D).
```{r}
#pdf(paste0("results/", experiment.title, "-9_heatmaps_log2FC-Add-vs-Combi.pdf"))
for (i in 1:length(levels(log2FC.matrix$magnitude.syn))){
  breaks <- c(seq(-6, -0.3,by=0.1),seq(0.3, 6,by=0.1))
  breaks <- append(breaks, -9,0)
  breaks <- append(breaks, 9)
  tmp <- log2FC.matrix[log2FC.matrix$magnitude.syn == 
                         levels(log2FC.matrix$magnitude.syn)[i],
                       c("Additive.logFC","Combinatorial.logFC")]
  h <- pheatmap(tmp,
                kmeans_k = 30, 
                cellwidth = 70, cellheight = 5,
                border_color = NA, 
                breaks=breaks,
                cluster_cols = F,
                show_rownames = F,
                color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(117),
                main = paste0("logFC expected vs. measured:\n",
                              levels(log2FC.matrix$magnitude.syn)[i]))
  print(h)
}
#dev.off()
```

## Enrichment analysis

### All comparisons: GSEA

#### 23.	Create a list containing the gene set groups of interest in the required format using ids2indices(), geneIds() and getGmt(). 
In this example these are manually curated gene set groups, saved in the “genesets” folder: disorder.gmt, behavior.gmt, connectivity.gmt, head.gmt, neural.gmt, postsynapse.gmt, presynapse.gmt and synapse.gmt.
```{r, warning=FALSE, message=FALSE}
gs.list <- list(
  "disorder" = ids2indices(geneIds(getGmt("genesets/disorder.gmt")), id = v$genes$Gene_name),
  "behavior" = ids2indices(geneIds(getGmt("genesets/behavior.gmt")), id=v$genes$Gene_name),
  "connectivity" = ids2indices(geneIds(getGmt("genesets/connectivity.gmt")), id=v$genes$Gene_name),
  "head" = ids2indices(geneIds(getGmt("genesets/head.gmt")), id=v$genes$Gene_name),
  "neural" = ids2indices(geneIds(getGmt("genesets/neural.gmt")), id=v$genes$Gene_name),
  "postsynapse" = ids2indices(geneIds(getGmt("genesets/postsynapse.gmt")), id=v$genes$Gene_name),
  "presynapse" = ids2indices(geneIds(getGmt("genesets/presynapse.gmt")), id=v$genes$Gene_name),
  "synapse" = ids2indices(geneIds(getGmt("genesets/synapse.gmt")), id=v$genes$Gene_name))
```

#### 24.	Optional: Create a custom color palette.
Create a custom color palette with at least as many colors as gene set groups to be tested.
```{r}
catcols=c("behavior" = "#8c3800",     #brown
          "disorder" = "#3C4347", 		#darkgrey 
          "connectivity" = "#e0a81c", #mustard
          "head" = "#702658", 			  #purple
          "neural" = "#004878", 		  #blue
          "postsynapse" = "#486030", 	#darkgreen
          "presynapse" = "#a8c018",		#lightgreen
          "synapse" = "#5c9340") 		  #green
# additional colors:
# "#6CA7AD", "#B72415", "#A88C05", "#E06C03", "#653C82", 
# "#0287AA", "#AD5A60", "#9B0420")
show_col(catcols)
```

#### 25.	Run gene set enrichment for all comparisons using camera.
Loop through all contrasts in the cont.matrix object. Perform enrichment and visualize as scatter (Fig. 7A) and bar plots (Fig. 7B) using the custom cameraplusplots() function.
```{r, warning=FALSE, message=FALSE, results="hide"}
camera.res.list <- list()
for (j in 1:length(colnames(cont.matrix))){
  print(paste0("Contrast: ", colnames(cont.matrix)[j]))
  #pdf(paste0("results/", experiment.title, "-10_GSEA-", colnames(cont.matrix)[j], "-plots.pdf"))
  camera.res <- cameraplusplots(contrast = cont.matrix[ ,j], 
                                genesetlist = gs.list, vobject = v, design = design, 
                                catcolors = catcols, title = paste0(colnames(cont.matrix)[j]))
  #dev.off()
  camera.res.list[[j]] <- camera.res
  names(camera.res.list)[j] <- colnames(cont.matrix)[j]
  write.csv(data.frame(camera.res), paste0("results/", experiment.title,
                                           "_GSEA-", colnames(cont.matrix)[j], ".csv"))
}
```

Plot legend (Fig. 7A).
```{r}
#pdf(paste0("results/", experiment.title, "-10_GSEA-plot-legend.pdf"))
plot(1,type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
legend("center", names(catcols), cex = 1.2, fill = catcols)
#dev.off()
```


### Specific gene subsets: ORA

#### 26. Create synergistic gene subsets to be analyzed.
Adjust FDR cutoff depending on the subtlety of the synergistic effect. In this example a cutoff of synergistic FDR < 1% was chosen.
```{r}
log2FC.matrix.sub <- subset(log2FC.matrix, Synergistic.FDR < 0.01)
```

Further possibly interesting subsets (commented out):
```{r}
  ##### Synergy genes with combinatorial FDR < 5%
  #log2FC.matrix.sub <- subset(log2FC.matrix, Combinatorial.FDR < 0.05)
  ##### Genes with synergistic Fold Change > 2 or < 0.5
  #log2FC.matrix.sub <- subset(log2FC.matrix, Synergistic.logFC > 1.5 | Synergistic.logFC < -1.5)
```

#### 27. Stratify chosen subset by synergy category using the custom stratify.by.syn.cat() function.
```{r}
syn.cat.list <- stratify.by.syn.cat(log2FC.matrix.sub)
```

#### 28.	Define reference genes as all genes analyzed. 
Since lowly expressed genes were filtered out at the beginning of the protocol, this can be interpreted as all expressed genes.
```{r}
allgenes <- as.character(y$genes$Gene_name)
```

#### 29.	Create a list of paths for the gene sets of interest. 
WebgestaltR, which is used for over-representation analysis, requires gene sets to be provided in a different format than camera.
```{r}
gs.list <- list(
  "disorder" = "genesets/disorder.gmt",
  "behavior" = "genesets/behavior.gmt",
  "connectivity" = "genesets/connectivity.gmt",
  "head" = "genesets/head.gmt",
  "neural" = "genesets/neural.gmt",
  "postsynapse" = "genesets/postsynapse.gmt",
  "presynapse" = "genesets/presynapse.gmt",
  "synapse" = "genesets/synapse.gmt")
```

#### 30.	Run over-representation analysis for “more up” and “more down” synergy categories using the WebGestaltR() function.
Loop through the “more.up” and “more.down” vectors in the previously created list object syn.cat.list.
Create a data frame for all results.
Visualize using the custom oraplot() function (Fig. 8A, B).
```{r, warning=FALSE, message=FALSE, results="hide"}
for (i in 3:4){
  tryCatch({
    ora <- data.frame(matrix(ncol = 11, nrow = 0))
    ora.list <- list()
    goi <- syn.cat.list[[i]]
    for (j in 1:length(gs.list)){
      tryCatch({
        ora.s <- WebGestaltR(enrichMethod = "ORA", 
                             organism = "hsapiens",
                             interestGene = goi,
                             interestGeneType = "genesymbol",
                             referenceGene = allgenes, 
                             referenceGeneType = "genesymbol",
                             enrichDatabase = "others",
                             enrichDatabaseFile = file.path(gs.list[j]), 
                             enrichDatabaseType = "genesymbol",
                             sigMethod = "top", topThr = 50, minNum = 3,
                             isOutput = F)
        ora.list[[j]] <- ora.s
        names(ora.list)[j] <- names(gs.list)[j]
        ora.list[[j]]$category <- names(ora.list[j])
        colnames(ora) <- names(ora.list[[1]])
        ora <- rbind.data.frame(ora, ora.list[[j]])
      }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    write.csv(ora, paste0("results/", experiment.title, "_ORA-", 
                          levels(log2FC.matrix$magnitude.syn)[i], ".csv"))
    
    g <- oraplot(ora, catcols, 
                 paste0(levels(log2FC.matrix$magnitude.syn)[i]))
    #pdf(paste0("results/", experiment.title, "-11_ORA-", levels(log2FC.matrix$magnitude.syn)[i], "-plots.pdf"))
    print(g)
    #dev.off()
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
```
