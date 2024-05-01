#Calculation
#source("/Users/dokada/Dropbox/analysis/2022.5/opp_entropy0416.R") #Run2024.4.27 in R 4.3.2
out_path <- "/Users/dokada/Desktop/work/opp_entropy0416/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Prepare data
library(TabulaMurisSenisData)
library(WebGestaltR)
library(SingleCellExperiment)
library(effsize)
set.seed(1000)
dataset <- c("FACS", "Droplet") #FACS, droplet
an_type <- c("Entropy", "Num_of_Expressed_Genes")
dbs <- c("geneontology_Biological_Process_noRedundant") #for GSEA
#dat <- NULL
#dat_colnames <- NULL
for(ds in 1:2){

    #set parameters
    rm(list=setdiff(ls(), c("out_path", "dataset", "an_type", "ds", "dbs")))
    gc()
    tmp_ds <- dataset[ds]

    #load data
    if(tmp_ds == "FACS"){
        raw_data <- TabulaMurisSenisFACS(tissues = "All", processedCounts = TRUE)
    }else{
        raw_data <- TabulaMurisSenisDroplet(tissues = "All", processedCounts = TRUE)
    }
    count_mat <- as.matrix(logcounts(raw_data$All))
    sample_annot <- colData(raw_data$All)
    gene_annot <- rowData(raw_data$All)

    #Association analysis
    id_all <- as.character(sample_annot$mouse.id)
    if(tmp_ds == "FACS"){
        age_all <-  sapply(id_all, function(chr){as.numeric(strsplit(chr,"_")[[1]][1])})
    }else{
        age_all <- sapply(id_all, function(chr){as.numeric(strsplit(chr,"-")[[1]][1])})
    }
    sex_all <- sample_annot$sex
    tissue_all <- sample_annot$tissue
    celltype_all <- sample_annot$cell_ontology_class
    stcia_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all, ".", id_all,".", age_all)
    stc_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all)
    stca_all <- paste0(sex_all, ".", tissue_all, ".", celltype_all, ".", age_all)

    for(an in 1:2){

        tmp_an <- an_type[an]

        #Calculate entropy
        div <- rep(NA, ncol(count_mat))
        if(tmp_an == "Entropy"){
            for(i in 1:ncol(count_mat)){
                vec <- count_mat[,i]
                p <- vec/sum(vec)
                logp <- log2(p)
                logp[logp==-Inf] <- 0
                div[i] <- - sum(p * logp)
                cat(i, "\n")
            }
        }else{
            for(i in 1:ncol(count_mat)){
                vec <- count_mat[,i]
                div[i] <- log2(sum(vec > 0))
                cat(i, "\n")
            }
        }

        #Asociation analysis to age
        unq_stc <- sort(unique(stc_all))
        pvals_ent <- estimates_ent <- rep(NA, length(unq_stc))
        names(pvals_ent) <- names(estimates_ent) <- unq_stc
        passed_idx <- NULL
        for(i in 1:length(unq_stc)){
            young_idx <- which((stc_all == unq_stc[i]) & (age_all == 3))
            old_idx <- which((stc_all == unq_stc[i]) & (age_all >= 18))
            if(length(young_idx) < 25 | length(old_idx) < 25){
                pvals_ent[i] <- NA
                estimates_ent[i] <- NA
                next
            }else{
                passed_idx <- c(passed_idx, i)
            }
            young_div <- div[young_idx]
            old_div <- div[old_idx]
            test_res_ent <- t.test(old_div, young_div)
            pvals_ent[i] <- test_res_ent$p.value
            cohend <- cohen.d(old_div, young_div) #ssv - nsv
            estimates_ent[i] <- cohend$estimate
            #cat(i, "\n")
        }
        unq_stc2 <- unq_stc[passed_idx]
        pvals_ent2 <- pvals_ent[passed_idx]
        estimates_ent2 <- estimates_ent[passed_idx]
        adj_pvals_ent2 <- p.adjust(pvals_ent2, method="bonferroni")

        #Volcano plot (Entropy)
        threshold <- 0.05 / length(pvals_ent2)
        png(paste0(out_path, tmp_ds, tmp_an, ".volc.png"), width=960, height=960)
        par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
        xlab = "Cohens's D (Old - Young)"
        ylab = "-log10(Adjusted p-values)"
        main <- paste0("", tmp_ds, ".", tmp_an, "" )
        cols <- rep("black", length(estimates_ent2))
        cols[pvals_ent2 < threshold & estimates_ent2 > 0] <- "red"
        cols[pvals_ent2 < threshold & estimates_ent2 < 0] <- "blue"
        plot(estimates_ent2,  -log10(pvals_ent2), col=cols, xlab="", ylab="", main="", cex=3, pch = 19, cex.axis=4)
        abline(h = -log10(threshold), col = "red", lty = 2, lwd=4)
        mtext(xlab, side=1, line=6, cex=4)
        mtext(ylab, side=2, line=6, cex=4)
        mtext(main, side=3, line=3, cex=4, adj=0)
        dev.off()

        #wariai
        tab_col <- table(cols)
        tab_col <-  tab_col/sum(tab_col) * 100
        write.csv(tab_col, paste0(out_path, tmp_ds, tmp_an, ".tab_col.csv"), row.names=TRUE)


        #data output
        dat <- data.frame(unq_stc2, pvals_ent2, estimates_ent2)
        write.csv(dat, paste0(out_path, tmp_ds, ".", tmp_an, "div_association.csv"), row.names=FALSE)
        #dat <- cbind(dat, unq_stc2, pvals_ent2, estimates_ent2, adj_pvals_ent2)
        #dat_colnames <- c(dat_colnames, paste0(tmp_ds, ".", tmp_an, c("Cell_subset", "P-value", "Cohen's D (Old - Young)", "Adjusted P-value")))
    }
}
#colnames(dat) <- dat_colnames
