#Calculation
#source("/Users/dokada/Dropbox/analysis/2022.5/opp_facsdrop0416.R") #Run2024.4.30 in R 4.3.2
out_path <- "/Users/dokada/Desktop/work/opp_facsdrop0416/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}


#Prepare data
set.seed(1000)
library(TabulaMurisSenisData)
library(WebGestaltR)
library(SingleCellExperiment)
library(vegan)
dbs <- c("geneontology_Biological_Process_noRedundant") #for GSEA
dataset <- c("FACS", "Droplet") #FACS, droplet
enrichResult_list <- list()
for(ds in 1:2){

    #set parameters
    rm(list=setdiff(ls(), c("out_path", "dataset", "an_type", "ds", "dbs", "enrichResult_list")))
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

    #DEG analysis based on t test
    unq_stc <- sort(unique(stc_all))
    used_genes <- rownames(count_mat)
    n_genes <- length(used_genes)
    scores <-  scores_p <- scores_y <- scores_o <- matrix(NA, nrow = length(unq_stc), ncol = n_genes)
    rownames(scores) <- rownames(scores_p) <- rownames(scores_y) <- rownames(scores_o) <- unq_stc
    colnames(scores) <- colnames(scores_p) <- colnames(scores_y) <- colnames(scores_o) <- used_genes
    passed_idx <- NULL
    for(i in 1:length(unq_stc)){
        young_idx <- which((stc_all == unq_stc[i]) & (age_all == 3))
        old_idx <- which((stc_all == unq_stc[i]) & (age_all >= 18))
        count_mat_young <- count_mat[used_genes, young_idx]
        count_mat_old <- count_mat[used_genes,old_idx]
        if(length(young_idx) < 25 | length(old_idx) < 25){
            scores[i,] <- NA
            scores_p[i,] <- NA
            scores_y[i,] <- NA
            scores_o[i,] <- NA
            next
        }else{
            passed_idx <- c(passed_idx, i)
        }

        #Calculation (DEG)
        for(j in 1:n_genes){
            vec_young <- count_mat_young[j,]
            vec_old <- count_mat_old[j,]

            #DEG
            scores[i, j] <- mean(vec_old) - mean(vec_young)
            if(sd(vec_young) == 0 | sd(vec_old) == 0){
                scores_p[i, j] <- NA
                scores_y[i, j] <- NA
                scores_o[i, j] <- NA
                next
            }else{
                res <- t.test(vec_young, vec_old)
                scores_p[i, j] <- res$p.value
                scores_y[i, j] <- mean(vec_young)
                scores_o[i, j] <- mean(vec_old)
            }
        }

        cat(i, "\n")
    }

    #Passed subset
    scores_pass <- scores[passed_idx,]
    scores_p_pass <- scores_p[passed_idx,]
    scores_y_pass <- scores_y[passed_idx,]
    scores_o_pass <- scores_o[passed_idx,]
    write.csv(scores_pass, paste0(out_path, tmp_ds, ".scores_pass.csv"))
    print(sum(is.na(scores_pass))) #No NA


    #score_sig行列の計算
    n_genes <- ncol(scores_pass)
    n_stcpass <- nrow(scores_pass)
    scores_sig <- matrix(NA, nrow=n_stcpass, ncol=n_genes)
    for(i in 1:n_stcpass){
        pval <- scores_p_pass[i,]
        signs <- sign(scores_pass[i,])
        adjp <- pval * nrow(scores_pass) * ncol(scores_pass)
        #adjp <- pval *  ncol(scores_pass)
        adjp[is.nan(adjp)] <- 1
        tmp_sig <- signs * ifelse(adjp < 0.05, 1, 0)
        scores_sig[i,] <-  tmp_sig
    }
    colnames(scores_sig) <- colnames(scores)
    rownames(scores_sig) <- rownames(scores_pass)

    #遺伝子ごとに集計
    sig_tab <- matrix(NA, nrow=n_genes, ncol=3)
    rownames(sig_tab) <- colnames(scores_pass)
    colnames(sig_tab) <- c("UP", "DOWN", "Not_change")
    for(i in 1:n_genes){
        tmp_sig <- scores_sig[,i]
        sig_tab[i,"UP"] <- length(which(tmp_sig == 1))
        sig_tab[i,"DOWN"] <- length(which(tmp_sig == -1))
        sig_tab[i,"Not_change"] <- length(which(tmp_sig == 0))
    }

    #histogram
    opp_score <- apply(sig_tab, 1, min)/nrow(scores_sig)
    png(paste0(out_path, tmp_ds, ".opphist.png"), width=960, height=960)
    par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
    xlab = "Opposite Effect Score"
    ylab = "Frequency"
    main = paste0(tmp_ds)
    hist(opp_score, cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", main="")
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4)
    dev.off()
    #GSEA for OPP
    gene_set <- names(opp_score[opp_score>0.05])
    bg_geneset <- names(opp_score)
    dbs <- c("geneontology_Biological_Process_noRedundant")
    enrichResult <- WebGestaltR(enrichMethod="ORA", organism="mmusculus", enrichDatabase=dbs, interestGene=gene_set, interestGeneType="genesymbol", 
    referenceGene=bg_geneset, referenceGeneType="genesymbol", outputDirectory=out_path, fdrThr=0.05, projectName=paste0("ORA_result_", tmp_ds))
    enrichResult_list[[ds]] <- enrichResult
    write.csv(enrichResult, paste0(out_path, tmp_ds, ".opp_ora.csv"), row.names=FALSE)

    #Opp common barplot
    n_subset <- nrow(scores_pass)
    n_top_opp <- 10
    opp_nums <- sig_tab[names(sort(opp_score[order(opp_score, decreasing=T)[1:n_top_opp]])),c("UP", "DOWN"), drop=F]
    write.csv( opp_nums, paste0(out_path, tmp_ds, ". opp_nums", ".csv"))
    png(paste0(out_path, tmp_ds, ".all.barplot.png"), width=960, height=960)
    par(mar=c(5,20,5,5)) #bot, left, top, roght
    par(oma=c(5,2,2,5))
    par(mgp=c(6,0,0)) #軸からラベル, メモリ、軸線
    main = paste0(tmp_ds)
    barplot(t(opp_nums), beside = TRUE, las = 1, horiz=T, col = c("red", "blue"), names.arg = rownames(opp_nums), cex.axis=3, cex.lab=3, cex.names=4)
    mtext(main, side=3, line=3, cex=4)
    mtext("Num of subsets", side=1, line=6, cex=4)
    dev.off()

    #plot of UP and Down
    set.seed(1000)
    n_genes <- nrow(sig_tab)
    cols <- rep("black", n_genes)
    cols[sig_tab[,"UP"] > 0 & sig_tab[,"DOWN"] == 0] <- "blue"
    cols[sig_tab[,"DOWN"] > 0 & sig_tab[,"UP"] == 0] <- "red"
    cols[sig_tab[,"UP"] > 0 & sig_tab[,"DOWN"] > 0] <- "purple"
    png(paste0(out_path, tmp_ds, ".up_down.png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    xlab = "Num of subsets (Up)"
    ylab = "Num of subsets (Down)"
    main <- paste0(tmp_ds)
    x_up <- sig_tab[,"UP"]
    y_down <- sig_tab[,"DOWN"]
    plot(jitter(x_up), jitter(y_down), col=cols, xlab="", ylab="", main="", cex=3, pch = 19, cex.axis=4)
    mtext(xlab, side=1, line=6, cex=4)
    mtext(ylab, side=2, line=6, cex=4)
    mtext(main, side=3, line=3, cex=4, adj=0)
    dev.off()


    #GenAgeのデータベースを使って、老化関連遺伝子の中での内訳を集計
    gene_age <- read.csv("/Users/dokada/Desktop/work/detarmin_data/genage_models_export.csv", header=TRUE)
    tab_geneage <- table(gene_age$Gene.Symbol)
    #Correct "GMFB"  "NUDT1"  "SOD3"  "G6PD"
    names(tab_geneage)[names(tab_geneage) == "SOD3"] <- "Sod3"
    names(tab_geneage)[names(tab_geneage) == "NUDT1"] <- "Nudt1"
    names(tab_geneage)[names(tab_geneage) == "GMFB"] <- "Gmfb"
    names(tab_geneage)[names(tab_geneage) == "G6PD"] <- "H6pd" #refer to https://www.ncbi.nlm.nih.gov/gene/100198
    used_symbols <- intersect(names(tab_geneage), rownames(sig_tab))
    #if(length(used_symbols) != length(names(tab_geneage))) stop("Error")
    pi_gene <- c("All", "GenAge")
    for(i in 1:2){
        if(i == 1){
            tmp_symbols <- rownames(sig_tab)
        }else{
            tmp_symbols <- used_symbols
        }
        up_common <- sum(sig_tab[tmp_symbols,"UP"] > 0 & sig_tab[tmp_symbols,"DOWN"] == 0)
        down_common <- sum(sig_tab[tmp_symbols,"DOWN"] > 0 & sig_tab[tmp_symbols,"UP"] == 0)
        opps <- sum(sig_tab[tmp_symbols,"UP"] > 0 & sig_tab[tmp_symbols,"DOWN"] > 0)
        stay_common <- sum(sig_tab[tmp_symbols,"UP"] == 0 & sig_tab[tmp_symbols,"DOWN"] == 0)
        png(paste0(out_path, tmp_ds, pi_gene[i], ".pie.png"), width=960, height=960)
        par(mar = c(4, 18, 4, 18)) ##bottom, left, top, right
        x <- c(up_common, down_common, opps, stay_common)
        wariai <- paste0(" (", round(x/sum(x)*100,1),"%)")
        lbls <- c("UP_only", "DOWN_only", "Opposite", "Not_change")
        main <- paste0(tmp_ds)
        pie(x, labels = paste0(lbls, wariai), col = c("red", "blue", "purple", "black"), cex=2.5, main = main, cex.main=4)
        dev.off()
        names(x) <- lbls
        write.csv(x, paste0(out_path, tmp_ds, pi_gene[i], ".pie.csv"))
    }

    #Opp common example plot
    opp_score_genage <-  opp_score[used_symbols]
    n_top_opp <- 10
    opp_nums <- sig_tab[names(sort(opp_score_genage[order(opp_score_genage, decreasing=T)[1:n_top_opp]])),c("UP", "DOWN"), drop=F]
    write.csv(opp_nums, paste0(out_path, tmp_ds, ". opp_nums", ".csv"))
    png(paste0(out_path, tmp_ds, ".genage.barplot.png"), width=960, height=960)
    par(mar=c(5,10,5,5)) #bot, left, top, roght
    par(oma=c(5,2,2,5))
    par(mgp=c(6,0,0)) #軸からラベル, メモリ、軸線
    main = paste0(tmp_ds)
    barplot(t(opp_nums), beside = TRUE, las = 1, horiz=T, col = c("red", "blue"), names.arg = rownames(opp_nums), cex.axis=3, cex.lab=3, cex.names=4)
    mtext(main, side=3, line=3, cex=4)
    mtext("Num of subsets", side=1, line=6, cex=4)
    dev.off()

    #Example plot
    out_dir <- paste0(out_path, tmp_ds, "_gene/")
    if(file.exists(out_dir)){
        unlink(out_dir, recursive=TRUE)
        dir.create(out_dir)
    }else{
        dir.create(out_dir)
    }
    genes <- names(opp_score_genage[opp_score_genage > 0])
    for(gene in genes){
        png(paste0(out_dir, gene, tmp_ds, ".opp_geneage_top.png"), width=960, height=960)
        par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
        xlab = "Young"
        ylab = "Old"
        main <- paste0(gene, " (", tmp_ds, ")" )
        cols <- rep("black", nrow(scores_pass))
        cols[scores_sig[, gene] == 1] <- "red"
        cols[scores_sig[, gene] == -1] <- "blue"
        plot(scores_y_pass[, gene], scores_o_pass[, gene], col=cols, xlab="", ylab="", main="", cex=3, pch = 19, cex.axis=4)
        abline(c(0, 1), col = "red", lwd=4)
        mtext(xlab, side=1, line=6, cex=4)
        mtext(ylab, side=2, line=6, cex=4)
        mtext(main, side=3, line=3, cex=4, adj=0)
        dev.off()
    }

    #output the sigtab
    opposite_effect_score <- opp_score
    sig_tab2 <- cbind(sig_tab, opposite_effect_score)
    write.csv(sig_tab2, paste0(out_path, tmp_ds, ".sig_tab.csv"))

}

#four gene is not uncluded in Droplet dataset
facs_sig_tab <- read.csv(paste0(out_path, "FACS.sig_tab.csv"), row.names=1)
drop_sig_tab <- read.csv(paste0(out_path, "Droplet.sig_tab.csv"), row.names=1)
gene_age <- read.csv("/Users/dokada/Desktop/work/detarmin_data/genage_models_export.csv", header=TRUE)
unq_geneage <- unique(gene_age$Gene.Symbol)
unq_geneage[unq_geneage == "SOD3"] <- "Sod3"
unq_geneage[unq_geneage == "NUDT1"] <- "Nudt1"
unq_geneage[unq_geneage == "GMFB"] <- "Gmfb"
unq_geneage[unq_geneage == "G6PD"] <- "H6pd" #refer to https://www.ncbi.nlm.nih.gov/gene/100198
facs_sig_tab_sub <- facs_sig_tab[unq_geneage,]
facs_score <- apply(facs_sig_tab_sub[,c("UP", "DOWN")], 1, min)
drop_sig_tab_sub <- drop_sig_tab[unq_geneage,]
drop_score <- apply(drop_sig_tab_sub[,c("UP", "DOWN")], 1, min)
opp_geneage <- unq_geneage[intersect(which(facs_score > 0) , which(drop_score > 0))]
write.csv(opp_geneage, paste0(out_path, "opp_geneage_and.csv"), row.names=FALSE)
opp_geneage <- unq_geneage[union(which(facs_score > 0) , which(drop_score > 0))]
write.csv(opp_geneage, paste0(out_path, "opp_geneage_or.csv"), row.names=FALSE)

#Enrichment analysis
en1 <- enrichResult_list[[1]]
en2 <- enrichResult_list[[2]]
com_geneset <- intersect(en1$geneSet, en2$geneSet)
write.csv(com_geneset, paste0(out_path, "common_geneset.csv"), row.names=FALSE)

##Visualization
for(ds in 1:length(enrichResult_list)){
    tmp <- enrichResult_list[[ds]]
    en <- tmp$enrichmentRatio
    names(en) <-tmp$description
    en <- sort(en)
    left <- ifelse(ds==1, 40, 60)
    png(paste0(out_path, dataset[ds], ".ora.png"), width=960, height=960)
    par(mar = c(9, left, 2, 2)) ##bottom, left, top, right
    xlab = "Num of genes"
    ylab = "Gene set"
    main = paste0(dataset[ds])
    barplot(en, beside = TRUE, las = 1, horiz=T, col = "blue", names.arg = names(en), cex.axis=3, cex.lab=3, cex.names=2)
    dev.off()
}