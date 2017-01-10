setwd("/Users/sfoltz/Desktop/lab/msi/")

#Libraries
library(ggplot2)
library(MASS)
library(glmnet)

#Load data matrix from folder
load("msi_data_matrix.Rdata")
load("univariate_analysis_dfs.Rdata")
load("glm_logistic_models.Rdata")

#Import data and get into proper data frame format
{
  #Read in relevant data -- variable name refers to specfic information derived from that data frame
  beta_values <- read.table("data/MMR.beta.values.panUnion", header=T, fill=T)
  pointMut_MSIsensor <- read.table("data/Sample_CT_MLH1meth_PointMut_MSI.txt", header=T)
  cadd_maf <- read.table("data/Somatic_pan8000_mmr_tier1.union.CADD.maf", header=T, sep="\t")
  clinical_MSI <- read.table("data/TCGA_sample_cancer_MMRstatus_SG_BRAF_panel_vital_MMRstatus_all.union.txt", header=T)

  #Fix beta_values lines with NA in them
  #Format starts as  gene transcript sample NA NA cantype [empty]
  #Format changes to gene transcript sample NA NA NA cantype
  beta_values[!complete.cases(beta_values),6:7] <- beta_values[!complete.cases(beta_values),5:6]

  #Totally inclusive list of all sample names from all four data files
  #This list consists of cases (>0 mutations) and controls (no mutations)
  all_sample_names <- data.frame(sample=c(substring(beta_values$sample,1,12), substring(pointMut_MSIsensor$sample,1,12), substring(cadd_maf$Tumor_Sample_Barcode,1,12), substring(clinical_MSI$Sample,1,12)), cancer_type=c(as.character(beta_values$cantype), as.character(pointMut_MSIsensor$cantype), as.character(cadd_maf$tumor_type), as.character(clinical_MSI$Cancer_Type)))
  sample_names <- unique(sort(all_sample_names$sample))
  for(name in sample_names){
    if(nrow(unique(all_sample_names[all_sample_names$sample==name,]))!=1){
      cts <- unique(as.character(all_sample_names[all_sample_names$sample==name,"cancer_type"]))
      ct <- cts[which(cts %in% c("COAD","READ"))]
      all_sample_names[all_sample_names$sample==name,"cancer_type"] <- ct
      rm(list=c("cts","ct"))
    }
  }
  rm(name)
  all_sample_names_cancer_types <- unique(all_sample_names)
  all_sample_names_cancer_types <- all_sample_names_cancer_types[order(all_sample_names_cancer_types$sample),]
  rm(all_sample_names)

  #Plot showing CADD score of different variant classes (missense, nonsense, nonstop, splice site mutations)
  if(FALSE){
    plot_df <- data.frame(cadd_maf$CADD_Score, cadd_maf$Variant_Classification)
    plot_df <- droplevels(plot_df[complete.cases(plot_df),])
    names(plot_df) <- c("CADD_Score","Variant_Classification")
    levels(plot_df$Variant_Classification) <- paste0(c("Missense\n(", "Nonsense\n(", "Nonstop\n(", "Splice Site\n("), summary(plot_df$Variant_Classification), c(")",")",")",")"))
    p <- ggplot(plot_df, aes(y=CADD_Score, x=Variant_Classification, fill=Variant_Classification)) + geom_violin(alpha=0.5) + geom_boxplot(width=0.05, fill="black", outlier.color=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x="Variant Class", y="CADD Score", title="CADD Score by Variant Class") + guides(fill=F)
    rm(plot_df)
    pdf("plots/cadd_score_by_variant_class.pdf",10,10)
    print(p)
    dev.off()
    rm(p)
  }

  #genes with methylation data
  methylation_genes <- sort(levels((beta_values$gene)))
  #genes with mutation data
  mutation_genes <- sort(levels(cadd_maf$Hugo_Symbol))

  #canonical cancer type (Y/N), point mutation rate, tcga_msi(Stable, Low, High, NA)
  canonical_cancer_type <- vector()
  point_mutation_rate <- vector()
  tcga_msi <- vector()
  msisensor <- vector()
  count <- 0
  for(name in sample_names){
    count <- count + 1
    canonical_cancer_type[count] <- all_sample_names_cancer_types[all_sample_names_cancer_types==name,"cancer_type"] %in% c("COAD","READ","STAD","UCEC","COADREAD")
    if(name %in% substr(pointMut_MSIsensor$sample,1,12)){
      point_mutation_rate[count] <- unique(pointMut_MSIsensor[substr(pointMut_MSIsensor$sample,1,12) == name,"Point_Mut"])
      msisensor[count] <- unique(pointMut_MSIsensor[substr(pointMut_MSIsensor$sample,1,12) == name,"MSI"])
    }
    else{
      point_mutation_rate[count] <- NA
      msisensor[count] <- NA
    }
    if(name %in% substr(clinical_MSI$Sample,1,12)){
      tcga_msi[count] <- as.character(clinical_MSI[clinical_MSI$Sample==name,"TCGA_panel"])
    }
    else{
      tcga_msi[count] <- NA
    }

  }
  tcga_msi[is.na(tcga_msi) | tcga_msi=="INDETERMINATE"] <- "Missing"
  tcga_msih <- c(tcga_msi==c("MSI-H"))

  #start to build data frame
  df <- data.frame(all_sample_names_cancer_types, canonical_cancer_type, point_mutation_rate, tcga_msi, tcga_msih, msisensor)
  rm(list=c("count","all_sample_names_cancer_types","canonical_cancer_type","point_mutation_rate","tcga_msi","tcga_msih", "msisensor"))
  first_beta_col <- ncol(df)+1

  #Create main data frame
  #Format: sample_name cancer_type point_mutation_rate tcga_msi beta_[geneName*] meth_[geneName*] cadd_[geneName^] mut_[geneName^] nonsense_[geneName^]
  #*Genes for which we have methylation data
  #^Genes for which we have mutation status
  #beta_* is the methylation beta value [0,1]
  #meth_* is a quantile category like low/medium/high or meth{1:5}
  #cadd_^ is the CADD score of the variant, 0 otherwise
  #mut_^ is a 0 or 1 for having any mutation there
  #nonsense_^ is a 0 or 1 for it being a nonsense (truncation) variant

  #beta_[geneName] values
  for(gene_name in methylation_genes){
    count <- 0
    beta_vector <- vector()
    vector_name <- paste0("beta_",gene_name)
    #get info for this gene
    this_gene_betas <- beta_values[beta_values$gene==gene_name,]
    this_sample <- substr(this_gene_betas$sample,1,12)
    for(name in sample_names){
      count <- count + 1
      #get info for this sample
      if(name %in% this_sample){
        beta_vector[count] <- mean(this_gene_betas[this_sample == name, "mean_beta"])
      }
      else{
        beta_vector[count] <- NA
      }
    }
    df <- data.frame(df, beta_vector)
    names(df)[ncol(df)] <- vector_name
    rm(list=c("count","beta_vector","gene_name","name","vector_name","this_gene_betas","this_sample"))
  }

  #Get 20% quantiles for each beta vector, store as meth_[geneName] in df
  if(TRUE){
    meth_df <- as.data.frame(apply(df[first_beta_col:ncol(df)], 2, function(x) cut(x, breaks=quantile(x, seq(0,1,0.2), na.rm=TRUE), include.lowest=TRUE,labels=paste0("meth",c(1:5)))))
    names(meth_df) <- paste0("meth_",methylation_genes)
    has_methylation_data <- apply(df[,first_beta_col:ncol(df)], 1, function(x) !all(is.na(x)))
    df <- data.frame(df, has_methylation_data, meth_df)
    rm(list=c("meth_df","has_methylation_data","first_beta_col"))
  }

  #cadd_[geneName], mut_[geneName], and nonsense_[geneName] for each sample
  cadd_df <- sample_names
  mut_df <- sample_names
  nonsense_df <- sample_names
  for(gene_name in mutation_genes){
    cadd_vector <- rep(0,nrow(df))
    mut_vector <- rep(0,nrow(df))
    nonsense_vector <- rep(0,nrow(df))
    this_gene_maf <- cadd_maf[cadd_maf$Hugo_Symbol==gene_name,c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","CADD_Score")]
    this_gene_samples <- substr(this_gene_maf$Tumor_Sample_Barcode,1,12)
    count <- 0
    if(nrow(this_gene_maf)>0){
      for(name in sample_names){
        count <- count + 1
        if(name %in% this_gene_samples){
          variant_class <- as.vector(this_gene_maf[this_gene_samples==name,"Variant_Classification"])
          cadd_scores <- this_gene_maf[this_gene_samples==name,"CADD_Score"]
          if(NaN %in% cadd_scores){
            variant_class <- variant_class[!is.na(cadd_scores)]
            cadd_scores <- cadd_scores[!is.na(cadd_scores)]
          }
          if(length(cadd_scores)>0){
            index <- which.max(cadd_scores)
            cadd_vector[count] <- cadd_scores[index]
            mut_vector[count] <- 1
            if(variant_class[index]=="Nonsense_Mutation"){
              nonsense_vector[count] <- 1
            }
            rm(index)
          }
          rm(list=c("variant_class","cadd_scores"))
        }
      }
    }
    rm(name)
    cadd_df <- data.frame(cadd_df, cadd_vector)
    names(cadd_df)[ncol(cadd_df)] <- paste0("cadd_",gene_name)
    mut_df <- data.frame(mut_df,mut_vector)
    names(mut_df)[ncol(mut_df)] <- paste0("mut_",gene_name)
    nonsense_df <- data.frame(nonsense_df,nonsense_vector)
    names(nonsense_df)[ncol(nonsense_df)] <- paste0("nonsense_",gene_name)
    rm(list=c("this_gene_maf","this_gene_samples","count","cadd_vector","mut_vector","nonsense_vector","gene_name"))
  }
  case <- apply(mut_df[,-1],1,function(x) max(x)==1)
  x <- data.frame(df, case, cadd_df[,-1], mut_df[,-1], nonsense_df[,-1])
  xcc <- x[complete.cases(x),]
  rownames(x) <- seq(length=nrow(x))
  tcga_sample_names <- sample_names
  rm(list=c("cadd_df","mut_df","nonsense_df","beta_values","cadd_maf","clinical_MSI","df","pointMut_MSIsensor","sample_names","case"))

  save(x, xcc, methylation_genes, mutation_genes, tcga_sample_names, file="msi_data_matrix.Rdata")
  rm(list=c("x", "methylation_genes", "mutation_genes", "tcga_sample_names", "xcc"))
}

#Dotplot of number of samples in each cancer type
{
  cts <- as.matrix(unique(sort(x$cancer_type)))
  num <- apply(cts,1,function(y) sum(x$cancer_type==y))
  plot_df <- data.frame(cts, num)
  p <- ggplot(plot_df, aes(x=num, y=reorder(cts,num))) + geom_point(size=3) + labs(x="Number of samples", y="Cancer type", title="Number of samples from each cancer type")
  pdf("plots/number_samples_by_cancer_type.pdf",10,10)
  print(p)
  dev.off()
  rm(list=c("cts","num","plot_df","p"))
}

#Dotplot with jitter of MSIsensor score by TCGA status
{
  num_msisensor_na <- sum(is.na(x$msisensor))
  plot_df <- x[!is.na(x$msisensor),c("tcga_msi","msisensor","cancer_type","canonical_cancer_type")]
  p <- ggplot(plot_df, aes(x=tcga_msi, y=msisensor, color=cancer_type, shape=canonical_cancer_type)) + geom_jitter() + labs(x="TCGA MSI status", y="MSIsensor score", title=paste0("MSIsensor score by TCGA MSI status\nNumber missing MSIsensor = ", num_msisensor_na), color="Cancer type", shape="COAD, READ,\nSTAD, or UCEC")
  pdf("plots/msisensor_by_tcga_status.pdf",10,10)
  print(p)
  dev.off()
  rm(list=c("num_msisensor_na","p","plot_df"))
}

#Dotplot with jitter of Point Mutation Rate score by TCGA status
{
  num_pmr_na <- sum(is.na(x$point_mutation_rate))
  plot_df <- x[!is.na(x$point_mutation_rate),c("tcga_msi","cancer_type","canonical_cancer_type","point_mutation_rate")]
  p <- ggplot(plot_df, aes(x=tcga_msi, y=point_mutation_rate, color=cancer_type, shape=canonical_cancer_type)) + geom_jitter() + labs(x="TCGA MSI status", y="Point mutation rate", title=paste0("Point Mutation Rate by TCGA MSI status\nNumber missing PMR = ", num_pmr_na), color="Cancer type", shape="COAD, READ,\nSTAD, or UCEC") + scale_y_log10()
  pdf("plots/point_mutation_rate_by_tcga_status.pdf",10,10)
  print(p)
  dev.off()
  rm(list=c("num_pmr_na","p","plot_df"))
}

#Univariate analysis of MSI-H vs. (MSI-L/MSS) status
{
  has_tcga_msi_status <- (xcc$tcga_msi %in% c("MSI-H","MSI-L","MSS"))
  xcc_tcga_msi <- droplevels(xcc[has_tcga_msi_status,])
  num_sim <- 1e4

  #betas
  log_reg_sim_beta <- vector()
  log_reg_p_beta <- vector()
  beta_genes <- paste0("beta_",methylation_genes)
  for(i in 1:length(beta_genes)){
    print(paste0("beta: ",i))
    simulated_z <- matrix(NA,num_sim,1)
    for(j in 1:num_sim){
      simulated_z[j] <- summary(glm(sample(xcc_tcga_msi$tcga_msih) ~ xcc_tcga_msi[,beta_genes[i]], family=binomial))$coefficients[2,3]
    }
    actual_z <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,beta_genes[i]], family=binomial))$coefficients[2,3]
    actual_p <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,beta_genes[i]], family=binomial))$coefficients[2,4]
    log_reg_sim_beta[i] <- 1 - mean(abs(actual_z) > abs(simulated_z))
    log_reg_p_beta[i] <- actual_p
  }

  #methylation quintiles
  log_reg_sim_meth <- vector()
  log_reg_p_meth <- vector()
  meth_genes <- paste0("meth_",methylation_genes)
  for(i in 1:length(meth_genes)){
    print(paste0("meth: ",i))
    simulated_z <- matrix(NA,num_sim,1)
    for(j in 1:num_sim){
      simulated_z[j] <- summary(glm(sample(xcc_tcga_msi$tcga_msih) ~ xcc_tcga_msi[,meth_genes[i]], family=binomial))$coefficients[5,3]
    }
    actual_z <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,meth_genes[i]], family=binomial))$coefficients[5,3]
    actual_p <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,meth_genes[i]], family=binomial))$coefficients[5,4]
    log_reg_sim_meth[i] <- 1 - mean(abs(actual_z) > abs(simulated_z))
    log_reg_p_meth[i] <- actual_p
  }

  #cadd scores
  log_reg_sim_cadd <- vector()
  log_reg_p_cadd <- vector()
  cadd_genes <- paste0("cadd_",mutation_genes)
  for(i in 1:length(cadd_genes)){
    print(paste0("cadd: ",i))
    simulated_z <- matrix(NA,num_sim,1)
    for(j in 1:num_sim){
      simulated_z[j] <- summary(glm(sample(xcc_tcga_msi$tcga_msih) ~ xcc_tcga_msi[,cadd_genes[i]], family=binomial))$coefficients[2,3]
    }
    actual_z <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,cadd_genes[i]], family=binomial))$coefficients[2,3]
    actual_p <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,cadd_genes[i]], family=binomial))$coefficients[2,4]
    log_reg_sim_cadd[i] <- 1 - mean(abs(actual_z) > abs(simulated_z))
    log_reg_p_cadd[i] <- actual_p
  }

  #mut events
  log_reg_sim_mut <- vector()
  log_reg_p_mut <- vector()
  mut_genes <- paste0("mut_",mutation_genes)
  for(i in 1:length(mut_genes)){
    print(paste0("mut: ",i))
    simulated_z <- matrix(NA,num_sim,1)
    for(j in 1:num_sim){
      simulated_z[j] <- summary(glm(sample(xcc_tcga_msi$tcga_msih) ~ xcc_tcga_msi[,mut_genes[i]], family=binomial))$coefficients[2,3]
    }
    actual_z <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,mut_genes[i]], family=binomial))$coefficients[2,3]
    actual_p <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,mut_genes[i]], family=binomial))$coefficients[2,4]
    log_reg_sim_mut[i] <- 1 - mean(abs(actual_z) > abs(simulated_z))
    log_reg_p_mut[i] <- actual_p
  }

  #nonsense events
  log_reg_sim_nonsense <- vector()
  log_reg_p_nonsense <- vector()
  nonsense_genes <- paste0("nonsense_",mutation_genes)
  for(i in 1:length(nonsense_genes)){
    print(paste0("nonsense: ",i))
    if(max(xcc_tcga_msi[,nonsense_genes[i]])==0){
      log_reg_sim_nonsense[i] <- NA
      log_reg_p_nonsense[i] <- NA
      next
    }
    else{
      simulated_z <- matrix(NA,num_sim,1)
      for(j in 1:num_sim){
        simulated_z[j] <- summary(glm(sample(xcc_tcga_msi$tcga_msih) ~ xcc_tcga_msi[,nonsense_genes[i]], family=binomial))$coefficients[2,3]
      }
      actual_z <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,nonsense_genes[i]], family=binomial))$coefficients[2,3]
      actual_p <- summary(glm(xcc_tcga_msi$tcga_msih ~ xcc_tcga_msi[,nonsense_genes[i]], family=binomial))$coefficients[2,4]
      log_reg_sim_nonsense[i] <- 1 - mean(abs(actual_z) > abs(simulated_z))
      log_reg_p_nonsense[i] <- actual_p
    }
  }
  univariate_analysis_methylation_df <- data.frame(methylation_genes, log_reg_sim_beta, log_reg_sim_meth, log_reg_p_beta, log_reg_p_meth)
  univariate_analysis_mutation_df <- data.frame(mutation_genes, log_reg_sim_cadd, log_reg_sim_mut, log_reg_sim_nonsense, log_reg_p_cadd, log_reg_p_mut, log_reg_p_nonsense)
  save(univariate_analysis_methylation_df, univariate_analysis_mutation_df, file="univariate_analysis_dfs.Rdata")
  rm(list=c("simulated_z","actual_p","actual_z","i","j","log_reg_sim_beta","log_reg_p_beta","log_reg_sim_meth","log_reg_p_meth","log_reg_sim_cadd","log_reg_p_cadd","log_reg_sim_mut","log_reg_p_mut","log_reg_sim_nonsense","log_reg_p_nonsense","beta_genes","meth_genes","cadd_genes","mut_genes","nonsense_genes","num_sim","has_tcga_msi_status","xcc_tcga_msi"))
}

#Build models with methylation Betas and mutation CADDs
{
  has_tcga_msi_data <- (xcc$tcga_msi %in% c("MSI-H","MSI-L","MSS"))
  xcc_tcga_msi <- subset(xcc, has_tcga_msi_data)
  df <- data.frame(xcc_tcga_msi[,c(1:7)], xcc_tcga_msi[,"beta" == substring(names(xcc_tcga_msi),1,4)], xcc_tcga_msi[,"cadd" == substring(names(xcc_tcga_msi),1,4)])

  #Variable string
  #paste0(c("point_mutation_rate",names(df)[8:ncol(df)])," + ", collapse='')

  #Use df in model building steps
  #GLM using logistic regression with TCGA MSI-H status as outcome
  glm.log.fit <- glm(tcga_msih ~ point_mutation_rate + beta_ERCC1 + beta_ERCC2 + beta_ERCC3 + beta_ERCC4 + beta_ERCC5 + beta_ERCC6 + beta_EXO1 + beta_LIG1 + beta_MGMT + beta_MLH1 + beta_MLH3 + beta_MSH2 + beta_MSH3 + beta_MSH4 + beta_MSH5 + beta_MSH6 + beta_PCNA + beta_PMS1 + beta_PMS2 + beta_POLD1 + beta_POLD3 + beta_POLE + beta_POLE2 + beta_POLE3 + beta_POLE4 + beta_POLK + beta_RFC1 + beta_RFC2 + beta_RFC3 + beta_RFC4 + beta_RFC5 + beta_RPA1 + beta_RPA2 + beta_RPA4 + cadd_ERCC1 + cadd_ERCC2 + cadd_ERCC3 + cadd_ERCC4 + cadd_ERCC5 + cadd_ERCC6 + cadd_EXO1 + cadd_LIG1 + cadd_MGMT + cadd_MLH1 + cadd_MLH3 + cadd_MSH2 + cadd_MSH3 + cadd_MSH4 + cadd_MSH5 + cadd_MSH6 + cadd_PCNA + cadd_PMS1 + cadd_PMS2 + cadd_POLD1 + cadd_POLD3 + cadd_POLE + cadd_POLE2 + cadd_POLE3 + cadd_POLE4 + cadd_POLK + cadd_RFC1 + cadd_RFC2 + cadd_RFC3 + cadd_RFC4 + cadd_RFC5 + cadd_RPA1 + cadd_RPA2 + cadd_RPA3 + cadd_RPA4, data=df, family=binomial)

  glm.log.fit.stepwise <- stepAIC(glm.log.fit, direction="both", step=1e4)
  table(predict(glm.log.fit.stepwise, newdata=df, type="response")>0.5, df$tcga_msih)

  #GLM using logistic regression with MSIsensor > 4 as outcome
  glm.log.msisensor.fit <- glm((msisensor > 4) ~ point_mutation_rate + beta_ERCC1 + beta_ERCC2 + beta_ERCC3 + beta_ERCC4 + beta_ERCC5 + beta_ERCC6 + beta_EXO1 + beta_LIG1 + beta_MGMT + beta_MLH1 + beta_MLH3 + beta_MSH2 + beta_MSH3 + beta_MSH4 + beta_MSH5 + beta_MSH6 + beta_PCNA + beta_PMS1 + beta_PMS2 + beta_POLD1 + beta_POLD3 + beta_POLE + beta_POLE2 + beta_POLE3 + beta_POLE4 + beta_POLK + beta_RFC1 + beta_RFC2 + beta_RFC3 + beta_RFC4 + beta_RFC5 + beta_RPA1 + beta_RPA2 + beta_RPA4 + cadd_ERCC1 + cadd_ERCC2 + cadd_ERCC3 + cadd_ERCC4 + cadd_ERCC5 + cadd_ERCC6 + cadd_EXO1 + cadd_LIG1 + cadd_MGMT + cadd_MLH1 + cadd_MLH3 + cadd_MSH2 + cadd_MSH3 + cadd_MSH4 + cadd_MSH5 + cadd_MSH6 + cadd_PCNA + cadd_PMS1 + cadd_PMS2 + cadd_POLD1 + cadd_POLD3 + cadd_POLE + cadd_POLE2 + cadd_POLE3 + cadd_POLE4 + cadd_POLK + cadd_RFC1 + cadd_RFC2 + cadd_RFC3 + cadd_RFC4 + cadd_RFC5 + cadd_RPA1 + cadd_RPA2 + cadd_RPA3 + cadd_RPA4, data=df, family=binomial)

  glm.log.msisensor.fit.stepwise <- stepAIC(glm.log.msisensor.fit, direction="both", step=1e4)
  table(predict(glm.log.msisensor.fit.stepwise, newdata=df, type="response")>0.5, df$msisensor>4)

  save(glm.log.fit, glm.log.fit.stepwise, glm.log.msisensor.fit, glm.log.msisensor.fit.stepwise, file="glm_logistic_models.Rdata")

  #Visualize using jittered dot plots
  predicted_tcga_msih <- predict(glm.log.fit.stepwise, newdata=df, type="response")
  predicted_msisensor <- predict(glm.log.msisensor.fit.stepwise, newdata=df, type="response")
  plot_df <- data.frame(predicted_tcga_msih, predicted_msisensor, cancer_type=droplevels(df$cancer_type), tcga_msih=df$tcga_msih, msisensor_high=df$msisensor>4)
  p <- ggplot(plot_df) + aes(x=tcga_msih, y=predicted_tcga_msih, color=cancer_type) + geom_jitter(width=0.3, height=0) + labs(x="TCGA MSI-H", y="Predicted MSI-H probability", title="Predicted MSI-H status", fill="Cancer type")
  pdf("plots/predicted_msih_status.pdf",10,10)
  print(p)
  dev.off()
  p <- ggplot(plot_df) + aes(x=predicted_tcga_msih, fill=tcga_msih) + geom_density(alpha=0.3) + labs(x="Predicted MSI-H probability", y="Density", title="Distribution of MSI-H prediction scores", fill="TCGA MSI-H")
  pdf("plots/distribution_msih_prediction_scores.pdf",10,10)
  print(p)
  dev.off()

  rm(list=c("df", "plot_df", "xcc_tcga_msi", "glm.log.fit", "glm.log.fit.stepwise", "glm.log.msisensor.fit", "glm.log.msisensor.fit.stepwise", "has_tcga_msi_data", "p", "predicted_msisensor", "predicted_tcga_msih"))
}

#TEST and TRAIN
{
  xcc_tcga_msi <- xcc[xcc$tcga_msi %in% c("MSI-H", "MSI-L", "MSS"),]
  train <- sample(c(TRUE,FALSE),size=nrow(xcc_tcga_msi),replace=TRUE)
  test <- !train

  #Refit model using training data, see error in test data
  glm.log.train.refit <- glm(glm.log.fit.stepwise$formula, data=xcc_tcga_msi[train,], family=binomial)

  #Predict outcomes of test data
  predicted <- predict(glm.log.train.refit, newdata=xcc_tcga_msi[test,], type="response")

  #Confusion matrix
  table(predicted > 0.5, xcc_tcga_msi$tcga_msih[test])
}

#Lasso model for variable selection
{
  xcc_tcga_msi <- xcc[xcc$tcga_msi %in% c("MSI-H","MSI-L","MSS"),]
  myX <- as.matrix(xcc_tcga_msi[,substring(names(xcc_tcga_msi),1,5) %in% c("point","beta_","cadd_")])
  myY <- as.numeric(xcc_tcga_msi$tcga_msih)
  variables_in_model <- vector()
  for(i in 1:1000){
    print(i)
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=10, family="binomial", type.measure="class", alpha=0.9)
    variables_in_model <- append(variables_in_model, row.names(coef(cvfit, s="lambda.min"))[coef(cvfit, s="lambda.min")[,1] != 0])
  }

  #FIND THE BEST LAMBDA
  #http://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results/173895#173895
  lambdas = NULL
  for (i in 1:1000)
  {
    print(i)
    fit <- cv.glmnet(myX,myY,nfolds=10, family="binomial", type.measure="class", alpha=0.9)
    errors = data.frame(fit$lambda,fit$cvm)
    lambdas <- rbind(lambdas,errors)
  }
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[,2], list(lambdas$fit.lambda), mean)

  # select the best one
  bestindex = which(lambdas[,2]==min(lambdas[,2]))
  bestlambda = lambdas[bestindex,1]

  # and now run glmnet once more with it
  glmnet.bestLambda.fit <- glmnet(myX,myY,lambda=bestlambda, family="binomial", alpha=0.9)
  best_coef_indicator <- coef(glmnet.bestLambda.fit)[,1] != 0
  coefs_in_best_model <- rownames(coef(glmnet.bestLambda.fit))[best_coef_indicator]

  #Plot to visualize how many times each variable appeared
  vars <- as.matrix(unique(sort(variables_in_model)))
  num <- apply(vars,1,function(y) sum(variables_in_model==y))
  plot_df <- data.frame(vars, num, best_coefs=(vars %in% coefs_in_best_model))
  p <- ggplot(plot_df, aes(x=num, y=reorder(vars,num), color=best_coefs)) + geom_point(size=3) + labs(x="Number of models", y="Model variable", title="Number of models each variable appeared in", color="Included in\n'best' model")
  pdf("plots/variables_in_model.pdf",10,10)
  print(p)
  dev.off()

  model_df <- xcc_tcga_msi[,names(xcc_tcga_msi) %in% c("tcga_msih", coefs_in_best_model)]
  glm.bestmodel.fit <- glm( tcga_msih ~ ., data=model_df, family=binomial)
  summary(glm.bestmodel.fit)

  plot_df <- data.frame(predicted=glm.bestmodel.fit$fitted.values, actual=xcc_tcga_msi$tcga_msih, cancer_type=xcc_tcga_msi$cancer_type)
  p <- ggplot(plot_df) + aes(x=actual, y=predicted, color=cancer_type) + geom_jitter(width=0.5, height=0) + labs(x="TCGA MSI-H status", y="Fitted probabilty of being MSI-H", title="Probabilty of MSI-H", color="Cancer type")
  pdf("plots/probability_of_msih_dots.pdf",10,10)
  print(p)
  dev.off()
  p <- ggplot(plot_df) + aes(x=actual, y=predicted) + geom_violin() + labs(x="TCGA MSI-H status", y="Fitted probabilty of being MSI-H", title="Probabilty of MSI-H")
  pdf("plots/probability_of_msih_violin.pdf",10,10)
  print(p)
  dev.off()

  #Calculate ROC curves for model and MSIsensor
  num_obs = length(msisensor_range)
  model_range <- sort(glm.bestmodel.fit$fitted.values)
  msisensor_range <- sort(xcc_tcga_msi$msisensor)
  model_roc <- matrix(NA,num_obs,3)
  msisensor_roc <- matrix(NA,num_obs,3)
  for(i in 1:num_obs){
    TT = sum(xcc_tcga_msi$tcga_msih & glm.bestmodel.fit$fitted.values >= model_range[i])
    TF = sum(xcc_tcga_msi$tcga_msih & !(glm.bestmodel.fit$fitted.values >= model_range[i]))
    FT = sum(!xcc_tcga_msi$tcga_msih & glm.bestmodel.fit$fitted.values >= model_range[i])
    FF = sum(!xcc_tcga_msi$tcga_msih & !(glm.bestmodel.fit$fitted.values >= model_range[i]))
    model_tab <- matrix(c(TT, FT, TF, FF),2,2)
    TT = sum(xcc_tcga_msi$tcga_msih & xcc_tcga_msi$msisensor > msisensor_range[i])
    TF = sum(xcc_tcga_msi$tcga_msih & !(xcc_tcga_msi$msisensor > msisensor_range[i]))
    FT = sum(!xcc_tcga_msi$tcga_msih & xcc_tcga_msi$msisensor > msisensor_range[i])
    FF = sum(!xcc_tcga_msi$tcga_msih & !(xcc_tcga_msi$msisensor > msisensor_range[i]))
    msisensor_tab <- matrix(c(TT, FT, TF, FF),2,2)
    #TPR, FPR for each type
    model_roc[i,] <- c(model_range[i], model_tab[1,1]/sum(model_tab[1,]), model_tab[2,1]/sum(model_tab[2,]))
    msisensor_roc[i,] <- c(msisensor_range[i], msisensor_tab[1,1]/sum(msisensor_tab[1,]), msisensor_tab[2,1]/sum(msisensor_tab[2,]))
  }
  plot_df <- rbind(data.frame(value = model_roc[,1], tpr = model_roc[,2], fpr = model_roc[,3], method=rep("Model",num_obs)), data.frame(value = msisensor_roc[,1], tpr = msisensor_roc[,2], fpr = msisensor_roc[,3], method=rep("MSIsensor",num_obs)))
  p <- ggplot(plot_df) + aes(x=fpr, y=tpr, color=method) + geom_line(size=1) + geom_abline(intercept=0, slope=1, linetype=2, alpha=0.5) + labs(x="False positive rate", y="True positive rate", title="ROC curve of regression model vs. MSIsensor", color="Method")
  pdf("plots/roc_curve_model_msisensor.pdf",10,10)
  print(p)
  dev.off()
}
