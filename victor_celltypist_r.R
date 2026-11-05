# Function for preprocessing celltypist output (probability_matrix) for Immune_All_Low.pkl V2
# The input should be the celltypist output: "probability_matrix.csv" read into R, where the first column are cell barcodes and the next 98 columns are cell probabilities

r_celltypist_immune <- function(prob_matrix = NULL, prob_threshold = 90){
  
  # Require packs
  require(stringr)
  
  # Lineages
  blin <- c("B cells", "Follicular B cells", "Proliferative germinal center B cells",
                    "Germinal center B cells", "Memory B cells", "Age-associated B cells",
                    "Naive B cells", "Transitional B cells", "Large pre-B cells",
                    "Small pre-B cells", "Pre-pro-B cells", "Pro-B cells",
                    "Cycling B cells", "Plasma cells", "Plasmablasts")
  tlin <- c("Cycling gamma-delta T cells", "Cycling T cells", "Double-negative thymocytes",
                            "Double-positive thymocytes", "ETP", "CD8a/a", "CD8a/b(entry)", "Follicular helper T cells",
                            "CRTAM+ gamma-delta T cells", "gamma-delta T cells", "MAIT cells", "Memory CD4+ cytotoxic T cells",
                            "NKT cells", "Regulatory T cells", "T(agonist)", "Tcm/Naive cytotoxic T cells",
                            "Tcm/Naive helper T cells", "Tem/Trm cytotoxic T cells", "Tem/Temra cytotoxic T cells",
                            "Trm cytotoxic T cells", "Tem/Effector helper T cells", "Tem/Effector helper T cells PD1+",
                            "Treg(diff)", "Type 1 helper T cells", "Type 17 helper T cells")
  
  # Probability
  prob <- prob_threshold
  
  # Process data
  df <- data.frame(
    ID = prob_matrix[,1]
  )

  # Incorporate columns for celltypist max probability celltype
  df[,"celltypist_max"] = apply(prob_matrix[,-1],1,function(y){
    colnames(prob_matrix)[-1][which.max(y)]
  })
  
  df[,"celltypist_max_prob"] = apply(prob_matrix[,-1],1,function(y){
    y[which.max(y)]
  })
  
  # Incorporate colums for celltypist celltypes with probability > prob
  df[,paste0("celltypist_prob",prob)] = apply(prob_matrix[,-1],1,function(y){
    paste(colnames(prob_matrix)[-1][y > (prob/100)],collapse="|")
  })
  
  df[,paste0("celltypist_prob",prob)][df[,paste0("celltypist_prob",prob)] == ""] <- "Unassigned"
  
  # Incorporate lineage and mixed columns for celltypist with probability > prob
  
  # Split each celltype label
  splits <- str_split(df[,paste0("celltypist_prob",prob)], pattern = "\\|")
  names(splits) <- df$ID
  
  # Define mixed cells
  mixing <- sapply(splits,function(y){
    if(length(y) <= 1){
      "No_mixed"
    }else if(length(y) > 1){
      if(sum(y %in% blin) == length(y)){
        "Mixed_B"
      }else if(sum(y %in% tlin) == length(y)){
        "Mixed_T"
      }else{
        "Mixed_other"
      }
    }
  })
  
  # Define Bcell or Tcell lineage
  lineage <- case_when(
    df[,paste0("celltypist_prob",prob)] %in% blin | mixing == "Mixed_B" ~ "B_lineage",
    df[,paste0("celltypist_prob",prob)] %in% tlin | mixing == "Mixed_T" ~ "T_lineage",
    TRUE ~ "Other"
  )
  
  # Incorporate
  df[,paste0("celltypist_",prob,"_mixing")] <- mixing
  df[,paste0("celltypist_",prob,"_lineage")] <- lineage
  
  # Final annotation
  df[,paste0("celltypist_",prob,"_final")] <- df[,paste0("celltypist_prob",prob)]
  df[,paste0("celltypist_",prob,"_final")][!grepl("No_mixed",df[,paste0("celltypist_",prob,"_mixing")])] <- df[,paste0("celltypist_",prob,"_mixing")][!grepl("No_mixed",df[,paste0("celltypist_",prob,"_mixing")])]
  
  # Specify unassigned cells
  df[,paste0("celltypist_",prob,"_mixing")][df[,paste0("celltypist_prob",prob)] == "Unassigned"] <- "Unassigned"
  df[,paste0("celltypist_",prob,"_lineage")][df[,paste0("celltypist_prob",prob)]== "Unassigned"] <- "Unassigned"
  df[,paste0("celltypist_",prob,"_final")][df[,paste0("celltypist_prob",prob)] == "Unassigned"] <- "Unassigned"
  

  # Output
  return(df)
}
