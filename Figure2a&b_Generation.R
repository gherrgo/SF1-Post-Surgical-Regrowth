###############################################################################################
# --------------------------------------------------------------------------------------------#
# Figure 2a-b generation - Kaplan-Meier survival curves
# This script provides survival curve generation as depicted in the publication.
# --------------------------------------------------------------------------------------------#
###############################################################################################
rm(list = ls()) %>% gc()

bioc_pkgs <- c('survival', 'survminer', 'ranger', 'minfi', 'ggplot2',
               'ggfortify', 'gridExtra', 'magrittr', 'gridExtra')
cran_pkgs <- c( 'dplyr', 'readxl', 'readr', 'xlsx')

f (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos='https://cloud.r-project.org')
}
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
install.packages(cran_pkgs, repos='https://cloud.r-project.org', dependencies = TRUE)
lapply(c(bioc_pkgs, cran_pkgs), library, character.only = TRUE)

############################################################
# 1. Necessary materials acquisition & event coding
############################################################
pd <- read.xlsx('./Materials/Supplementary_File_S1.xlsx', sheetName = 'Sheet1')
pd$event <- ifelse(pd$Postoperative.Tumor.Volume < pd$Postoperative.Volume, 1, 0)
pd$event <- pd$Reintervention

###############################################################################################
# 2. Collapse the dataset by 'Record.id' to capture the first event time or the last time point
###############################################################################################
pd_collapsed <- pd %>%
  group_by(Record.id) %>%
  arrange(Months.after.1..postoperative.scan) %>%
  # Get the first event time where event == 1
  mutate(event_time = ifelse(event == 1, Months.after.1..postoperative.scan, NA)) %>%
  # If no event occurred, take the last observation's time
  summarise(event_time = ifelse(any(event == 1), min(event_time, na.rm = TRUE), max(Months.after.1..postoperative.scan)),
            event = ifelse(any(event == 1), 1, 0), # Mark as event if any event happened
            Kcluster_5 = dplyr::first(Kcluster_5)) %>% 
  ungroup()

###############################################################################################
# 3. Kaplan-Meier Survival model fitting & generation of Figure 2a
###############################################################################################
fit <- survfit(Surv(event_time, event) ~ Kcluster_5, data = pd_collapsed)
cluster_levels <- levels(factor(pd_collapsed$Kcluster_5))
num_levels <- length(cluster_levels)

cox_model <- coxph(Surv(event_time, event) ~ factor(Kcluster_5), data = pd_collapsed)
cox_summary <- summary(cox_model)

clusters_present <- sort(unique(pd_collapsed$Kcluster_5))  # Ensure all present clusters are included

# ----- Extract Hazard Ratio (HR) and relevant p-values from Cox summary
hr_values <- exp(cox_summary$coefficients[, "coef"])
lower_ci      <- cox_summary$conf.int[, "lower .95"]
upper_ci      <- cox_summary$conf.int[, "upper .95"]
p_values <- cox_summary$coefficients[, "Pr(>|z|)"]

# ----- Create a full vector for all clusters (initialize with default values)
full_hr_values <- rep(1, length(clusters_present)) # Start with HR = 1 for all
full_lower_ci     <- rep(1, length(clusters_present))  # CI for reference
full_upper_ci     <- rep(1, length(clusters_present))
full_p_values <- rep(1, length(clusters_present))  # Start with p = 1 for all

# ----- Fill in HR and p-values for non-reference clusters
non_ref_clusters <- as.numeric(gsub("factor\\(Kcluster_5\\)", "", rownames(cox_summary$coefficients)))
full_hr_values[non_ref_clusters]    <- hr_values
full_lower_ci[non_ref_clusters]     <- lower_ci
full_upper_ci[non_ref_clusters]     <- upper_ci
full_p_values[non_ref_clusters]     <- p_values

# ----- Add HR = 1, CI = 1, p = NA for the reference group
hr_values  <- c(1, hr_values_raw)
p_values   <- c(NA, p_values_raw)
lower_ci   <- c(1, lower_ci_raw)
upper_ci   <- c(1, upper_ci_raw)

# ----- Create the HR dataframe
hr_df <- data.frame(
  Cluster = clusters_present, HazardRatio = hr_values,
  LowerCI = lower_ci, UpperCI = upper_ci,
  Pvalue = p_values
)

# ----- Generation of the kaplan-meier survival curve plot
km_plot <- ggsurvplot(
  fit,
  data = pd_collapsed,
  risk.table = "nrisk_cumcensor",                         # Include risk table
  risk.table.height = 0.25,
  palette = c("red", "orange", "yellow2", "blue", "green2"),  # Custom colors
  xlab = "Months after surgery",             # X-axis label
  ylab = "Progression-Free Survival",        # Y-axis label
  break.time.by = 24,                        # Time intervals (24 months)
  title = "Progression-Free Survival",  # Title
  legend.title = "Clusters",                 # Legend title
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")
)
km_plot
# ----- Extraction of the main plot and risk table separately
main_plot <- km_plot$plot       # Extract the survival plot
risk_table_plot <- km_plot$table  # Extract the risk table plot (with nrisk_cumcensor)

# ----- Add annotations for hazard ratios and p-values to the survival plot
main_plot_with_annotations <- main_plot + 
  annotate("text", x = 130, y = 0.8, 
           label = paste("HR (Cluster 1):", 
                         round(hr_df$HazardRatio[1], 2),
                         "[95% CI:",
                         round(hr_df$LowerCI[1], 2), "-", round(hr_df$UpperCI[1], 2), "]",
                         ", p =", round(hr_df$Pvalue[1], 3)), 
           color = "black", size = 4) +
  
  annotate("text", x = 130, y = 0.75, 
           label = paste("HR (Cluster 2):", 
                         round(hr_df$HazardRatio[2], 2),
                         "[95% CI:",
                         round(hr_df$LowerCI[2], 2), "-", round(hr_df$UpperCI[2], 2), "]",
                         ", p =", round(hr_df$Pvalue[2], 3)), 
           color = "black", size = 4) +
  
  annotate("text", x = 130, y = 0.7, 
           label = paste("HR (Cluster 3):", 
                         round(hr_df$HazardRatio[3], 2),
                         "[95% CI:",
                         round(hr_df$LowerCI[3], 2), "-", round(hr_df$UpperCI[3], 2), "]",
                         ", p =", round(hr_df$Pvalue[3], 3)), 
           color = "black", size = 4) +
  
  annotate("text", x = 130, y = 0.65, 
           label = paste("HR (Cluster 4):", 
                         round(hr_df$HazardRatio[4], 2),
                         "[95% CI:",
                         round(hr_df$LowerCI[4], 2), "-", round(hr_df$UpperCI[4], 2), "]",
                         ", p =", round(hr_df$Pvalue[4], 3)), 
           color = "black", size = 4) +
  
  annotate("text", x = 130, y = 0.6, 
           label = paste("HR (Cluster 5):", 
                         round(hr_df$HazardRatio[5], 2),
                         "[95% CI:",
                         round(hr_df$LowerCI[5], 2), "-", round(hr_df$UpperCI[5], 2), "]",
                         ", p =", round(hr_df$Pvalue[5], 3)), 
           color = "black", size = 4)
# ----- Use gridExtra to combine the main plot with annotations and the risk table plot
grid.arrange(main_plot_with_annotations, risk_table_plot, ncol = 1, heights = c(4, 1.5))

pdf('./Figures/Figure2a.pdf', height = 8, width = 10)
grid.arrange(main_plot_with_annotations, risk_table_plot, ncol = 1, heights = c(4, 1.2))
dev.off()

###############################################################################################
# 3. Exclusion of k4 and generation of Figure 2b (additional kaplan-meier survival curve)
###############################################################################################

# ----- Exclusion of k4 
pd_collapsed_filtered <- pd_collapsed %>%
  filter(Kcluster_5 != 4)

# ----- Get cluster levels after filtering
clusters_present <- levels(factor(pd_collapsed_filtered$Kcluster_5))
legend_labels <- paste("Cluster", clusters_present)

# ----- Adjust color palette (based on number of remaining clusters)
palette_used <- c("red", "orange", "yellow2", "green2")[1:length(clusters_present)]

# ----- Fit survival and Cox models
surv_fit <- survfit(Surv(event_time, event) ~ Kcluster_5, data = pd_collapsed_filtered)
cox_model <- coxph(Surv(event_time, event) ~ factor(Kcluster_5), data = pd_collapsed_filtered)
cox_summary <- summary(cox_model)

# ----- Extract HR, CI, and p-values (excluding reference)
hr_values_raw   <- exp(cox_summary$coefficients[, "coef"])
lower_ci_raw    <- cox_summary$conf.int[, "lower .95"]
upper_ci_raw    <- cox_summary$conf.int[, "upper .95"]
p_values_raw    <- cox_summary$coefficients[, "Pr(>|z|)"]

# ----- Add HR = 1, CI = 1, p = NA for reference group
hr_values  <- c(1, hr_values_raw)
lower_ci   <- c(1, lower_ci_raw)
upper_ci   <- c(1, upper_ci_raw)
p_values   <- c(NA, p_values_raw)

# ----- Ensure all vectors match the number of clusters
hr_values  <- hr_values[seq_along(clusters_present)]
lower_ci   <- lower_ci[seq_along(clusters_present)]
upper_ci   <- upper_ci[seq_along(clusters_present)]
p_values   <- p_values[seq_along(clusters_present)]

# ----- Create the HR dataframe
hr_df <- data.frame(
  Cluster = clusters_present, HazardRatio = hr_values,
  LowerCI = lower_ci, UpperCI = upper_ci,
  Pvalue = p_values
)

# ----- Create the Kaplan-Meier plot
km_plot <- ggsurvplot(
  surv_fit,
  data = pd_collapsed_filtered,  
  risk.table = "nrisk_cumcensor",  
  risk.table.height = 0.25,
  palette = palette_used,  
  xlab = "Months after surgery",  
  ylab = "Progression-Free Survival",  
  break.time.by = 24,  
  title = "Progression-Free Survival",  
  subtitle = "SF1 positive PitNETs",  
  legend.title = "Clusters",  
  legend.labs = legend_labels
)

# ----- Extract main plot and risk table
main_plot <- km_plot$plot
risk_table_plot <- km_plot$table

# ----- Add annotations for HR + 95% CI and p-values
main_plot_with_annotations <- main_plot
for (i in seq_along(clusters_present)) {
  main_plot_with_annotations <- main_plot_with_annotations + 
    annotate(
      "text", x = 130, y = 0.8 - (i - 1) * 0.05, 
      label = paste0(
        "HR (", legend_labels[i], "): ",
        round(hr_df$HazardRatio[i], 2),
        " [95% CI: ",
        round(hr_df$LowerCI[i], 2), "–", round(hr_df$UpperCI[i], 2), "]",
        ifelse(is.na(hr_df$Pvalue[i]), "", paste0(", p = ", round(hr_df$Pvalue[i], 3)))
      ),
      color = "black", size = 4
    )
}

# ----- Arrange the plots
final_plot <- grid.arrange(main_plot_with_annotations, risk_table_plot, ncol = 1, heights = c(4, 1))

# ----- Save to PDF
pdf('./Figures/Figure2b.pdf', height = 8, width = 10)
grid.arrange(main_plot_with_annotations, risk_table_plot, ncol = 1, heights = c(4, 1))
dev.off()
