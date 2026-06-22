# =====================================================================
# run_cluster_analysis()
#
# Wraps the original clustering script into a single reusable function.
# Performs:
#   1. (optional) within-subject scaling
#   2. clustering via "regular" (raw variables) or "umap" (UMAP + kmeans)
#   3. an elbow plot (fviz_nbclust) and a dendrogram (fviz_dend)
#   4. group-wise comparison plots (each cluster vs. the rest), with
#      t-test p-values and Cohen's d annotated on the bars
#   5. a combined grid of all group plots
#
# Returns a list containing the resulting dataset AND all figures, so
# nothing is silently left in the global environment.
# =====================================================================

required_pkgs <- c("dplyr", "ggplot2", "ggpubr", "gridExtra", "reshape2",
                   "Rmisc", "psych", "factoextra", "umap")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(paste0("Package '", p, "' is required but not installed. ",
                "Install it with install.packages('", p, "')"))
  }
}))

library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(Rmisc)
library(psych)
library(factoextra)

#' Run cluster analysis on a set of variables
#'
#' @param data         data.frame containing id_var and the clustering variables
#'                      (e.g. already filtered to baseline/T1, as in the original script)
#' @param vars         character vector of column names to cluster on
#'                      (e.g. c("d_proc","d_att_ctrl","d_mem","d_ver_wm","d_vis_wm","d_visuom"))
#' @param id_var       name of the subject/id column (default "subj")
#' @param scaling      "unscaled" (default) or "scaled" (row-centers each subject's scores)
#' @param method       "umap" (default) or "regular" (cluster directly on `vars`)
#' @param k            number of clusters, between 2 and 6 (default 4)
#' @param seed         random seed for reproducibility (default 123)
#' @param colors       fill colors used for each cluster's bar plot
#' @param save_figures logical; if TRUE, also writes .pdf figures to `figures_dir`
#' @param figures_dir  directory for saved figures (created if it doesn't exist)
#' @param save_data    logical; if TRUE, also writes the resulting dataset as .csv to `data_dir`
#' @param data_dir     directory for the saved .csv (created if it doesn't exist)
#' @param file_tag     prefix used in output filenames
#'
#' @return a list with:
#'   data          - data.frame: id_var + vars + groups (cluster assignment)
#'   cluster_model - the kmeans model object
#'   nbclust_plot  - elbow plot (ggplot) used to help choose k
#'   dendrogram    - hierarchical clustering dendrogram (ggplot)
#'   group_plots   - named list of per-cluster bar plots ("group_1", "group_2", ...)
#'   combined_plot - grid-arranged combination of all group plots (a grob; draw with grid::grid.draw())
#'   stats         - list with p_values and cohen_d per group/variable
#'
#' @examples
#' \dontrun{
#'   dat <- read.csv2("et_domains_normed_peds_COGCOR_age_2026-05-12.csv",
#'                     stringsAsFactors = FALSE, sep = ";")
#'   dat$test <- lubridate::ymd(dat$test)
#'   data_T1 <- dat[dat$group == 1, ]
#'
#'   res <- run_cluster_analysis(
#'     data    = data_T1,
#'     vars    = c("d_proc","d_att_ctrl","d_mem","d_ver_wm","d_vis_wm","d_visuom"),
#'     scaling = "unscaled",
#'     method  = "umap",
#'     k       = 4,
#'     save_figures = TRUE,
#'     save_data    = TRUE
#'   )
#'
#'   res$data           # clustered dataset (subj + domains + groups)
#'   res$nbclust_plot    # elbow plot
#'   res$dendrogram      # dendrogram
#'   res$group_plots$group_1   # bar plot for cluster 1 vs rest
#'   grid::grid.draw(res$combined_plot)  # all group plots together
#' }
run_cluster_analysis <- function(data,
                                 vars,
                                 id_var = "subj",
                                 scaling = c("unscaled", "scaled"),
                                 method = c("umap", "regular"),
                                 k = 4,
                                 seed = 123,
                                 colors = c("skyblue", "palegreen", "orange",
                                            "tomato", "purple", "gold"),
                                 save_figures = FALSE,
                                 figures_dir = "figures",
                                 save_data = FALSE,
                                 data_dir = "databases",
                                 file_tag = "cluster_analysis") {
  
  scaling <- match.arg(scaling)
  method  <- match.arg(method)
  
  if (!all(vars %in% colnames(data))) {
    stop("Some `vars` are not present in `data`: ",
         paste(setdiff(vars, colnames(data)), collapse = ", "))
  }
  if (!id_var %in% colnames(data)) {
    stop("`id_var` ('", id_var, "') is not present in `data`.")
  }
  if (k < 2 || k > 6) {
    stop("k must be between 2 and 6 (the combined-plot layout only supports this range).")
  }
  
  subj <- data[[id_var]]
  data_sel <- dplyr::select(data, dplyr::all_of(vars))
  
  # ---- optional within-subject scaling ----
  if (scaling == "scaled") {
    data_sel <- data_sel - rowMeans(data_sel)
  }
  
  # ============================ CLUSTERING ============================
  if (method == "regular") {
    
    set.seed(seed)
    nbclust_plot <- factoextra::fviz_nbclust(data_sel, kmeans, method = "wss")
    
    set.seed(seed)
    hc <- data_sel %>%
      scale() %>%
      dist(method = "manhattan") %>%
      hclust(method = "ward.D2")
    
    set.seed(seed)
    dend_plot <- factoextra::fviz_dend(hc, k = k, cex = 0.5,
                                       color_labels_by_k = TRUE, rect = TRUE)
    
    set.seed(seed)
    clustering_model <- kmeans(as.matrix(data_sel), k, nstart = 100)
    data_sel$groups <- clustering_model$cluster
    
  } else { # method == "umap"
    
    set.seed(seed)
    umap_result <- umap::umap(data_sel, config = umap::umap.defaults)
    umap_embedding <- umap_result$layout
    
    set.seed(seed)
    nbclust_plot <- factoextra::fviz_nbclust(umap_embedding, kmeans, method = "wss")
    
    set.seed(seed)
    hc <- umap_embedding %>%
      scale() %>%
      dist(method = "manhattan") %>%
      hclust(method = "ward.D2")
    
    set.seed(seed)
    dend_plot <- factoextra::fviz_dend(hc, k = k, cex = 0.5,
                                       color_labels_by_k = TRUE, rect = TRUE)
    
    set.seed(seed)
    clustering_model <- kmeans(umap_embedding, k, nstart = 100)
    data_sel$groups <- clustering_model$cluster
  }
  
  # ======================= GROUP-WISE COMPARISON =======================
  group_plots <- list()
  stats_p <- list()
  stats_d <- list()
  
  if (save_figures && !dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
  
  for (z in seq_len(k)) {
    
    grp_data <- data_sel
    grp_data$groups <- as.numeric(grp_data$groups == z)
    
    if (sum(grp_data$groups) <= 1) next  # need >1 observation in the target cluster
    
    data_long <- reshape2::melt(grp_data, id = "groups")
    data_long_m <- Rmisc::summarySE(data_long, measurevar = "value",
                                    groupvars = c("variable", "groups"))
    
    p_vals <- sapply(vars, function(v) {
      round(t.test(grp_data[grp_data$groups == 0, v],
                   grp_data[grp_data$groups == 1, v],
                   alternative = "two.sided")$p.value, 3)
    })
    
    d_vals <- sapply(vars, function(v) {
      round(psych::cohen.d(grp_data[[v]], group = grp_data$groups)$cohen.d[2], 2)
    })
    
    stats_p[[paste0("group_", z)]] <- p_vals
    stats_d[[paste0("group_", z)]] <- d_vals
    
    p_coding <- ifelse(p_vals < .001, "***",
                       ifelse(p_vals < .01,  "**",
                              ifelse(p_vals < .05,  "*", "")))
    
    labels <- paste0("d = ", d_vals)
    col_z <- colors[((z - 1) %% length(colors)) + 1]
    
    plot <- ggplot(data_long_m, aes(x = variable, y = value, fill = factor(groups))) +
      geom_errorbar(aes(ymin = value - se, ymax = value + se, width = 0.5),
                    position = position_dodge(width = 0.9), colour = "grey") +
      geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
      scale_fill_manual(values = c("grey", col_z)) +
      ggtitle(paste0("Group ", z, " (n = ", sum(grp_data$groups == 1), ")")) +
      xlab("Domain") + ylab("Performance (z-score)") +
      annotate("text", x = seq_along(vars), y = -0.5, label = labels, size = 2) +
      annotate("text", x = seq_along(vars), y = 0, label = p_coding, size = 5) +
      theme_pubr() +
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 9),
            plot.title = element_text(hjust = 0.5),
            legend.title = element_blank()) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2))
    
    group_plots[[paste0("group_", z)]] <- plot
    
    if (save_figures) {
      grDevices::cairo_pdf(file.path(figures_dir, paste0("plot_group_", z, ".pdf")))
      print(plot)
      grDevices::dev.off()
    }
  }
  
  # ---- combined grid of all group plots ----
  nrow_lookup <- c(`2` = 2, `3` = 2, `4` = 2, `5` = 3, `6` = 3)
  combined_plot <- gridExtra::arrangeGrob(grobs = group_plots,
                                          nrow = nrow_lookup[as.character(k)])
  
  if (save_figures) {
    grDevices::cairo_pdf(file.path(figures_dir,
                                   paste0("plot_", scaling, "_", method, "_clusters_", k, "_all.pdf")))
    grid::grid.draw(combined_plot)
    grDevices::dev.off()
  }
  
  # ====================== FINAL DATASET ======================
  result_data <- data.frame(setNames(list(subj), id_var), data_sel,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  if (save_data) {
    if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
    write.csv2(result_data,
               file.path(data_dir, paste0(file_tag, "_", scaling, "_clusters_", k, "_all.csv")),
               row.names = FALSE)
  }
  
  list(
    data          = result_data,
    cluster_model = clustering_model,
    nbclust_plot  = nbclust_plot,
    dendrogram    = dend_plot,
    group_plots   = group_plots,
    combined_plot = combined_plot,
    stats         = list(p_values = stats_p, cohen_d = stats_d)
  )
}

# =====================================================================
# Example usage (mirrors the original script's workflow)
# =====================================================================
#
# library(lubridate)
#
# dat <- read.csv2("et_domains_normed_peds_COGCOR_age_2026-05-12.csv",
#                   stringsAsFactors = FALSE, sep = ";")
# dat$test <- ymd(dat$test)
# data_T1 <- dat[dat$group == 1, ]
# write.csv2(data_T1, "COGCOR_BASELINE_DATA_2026.csv", row.names = FALSE)
#
# domain_vars <- c("d_proc", "d_att_ctrl", "d_mem", "d_ver_wm", "d_vis_wm", "d_visuom")
#
# res <- run_cluster_analysis(
#   data         = data_T1,
#   vars         = domain_vars,
#   id_var       = "subj",
#   scaling      = "unscaled",
#   method       = "umap",
#   k            = 4,
#   save_figures = TRUE,
#   figures_dir  = "figures",
#   save_data    = TRUE,
#   data_dir     = "databases",
#   file_tag     = "COGCO_ET_DATABASE"
# )
#
# # the clustered dataset:
# head(res$data)
#
# # the figures:
# res$nbclust_plot
# res$dendrogram
# res$group_plots$group_1
# grid::grid.draw(res$combined_plot)