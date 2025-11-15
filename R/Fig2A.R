library(ComplexHeatmap)

title <- "glioma_baseline_responder_Vs_nonresponder_logfc_1_pval_0.05_Top_50"

load(paste0(title, "_data.Rdata"))

top_annot <- HeatmapAnnotation(df = res_data[["FinalPheno"]], col = list(
  group = c("nonresponder" = "#FC8D62", "responder" = "#66C2A5",  
            "NA" = "#CCCCCC")),
  gp = gpar(col = "white", lwd = 1),
  border = FALSE,
  annotation_legend_param = list(border = TRUE))

hm <- Heatmap(res_data[["FinalMat"]],
              name = "scale",
              row_title_side = "right",
              top_annotation = top_annot,
              show_row_names = F,
              row_names_gp = gpar(fontsize = 8),
              show_column_names = F, 
              show_heatmap_legend = TRUE,
              cluster_rows = T,
              cluster_columns = T,
              row_split = res_data[["resSig"]]$Dir,
              column_split = res_data[["FinalPheno"]]$group,
              cluster_row_slices = FALSE)

pdf(paste0(title, "_HeatMap.pdf"), 
    height = 7, width = 9)
print (hm)
dev.off()
