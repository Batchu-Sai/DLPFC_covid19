# Modify the functions form the SptialLIBD package to remove some superfluous labeling

vis_clus_p2 <- function (spe, d, clustervar, sampleid, colors, spatial, title, 
                         image_id = "lowres", alpha = 1, point_size = 2) 
{
  pxl_row_in_fullres <- pxl_col_in_fullres <- key <- NULL
  if (clustervar %in% c("layer_guess", "layer_guess_reordered", 
                        "layer_guess_reordered_short", "spatialLIBD")) {
    title <- gsub(clustervar, "LIBD Layers", title)
  }
  img <- SpatialExperiment::imgRaster(spe, sample_id = sampleid, 
                                      image_id = image_id)
  p <- ggplot(d, aes(x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, 
                                                                              sample_id = sampleid, image_id = image_id), y = pxl_row_in_fullres * 
                       SpatialExperiment::scaleFactors(spe, sample_id = sampleid, 
                                                       image_id = image_id), fill = factor(!!sym(clustervar)), 
                     key = key))
  if (spatial) {
    grob <- grid::rasterGrob(img, width = grid::unit(1, "npc"), 
                             height = grid::unit(1, "npc"))
    p <- p + geom_spatial(data = tibble::tibble(grob = list(grob)), 
                          aes(grob = grob), x = 0.5, y = 0.5)
  }
  p <- p + geom_point(shape = 21, size = point_size, stroke = 0, 
                      colour = "transparent", alpha = alpha) + coord_cartesian(expand = FALSE) + 
    scale_fill_manual(values = colors) + xlim(0, ncol(img)) + 
    ylim(nrow(img), 0) + xlab("") + ylab("") + labs(fill = NULL) + 
    guides(fill = guide_legend(override.aes = list(size = 8)),
           color = guide_legend(override.aes = list(size = 8)),
    ) + 
    ggtitle(title) + 
    theme_set(theme_bw(base_size = 20)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.text = element_text(size = 8),
          #legend.position="none",
          #plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  return(p)
}

vis_clus2 <- function (spe, sampleid, clustervar, colors = c("#b2df8a", "#e41a1c", 
                                                             "#377eb8", "#4daf4a", "#ff7f00", "gold", "#a65628", "#999999", 
                                                             "black", "grey", "white", "purple"), spatial = TRUE, image_id = "lowres", 
                       alpha = 1, point_size = 2, ...) 
{
  spe_sub <- spe[, spe$sample_id == sampleid]
  d <- as.data.frame(cbind(colData(spe_sub), SpatialExperiment::spatialCoords(spe_sub)), 
                     optional = TRUE)
  vis_clus_p2(spe = spe_sub, d = d, 
              clustervar = clustervar, 
              sampleid = sampleid, 
              spatial = spatial, 
              title = '', 
              colors = get_colors(colors, d[, clustervar]), 
              image_id = image_id, 
              alpha = alpha, 
              point_size = point_size)
}

vis_gene2 <- function (spe, sampleid, geneid = "SCGB2A2; ENSG00000110484", 
                       spatial = TRUE, assayname = "logcounts", minCount = 0, viridis = TRUE, 
                       image_id = "lowres", alpha = 1, cont_colors = if (viridis) viridisLite::viridis(21) else c("aquamarine4", 
                                                                                                                  "springgreen", "goldenrod", "red"), point_size = 2, ...) 
{
  spe_sub <- spe[, spe$sample_id == sampleid]
  d <- as.data.frame(cbind(colData(spe_sub), SpatialExperiment::spatialCoords(spe_sub)), 
                     optional = TRUE)
  if (geneid %in% colnames(colData(spe_sub))) {
    d$COUNT <- colData(spe_sub)[[geneid]]
  }
  else if (geneid %in% rowData(spe_sub)$gene_search) {
    d$COUNT <- assays(spe_sub)[[assayname]][which(rowData(spe_sub)$gene_search == 
                                                    geneid), ]
  }
  else if (geneid %in% rownames(spe_sub)) {
    d$COUNT <- assays(spe_sub)[[assayname]][which(rownames(spe_sub) == 
                                                    geneid), ]
  }
  else {
    stop("Could not find the 'geneid' ", geneid, call. = FALSE)
  }
  d$COUNT[d$COUNT <= minCount] <- NA
  geneid2 <- gsub(";.*", "", geneid)
  p <- vis_gene_p(spe = spe_sub, d = d, 
                  sampleid = sampleid, 
                  spatial = spatial, 
                  title = paste('', ...), 
                  viridis = viridis, 
                  image_id = image_id, 
                  alpha = alpha, 
                  cont_colors = cont_colors, 
                  point_size = point_size)
  
}
