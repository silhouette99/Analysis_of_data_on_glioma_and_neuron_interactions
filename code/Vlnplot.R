library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)

individual_vln <- function(obj,
                           features,
                           group,
                           assay = 'RNA',
                           my_comparation = NULL,
                           layer = 'data',
                           raster = NULL,
                           method = 'wilcox.test',
                           ylims = NULL,
                           pt.size = 1,
                           alpha = 0.1,
                           dis_num = 'Mean',
                           lines = F) {
  if (length(intersect(group, colnames(obj@meta.data))) == 0) {
    stop('Missing this group')
  }
  
  if (length(intersect(features, c(
    colnames(obj@meta.data), rownames(GetAssay(obj, assay = assay))
  ))) == 0) {
    stop('Missing this feature')
  }
  
  if (length(intersect(features, c(colnames(obj@meta.data)))) == 0) {
    obj@meta.data[, gsub(pattern = '-',
                         replacement = '_',
                         paste0(features, '_RNA'))] <-
      GetAssayData(obj, layer = layer, assay = assay)[features, ]
    features <-
      gsub(pattern = '-',
           replacement = '_',
           paste0(features, '_RNA'))
  }
  
  if (is.null(ylims) | length(ylims) != 2 | !is.numeric(ylims)) {
    ylim_ = NULL
  } else{
    ylim_ <- ylim(ylims)
  }
  
  if (!is.null(my_comparation)) {
    stat_com <-
      stat_compare_means(comparisons = my_comparation, method = method)
  } else{
    stat_com <- NULL
  }
  
  
  if(is.null(dis_num)){
    labes <- NULL
  }else{
    if(dis_num == 'Mean'){
      funcs <- mean
    }else if(dis_num == 'Median'){
      funcs <- median
    }
    
    num_df <- obj@meta.data %>%
      dplyr::group_by(!!rlang::sym(group)) %>%
      dplyr::summarise(
        num = funcs(!!rlang::sym(features)),
        max_val = max(!!rlang::sym(features), na.rm = TRUE),
        q75 = quantile(!!rlang::sym(features),probs = 0.75, na.rm = TRUE)
      ) %>%
      dplyr::rename(group = !!rlang::sym(group))
    
    labes <- geom_label(
      data = num_df,
      aes(x = group, y = 1.5*q75, label = round(num, 3)),
      fill = 'grey98',
      inherit.aes = FALSE,
      vjust = 0,
      size = 4
    )
    
  }
  
  
  
  
  p <- VlnPlot(
    obj,
    features = features,
    assay = assay,
    alpha = alpha,
    raster = raster,
    group.by = group,
    pt.size = pt.size
  ) +
    NoLegend() + theme_pubclean() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        colour = 'black',
        size = 12
      ),
      axis.text.y = element_text(
        hjust = 1,
        vjust = 1,
        colour = 'black',
        size = 12
      )
    ) +
    stat_com +
    ylim_ + geom_boxplot(
      data = obj@meta.data,
      width = 0.1,
      aes(
        x = !!sym(group),
        y = !!sym(features),
        fill = !!sym(group)
      ),
      color = 'black',
      fill = 'white',
      # 显式定义 fill
      inherit.aes = FALSE  # 禁止继承 VlnPlot 的 aes 映射
    ) + ylab(label = 'Normalized Exp') + labes
  
  if(lines == T){
    median_data <- obj@meta.data %>%
      group_by(!!rlang::sym(group)) %>%  # 替换为您的分组变量
      summarise(mean_value = funcs(!!rlang::sym(features), na.rm = TRUE))
    
    p <- p + geom_smooth(data = median_data,
                    aes(x = as.numeric(factor(!!rlang::sym(group))), 
                        y = mean_value,fill = NULL),
                    method = "loess",  # 使用局部回归拟合曲线
                    formula = y ~ x,
                    se = FALSE,        # 不显示置信区间
                    color = "black",
                    linewidth = 1.2) +
      theme_minimal()
  }
  
  
  
  return(p)
}
