sesameQC_calcStats <- function(sdf, funs = NULL) {
  if (is.null(funs)) {
    funs <- c(
      sesameQC_calcStats_detection,
      sesameQC_calcStats_intensity,
      sesameQC_calcStats_numProbes,
      sesameQC_calcStats_channel,
      sesameQC_calcStats_dyeBias,
      sesameQC_calcStats_betas)
  }
  if (!is(funs, "list")) { funs <- c(funs) }
  qc <- new("sesameQC")
  
  for (func in funs) {
    if (is.character(func)) {
      func <- get(paste0("sesameQC_calcStats_",func))
      stopifnot(is(func, "function"))
    }
    qc <- func(sdf, qc = qc)
  }
  qc
}



sesameQC_calcStats_detection <- function(sdf, qc = NULL) {
  
  g1 <- .setGroup_detection()
  group_nm <- names(g1)[1]
  if (is.null(qc)) { s <- list(); g <- list()
  } else { s <- qc@stat; g <- qc@group }
  if (group_nm %in% names(g)) { return(qc); }
  g[[group_nm]] <- g1[[group_nm]]
  
  pvals0 <- sesame::pOOBAH(sdf, return.pval = TRUE)
  pvals <- na.omit(pvals0)
  s$num_dtna <- sum(is.na(pvals0))
  s$frac_dtna <- s$num_dtna / length(pvals0)
  s$num_dt <- sum(pvals <= 0.05)
  s$frac_dt <- s$num_dt / length(pvals)
  idx_mk <- !is.na(pvals0) & !sdf$mask
  s$num_dt_mk <- sum(pvals0[idx_mk] <= 0.05)
  s$frac_dt_mk <- s$num_dt_mk / sum(idx_mk)
  for (pt in c('cg','ch','rs')) {
    p1 <- pvals[grep(paste0('^', pt), names(pvals))]
    s[[paste0('num_dt_', pt)]] <- sum(p1 <= 0.05)
    s[[paste0('frac_dt_', pt)]] <- sum(p1 <= 0.05) / length(p1)
  }
  new("sesameQC", stat=s, group=g)
}



sesameQC_calcStats_intensity <- function(sdf, qc = NULL) {
  
  g1 <- .setGroup_intensity()
  group_nm <- names(g1)[1]
  if (is.null(qc)) { s <- list(); g <- list()
  } else { s <- qc@stat; g <- qc@group }
  if (group_nm %in% names(g)) { return(qc); }
  g[[group_nm]] <- g1[[group_nm]]
  
  dG <- InfIG(sdf); dR <- InfIR(sdf); d2 <- InfII(sdf)
  s$mean_intensity <- sesame::meanIntensity(sdf) # excluding type-I out-of-band
  s$mean_intensity_MU <- mean(sesame::totalIntensities(sdf), na.rm=TRUE) # M + U
  s$mean_ii <- mean(c(d2$UG,d2$UR), na.rm = TRUE)
  s$mean_inb_grn <- mean(c(dG$MG, dG$UG), na.rm = TRUE)
  s$mean_inb_red <- mean(c(dR$MR, dR$UR), na.rm = TRUE)
  s$mean_oob_grn <- mean(c(dR$MG, dR$UG), na.rm = TRUE)
  s$mean_oob_red <- mean(c(dG$MR, dG$UR), na.rm = TRUE)
  mu <- sesame::signalMU(sdf)
  s$na_intensity_M <- sum(is.na(mu$M))
  s$na_intensity_U <- sum(is.na(mu$U))
  s$na_intensity_ig <- sum(is.na(c(dG$MG, dG$MR, dG$UG, dG$UR)))
  s$na_intensity_ir <- sum(is.na(c(dR$MG, dR$MR, dR$UG, dR$UR)))
  s$na_intensity_ii <- sum(is.na(c(d2$UG, d2$UR)))
  
  new("sesameQC", stat=s, group=g)
}







sesameQC_calcStats_numProbes <- function(sdf, qc = NULL) {
  
  g1 <- .setGroup_numProbes()
  group_nm <- names(g1)[1]
  if (is.null(qc)) { s <- list(); g <- list()
  } else { s <- qc@stat; g <- qc@group }
  if (group_nm %in% names(g)) { return(qc); }
  g[[group_nm]] <- g1[[group_nm]]
  
  s$num_probes <- nrow(sdf)
  s$num_probes_II <- nrow(InfII(sdf))
  s$num_probes_IR <- nrow(InfIR(sdf))
  s$num_probes_IG <- nrow(InfIG(sdf))
  s$num_probes_cg <- sum(startsWith(sdf$Probe_ID,"cg"))
  s$num_probes_ch <- sum(startsWith(sdf$Probe_ID,"ch"))
  s$num_probes_rs <- sum(startsWith(sdf$Probe_ID,"rs"))
  new("sesameQC", stat=s, group=g)
}



sesameQC_calcStats_channel <- function(sdf, qc = NULL) {
  
  g1 <- .setGroup_channel()
  group_nm <- names(g1)[1]
  if (is.null(qc)) { s <- list(); g <- list()
  } else { s <- qc@stat; g <- qc@group }
  if (group_nm %in% names(g)) { return(qc); }
  g[[group_nm]] <- g1[[group_nm]]
  
  res <- sesame::inferInfiniumIChannel(sdf, summary = TRUE)
  for (nm in names(res)) {
    s[[paste0('InfI_switch_', nm)]] <- unname(res[nm])
  }
  new("sesameQC", stat=s, group=g)
}


sesameQC_calcStats_dyeBias <- function(sdf, qc = NULL) {
  
  g1 <- .setGroup_dyeBias()
  group_nm <- names(g1)[1]
  if (is.null(qc)) { s <- list(); g <- list()
  } else { s <- qc@stat; g <- qc@group }
  if (group_nm %in% names(g)) { return(qc); }
  g[[group_nm]] <- g1[[group_nm]]
  
  t1 <- InfI(sdf)
  intens <- sesame::totalIntensities(sdf)
  s$medR <- median(sort(intens[t1[t1$col == "R", "Probe_ID"]]))
  s$medG <- median(sort(intens[t1[t1$col == "G", "Probe_ID"]]))
  s$topR <- median(tail(sort(intens[t1[t1$col == "R", "Probe_ID"]]), n=20))
  s$topG <- median(tail(sort(intens[t1[t1$col == "G", "Probe_ID"]]), n=20))
  s$RGratio <- s$medR / s$medG
  s$RGdistort <- (s$topR / s$topG) / (s$medR / s$medG)
  
  new("sesameQC", stat=s, group=g)
}


sesameQC_calcStats_betas <- function(sdf, qc = NULL) {
  
  g1 <- .setGroup_betas()
  group_nm <- names(g1)[1]
  if (is.null(qc)) { s <- list(); g <- list()
  } else { s <- qc@stat; g <- qc@group }
  if (group_nm %in% names(g)) { return(qc); }
  g[[group_nm]] <- g1[[group_nm]]
  
  betas <- sesame::getBetas(sesame::pOOBAH(sesame::noob(dyeBiasNL(sdf))))
  s$mean_beta <- mean(betas, na.rm = TRUE)
  s$median_beta <- median(betas, na.rm = TRUE)
  s$frac_unmeth <- sum(betas < 0.3, na.rm = TRUE)/sum(!is.na(betas))
  s$frac_meth <- sum(betas > 0.7, na.rm = TRUE)/sum(!is.na(betas))
  s$num_na <- sum(is.na(betas))
  s$frac_na <- sum(is.na(betas)) / length(betas)
  
  for (pt in c('cg','ch','rs')) {
    b1 <- betas[grep(paste0('^', pt), names(betas))]
    s[[paste0('mean_beta_', pt)]] <- mean(b1, na.rm = TRUE)
    s[[paste0('median_beta_', pt)]] <- median(b1, na.rm = TRUE)
    s[[paste0('frac_unmeth_', pt)]] <-
      sum(b1 < 0.3, na.rm = TRUE) / sum(!is.na(b1))
    s[[paste0('frac_meth_', pt)]] <-
      sum(b1 > 0.7, na.rm = TRUE) / sum(!is.na(b1))
    s[[paste0('num_na_', pt)]] <- sum(is.na(b1))
    s[[paste0('frac_na_', pt)]] <- sum(is.na(b1)) / length(b1)
  }
  new("sesameQC", stat=s, group=g)
}


.setGroup_detection <- function() {
  list("Detection" = c(
    num_dtna    = "N. Probes w/ Missing Raw Intensity  ",
    frac_dtna   = "% Probes w/ Missing Raw Intensity   ",
    num_dt      = "N. Probes w/ Detection Success      ",
    frac_dt     = "% Detection Success                 ",
    num_dt_mk   = "N. Detection Succ. (after masking)  ",
    frac_dt_mk  = "% Detection Succ. (after masking)   ",
    num_dt_cg   = "N. Probes w/ Detection Success (cg) ",
    frac_dt_cg  = "% Detection Success (cg)            ",
    num_dt_ch   = "N. Probes w/ Detection Success (ch) ",
    frac_dt_ch  = "% Detection Success (ch)            ",
    num_dt_rs   = "N. Probes w/ Detection Success (rs) ",
    frac_dt_rs  = "% Detection Success (rs)            "))
}




.setGroup_intensity <- function() {
  list("Signal Intensity" = c(
    mean_intensity    = "Mean sig. intensity         ",
    mean_intensity_MU = "Mean sig. intensity (M+U)   ",
    mean_ii           = "Mean sig. intensity (Inf.II)",
    mean_inb_grn      = "Mean sig. intens.(I.Grn IB) ",
    mean_inb_red      = "Mean sig. intens.(I.Red IB) ",
    mean_oob_grn      = "Mean sig. intens.(I.Grn OOB)",
    mean_oob_red      = "Mean sig. intens.(I.Red OOB)",
    na_intensity_M    = "N. NA in M (all probes)     ",
    na_intensity_U    = "N. NA in U (all probes)     ",
    na_intensity_ig   = "N. NA in raw intensity (IG) ",
    na_intensity_ir   = "N. NA in raw intensity (IR) ",
    na_intensity_ii   = "N. NA in raw intensity (II) "))
}




noMasked <- function(sdf) { # filter masked probes
  sdf[!sdf$mask,,drop=FALSE]
}

InfIR <- function(sdf) {
  sdf[sdf$col == "R",,drop=FALSE]
}

InfIG <- function(sdf) {
  sdf[sdf$col == "G",,drop=FALSE]
}

InfI <- function(sdf) {
  sdf[sdf$col != "2",,drop=FALSE]
}

InfII <- function(sdf) {
  sdf[sdf$col == "2",,drop=FALSE]
}

oobG <- function(sdf) {
  dR <- InfIR(sdf)
  c(dR$MG, dR$UG)
}

oobR <- function(sdf) {
  dG <- InfIG(sdf)
  c(dG$MR, dG$UR)
}

.setGroup_numProbes <- function() {
  list("Number of Probes" = c(
    num_probes    = "N. Probes         ",
    num_probes_II = "N. Inf.-II Probes ",
    num_probes_IR = "N. Inf.-I (Red)   ",
    num_probes_IG = "N. Inf.-I (Grn)   ",
    num_probes_cg = "N. Probes (CG)    ",
    num_probes_ch = "N. Probes (CH)    ",
    num_probes_rs = "N. Probes (RS)    "))
}

.setGroup_channel <- function() {
  list("Color Channel" = c(
    InfI_switch_R2R = "N. Inf.I Probes Red -> Red ",
    InfI_switch_G2G = "N. Inf.I Probes Grn -> Grn ",
    InfI_switch_R2G = "N. Inf.I Probes Red -> Grn ",
    InfI_switch_G2R = "N. Inf.I Probes Grn -> Red "))
}

.setGroup_dyeBias <- function() {
  list("Dye Bias" = c(
    medR      = "Median Inf.I Intens. Red           ",
    medG      = "Median Inf.I Intens. Grn           ",
    topR      = "Median of Top 20 Inf.I Intens. Red ",
    topG      = "Median of Top 20 Inf.I Intens. Grn ",
    RGratio   = "Ratio of Red-to-Grn median Intens. ",
    RGdistort = "Ratio of Top vs. Global R/G Ratios "))
}

.setGroup_betas <- function() {
  list("Beta Value" = c(
    mean_beta      = "Mean Beta           ",
    median_beta    = "Median Beta         ",
    frac_unmeth    = "% Beta < 0.3        ",
    frac_meth      = "% Beta > 0.7        ",
    num_na         = "N. is.na(Beta)      ",
    frac_na        = "% is.na(Beta)       ",
    mean_beta_cg   = "Mean Beta (CG)      ",
    median_beta_cg = "Median Beta (CG)    ",
    frac_unmeth_cg = "% Beta < 0.3 (CG)   ",
    frac_meth_cg   = "% Beta > 0.7 (CG)   ",
    num_na_cg      = "N. is.na(Beta) (CG) ",
    frac_na_cg     = "% is.na(Beta) (CG)  ",
    mean_beta_ch   = "Mean Beta (CH)      ",
    median_beta_ch = "Median Beta (CH)    ",
    frac_unmeth_ch = "% Beta < 0.3 (CH)   ",
    frac_meth_ch   = "% Beta > 0.7 (CH)   ",
    num_na_ch      = "N. is.na(Beta) (CH) ",
    frac_na_ch     = "% is.na(Beta) (CH)  ",
    mean_beta_rs   = "Mean Beta (RS)      ",
    median_beta_rs = "Median Beta (RS)    ",
    frac_unmeth_rs = "% Beta < 0.3 (RS)   ",
    frac_meth_rs   = "% Beta > 0.7 (RS)   ",
    num_na_rs      = "N. is.na(Beta) (RS) ",
    frac_na_rs     = "% is.na(Beta) (RS)  "))
}



dyeBiasNL <- function(sdf, mask = TRUE, verbose = FALSE) {
  
  stopifnot(is(sdf, "SigDF"))
  rgdistort <- sesameQC_calcStats(sdf, "dyeBias")@stat$RGdistort
  if (is.na(rgdistort) || rgdistort >10) {
    return(maskIG(sdf)); }
  
  ## we use all Inf-I probes so we capture the entire support range
  if (mask) { dG <- InfIG(sdf); dR <- InfIR(sdf)
  } else { dG <- InfIG(noMasked(sdf)); dR <- InfIR(noMasked(sdf)) }
  IG0 <- c(dG$MG, dG$UG); IR0 <- c(dR$MR, dR$UR)
  
  maxIG <- max(IG0, na.rm = TRUE); minIG <- min(IG0, na.rm = TRUE)
  maxIR <- max(IR0, na.rm = TRUE); minIR <- min(IR0, na.rm = TRUE)
  
  if (maxIG <= 0 || maxIR <= 0) { return(sdf); }
  IR1 <- sort(as.numeric(IR0))
  IR2 <- sort(as.vector(normalize.quantiles.use.target_new(
    matrix(IR1), as.vector(IG0))))
  IRmid <- (IR1 + IR2) / 2.0
  maxIRmid <- max(IRmid); minIRmid <- min(IRmid)
  fitfunRed <- function(data) {
    insupp    <- data <= maxIR & data >= minIR & (!is.na(data))
    oversupp  <- data > maxIR & (!is.na(data))
    undersupp <- data < minIR & (!is.na(data))
    data[insupp] <- approx(x=IR1, y=IRmid, xout=data[insupp], ties=mean)$y
    data[oversupp]  <- data[oversupp] - maxIR + maxIRmid
    data[undersupp] <- minIRmid/ minIR * data[undersupp]
    data
  }
  
  IG1 <- sort(as.numeric(IG0))
  IG2 <- sort(as.vector(normalize.quantiles.use.target_new(
    matrix(IG1), as.vector(IR0))))
  IGmid <- (IG1 + IG2) / 2.0
  maxIGmid <- max(IGmid); minIGmid <- min(IGmid)
  fitfunGrn <- function(data) {
    insupp    <- data <= maxIG & data >= minIG & (!is.na(data))
    oversupp  <- data > maxIG & (!is.na(data))
    undersupp <- data < minIG & (!is.na(data))
    data[insupp] <- approx(x=IG1, y=IGmid, xout=data[insupp], ties=mean)$y
    data[oversupp]  <- data[oversupp] - maxIG + maxIGmid
    data[undersupp] <- minIGmid/ minIG * data[undersupp]
    data
  }
  
  sdf$MR <- fitfunRed(sdf$MR); sdf$UR <- fitfunRed(sdf$UR)
  sdf$MG <- fitfunGrn(sdf$MG); sdf$UG <- fitfunGrn(sdf$UG)
  sdf
}







# normalize.quantiles.use.target <-
#   function (x,
#             target,
#             copy = TRUE,
#             subset = NULL)
#   {
#     if (!is.matrix(x)) {
#       stop("This function expects supplied argument to be matrix")
#     }
#     if (!is.numeric(x)) {
#       stop("Supplied argument should be a numeric matrix")
#     }
#     rows <- dim(x)[1]
#     cols <- dim(x)[2]
#     if (is.integer(x)) {
#       x <- matrix(as.double(x), rows, cols)
#     }
#     if (!is.vector(target)) {
#       stop("This function expects target to be vector")
#     }
#     if (!is.numeric(target)) {
#       stop("Supplied target argument should be a numeric vector")
#     }
#     if (is.integer(target)) {
#       target <- as.double(target)
#     }
#     if (is.null(subset)) {
#       return(.Call("R_qnorm_using_target", x, target, copy,
#                    PACKAGE = "preprocessCore"))
#     } else {
#       if (length(subset) != rows) {
#         stop("subset should have same length as nrows(x)")
#       }
#       subset <- as.integer(subset)
#       return(
#         .Call(
#           "R_qnorm_using_target_via_subset",
#           x,
#           subset,
#           target,
#           copy,
#           PACKAGE = "preprocessCore"
#         )
#       )
#     }
#   }


normalize.quantiles.use.target_new <- function(x, target, copy = TRUE, subset = NULL) {
  # 输入验证
  if (!is.matrix(x)) {
    stop("This function expects supplied argument to be a matrix")
  }
  if (!is.numeric(x)) {
    stop("Supplied argument should be a numeric matrix")
  }
  if (!is.vector(target)) {
    stop("This function expects target to be a vector")
  }
  if (!is.numeric(target)) {
    stop("Supplied target argument should be a numeric vector")
  }
  if (!is.null(subset) && length(subset) != nrow(x)) {
    stop("subset should have same length as nrows(x)")
  }
  
  # 如果需要拷贝，创建数据的副本
  if (copy) {
    result <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
  } else {
    result <- x
  }
  
  # 将target转换为数值向量
  target <- as.numeric(target)
  
  # 对target进行排序
  sorted_target <- sort(target)
  
  if (is.null(subset)) {
    # 对整个矩阵进行归一化
    for (j in 1:ncol(result)) {
      # 获取当前列的非NA值
      col_data <- result[, j]
      na_mask <- is.na(col_data)
      valid_data <- col_data[!na_mask]
      
      if (length(valid_data) > 0) {
        # 对当前列的值进行排序并获取排序索引
        sorted_indices <- order(valid_data)
        
        # 创建与目标分布匹配的归一化值
        normalized_values <- rep(NA, length(valid_data))
        
        # 将排序后的值映射到目标分布
        # 使用线性插值来处理长度不匹配的情况
        if (length(valid_data) == length(sorted_target)) {
          normalized_values[sorted_indices] <- sorted_target
        } else {
          # 创建插值函数
          source_indices <- seq_along(valid_data)
          target_indices <- seq(1, length(sorted_target), 
                                length.out = length(valid_data))
          
          # 对目标值进行插值
          interpolated_target <- approx(seq_along(sorted_target), sorted_target, 
                                        target_indices, method = "linear")$y
          
          normalized_values[sorted_indices] <- interpolated_target
        }
        
        # 将归一化后的值放回原位置
        result[!na_mask, j] <- normalized_values
      }
    }
  } else {
    # 只对子集进行归一化
    subset <- as.logical(subset)
    if (sum(subset) == 0) {
      warning("No elements selected by subset")
      return(result)
    }
    
    # 提取子集数据
    subset_data <- result[subset, , drop = FALSE]
    
    # 对子集数据进行归一化
    for (j in 1:ncol(subset_data)) {
      # 获取当前列的非NA值
      col_data <- subset_data[, j]
      na_mask <- is.na(col_data)
      valid_data <- col_data[!na_mask]
      
      if (length(valid_data) > 0) {
        # 对当前列的值进行排序并获取排序索引
        sorted_indices <- order(valid_data)
        
        # 创建与目标分布匹配的归一化值
        normalized_values <- rep(NA, length(valid_data))
        
        # 将排序后的值映射到目标分布
        if (length(valid_data) == length(sorted_target)) {
          normalized_values[sorted_indices] <- sorted_target
        } else {
          # 创建插值函数
          source_indices <- seq_along(valid_data)
          target_indices <- seq(1, length(sorted_target), 
                                length.out = length(valid_data))
          
          # 对目标值进行插值
          interpolated_target <- approx(seq_along(sorted_target), sorted_target, 
                                        target_indices, method = "linear")$y
          
          normalized_values[sorted_indices] <- interpolated_target
        }
        
        # 将归一化后的值放回子集
        subset_data[!na_mask, j] <- normalized_values
      }
    }
    
    # 将处理后的子集数据放回原矩阵
    result[subset, ] <- subset_data
  }
  
  return(result)
}





dyeBiasCorrTypeINorm <- function (sdf, mask = TRUE, verbose = FALSE) 
{
  stopifnot(is(sdf, "SigDF"))
  rgdistort <- sesameQC_calcStats(sdf, "dyeBias")@stat$RGdistort
  if (is.na(rgdistort) || rgdistort > 10) {
    return(maskIG(sdf))
  }
  if (mask) {
    dG <- InfIG(sdf)
    dR <- InfIR(sdf)
  }
  else {
    dG <- InfIG(noMasked(sdf))
    dR <- InfIR(noMasked(sdf))
  }
  IG0 <- c(dG$MG, dG$UG)
  IR0 <- c(dR$MR, dR$UR)
  maxIG <- max(IG0, na.rm = TRUE)
  minIG <- min(IG0, na.rm = TRUE)
  maxIR <- max(IR0, na.rm = TRUE)
  minIR <- min(IR0, na.rm = TRUE)
  if (maxIG <= 0 || maxIR <= 0) {
    return(sdf)
  }
  IR1 <- sort(as.numeric(IR0))
  IR2 <- sort(as.vector(normalize.quantiles.use.target_new(matrix(IR1), 
                                                           as.vector(IG0))))
  IRmid <- (IR1 + IR2)/2
  maxIRmid <- max(IRmid)
  minIRmid <- min(IRmid)
  fitfunRed <- function(data) {
    insupp <- data <= maxIR & data >= minIR & (!is.na(data))
    oversupp <- data > maxIR & (!is.na(data))
    undersupp <- data < minIR & (!is.na(data))
    data[insupp] <- approx(x = IR1, y = IRmid, xout = data[insupp], 
                           ties = mean)$y
    data[oversupp] <- data[oversupp] - maxIR + maxIRmid
    data[undersupp] <- minIRmid/minIR * data[undersupp]
    data
  }
  IG1 <- sort(as.numeric(IG0))
  IG2 <- sort(as.vector(normalize.quantiles.use.target_new(matrix(IG1), 
                                                           as.vector(IR0))))
  IGmid <- (IG1 + IG2)/2
  maxIGmid <- max(IGmid)
  minIGmid <- min(IGmid)
  fitfunGrn <- function(data) {
    insupp <- data <= maxIG & data >= minIG & (!is.na(data))
    oversupp <- data > maxIG & (!is.na(data))
    undersupp <- data < minIG & (!is.na(data))
    data[insupp] <- approx(x = IG1, y = IGmid, xout = data[insupp], 
                           ties = mean)$y
    data[oversupp] <- data[oversupp] - maxIG + maxIGmid
    data[undersupp] <- minIGmid/minIG * data[undersupp]
    data
  }
  sdf$MR <- fitfunRed(sdf$MR)
  sdf$UR <- fitfunRed(sdf$UR)
  sdf$MG <- fitfunGrn(sdf$MG)
  sdf$UG <- fitfunGrn(sdf$UG)
  sdf
}








