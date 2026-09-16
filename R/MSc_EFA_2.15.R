# Link biomechanical measures to EFA shape positions and save the plots.
allo_targets <- c("shapePC1","shapePC2","wing_curvature")
allo_preds   <- BIO_VARS[BIO_VARS %in% names(sp)]

# H1. Correlations: biomechanical metrics ↔ shape PC positions
# H1. Correlations: biomechanical metrics ↔ shape PC positions
# Check which measurements rise or fall with each shape axis.
bio_mat   <- sp[, c(allo_preds,"shapePC1","shapePC2")]
corr_mat  <- cor(bio_mat, use="pairwise.complete.obs")
save_csv(as.data.frame(corr_mat), "EFA_biomech_shapePC_correlations")
cat("  Correlation (shapePC1):\n")
print(round(corr_mat["shapePC1", allo_preds], 3))
cat("  Correlation (shapePC2):\n")
print(round(corr_mat["shapePC2", allo_preds], 3))

# H2. Linear models: each biomechanical predictor → shapePC1 and shapePC2
# H2. Linear models: each biomechanical predictor → shapePC1 and shapePC2
# Test one predictor at a time against each shape axis.
lm_results <- do.call(rbind, lapply(allo_preds, function(pred) {
  do.call(rbind, lapply(c("shapePC1","shapePC2"), function(resp) {
    d <- sp[,c(resp,pred)]; d <- d[complete.cases(d),]
    if (nrow(d) < 5) return(NULL)
    m <- lm(as.formula(paste(resp,"~",pred)), data=d)
    s <- summary(m)
    data.frame(Response=resp, Predictor=pred,
               R2=round(s$r.squared,4),
               Adj_R2=round(s$adj.r.squared,4),
               F_stat=round(s$fstatistic[1],3),
               p_value=round(pf(s$fstatistic[1],s$fstatistic[2],
                                s$fstatistic[3],lower.tail=FALSE),6))
  }))
}))
save_csv(lm_results, "EFA_lm_biomech_predictors_shapePC")

# H3. Multiple regression: all biomechanical metrics → shapePC1 and shapePC2
# H3. Multiple regression: all biomechanical metrics → shapePC1 and shapePC2
# Test all measurements together for each shape axis.
for (resp in c("shapePC1","shapePC2")) {
  d_multi <- sp[, c(resp, allo_preds)] %>% drop_na()
  if (nrow(d_multi) < 10) next
  m_multi <- lm(as.formula(paste(resp, "~", paste(allo_preds, collapse="+"))), data=d_multi)
  s <- summary(m_multi)
  coef_df <- as.data.frame(coef(s))
  coef_df$term <- rownames(coef_df)
  coef_df$Response <- resp
  save_csv(coef_df, paste0("EFA_multiple_regression_", resp))
  cat(sprintf("  Multiple regression %s: R2=%.3f, adj-R2=%.3f, F-p=%.4f\n",
              resp, s$r.squared, s$adj.r.squared,
              pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)))
}

# H4. RandomForest: predict shapePC1 and shapePC2 from biomechanical metrics
# H4. RandomForest: predict shapePC1 and shapePC2 from biomechanical metrics
# Use a tree-based model to see which measurements matter most.
# (clade-level analysis)
rf_results <- list()
for (resp in c("shapePC1","shapePC2")) {
  d_rf <- sp[, c(resp, allo_preds)] %>% drop_na()
  if (nrow(d_rf) < 20) next
  rf <- tryCatch(
    randomForest(as.formula(paste(resp,"~",paste(allo_preds,collapse="+"))),
                 data=d_rf, ntree=500, importance=TRUE),
    error=function(e) NULL
  )
  if (!is.null(rf)) {
    rf_results[[resp]] <- rf
    imp_df <- as.data.frame(importance(rf))
    imp_df$Variable <- rownames(imp_df)
    imp_df$Response <- resp
    save_csv(imp_df, paste0("EFA_rf_importance_", resp))

    # Variable importance plot
    p_imp <- ggplot(imp_df %>% slice_max(`%IncMSE`, n=10),
                    aes(reorder(Variable,`%IncMSE`), `%IncMSE`)) +
      geom_col(fill="#2980B9", alpha=0.85) +
      coord_flip() +
      labs(title=paste("RF variable importance for", resp),
           subtitle="Top 10 biomechanical predictors of EFA shape position",
           x="Biomechanical variable", y="% Increase in MSE") +
      theme_bw(base_size=12)
    save_plot(p_imp, paste0("EFA_rf_importance_", resp), w=9, h=6)

    # Predicted vs actual
    d_rf$predicted <- predict(rf)
    r2 <- round(cor(d_rf[[resp]], d_rf$predicted, use="complete.obs")^2, 3)
    p_pva <- ggplot(d_rf, aes(.data[[resp]], predicted)) +
      geom_point(alpha=0.6, colour="#2980B9", size=2) +
      geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey40") +
      geom_smooth(method="lm", se=FALSE, colour="black", linetype="dotted") +
      labs(title=paste("RF predicted vs actual:", resp),
           subtitle=paste0("R² = ", r2),
           x=paste("Actual", resp), y=paste("Predicted", resp)) +
      theme_bw(base_size=12)
    save_plot(p_pva, paste0("EFA_rf_predicted_vs_actual_", resp), w=8, h=7)
  }
}

# H5. Wing curvature gap by clade (Msc_allo_3 adaptation)
# H5. Wing curvature gap by clade (Msc_allo_3 adaptation)
# Compare observed curvature with the value predicted from the other measures.
if ("wing_curvature" %in% names(sp) && !is.null(rf_results$shapePC2)) {
  d_curv <- sp[, c("shapePC2","wing_curvature","clade",allo_preds)] %>% drop_na()
  rf_c <- tryCatch(
    randomForest(wing_curvature ~ ., data=d_curv %>% select(-clade,-shapePC2),
                 ntree=500, importance=TRUE),
    error=function(e) NULL
  )
  if (!is.null(rf_c)) {
    d_curv$curv_pred <- predict(rf_c, newdata=d_curv)
    d_curv$curv_gap  <- d_curv$curv_pred - d_curv$wing_curvature
    p_gap <- ggplot(d_curv, aes(fct_infreq(as.factor(clade)), curv_gap)) +
      geom_boxplot(aes(fill=clade), alpha=0.55, outlier.shape=NA) +
      geom_jitter(size=1.2, alpha=0.45, width=0.2) +
      geom_hline(yintercept=0, linetype="dashed", colour="red") +
      labs(title="Curvature gap by clade (EFA space)",
           subtitle="Positive = skeleton predicts more curvature than observed",
           x="Clade", y="Curvature gap (predicted − actual)") +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none")
    save_plot(p_gap, "EFA_curvature_gap_by_clade", w=12, h=7)
    gap_summary <- d_curv %>% group_by(clade) %>%
      summarise(n=n(), mean_gap=round(mean(curv_gap,na.rm=TRUE),5),
                sd_gap=round(sd(curv_gap,na.rm=TRUE),5), .groups="drop")
    save_csv(gap_summary, "EFA_curvature_gap_summary")
  }
}

# H6. Scatter plots: biomechanical metrics vs shapePC1/shapePC2 coloured by clade
# H6. Scatter plots: biomechanical metrics vs shapePC1/shapePC2 coloured by clade
# Save simple plots so each relationship can be inspected by eye.
for (pred in allo_preds) {
  for (resp in c("shapePC1","shapePC2")) {
    d_plot <- sp %>% filter(!is.na(.data[[pred]]), !is.na(.data[[resp]]))
    if (nrow(d_plot) < 5) next
    r_val <- round(cor(d_plot[[pred]], d_plot[[resp]], use="complete.obs"), 3)
    p_sc <- ggplot(d_plot, aes(.data[[pred]], .data[[resp]], colour=clade)) +
      geom_point(size=2, alpha=0.75) +
      geom_smooth(method="lm", se=TRUE, colour="black", linewidth=0.8, alpha=0.15) +
      labs(title=paste(pred, "vs", resp),
           subtitle=paste0("Pearson r = ", r_val),
           x=pred, y=if(resp=="shapePC1") efa_lab1 else efa_lab2,
           colour="Clade") +
      theme_bw(base_size=11)
    save_plot(p_sc, paste0("EFA_scatter_", pred, "_vs_", resp), w=9, h=6)
  }
}
