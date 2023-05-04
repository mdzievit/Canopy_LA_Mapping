qtl.analysis <- function (crossobj = crossobj, trait = "pred", step, method, 
                          threshold, distance, cofactors, window.size = 50, print = FALSE,
                          write = FALSE) 
{
  dir.create("qtl_analysis_reports", showWarnings = F)
  QTL.result <- NULL
  trait <- tolower(trait)
  crossobj <- calc.genoprob(crossobj, step = step)
  mark.names <- c()
  for (i in 1:length(crossobj$geno)) {
    mark.names <- c(mark.names, colnames(crossobj$geno[[i]]$data))
  }
  probs1 <- NULL
  probs2 <- NULL
  for (i in 1:nchr(crossobj)) {
    names <- dimnames(crossobj$geno[[i]]$prob)[[2]]
    sel.mark <- which(names %in% mark.names == FALSE)
    names[sel.mark] <- paste(i, names[sel.mark], sep = "_")
    dimnames(crossobj$geno[[i]]$prob)[[2]] <- names
    names(attributes(crossobj$geno[[i]]$prob)$map) <- names
    p1 <- crossobj$geno[[i]]$prob
    names(p1) <- dimnames(p1)[[2]]
    probs1 <- cbind(probs1, p1[, , 1])
    if (class(crossobj)[1] == "f2" | class(crossobj)[1] == 
        "4way") {
      probs2 <- cbind(probs2, p1[, , 3])
    }
    if (class(crossobj)[1] == "dh" | class(crossobj)[1] == 
        "bc" | class(crossobj)[1] == "riself" | class(crossobj)[1] == 
        "ri4self" | class(crossobj)[1] == "ri8self" | class(crossobj)[1] == 
        "risib" | class(crossobj)[1] == "ri4sib" | class(crossobj)[1] == 
        "ri8sib") {
      probs2 <- cbind(probs2, p1[, , 2])
    }
  }
  additive <- probs2 - probs1
  new.additive <- additive
  P.data <- crossobj$pheno
  Gen <- "id"
  Trait <- paste(trait)
  ENV <- as.factor(rep("env", nrow(P.data)))
  GEN <- as.factor(P.data[, Gen])
  MEAN <- as.numeric(as.matrix(P.data[, Trait]))
  all.means <- data.frame(ENV, GEN, MEAN, stringsAsFactors = FALSE)
  b <- matrix(, 0, 2)
  for (i in 1:nchr(crossobj)) {
    a <- paste("crossobj$geno$'", i, "'$prob", sep = "")
    mp <- attributes(eval(parse(text = a)))$map
    mp <- cbind(rep(i, length(mp)), as.matrix(mp))
    b <- rbind(b, mp)
  }
  if (method == "SIM") {
    p.values <- NULL
    fixeff <- NULL
    for (i in 1:dim(new.additive)[2]) {
      marker <- new.additive[, i]
      GE.data <- data.frame(all.means, marker)
      CS.Model <- lm(MEAN ~ marker, data = GE.data)
      fstats <- as.vector(summary(CS.Model)$fstatistic)
      p.value <- 1 - pf(fstats[1], fstats[2], fstats[3])
      p.value <- data.frame(rownames(b)[i], b[i, 1], b[i, 
                                                       2], p.value, stringsAsFactors = FALSE)
      p.values <- rbind(p.values, p.value)
      f <- CS.Model$coefficients[2]
    }
    fixeff <- cbind(fixeff, f)
  }
  if (method == "CIM") {
    cofactor.list <- NULL
    cofactor.pos <- NULL
    cofactor.list <- as.matrix(new.additive[, match(cofactors, 
                                                    colnames(new.additive))])
    cofactor.pos <- as.matrix(b[match(cofactors, row.names(b)), 
                                ])
    if (ncol(cofactor.pos) == 1) {
      cofactor.pos <- t(cofactor.pos)
    }
    cofactor.win.f <- rep(0, dim(b)[1])
    for (j in 1:length(cofactors)) {
      c <- b[, 1] == as.numeric(as.character(cofactor.pos[j, 
                                                          1]))
      c[c == FALSE] <- 0
      c[c == TRUE] <- 1
      win <- (b[, 2] > (as.numeric(as.character(cofactor.pos[j, 
                                                             2]))[1] - 0.5 * window.size) & b[, 1] < (as.numeric(as.character(cofactor.pos[j, 
                                                                                                                                           2]))[1] + 0.5 * window.size))
      win[win == FALSE] <- 0
      win[win == TRUE] <- 1
      cofactor.win <- c * win
      cofactor.win[cofactor.win == 1] <- j
      cofactor.win.f <- cofactor.win.f + cofactor.win
    }
    p.values <- NULL
    fixeff <- NULL
    for (i in 1:dim(new.additive)[2]) {
      marker <- new.additive[, i]
      if (unlist(cofactor.win.f[i]) > 0) {
        new.list <- cofactor.list[, -c(unlist(cofactor.win.f[i]))]
        number.cofactors <- length(cofactors) - 1
      }
      if (unlist(cofactor.win.f[i]) == 0) {
        new.list <- cofactor.list
        number.cofactors <- length(cofactors)
      }
      GE.data <- data.frame(all.means, marker)
      if (length(levels(ENV)) == 1 & number.cofactors > 
          0) {
        CS.Model <- lm(MEAN ~ new.list + marker, data = GE.data)
        fstats <- c(anova(CS.Model)[2, 4], anova(CS.Model)[2, 
                                                           1], anova(CS.Model)[3, 1])
      }
      if (length(levels(ENV)) == 1 & number.cofactors == 
          0) {
        CS.Model <- lm(MEAN ~ +marker, data = GE.data)
        fstats <- c(anova(CS.Model)[1, 4], anova(CS.Model)[1, 
                                                           1], anova(CS.Model)[2, 1])
      }
      p.value <- 1 - pf(fstats[1], fstats[2], fstats[3])
      p.value <- data.frame(rownames(b)[i], b[i, 1], b[i, 
                                                       2], p.value, stringsAsFactors = FALSE)
      p.values <- rbind(p.values, p.value)
      if (length(levels(ENV)) == 1) {
        f <- CS.Model$coefficients[ncol(new.list) + 
                                     2]
      }
      fixeff <- cbind(fixeff, f)
    }
  }
  b
  names(p.values) <- c("marker", "Chr", "Pos", "-log10(p)")
  outem <- p.values
  if (threshold > 0) {
    threshold.f <- threshold
  }
  if (threshold == "Li&Ji") {
    scores <- NULL
    for (i in 1:nchr(crossobj)) {
      a <- paste("crossobj$geno$'", i, "'$data", sep = "")
      s <- eval(parse(text = a))
      scores <- cbind(scores, s)
    }
    scores[scores == 1] <- 0
    average <- NULL
    for (i in 1:dim(scores)[2]) {
      a <- mean(scores[, i], na.rm = TRUE)
      average <- cbind(average, a)
    }
    for (i in 1:dim(scores)[1]) {
      for (j in 1:dim(scores)[2]) {
        if (is.na(scores[i, j]) == TRUE) {
          scores[i, j] <- average[1, j]
        }
      }
    }
    pca.analysis1 <- prcomp(scores, scale = TRUE)
    meff1 <- nind(crossobj) - 1
    ngeno <- nind(crossobj)
    lambda <- summary(pca.analysis1)$importance[2, ]
    lambda <- lambda[1:meff1]
    TM <- NULL
    x <- 10
    while (x > 0.9792895) {
      meff <- length(lambda)
      sum.lambda <- sum(lambda, na.rm = TRUE)
      sum.lambda2 <- sum(lambda^2, na.rm = TRUE)
      neff <- ((ngeno + 1) * (sum.lambda^2))/(((ngeno - 
                                                  1) * sum.lambda2) - (sum.lambda^2))
      mu <- ((sqrt(neff - 1) + sqrt(meff))^2)/(neff)
      sigma <- ((sqrt(neff - 1) + sqrt(meff))/neff) * 
        (((1/(sqrt(neff - 1))) + (1/sqrt(meff)))^(1/3))
      L1 <- (meff * lambda[1])/sum.lambda
      x <- (L1 - mu)/sigma
      TM <- c(TM, x)
      lambda <- lambda[2:length(lambda)]
    }
    n.signif <- length(TM) - 1
    alpha.p <- 1 - ((1 - 0.05)^(1/n.signif))
    threshold.f <- (-log10(alpha.p))
  }
  pot.qtl <- outem[which(outem$"-log10(p)" < threshold.f), 
                   ]
  res.qtl <- c()
  if (nrow(pot.qtl) > 0) {
    pot.qtl$select <- 1
    pot.qtl$eval <- 0
    for (chr in unique(pot.qtl$Chr)) {
      t.pot.qtl <- pot.qtl[pot.qtl$Chr == chr, ]
      while (sum(t.pot.qtl$eval) < nrow(t.pot.qtl)) {
        min.p <- min(t.pot.qtl$"-log10(p)"[which(t.pot.qtl$eval == 
                                                   0)])
        sel.row <- which(t.pot.qtl$"-log10(p)" == min.p & 
                           t.pot.qtl$eval == 0)[1]
        d <- abs(t.pot.qtl$Pos - t.pot.qtl$Pos[sel.row])
        t.pot.qtl$select[d <= distance] <- 0
        t.pot.qtl$eval[d <= distance] <- 1
        t.pot.qtl$select[sel.row] <- 1
      }
      t.pot.qtl <- t.pot.qtl[which(t.pot.qtl$select == 
                                     1), ]
      res.qtl <- rbind(res.qtl, t.pot.qtl)
    }
    res.qtl$select <- NULL
    res.qtl$eval <- NULL
  }
  if (max(-log10(outem[, 4]), na.rm = TRUE) == "Inf") {
    max <- 10
  }
  if (max(-log10(outem[, 4]), na.rm = TRUE) != "Inf") {
    max <- (max(-log10(outem[, 4])) + 0.05)
  }
  bw.sel <- function(y, X, alpha = 0) {
    XXX <- as.matrix(X)
    pval <- rep(0, ncol(XXX))
    names(pval) <- colnames(XXX)
    R.vec <- c()
    test <- 1
    while (test == 1) {
      toss <- which(pval == max(pval) & pval > alpha)[1]
      sel <- setdiff(c(1:ncol(XXX)), toss)
      if (length(sel) == 0) {
        break
      }
      nms <- colnames(XXX)[sel]
      XXX <- as.matrix(XXX[, sel])
      colnames(XXX) <- nms
      CS.Model <- lm(y ~ XXX)
      pval <- summary(CS.Model)$coefficients[, 4]
      pval <- pval[2:length(pval)]
      test <- ((sum(pval >= alpha) > 0) * 1)
    }
    X <- XXX
    CS.Model <- lm(y ~ X)
    coef <- coefficients(CS.Model)
    coef <- as.vector(coef[2:length(coef)])
    qtl.names <- names(coef)
    qtl.names <- colnames(X)
    Rsq <- summary(CS.Model)$r.squared
    for (i in 1:ncol(X)) {
      toss <- i
      if (ncol(X) > 1) {
        sel <- setdiff(c(1:ncol(X)), toss)
      }
      else {
        sel <- toss
      }
      XXX <- as.matrix(X[, sel])
      CS.Model <- lm(y ~ XXX)
      r <- summary(CS.Model)$r.squared
      R.vec <- c(R.vec, Rsq - r)
    }
    if (ncol(X) == 1) {
      qtl.names <- colnames(X)
      R.vec <- r
    }
    return(list(qtl.names, R.vec, coef))
  }
  if (nrow(pot.qtl) > 0) {
    new.additive.final <- as.matrix(new.additive[, colnames(new.additive) %in% 
                                                   res.qtl$marker])
    colnames(new.additive.final) <- colnames(new.additive)[colnames(new.additive) %in% 
                                                             res.qtl$marker]
    BW <- bw.sel(MEAN, new.additive.final, alpha = 0.05)
    qtl.names <- BW[[1]]
    res.qtl <- res.qtl[match(qtl.names, res.qtl$marker), 
                       ]
    m.eff <- BW[[3]]
    Rsq <- BW[[2]]
  }
  if (nrow(pot.qtl) == 0) {
    m.eff <- NA
    Rsq <- NA
  }
  QTL.result$all <- outem
  if (method == "CIM") {
    QTL.result$selected <- cbind(res.qtl, m.eff, Rsq)
  }
  if (method == "SIM") {
    QTL.result$selected <- res.qtl
  }
  m.eff <- NULL
  QTL.result$all$"-log10(p)" <- (-log10(QTL.result$all$"-log10(p)"))
  if (nrow(pot.qtl) > 0) {
    QTL.result$selected$"-log10(p)" <- (-log10(QTL.result$selected$"-log10(p)"))
  }
  QTL.result$interaction <- NULL
  if (step == 0) {
    method <- "SMA"
  }
  if(write) {
  filename.1 <- paste("qtl_analysis_reports/QTL_summary_", 
                      trait, "_", method, "_", ".txt", sep = "")
  write.table(QTL.result$all, file = filename.1, sep = ",", 
              eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)
  filename.2 <- paste("qtl_analysis_reports/QTL_selected_", 
                      trait, "_", method, "_", ".txt", sep = "")
  write.table(QTL.result$selected, file = filename.2, sep = ",", 
              eol = "\n", na = "-", dec = ".", col.names = T, row.names = F)
  }
  if(print){
    print(QTL.result$all)
    print(QTL.result$selected)
    print(xyplot(-log10(outem[, 4]) ~ outem[, 3] | factor(outem[, 
                                                                2]), type = "l", layout = c(nchr(crossobj), 1), col = "red", 
                 xlab = "Chromosome position", ylab = "-log10(P)", main = paste("QTL mapping", 
                                                                                method, sep = ""), scales = list(x = "free"), lwd = 3, 
                 panel = function(x, y, ...) {
                   panel.abline(h = -log10(threshold.f), lty = 2)
                   llines(x, y, col = "red", lwd = 2)
                 }))
  }
  QTL.result$threshold <- threshold.f
  QTL.result
}
