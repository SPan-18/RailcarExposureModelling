# print predictive information criteria

LL <- post.loglik
Dev <- -2*colSums(LL)
dic <- list(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)
waic <- WAIC(LL)

cat("\n \t PREDICTIVE INFORMATION CRITERIA \n")
cat("Deviance Information Criteria (Dbar = mean(deviance), pV = var(deviance)/2) \n")
cat("DIC = ", dic$DIC, "\t Dbar = ", dic$Dbar, "\t pV = ", dic$pV, "\n")
cat("Widely Applicable Information Criteria (Watanabe, 2010) \n")
cat("WAIC = ", waic$WAIC, "\t lppd = ", waic$lppd, "\t pWAIC = ", 
    waic$pWAIC, "or", waic$pWAIC1, "\n")
cat("(lppd = log pointwise predictive density,\n pWAIC = effective number of parameters)\n")
