library(RegionScan)
data("REGIONinfo","SNPinfo","phenocov","geno")

results<-regscan(phenocov = phenocov, pheno="sim_QT", REGIONinfo=REGIONinfo,
                 geno_type="D", pheno_type="C", data = geno, SNPinfo = SNPinfo )

# locus plot 
LocusPlot(chr=16,pheno="sim_QT",regscanout=results,regionlist=c(1:15),outname="output/RegionScan/LocusPlot_region1_15",
          region_tests=c("Wald.p","PC80.p","MLCB.p","SKATO.p"))

regscan <- function (phenocov = NULL, pheno, REGIONinfo, geno_type, pheno_type, 
          covlist = NULL, data = NULL, SNPinfo = NULL, vcfname = NULL, 
          machr2 = NULL, qcinput = NULL, info_score = NULL, multiallelic = FALSE, 
          multial_nmaxalleles = 2, mafcut = 0.05, rcut = 0.99, firthreg = FALSE, 
          MLCheatmap = FALSE, regionlist = NULL, alltests = FALSE, 
          edgecut = 0.5, tol = 1e-16, qcmachr2 = NULL, covout = FALSE, 
          LDpruning = TRUE, singleSNPall = FALSE, SKAT_kernel = "linear.weighted", 
          SKAT_weights = NULL, SKAT_weights_beta = c(1, 25), verbose = FALSE, 
          debug = FALSE, parallel = FALSE, path=NULL) 
{
  if (is.null(pheno_type)) {
    return("Please specify if the phenotype is continuous (pheno_type='C') or dichotomous (pheno_type='D')")
  }
  if (is.null(geno_type)) {
    return("Please specify if the genotypes are in allele dosage format (geno_type='D') or genotype format (geno_type='G')")
  }
  if (pheno_type == "C") 
    family = "gaussian"
  if (pheno_type == "D") 
    family = "binomial"
  options(warn = -1)
  message("*********************************")
  message("Input parameters")
  message(paste("pheno:", pheno, sep = ""))
  message(toString(paste("covlist:", toString(covlist))))
  message(paste("mafcut:", mafcut, sep = ""))
  message(paste("edgecut:", edgecut, sep = ""))
  message(paste("rcut:", rcut, sep = ""))
  message(paste("geno_type:", geno_type, sep = ""))
  message(paste("pheno_type:", pheno_type, sep = ""))
  message(paste("multiallelic:", multiallelic, sep = ""))
  message(paste("firthreg:", firthreg, sep = ""))
  message(paste("covout:", covout, sep = ""))
  message(paste("LDpruning:", LDpruning, sep = ""))
  message(paste("alltests:", alltests, sep = ""))
  message(paste("tol:", tol, sep = ""))
  message(paste("singleSNPall:", singleSNPall, sep = ""))
  message(paste("MLCheatmap:", MLCheatmap, sep = ""))
  message(paste("debug:", debug, sep = ""))
  message("*********************************")
  if (is.null(vcfname) && is.null(data)) {
    message(cat("ERROR: vcfname or data is required"))
    break
  }
  if (is.null(SNPinfo)) {
    message(cat("ERROR: SNPinfo is required"))
    break
  }
  if (is.null(REGIONinfo)) {
    message(cat("ERROR: REGIONinfo is required"))
    break
  }
  if (is.null(phenocov)) {
    message(cat("ERROR: phenocov is required"))
    break
  }
  if (is.null(pheno)) {
    message(cat("ERROR: pheno is required"))
    break
  }
  if (is.null(pheno_type)) {
    message(cat("ERROR: pheno_type is required"))
    break
  }
  if (is.null(geno_type)) {
    message(cat("ERROR: geno_type is required"))
    break
  }
  if (!is.null(regionlist)) {
    REGIONinfo <- subset(REGIONinfo, region %in% regionlist)
  }
  inputs <- lapply(REGIONinfo$region, function(regionc) {
    tryCatch({
      regioncur <- REGIONinfo[which(REGIONinfo$region == 
                                      regionc), ]
      region <- as.character(regioncur["region"])
      chr <- as.numeric(regioncur["chr"])
      start <- as.numeric(regioncur["start.bp"])
      end <- as.numeric(regioncur["end.bp"])
      SNPinfo_region <- subset(SNPinfo, SNPinfo$bp >= 
                                 start & SNPinfo$bp <= end)
      datapcov <- phenocov[, c(pheno, covlist), drop = FALSE]
      data_region <- na.omit(cbind(datapcov, data[, SNPinfo_region$variant, 
                                                  drop = F]))
      return(regionc = list(SNPinfo_region = SNPinfo_region, 
                            data_region = data_region, regioninfo = regioncur))
    }, error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    })
  })
  rm(data)
  rm(REGIONinfo)
  if (isTRUE(parallel)) {
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    foreach::getDoParRegistered()
    foreach::getDoParWorkers()
  }
  rscan <- foreach::foreach(regionc = inputs, .combine = "rbind", 
                            .multicombine = T, .verbose = verbose, .packages = c("igraph", 
                                                                                 "rms", "SKAT", "foreach", "COMBAT", "brglm2", "ggplot2", 
                                                                                 "reshape2", "doSNOW", "karyoploteR")) %dopar% {
                                                                                   regionout <- binout <- snpout <- outremoved <- filterout <- outsingleSNPall <- prunout <- assignout <- covoutput <- covarout <- NULL
                                                                                   regioninfo <- regionc$regioninfo
                                                                                   message(paste("Region", regioninfo$region, sep = ":"))
                                                                                   region <- as.character(regioninfo["region"])
                                                                                   chr <- as.numeric(regioninfo["chr"])
                                                                                   start <- as.numeric(regioninfo["start.bp"])
                                                                                   end <- as.numeric(regioninfo["end.bp"])
                                                                                   SNPinfo_region <- regionc[["SNPinfo_region"]]
                                                                                   data_region <- regionc[["data_region"]]
                                                                                   tryCatch({
                                                                                     if (!is.null(SNPinfo_region)) {
                                                                                       processout <- RegionScan:::processing(data = data_region, 
                                                                                                                             pheno = pheno, covlist = covlist, SNPinfo = SNPinfo_region, 
                                                                                                                             mafcut = mafcut, multiallelic = multiallelic)
                                                                                       nSNPs <- nSNPs.kept <- length(processout$SNPinfo$variant)
                                                                                       if (isTRUE(debug)) {
                                                                                         save(processout, pheno, covlist, file = paste("chr", 
                                                                                                                                       chr, "_", region, "_checking.Rdata", sep = ""))
                                                                                       }
                                                                                       if (isTRUE(nSNPs > 1)) {
                                                                                         clustout <- RegionScan:::clustering(data = processout$data[, 
                                                                                                                                                    processout$bialSNP], edgecut = edgecut)
                                                                                         if (length(processout$multialSNP.kept) > 0) {
                                                                                           codechange.bial <- RegionScan:::recoding(binlist = clustout$binlist, 
                                                                                                                                    data = processout$data)
                                                                                           assignout <- RegionScan:::assigning(binvector = clustout$binvector, 
                                                                                                                               codechange = codechange.bial$newcode, 
                                                                                                                               multialSNP.kept = processout$multialSNP.kept, 
                                                                                                                               data = processout$data, edgecut = edgecut)
                                                                                           if (isTRUE(LDpruning)) {
                                                                                             prunout <- RegionScan:::pruning(binvector = assignout, 
                                                                                                                             data = processout$data, rcut = rcut, 
                                                                                                                             multialSNP.setA = processout$multialSNP.setA, 
                                                                                                                             multialSNP.setB = processout$multialSNP.setB)
                                                                                             aliasout <- RegionScan:::aliasing(pheno = pheno, 
                                                                                                                               data = processout$data, binvector = prunout$binvector)
                                                                                           }
                                                                                           else {
                                                                                             aliasout <- RegionScan:::aliasing(pheno = pheno, 
                                                                                                                               data = processout$data, binvector = assignout)
                                                                                           }
                                                                                           codechange.MLC <- codechange.bial$newcode[names(codechange.bial$newcode) %in% 
                                                                                                                                       names(aliasout$final)]
                                                                                           mual <- setdiff(names(aliasout$final), names(codechange.bial$newcode))
                                                                                           binvector.LC <- rep(1, length(codechange.MLC))
                                                                                           names(binvector.LC) <- names(codechange.MLC)
                                                                                           codechange.LC <- RegionScan:::recoding(binlist = list(names(binvector.LC)), 
                                                                                                                                  data = processout$data)$newcode
                                                                                           if (length(mual) > 0) {
                                                                                             codechange.mual <- rep(0, length(mual))
                                                                                             names(codechange.mual) <- mual
                                                                                             codechange.MLC <- c(codechange.MLC, codechange.mual)
                                                                                             codechange.MLC <- codechange.MLC[names(aliasout$final)]
                                                                                             codechange.LC <- c(codechange.LC, codechange.mual)
                                                                                             codechange.LC <- codechange.LC[names(aliasout$final)]
                                                                                             binvector.LC <- rep(1, length(codechange.MLC))
                                                                                           }
                                                                                         }
                                                                                         else {
                                                                                           if (isTRUE(LDpruning)) {
                                                                                             prunout <- RegionScan:::pruning(binvector = clustout$binvector, 
                                                                                                                             data = processout$data, rcut = rcut)
                                                                                             aliasout <- RegionScan:::aliasing(pheno = pheno, 
                                                                                                                               data = processout$data, binvector = prunout$binvector)
                                                                                           }
                                                                                           else {
                                                                                             aliasout <- RegionScan:::aliasing(pheno = pheno, 
                                                                                                                               data = processout$data, binvector = clustout$binvector)
                                                                                           }
                                                                                           binlist.MLC <- lapply(unique(aliasout$final), 
                                                                                                                 function(bin) {
                                                                                                                   names(which(aliasout$final == bin))
                                                                                                                 })
                                                                                           codechange.MLC <- RegionScan:::recoding(binlist.MLC, 
                                                                                                                                   processout$data)$newcode
                                                                                           binlist.LC <- list(unlist(names(aliasout$final)))
                                                                                           binvector.LC <- rep(1, length(unlist(aliasout$final)))
                                                                                           names(binvector.LC) <- unlist(names(aliasout$final))
                                                                                           codechange.LC <- RegionScan:::recoding(binlist.LC, 
                                                                                                                                  processout$data)$newcode
                                                                                         }
                                                                                         nSNPs.kept <- length(aliasout$final)
                                                                                         if (isFALSE(singleSNPall)) {
                                                                                           sgout <- data.frame(do.call("rbind", lapply(names(aliasout$final), 
                                                                                                                                       function(x) {
                                                                                                                                         binvector <- 1
                                                                                                                                         names(binvector) <- x
                                                                                                                                         sglmout <- RegionScan:::glmfit(data = processout$data, 
                                                                                                                                                                        pheno = pheno, binvector = binvector, 
                                                                                                                                                                        covlist = covlist, family = family, 
                                                                                                                                                                        tol = tol, firthreg = firthreg)
                                                                                                                                         c(variant = x, sglm.beta = unname(sglmout$beta_g), 
                                                                                                                                           sglm.se = unname(sglmout$beta_SE), 
                                                                                                                                           sglm.pvalue = sglmout$p_g)
                                                                                                                                       })))
                                                                                         }
                                                                                         else {
                                                                                           if (length(assignout) > 0) {
                                                                                             binvectorsg <- assignout
                                                                                           }
                                                                                           else {
                                                                                             binvectorsg <- clustout$binvector
                                                                                           }
                                                                                           sgout <- data.frame(do.call("rbind", lapply(names(binvectorsg), 
                                                                                                                                       function(x) {
                                                                                                                                         binvector <- 1
                                                                                                                                         names(binvector) <- x
                                                                                                                                         sglmout <- RegionScan:::glmfit(data = processout$data, 
                                                                                                                                                                        pheno = pheno, binvector = binvector, 
                                                                                                                                                                        covlist = covlist, family = family, 
                                                                                                                                                                        tol = tol, firthreg = firthreg)
                                                                                                                                         c(variant = x, sglm.beta = unname(sglmout$beta_g), 
                                                                                                                                           sglm.se = unname(sglmout$beta_SE), 
                                                                                                                                           sglm.pvalue = sglmout$p_g)
                                                                                                                                       })))
                                                                                           outsingleSNPall <- data.frame(cbind(chr = chr, 
                                                                                                                               region = region, start.bp = start, end.bp = end, 
                                                                                                                               bin = binvectorsg, processout$SNPinfo[match(names(binvectorsg), 
                                                                                                                                                                           processout$SNPinfo$variant), c("bp", 
                                                                                                                                                                                                          "multiallelic", "ref", "alt", "maf")], 
                                                                                                                               sgout))
                                                                                         }
                                                                                         glmout <- RegionScan:::glmfit(data = processout$data, 
                                                                                                                       pheno = pheno, covlist = covlist, binvector = aliasout$final, 
                                                                                                                       family = family, tol = tol, firthreg = firthreg)
                                                                                         MLCBout <- RegionScan:::MLC(beta = glmout$beta_g, 
                                                                                                                     sigmainv = glmout$invcov_g, binvector = aliasout$final, 
                                                                                                                     codechange = codechange.MLC, tol = tol)
                                                                                         LCBout <- RegionScan:::MLC(beta = glmout$beta_g, 
                                                                                                                    sigmainv = glmout$invcov_g, binvector = binvector.LC, 
                                                                                                                    codechange = codechange.LC, tol = tol)
                                                                                         Waldout <- RegionScan:::Wald(glmout$beta_g, 
                                                                                                                      glmout$invcov_g)
                                                                                         PC80out <- RegionScan:::PC80(data = processout$data, 
                                                                                                                      pheno = pheno, covlist = covlist, binvector = aliasout$final, 
                                                                                                                      family = family, tol = tol)
                                                                                         SKATout <- RegionScan:::SKATp(data = processout$data, 
                                                                                                                       pheno = pheno, covlist = covlist, pheno_type = pheno_type, 
                                                                                                                       geno_type = geno_type, binvector = aliasout$final, 
                                                                                                                       SKAT_kernel = SKAT_kernel, SKAT_weights = SKAT_weights, 
                                                                                                                       SKAT_weights_beta = SKAT_weights_beta)
                                                                                         regionout <- data.frame(chr = chr, region = region, 
                                                                                                                 start.bp = start, end.bp = end, nSNPs = nSNPs, 
                                                                                                                 nSNPs.kept = nSNPs.kept, maxVIF = max(glmout$vif), 
                                                                                                                 Wald = Waldout$stat, Wald.df = Waldout$df, 
                                                                                                                 Wald.p = Waldout$pvalue, PC80 = PC80out$stat, 
                                                                                                                 PC80.df = PC80out$df, PC80.p = PC80out$pvalue, 
                                                                                                                 MLCB = MLCBout$stat, MLCB.df = MLCBout$df, 
                                                                                                                 MLCB.p = MLCBout$pvalue, LCB = LCBout$stat, 
                                                                                                                 LCB.df = LCBout$df, LCB.p = LCBout$pvalue, 
                                                                                                                 SKAT.p = SKATout$SKAT.p, SKATO.p = SKATout$SKATO.p)
                                                                                         if (isTRUE(nSNPs.kept > 1)) {
                                                                                           sgsel <- subset(sgout, variant %in% names(aliasout$final))
                                                                                           corsel <- cor(processout$data[, sgsel$variant, 
                                                                                                                         drop = F])
                                                                                           pval <- as.numeric(as.character(sgsel$sglm.pvalue))
                                                                                           names(pval) <- sgsel$variant
                                                                                           out_simes <- RegionScan:::ext_simes_v2(pval, 
                                                                                                                                  corsel)
                                                                                           out_simpleM <- RegionScan:::simpleM_v2(pval, 
                                                                                                                                  corsel)
                                                                                           out_gates <- RegionScan:::gates_v2(pval, 
                                                                                                                              corsel)
                                                                                           minp_uncorected <- min(pval)
                                                                                           regionout <- c(regionout, simes.p = out_simes$p, 
                                                                                                          simpleM.df = out_simpleM$simpleM.df, simpleM.p = out_simpleM$simpleM.p, 
                                                                                                          GATES.p = out_gates, single_Wald.p = minp_uncorected)
                                                                                         }
                                                                                         else {
                                                                                           regionout$maxVIF <- NA
                                                                                           regionout <- c(regionout, simes.p = glmout$p_g, 
                                                                                                          simpleM.df = 1, simpleM.p = glmout$p_g, 
                                                                                                          GATES.p = glmout$p_g, single_Wald.p = glmout$p_g)
                                                                                         }
                                                                                         if (length(assignout) > 0) {
                                                                                           bin.bfP <- lapply(unique(assignout), function(x) {
                                                                                             names(which(assignout == x))
                                                                                           })
                                                                                         }
                                                                                         else {
                                                                                           bin.bfP <- lapply(unique(clustout$binvector), 
                                                                                                             function(x) {
                                                                                                               names(which(clustout$binvector == x))
                                                                                                             })
                                                                                         }
                                                                                         binstart.bfP.bp <- unlist(lapply(bin.bfP, 
                                                                                                                          function(x) {
                                                                                                                            min(as.numeric(as.character(subset(processout$SNPinfo, 
                                                                                                                                                               variant %in% x)$bp)))
                                                                                                                          }))
                                                                                         binend.bfP.bp <- unlist(lapply(bin.bfP, function(x) {
                                                                                           max(as.numeric(as.character(subset(processout$SNPinfo, 
                                                                                                                              variant %in% x)$bp)))
                                                                                         }))
                                                                                         bin.afP <- lapply(unique(aliasout$final), 
                                                                                                           function(x) {
                                                                                                             names(which(aliasout$final == x))
                                                                                                           })
                                                                                         binstart.afP.bp <- unlist(lapply(bin.afP, 
                                                                                                                          function(x) {
                                                                                                                            min(as.numeric(as.character(subset(processout$SNPinfo, 
                                                                                                                                                               variant %in% x)$bp)))
                                                                                                                          }))
                                                                                         binend.afP.bp <- unlist(lapply(bin.afP, function(x) {
                                                                                           max(as.numeric(as.character(subset(processout$SNPinfo, 
                                                                                                                              variant %in% x)$bp)))
                                                                                         }))
                                                                                         binsize.bfP <- unlist(lapply(bin.bfP, length))
                                                                                         binsize.afP <- unlist(lapply(bin.afP, length))
                                                                                         binout <- data.frame(cbind(chr = chr, region = region, 
                                                                                                                    start.bp = start, end.bp = end, binstart.bfP.bp = binstart.bfP.bp, 
                                                                                                                    binend.bfP.bp = unname(binend.bfP.bp), binstart.afP.bp = binstart.afP.bp, 
                                                                                                                    binend.afP.bp = unname(binend.afP.bp), binsize.bfP = unname(binsize.bfP), 
                                                                                                                    binsize.afP = unname(binsize.afP), MLCBout$deltabin))
                                                                                         if (isTRUE(alltests)) {
                                                                                           MLCZout <- RegionScan:::MLC(Z = glmout$Z_g, 
                                                                                                                       invcor = glmout$invcor_g, binvector = aliasout$final, 
                                                                                                                       codechange = codechange.MLC, tol = tol)
                                                                                           LCZout <- RegionScan:::MLC(Z = glmout$Z_g, 
                                                                                                                      invcor = glmout$invcor_g, binvector = binvector.LC, 
                                                                                                                      codechange = codechange.LC, tol = tol)
                                                                                           regionout <- c(regionout, MLCZ = MLCZout$stat, 
                                                                                                          MLCZ.df = MLCZout$df, MLCZ.p = MLCZout$pvalue, 
                                                                                                          LCZ = LCZout$stat, LCZ.df = LCZout$df, 
                                                                                                          LCZ.p = LCZout$pvalue)
                                                                                           binout <- c(binout, MLCZout$deltabin)
                                                                                         }
                                                                                         snpout <- data.frame(chr = chr, region = region, 
                                                                                                              start.bp = start, end.bp = end, bin = aliasout$final, 
                                                                                                              processout$SNPinfo[match(names(aliasout$final), 
                                                                                                                                       processout$SNPinfo$variant), c("bp", "multiallelic", 
                                                                                                                                                                      "ref", "alt", "maf")], MLC.codechange = codechange.MLC[names(glmout$beta_g)], 
                                                                                                              LC.codechange = codechange.LC[names(glmout$beta_g)], 
                                                                                                              subset(sgout, variant %in% names(aliasout$final)), 
                                                                                                              mglm.vif = glmout$vif, mglm.beta = glmout$beta_g, 
                                                                                                              mglm.se = glmout$beta_SE, mglm.pvalue = glmout$p_g)
                                                                                         allremoved <- filterout <- NULL
                                                                                         if (length(processout$removed.mafcut) > 0) {
                                                                                           allremoved <- data.frame(cbind(variant = processout$removed.mafcut, 
                                                                                                                          bin = "NA", reason = "mafcut"))
                                                                                         }
                                                                                         if (length(processout$removed.multiallelic) > 
                                                                                             0) {
                                                                                           allremoved <- data.frame(rbind(allremoved, 
                                                                                                                          cbind(variant = processout$removed.multiallelic, 
                                                                                                                                bin = "NA", reason = "multial")))
                                                                                         }
                                                                                         if (length(names(prunout$removed)) > 0) {
                                                                                           allremoved <- data.frame(rbind(allremoved, 
                                                                                                                          cbind(variant = names(prunout$removed), 
                                                                                                                                bin = prunout$removed, reason = "rcut")))
                                                                                         }
                                                                                         if (length(names(aliasout$aliasout)) > 0) {
                                                                                           allremoved <- data.frame(rbind(allremoved, 
                                                                                                                          cbind(variant = names(aliasout$aliasout), 
                                                                                                                                bin = aliasout$aliasout, reason = "alias")))
                                                                                         }
                                                                                         if (!is.null(allremoved)) {
                                                                                           bp <- SNPinfo_region[match(allremoved$variant, 
                                                                                                                      SNPinfo_region$variant), ]$bp
                                                                                           names(bp) <- SNPinfo_region[match(allremoved$variant, 
                                                                                                                             SNPinfo_region$variant), ]$variant
                                                                                           filterout <- cbind(chr = chr, region = region, 
                                                                                                              start = start, end = end, bp = bp, allremoved)
                                                                                         }
                                                                                         if (!is.null(filterout) && nrow(filterout) > 
                                                                                             1) {
                                                                                           filterout <- filterout[order(filterout$bin, 
                                                                                                                        filterout$bp), ]
                                                                                         }
                                                                                         if (nrow(binout) > 1) {
                                                                                           binout <- binout[order(binout$bin), ]
                                                                                         }
                                                                                         if (nrow(snpout) > 1) {
                                                                                           snpout <- snpout[order(snpout$bin, snpout$bp), 
                                                                                           ]
                                                                                         }
                                                                                         if (!is.null(outsingleSNPall)) {
                                                                                           outsingleSNPall <- outsingleSNPall[order(outsingleSNPall$bin, 
                                                                                                                                    outsingleSNPall$bp), ]
                                                                                         }
                                                                                         if (isTRUE(MLCheatmap)) {
                                                                                           if (length(assignout) > 0) {
                                                                                             binvector <- assignout
                                                                                           }
                                                                                           else {
                                                                                             binvector <- clustout$binvector
                                                                                           }
                                                                                           nbins <- max(binvector)
                                                                                           colvect <- (scales::hue_pal())(nbins)
                                                                                           outfile <- paste(path,"MLC_heatmap_chr", chr, 
                                                                                                            ".region", region, sep = "")
                                                                                           tmp0 <- subset(processout$SNPinfo, variant %in% 
                                                                                                            names(binvector))
                                                                                           tmp0 <- tmp0[order(tmp0$bp), ]
                                                                                           cormat1 <- apply(cor(processout$data[, tmp0$variant, 
                                                                                                                                drop = F]), 2, function(x) {
                                                                                                                                  round(x, 2)
                                                                                                                                })
                                                                                           tdata1 <- processout$data[, tmp0$variant, 
                                                                                                                     drop = F]
                                                                                           tdata1[, names(which(codechange.MLC == 1))] <- 2 - 
                                                                                             tdata1[, names(which(codechange.MLC == 
                                                                                                                    1))]
                                                                                           cormat <- apply(cor(tdata1), 2, function(x) {
                                                                                             round(x, 2)
                                                                                           })
                                                                                           tmp <- subset(processout$SNPinfo, variant %in% 
                                                                                                           names(aliasout$final))
                                                                                           tmp <- tmp[order(tmp$bp), ]
                                                                                           tdata2 <- processout$data[, names(aliasout$final)]
                                                                                           tdata2[, names(which(codechange.MLC == 1))] <- 2 - 
                                                                                             tdata2[, names(which(codechange.MLC == 
                                                                                                                    1))]
                                                                                           cormat3 <- apply(cor(tdata2), 2, function(x) {
                                                                                             round(x, 2)
                                                                                           })
                                                                                           RegionScan:::RegionHeatmap(cormat1, type_mat = "Correlation,\n", 
                                                                                                                      ptitle = "Within region correlation (before pruning & recoding),\n SNPs ordered by pos", 
                                                                                                                      colvect = colvect, clustorder = F)
                                                                                           ggplot2::ggsave(paste(outfile, "_before_pruning_and_recoding_ordered_by_pos.pdf", 
                                                                                                                 sep = ""))
                                                                                           RegionScan:::RegionHeatmap(cormat1[names(binvector), 
                                                                                                                              names(binvector)], binvector = binvector, 
                                                                                                                      type_mat = "Correlation", ptitle = "Within region correlation (before pruning & recoding),\n SNPs ordered by LDbin", 
                                                                                                                      colvect = colvect, clustorder = F)
                                                                                           ggplot2::ggsave(paste(outfile, "_before_pruning_and_recoding_ordered_by_LDbin.pdf", 
                                                                                                                 sep = ""))
                                                                                           RegionScan:::RegionHeatmap(cormat1[tmp$variant, 
                                                                                                                              tmp$variant], type_mat = "Correlation", 
                                                                                                                      ptitle = "Within region correlation (after pruning & before recoding),\n SNPs ordered by pos", 
                                                                                                                      colvect = colvect, clustorder = F)
                                                                                           ggplot2::ggsave(paste(outfile, "_after_pruning_and_before_recoding_ordered_by_pos.pdf", 
                                                                                                                 sep = ""))
                                                                                           RegionScan:::RegionHeatmap(cormat1[names(aliasout$final), 
                                                                                                                              names(aliasout$final)], type_mat = "Correlation", 
                                                                                                                      ptitle = "Within region correlation (after pruning & before recoding),\n SNPs ordered by LDbin", 
                                                                                                                      binvector = aliasout$final, colvect = colvect, 
                                                                                                                      clustorder = F)
                                                                                           ggplot2::ggsave(paste(outfile, "_after_pruning_and_before_recoding_ordered_by_LDbin.pdf", 
                                                                                                                 sep = ""))
                                                                                           RegionScan:::RegionHeatmap(cormat3[tmp$variant, 
                                                                                                                              tmp$variant], type_mat = "Correlation", 
                                                                                                                      ptitle = "Within region correlation (after pruning & recoding),\n SNPs ordered by pos", 
                                                                                                                      colvect = colvect, clustorder = F)
                                                                                           ggplot2::ggsave(paste(outfile, "_after_pruning_and_after_recoding_ordered_by_pos.pdf", 
                                                                                                                 sep = ""))
                                                                                           RegionScan:::RegionHeatmap(cormat3, type_mat = "Correlation", 
                                                                                                                      ptitle = "Within region correlation (after pruning & recoding),\n SNPs ordered by LDbin", 
                                                                                                                      binvector = aliasout$final, colvect = colvect, 
                                                                                                                      clustorder = F)
                                                                                           ggplot2::ggsave(paste(outfile, "_after_pruning_and_after_recoding_ordered_by_LDbin.pdf", 
                                                                                                                 sep = ""))
                                                                                         }
                                                                                       }
                                                                                       else {
                                                                                         binvector <- 1
                                                                                         names(binvector) <- processout$SNPinfo[, "variant"]
                                                                                         mglmtime <- glmout <- RegionScan:::glmfit(data = processout$data, 
                                                                                                                                   pheno, binvector = binvector, covlist = covlist, 
                                                                                                                                   family = family, tol = tol, firthreg = firthreg)
                                                                                         MLCBout <- RegionScan:::MLC(beta = glmout$beta_g, 
                                                                                                                     sigmainv = glmout$invcov_g, binvector = binvector, 
                                                                                                                     codechange = 0, tol = tol)
                                                                                         PC80out <- Waldout <- LCBout <- MLCBout
                                                                                         SKATout <- RegionScan:::SKATp(data = processout$data, 
                                                                                                                       pheno = pheno, covlist = covlist, pheno_type = pheno_type, 
                                                                                                                       geno_type = geno_type, binvector = binvector, 
                                                                                                                       SKAT_kernel = SKAT_kernel, SKAT_weights = SKAT_weights, 
                                                                                                                       SKAT_weights_beta = SKAT_weights_beta)
                                                                                         sgout <- data.frame(variant = processout$SNPinfo[, 
                                                                                                                                          "variant"], sglm.beta = unname(glmout$beta_g), 
                                                                                                             sglm.se = unname(glmout$beta_SE), sglm.pvalue = unname(glmout$p_g))
                                                                                         binout <- data.frame(cbind(chr = chr, region = region, 
                                                                                                                    start.bp = start, end.bp = end, binstart.bfP.bp = unname(unique(processout$SNPinfo["bp"])), 
                                                                                                                    binend.bfP.bp = unname(unique(processout$SNPinfo["bp"])), 
                                                                                                                    binstart.afP.bp = unname(unique(processout$SNPinfo["bp"])), 
                                                                                                                    binend.afP.bp = unname(unique(processout$SNPinfo["bp"])), 
                                                                                                                    binsize.bfP = 1, binsize.afP = 1, unique(MLCBout$deltabin)))
                                                                                         snpout <- data.frame(chr = chr, region = region, 
                                                                                                              start.bp = start, end.bp = end, bin = 1, 
                                                                                                              processout$SNPinfo[, c("bp", "multiallelic", 
                                                                                                                                     "ref", "alt", "maf")], MLC.codechange = 0, 
                                                                                                              LC.codechange = 0, sgout, mglm.vif = "NA", 
                                                                                                              mglm.beta = unname(unlist(glmout$beta_g)), 
                                                                                                              mglm.se = glmout$beta_SE, mglm.pvalue = glmout$p_g)
                                                                                         regionout <- data.frame(chr = chr, region = region, 
                                                                                                                 start.bp = start, end.bp = end, nSNPs = nSNPs, 
                                                                                                                 nSNPs.kept = nSNPs.kept, maxVIF = "NA", 
                                                                                                                 Wald = Waldout$stat, Wald.df = Waldout$df, 
                                                                                                                 Wald.p = Waldout$pvalue, PC80 = PC80out$stat, 
                                                                                                                 PC80.df = PC80out$df, PC80.p = PC80out$pvalue, 
                                                                                                                 MLCB = MLCBout$stat, MLCB.df = MLCBout$df, 
                                                                                                                 MLCB.p = MLCBout$pvalue, LCB = LCBout$stat, 
                                                                                                                 LCB.df = LCBout$df, LCB.p = LCBout$pvalue, 
                                                                                                                 SKAT.p = SKATout$SKAT.p, SKATO.p = SKATout$SKATO.p, 
                                                                                                                 simes.p = glmout$p_g, simpleM.df = 1, simpleM.p = glmout$p_g, 
                                                                                                                 GATES.p = glmout$p_g, single_Wald.p = glmout$p_g)
                                                                                         if (isTRUE(alltests)) {
                                                                                           MLCZout <- LCZout <- RegionScan:::MLC(Z = glmout$Z_g, 
                                                                                                                                 invcor = glmout$invcor_g, binvector = binvector, 
                                                                                                                                 codechange = 0, tol = tol)
                                                                                           regionout <- c(regionout, MLCZ = MLCZout$stat, 
                                                                                                          MLCZ.df = MLCZout$df, MLCZ.p = MLCZout$pvalue, 
                                                                                                          LCZ = LCZout$stat, LCZ.df = LCZout$df, 
                                                                                                          LCZ.p = LCZout$pvalue)
                                                                                           binout <- c(binout, MLCZout$deltabin)
                                                                                         }
                                                                                       }
                                                                                       if (isTRUE(covout) && length(covlist) > 0) {
                                                                                         covoutput <- data.frame(cbind(chr = chr, region = region, 
                                                                                                                       start.bp = start, end.bp = end, covlist = covlist, 
                                                                                                                       mglm_beta = glmout$beta_covlist, mglm_pvalue = glmout$p_covlist))
                                                                                       }
                                                                                     }
                                                                                   }, error = function(e) {
                                                                                     cat("ERROR :", conditionMessage(e), "\n")
                                                                                   })
                                                                                   return(list(regionout = regionout, binout = binout, 
                                                                                               snpout = snpout, filterout = filterout, outsingleSNPall = outsingleSNPall, 
                                                                                               covarout = covoutput))
                                                                                 }
  if (length(rscan) > 6) {
    regionout <- as.data.frame(do.call("rbind", lapply(rscan[, 
                                                             "regionout"], function(x) {
                                                               unlist(x)
                                                             })))
    binout <- as.data.frame(do.call("rbind", rscan[, "binout"]))
    snpout <- as.data.frame(do.call("rbind", rscan[, "snpout"]))
    filterout <- as.data.frame(do.call("rbind", rscan[, 
                                                      "filterout"]))
    outsingleSNPall <- as.data.frame(do.call("rbind", rscan[, 
                                                            "outsingleSNPall"]))
    covarout <- as.data.frame(do.call("rbind", rscan[, "covarout"]))
  }
  else {
    regionout <- as.data.frame(do.call("rbind", lapply(rscan["regionout"], 
                                                       function(x) {
                                                         unlist(x)
                                                       })))
    binout <- as.data.frame(do.call("rbind", rscan["binout"]))
    snpout <- as.data.frame(do.call("rbind", rscan["snpout"]))
    filterout <- as.data.frame(do.call("rbind", rscan["filterout"]))
    outsingleSNPall <- as.data.frame(do.call("rbind", rscan["outsingleSNPall"]))
    covarout <- as.data.frame(do.call("rbind", rscan["covarout"]))
  }
  regionout$chr <- as.numeric(as.character(regionout$chr))
  if (isTRUE(parallel)) {
    parallel::stopCluster(cl = my.cluster)
  }
  rownames(regionout) <- rownames(binout) <- rownames(snpout) <- rownames(filterout) <- rownames(outsingleSNPall) <- rownames(covarout) <- NULL
  return(list(regionout = regionout, binout = binout, snpout = snpout, 
              filterout = filterout, outsingleSNPall = outsingleSNPall, 
              covout = covarout))
}


# LD heatmaps 
regscan(phenocov = phenocov, pheno="sim_QT", REGIONinfo=REGIONinfo,
        geno_type="D", pheno_type="C", data = geno, SNPinfo = SNPinfo,
        MLCheatmap = TRUE, regionlist ="8",path="output/RegionScan/")


MLCbinsnpPlot <- function (rscanout, chr_, region_,path=NULL) # MLCbinPlot with path
{
  snpout <- rscanout$snpout[, c("chr", "region", "bin", "variant", 
                                "bp")]
  filterout <- rscanout$filterout[, c("chr", "region", "bin", 
                                      "variant", "bp")]
  binout <- rscanout$binout[, c("chr", "region", "start.bp", 
                                "end.bp", "bin", "binstart.afP.bp", "binend.afP.bp", 
                                "binstart.bfP.bp", "binend.bfP.bp", "binsize.bfP", "binsize.afP")]
  snpout <- subset(snpout, chr == chr_ & region == region_)
  filterout <- na.omit(subset(filterout, chr == chr_ & region == 
                                region_))
  binout <- subset(binout, chr == chr_ & region == region_)
  allout <- unique(cbind(rbind(snpout, filterout), kept = c(rep(1, 
                                                                nrow(snpout)), rep(0, nrow(filterout)))))
  nbins <- max(as.numeric(as.character(snpout$bin)))
  colvect <- (scales::hue_pal())(nbins)
  binout <- binout[order(as.numeric(as.character(binout$bin))), 
  ]
  binout$binstart.afP.bp <- as.numeric(as.character(binout$binstart.afP.bp))
  binout$binend.afP.bp <- as.numeric(as.character(binout$binend.afP.bp))
  binout$binstart.bfP.bp <- as.numeric(as.character(binout$binstart.bfP.bp))
  binout$binend.bfP.bp <- as.numeric(as.character(binout$binend.bfP.bp))
  binout$bin <- factor(binout$bin, level = unique(binout$bin), 
                       order = TRUE)
  nSNP.bfP <- sum(as.numeric(as.character(binout$binsize.bfP)))
  nSNP.afP <- sum(as.numeric(as.character(binout$binsize.afP)))
  snpout$bin <- factor(snpout$bin, level = unique(snpout$bin), 
                       order = TRUE)
  snpout$bp <- as.numeric(as.character(snpout$bp))
  outname <- paste(path,"MLC_LDbin_SNPpos_chr", chr_, "_region", 
                   region_, sep = "")
  ggplot2::ggplot() + ggplot2::geom_segment(data = binout, 
                                            aes(x = binstart.afP.bp, y = bin, xend = binend.afP.bp, 
                                                yend = bin, color = bin)) + ggplot2::geom_point(data = snpout, 
                                                                                                aes(x = bp, y = bin, color = bin)) + ggplot2::ggtitle(paste("region: ", 
                                                                                                                                                            region_, ", after pruning : ", nSNP.afP, " SNP(s)", 
                                                                                                                                                            sep = "")) + ggplot2::xlab("pos") + ggplot2::ylab("MLC LD bin#") + 
    ggplot2::scale_color_manual(values = colvect, labels = paste("Bin", 
                                                                 binout$bin, ": ", binout$binsize.afP, " SNP(s)", 
                                                                 sep = ""))
  ggplot2::ggsave(paste(outname, "_after_pruning.pdf", sep = ""), 
                  width = 15, height = 10, units = "cm")
  ggplot2::ggplot() + ggplot2::geom_segment(data = binout, 
                                            aes(x = binstart.bfP.bp, y = bin, xend = binend.bfP.bp, 
                                                yend = bin, color = bin)) + ggplot2::geom_point(data = subset(allout, 
                                                                                                              kept == 0), aes(x = bp, y = bin), color = "darkgrey") + 
    ggplot2::geom_point(data = subset(allout, kept == 1), 
                        aes(x = bp, y = bin, color = bin)) + ggplot2::ggtitle(paste("region", 
                                                                                    region_, ", after (before) pruning : ", nSNP.afP, " (", 
                                                                                    nSNP.bfP, ") ", " SNPs", sep = "")) + ggplot2::xlab("pos") + 
    ggplot2::ylab("MLC LD bin#") + ggplot2::scale_color_manual(values = colvect, 
                                                               labels = paste("Bin", binout$bin, ": ", binout$binsize.afP, 
                                                                              " (", binout$binsize.bfP, ") SNP(s)", sep = ""))
  ggplot2::ggsave(paste(outname, "_before_and_after_pruning.pdf", 
                        sep = ""), width = 15, height = 10, units = "cm")
}


# SNP LD bin position
MLCbinsnpPlot(rscanout = results, chr_=16, region_ =8, path="output/RegionScan/")

