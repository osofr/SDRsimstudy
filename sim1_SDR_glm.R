`%+%` <- function(a, b) paste0(a, b)
require("data.table")
require("simcausal")
require("foreach")
require("doParallel")
require("magrittr")
options(simcausal.verbose = FALSE)

# nsims <- 50
nsims <- 1000
nsamp <- 500
tvals <- 2 # tvals <- 1

# g model (same for both time-points)
gform.c <- "A ~ L + Atm1"
gform.i <- "A ~ L"

# Q models
# Qforms.c <- c("Qkplus1 ~ L + A", "Qkplus1 ~ A0L1 + L0A0 + L + A + Atm1 + Ltm1")
# Qforms.c <- c("Qkplus1 ~ L",
#               "Qkplus1 ~ A0L1 + L0A0 + L + A + Atm1 + Ltm1")
Qforms.c <- c("Qkplus1 ~ L",
              "Qkplus1 ~ A0L1 + L0A0 + L + A + Atm1 + Ltm1",
              "Qkplus1 ~ A0L1 + L0A0 + L + A + Atm1 + Ltm1")
Qinteract.c <- list(c("A", "L"), c("A", "Ltm1"))
# Qforms.c <- c("Qkplus1 ~ L0A1 + A0L1 + L1A1 + L0A0 + L + A + Atm1 + Ltm1", "Qkplus1 ~ L + A")
# Qforms.i <- c("Qkplus1 ~ L + A", "Qkplus1 ~ L + A + Atm1 + Ltm1")
# Qforms.i <- c("Qkplus1 ~ L", "Qkplus1 ~ Atm1")
Qforms.i <- c("Qkplus1 ~ L", "Qkplus1 ~ Atm1", "Qkplus1 ~ Atm1")
Qinteract.i <- NULL

registerDoParallel(cores=detectCores())
## ----------------------------------------------------------------------------
## simcausal DAG
## ----------------------------------------------------------------------------
simDAG_t1 = function(rndseed) {
  require("simcausal")
  D <- DAG.empty() +
      node("L", t = 0, distr = "rnorm") +
      node("A", t = 0, distr = "rbern", prob = plogis(L[0])) +
      node("Y", t = 0, distr = "rconst", const = 0) +
      node("L", t = 1, distr = "rnorm") +
      node("A", t = 1, distr = "rbern", prob = plogis(L[1] + A[0])) +
      node("Y", t = 1, distr = "rbern", prob = plogis(L[0] * A[1] + A[0]*L[1] + L[1]*A[1]))

  Dset <-
      set.DAG(D) +
      action("A1is1", nodes =
          c(node("A", t = 1, distr = "rbern", prob = 1),
              node("Y", t = 1, distr = "rconst", const = plogis(L[0]*A[1] + A[0]*L[1] + L[1]*A[1])))
          ) +
      action("allAare1", nodes =
          c(node("A", t = c(0,1), distr = "rbern", prob = 1),
              node("Y", t = 1, distr = "rconst", const = plogis(L[0]*A[1] + A[0]*L[1] + L[1]*A[1])))
          )
  return(Dset)
}

simDAG_t2 = function(rndseed) {
  require("simcausal")
  D <- DAG.empty() +
      node("L", t = 0, distr = "rnorm") +
      node("A", t = 0, distr = "rbern", prob = plogis(L[0])) +
      node("Y", t = 0, distr = "rconst", const = 0) +
      node("L", t = 1, distr = "rnorm") +
      node("A", t = 1, distr = "rbern", prob = plogis(L[1] + A[0])) +
      node("Y", t = 1, distr = "rconst", const = 0) +
      node("L", t = 2, distr = "rnorm", mean = L[0]*A[1] + A[0]*L[1] + L[1]*A[1]) +
      node("A", t = 2, distr = "rbern", prob = plogis(L[2] + A[1])) +
      node("Y", t = 2, distr = "rbern", prob = plogis(L[1]*A[2] + A[1]*L[2] + L[2]*A[2]))

  Dset <-
      set.DAG(D) +
      action("A1A2are1", nodes =
          c(node("A", t = c(1,2), distr = "rbern", prob = 1),
              node("Y", t = 2, distr = "rconst", const = plogis(L[1]*A[2] + A[1]*L[2] + L[2]*A[2])))
          ) +
      action("allAare1", nodes =
          c(node("A", t = c(0,1,2), distr = "rbern", prob = 1),
              node("Y", t = 2, distr = "rconst", const = plogis(L[1]*A[2] + A[1]*L[2] + L[2]*A[2])))
          )
  return(Dset)
}

eval_true_EY_WA <- function(dag_obj, action_name = "allAare1", tval = 2, ntest = 1e4, mc.repl = 1000, rndseed = 12345, rndseed.reset.node = "Y_0") {
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  require("doParallel")
  require("foreach")
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  ex <- c(ex, "nfull", "simDAG_t1", "simDAG_t2")
  Dset <- dag_obj
  dat0_repl <- foreach(t.counter = (1 : mc.repl), .packages = c("simcausal"), .export = ex, .combine = "+") %dopar% {
    dat <- data.table::data.table(sim(Dset, actions = action_name, n = ntest, rndseed = rndseed, rndseed.reset.node = rndseed.reset.node)[[1]])
    dat[[c("Y_" %+% tval)]]
  }
  P0_EY_WA <- dat0_repl / mc.repl
  dat0 <- data.table::data.table(sim(Dset, actions = action_name, n = ntest, rndseed = rndseed)[[1]])
  test.dat <- dat0[, c("L_0","A_0"), with = FALSE]
  setnames(test.dat, c("L", "A"))
  test.dat[, ("ID") := seq.int(nrow(test.dat))]
  test.dat[, ("t") := 0]
  return(list(test.dat = test.dat, P0_EY_WA = P0_EY_WA))
}

dag_obj <- do.call("simDAG_t" %+% tvals, args = list())
# dag_obj <- simDAG_t1(); dag_obj <- simDAG_t2()
truth_psi0 <- mean(sim(dag_obj, actions = "allAare1", n = 2000000)[[1]][["Y_"%+%tvals]])
# truth_EY_WA <- eval_true_EY_WA(dag_obj, ntest = 1e4, tval = tvals, mc.repl = 100, rndseed = 12345)
truth_EY_WA <- eval_true_EY_WA(dag_obj, ntest = 1e4, tval = tvals, mc.repl = 20000, rndseed = 12345)
truth_EY_WA[["truth_psi0"]] <- truth_psi0
truth_EY_WA[["test.dat"]]
print(truth_EY_WA[["test.dat"]])


sdr_tmle_gcomp_1sim <- function(dat, truth, tvals = 2, Qforms, gforms, Qinteract, ...) {
  library("data.table")
  setDTthreads(1)
  library("xgboost")
  # library("gridisl")
  library("stremr")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)

  # truth_psi0 <- psi0_byt[tvals]
  truth_psi0 <- truth[["truth_psi0"]]
  test.dat <- truth[["test.dat"]]
  test.dat <- data.table(test.dat)
  P0_EY_WA <- truth[["P0_EY_WA"]]

  ID <- "ID"
  t <- "t"
  TRT <- "A"
  outcome <- "Y"
  # covars <- c("L", "Ltm1", "Atm1", "L0A1", "A0L1", "L1A1")
  covars <- c("L", "Ltm1", "Atm1", "A0L1", "L0A0")

  # Define counterfactual (intervention) node
  # dat[, ("A.star") := A][t == 1, ("A.star") := 1]
  dat[, ("A.star") := 1]

  # Define lag nodes and interactions.
  setkeyv(dat, cols = c("ID", "t"))
  lagnodes <- c("L", "A")
  newVarnames <- paste0(lagnodes, "tm1")
  dat[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # dat[, ("L0A1") := A * Ltm1]
  dat[, ("A0L1") := Atm1 * L]
  # dat[, ("L1A1") := A * L]
  dat[, ("L0A0") := Atm1 * Ltm1]
  # Dummies for censoring and monitoring
  dat[, ("C") := 0][, ("N") := 1]

  # Propensity score models for Treatment, Censoring & Monitoring:
  gform_TRT <- gforms
  # stratify_TRT <- list(A=c("t == 0L","t == 1L"))  # Separate models for t=0 and t=1
  stratify_TRT <- list(A=c("t == 0L","t == 1L", "t == 2L"))  # Separate models for t=0, t=1 and t=2
  gform_CENS <- "C ~ 1"
  gform_MONITOR <- "N ~ 1"

  # IMPORT DATA
  OData <- stremr::importData(dat, ID = ID, t = t, covars = covars, CENS = "C", TRT = TRT, MONITOR = "N", OUTCOME = "Y")
  OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)

  # FIT PROPENSITY SCORES WITH speedglm
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                          estimator = "speedglm__glm",
                          family = "quasibinomial", # estimator = "xgboost__gbm",
                          fit_method = "none", # fit_method = "cv",
                          # fold_column = "fold_ID",
                          )

  analysis <- list(intervened_TRT = c("A.star"),
                   stratifyQ_by_rule = TRUE) %>%
                purrr::cross_d() %>%
                dplyr::arrange(stratifyQ_by_rule)
  # params <- gridisl::defModel(estimator = "xgboost__gbm",
  #                             family = "quasibinomial",
  #                             nthread = 2,
  #                             nrounds = 100,
  #                             early_stopping_rounds = 2,
  #                             interactions = interactions)
  paramsGLM <-
      gridisl::defModel(estimator = "xgboost__glm",
                      family = "quasibinomial",
                      nthread = 1,
                      nrounds = 200,
                      # max_delta_step = 2,
                      interactions = Qinteract
                      )

  eval_est_perf <- function(est_obj, est = "sdr") {
    est_name <- attributes(est_obj)$estimator_short
    psi_n <- (1 - est_obj[["estimates"]][["St."%+%est_name]])
    est_bias <- (truth_psi0-psi_n)
    fW_fit <- est_obj[["estimates"]][["fW_fit"]][[1]]
    preds_fW <- gridisl::predict_SL(fW_fit, test.dat)
    est_MSE <- mean((P0_EY_WA - preds_fW[["preds"]])^2)
    est_SE <- est_obj[["estimates"]][["SE."%+%est_name]]
    est_CI <- c(psi_n - abs(qnorm(0.025))*est_SE, psi_n + abs(qnorm(0.025))*est_SE)
    CIlen <- abs(est_CI[2]-est_CI[1])
    est_cover <- as.integer((truth_psi0 <= est_CI[2]) && (truth_psi0 >= est_CI[1]))
    # print("psi_n: "); print(psi_n); print("est_CI: "); print(est_CI)
    res <- c(est_cover, est_MSE, est_bias, psi_n, CIlen)
    names(res) <- c(est%+%"_cover", est%+%"_mse", est%+%"_bias", est%+%"_psi_n", est%+%"_CIlen")
    fW_n <- preds_fW[["preds"]]
    return(list(res = res, fW_n = fW_n, est_CI = est_CI))
  }

  # IPW <-  analysis %>%
  #         dplyr::distinct(intervened_TRT) %>%
  #         dplyr::mutate(wts_data = purrr::map(intervened_TRT, getIPWeights, OData = OData)) %>%
  #         dplyr::mutate(NPMSM = purrr::map(wts_data, ~ survNPMSM(wts_data = .x, OData = OData))) %>%
  #         dplyr::mutate(NPMSM = purrr::map(NPMSM, "estimates"))
  # St.IPW <- IPW %>% dplyr::select(NPMSM) %>% tidyr::unnest()
  # ipw_bias <- (data.table::as.data.table(St.IPW)[time == tvals, ][["St.NPMSM"]] - truth_psi0)

  ## direct IPW w/ GLM g
  IPW <-  analysis %>%
          dplyr::distinct(intervened_TRT) %>%
          dplyr::mutate(wts_data = purrr::map(intervened_TRT, getIPWeights, OData = OData)) %>%
          dplyr::mutate(IPW = purrr::map(wts_data, ~ survDirectIPW(wts_data = .x, OData = OData))) %>%
          dplyr::mutate(IPW = purrr::map(IPW, "estimates"))
  St.IPW <- IPW %>% dplyr::select(IPW) %>% tidyr::unnest()
  ipw_psi_n <- 1 - data.table::as.data.table(St.IPW)[time == tvals, ][["St.DirectIPW"]]
  ipw_bias <- (truth_psi0-ipw_psi_n)
  ipw_res <- c(ipw_bias, ipw_psi_n)
  names(ipw_res) <- c("ipw_bias", "ipw_psi_n")
  print("ipw_res"); print(ipw_res)

  #GCOMP
  # cat("fitting SeqGcomp\n")
  gcomp_est <- fitSeqGcomp(OData, tvals = tvals,
                           intervened_TRT = "A.star", Qforms = Qforms,
                           # models = params,
                           stratifyQ_by_rule = TRUE,
                           estimator = "speedglm__glm",
                           fit_method = "none",
                           # fit_method = "cv",
                           # fold_column = "fold_ID",
                           return_fW = TRUE,
                           interactions = Qinteract)
  gcomp_res <- eval_est_perf(gcomp_est, est = "gcomp")
  print("gcomp_res"); print(gcomp_res[["res"]])
  # psi_n_gcomp <- gcomp_est[["estimates"]][["St.GCOMP"]]
  # gcomp_bias <- (truth_psi0-psi_n_gcomp)
  # fW_fit <- gcomp_est[["estimates"]][["fW_fit"]][[1]]
  # preds_fW <- gridisl::predict_SL(fW_fit, test.dat)
  # gcomp_MSE <- mean((P0_EY_WA - preds_fW[["preds"]])^2)

  #TMLE
  tmle_est <- fitTMLE(OData, tvals = tvals,
                      intervened_TRT = "A.star", Qforms = Qforms,
                      stratifyQ_by_rule = TRUE,
                      estimator = "speedglm__glm",
                      fit_method = "none",
                      # fit_method = "cv",
                      # models = params,
                      # fold_column = "fold_ID",
                      return_fW = TRUE,
                      interactions = Qinteract)
  tmle_res <- eval_est_perf(tmle_est, est = "tmle")
  print("tmle_res"); print(tmle_res[["res"]])
  # psi_n_tmle <- tmle_est[["estimates"]][["St.TMLE"]]
  # tmle_bias <- (truth_psi0-psi_n_tmle)
  # fW_fit <- tmle_est[["estimates"]][["fW_fit"]][[1]]
  # preds_fW <- gridisl::predict_SL(fW_fit, test.dat)
  # tmle_MSE <- mean((P0_EY_WA - preds_fW[["preds"]])^2)
  # tmle_SE <- tmle_est[["estimates"]][["SE.TMLE"]]
  # tmle_CI <- c(psi_n_tmle - abs(qnorm(0.025))*tmle_SE, psi_n_tmle + abs(qnorm(0.025))*tmle_SE)
  # tmle_cover <- as.integer((truth_psi0 <= tmle_CI[2]) && (truth_psi0 >= tmle_CI[1]))

  # SDR
  # cat("fitting SeqDR\n")
  SDR_est <- stremr:::fitSDR(OData, tvals = tvals,
                      intervened_TRT = "A.star", Qforms = Qforms,
                      stratifyQ_by_rule = TRUE,
                      estimator = "speedglm__glm",
                      fit_method = "none",
                      # fit_method = "cv",
                      # models = params,
                      fold_column = "fold_ID",
                      return_fW = TRUE,
                      interactions = Qinteract)
  sdr_res <- eval_est_perf(SDR_est, est = "sdr")
  print("sdr_res"); print(sdr_res[["res"]])
  # psi_n_sdr <- SDR_est[["estimates"]][["St.SDR"]]
  # sdr_bias <- (truth_psi0-psi_n_sdr)
  # fW_fit <- SDR_est[["estimates"]][["fW_fit"]][[1]]
  # preds_fW <- gridisl::predict_SL(fW_fit, test.dat)
  # sdr_MSE <- mean((P0_EY_WA - preds_fW[["preds"]])^2)
  # sdr_SE <- SDR_est[["estimates"]][["SE.SDR"]]
  # sdr_CI <- c(psi_n_sdr - abs(qnorm(0.025))*sdr_SE, psi_n_sdr + abs(qnorm(0.025))*sdr_SE)
  # sdr_cover <- as.integer((truth_psi0 <= sdr_CI[2]) && (truth_psi0 >= sdr_CI[1]))

  DR_trans_est <- stremr:::fitSDR(OData, tvals = tvals,
                           intervened_TRT = "A.star", Qforms = Qforms,
                           stratifyQ_by_rule = TRUE,
                           fit_method = "none",
                           # estimator = "speedglm__glm",
                           models = paramsGLM,
                           fold_column = "fold_ID",
                           return_fW = TRUE,
                           use_DR_transform = TRUE # stabilize = FALSE,
                           # interactions = Qinteract
                          )
  drtrans_res <- eval_est_perf(DR_trans_est, est = "drtrans")
  print("drtrans_res"); print(drtrans_res[["res"]])

  # cat("sdr_bias: ", sdr_bias, "\n")
  # cat("tmle_bias: ", tmle_bias, "\n")
  # cat("gcomp_bias: ", gcomp_bias, "\n")
  # cat("ipw_bias: ", ipw_bias, "\n")
  # cat("sdr_MSE: ", sdr_MSE, "\n")
  # cat("tmle_MSE: ", tmle_MSE, "\n")
  # cat("gcomp_MSE: ", gcomp_MSE, "\n")
  # cat("sdr_CI: ", sdr_CI, "\n")
  # cat("tmle_CI: ", tmle_CI, "\n")

  # return(c(sdr_cover = sdr_cover, tmle_cover = tmle_cover,
  #          sdr_MSE = sdr_MSE, tmle_MSE = tmle_MSE, gcomp_MSE = gcomp_MSE,
  #          sdr_bias = sdr_bias, tmle_bias = tmle_bias, gcomp_bias = gcomp_bias, ipw_bias = ipw_bias))
  res_sim <- c(sdr_res[["res"]], drtrans_res[["res"]], tmle_res[["res"]], gcomp_res[["res"]], ipw_res)
  ## uncommment to add gcomp.glm and br:
  # res_sim <- c(sdr_res[["res"]], drtrans_res[["res"]], tmle_res[["res"]], br_res[["res"]], gcomp_res[["res"]], gcompGLM_res[["res"]], ipw_res)
  print("sim res: "); print(res_sim)

  fW_n <- list(
    sdr_fW_n = sdr_res[["fW_n"]],
    drtran_fW_n = drtrans_res[["fW_n"]],
    tmle_fW_n = tmle_res[["fW_n"]],
    # br_fW_n = br_res[["fW_n"]],
    gcomp_fW_n = gcomp_res[["fW_n"]]
    # gcompGLM_fW_n = gcompGLM_res[["fW_n"]]
    )

  # return(res_sim)
  return(list(res_sim = res_sim, fW_n = fW_n))
}

Qforms <- list(Qc = Qforms.c, Qi = Qforms.i)
gforms <- list(gc = gform.c, gi = gform.i)
Qinteract <- list(Qc = Qinteract.c, Qi = Qinteract.i)
sim_scens <-
    list(Q = c("Qc", "Qi"), g = c("gc", "gi")) %>%
    purrr::cross_d()

sim_scens <- sim_scens %>%
                # rbind(sim_scens,
                #    tibble::tibble(Q = "QSDR", g = "gSDR")) %>%
  dplyr::mutate(Qforms = purrr::map(Q, ~ Qforms[[.x]])) %>%
  dplyr::mutate(gforms = purrr::map(g, ~ gforms[[.x]])) %>%
  dplyr::mutate(Qinteract = purrr::map(Q, ~ Qinteract[[.x]])) %>%
  dplyr::mutate(tvals = tvals) %>%
  tidyr::unite(scen, Q, g, sep = ".", remove = FALSE)

tmp <- foreach::foreach(i = 1:nsims) %dopar% {
  cat("\nsim: ", i, "\n")
  print("sim: " %+% i)
  Dset <- dag_obj
  dat <- data.table(sim(Dset, n = nsamp, wide = FALSE))

  res <- sim_scens %>%
    dplyr::mutate(res = purrr::pmap(sim_scens, sdr_tmle_gcomp_1sim,
                                    dat = dat,
                                    truth = truth_EY_WA))
  res <- res %>%
    dplyr::mutate(sim_res = purrr::map(res, "res_sim")) %>%
    dplyr::mutate(sim_res = purrr::map2(sim_res, scen,
                    function(.x, .y) {names(.x) <- paste0(paste(.y, ".", sep=""), names(.x)); .x})
                  ) %>%
    dplyr::mutate(fW_res = purrr::map(res, "fW_n")) %>%
    dplyr::mutate(fW_res = purrr::map2(fW_res, scen,
                    function(.x, .y) {names(.x) <- paste0(paste(.y, ".", sep=""), names(.x)); .x})
                  ) %>%
    dplyr::mutate(fW_res = purrr::map(fW_res, ~ as.data.table(.x)))
    # res %>% dplyr::select(fW_res) %>% tidyr::unnest() %>% data.table

  fW_DT = do.call(cbind, res[["fW_res"]])
  fW_DT[, ("ID") := 1:nrow(fW_DT)]
  list(res = unlist(res[["sim_res"]]), fW = fW_DT)
}

point_res <- lapply(tmp, '[[', "res")
out = do.call(cbind, point_res)
out_cover <- out[grep("cover", rownames(out)),]
out_MSE <- out[grep("mse", rownames(out)),]
out_bias <- out[grep("bias", rownames(out)),]
out_psi_n <- out[grep("psi_n", rownames(out)),]
out_CIlen <- out[grep("CIlen", rownames(out)),]

results_cov = data.frame(round(abs(rowMeans(out_cover)),3))
colnames(results_cov)=c('Coverage prob.')

results_CIlen = data.frame(round(abs(rowMeans(out_CIlen)),3))
colnames(results_CIlen)=c('Mean CI length')

results_MSE_mean = data.frame(round(abs(rowMeans(out_MSE)),6))
colnames(results_MSE_mean)=c('Mean MSE')
results_MSE = data.frame(round(abs(rowMeans(out_MSE))/min(abs(rowMeans(out_MSE))),2))
colnames(results_MSE)=c('Relative MSE')

results_bias_mean = data.frame(round(abs(rowMeans(out_bias)),6))
colnames(results_bias_mean)=c('Mean Bias')
results_bias = data.frame(round(abs(rowMeans(out_bias))/min(abs(rowMeans(out_bias))),2))
colnames(results_bias)=c('Relative Bias')

results_var = data.frame(apply(out_psi_n, 1, var))
colnames(results_var)=c('Variance')
results_rel_var = results_var / min(results_var)
colnames(results_rel_var)=c('Relative Variance')

fW_res <- data.table::rbindlist(lapply(tmp, '[[', "fW"))
data.table::setkeyv(fW_res, "ID")
fW_var <- fW_res[, lapply(.SD, stats::var), by = ID]
# fW_var <- fW_res[, lapply(.SD, var), by = ID]
fW_var[, "ID" := NULL]
mean_var_fW <- data.frame(round(unlist(fW_var[, lapply(.SD, mean)]),6))
colnames(mean_var_fW) <-  "Mean Variance fW"

res_name <- "results for E[Y_d(t)] for \bar{A}=1, t=" %+% tvals %+% ", N = " %+% nsamp %+% ", nsims = " %+% nsims
print(res_name)
print(results_bias_mean)
print(results_var)
print(results_MSE_mean)
print(results_bias)
print(results_rel_var)
print(results_cov)
print(results_CIlen)
print(results_MSE)
print(mean_var_fW)

results <- list(res_name = res_name, nsamp = nsamp, nsims = nsims,
                results_cov = results_cov, results_CIlen = results_CIlen,
                results_bias_mean = results_bias_mean,
                results_bias = results_bias,
                results_var = results_var,
                results_rel_var = results_rel_var,
                results_MSE_mean = results_MSE_mean,
                results_MSE = results_MSE,
                mean_var_fW = mean_var_fW)
save(results, file = "sim1_results_t" %+% tvals %+% ".Rd")


