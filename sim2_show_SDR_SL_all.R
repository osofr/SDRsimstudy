`%+%` <- function(a, b) paste0(a, b)
require("data.table")
require("simcausal")
require("foreach")
require("doParallel")
require("magrittr")
options(simcausal.verbose = FALSE)
# nsims <- 2
nsims <- 1000
# nsamp <- 1000
# nsamp <- 2000
nsamp <- 5000
# nsamp <- 10000
tvals <- 4

registerDoParallel(cores=detectCores())

# source("data.gen_v1.R")
# source("data.gen_v3.R")
source("data.gen_v4.R")

## g model (same across all time-points):
gform.c <- c(rep.int("A ~ hiRsk + L + Atm1 + Ytm1", tvals), "A ~ hiRsk + Atm1 + Ytm1")
gform.i <- "A ~ 1"

## Q models:
Qforms.c <- c("Qkplus1 ~ L + W",
              "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
              "Qkplus1 ~ L + Ltm1 + Ltm2 + A + Atm1 + Atm2 + Ytm1 + W",
              "Qkplus1 ~ L + Ltm1 + Ltm2 + A + Atm1 + Atm2 + Ytm1 + W",
              "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W"
              )
Qinteract.c <- NULL
## v1: SDR-methods do not save bias/MSE under this misspecification (almost as bad as TMLE):
Qforms.i <- c("Qkplus1 ~ L + W", "Qkplus1 ~ W", "Qkplus1 ~ W", "Qkplus1 ~ W", "Qkplus1 ~ W")
## v2: Even though Q is wrong, SDR methods still give a big reduction in bias/MSE:
# Qforms.i <- c("Qkplus1 ~ L + W", "Qkplus1 ~ Ltm1", "Qkplus1 ~ Ltm1", "Qkplus1 ~ Ltm1", "Qkplus1 ~ Ltm1")

Qinteract.i <- NULL

## g and Q model to demonstrate SDR property (only last g is correctly specified)
gform.SDR <- c(rep.int("A ~ 1", tvals-1), "A ~ 1", "A ~ hiRsk + Atm1 + Ytm1")
## v1: TMLE is now biased for marginal parameter, but reduction in MSE is not as crazy for SDR (as with v2):
Qforms.SDR <- c("Qkplus1 ~ L + W",
                "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
                "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
                "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
                "Qkplus1 ~ W")
## v2: TMLE is still unbiased for marginal parameter, but crazy reduction in MSE for SDR:
# Qforms.SDR <- c("Qkplus1 ~ L + W",
#                 "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
#                 "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
#                 "Qkplus1 ~ L + Ltm1 + A + Atm1 + Ytm1 + W",
#                 "Qkplus1 ~ Ltm1")
Qinteract.SDR <- NULL

eval_true_EY_WA <- function(dag_obj, action_name = "allAare1", tval = 4, ntest = 1e4, mc.repl = 1000, rndseed = 12345, rndseed.reset.node = "Y_0") {
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  require("doParallel")
  require("foreach")
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  ex <- c(ex, "nfull", "simDAG_t1", "simDAG_t2")
  Dset <- dag_obj
  dat0_repl <- foreach(t.counter = (1 : mc.repl), .packages = c("simcausal"), .export = ex, .combine = "+") %dopar% {
    dat <- data.table::data.table(sim(Dset, actions = action_name, LTCF = "Y", n = ntest, rndseed = rndseed, rndseed.reset.node = rndseed.reset.node)[[1]])
    dat[[c("Y_" %+% tval)]]
  }
  P0_EY_WA <- dat0_repl / mc.repl
  dat0 <- data.table::data.table(sim(Dset, actions = action_name, n = ntest, rndseed = rndseed)[[1]])
  test.dat <- dat0[, c("W_0", "L_0","A_0"), with = FALSE]
  setnames(test.dat, c("W", "L", "A"))
  test.dat[, ("ID") := seq.int(nrow(test.dat))]
  test.dat[, ("t") := 0]
  return(list(test.dat = test.dat, P0_EY_WA = P0_EY_WA))
}

dag_obj <- do.call("simDAG_t" %+% tvals, args = list()) # dag_obj <- simDAG_t1(); dag_obj <- simDAG_t2()

## option 1 to evaluate psi0:
# truth_psi0 <- mean(sim(dag_obj, actions = "allAare1", LTCF = "Y", n = 1000000, rndseed = 12345)[[1]][["Y_"%+%tvals]])
## option 2 to evaluate psi0:
# X2 <- sim(DAG = dag_obj, actions = "allAare1", n = 1000000, rndseed = 12345)
X2 <- sim(DAG = dag_obj, actions = "allAare1", n = 2000000, rndseed = 12345)
D2 <- set.targetE(dag_obj, outcome="Y", t=tvals, param="allAare1")
truth_psi0 <- eval.target(DAG = D2, data = X2)$res
print("truth: "); print(truth_psi0)

# truth_EY_WA <- eval_true_EY_WA(dag_obj, ntest = 1e4, tval = tvals, mc.repl = 200, rndseed = 12345)
truth_EY_WA <- eval_true_EY_WA(dag_obj, ntest = 1e4, tval = tvals, mc.repl = 20000, rndseed = 12345)
truth_EY_WA[["P0_EY_WA"]]
print(var(truth_EY_WA[["P0_EY_WA"]]))
truth_EY_WA[["truth_psi0"]] <- truth_psi0
print(truth_EY_WA[["test.dat"]])
save(truth_EY_WA, file = "./truth_EY_WA_t" %+% tvals %+% ".Rd")
# load("./truth_EY_WA_t" %+% tvals %+% ".Rd")

sdr_tmle_gcomp_1sim <- function(dat, truth, tvals = 5, Qforms, gforms, Qinteract, ...) {
  library("data.table")
  setDTthreads(1)
  library("xgboost")
  # library("gridisl")
  library("stremr")
  options(stremr.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  options(gridisl.verbose = FALSE)
  # options(gridisl.verbose = TRUE)

  tvector <- 0:tvals
  truth_psi0 <- truth[["truth_psi0"]]
  test.dat <- truth[["test.dat"]]
  test.dat <- data.table(test.dat)
  P0_EY_WA <- truth[["P0_EY_WA"]]

  ID <- "ID"
  t <- "t"
  TRT <- "A"
  outcome <- "Y"
  covars <- c("L", "Ltm1", "Ltm2", "Atm1", "Atm2", "Ytm1")
  # covars <- c("L", "Ltm1", "Atm1", "A0L1", "L0A0", "Ytm1")

  # Define counterfactual (intervention) node
  # dat[, ("A.star") := A][t == 1, ("A.star") := 1]
  dat[, ("A.star") := 1]

  # Define lag nodes and interactions.
  setkeyv(dat, cols = c("ID", "t"))
  lag1nodes <- c("L", "A", "Y")
  lag2tnodes <- c("L", "A")
  newlag1names <- paste0(lag1nodes, "tm1")
  newlag2names <- paste0(lag2tnodes, "tm2")
  dat[, (newlag1names) := shift(.SD, n=1L, fill=1L, type="lag"), by=ID, .SDcols=(lag1nodes)]
  dat[, (newlag2names) := shift(.SD, n=2L, fill=1L, type="lag"), by=ID, .SDcols=(lag2tnodes)]
  # dat[, ("L0A1") := A * Ltm1]
  # dat[, ("A0L1") := Atm1 * L]
  # dat[, ("L1A1") := A * L]
  # dat[, ("L0A0") := Atm1 * Ltm1]
  # Dummies for censoring and monitoring
  dat[, ("C") := 0][, ("N") := 1]

  gform_CENS <- "C ~ 1"
  gform_MONITOR <- "N ~ 1"

  # IMPORT DATA
  OData <- stremr::importData(dat, ID = ID, t = t, covars = covars, CENS = "C", TRT = TRT, MONITOR = "N", OUTCOME = "Y", remove_extra_rows = FALSE)
  OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)

  # Time-specific Propensity score models for Treatment (A):
  params_g <- gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial")
  # params_g <- gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial", interactions = list(c("Atm1","L")))
  gform_TRT <- gforms
  if (length(gform_TRT)==1L) gform_TRT <- rep.int(gform_TRT, length(tvector))

  reg_TRT <- unlist(lapply(seq_along(tvector), function(t_idx) {
    define_single_regression(OData, gform_TRT[t_idx], stratify = list(A = "(t == " %+% tvector[t_idx] %+% ") & (Atm1 == 1L)"), models = params_g)
    }))
  reg_TRT <- c(
    reg_TRT,
    define_single_regression(OData, "A ~ Atm1", stratify = list(A = "(t > 0L) & (Atm1 == 0L)"), models = params_g))

  # FIT PROPENSITY SCORES WITH speedglm
  OData <- fitPropensity(OData, gform_CENS = gform_CENS,
                          reg_TRT = reg_TRT,
                          gform_MONITOR = gform_MONITOR,
                          estimator = "speedglm__glm",
                          family = "quasibinomial", # estimator = "xgboost__gbm",
                          fit_method = "none"
                          # gform_TRT = gform_TRT,
                          # stratify_TRT = stratify_TRT,
                          # fit_method = "cv", # fold_column = "fold_ID",
                          )

  analysis <- list(intervened_TRT = c("A.star"),
                   stratifyQ_by_rule = TRUE) %>%
                purrr::cross_d() %>%
                dplyr::arrange(stratifyQ_by_rule)

  # params <- gridisl::defModel(estimator = "speedglm__glm", family = "gaussian")
  params <-
    gridisl::defModel(estimator = "xgboost__gbm", family = "quasibinomial",
                      # nrounds = 250,
                      # nrounds = 350,
                      # nrounds = 200,
                      nthread = 1,
                      colsample_bytree = 1,
                      # gamma = 0, # gamma = 1, # gamma = 3,
                      # min_child_weight = 20, # min_child_weight = 8, # min_child_weight = 3,
                      # max_delta_step = 5,
                      subsample = .5,
                      # alpha = .75,
                      # lambda = 2,
                      remove_const_cols = TRUE,
                      interactions = Qinteract,
                      param_grid = list(
                        nrounds = c(50, 200, 250),
                        # lambda = c(0,2),
                        # gamma = c(0,1),
                        min_child_weight = c(5, 10, 15), # min_child_weight = c(15, 18, 20), # min_child_weight = c(20,30),
                        max_depth = c(2,6), # max_depth = c(2,6,10,15), # max_depth = c(6,8,10), # max_depth = c(3),
                        learning_rate = c(.05) # learning_rate = c(.1, .2),
                        )
                      ) +
    gridisl::defModel(estimator = "xgboost__glm",
                      family = "quasibinomial",
                      nthread = 1,
                      nrounds = 200,
                      max_delta_step = 2,
                      interactions = Qinteract
                      )

  paramsGLM <-
      gridisl::defModel(estimator = "xgboost__glm",
                      family = "quasibinomial",
                      nthread = 1,
                      nrounds = 200,
                      max_delta_step = 2,
                      interactions = NULL
                      )

  # reg_Q <- list(paramsGLM, paramsGLM, paramsGLM, paramsGLM, params)
  reg_Q <- list(paramsGLM, paramsGLM, paramsGLM, params, params)

  print(OData$modelfit.gA$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$get.fits()[[1]]$get_best_models()[[1]])
  print(OData$modelfit.gA$getPsAsW.models()[[2]]$getPsAsW.models()[[1]]$get.fits()[[1]]$get_best_models()[[1]])
  print(OData$modelfit.gA$getPsAsW.models()[[3]]$getPsAsW.models()[[1]]$get.fits()[[1]]$get_best_models()[[1]])
  print(OData$modelfit.gA$getPsAsW.models()[[4]]$getPsAsW.models()[[1]]$get.fits()[[1]]$get_best_models()[[1]])
  print(OData$modelfit.gA$getPsAsW.models()[[5]]$getPsAsW.models()[[1]]$get.fits()[[1]]$get_best_models()[[1]])

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

  ## hazard IPW (wrong for none time-to-event outcome) w/ GLM g
  # hazardIPW <-  analysis %>%
  #         dplyr::distinct(intervened_TRT) %>%
  #         dplyr::mutate(wts_data = purrr::map(intervened_TRT, getIPWeights, OData = OData)) %>%
  #         dplyr::mutate(NPMSM = purrr::map(wts_data, ~ survNPMSM(wts_data = .x, OData = OData))) %>%
  #         dplyr::mutate(NPMSM = purrr::map(NPMSM, "estimates"))
  # St.hazardIPW <- hazardIPW %>% dplyr::select(NPMSM) %>% tidyr::unnest()
  # hazardipw_bias <- ((1 - data.table::as.data.table(St.hazardIPW)[time == tvals, ][["St.NPMSM"]]) - truth_psi0)

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

  ## SDR w/ SL Q and split-specific CV for targeting (using xgboost)
  SDR_est <- stremr:::fitSDR(OData, tvals = tvals,
                      intervened_TRT = "A.star", Qforms = Qforms,
                      stratifyQ_by_rule = TRUE,
                      fit_method = "cv",
                      models = params,
                      return_fW = TRUE,
                      reg_Q = reg_Q)
  sdr_res <- eval_est_perf(SDR_est, est = "sdr")
  print("sdr_res"); print(sdr_res[["res"]]); cat("SDR CI: ", sdr_res[["est_CI"]])

  ## DR-transform w/ SL Q
  DR_trans_est <- stremr:::fitSDR(OData, tvals = tvals,
                           intervened_TRT = "A.star", Qforms = Qforms,
                           stratifyQ_by_rule = TRUE,
                           fit_method = "cv",
                           models = params,
                           return_fW = TRUE,
                           use_DR_transform = TRUE, # stabilize = FALSE,
                           reg_Q = reg_Q
                          )
  drtrans_res <- eval_est_perf(DR_trans_est, est = "drtrans")
  print("drtrans_res"); print(drtrans_res[["res"]]); cat("DR Transform CI: ", drtrans_res[["est_CI"]])

  ## TMLE w/ SL Q
  tmle_est <- fitTMLE(OData, tvals = tvals,
                      intervened_TRT = "A.star", Qforms = Qforms,
                      stratifyQ_by_rule = TRUE,
                      fit_method = "cv",
                      models = params,
                      return_fW = TRUE,
                      reg_Q = reg_Q)
  tmle_res <- eval_est_perf(tmle_est, est = "tmle")
  print("tmle_res"); print(tmle_res[["res"]]); cat("TMLE CI: ", tmle_res[["est_CI"]])

  ## GCOMP w/ SL Q
  gcomp_est <- fitSeqGcomp(OData, tvals = tvals,
                           intervened_TRT = "A.star", Qforms = Qforms,
                           stratifyQ_by_rule = TRUE,
                           fit_method = "cv",
                           models = params,
                           return_fW = TRUE,
                           reg_Q = reg_Q)
  gcomp_res <- eval_est_perf(gcomp_est, est = "gcomp")
  print("gcomp_res"); print(gcomp_res[["res"]])
  # gcompSL_fW_n <- gcomp_res[["fW_n"]]
  # test.dat[, ("fWgcompSL") := gcompSL_fW_n]
  # gcomp_res <- NULL

  ## BR (Band & Robins DR estimator, TMLE w/ GLM Q, main terms logistic regression)
  # br_est <- fitTMLE(OData, tvals = tvals,
  #                   intervened_TRT = "A.star", Qforms = Qforms,
  #                   stratifyQ_by_rule = TRUE,
  #                   estimator = "speedglm__glm",
  #                   fit_method = "none",
  #                   return_fW = TRUE)
  # br_res <- eval_est_perf(br_est, est = "br")
  # print("br_res"); print(br_res[["res"]])

  ## GCOMP w/ GLM Q
  # gcompGLM_est <- fitSeqGcomp(OData, tvals = tvals,
  #                          intervened_TRT = "A.star", Qforms = Qforms,
  #                          stratifyQ_by_rule = TRUE,
  #                          estimator = "speedglm__glm",
  #                          fit_method = "none",
  #                          return_fW = TRUE)
  # gcompGLM_res <- eval_est_perf(gcompGLM_est, est = "gcompGLM")
  # print("gcompGLM_res"); print(gcompGLM_res[["res"]])

  # gcompGLM_fW_n <- gcompGLM_res[["fW_n"]]
  # test.dat[, ("fWgcompGLM") := gcompGLM_fW_n]
  # test.dat[, ("Ytrue") := P0_EY_WA]
  # round(test.dat[, fWgcompGLM - Y1], 3)
  # mean(round(test.dat[, fWgcompGLM - Y1], 3))
  # mean(test.dat[, Y1])
  # mean(test.dat[, fWgcompGLM])
  # mean(test.dat[, fWgcompSL])
  # test.dat[1:100, ]
  # var(test.dat[, fWgcompSL])
  # mean(abs(test.dat[, fWgcompSL]-test.dat[, Ytrue]))
  # mean(abs(test.dat[, fWgcompGLM]-test.dat[, Ytrue]))

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

## ------------------------------------------
## Define simulation scenarios:
## ------------------------------------------
Qforms <- list(Qc = Qforms.c, Qi = Qforms.i, QSDR = Qforms.SDR)
gforms <- list(gc = gform.c, gi = gform.i, gSDR = gform.SDR)
Qinteract <- list(Qc = Qinteract.c, Qi = Qinteract.i, QSDR = Qinteract.SDR)

sim_scens <-
    list(Q = c("Qc", "Qi"), g = c("gc", "gi")) %>%
    # list(Q = c("Qc"), g = c("gc")) %>%
    purrr::cross_d()

sim_scens <- rbind(sim_scens,
                   tibble::tibble(Q = "QSDR", g = "gSDR")) %>%
  dplyr::mutate(Qforms = purrr::map(Q, ~ Qforms[[.x]])) %>%
  dplyr::mutate(gforms = purrr::map(g, ~ gforms[[.x]])) %>%
  dplyr::mutate(Qinteract = purrr::map(Q, ~ Qinteract[[.x]])) %>%
  dplyr::mutate(tvals = tvals) %>%
  tidyr::unite(scen, Q, g, sep = ".", remove = FALSE)

# sim_scens <- sim_scens[5, ]

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
save(results, file = "sim2_results_t" %+% tvals %+% ".Rd")


