## ----------------------------------------------------------------------------
## simcausal DAG
## ----------------------------------------------------------------------------
## Standard of care is A(t)=1
## Experimental treatment is A(t)=0
## We want to evaluate mean survival for staying on standard of care for entire follow-up (\bar{A}(t)=1)
## Once people switch to experimental treatment A(t-1)=0, they stay on it till end of FUP A(t)=0
## * Confounding *:
## * Subjects have a high probability of switching to A to0 (P(A(t)=0)=0.75) only when the probability of adverse outcome (Y=1) is >= 0.8;
## * Otherwise the probability of going on experimental treatment is much lower (0.2)

simDAG_t4 = function(rndseed) {
  require("simcausal")
  D <- DAG.empty() +
      node("UL",    t = 0, distr = "rnorm") +
      node("W",     t = 0, distr = "rnorm", mean = 0) +
      node("L",     t = 0, distr = "rconst", const = abs(UL[t])) +
      node("phiRsk", t = 0, distr = "rconst", const = L[0]) +
      node("hiRsk", t = 0, distr = "rconst", const = plogis(L[0]) > 0.8) +
      node("A",     t = 0, distr = "rbern", prob = plogis(L[0])) +
      node("Y",     t = 0, distr = "rconst", const = 0) +

      node("UL",t = 1, distr = "rnorm") +
      node("L", t = 1, distr = "rconst", const = abs(UL[t])) +
      node("phiRsk", t = 1, distr = "rconst", const = -2 + 0.5*L[t-1] + 0.5*2*L[t]) +
      node("hiRsk", t = 1, distr = "rconst", const = plogis(-2 + 0.5*L[t-1] + 0.5*2*L[t]) > 0.9) +
      node("A", t = 1, distr = "rbern", prob = A[t-1]*plogis(1.7 - 2.0*hiRsk[t])) +
      node("Y", t = 1, distr = "rbern", prob = plogis(-3 + 0.5*L[t-1]*A[t] + 0.5*A[t-1]*L[t] + 0.5*L[t]*A[t])) +

      node("UL",t = 2, distr = "rnorm", mean = A[0]*L[1] + L[1]*A[1]) +
      node("L", t = 2, distr = "rconst", const = abs(UL[t])) +
      node("phiRsk", t = 2, distr = "rconst", const = -2 + 0.5*L[t-1] + 0.5*2*L[t]) +
      node("hiRsk", t = 2, distr = "rconst", const = plogis(-2 + 0.5*L[t-1] + 0.5*2*L[t]) > 0.85) +
      node("A", t = 2, distr = "rbern", prob = A[t-1]*plogis(1.7 - 2.0*hiRsk[t])) +
      node("Y", t = 2, distr = "rbern", prob = plogis(-3*Y[t-1] + 0.5*L[t-1]*A[t] + 0.5*A[t-1]*L[t] + 0.5*L[t]*A[t])) +

      node("UL",t = 3, distr = "rnorm", mean = L[1] * A[t-1] + A[1]*L[t-1] + L[t-1]*A[t-1]) +
      node("L", t = 3, distr = "rconst", const = abs(UL[t])) +
      node("phiRsk", t = 3, distr = "rconst", const = -2 + 0.5*L[t-1] + 0.5*2*L[t]) +
      node("hiRsk", t = 3, distr = "rconst", const = plogis(-2 + 0.5*L[t-1] + 0.5*2*L[t]) > 0.80) +
      node("A", t = 3, distr = "rbern", prob = A[t-1]*plogis(1.7 - 2.0*hiRsk[t])) +
      node("Y", t = 3, distr = "rbern", prob = plogis(-3*Y[t-1] + 0.5*L[t-1]*A[t] + 0.5*A[t-1]*L[t] + 0.5*L[t]*A[t])) +

      node("UL",t = 4, distr = "rnorm", mean = L[1] * A[t-1] + A[1]*L[t-1] + L[t-1]*A[t-1]) +
      node("L", t = 4, distr = "rconst", const = abs(UL[t])) +
      # 0.5*W[0] +
      # 0.1*W[0]*L[t]*(L[t]<3.4)
      # *(L[t-1]>=3.4)
      # (-0.3*W[0]*(W[0] < 0) + 0.9*W[0]*(W[0] >= 0))
      node("phiRsk", t = 4, distr = "rconst", const = -1 + 0.25*L[t-1] + 0.25*2*L[t] - 0.1*L[t]*L[t-1] + 1.5*W[0]*L[t-1]) +
      # node("phiRsk", t = 4, distr = "rconst", const = -1 + 0.5*L[t-1] + 0.5*2*L[t] + 1.5*W[0]) +
      node("hiRsk", t = 4, distr = "rconst", const = plogis(phiRsk[t]) > 0.80) +
      # node("hiRsk", t = 4, distr = "rconst", const = plogis(-1 + 0.5*L[t-1] + 0.5*2*L[t] +  1.5*W[0]) > 0.80) +
      node("A", t = 4, distr = "rbern", prob = A[t-1]*plogis(2 - 2.0*hiRsk[t])) +
      node("Y", t = 4, distr = "rbern", prob = plogis(-1*Y[t-1] + A[t] + phiRsk[t]*A[t] + 0.20*A[t-1]*L[t]) )
      # node("Y", t = 4, distr = "rbern", prob = plogis(-1*Y[t-1] + 0.5*L[t-1]*A[t] + 0.5*A[t-1]*L[t] + 0.5*L[t]*A[t] + 1.5*W[0]*A[t]))

  Dset <-
      set.DAG(D) +
      action("noIntervention", nodes = node("L", t = 0, distr = "rnorm")) +
      action("allAare1", nodes =
          c(node("A", t = c(0:4), distr = "rbern", prob = 1),
              node("Y", t = 4, distr = "rconst", const = plogis(-1*Y[t-1] + A[t] + phiRsk[t]*A[t] + 0.20*A[t-1]*L[t]) ))
              # node("Y", t = 4, distr = "rconst", const = plogis(-1*Y[t-1] + 0.5*L[t-1]*A[t] + 0.5*A[t-1]*L[t] + 0.5*L[t]*A[t] + 1.5*W[0]*A[t]) ))
          )

  # dat <- data.table(sim(Dset, n = 10000, wide = FALSE))
  # dat[, ("Atm1") := shift(.SD, n=1L, fill=1L, type="lag"), by=ID, .SDcols="A"]
  # # dat[Atm1==1, A]
  # dat[, all(A==0), by = ID][,sum(V1)]
  # dat[, all(A==1), by = ID][,sum(V1)]
  # dat[, sum(A==0, na.rm = TRUE), by = t]
  # dat[, sum(A==1, na.rm = TRUE), by = t]

  # X2 <- sim(DAG = Dset, actions = "noIntervention", n = 10000, rndseed = 12345)
  # # # # # # X2 <- sim(DAG = Dset, actions = "noIntervention", n = 1000000, rndseed = 12345)
  # D2 <- set.targetE(Dset, outcome="Y", t=1:4, param="noIntervention")
  # psi0_byt <- eval.target(DAG = D2, data = X2)$res
  # print("observed event prob: "); print(psi0_byt)

  # X2 <- sim(DAG = Dset, actions = "allAare1", n = 10000, rndseed = 12345, wide = FALSE)
  # dat <- data.table(X2[[1]])
  # tval <- 4
  # plot(density(dat[t==tval, L]))
  # round(dat[t==tval, L], 2)
  # summary(round(dat[t==tval, L], 2))
  # round(dat[t==4, L], 2)*round(dat[t==3, L], 2)
  # summary(round(dat[t==4, L], 2)*round(dat[t==3, L], 2))
  # plot(density(round(dat[t==2, L], 2)*round(dat[t==1, L], 2)))

  # summary(dat[t==tval, L])
  # round(dat[t==tval, ][["phiRsk"]], 3)
  # summary(dat[t==tval, ][["phiRsk"]])
  # round(dat[t==tval, ][["phiRsk"]],4)
  # round(plogis(dat[t==tval, ][["phiRsk"]]),4)
  # plot(density(dat[t==tval, ][["phiRsk"]]))
  # plot(density(plogis(dat[t==tval, ][["phiRsk"]])))
  # hist(plogis(dat[t==tval, ][["phiRsk"]]))
  # mean(dat[t==tval & plogis(phiRsk) > 0.8, ][["A"]])
  # mean(dat[t==tval & plogis(phiRsk) > 0.8, ][["Y"]])
  # mean(dat[t==tval & plogis(phiRsk) < 0.6, ][["A"]])
  # mean(dat[t==tval & plogis(phiRsk) < 0.6, ][["Y"]])
  # mean(dat[t==tval, ][["Y"]])

  return(Dset)
}
