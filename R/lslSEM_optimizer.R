


.sem_estep <- function(Cyc, ey, Sg, mu, Ld, Ps, Bt, Ph, nu, ap) {
  Sgiv <- solve(Sg)
  IBtiv <- solve(diag(1, dim(Ld)[2]) - Bt)
  Xi <- IBtiv %*% Ph %*% t(IBtiv)
  KM <- Xi %*% t(Ld) %*% Sgiv
  JM <- IBtiv %*% ap - KM %*% mu
  eet <- JM + KM %*% ey
  Cyet <- ey %*% t(JM) + (Cyc + tcrossprod(ey)) %*% t(KM)
  Cet <-
    Xi - Xi %*% t(Ld) %*% t(KM) + JM %*% t(JM) + JM %*% t(ey) %*% t(KM) + KM %*% ey %*% t(JM) +
    KM %*% (Cyc + tcrossprod(ey)) %*% t(KM)
  mis_moment <- list(Cet = Cet, Cyet = Cyet, eet = eet)
  return(mis_moment)
}



.sem_mstep_Ld <-
  function(Ldp,
           Ld,
           Ps,
           nu,
           Cyet,
           Cet,
           ey,
           eet,
           type,
           gm,
           dt,
           P,
           Psp,
           Psp_type) {
    if (all(!.is_est(Ldp))) {
      Ld <- Ld
    } else {
      if (Psp_type == 1) {
        Psiv <- matrix(0, P, P)
      } else if (Psp_type == 2) {
        Psiv <- diag(1 / diag(Ps))
      } else if (Psp_type == 3) {
        Psiv <- solve(Ps)
      } else {
        merror_idc = .is_est(diag(Psp))
        if (all(merror_idc)) {
          Psiv <- solve(Ps)
        } else {
          Psiv <- matrix(0, P, P)
          Psiv[merror_idc, merror_idc] <-
            solve(Ps[merror_idc, merror_idc])
        }
      }
      for (p in 1:P) {
        idx <- which(.is_est(Ldp[p,]))
        for (j in idx) {
          CRsum <- 0
          if (Psp_type != 1) {
            for (k in (1:P)[(1:P) != p]) {
              CRsum <- CRsum + (Psiv[p, k] / Psiv[p, p]) *
                (Cyet[k, j] - nu[k, 1] * eet[j, 1] - Ld[k, , drop = FALSE] %*% Cet[, j, drop = FALSE])
            }
          }
          ldq <-
            ((Cyet[p, j] - nu[p, 1] * eet[j, 1] - Ld[p, -j, drop = FALSE] %*% Cet[-j, j, drop = FALSE]) + CRsum) / Cet[j, j]
          if (.is_one(Ldp[p, j])) {
            Ld[p, j] <- ldq
          } else {
            wq <- 1 / (Psiv[p, p] * Cet[j, j])
            Ld[p, j] <- .sem_threshold(ldq, type, wq, gm, dt)
          }
        }
      }
    }
    return(Ld)
  }


.sem_mstep_Ps <-
  function(Psp,
           Ps,
           Ld,
           nu,
           Cyc,
           Cyet,
           Cet,
           ey,
           eet,
           type,
           gm,
           dt,
           P,
           Psp_type) {
    Cep <-
      Cyc + tcrossprod(ey) - nu %*% t(ey) - Ld %*% t(Cyet) - ey %*% t(nu) + nu %*% t(nu) +
      Ld %*% eet %*% t(nu) - Cyet %*% t(Ld) + nu %*% t(eet) %*% t(Ld) + Ld %*% Cet %*% t(Ld)
    if (Psp_type == 1) {
      Ps <- Ps
    } else if (Psp_type == 2) {
      Ps <- diag(diag(Cep))
    } else if (Psp_type == 3) {
      Ps <- Cep
    } else {
      merror_idc <- .is_est(diag(Psp))
      for (p in 1:P) {
        if (merror_idc[p]) {
          idx <- which(.is_est(Psp[, p]))
          yidxp <- which((.is_est(diag(Psp[-p,-p, drop = FALSE]))))
          Psivp <- matrix(0, P - 1, P - 1)
          Psivp[yidxp, yidxp] <- solve((Ps[-p,-p, drop = FALSE])[yidxp, yidxp, drop = FALSE])
          for (j in idx[idx > p]) {
            U <- Psivp %*% Cep[-p, p, drop = FALSE]
            V <- Psivp %*% Cep[-p, -p, drop = FALSE] %*% Psivp
            Vsum = 0
            for (k in idx) {
              if (k < p) {
                Vk <- V[(j - 1), k]
              }
              else if (k == p) {
                Vk <- 0
              }
              else if (k == j) {
                Vk <- 0
              }
              else {
                Vk <- V[(j - 1), (k - 1)]
              }
              Vsum <- Vsum + Ps[k, p] * Vk
            }
            psq <- (U[(j - 1)] - Vsum) / V[(j - 1), (j - 1)]
            vpsp <-
              Ps[p, p] - Ps[p, -p, drop = FALSE] %*% Psivp %*% Ps[-p, p, drop = FALSE]
            if (.is_one(Psp[p, j])) {
              Ps[p, j] <- psq
              Ps[j, p] <- Ps[p, j]
            } else if (is.na(Psp[p, j])) {
              wq <- 1 / ((1 / vpsp) * V[(j - 1), (j - 1)])
              Ps[p, j] <- .sem_threshold(psq, type, wq, gm, dt)
              Ps[j, p] <- Ps[p, j]
            } else {
              
            }
          }
          U <- Psivp %*% matrix(Cep[-p, p, drop = FALSE], (P - 1), 1)
          V <- Psivp %*% Cep[-p, -p, drop = FALSE] %*% Psivp
          Vsum <- 0
          Usum <- 0
          for (k in idx) {
            if (k < p) {
              Uk <- U[k]
            }
            else if (k == p) {
              Uk <- 0
            }
            else {
              Uk <- U[(k - 1)]
            }
            Usum = Usum + Ps[k, p] * Uk
            for (l in idx) {
              if (l < p & k < p) {
                Vkl <- V[l, k]
              }
              else if (k == p | l == p) {
                Vkl <- 0
              }
              else if (l < p & k > p) {
                Vkl <- V[l, (k - 1)]
              }
              else if (l > p & k < p) {
                Vkl <- V[(l - 1), k]
              }
              else {
                Vkl <- V[(l - 1), (k - 1)]
              }
              Vsum <- Vsum + Ps[k, p] * Ps[l, p] * Vkl
            }
          }
          yidxp <- which((.is_est(diag(Psp[-p, -p, drop = FALSE]))))
          Psivp <- matrix(0, P - 1, P - 1)
          Psivp[yidxp, yidxp] <- solve((Ps[-p, -p, drop = FALSE])[yidxp, yidxp, drop = FALSE])
          Ps[p, p] <- Cep[p, p] - 2 * Usum + Vsum + Ps[p, -p, drop = FALSE] %*% Psivp %*% Ps[-p, p, drop = FALSE]
        }
      }
    }
    return(Ps)
  }


.sem_mstep_Bt <-
  function(Btp,
           Bt,
           Ph,
           ap,
           Cet,
           eet,
           type,
           gm,
           dt,
           M,
           Php_type) {
    if (all(!.is_est(Btp))) {
      Bt <- Bt
    } else {
      if (Php_type == 1) {
        Phiv <- matrix(0, M, M)
      } else if (Php_type == 2) {
        Phiv <- diag(1 / diag(Ph))
      } else if (Php_type == 3) {
        Phiv <- solve(Ph)
      } else {
        Phiv <- solve(Ph)
      }
      for (m in 1:M) {
        idx <- which(.is_est(Btp[m,]))
        for (j in idx) {
          CRsum <- 0
          if (Php_type != 1) {
            for (k in (1:M)[(1:M) != m]) {
              CRsum <-
                CRsum + (Phiv[m, k] / Phiv[m, m]) * (Cet[k, j] - ap[k, 1] * eet[j, 1] - Bt[k, , drop = FALSE] %*% Cet[, j, drop = FALSE])
            }
          }
          btq <-
            ((Cet[m, j] - ap[m, 1] * eet[j, 1] - Bt[m, -j] %*% Cet[-j, j, drop = FALSE]) + CRsum) / Cet[j, j]
          if (.is_one(Btp[m, j])) {
            Bt[m, j] <- btq
          } else {
            wq <- 1 / (Phiv[m, m] * Cet[j, j])
            Bt[m, j] <- .sem_threshold(btq, type, wq, gm, dt)
          }
        }
      }
    }
    return(Bt)
  }


.sem_mstep_Ph <-
  function(Php,
           Ph,
           Bt,
           ap,
           Cet,
           eet,
           type,
           gm,
           dt,
           M,
           Php_type) {
    Czt <-
      Cet - ap %*% t(eet) - Bt %*% Cet - eet %*% t(ap) + ap %*% t(ap) + Bt %*% eet %*% t(ap) -
      Cet %*% t(Bt) + ap %*% t(eet) %*% t(Bt) + Bt %*% Cet %*% t(Bt)
    if (Php_type == 1) {
      Ph <- Ph
    } else if (Php_type == 2) {
      if (all(dim(Czt) == c(1, 1))) {
        Ph <- Czt
      } else {
        Ph <- diag(diag(Czt))
      }
    } else if (Php_type == 3) {
      Ph <- Czt
    } else {
      for (m in 1:M) {
        idx <- which(.is_est(Php[, m]))
        for (j in idx[idx > m]) {
          Phivm <- solve(Ph[-m, -m])
          U <- Phivm %*% Czt[-m, m, drop = FALSE]
          V <- Phivm %*% Czt[-m, -m, drop = FALSE] %*% Phivm
          Vsum <- 0
          for (k in idx) {
            if (k < m) {
              Vk <- V[(j - 1), k]
            }
            else if (k == m) {
              Vk <- 0
            }
            else if (k == j) {
              Vk <- 0
            }
            else {
              Vk <- V[(j - 1), (k - 1)]
            }
            Vsum <- Vsum + Ph[k, m] * Vk
          }
          
          vphm <-
            Ph[m, m] - Ph[m, -m, drop = FALSE] %*% Ph[-m, -m, drop = FALSE] %*% Ph[-m, m, drop = FALSE]
          phq <- (U[(j - 1)] - Vsum) / V[(j - 1), (j - 1)]
          if (.is_one(Php[m, j])) {
            Ph[m, j] <- phq
            Ph[j, m] <- Ph[m, j]
          } else {
            wq <- 1 / ((1 / vphm) * V[(j - 1), (j - 1)])
            Ph[m, j] <- .sem_threshold(phq, type, wq, gm, dt)
            Ph[j, m] <- Ph[m, j]
          }
        }
        Phivm <- solve(Ph[-m, -m, drop = FALSE])
        U <- Phivm %*% Czt[-m, m, drop = FALSE]
        V <- Phivm %*% Czt[-m, -m, drop = FALSE] %*% Phivm
        Vsum <- 0
        Usum <- 0
        for (k in idx) {
          if (k < m) {
            Uk <- U[k]
          }
          else if (k == m) {
            Uk <- 0
          }
          else {
            Uk <- U[(k - 1)]
          }
          Usum <- Usum + Ph[k, m] * Uk
          for (l in idx) {
            if (l < m & k < m) {
              Vkl <- V[l, k]
            }
            else if (k == m | l == m) {
              Vkl <- 0
            }
            else if (l < m & k > m) {
              Vkl <- V[l, (k - 1)]
            }
            else if (l > m & k < m) {
              Vkl <- V[(l - 1), k]
            }
            else {
              Vkl <- V[(l - 1), (k - 1)]
            }
            Vsum <- Vsum + Ph[k, m] * Ph[l, m] * Vkl
          }
        }
        Ph[m, m] <- Czt[m, m] - 2 * Usum + Vsum + Ph[m, -m, drop = FALSE] %*% solve(Ph[-m, -m, drop = FALSE]) %*% Ph[-m, m, drop = FALSE]
      }
    }
    return(Ph)
  }



.sem_mstep_nu <-
  function(nup,
           nu,
           Ld,
           Ps,
           ey,
           eet,
           type,
           gm,
           dt,
           P,
           Psp,
           Psp_type) {
    if (all(!.is_est(nup))) {
      nu <- nu
    } else {
      if (all(nup == 1)) {
        nu <- ey
      } else {
        if (Psp_type == 1) {
          Psiv <- matrix(0, P, P)
        } else if (Psp_type == 2) {
          Psiv <- diag(1 / diag(Ps))
        } else if (Psp_type == 3) {
          Psiv <- solve(Ps)
        } else {
          merror_idc = .is_est(diag(Psp))
          if (all(merror_idc)) {
            Psiv <- solve(Ps)
          } else {
            Psiv <- matrix(0, P, P)
            Psiv[merror_idc, merror_idc] <-
              solve(Ps[merror_idc, merror_idc])
          }
          idx <- which(.is_est(nup[, 1]))
          for (p in idx) {
            CRsum <- 0
            for (k in (1:P)[(1:P) != p]) {
              CRsum <-
                CRsum + (Psiv[p, k] / Psiv[p, p]) * (ey[k, 1] - nu[k, 1] - Ld[k, , drop = FALSE] %*% eet)
            }
            nuq <- ((ey[p, 1] - Ld[p,] %*% eet) + CRsum)
            if (.is_one(nup[p, 1])) {
              nu[p, 1] <- nuq
            } else {
              wq <- 1 / (Psiv[p, p])
              nu[p, 1] <- .sem_threshold(nuq, type, wq, gm, dt)
            }
          }
        }
      }
    }
    return(nu)
  }


.sem_mstep_ap = function(app, ap, Bt, Ph, eet, type, gm, dt, M, Php_type) {
  if (all(!.is_est(app))) {
    ap <- ap
  } else {
    if (Php_type == 1) {
      Phiv <- matrix(0, M, M)
    } else if (Php_type == 2) {
      Phiv <- diag(1 / diag(Ph))
    } else if (Php_type == 3) {
      Phiv <- solve(Ph)
    } else {
      Phiv <- solve(Ph)
    }
    idx <- which(.is_est(app[, 1]))
    for (m in idx) {
      CRsum <- 0
      for (k in (1:M)[(1:M) != m]) {
        CRsum <-
          CRsum + (Phiv[m, k] / Phiv[m, m]) * (eet[k, 1] - ap[k, 1] - matrix(Bt[k, , drop = FALSE], 1, M) %*% eet)
      }
      apq <-
        ((eet[m, 1] - matrix(Bt[m, , drop = FALSE], 1, M) %*% eet) + CRsum)
      if (.is_one(app[m, 1])) {
        ap[m, 1] <- apq
      } else {
        wq <- 1 / (Phiv[m, m])
        ap[m, 1] <- .sem_threshold(apq, type, wq, gm, dt)
      }
    }
  }
  return(ap)
}


.sem_ecm <-
  function(pattern,
           value,
           penalty_fit,
           control,
           info,
           obs_moment) {
    N <- info$N
    P <- info$P
    M <- info$M
    Qall <- info$Qall
    
    Ldp <- pattern$lambda
    Psp <- pattern$psi
    Btp <- pattern$beta
    Php <- pattern$phi
    nup <- pattern$nu
    app <- pattern$alpha
    
    if (all(!.is_est(Psp))) {
      Psp_type = 1
    } else if (all(diag(.is_one(Psp))) &
               all(!.is_est(Psp[lower.tri(Psp)]))) {
      Psp_type = 2
    } else if (all(.is_one(Psp))) {
      Psp_type = 3
    } else {
      Psp_type = 4
    }
    
    if (all(!.is_est(Php))) {
      Php_type = 1
    } else if (all(diag(.is_one(Php))) &
               all(!.is_est(Php[lower.tri(Php)]))) {
      Php_type = 2
    } else if (all(.is_one(Php))) {
      Php_type = 3
    } else {
      Php_type = 4
    }
    
    Ld <- value$lambda
    Ps <- value$psi
    Bt <- value$beta
    Ph <- value$phi
    
    if (all(nup == 1) & all(!.is_est(app)) & all(value$ap == 0)) {
      intercept_type = 1
      nu <- obs_moment$mu
      ap <- value$alpha
    } else {
      intercept_type = 2
      nu <- value$nu
      ap <- value$alpha
    }
    
    type <- penalty_fit$type
    gm <- penalty_fit$gm
    dt <- penalty_fit$dt
    
    maxit <- control$maxit
    epsilon <- control$epsilon
    
    Cyc <- obs_moment$Sg
    ey <- obs_moment$mu
    
    theta <-
      .sem_theta_cal(Ldp, Psp, Btp, Php, nup, app, Ld, Ps, Bt, Ph, nu, ap)
    thetap <-
      .sem_theta_cal(Ldp, Psp, Btp, Php, nup, app, Ldp, Psp, Btp, Php, nup, app)
    implied_moment <-
      .sem_implied_moment_cal(Ld, Ps, Bt, Ph, nu, ap)
    Sg <- implied_moment$Sg
    mu <- implied_moment$mu
    dpl <- Inf
    
    for (it in 1:maxit) {
      mis_moment <- .sem_estep(Cyc, ey, Sg, mu, Ld, Ps, Bt, Ph, nu, ap)
      Cet <- mis_moment$Cet
      Cyet <- mis_moment$Cyet
      eet <- mis_moment$eet
      if (intercept_type == 2) {
        nu <-
          .sem_mstep_nu(nup, nu, Ld, Ps, ey, eet, type, gm, dt, P, Psp, Psp_type)
        ap <-
          .sem_mstep_ap(app, ap, Bt, Ph, eet, type, gm, dt, M, Php_type)
      }
      Ld <-
        .sem_mstep_Ld(Ldp, Ld, Ps, nu, Cyet, Cet, ey, eet, type, gm, dt, P, Psp, Psp_type)
      Bt <-
        .sem_mstep_Bt(Btp, Bt, Ph, ap, Cet, eet, type, gm, dt, M, Php_type)
      Ps <-
        .sem_mstep_Ps(Psp, Ps, Ld, nu, Cyc, Cyet, Cet, ey, eet, type, gm, dt, P, Psp_type)
      Ph <-
        .sem_mstep_Ph(Php, Ph, Bt, ap, Cet, eet, type, gm, dt, M, Php_type)
      
      theta_new <-
        .sem_theta_cal(Ldp, Psp, Btp, Php, nup, app, Ld, Ps, Bt, Ph, nu, ap)
      thetap_new <-
        .sem_theta_cal(Ldp, Psp, Btp, Php, nup, app, Ldp, Psp, Btp, Php, nup, app)
      implied_moment <-
        .sem_implied_moment_cal(Ld, Ps, Bt, Ph, nu, ap)
      Sg <- implied_moment$Sg
      mu <- implied_moment$mu
      dml <- .sem_dml_cal(Cyc, ey, Sg, mu)
      rpl <- .sem_rpl_cal(theta_new, thetap_new, type, gm, dt)
      dpl_new <- dml + rpl
      if (mean(abs(theta_new - theta)) < epsilon) {
        theta <- theta_new
        dpl <- dpl_new
        break
      } else {
        theta <- theta_new
        dpl <- dpl_new
      }
    }
    Q <-
      sum(theta[is.na(thetap)] != 0, na.rm = T) + sum(thetap, na.rm = T)
    df <- P * (P + 3) / 2 - Q
    theta_mat <-
      list(
        Ld = Ld,
        Ps = Ps,
        Bt = Bt,
        Ph = Ph,
        nu = nu,
        ap = ap
      )
    rst_ecm <-
      list(
        dpl = dpl,
        dml = dml,
        Q = Q,
        df = df,
        it = it,
        theta = theta,
        theta_mat = theta_mat,
        implied_moment = implied_moment
      )
    return(rst_ecm)
  }
