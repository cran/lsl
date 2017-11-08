

.sem_theta_cal <- function(Ldp, Psp, Btp, Php, nup, app, 
                           Ld, Ps, Bt, Ph, nu, ap) {
  theta <- c(.vec(Ld)[.is_est(.vec(Ldp))], .vech(Ps)[.is_est(.vech(Psp))],
             .vec(Bt)[.is_est(.vec(Btp))], .vech(Ph)[.is_est(.vech(Php))],
             .vec(nu)[.is_est(.vec(nup))], .vec(ap)[.is_est(.vec(app))])
  return(theta)
}

.sem_theta_names <- function(Ldp, Psp, Btp, Php, nup, app) {
  if (prod(!.is_est(Ldp)) == 1) {
    Ld_names <- NULL
  } else {
    Ldp_idx <- which(.is_est(Ldp), arr.ind=TRUE)
    Ld_names <- paste("lambda", "[", Ldp_idx[, 1], ",", Ldp_idx[, 2], "]", sep="")    
  }
  if (prod(!.is_est(Psp)) == 1) {
    Ps_names <- NULL
  } else {
    Psp_idx <- which(.is_est(.ltri(Psp)), arr.ind=TRUE)
    Ps_names <- paste("psi", "[", Psp_idx[, 1], ",", Psp_idx[, 2], "]", sep="")    
  }
  if (prod(!.is_est(Btp)) == 1) {
    Bt_names <- NULL
  } else {
    Btp_idx <- which(.is_est(Btp), arr.ind=TRUE)
    Bt_names <- paste("beta", "[", Btp_idx[, 1], ",", Btp_idx[, 2], "]", sep="")    
  }
  if (prod(!.is_est(Php)) == 1) {
    Ph_names <- NULL
  } else {
    Php_idx <- which(.is_est(.ltri(Php)), arr.ind=TRUE)
    Ph_names <- paste("phi", "[", Php_idx[, 1], ",", Php_idx[, 2], "]", sep="")    
  }
  if (prod(!.is_est(nup)) == 1) {
    nu_names <- NULL
  } else {
    nup_idx <- which(.is_est(nup), arr.ind=TRUE)
    nu_names <- paste("nu", "[", nup_idx[, 1], ",", nup_idx[, 2], "]", sep="")  
  }  
  if (prod(!.is_est(app)) == 1) {
    ap_names <- NULL
  } else {
    app_idx <- which(.is_est(app), arr.ind=TRUE)
    ap_names <- paste("alpha", "[", app_idx[, 1], ",", app_idx[, 2], "]", sep="")  
  }  
  theta_names <- c(Ld_names, Ps_names, Bt_names, Ph_names, nu_names, ap_names)
  return(theta_names)
}


.sem_implied_moment_cal <- function(Ld, Ps, Bt, Ph, nu, ap) {
  IBtiv <- solve(diag(1, dim(Ld)[2]) - Bt)
  Xi <- IBtiv %*% Ph %*% t(IBtiv)
  Sg <- Ld %*% Xi %*% t(Ld) + Ps
  mu <- nu + Ld %*% IBtiv %*% ap
  implied_moment <- list(Sg = Sg, mu = mu)
  return(implied_moment)
}


.sem_dml_cal <- function(Cyc, ey, Sg, mu) {
  Sgiv <- solve(Sg)
  dml <- (-log(det(Cyc %*% Sgiv)) + sum(diag(Cyc %*% Sgiv)) - dim(Sg)[1] + t(ey - mu) %*% Sgiv %*% (ey - mu))
  return(dml)
}


.sem_rpl_cal <- function(theta, thetap, type, gm, dt) {
  gm <- gm
  dt <- dt
  if (sum(is.na(thetap)) > 0) {
    theta_pl <- c(theta[is.na(thetap)])
    if (type == 'l2') {
      rpl <- gm * sum((theta_pl)^2)
    } else if (type == 'l1') {
      rpl <- gm * sum(abs(theta_pl))
    } else if (type == 'scad') {
      rpl <- sum(gm * (abs(theta_pl) * (abs(theta_pl) <= gm))) + 
        sum((((gm * dt * abs(theta_pl) - 0.5 * (abs(theta_pl)^2 + gm^2))/(dt - 1)) * (abs(theta_pl) > gm & abs(theta_pl) <= (gm * dt)))) + 
        sum((((gm^2) * (dt^2 - 1) / (2 * (dt - 1))) * (abs(theta_pl) > (gm * dt))))
    } else if (type == 'mcp') {
      rpl <- sum(((gm * abs(theta_pl) - (abs(theta_pl)^2) / (2 * dt)) * (abs(theta_pl) <= (gm * dt)))) + 
        sum(((0.5 * dt * gm^2) * (abs(theta_pl) > (gm * dt))))
    } else {rpl <- 0}
  } else {rpl <- 0}
  return(rpl)
}


.sem_threshold <- function(ttq, type, wq, gm, dt) {
  if (type == "l2") {
    ttq <- ttq / (1 + gm * wq)
  } else if (type == "l1") {
    ttq <- sign(ttq) * max(abs(ttq) - 0.5 * gm * wq, 0)
  } else if (type == "scad") {
    if (abs(ttq) <= 0.5 * gm * (1+wq)) {
      ttq <- sign(ttq) * max(abs(ttq) - 0.5 *  gm * wq, 0)
    } else if (abs(ttq) > 0.5 * gm * (1 + wq) & abs(ttq) <= 0.5 * gm * dt) {
      ttq <- sign(ttq) * max(abs(ttq) - 0.5 * gm * wq * dt / (dt - 1), 0) / (1 - (wq / (dt - 1)))
    } else {}
  } else if (type == "mcp") {
    if (abs(ttq) <= gm * dt) {
      ttq <- sign(ttq) * max(abs(ttq) - 0.5 * gm * wq, 0) / (1 - (wq / (2 * dt)))
    } else {}
  } else {}
  return(ttq)
}
