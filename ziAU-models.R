library(sp)
library(GWmodel)
library(caret)

set.seed(123)

## DATA IMPORT
DB=read.csv("https://raw.githubusercontent.com/ProfNascimento/ziAU/refs/heads/main/Cobalt.csv",sep=";")

#-------------------------------------------------------------------------------#
# (HOLD OUT)
# (PCA) --Scaled--

fit.pca <- prcomp(DB[,-c(1:4,18)], scale. = T, center = T)

# SET + FULL COMPONENTES
DB0=cbind(V1=DB$Co/10,                 # TRANSFORMAR Y
          fit.pca$x,                        # PCs (Nuevos Xs)
          scale(DB[,-c(1:4,18)]),   # Elementos Geoquimicos Normalizados
          scale(DB[,1:2]),          # Lat,Lon Normalizados
          DB[,4])                   # Zona

DB0=as.data.frame(DB0)
str(DB0)

indices <- caret::createDataPartition(DB0$V1, p = 0.7, list = FALSE)
DB_train <- DB0[indices, ]
DB_test <- DB0[-indices, ]

#------------------------------------------------------------------------------#

# GWR-ZIAU MODEL
#######################################################################
#  GWR–ZIAU  — PART 1: utilidades, densidades y familias
#######################################################################

## ---------- utilidades básicas  ----
if (!exists("expit", mode = "function"))
  expit  <- function(x) plogis(x)

if (!exists("clip01", mode = "function"))
  clip01 <- function(v) pmin(pmax(v, 1e-6), 1 - 1e-6)

## ---------- make.link() y print.family  --------------------
if (!exists("family", mode = "function"))
  family <- function(object, ...) UseMethod("family")

if (!exists("print.family", mode = "function")) {
  print.family <- function(x, ...){
    cat("\nLink:", x$link, "\n\n"); invisible(x)
  }
}

if (!exists("make.link", mode = "function")) {
  make.link <- function(link){
    switch(link,
           "logit" = {lf<-function(mu).Call(C_logit_link,mu)
           li<-function(et).Call(C_logit_linkinv,et)
           me<-function(et).Call(C_logit_mu_eta,et); ve<-TRUE},
           "probit"= {lf<-qnorm; li<-function(et){
             thr<- -qnorm(.Machine$double.eps)
             pnorm(pmax(pmin(et,thr),-thr))}
           me<-function(et)pmax(dnorm(et),.Machine$double.eps); ve<-TRUE},
           "cloglog"={lf<-function(mu)log(-log1p(-mu))
           li<-function(et)pmax(pmin(-expm1(-exp(pmin(et,700))),
                                     1-.Machine$double.eps),
                                .Machine$double.eps)
           me<-function(et)pmax(exp(et)*exp(-exp(et)),
                                .Machine$double.eps); ve<-TRUE},
           "identity"={lf<-identity; li<-identity; me<-function(et)1; ve<-TRUE},
           "log"={lf<-log; li<-function(et)pmax(exp(et),.Machine$double.eps)
           me<-function(et)pmax(exp(et),.Machine$double.eps); ve<-TRUE},
           "inverse"={lf<-function(mu)1/mu; li<-function(et)1/et
           me<-function(et)-1/et^2; ve<-function(et)all(et!=0)},
           stop(gettextf("%s link not recognised", sQuote(link)), domain=NA))
    environment(lf)<-environment(li)<-environment(me)<-environment(ve)<-
      asNamespace("stats")
    structure(list(linkfun=lf,linkinv=li,mu.eta=me,valideta=ve,name=link),
              class="link-glm")
  }
}

#######################################################################
#  Densidades: Alpha-Unit (AU) y Zero-Inflated Alpha-Unit (ZIAU)
#######################################################################

# --- Densidad AU: soporte (0,1). Devuelve densidad o log-densidad ----
dau <- function(y, alpha, log = FALSE){
  stopifnot(all(alpha > 0))
  y_ok <- (y > 0 & y < 1)
  # f_AU(y; alpha) = 2/(y*alpha) * (log(y)/alpha)^2 * phi(log(y)/alpha)
  logfy <- rep(-Inf, length(y))
  if (any(y_ok)) {
    z <- log(y[y_ok]) / alpha
    logfy[y_ok] <- log(2) - log(alpha) - log(y[y_ok]) + 2*log(abs(log(y[y_ok]))/alpha) +
      dnorm(z, log = TRUE)
  }
  if (log) logfy else exp(logfy)
}

# --- Densidad ZIAU inflada en 0: mezcla en y=0 con prob. nu -----------
dziau <- function(y, alpha, nu, log = FALSE){
  stopifnot(all(alpha > 0), all(nu >= 0), all(nu < 1))
  # Si y==0 → masa puntual nu
  # Si y∈(0,1) → (1-nu)*AU(y; alpha)
  is0   <- (y == 0)
  is01  <- (y > 0 & y < 1)
  out_l <- rep(-Inf, length(y))
  if (any(is0))  out_l[is0]  <- log(pmax(nu[is0], .Machine$double.eps))
  if (any(is01)) out_l[is01] <- log1p(-nu[is01]) + dau(y[is01], alpha[is01], log = TRUE)
  # (y>=1) o (y<0) → densidad 0
  if (log) out_l else exp(out_l)
}

#######################################################################
#  Familias: AU y ZIAU 
#######################################################################

# Familia AU (solo para referencia; α es el parámetro de forma)
au <- function(link.alpha = "log"){
  la <- make.link(link.alpha)
  structure(list(
    family       = "au",
    link         = link.alpha,
    linkfun      = la$linkfun,   # para α
    linkinv      = la$linkinv,
    alpha.eta    = la$mu.eta,    # dα/dη (con η = link(α))
    # Funciones auxiliares opcionales (no usadas por el motor por ahora)
    logdens      = function(y, alpha) dau(y, alpha, log = TRUE)
  ), class = "family")
}

# Familia ZIAU: ν en (0,1) (global), α>0 (local)
ziau <- function(link.alpha = "log", link.nu = "logit"){
  la <- make.link(link.alpha)  # α
  ln <- make.link(link.nu)     # ν
  structure(list(
    family         = "ziau",
    link.alpha     = link.alpha,
    link.nu        = link.nu,
    alpha.linkfun  = la$linkfun,
    alpha.linkinv  = la$linkinv,
    alpha.eta      = la$mu.eta,   # dα/dη
    nu.linkfun     = ln$linkfun,
    nu.linkinv     = ln$linkinv,
    nu.eta         = ln$mu.eta,   # dν/dη
    # Densidad (log) para usar en la verosimilitud:
    logdens        = function(y, alpha, nu) dziau(y, alpha, nu, log = TRUE)
    # Nota: la varianza/esperanza cerradas no son necesarias aquí; el motor
    #       GWR-ZIAU manejará pesos específicos para α y ν.
  ), class = "family")
}

#######################################################################

# Requiere PARTE 1: expit, clip01, make.link, dau(), dziau() (si tienes
# dziau solo para c=0, abajo incluyo una versión generalizada opcional).

## -------------------------------------------------------------------
##  utilitarios de álgebra local
## -------------------------------------------------------------------
Ci_local <- function(X, w, ridge = 1e-6){
  XtW  <- sweep(t(X), 2, w, "*")             # (p × n)
  XtWX <- XtW %*% X + diag(ridge, ncol(X))   # (p × p)
  Ci   <- solve(XtWX, XtW)                   # (p × n)   (X'WX)^(-1) X'W
  hii  <- rowSums(X * t(Ci))                 # diag(S)≈leverages
  list(Ci = Ci, hii = hii)
}
Ci_mat <- function(X, w, ridge = 1e-6) Ci_local(X, w, ridge)$Ci

gw_reg <- function(X, y, w, hatmatrix = FALSE, focus = 1L, ridge = 1e-6){
  XtW  <- sweep(t(X), 2, w, "*")
  XtWX <- XtW %*% X + diag(ridge, ncol(X))
  beta <- as.vector(solve(XtWX, XtW %*% y))
  if(!hatmatrix) return(list(beta))
  Ci   <- solve(XtWX, XtW)
  hii  <- sum(X[focus,] * Ci[,focus])
  list(beta, Ci = Ci, hii = hii)
}
gw_fitted <- function(X, B) rowSums(X * B)

## -------------------------------------------------------------------
##  Momentos cerrados de AU y momentos de ZIAU (inflación en c=0/1)
## -------------------------------------------------------------------
# μ_r = E_AU(Y^r; α) — forma cerrada
.au_mr_closed <- function(alpha, r = 1){
  stopifnot(all(alpha > 0), length(r) == 1, r > 0)
  ar   <- r * alpha
  term <- ((ar^2 + 1) * (1 - pnorm(ar)) - ar * dnorm(ar))
  mu_r <- 2 * exp(0.5 * ar^2) * term
  pmin(pmax(mu_r, 0), 1)  # guarda numérica al soporte [0,1]
}
.au_m1 <- function(alpha) .au_mr_closed(alpha, r = 1)
.au_m2 <- function(alpha) .au_mr_closed(alpha, r = 2)

# Media y varianza de ZIAU para inflación en c ∈ {0,1}
.ziau_moments <- function(alpha, nu, c_inflate = 0){
  m1 <- nu * c_inflate + (1 - nu) * .au_m1(alpha)
  m2 <- nu * (c_inflate^2) + (1 - nu) * .au_m2(alpha)
  v  <- pmax(m2 - m1^2, 0)
  list(mean = m1, var = v, m2 = m2)
}

# (Opcional) Generalización de dziau a c ∈ {0,1}, por si tu PARTE 1 solo trae c=0
dziau <- function(y, alpha, nu, c = 0, log = FALSE){
  stopifnot(all(alpha > 0), all(nu >= 0 & nu < 1), c %in% c(0,1))
  is_c  <- (y == c)
  is01  <- (y > 0 & y < 1)
  out_l <- rep(-Inf, length(y))
  if (any(is_c))  out_l[is_c]  <- log(pmax(nu[is_c], .Machine$double.eps))
  if (any(is01))  out_l[is01] <- log1p(-nu[is01]) + dau(y[is01], alpha[is01], log = TRUE)
  if (log) out_l else exp(out_l)
}

## -------------------------------------------------------------------
##  núcleo local para α (NR en β_α(s) con enlace log)
##  Maximiza: sum_i w_i [ -3 log α_i - 0.5 (log(y_i)/α_i)^2 ],  y_i∈(0,1)
##  con: α_i = exp(η_i), η_i = X_i β(s)
##  Derivadas cerradas:
##   s_i = ∂ℓ/∂η_i = -3 + (log y_i)^2 / α_i^2
##   d s_i / dη_i  = -2 (log y_i)^2 / α_i^2
## -------------------------------------------------------------------
.local_alpha_fit <- function(X, y, w, beta0 = NULL,
                             tol = 1e-6, maxiter = 100, ridge = 1e-6){
  n  <- nrow(X); p <- ncol(X)
  is01 <- (y > 0 & y < 1)
  if(!any(is01))
    return(list(beta = rep(NA_real_, p),
                se   = rep(NA_real_, p),
                conv = FALSE,
                iter = 0))
  
  # submuestras informativas para α
  Xs <- X[is01, , drop = FALSE]
  ws <- w[is01]
  t  <- log(y[is01])
  
  # inicial: α0^2 = mean_w(t^2)/3  ⇒ η0 = log α0
  if(is.null(beta0)){
    a0 <- sqrt( sum(ws * t^2) / (3 * sum(ws)) )
    eta0 <- log(pmax(a0, 1e-6))
    zs <- rep(eta0, length(ws))
    XtW <- sweep(t(Xs), 2, ws, "*")
    beta  <- as.vector(solve(XtW %*% Xs + diag(ridge, p), XtW %*% zs))
  } else beta <- beta0
  
  it <- 0L; conv <- FALSE
  repeat{
    eta   <- as.vector(Xs %*% beta)
    alpha <- exp(eta)
    
    s     <- -3 + (t^2)/(alpha^2)                 # score en η
    ds    <- -2 * (t^2)/(alpha^2)                 # derivada del score
    
    g     <- t(Xs) %*% (ws * s)                   # gradiente (p×1)
    H     <- t(Xs) %*% (Xs * as.numeric(ws * ds)) # Hessiano (p×p)
    
    # Regularizar e iterar Newton (maximización): β_new = β - H^{-1} g
    Hreg  <- H - diag(diag(H)) + diag(diag(H) - ridge, p)
    step  <- tryCatch(solve(Hreg, g), error = function(e) rep(0, p))
    beta_new <- beta - as.vector(step)
    
    den  <- max(1, sqrt(sum(beta^2)))
    if(sqrt(sum((beta_new - beta)^2)) / den < tol){ conv <- TRUE; beta <- beta_new; break }
    beta <- beta_new
    it   <- it + 1L
    if(it >= maxiter) break
  }
  
  # SE(β): var ≈ [ -H(β̂) ]^{-1}
  eta_hat   <- as.vector(Xs %*% beta); alpha_hat <- exp(eta_hat)
  ds_hat    <- -2 * (t^2) / (alpha_hat^2)
  H_hat     <- t(Xs) %*% (Xs * as.numeric(ws * ds_hat))
  Info_pos  <- -(H_hat) + diag(ridge, p)
  Vbeta     <- tryCatch(solve(Info_pos), error=function(e) diag(NA_real_, p))
  se        <- sqrt(pmax(diag(Vbeta), 0))
  
  list(beta = beta, se = se, conv = conv, iter = it)
}

## -------------------------------------------------------------------
##  ν global (logit) mediante GLM binomial estándar
## -------------------------------------------------------------------
.fit_nu_global <- function(X_nu, y, link = "logit"){
  z <- as.integer(y == 0)    # 1 si cero “inflado”, 0 si (0,1)
  fam <- binomial(link = link)
  fit <- stats::glm.fit(x = X_nu, y = z, family = fam)
  beta <- coef(fit); beta[is.na(beta)] <- 0
  nu   <- clip01(fam$linkinv(drop(X_nu %*% beta)))
  list(beta = beta, nu = nu, fitted = fit)
}

## -------------------------------------------------------------------
##  MOTOR: GWR–ZIAU  (α local con kernel, ν global)
## -------------------------------------------------------------------
gwr.ziau <- function(y, x,                        # diseño para α (local)
                     regression.points,
                     W1.mat, W2.mat,             # pesos α y predicción
                     hatmatrix   = TRUE,
                     tol         = 1e-5,
                     maxiter     = 200,
                     link.alpha  = "log",
                     link.nu     = "logit",
                     formula.nu  = NULL,         # diseño para ν global (opcional)
                     dp.locat,
                     p4s         = NA_character_,
                     nu.cut      = 0.5,          # umbral para clasificación a c
                     ridge       = 1e-6,
                     c_inflate   = 0)            # 0 ó 1 (por defecto 0)
{
  eps <- 1e-8
  
  ## --- coordenadas --------------------------------------------------
  if (missing(dp.locat) || is.null(dp.locat)) {
    if (is(regression.points, "Spatial"))
      dp.locat <- coordinates(regression.points)
    else stop("'dp.locat' ausente y no deducible")
  }
  if (is(regression.points, "Spatial"))
    p4s <- proj4string(regression.points)
  
  n <- length(y); k <- ncol(x)
  is01 <- (y > 0 & y < 1)
  
  ## --- diseño para ν (global) --------------------------------------
  X_nu <- if(!is.null(formula.nu)){
    mf.nu <- model.frame(formula.nu, data = as.data.frame(x))
    model.matrix(attr(mf.nu, "terms"), mf.nu)
  } else {
    x
  }
  
  ## --- ajuste global de ν ------------------------------------------
  nu_fit   <- .fit_nu_global(X_nu, y, link = link.nu)
  nu       <- nu_fit$nu
  beta_nu  <- nu_fit$beta
  
  ## --- ajuste local de α (uno por punto i) -------------------------
  Balpha <- matrix(NA_real_, n, k)
  SEa    <- matrix(NA_real_, n, k)
  alpha  <- rep(NA_real_, n)
  hii    <- numeric(n)
  
  for(i in 1:n){
    w_i <- W1.mat[, i] * as.numeric(is01)    # solo y∈(0,1) aporta a α
    if(sum(w_i) <= 0){
      Balpha[i, ] <- NA_real_; SEa[i, ] <- NA_real_; alpha[i] <- NA_real_
      next
    }
    fit_i <- .local_alpha_fit(X = x, y = y, w = w_i,
                              beta0 = NULL, tol = tol,
                              maxiter = maxiter, ridge = ridge)
    Balpha[i, ] <- fit_i$beta
    SEa[i, ]    <- fit_i$se
    # α̂ EN EL PUNTO i (NO sobre toda la muestra):
    alpha[i]    <- exp(drop(x[i, ] %*% fit_i$beta))
    # leverage proxy en i
    Ci_i        <- Ci_mat(x, w_i, ridge = ridge)
    hii[i]      <- sum(x[i, ] * Ci_i[, i])
  }
  
  ## --- predicción, residuales, log-lik (FÓRMULAS CERRADAS) ----------
  mom   <- .ziau_moments(alpha, nu, c_inflate = c_inflate)
  Ey    <- mom$mean
  yhat  <- ifelse(nu >= nu.cut, c_inflate, Ey)  # si cortas por ν, respeta c
  resid <- y - yhat
  ll    <- sum(dziau(y, alpha = alpha, nu = nu, c = c_inflate, log = TRUE))
  
  ## --- diagnóstico AIC/AICc (traza S aprox por suma h_ii) -----------
  trS <- sum(pmax(pmin(hii, 1), 0), na.rm = TRUE)
  AIC  <- -2 * ll + 2 * trS
  AICc <- AIC + 2 * trS * (trS + 1) / pmax(n - trS - 1, 1)
  
  ## --- salida tipo SpatialPointsDataFrame ---------------------------
  df <- data.frame(
    y     = y,
    Ey    = Ey,
    yhat  = yhat,
    resid = resid,
    alpha = alpha,
    nu    = nu
  )
  colnames(Balpha) <- colnames(x)
  colnames(SEa)    <- paste0(colnames(x), "_SE")
  df <- cbind(df, Balpha, SEa)
  
  SDF <- SpatialPointsDataFrame(dp.locat, df, proj4string = CRS(p4s))
  
  list(
    SDF           = SDF,
    glms          = list(beta.alpha = Balpha, beta.nu = beta_nu),
    GW.diagnostic = list(gw.deviance = -2 * ll,
                         AIC = AIC, AICc = AICc,
                         n.par = trS),
    alpha         = alpha,
    nu            = nu,
    nu.cut        = nu.cut,
    c_inflate     = c_inflate,
    logLik        = ll
  )
}



#######################################################################
#  GWR–ZIAU — PARTE 3: versiones .wt, CV/AICc, bw, y wrapper básico
##############################

## ---------------------------------------------------------------
##  Pesos kernel (define solo si no existe en el entorno)
## ---------------------------------------------------------------
if (!exists("gw.weight", mode = "function")) {
  gw.weight <- function(d, bw, kernel = "bisquare", adaptive = FALSE){
    if (adaptive){ kw <- (rank(d, ties.method="first") - 1L) / bw }
    else          { kw <- d / bw }
    switch(kernel,
           gauss     = exp(-.5 * kw^2),
           bisquare  = (abs(kw) < 1) * (1 - kw^2)^2,
           tricube   = (abs(kw) < 1) * (1 - abs(kw)^3)^3,
           stop("kernel desconocido"))
  }
}

## -------------------------------------------------------------------
##  .WT para ZIAU (rápido para AICc/CV): una pasada IRLS en α local
## -------------------------------------------------------------------
# Devuelve:
#   - wt2  : pesos efectivos para α (≈ -ds/dη en (0,1); 0 en y==0 o 1)
#   - llik : log-verosimilitud total con α̂ y ν̂ (ν global)
#   - Ey   : E[Y] del modelo (cerrado)
#   - alpha, nu : parámetros estimados (para debug/inspección)
gwr.ziau.wt <- function(y, X, bw, W,                 # W: n×n (columnas = focos)
                        kernel = "bisquare",
                        adaptive = FALSE,
                        link.alpha = "log",
                        link.nu    = "logit",
                        X_nu       = NULL,           # diseño para ν (global)
                        c_inflate  = 0,
                        tol        = 1e-6)
{
  n <- length(y)
  is01 <- (y > 0 & y < 1)
  
  # ---- ν global (una vez) ------------------------------------------
  if (is.null(X_nu)) X_nu <- X
  nu <- .fit_nu_global(X_nu, y, link = link.nu)$nu
  
  # ---- inicial para α: α0^2 = mean(t^2)/3, con t = log y ------------
  t      <- log(pmax(y, 1e-12))
  a0     <- rep(NA_real_, n)
  a0[is01]  <- sqrt( mean( t[is01]^2 ) / 3 )
  a0[!is01] <- median(a0[is01], na.rm = TRUE)
  eta0   <- log(pmax(a0, 1e-6))
  
  # ---- una iteración IRLS equivalente a NR en η:
  #      s  = -3 + t^2/α^2;  ds = -2 t^2/α^2;  z = η + s/ds; w_eff = -ds
  z_all <- eta0                    # ← FIX: antes era rep(eta0, n)
  w_eff <- rep(0, n)
  if (any(is01)) {
    s        <- -3 + (t[is01]^2)/(a0[is01]^2)
    ds       <- -2 * (t[is01]^2)/(a0[is01]^2)
    z_all[is01] <- eta0[is01] + s/ds
    w_eff[is01] <- -ds
  }
  
  # ---- β_α(s) local (WLS) para cada foco i -------------------------
  Balpha <- matrix(NA_real_, n, ncol(X))
  for(i in 1:n){
    wi <- W[, i] * w_eff
    if(sum(wi) <= 0){ Balpha[i,] <- 0; next }
    Balpha[i,] <- gw_reg(X, z_all, wi, hatmatrix = FALSE, focus = i)[[1]]
  }
  
  # α̂_i = exp(x_i' β̂_i)  (vector completo)
  eta_hat <- as.vector(rowSums(X * Balpha))
  alpha   <- exp(eta_hat)
  
  # ---- pesos para hat-matrix/AICc
  wt2 <- w_eff
  
  # ---- E[Y] (cerrado) y log-verosimilitud total --------------------
  Ey    <- nu * c_inflate + (1 - nu) * .au_m1(alpha)
  llik  <- sum(dziau(y, alpha = alpha, nu = nu, c = c_inflate, log = TRUE))
  
  list(wt2 = wt2, llik = llik, Ey = Ey, alpha = alpha, nu = nu)
}

## -------------------------------------------------------------------
##  Contribution CV (LOOCV fast, ν global)
## -------------------------------------------------------------------
ggwr.ziau.cv.contrib <- function(bw, X, Y,
                                 kernel="bisquare", adaptive=FALSE,
                                 dp.locat, p=2, theta=0, longlat=FALSE,
                                 dMat=NULL, link.alpha="log", link.nu="logit",
                                 X_nu=NULL, c_inflate = 0)
{
  n <- nrow(dp.locat)
  is01 <- (Y > 0 & Y < 1)
  
  # ν global una vez
  if (is.null(X_nu)) X_nu <- X
  nu <- .fit_nu_global(X_nu, Y, link = link.nu)$nu
  
  # Inicial para α (como en .wt)
  t      <- log(pmax(Y, 1e-12))
  a0     <- rep(NA_real_, n)
  a0[is01]  <- sqrt( mean(t[is01]^2) / 3 )
  a0[!is01] <- median(a0[is01], na.rm = TRUE)
  eta0   <- log(pmax(a0, 1e-6))
  
  z_all <- eta0                    # ← FIX: antes era rep(eta0, n)
  w_eff <- rep(0, n)
  if (any(is01)) {
    s        <- -3 + (t[is01]^2)/(a0[is01]^2)
    ds       <- -2 * (t[is01]^2)/(a0[is01]^2)
    z_all[is01] <- eta0[is01] + s/ds
    w_eff[is01] <- -ds
  }
  
  CV <- numeric(n)
  for(i in 1:n){
    dist <- if(is.null(dMat)) gw.dist(dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    else dMat[, i]
    wi   <- gw.weight(dist, bw, kernel, adaptive)
    wi[i] <- 0  # LOO
    if(sum(wi)==0){ CV[i] <- NA; next }
    # fit local en η (una pasada)
    beta_i  <- gw_reg(X, z_all, wi * w_eff, hatmatrix = FALSE, focus = i)[[1]]
    eta_i   <- drop(X[i,] %*% beta_i)
    alpha_i <- exp(eta_i)
    Ey_i    <- nu[i] * c_inflate + (1 - nu[i]) * .au_m1(alpha_i)
    CV[i]   <- Y[i] - Ey_i
  }
  CV
}

## -------------------------------------------------------------------
##  Objective Genéric CV (optimize)
## -------------------------------------------------------------------
ggwr.ziau.cv <- function(bw, X, Y,
                         kernel="bisquare", adaptive=FALSE,
                         dp.locat, p=2, theta=0, longlat=FALSE,
                         dMat=NULL, link.alpha="log", link.nu="logit",
                         X_nu=NULL, c_inflate = 0){
  cv_vec <- ggwr.ziau.cv.contrib(bw, X, Y, kernel, adaptive,
                                 dp.locat, p, theta, longlat,
                                 dMat, link.alpha, link.nu,
                                 X_nu, c_inflate)
  if(!all(is.finite(cv_vec))) return(.Machine$double.xmax)
  sum(cv_vec^2)
}

## -------------------------------------------------------------------
##  AICc (optimize) using gwr.ziau.wt
## -------------------------------------------------------------------
ggwr.ziau.aic <- function(bw, X, Y,
                          kernel="bisquare", adaptive=FALSE,
                          dp.locat, p=2, theta=0, longlat=FALSE,
                          dMat=NULL, link.alpha="log", link.nu="logit",
                          X_nu=NULL, c_inflate = 0)
{
  n <- nrow(dp.locat)
  
  # Matriz de pesos W (n×n) para α local
  W <- matrix(0, n, n)
  for(i in 1:n){
    dist <- if(is.null(dMat)) gw.dist(dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    else dMat[, i]
    W[, i] <- gw.weight(dist, bw, kernel, adaptive)
  }
  
  # Fit rápido (.wt)
  res  <- gwr.ziau.wt(Y, X, bw, W, kernel, adaptive,
                      link.alpha, link.nu, X_nu, c_inflate)
  wt2  <- res$wt2
  llik <- res$llik
  
  # traza S ~ suma de leverages con pesos W[,i]*wt2
  hii <- numeric(n)
  for(i in 1:n){
    Ci <- Ci_mat(X, W[, i] * wt2)
    hii[i] <- sum(X[i, ] * Ci[, i])
  }
  trS <- sum(hii, na.rm = TRUE)
  
  AIC  <- -2 * llik + 2 * trS
  AICc <- AIC + 2 * trS * (trS + 1) / max(n - trS - 1, 1)
  if(!is.finite(AICc)) AICc <- .Machine$double.xmax
  AICc
}

## -------------------------------------------------------------------
##  Band for the ZIAU
## -------------------------------------------------------------------
bw.ggwr.ziau <- function(formula, data, approach = "CV",
                         kernel = "bisquare", adaptive = FALSE,
                         p = 2, theta = 0, longlat = FALSE, dMat = NULL,
                         link.alpha = "log", link.nu = "logit",
                         formula.nu = NULL, c_inflate = 0)
{
  if (!is(data, "Spatial"))
    stop("‘data’ debe ser Spatial*DataFrame")
  
  loc <- coordinates(data)
  mf  <- model.frame(formula, data = as(data, "data.frame"))
  X   <- model.matrix(attr(mf, "terms"), mf)
  Y   <- model.response(mf)
  
  # diseño para ν global (si se especifica formula.nu)
  X_nu <- if(!is.null(formula.nu)){
    mf.nu <- model.frame(formula.nu, data = as.data.frame(data))
    model.matrix(attr(mf.nu, "terms"), mf.nu)
  } else X
  
  # rango inicial heurístico
  upper <- if(is.null(dMat)) max(gw.dist(loc, focus = 1, p = p, theta = theta, longlat = longlat))
  else max(dMat)
  lower <- upper / 50
  
  score <- if (tolower(approach) == "cv") ggwr.ziau.cv else ggwr.ziau.aic
  optimize(score, c(lower, upper), tol = 0.5,
           X = X, Y = Y,
           kernel = kernel, adaptive = adaptive,
           dp.locat = loc, p = p, theta = theta, longlat = longlat,
           dMat = dMat, link.alpha = link.alpha, link.nu = link.nu,
           X_nu = X_nu, c_inflate = c_inflate)$minimum
}

## -------------------------------------------------------------------
##  Wrapper principal: ggwr.ziau.basic 
## -------------------------------------------------------------------
ggwr.ziau.basic <- function(formula, data,                 # especificación
                            regression.points,
                            bw,
                            kernel      = "bisquare",
                            adaptive    = FALSE,
                            cv          = TRUE,
                            tol         = 1e-5,
                            maxiter     = 200,
                            p           = 2,
                            theta       = 0,
                            longlat     = FALSE,
                            dMat        = NULL,
                            dMat1       = NULL,
                            link.alpha  = "log",
                            link.nu     = "logit",
                            formula.nu  = NULL,           # diseño ν global (opcional)
                            nu.cut      = 0.5,
                            c_inflate   = 0)              # 0 ó 1 (default 0)
{
  timings   <- list(start = Sys.time())
  this.call <- match.call()
  
  if (!is(data, "Spatial"))
    stop("‘data’ debe ser un Spatial*DataFrame")
  
  ## ---------- datos y diseño ----------------------------------------
  p4s    <- proj4string(data)
  dp.loc <- coordinates(data)
  mf     <- model.frame(formula, data = as(data, "data.frame"))
  y      <- model.response(mf)
  X      <- model.matrix(attr(mf, "terms"), mf)
  n      <- nrow(X)
  
  ## ---------- diseño para ν (global) --------------------------------
  X_nu <- if(!is.null(formula.nu)){
    mf.nu <- model.frame(formula.nu, data = as.data.frame(data))
    model.matrix(attr(mf.nu, "terms"), mf.nu)
  } else X
  
  ## ---------- puntos de regresión -----------------------------------
  if (missing(regression.points)) {
    rp.loc    <- dp.loc
    hatmatrix <- TRUE
  } else {
    if (!is(regression.points, "Spatial"))
      stop("‘regression.points’ debe ser Spatial*")
    rp.loc    <- coordinates(regression.points)
    hatmatrix <- FALSE
  }
  nrp <- nrow(rp.loc)
  
  ## ---------- distancias (si no se pasan) ---------------------------
  if (is.null(dMat))
    dMat  <- gw.dist(dp.loc, rp.loc, p = p, theta = theta, longlat = longlat)
  if (is.null(dMat1))
    dMat1 <- gw.dist(dp.loc, dp.loc,  p = p, theta = theta, longlat = longlat)
  
  ## ---------- matrices de pesos -------------------------------------
  W1 <- matrix(0, n, n)
  for (i in seq_len(n))
    W1[, i] <- gw.weight(dMat1[, i], bw, kernel, adaptive)
  
  W2 <- matrix(0, n, nrp)
  for (i in seq_len(nrp))
    W2[, i] <- gw.weight(dMat[,  i], bw, kernel, adaptive)
  
  ## ---------- motor principal (PARTE 2) ------------------------------
  
  mod <- gwr.ziau(y, X,
                  regression.points = data,
                  W1.mat = W1, W2.mat = W2,
                  hatmatrix = hatmatrix,
                  tol = tol, maxiter = maxiter,
                  link.alpha = link.alpha, link.nu = link.nu,
                  formula.nu = NULL,
                  dp.locat = dp.loc, p4s = p4s,
                  nu.cut = nu.cut, c_inflate = c_inflate)
  
  ## ---------- actualizar ν con formula.nu (si se especifica) --------
  if (!is.null(formula.nu)) {
    nu_fit <- .fit_nu_global(X_nu, y, link = link.nu)
    mod$nu <- nu_fit$nu
    # recomputar Ey, yhat, resid con c_inflate
    df <- as.data.frame(mod$SDF)
    Ey <- nu_fit$nu * c_inflate + (1 - nu_fit$nu) * .au_m1(df$alpha)
    df$Ey    <- Ey
    df$yhat  <- ifelse(nu_fit$nu >= nu.cut, c_inflate, Ey)
    df$resid <- df$y - df$yhat
    mod$SDF@data <- df
    mod$glms$beta.nu <- nu_fit$beta
  }
  
  ## ---------- CV global (solo hat-matrix) ----------------------------
  CV <- if (cv && hatmatrix)
    ggwr.ziau.cv.contrib(bw, X, y, kernel, adaptive,
                         dp.locat = dp.loc,
                         p = p, theta = theta, longlat = longlat,
                         dMat = dMat1, link.alpha = link.alpha,
                         link.nu = link.nu, X_nu = X_nu, c_inflate = c_inflate)
  else numeric(n)
  
  timings$stop <- Sys.time()
  
  structure(
    list(GW.arguments = list(formula   = formula,
                             bw        = bw,
                             kernel    = kernel,
                             adaptive  = adaptive,
                             family    = "ziau",
                             nu.cut    = nu.cut,
                             c_inflate = c_inflate,
                             link.alpha= link.alpha,
                             link.nu   = link.nu),
         SDF            = mod$SDF,
         glms           = mod$glms,
         GW.diagnostic  = mod$GW.diagnostic,
         CV             = CV,
         timings        = timings,
         this.call      = this.call),
    class = "ggwrm"
  )
}

#######################################################################
#  MÉTODOS S3 (INFOVIS)
#######################################################################
## ---------- print.ggwrm --------------------------------------------
print.ggwrm <- function(x, ...){
  cat("\n─ GWR-", toupper(x$GW.arguments$family),
      "  (kernel:", x$GW.arguments$kernel, ")\n", sep="")
  cat(" Fórmula : "); print(x$GW.arguments$formula)
  cat(" N puntos:", nrow(x$SDF), 
      "   Bandwidth:", format(x$GW.arguments$bw, digits = 5),
      if(x$GW.arguments$adaptive) "(adaptativo)" else "(fijo)", "\n")
  cat(" Tiempo  :", round(difftime(x$timings$stop,
                                   x$timings$start,
                                   units="secs"), 2), "seg\n")
  cat(" --- Diagnóstico global ---\n")
  print(x$GW.diagnostic)
  invisible(x)
}

## ---------- as.data.frame.ggwrm ------------------------------------
as.data.frame.ggwrm <- function(x, row.names = NULL, optional = FALSE, ...){
  df <- as.data.frame(x$SDF)
  attr(df, "coords") <- coordinates(x$SDF)   # por si se necesita luego
  df
}

## ---------- plot.ggwrm  (mapa rápido con ggplot2) -------------------

plot.ggwrm <- function(x,
                       what = c("beta","yhat","resid","mu","nu"),
                       param = 1,
                       palette = scales::viridis_pal(),
                       size = 2,
                       ...){
  what <- match.arg(what)
  df   <- as.data.frame(x)
  crd  <- attr(df,"coords")
  df$X <- crd[,1]; df$Y <- crd[,2]
  
  if(what == "beta"){
    betas <- grep("^Bmu\\.|^Bnu\\.|^\\(Intercept\\)$|^[^_]+$", names(df),
                  value = TRUE)
    if(param < 1 || param > length(betas))
      stop("param fuera de rango (hay ", length(betas), " coeficientes)")
    zvar <- betas[param]
    ttl  <- paste0("Coef. local: ", zvar)
  } else if(what == "yhat"){
    zvar <- if("yhat" %in% names(df)) "yhat" else "Ey"
    ttl  <- "ŷ / E[ŷ]"
  } else if(what == "resid"){
    zvar <- "resid"; ttl <- "Residuales"
  } else {     # mu o nu
    zvar <- what; ttl <- paste0("ĥ ", what)
    if(!zvar %in% names(df))
      stop("La columna '", zvar,"' no existe en SDF")
  }
  
  ggplot(df, aes(X, Y, colour = .data[[zvar]])) +
    geom_point(size = size) +
    coord_equal() +
    scale_colour_gradientn(colours = palette(256)) +
    labs(title = ttl, colour = "") +
    theme_bw()
}

#------------------------------------------------------------------------------#
# Function Predictions 
## ===============================================================
##  MÉTODO S3: predict.ggwrm  (GWR–ZIAU, MLE local, sin nu.cut)
## ===============================================================

# --- helpers numéricos (locales al archivo; no exportar) ---
.au_m1_stable <- function(a){
  ltail <- pnorm(-a, log.p = TRUE)
  t1    <- (a^2 + 1) * exp(ltail + 0.5 * a^2)
  t2    <- a / sqrt(2*pi)
  out   <- 2 * (t1 - t2)
  pmin(pmax(out, 0), 1)
}
.au_logpdf <- function(y, alpha, eps=1e-12){
  y  <- pmin(pmax(y, eps), 1 - eps)
  z  <- log(y) / alpha
  val <- log(2) - 0.5*log(2*pi) + 2*log(abs(log(y))) - 3*log(alpha) - 0.5*z^2 - log(y)
  val[!is.finite(val)] <- -Inf
  val
}
.au_m1_inv <- function(target, a_lo=1e-3, a_hi=25){
  target <- pmin(pmax(target, 1e-10), 1 - 1e-10)
  f <- function(a) .au_m1_stable(a) - target
  out <- try(uniroot(f, c(a_lo, a_hi))$root, silent=TRUE)
  if (inherits(out, "try-error")) {
    aa <- exp(seq(log(a_lo+1e-6), log(a_hi), length.out=200))
    aa[which.min(abs(.au_m1_stable(aa) - target))]
  } else out
}
.kernel_w <- function(d, bw, kernel=c("bisquare","gaussian")){
  kernel <- match.arg(kernel)
  if (kernel == "bisquare"){
    w <- (1 - (d/bw)^2); w[d>bw] <- 0; pmax(w,0)^2
  } else {
    exp(-(d^2)/(2*bw^2))
  }
}
.dist2d <- function(A, B){
  dx <- outer(A[,1], B[,1], "-")
  dy <- outer(A[,2], B[,2], "-")
  sqrt(dx^2 + dy^2)
}

# --- helpers para ν (extraer coeficientes y alinear diseño) ---
.get_coef_nu <- function(obj){
  if (!is.null(obj$coef_nu)) return(obj$coef_nu)
  if (!is.null(obj$coef.nu)) return(obj$coef.nu)
  if (!is.null(obj$nu_hat))  return(NULL)  # ν constante
  tries <- list(
    try(obj$nu$coef, silent=TRUE),
    try(coef(obj$nu), silent=TRUE),
    try(coef(obj$fit_nu), silent=TRUE),
    try(coef(obj$mod_nu), silent=TRUE),
    try(obj$gamma, silent=TRUE),
    try(obj$coefNu, silent=TRUE)
  )
  for (x in tries) if (!inherits(x, "try-error") && !is.null(x)) return(x)
  if (!is.null(names(obj))){
    for (nm in names(obj)){
      xo <- obj[[nm]]
      if (is.numeric(xo) && length(xo) >= 2) return(xo)
      if (inherits(xo, "glm")) return(coef(xo))
      if (is.list(xo) && !is.null(xo$coef)) return(xo$coef)
    }
  }
  NULL
}
.prep_nu_hat <- function(coef_nu, form_nu, data_frame){
  Z <- model.matrix(form_nu, data_frame)
  b <- as.numeric(coef_nu); nm_b <- names(coef_nu)
  if (is.null(nm_b) || any(nm_b=="")){
    if (length(b) == ncol(Z)) nm_b <- colnames(Z) else
      stop("coef_nu sin nombres y longitud distinta a ncol(Z).")
  }
  names(b) <- nm_b
  if ("(Intercept)" %in% names(b) && !("(Intercept)" %in% colnames(Z)))
    Z <- cbind("(Intercept)" = 1, Z)
  miss <- setdiff(names(b), colnames(Z))
  if (length(miss)){
    for (mc in miss){
      Z <- cbind(Z, setNames(matrix(if (mc=="(Intercept)") 1 else 0,
                                    nrow(Z), 1), mc))
    }
  }
  Z <- Z[, names(b), drop=FALSE]
  drop(plogis(Z %*% b))
}

# --- objetivo local (ZIAU, c=0): -sum w_i * log L_i(beta) ---
.make_local_nll <- function(link_a="log"){
  function(beta, X, y, nu, w){
    eta   <- as.vector(X %*% beta)
    eta   <- pmin(pmax(eta, -30), 30)
    alpha <- if (link_a=="log")   exp(eta) else
      if (link_a=="logit") plogis(eta) else
        if (link_a=="identity") pmax(eta,1e-8) else exp(eta)
    is_zero <- (y <= 0); is_pos <- (y > 0) & (y < 1)
    ll <- numeric(length(y))
    ll[is_zero] <- log(nu[is_zero])
    if (any(is_pos)){
      lpdf <- .au_logpdf(y[is_pos], alpha[is_pos])
      ll[is_pos] <- log1p(-nu[is_pos]) + lpdf
    }
    if (any(!is_zero & !is_pos)){   # bordes (p.ej., y≈1)
      yy <- pmin(pmax(y[!is_zero & !is_pos], 1e-12), 1 - 1e-12)
      lpdf <- .au_logpdf(yy, alpha[!is_zero & !is_pos])
      ll[!is_zero & !is_pos] <- log1p(-nu[!is_zero & !is_pos]) + lpdf
    }
    nll <- -sum(w * ll); if (!is.finite(nll)) nll <- 1e12; nll
  }
}

# ===============================================================
predict.ggwrm <- function(object,
                          newdata,
                          newcoords = NULL,
                          type = c("mean","alpha","eta","prob_zero"),
                          c_inflate = 0,
                          return_coefs = FALSE,
                          metrics = TRUE,
                          ...){
  type <- match.arg(type)
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  
  # --- recuperar hiperparámetros y fórmulas del objeto ---
  GA       <- object$GW.arguments %||% list()
  form_mu  <- attr(object, "form_alpha") %||% GA$formula
  form_nu  <- attr(object, "form_nu")    %||% GA$formula.nu
  if (is.null(form_mu)) stop("No encuentro la fórmula de alpha (form_mu).")
  if (is.null(form_nu)) stop("No encuentro la fórmula de nu (form_nu).")
  
  kernel   <- GA$kernel %||% "bisquare"
  bw       <- GA$bw     %||% stop("No encuentro 'bw' en el objeto.")
  link_a   <- GA$link.alpha %||% "log"
  
  # --- TRAIN: SDF o atributo guardado ---
  sp_train <- attr(object, "sp_train")
  if (is.null(sp_train)) {
    if (inherits(object$SDF, "SpatialPointsDataFrame")) {
      sp_train <- object$SDF
    } else stop("No encuentro el conjunto de entrenamiento (sp_train/SDF) dentro del objeto.")
  }
  
  # TRAIN: y, X, coords
  X_tr <- model.matrix(form_mu, sp_train@data)
  resp <- all.vars(form_mu)[1]
  if (!resp %in% names(sp_train@data)) stop("La variable respuesta '", resp,"' no está en TRAIN.")
  y_tr <- sp_train@data[[resp]]
  ctr  <- coordinates(sp_train)
  
  # --- TEST: SpatialPointsDataFrame o data.frame + coords ---
  if (inherits(newdata, "SpatialPointsDataFrame")) {
    sp_test <- newdata
  } else if (is.data.frame(newdata)) {
    if (is.null(newcoords)) {
      if (!all(c("LON","LAT") %in% names(newdata)))
        stop("Si 'newdata' es data.frame, pasa 'newcoords' o incluye columnas 'LON' y 'LAT'.")
      newcoords <- as.matrix(newdata[, c("LON","LAT")])
    }
    if (!resp %in% names(newdata)) newdata[[resp]] <- NA_real_
    sp_test <- sp::SpatialPointsDataFrame(
      coords = data.frame(LON=newcoords[,1], LAT=newcoords[,2]),
      data   = newdata,
      proj4string = sp_train@proj4string
    )
  } else stop("'newdata' debe ser SpatialPointsDataFrame o data.frame.")
  
  # TEST: matrices
  X_te <- model.matrix(form_mu, sp_test@data)
  cte  <- coordinates(sp_test)
  y_te <- sp_test@data[[resp]]
  
  # --- modelo global de ν (sin nu.cut) ---
  coef_nu  <- .get_coef_nu(object)
  nu_const <- object$nu_hat %||% NULL
  
  if (!is.null(coef_nu)) {
    nu_tr <- try(.prep_nu_hat(coef_nu, form_nu, sp_train@data), silent=TRUE)
    if (inherits(nu_tr, "try-error")) {
      dat_tr <- sp_train@data; dat_tr$..z0.. <- as.integer(dat_tr[[resp]] <= 0)
      mod_nu <- glm(update(form_nu, ..z0.. ~ .), data = dat_tr, family = binomial())
      coef_nu <- coef(mod_nu)
      nu_tr   <- .prep_nu_hat(coef_nu, form_nu, sp_train@data)
    }
  } else if (!is.null(nu_const)) {
    nu_tr <- rep(nu_const, nrow(X_tr))
  } else {
    dat_tr <- sp_train@data; dat_tr$..z0.. <- as.integer(dat_tr[[resp]] <= 0)
    mod_nu <- glm(update(form_nu, ..z0.. ~ .), data = dat_tr, family = binomial())
    coef_nu <- coef(mod_nu)
    nu_tr   <- .prep_nu_hat(coef_nu, form_nu, sp_train@data)
  }
  nu_tr <- pmin(pmax(nu_tr, 1e-12), 1 - 1e-12)   # clip estable (sin nu.cut)
  
  # --- MLE local por punto de TEST ---
  D     <- .dist2d(ctr, cte)               # n_train x n_test
  p     <- ncol(X_tr)
  B_te  <- matrix(NA_real_, nrow = nrow(X_te), ncol = p)
  colnames(B_te) <- colnames(X_tr)
  .local_nll <- .make_local_nll(link_a)
  
  for (j in seq_len(nrow(X_te))){
    d  <- D[, j]
    wj <- .kernel_w(d, bw, kernel)
    if (sum(wj) == 0) wj[which.min(d)] <- 1
    
    # warm-start: alpha0 desde inversión de m_AU de la media local
    num <- sum(wj * y_tr, na.rm=TRUE)
    den <- sum(wj * (1 - nu_tr), na.rm=TRUE)
    t_j <- if (den > 0 && is.finite(num/den)) num / den else mean(y_tr[y_tr > 0], na.rm=TRUE)
    t_j <- ifelse(is.finite(t_j), pmin(pmax(t_j, 1e-6), 1 - 1e-6), 0.1)
    alpha0 <- .au_m1_inv(t_j)
    beta0  <- rep(0, p)
    beta0[1] <- if (link_a=="log") log(alpha0) else
      if (link_a=="logit") qlogis(alpha0) else alpha0
    
    fit_j <- try(
      optim(par = beta0, fn = .local_nll,
            X = X_tr, y = y_tr, nu = nu_tr, w = wj,
            method = "BFGS", control = list(maxit = 300, reltol = 1e-8)),
      silent = TRUE
    )
    if (inherits(fit_j, "try-error") || !is.finite(fit_j$value)){
      fit_j <- try(
        optim(par = beta0, fn = .local_nll,
              X = X_tr, y = y_tr, nu = nu_tr, w = wj,
              method = "Nelder-Mead", control = list(maxit = 600, reltol = 1e-8)),
        silent = TRUE
      )
    }
    B_te[j, ] <- if (inherits(fit_j, "try-error")) beta0 else fit_j$par
  }
  
  # --- predicción en TEST ---
  eta_te   <- rowSums(X_te * B_te)
  eta_te   <- pmin(pmax(eta_te, -30), 30)
  alpha_te <- if (link_a=="log")   exp(eta_te) else
    if (link_a=="logit") plogis(eta_te) else
      if (link_a=="identity") pmax(eta_te,1e-8) else exp(eta_te)
  
  if (!is.null(coef_nu)) {
    nu_te <- .prep_nu_hat(coef_nu, form_nu, sp_test@data)
  } else {
    nu_te <- rep(nu_const, nrow(X_te))
  }
  nu_te <- pmin(pmax(nu_te, 1e-12), 1 - 1e-12)   # sin umbralizar
  
  mu_AU <- .au_m1_stable(alpha_te)
  
  pred <- switch(type,
                 "mean"      = if (c_inflate==0) (1 - nu_te) * mu_AU else (1 - nu_te) * mu_AU + nu_te * 1,
                 "alpha"     = alpha_te,
                 "eta"       = eta_te,
                 "prob_zero" = if (c_inflate==0) nu_te else 0
  )
  
  # >>> incluye y_obs siempre (si no existe en newdata, quedará NA)
  out <- data.frame(
    y_obs     = y_te,
    pred      = pred,
    alpha_hat = alpha_te,
    eta_alpha = eta_te,
    nu_hat    = nu_te
  )
  
  # métricas (solo si newdata trae y observada)
  if (metrics && !all(is.na(y_te))) {
    ok <- is.finite(y_te) & is.finite(out$pred)
    yobs <- y_te[ok]; yprd <- out$pred[ok]
    MSE  <- mean((yobs - yprd)^2)
    MAE  <- mean(abs(yobs - yprd))
    RMSE <- sqrt(MSE)
    attr(out, "metrics") <- data.frame(MSE=MSE, MAE=MAE, RMSE=RMSE)
  }
  
  if (return_coefs) attr(out, "betas_locales") <- B_te
  out
}

# ---- registrar método S3 ----
registerS3method("predict", "ggwrm", predict.ggwrm)


# -------- registrar métodos S3 (por si el script se 'sourcea' a mitad) -----
registerS3method("print",       "ggwrm", print.ggwrm)
registerS3method("as.data.frame","ggwrm", as.data.frame.ggwrm)
registerS3method("plot",        "ggwrm", plot.ggwrm)

#--------------------------------------------------------#
#------------MODELS (ESPECIFICATIONS)--------------------#
# GWR-ZIAU  #1 NU= FULL PCA; Alpha = LAT,LON | Link: Logit (nu); log(alpha)
# GWR-ZIAU  #2 NU= LAT, LON; Alpha = FULL PCA| Link: logit (nu); log(alpha)
# GWR-ZIAU  #4 NU= LAT, LON; Alpha = FULL PCA| Link: probit (nu); log(alpha)
# GWR-ZIAU  #5 NU= FULL PCA; Alpha = FULL PCA| Link: logit (nu); log(alpha)

#------------------------------------------------------------------------------#
# TRAINING MODELO 
#-----------------------------------------------------------------------------#
# CREACION SPATIALPOINTS DATAFRAME

#GWR-ZIAU FULL #1 
spCob <- SpatialPointsDataFrame(
  coords      = data.frame(LON = DB_train$LON, LAT = DB_train$LAT),
  data        = data.frame(
    y     = DB_train$V1,
    LON  = DB_train$LON,
    LAT = DB_train$LAT,
    DB_train[, paste0("PC", 1:40)]
  ),
  proj4string = CRS("")
)

# GWR-ZIAU FULL #2| GWR-ZIAU  #3|  GWR-ZIAU  #4| 
spCob2 <- SpatialPointsDataFrame(
  coords = data.frame(LON = DB_train$LON, LAT = DB_train$LAT),
  data   = data.frame(
    y    = DB_train$V1,
    LON = DB_train$LON,
    LAT = DB_train$LAT,
    DB_train[, paste0("PC", 1:40), drop = FALSE]  
  ),
  proj4string = CRS("")   
)

###############################################################################
# FORMULAS
# ALPHA
form_alpha_coordenadas <- y ~ LON + LAT
form_alpha_FULL_PCA <- as.formula(paste("y ~", paste(paste0("PC", 1:40), collapse = " + ")))

# NU
form_nu_coordenadas <- ~ LON + LAT
form_nu_FULL_PCA    <- as.formula(paste("~", paste(paste0("PC", 1:40), collapse = " + ")))


################################################################################
# SPATIAL GRID
# GWR-ZIAU.1 
set.seed(123)
bw_cob1 <- bw.ggwr.ziau(
  formula    = form_alpha_coordenadas,
  data       = spCob,
  approach   = "CV",
  kernel     = "gaussian", 
  adaptive   = FALSE,
  link.alpha = "log",
  link.nu    = "logit",
  formula.nu = form_nu_FULL_PCA,
  c_inflate  = 0
)

# GWR-ZIAU.2
set.seed(123)
bw_cob3 <- bw.ggwr.ziau(
  formula    = form_alpha_FULL_PCA,
  data       = spCob2,
  approach   = "CV",
  kernel     = "gaussian",
  adaptive   = FALSE,
  link.alpha = "log",
  link.nu    = "logit",
  formula.nu = form_nu_coordenadas,   
  c_inflate  = 0
)

# GWR-ZIAU.3
set.seed(123)
bw_cob4 <- bw.ggwr.ziau(
  formula    = form_alpha_FULL_PCA,
  data       = spCob2,
  approach   = "CV",
  kernel     = "gaussian",
  adaptive   = FALSE,
  link.alpha = "log",
  link.nu    = "probit",
  formula.nu = form_nu_coordenadas,
  c_inflate  = 0
)

# GWR-ZIAU.4
set.seed(123)
bw_cob5 <- bw.ggwr.ziau(
  formula    = form_alpha_FULL_PCA,
  data       = spCob2,
  approach   = "CV",
  kernel     = "gaussian",
  adaptive   = FALSE,
  link.alpha = "log",
  link.nu    = "logit",
  formula.nu = form_nu_FULL_PCA,
  c_inflate  = 0
)

###############################################################################
## FITTING MODELS
# GWR-ZIAU.1 
fit_cob <- ggwr.ziau.basic(
  formula     = form_alpha_coordenadas,
  data        = spCob,
  bw          = bw_cob1,
  kernel      = "gaussian",
  adaptive    = FALSE,
  cv          = TRUE,
  link.alpha  = "log",
  link.nu     = "logit", 
  formula.nu  = form_nu_FULL_PCA,
  nu.cut      = 0.4676828,
  c_inflate   = 0
)
print(fit_cob)

# GWR-ZIAU.2 
fit_cob3 <- ggwr.ziau.basic(
  formula     = form_alpha_FULL_PCA,
  data        = spCob2,
  bw          = bw_cob3,
  kernel      = "gaussian",
  adaptive    = FALSE,
  cv          = TRUE,
  link.alpha  = "log",
  link.nu     = "logit",
  formula.nu  = form_nu_coordenadas,
  nu.cut      = 0.4676828,
  c_inflate   = 0
)
print(fit_cob3)

# GWR-ZIAU.3
fit_cob4 <- ggwr.ziau.basic(
  formula     = form_alpha_FULL_PCA,
  data        = spCob2,
  bw          = bw_cob4,
  kernel      = "gaussian",
  adaptive    = FALSE,
  cv          = TRUE,
  link.alpha  = "log",
  link.nu     = "probit",
  formula.nu  = form_nu_coordenadas,
  nu.cut      = 0.3534911,   
  c_inflate   = 0
)
print(fit_cob4)

# GWR-ZIAU.4
fit_cob5 <- ggwr.ziau.basic(
  formula     = form_alpha_FULL_PCA,
  data        = spCob2,
  bw          = bw_cob5,
  kernel      = "gaussian",
  adaptive    = FALSE,
  cv          = TRUE,
  link.alpha  = "log",
  link.nu     = "logit",
  formula.nu  = form_nu_FULL_PCA,
  nu.cut      = 0.4676828,   
  c_inflate   = 0
)
print(fit_cob5)