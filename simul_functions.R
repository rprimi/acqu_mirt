############### sim_trait_acq

sim_trait_acq <- function(
      
      N=10000,   # number of subjects
      mean_a1,   # mean of the discrimination distribution
      sd_a1,     # standard deviation of the discrimination distribution
      sd_acq,   # variance of the acquiescence factor
      n_pw,      # number of positively worded items
      n_nw,      # number of negatively worded items
      interc_pw,  # matrix with positively worded items itntercepts to draw from
      interc_nw  # matrix with negatively worded items itntercepts to draw from
  )      {
    
    # rescales discrimination parameters to lognorm 
    log_m <- log(mean_a1^2 / sqrt(sd_a1^2 + mean_a1^2))
    log_sd <- sqrt(log(1 + (sd_a1^2 / mean_a1^2)))
    
    # discrimination parameters for trait and acq factors
    n_it = n_pw + n_nw
    
    a1 <-rlnorm(n_it, meanlog = log_m, sdlog = log_sd)
    
    # fixed weights for balanced scales
    a2 <- c(rep(1, n_pw), rep(-1, n_nw))
    
    # fixed weights that compensates for unbalanced scales
    w = n_it / 2
    w_p = w / n_pw
    w_n = w / n_nw
    a3 <- c(rep(w_p, n_pw), rep(-w_n, n_nw))
    
    
    A <- cbind(a1, a2, a3)
    
    # intercepts
    interc_pw <-  interc_pw %>% sample_n(n_pw)
    interc_nw <-  interc_nw %>% sample_n(n_nw)
    D <- rbind(interc_pw, interc_nw)
    
    # generates theta (note standard deviation of acq is relatedt to its discrimination)
    
    trait_t <- rnorm(N, mean = 0, sd = 1)
    acq_t <- rnorm(N, mean = 0, sd = sd_acq)
    theta <- as.matrix(cbind(trait_t,  acq_t))
    
    D<-as.matrix(D)
    dimnames(D) <- NULL 
    
    # theta <- as.matrix(theta)
    # dimnames(theta) <- NULL 
    
    dt <- simdata(a=A[, 1:2], d=D, N=N, itemtype="graded", Theta = theta)
    
    items = cbind(A, D) %>% as.data.frame
    items$coditem <-dimnames(dt)[[2]]
    
    
    return( list(
        data = dt,
        persons = theta,
        items = items,
        acq_sd = sd_acq
    )
    )
  
}

##############  classical_scores

classical_scores <- function(obj, lik = 5){
    
    # data was generated reversed
    scr =  apply(obj$data, MARGIN = 1, mean, na.rm=TRUE)
    
    # reverse back items to calculate acq score
    dt_tmp =  obj$data %>% 
        as.data.frame %>% 
        mutate_if(obj$items[, "a2"] ==-1, funs((lik+1)-.))
    
    acq = apply(dt_tmp, MARGIN = 1, mean, na.rm=TRUE)
    
    acq2 = apply(
        dt_tmp, MARGIN = 1, function(x) { 
        as.numeric(t(x)) %*% as.numeric(abs(obj$items[, "a3"])) / length(x) 
        }
      )
    
    # compute scores corrected for aaquiescence
    
    scr_rec <- dt_tmp %>%  apply(MARGIN = 2, function(x) { (x - acq)} ) %>% as.data.frame %>%
        mutate_if(obj$items[, "a2"] ==-1, funs(-1*.)) %>% 
        apply(MARGIN = 1, mean, na.rm=TRUE)
    
    scr_rec2 <- dt_tmp %>%  apply(MARGIN = 2, function(x) { (x - acq2)} ) %>% as.data.frame %>%
        mutate_if(obj$items[, "a2"] ==-1, funs(-1*.)) %>%
        apply(MARGIN = 1, mean, na.rm=TRUE)
    
    return(cbind(scr, scr_rec2, scr_rec,acq, acq2))
    
}


############## PDIF scores



pdif_scores <- function(obj) {
    
    # 1. Calibrate parameters with all items
    
    grsm <- obj$data %>% 
      mirt(
        1, 
        itemtype = 'graded',
        TOL = .001
      )
    
    scr_1d <-  fscores(
      grsm, 
      method = "EAP", 
      full.scores.SE = TRUE
      )
    
    # 2. Recover item parameter
    
    prms0 <- mod2values(grsm)
    
    
    # 3. Scores from negatively phrased items 
    
    v <- obj$items %>% filter(a2 ==-1)
    dt_nw <- obj$data  %>% 
      as.data.frame %>% 
      select(v$coditem) 
    
    # 3.1 recover df of parameters from the subeset
    prms_nw <- dt_nw %>%
      mirt(
        1, 
        itemtype = 'graded',
        pars = 'values',
        TOL = .001
        )
    
    # 3.2 get corerresponding parameters
    calibrated <- prms0 %>% 
      filter(item %in% c(v$coditem, "GROUP"))
    
    # 3.3 replace parameters in df with calibrated parameters and fix it
    
    prms_nw$value <-  calibrated$value
    prms_nw$est <-FALSE
    
    # 3.4 calibrate q\with fixed parameters
    
    calibr_nw <- 
      dt_nw %>% 
      mirt(
        1, 
        itemtype = 'graded',
        pars = prms_nw,
        TOL = .001
        )
    
    # 3.5 scores negative worded
    scr_nw <-  fscores(
      calibr_nw , 
      method = "EAP", 
      full.scores.SE = TRUE
      )
    
    
    # 4. Scores from positively phrased items 
    
    v <- obj$items %>% filter(a2 ==1)
    dt_pw <- obj$data  %>% as.data.frame %>% select(v$coditem) 
    
    # 4.1 recover df of parameters from the subeset
    prms_pw <- dt_pw %>%
      mirt(
        1, 
        itemtype = 'graded',
        pars = 'values',
        TOL = .001
        )
    
    # 4.2 get corerresponding parameters
    calibrated <- prms0 %>% 
      filter(item %in% c(v$coditem, "GROUP"))
    
    # 4.3 replace parameters in df with calibrated parameters and fix it
    
    prms_pw$value <-  calibrated$value
    prms_pw$est <-FALSE
    
    # 4.4 calibrate q\with fixed parameters
    
    calibr_pw <-dt_pw %>% 
      mirt(
        1, 
        itemtype = 'graded',
        pars = prms_pw,
        TOL = .001
        )
    
    # 4.5 scores negative worded
    scr_pw <-  fscores(calibr_pw, method = "EAP", full.scores.SE = TRUE)
    
    # Calculates PDIF
    scrs <- cbind(
      f1_1d = scr_1d[ , 1],
      f1_pos = scr_pw[ ,1], 
      f1_neg = scr_nw[ ,1], 
      pdif_mirt = (scr_pw[ ,1] - scr_nw[ ,1]),
      f1_1d_se = scr_1d[ , 1],
      f1_pos_se = scr_pw[ ,2], 
      f1_neg_se = scr_nw[ ,2]
      ) %>% 
      as.data.frame()
    
    return(
      list(
        scrs = scrs,
        mirt_obj = grsm)
    )
    
}


############## RI scores


ri_scores <- function(obj, a = "a2", cov_ta = TRUE) {
    
    # recover df with parameters
    prms  <-obj$data %>%
      mirt(
        2, 
        itemtype = 'graded',
        pars = 'values',
        TOL = .001
        )
    
    # fix acq factor discriminations and estimate variance of f2
    
    prms[prms$name=="a2", ]$value <- obj$items[ , a] 
    prms[prms$name=="a2", ]$est <- FALSE 
    prms[prms$name=="COV_22", ]$est <- TRUE 
    prms[prms$name=="COV_21", ]$est <- cov_ta 
    prms[prms$name=="a1", ]$value <- abs(prms[prms$name=="a1", ]$value)
    
    # calibrate item parameters
    ri <- obj$data %>% 
      mirt(
        2, 
        pars = prms, 
        itemtype = 'graded', 
        TOL = .001
        )
    
    scrs <- fscores(
      ri, 
      method = "EAP", 
      full.scores.SE = TRUE
      )
    
    return(
      list(
        scrs = scrs,
        mirt_obj = ri
      )
    )
  
}