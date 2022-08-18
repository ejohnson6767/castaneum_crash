library(tidyverse)
library(gtools)
library(doParallel)
library(doRNG)
library(foreach)
library(doFuture)

# transform parameters between constrained scales (good for parameter interpretation) and unconstrained scales (good for optimization)
lower_constrain <- function(uncon_var, lower = 0) {exp(uncon_var) + lower}
lower_unconstrain <- function(con_var, lower = 0) {log(con_var - lower)}
lower_upper_unconstrain <- function(con_var, lower = 0,  upper = 1) {
  if(lower >= upper) stop("lower constraint is >= upper constraint")
  logit((con_var - lower)/(upper - lower))}
lower_upper_constrain <- function(uncon_var, lower = 0,  upper = 1) {
  if(lower >= upper) stop("lower constraint is >= upper constraint")
  lower + (upper - lower) * inv.logit(x = uncon_var)}


#lower_upper_unconstrain(0.9999999999, -1, 1)
lower_upper_constrain(Inf, -5, 100)
# constrain_pars <- function(named_pars)
# {
#   con_par_list <- list()
#   for(i in 1:length(named_pars))
#   {
#     if(names(named_pars)[[i]] == "high_zt_survival_prob_uncon")
#     {
#       con_par_list[[i]] <- lower_upper_constrain(named_pars[[i]])
#       names(con_par_list)[[i]] <- sub(pattern = "_uncon", replacement = "", x = names(named_pars)[[i]])
#     } else if(names(named_pars)[[i]] == "mu_z_star1_uncon")
#     {
#       con_par_list[[i]] <- named_pars[[i]]
#       names(con_par_list)[[i]] <- sub(pattern = "_uncon", replacement = "", x = names(named_pars)[[i]])
#     } else{
#       con_par_list[[i]] <- lower_constrain(named_pars[[i]])
#       names(con_par_list)[[i]] <- sub(pattern = "_uncon", replacement = "", x = names(named_pars)[[i]])
#       
#     }
#   }
#   return(con_par_list)
# }


constrain_pars <- function(named_pars)
{
  con_par_list <- list()
  for(i in 1:length(named_pars))
  {
    if(names(named_pars)[[i]] == "high_zt_survival_prob_uncon")
    {
      con_par_list[[i]] <- lower_upper_constrain(named_pars[[i]])
      names(con_par_list)[[i]] <- sub(pattern = "_uncon", replacement = "", x = names(named_pars)[[i]])
    } else{
      con_par_list[[i]] <- lower_constrain(named_pars[[i]])
      names(con_par_list)[[i]] <- sub(pattern = "_uncon", replacement = "", x = names(named_pars)[[i]])
      
    }
  }
  return(con_par_list)
}


unconstrain_pars <- function(named_pars)
{
  uncon_par_list <- list()
  for(i in 1:length(named_pars))
  {
    if(names(named_pars)[[i]] == "high_zt_survival_prob")
    {
      uncon_par_list[[i]] <- lower_upper_unconstrain(named_pars[[i]])
      names(uncon_par_list)[[i]] <- paste0(names(named_pars)[[i]],"_uncon")
    } else if(names(named_pars)[[i]] == "mu_z_star1_uncon")
    {
      con_par_list[[i]] <- named_pars[[i]]
      names(con_par_list)[[i]] <- sub(pattern = "_uncon", replacement = "", x = names(named_pars)[[i]])
    } else
    {
      uncon_par_list[[i]] <- lower_unconstrain(named_pars[[i]])
      names(uncon_par_list)[[i]] <- paste0(names(named_pars)[[i]],"_uncon")
      
    }
  }
  return(uncon_par_list)
}

# grepl( "mu_alpha", "mu_alpha_uncon", fixed = TRUE)
# nameVar <- "a name"
# sub(pattern = "_uncon", replacement = "", x = "mu_alpha_uncon")
# l <- list(69)
# names(l[[1]]) <- nameVar
# l
# v <- c(1,2,3)
# v[[1]]

range <- c(0, Inf)
m = 1; s =1
n <- 10

rnormt <- function(n, range, m, s) {
  a <- pnorm(min(range), mean = m, sd = s)
  b <- pnorm(max(range), mean = m, sd = s)
  u <- runif(n, a, b)
  res <- qnorm(u, mean = m, sd = s)
  res[which(!is.finite(res))] <- 0
  res
}

try_optim <- function(nll_fun, par_init_fun, dat, max_tries = 10, ...)
{
  j <- 1
  while(j <= max_tries)
  {
    tryCatch({
      fit <- optim(nll_fun, par = par_init_fun(), dat = dat, ...)
      return(fit)
    },error=function(e){
      print(e$message)
      j <<- j + 1
    },finally={})
  }
  stop(paste0(max_tries, " / ", max_tries, " attempts ended with error. Try altering the initalization function"))
}

multi_mle <- function(nll_fun, par_init_fun, dat, max_mle_attempts = 5, max_tries = 1000, max_non_convergence = 5, num_cores = 1, ...)
{
  non_convergence_counter <- 0
  fit_list <- list()
  for(i in 1:max_mle_attempts)
  {
    tryCatch({
      cur_fit <- try_optim(nll_fun, par_init_fun, dat, max_tries,...)
      fit_list <- append(fit_list, list(cur_fit))
      if(cur_fit$convergence != 0){
        non_convergence_counter <- non_convergence_counter + 1
      }
    },error=function(e){
      stop(paste0("On iter ", i, " / ", max_mle_attempts, " mle attempts : ", e))
    },finally={})
  }
  
  if(non_convergence_counter != 0)
  {
    warning(paste0(non_convergence_counter, " / ", max_mle_attempts, " mle attempts initiated but did not converge before the maximum number of optim iterations \n "))
  }
  #parallel::stopCluster(cl)
  return(fit_list)
}


multi_mle_parallel <- function(nll_fun, par_init_fun, dat, max_mle_attempts = 5, max_tries = 1000, max_non_convergence = 5, num_cores = 1, ...)
{
  if(num_cores == 1)
  {
    non_convergence_counter <- 0
    fit_list <- list()
    for(i in 1:max_mle_attempts)
    {
      tryCatch({
        cur_fit <- try_optim(nll_fun, par_init_fun, dat, max_tries,...)
        fit_list <- append(fit_list, list(cur_fit))
        if(cur_fit$convergence != 0){
          non_convergence_counter <- non_convergence_counter + 1
        }
      },error=function(e){
        stop(paste0("On iter ", i, " / ", max_mle_attempts, " mle attempts : ", e))
      },finally={})
    }
    
    if(non_convergence_counter != 0)
    {
      warning(paste0(non_convergence_counter, " / ", max_mle_attempts, " mle attempts initiated but did not converge before the maximum number of optim iterations \n "))
    }
    #parallel::stopCluster(cl)
    return(fit_list)
  } else
  {
    model_name <- as.character(substitute(nll_fun))
    model_num <- substr(model_name,4,5)
    if(is.na(as.numeric(model_num)))
    {
      model_num <- substr(model_name,4,4)
    }
    registerDoFuture()
    cl <- makeCluster(num_cores)
    plan(cluster, workers=cl)
    
    non_convergence_counter <- 0
    fit_list <- foreach(i = 1:max_mle_attempts,
                        .noexport = c(model_name),
                        .packages = c('Rcpp', 'here')) %dorng%
      {
        if(any( c(1:6) == model_num))  {
          Rcpp::sourceCpp(here("beetle_models_sourceCpp.cpp"))
        } else if(any( c(7:12) == model_num))  {
          Rcpp::sourceCpp(here("beetle_models_dd_cannibalism_sourceCpp.cpp"))
        } else if(any( c(13:16) == model_num))  {
          Rcpp::sourceCpp(here("beetle_models_no_crash_sourceCpp.cpp"))
        } else if(any( c(17:22) == model_num))  {
          Rcpp::sourceCpp(here("beetle_models_dd_cannibalism_dd_threshold_sourceCpp.cpp"))
        }
        model_name <- paste0("mod", model_num, "LlCpp")
        
        cur_fit <- try_optim(get(model_name), par_init_fun, dat, max_tries,...)
        cur_fit
      }
    
    stopCluster(cl)
    return(fit_list)
  }
}


robust_mle <- function(nll_fun, par_init_fun, dat, max_mle_attempts = 5, max_tries = 1000, ...)
{
  fit_list <- multi_mle(nll_fun, par_init_fun,dat, max_mle_attempts, max_tries, ...)
  fit_vals <- map_dbl(fit_list, function(x){x[["value"]]})
  return(fit_list[[first(which(fit_vals == min(fit_vals)))]])
}

select_mle <- function(fit_list)
{
  fit_vals <- map_dbl(fit_list, function(x){x[["value"]]})
  return(fit_list[[first(which(fit_vals == min(fit_vals)))]])
}

# dat <- dat_A_crash
# mod_par_test_fun <- function(){c(mu_alpha_uncon = lower_unconstrain(1),
#                                  beta0_uncon = lower_unconstrain(0.05),
#                                  eta_alpha_uncon = lower_unconstrain(20),
#                                  high_zt_survival_prob_uncon = lower_upper_unconstrain(0.05),
#                                  thresh_zt_range_uncon = lower_unconstrain(30),
#                                  mu_z_star_uncon = lower_unconstrain(500),
#                                  sigma_z_star_uncon = lower_unconstrain(10))}


# k_fold_CV <- function(model_list, par_init_fun, dat, max_mle_attempts = 5, max_tries = 1000, k = 2, delta = 4, M =20, ztStep = 1, ftStep = 1, lowZtSurvivalProb = 0.95, ...){
#   # print warnings in real time, messages automatically print in real time
#   options(warn = 1)
#   # prep return data
#   num_models <- length(model_list)
#   dat_len <- k * num_models
#   #cv_return_dat <- data.frame(model_name = as.character(dat_len), fold = as.integer(), log_lik = as.numeric, convergence = as.integer())
#   cv_return_dat <- data.frame(matrix(NA, nrow = dat_len, ncol = 5))
#   colnames(cv_return_dat) <- c("model_name", "fold", "log_lik", "convergence", "init_error", "minutes", "iterations")
#   
#   # combine general (i.e. cross-species) control parameters (for optim) with species-species control parameter(i.e. parscale)
#   for(j in 1:num_models)
#   {
#     if("control" %in% names(list(...)))
#     {
#       (model_list[[j]])[["control"]] <- c((list(...))[["control"]], (model_list[[j]])[["parscale"]])
#     } else
#     {
#       (model_list[[j]])[["control"]] <- (model_list[[j]])[["parscale"]]
#     }
#     control
#   }
#   
#   # take "control" away from the additional arguments
#   additional_args <- (list(...))[which(names(list(...)) != "control")]
#   
#   # split data into k subsets
#   dat$fold <- sample(rep( 1:k, length.out = NROW(dat)))
#   
#   # loop over each fold
#   for(i in 1:k)
#   {
#     train_dat <- dat %>% filter(fold != i)
#     valid_dat <- dat %>% filter(fold == i)
#     
#     # put data into list (preferred input of cpp log lik functions)
#     train_dat_list <- list(N = NROW(train_dat),
#                            nt = train_dat$nt,
#                            ntp1 = train_dat$ntp1,
#                            oviPeriod = train_dat$oviposition_days,
#                            delta = delta,
#                            M = M,
#                            ztStep = ztStep,
#                            ftStep = ftStep,
#                            lowZtSurvivalProb = lowZtSurvivalProb)
#     
#     
#     valid_dat_list <- list(N = NROW(valid_dat),
#                            nt = valid_dat$nt,
#                            ntp1 = valid_dat$ntp1,
#                            oviPeriod = valid_dat$oviposition_days,
#                            delta = delta,
#                            M = M,
#                            ztStep = ztStep,
#                            ftStep = ftStep,
#                            lowZtSurvivalProb = lowZtSurvivalProb)
#     
#     # loop over each model
#     for(j in 1:num_models)
#     {
#       # print progress
#       cat("fold ", i, ", model ", j, " : ~ ", round(100*(((i-1)*num_models + j)/dat_len)), "% complete. \n", sep = "")
#       
#       # return a data frame row
#       cv_return_dat[(i-1)*num_models + j,] <-
#         withCallingHandlers({
#           withRestarts({
#             cur_mle <- robust_mle(model_list[[j]], par_init_fun, train_dat_list, max_mle_attempts, max_tries, additional_args ) # get mle on training data
#             valid_ll <- model_list[[j]](cur_mle$par, valid_dat_list) # get LL on validation data
#             list(model_name = names(model_list)[[j]], 
#               fold = i, 
#               log_lik = valid_ll, 
#               convergence = ifelse(cur_mle$convergence == 0, 1, 0),
#               init_error = 0) # return data
#           },
#           muffle_stop = function(){
#             list(model_name = names(model_list)[[j]], 
#               fold = i, 
#               log_lik = NA, 
#               convergence = NA,
#               init_error = 1) # return row anyways if optim function can't initialize
#           }
#           )
#         },
#         warning=function(w){
#           warning(paste0("With model = ", names(model_list)[[j]], ", fold = ", k, " : ", w$message))
#           rlang::cnd_muffle(w) # prevents warning from printing again on function exit
#         },
#         error=function(e){
#           message(paste0("Muffled Error : With model = ", names(model_list)[[j]], ", fold = ", k, " : ", e$message))
#           invokeRestart("muffle_stop") # skips the the the withRestarts block code; goes to muffle_stop
#         })
#       
#     }
#   }
#   
#   # retrieve data frame with model_name, and log lik
#   return(cv_return_dat)
# }

k_fold_CV <- function(model_list, dat, max_mle_attempts = 5, max_tries = 1000, k = 2, delta = 4, M =20, ztStep = 1, ztIntervals = NULL, ftStep = 1,  ftIntervals = NULL, lowZtSurvivalProb = 0.95, ...){
  # print warnings in real time, messages automatically print in real time
  options(warn = 1)
  # prep return data
  num_models <- length(model_list)
  dat_len <- k * num_models
  #cv_return_dat <- data.frame(model_name = as.character(dat_len), fold = as.integer(), log_lik = as.numeric, convergence = as.integer())
  cv_return_dat <- data.frame(matrix(NA, nrow = dat_len, ncol = 7))
  colnames(cv_return_dat) <- c("model_name", "fold", "neg_log_lik", "convergence", "init_error", "minutes", "iterations")
  
  # combine general (i.e. cross-species) control parameters (for optim) with species-species control parameter(i.e. parscale)
  for(j in 1:num_models)
  {
    if("control" %in% names(list(...)))
    {
      model_list[[j]][["control"]] <- c((list(...))[["control"]], (model_list[[j]])["parscale"])
    } else
    {
      model_list[[j]][["control"]] <- (model_list[[j]])["parscale"]
    }
  }
  
  # take "control" away from the ellipsis argument
  additional_args <- (list(...))[which(names(list(...)) != "control")]
  
  # split data into k subsets
  dat$fold <- sample(rep( 1:k, length.out = NROW(dat)))
  
  # loop over each fold
  for(i in 1:k)
  {
    train_dat <- dat %>% filter(fold != i)
    valid_dat <- dat %>% filter(fold == i)
    
    # put data into list (preferred input of cpp log lik functions)
    train_dat_list <- list(N = NROW(train_dat),
                           nt = train_dat$nt,
                           ntp1 = train_dat$ntp1,
                           oviPeriod = train_dat$oviposition_days,
                           delta = delta,
                           M = M,
                           ztStep = ztStep,
                           ztIntervals = ztIntervals,
                           ftStep = ftStep,
                           ftIntervals = ftIntervals,
                           lowZtSurvivalProb = lowZtSurvivalProb)
    
    
    valid_dat_list <- list(N = NROW(valid_dat),
                           nt = valid_dat$nt,
                           ntp1 = valid_dat$ntp1,
                           oviPeriod = valid_dat$oviposition_days,
                           delta = delta,
                           M = M,
                           ztStep = ztStep,
                           ztIntervals = ztIntervals,
                           ftStep = ftStep,
                           ftIntervals = ftIntervals,
                           lowZtSurvivalProb = lowZtSurvivalProb)
    
    # loop over each model
    for(j in 1:num_models)
    {
      # print progress
      cat("fold ", i, ", model ", j, " : ~ ", round(100*(((i-1)*num_models + j)/dat_len)), "% complete. \n", sep = "")
      
      t1 <- Sys.time() # start stopwatch
      
      # return a data frame row
      cv_return_dat[(i-1)*num_models + j,] <-
        withCallingHandlers({
          withRestarts({
            
            cur_mle <- robust_mle((model_list[[j]])[["fun"]],
                                  (model_list[[j]])[["par_init_fun"]],
                                  train_dat_list,
                                  max_mle_attempts,
                                  max_tries,
                                  control = (model_list[[j]])[["control"]],
                                  unlist(additional_args) ) # get mle on training data
            valid_ll <- model_list[[j]][["fun"]](cur_mle$par, valid_dat_list) # get LL on validation data
            t2 <- Sys.time()# end stopwatch
            list(model_name = names(model_list)[[j]], 
                 fold = i, 
                 neg_log_lik = valid_ll, 
                 convergence = ifelse(cur_mle$convergence == 0, 1, 0),
                 init_error = 0,
                 minutes = difftime(t2, t1, units = "mins")[[1]], 
                 iterations = cur_mle$counts[["function"]]) # return data
          },
          muffle_stop = function(){
            t2 <- Sys.time()
            list(model_name = names(model_list)[[j]], 
                 fold = i, 
                 neg_log_lik = NA, 
                 convergence = NA,
                 init_error = 1,
                 minutes = difftime(t2, t1, units = "mins")[[1]], 
                 iterations = NA) # return row anyways if optim function can't initialize
          }
          )
        },
        warning=function(w){
          warning(paste0("With model = ", names(model_list)[[j]], ", fold = ", k, " : ", w$message))
          rlang::cnd_muffle(w) # prevents warning from printing again on function exit
        },
        error=function(e){
          message(paste0("Muffled Error : With model = ", names(model_list)[[j]], ", fold = ", k, " : ", e$message))
          invokeRestart("muffle_stop") # skips the the the withRestarts block code; goes to muffle_stop
        })
      
    }
  }
  
  # retrieve data frame with model_name, and log lik
  return(cv_return_dat)
}


# 
# 
# 
# dat %>%
#   filter(temp > 25) %>%
#   ggplot( aes(x = nt, y = ntp1, color = species)) + 
#   geom_point() + 
#   facet_wrap(~exp_name, scales = "free")
# 
# 
# dat_A_crash <- dat %>%
#   filter(temp > 25 & species == "A") %>%
#   filter((exp_name == "castaneum_chaos") |
#            (exp_name == "MVM_priors"))  %>%
#   select(nt, ntp1, oviposition_days,exp_name) %>%
#   na.omit()
# 
# # select only castaneum chaos
# 
# dat_A_crash  <- dat_A_crash %>%
#   filter(exp_name == "castaneum_chaos")
# 
# model_list <- list(
#   mod1 = list(
#     fun = mod1LlCpp,
#     par_init_fun = function(i) {c(mu_alpha_uncon = lower_unconstrain(19),
#                                   beta0_uncon = lower_unconstrain(0.05),
#                                   eta_alpha_uncon = lower_unconstrain(20),
#                                   high_zt_survival_prob_uncon = lower_upper_unconstrain(0.05),
#                                   thresh_zt_range_uncon = lower_unconstrain(30),
#                                   mu_z_star_uncon = lower_unconstrain(500))},
#     parscale = abs(c(mu_alpha_uncon = lower_unconstrain(19),
#                      beta0_uncon = lower_unconstrain(0.05),
#                      eta_alpha_uncon = lower_unconstrain(20),
#                      high_zt_survival_prob_uncon = lower_upper_unconstrain(0.05),
#                      thresh_zt_range_uncon = lower_unconstrain(30),
#                      mu_z_star_uncon = lower_unconstrain(500)))
#   ),
#   mod2 = list(
#     fun = mod1LlCpp,
#     par_init_fun = function(i) {c(mu_alpha_uncon = lower_unconstrain(19),
#                                   beta0_uncon = lower_unconstrain(0.05),
#                                   eta_alpha_uncon = lower_unconstrain(20),
#                                   high_zt_survival_prob_uncon = lower_upper_unconstrain(0.05),
#                                   thresh_zt_range_uncon = lower_unconstrain(30),
#                                   mu_z_star_uncon = lower_unconstrain(500))},
#     parscale = abs(c(mu_alpha_uncon = lower_unconstrain(19),
#                      beta0_uncon = lower_unconstrain(0.05),
#                      eta_alpha_uncon = lower_unconstrain(20),
#                      high_zt_survival_prob_uncon = lower_upper_unconstrain(0.05),
#                      thresh_zt_range_uncon = lower_unconstrain(30),
#                      mu_z_star_uncon = lower_unconstrain(500)))
#   )
# )
# 
# model_list[[1]][["init_par_fun"]]()
# 
# max_mle_attempts = 5
# max_tries = 2
# k = 2
# delta = 4
# M =20
# ztStep = 1
# ftStep = 1
# lowZtSurvivalProb = 0.95
# #debug(k_fold_CV)
# cv_res <- k_fold_CV(model_list, dat = dat_A_crash, max_mle_attempts = 5, max_tries = 10, k = 2, delta = 4, M = 20, ztStep = 1, ftStep = 1, lowZtSurvivalProb = 0.95, control = list(maxit = 10))
# # str(cv_res)
# 
# train_dat <- dat_A_crash
# train_dat_list <- list(N = NROW(train_dat),
#                        nt = train_dat$nt,
#                        ntp1 = train_dat$ntp1,
#                        oviPeriod = train_dat$oviposition_days,
#                        delta = delta,
#                        M = M,
#                        ztStep = ztStep,
#                        ftStep = ftStep,
#                        lowZtSurvivalProb = lowZtSurvivalProb)
# 
# cur_mle <- robust_mle((model_list[[1]])[["fun"]],
#                       (model_list[[1]])[["par_init_fun"]],
#                       train_dat_list,
#                       max_mle_attempts,
#                       max_tries, control = list(maxit = 5))
# l <- list(control = list(maxit = 5))
# 
# model_list[[1]][["control"]] <- c((l)[["control"]], (model_list[[1]])["parscale"])
