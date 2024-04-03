# RMST calculation and comparison for ADT group ----------
d_rmst <- DF_AA %>% 
  select(FU, DEATH, 
         TTSRE, TTSRE_PRTB, TTSRE_CF, TTSRE_SCC, TTSRE_STB, 
         SRE, SRE_PRTB, SRE_CF, SRE_SCC, SRE_STB, BMA_WITHIN, weight) %>% 
  arrange(BMA_WITHIN) %>% 
  as.data.frame() 

d_rmst$BMA_WITHIN <- as.factor(d_rmst$BMA_WITHIN)

# restricted mean survival time for recurrence-free survival ----------
# time to palliative radiation to bone
armst_ttsre_prtb <- function(tau) {
  j = length(unique(d_rmst$BMA_WITHIN))
  
  rmst <- rep(9999, length(1:j))
  groupval <- rep(9999, length(1:j))
  rmst_var <- rep(9999, length(1:j))
  rmst_se <- rep(9999, length(1:j))
  
  for (i in 1:j){
    groupval[i] <- (levels(d_rmst$BMA_WITHIN)[i])
    dat_group <- d_rmst[which(d_rmst$BMA_WITHIN == (groupval[i])),]
    
    #--- AKM ---
    # Based on 'adjusted.KM' function from {IPWsurvival} package
    # Author: F. Le Borgne and Y. Foucher
    tj <- c(0, sort(unique(dat_group$TTSRE_PRTB[dat_group$SRE_PRTB == 1])))
    dj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_PRTB == x & dat_group$SRE_PRTB == 1])})
    yj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_PRTB >= x])})
    st <- cumprod(1 - (dj/yj))
    m <- sapply(tj, function(x) {sum((dat_group$weight[dat_group$TTSRE_PRTB >= x])^2)})
    mj <- ((yj^2)/m)
    #ft <- data.frame(time = tj, n_risk = yj, n_event = dj, survival = st, variable = i, m = mj)
    ft <- data.frame(tj, yj, dj, st, i, mj)
    
    #--- RMST ---
    # Based on 'rmst1 function' from {survRM2} package
    # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
    rtime <- ft$tj <= tau
    tj_r <- sort(c(ft$tj[rtime], tau))
    st_r <- ft$st[rtime]
    yj_r <- ft$yj[rtime]
    dj_r <- ft$dj[rtime]
    time_diff <- diff(c(0, tj_r))
    areas <- time_diff * c(1, st_r)
    rmst[i] <- sum(areas)
    
    mj_r <- ft$mj[rtime]
    var_r <- ifelse((yj_r-dj_r) == 0, 0, dj_r / (mj_r * (yj_r - dj_r)))
    #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
    var_r <- c(var_r, 0)
    rmst_var[i] <- sum(cumsum(rev(areas[-1])) ^ 2 * rev(var_r)[-1])
    rmst_se[i] <- sqrt(rmst_var[i])
  }
  
  # --- Compare RMST between groups and compile output---
  
  output <- tibble(rmst = 0, 
                   lcl = 0, 
                   ucl = 0)
  
  output$rmst <- rmst[2] - rmst[1] # no = 1, yes = 2
  output$lcl <- rmst[2] - rmst[1] - qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  output$ucl <- rmst[2] - rmst[1] + qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  
  print(output)
}

data_rmst <- data.frame()
for(i in seq(from = 0, to = 63.9, by = 0.5)) {
  data_rmst = rbind(data_rmst, armst_ttsre_prtb(i))
}

tau_vec <- as.data.frame(seq(from = 0, to = 63.9, by = 0.5))

plot_rmst <- cbind(tau_vec, data_rmst)
colnames(plot_rmst)[1] <- "tau"
colnames(plot_rmst)[2] <- "rmst"
colnames(plot_rmst)[3] <- "lcl"
colnames(plot_rmst)[4] <- "ucl"

# plot RMST difference for IPW-adjusted model for TTSRE
rmst_plot_ttsre_prtb <- plot_rmst %>% 
  ggplot(data = ., aes(x = tau)) +
  geom_hline(aes(yintercept = 0), color = 'grey40', size = 0.7) +
  geom_line(aes(y = rmst), size = 0.7, color = "#39568CFF") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.4, fill = "#39568CFF") +
  scale_y_continuous(limits = c(-5, 12), breaks = seq(from = -5, to = 12, by = 1)) +
  scale_x_continuous(limits = c(0, 66), breaks = seq(from = 0, to = 66, by = 12)) +
  labs(x = "Months since randomization", 
       y = "RMST difference (95% CI), month", 
       title = 'Time to palliative radiation to bone') +
  theme_jikei() 

# time to clinical fracture
armst_ttsre_cf <- function(tau) {
  j = length(unique(d_rmst$BMA_WITHIN))
  
  rmst <- rep(9999, length(1:j))
  groupval <- rep(9999, length(1:j))
  rmst_var <- rep(9999, length(1:j))
  rmst_se <- rep(9999, length(1:j))
  
  for (i in 1:j){
    groupval[i] <- (levels(d_rmst$BMA_WITHIN)[i])
    dat_group <- d_rmst[which(d_rmst$BMA_WITHIN == (groupval[i])),]
    
    #--- AKM ---
    # Based on 'adjusted.KM' function from {IPWsurvival} package
    # Author: F. Le Borgne and Y. Foucher
    tj <- c(0, sort(unique(dat_group$TTSRE_CF[dat_group$SRE_CF == 1])))
    dj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_CF == x & dat_group$SRE_CF == 1])})
    yj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_CF >= x])})
    st <- cumprod(1 - (dj/yj))
    m <- sapply(tj, function(x) {sum((dat_group$weight[dat_group$TTSRE_CF >= x])^2)})
    mj <- ((yj^2)/m)
    #ft <- data.frame(time = tj, n_risk = yj, n_event = dj, survival = st, variable = i, m = mj)
    ft <- data.frame(tj, yj, dj, st, i, mj)
    
    #--- RMST ---
    # Based on 'rmst1 function' from {survRM2} package
    # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
    rtime <- ft$tj <= tau
    tj_r <- sort(c(ft$tj[rtime], tau))
    st_r <- ft$st[rtime]
    yj_r <- ft$yj[rtime]
    dj_r <- ft$dj[rtime]
    time_diff <- diff(c(0, tj_r))
    areas <- time_diff * c(1, st_r)
    rmst[i] <- sum(areas)
    
    mj_r <- ft$mj[rtime]
    var_r <- ifelse((yj_r-dj_r) == 0, 0, dj_r / (mj_r * (yj_r - dj_r)))
    #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
    var_r <- c(var_r, 0)
    rmst_var[i] <- sum(cumsum(rev(areas[-1])) ^ 2 * rev(var_r)[-1])
    rmst_se[i] <- sqrt(rmst_var[i])
  }
  
  # --- Compare RMST between groups and compile output---
  
  output <- tibble(rmst = 0, 
                   lcl = 0, 
                   ucl = 0)
  
  output$rmst <- rmst[2] - rmst[1] # no = 1, yes = 2
  output$lcl <- rmst[2] - rmst[1] - qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  output$ucl <- rmst[2] - rmst[1] + qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  
  print(output)
}

data_rmst <- data.frame()
for(i in seq(from = 0, to = 63.9, by = 0.5)) {
  data_rmst = rbind(data_rmst, armst_ttsre_cf(i))
}

tau_vec <- as.data.frame(seq(from = 0, to = 63.9, by = 0.5))

plot_rmst <- cbind(tau_vec, data_rmst)
colnames(plot_rmst)[1] <- "tau"
colnames(plot_rmst)[2] <- "rmst"
colnames(plot_rmst)[3] <- "lcl"
colnames(plot_rmst)[4] <- "ucl"

# plot RMST difference for IPW-adjusted model for TTSRE
rmst_plot_ttsre_cf <- plot_rmst %>% 
  ggplot(data = ., aes(x = tau)) +
  geom_hline(aes(yintercept = 0), color = 'grey40', size = 0.7) +
  geom_line(aes(y = rmst), size = 0.7, color = "#39568CFF") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.4, fill = "#39568CFF") +
  scale_y_continuous(limits = c(-5, 12), breaks = seq(from = -5, to = 12, by = 1)) +
  scale_x_continuous(limits = c(0, 66), breaks = seq(from = 0, to = 66, by = 12)) +
  labs(x = "Months since randomization", 
       y = "RMST difference (95% CI), month", 
       title = 'Time to clinical fracture') +
  theme_jikei() 

# time to spinal cord compression
armst_ttsre_scc <- function(tau) {
  j = length(unique(d_rmst$BMA_WITHIN))
  
  rmst <- rep(9999, length(1:j))
  groupval <- rep(9999, length(1:j))
  rmst_var <- rep(9999, length(1:j))
  rmst_se <- rep(9999, length(1:j))
  
  for (i in 1:j){
    groupval[i] <- (levels(d_rmst$BMA_WITHIN)[i])
    dat_group <- d_rmst[which(d_rmst$BMA_WITHIN == (groupval[i])),]
    
    #--- AKM ---
    # Based on 'adjusted.KM' function from {IPWsurvival} package
    # Author: F. Le Borgne and Y. Foucher
    tj <- c(0, sort(unique(dat_group$TTSRE_SCC[dat_group$SRE_SCC == 1])))
    dj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_SCC == x & dat_group$SRE_SCC == 1])})
    yj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_SCC >= x])})
    st <- cumprod(1 - (dj/yj))
    m <- sapply(tj, function(x) {sum((dat_group$weight[dat_group$TTSRE_SCC >= x])^2)})
    mj <- ((yj^2)/m)
    #ft <- data.frame(time = tj, n_risk = yj, n_event = dj, survival = st, variable = i, m = mj)
    ft <- data.frame(tj, yj, dj, st, i, mj)
    
    #--- RMST ---
    # Based on 'rmst1 function' from {survRM2} package
    # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
    rtime <- ft$tj <= tau
    tj_r <- sort(c(ft$tj[rtime], tau))
    st_r <- ft$st[rtime]
    yj_r <- ft$yj[rtime]
    dj_r <- ft$dj[rtime]
    time_diff <- diff(c(0, tj_r))
    areas <- time_diff * c(1, st_r)
    rmst[i] <- sum(areas)
    
    mj_r <- ft$mj[rtime]
    var_r <- ifelse((yj_r-dj_r) == 0, 0, dj_r / (mj_r * (yj_r - dj_r)))
    #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
    var_r <- c(var_r, 0)
    rmst_var[i] <- sum(cumsum(rev(areas[-1])) ^ 2 * rev(var_r)[-1])
    rmst_se[i] <- sqrt(rmst_var[i])
  }
  
  # --- Compare RMST between groups and compile output---
  
  output <- tibble(rmst = 0, 
                   lcl = 0, 
                   ucl = 0)
  
  output$rmst <- rmst[2] - rmst[1] # no = 1, yes = 2
  output$lcl <- rmst[2] - rmst[1] - qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  output$ucl <- rmst[2] - rmst[1] + qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  
  print(output)
}

data_rmst <- data.frame()
for(i in seq(from = 0, to = 63.9, by = 0.5)) {
  data_rmst = rbind(data_rmst, armst_ttsre_scc(i))
}

tau_vec <- as.data.frame(seq(from = 0, to = 63.9, by = 0.5))

plot_rmst <- cbind(tau_vec, data_rmst)
colnames(plot_rmst)[1] <- "tau"
colnames(plot_rmst)[2] <- "rmst"
colnames(plot_rmst)[3] <- "lcl"
colnames(plot_rmst)[4] <- "ucl"

# plot RMST difference for IPW-adjusted model for TTSRE
rmst_plot_ttsre_scc <- plot_rmst %>% 
  ggplot(data = ., aes(x = tau)) +
  geom_hline(aes(yintercept = 0), color = 'grey40', size = 0.7) +
  geom_line(aes(y = rmst), size = 0.7, color = "#39568CFF") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.4, fill = "#39568CFF") +
  scale_y_continuous(limits = c(-5, 12), breaks = seq(from = -5, to = 12, by = 1)) +
  scale_x_continuous(limits = c(0, 66), breaks = seq(from = 0, to = 66, by = 12)) +
  labs(x = "Months since randomization", 
       y = "RMST difference (95% CI), month", 
       title = 'Time to spinal cord compression') +
  theme_jikei() 

# time to spinal cord compression
armst_ttsre_stb <- function(tau) {
  j = length(unique(d_rmst$BMA_WITHIN))
  
  rmst <- rep(9999, length(1:j))
  groupval <- rep(9999, length(1:j))
  rmst_var <- rep(9999, length(1:j))
  rmst_se <- rep(9999, length(1:j))
  
  for (i in 1:j){
    groupval[i] <- (levels(d_rmst$BMA_WITHIN)[i])
    dat_group <- d_rmst[which(d_rmst$BMA_WITHIN == (groupval[i])),]
    
    #--- AKM ---
    # Based on 'adjusted.KM' function from {IPWsurvival} package
    # Author: F. Le Borgne and Y. Foucher
    tj <- c(0, sort(unique(dat_group$TTSRE_STB[dat_group$SRE_STB == 1])))
    dj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_STB == x & dat_group$SRE_STB == 1])})
    yj <- sapply(tj, function(x) {sum(dat_group$weight[dat_group$TTSRE_STB >= x])})
    st <- cumprod(1 - (dj/yj))
    m <- sapply(tj, function(x) {sum((dat_group$weight[dat_group$TTSRE_STB >= x])^2)})
    mj <- ((yj^2)/m)
    #ft <- data.frame(time = tj, n_risk = yj, n_event = dj, survival = st, variable = i, m = mj)
    ft <- data.frame(tj, yj, dj, st, i, mj)
    
    #--- RMST ---
    # Based on 'rmst1 function' from {survRM2} package
    # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
    rtime <- ft$tj <= tau
    tj_r <- sort(c(ft$tj[rtime], tau))
    st_r <- ft$st[rtime]
    yj_r <- ft$yj[rtime]
    dj_r <- ft$dj[rtime]
    time_diff <- diff(c(0, tj_r))
    areas <- time_diff * c(1, st_r)
    rmst[i] <- sum(areas)
    
    mj_r <- ft$mj[rtime]
    var_r <- ifelse((yj_r-dj_r) == 0, 0, dj_r / (mj_r * (yj_r - dj_r)))
    #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
    var_r <- c(var_r, 0)
    rmst_var[i] <- sum(cumsum(rev(areas[-1])) ^ 2 * rev(var_r)[-1])
    rmst_se[i] <- sqrt(rmst_var[i])
  }
  
  # --- Compare RMST between groups and compile output---
  
  output <- tibble(rmst = 0, 
                   lcl = 0, 
                   ucl = 0)
  
  output$rmst <- rmst[2] - rmst[1] # no = 1, yes = 2
  output$lcl <- rmst[2] - rmst[1] - qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  output$ucl <- rmst[2] - rmst[1] + qnorm(1 - 0.05 / 2) * sqrt(rmst_var[2] + rmst_var[1]) # no = 1, yes = 2
  
  print(output)
}

data_rmst <- data.frame()
for(i in seq(from = 0, to = 63.9, by = 0.5)) {
  data_rmst = rbind(data_rmst, armst_ttsre_stb(i))
}

tau_vec <- as.data.frame(seq(from = 0, to = 63.9, by = 0.5))

plot_rmst <- cbind(tau_vec, data_rmst)
colnames(plot_rmst)[1] <- "tau"
colnames(plot_rmst)[2] <- "rmst"
colnames(plot_rmst)[3] <- "lcl"
colnames(plot_rmst)[4] <- "ucl"

# plot RMST difference for IPW-adjusted model for TTSRE
rmst_plot_ttsre_stb <- plot_rmst %>% 
  ggplot(data = ., aes(x = tau)) +
  geom_hline(aes(yintercept = 0), color = 'grey40', size = 0.7) +
  geom_line(aes(y = rmst), size = 0.7, color = "#39568CFF") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.4, fill = "#39568CFF") +
  scale_y_continuous(limits = c(-5, 12), breaks = seq(from = -5, to = 12, by = 1)) +
  scale_x_continuous(limits = c(0, 66), breaks = seq(from = 0, to = 66, by = 12)) +
  labs(x = "Months since randomization", 
       y = "RMST difference (95% CI), month", 
       title = 'Time to surgery to bone') +
  theme_jikei() 

# merge figures
merge_aa <- ggarrange(rmst_plot_ttsre_prtb, rmst_plot_ttsre_cf,
                      rmst_plot_ttsre_scc, rmst_plot_ttsre_stb, 
                      nrow = 2, ncol = 2)

ggsave("rmstd_aap_detail.pdf", merge_aa,
       units = "in", width = 8, height = 8)

# summary statistics
armst_ttsre_prtb(63.9)
armst_ttsre_cf(63.9)
armst_ttsre_scc(63.9)
armst_ttsre_stb(63.9)