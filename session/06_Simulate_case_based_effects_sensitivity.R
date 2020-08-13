### DTQ analysis 2020-05-20
library(DTQbp)
library(ggplot2)
library(parallel)
library(dplyr)
library(future)
library(future.apply)
###
compile.sumdata <- function(RES, scenarios){
  sdfs <- purrr::map(RES, 2)
  
  scen_sum <- function(d){
    
    calc_pui <- function(x){
      
      if(class(x) =="logical"){
        pest <- mean(x, na.rm = T)
        se <- sqrt(pest*(1-pest)/length(x))
        iest <- c(pest - 1.959964*se, pest + 1.959964*se)
        iest <- pmax(0, iest)
      } else {
        pest <- median(x, na.rm = T)
        iest <- quantile(x, probs = c(.025, .975), names = F, na.rm = T)
      }
      return(data.frame(est = pest, est_l = iest[1], est_u = iest[2]))
    }
    
    est_names <- colnames(d)
    
    tmp <- lapply(est_names, function(xn){ 
      d <- as.data.frame(d)
      a <- calc_pui(d[,xn])
      a$est_name <- xn
      a
    }) %>% bind_rows()
    
  }
  
  df <- lapply(1:length(sdfs), function(k){
    out <- scen_sum(sdfs[[k]])
    out$scenario <- k
    out
  }) %>% bind_rows()%>%
    left_join(scenarios) %>%
    arrange(est_name) 
  return(df)
}

gmean2para_incu <- function(x){
  varx <- 2.84*1.93*1.93 # DTQ158
  scl <- varx/x
  shp <- x/scl
  return(list(incu_shape=shp, incu_scale=scl))
}

gmean2para_delay <- function(x){
  
  off <- -1
  x <- x - off
  varx <- 1.08*5.58*5.58  # DTQ158
  scl <- varx/x
  shp <- x/scl
  
  return(list(delay_shape=shp, delay_scale=scl, delay_off=off))
}

gmean2para_delay(2.5)

### range of one-way sensitivity analysis
para.range <- list(
  mean_delay = c(0, 20.51),
  mean_incu = c(1.06, 13.45),
  r0community = c(2, 3),
  disp.com = c(19.2, 48.8),
  prop.asym = c(0.04, .41),
  relR.asym = c(0, 1),
  prop.presym = c(0.39, 0.67),
  si_omega = c(1.87, 3.79),
  detect_sen = c(0.75, 1),
  prop.ascertain = c(0.75, 1),
  num.initial.cases = c(10, 500)
)

### Baseline setting
B0 <- list(
  num.initial.cases = 100,
  cap_gen = 20, 
  cap_max_days = Inf,
  cap_cases = 10000,

  r0isolated = 0,
  disp.iso = 1,
  
  r0community = 2.5, ## 
  relR.asym = .5, prop.asym = .15,
  
  incu_shape = 2.84, incu_scale = 1.93,
  delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158
  
  prop.presym = 0.524, si_omega = 2.689,
  disp.com = 19.196,
  
  detect_sen=0.95,
  prop.ascertain = .9,
  quarant.days = 14,
  quarant.retro.days = 4
)
N.iter <- 100


plan(multiprocess, workers = 6)

B1 <- B0
sims0 <- scenario_sim(save.case_data = F, n.sim = N.iter,
                     num.initial.cases = B1$num.initial.cases,
                     cap_gen = B1$cap_gen,
                     cap_max_days = B1$cap_max_days,
                     cap_cases = B1$cap_cases,
                     
                     r0isolated = B1$r0isolated,
                     r0community = B1$r0community,
                     disp.iso = B1$disp.iso,
                     disp.com = B1$disp.com,
                     
                     prop.asym = B1$prop.asym,
                     relR.asym = B1$relR.asym,
                     prop.presym = B1$prop.presym,
                     si_omega = B1$si_omega,
                     
                     incu_shape = B1$incu_shape, incu_scale = B1$incu_scale,
                     delay_shape = B1$delay_shape, delay_scale = B1$delay_scale,
                     delay_off = B1$delay_off,
                     detect_sen = B1$detect_sen,
                     prop.ascertain = B1$prop.ascertain,
                     quarant.days = B1$quarant.days,
                     quarant.retro.days = B1$quarant.retro.days)

res0 <- as.data.frame(t(apply(sims0$sum_data,2, mean, na.rm = T)))%>%
  mutate(type = "baseline", variable = "baseline")



res <- future_lapply(1:length(para.range), function(i){
    
  xv <- names(para.range)[i]
  
  if(xv == "mean_incu"){
    para.lu <- lapply(para.range[[xv]],gmean2para_incu)
  } else if (xv == "mean_delay"){
    para.lu <- lapply(para.range[[xv]],gmean2para_delay)
  } else {
    para.lu <- lapply(para.range[[xv]], function(x){
      out <- list(x)
      names(out) <- xv
      out
    })
  }
   
  sub.res <- lapply(para.lu, function(X){
    
    B1 <- B0
    for(vv in names(X)) B1[[vv]] <- X[[vv]]
    
    sims <- scenario_sim(save.case_data = F, n.sim = N.iter,
                         num.initial.cases = B1$num.initial.cases,
                         cap_gen = B1$cap_gen,
                         cap_max_days = B1$cap_max_days,
                         cap_cases = B1$cap_cases,
                                      
                         r0isolated = B1$r0isolated,
                         r0community = B1$r0community,
                         disp.iso = B1$disp.iso,
                         disp.com = B1$disp.com,
                        
                         prop.asym = B1$prop.asym,
                         relR.asym = B1$relR.asym,
                         prop.presym = B1$prop.presym,
                         si_omega = B1$si_omega,
                        
                         incu_shape = B1$incu_shape, incu_scale = B1$incu_scale,
                         delay_shape = B1$delay_shape, delay_scale = B1$delay_scale,
                         delay_off = B1$delay_off,
                         detect_sen = B1$detect_sen,
                         prop.ascertain = B1$prop.ascertain,
                         quarant.days = B1$quarant.days,
                         quarant.retro.days = B1$quarant.retro.days)
    
    
    out <- as.data.frame(t(apply(sims$sum_data,2, mean, na.rm = T)))
    out
  }) %>% bind_rows() %>%
    mutate(type = c("lower", "upper"), variable = xv)
  sub.res
  
}) 
plan(sequential) 

df <- bind_rows(res, res0)
  
saveRDS(df, "mid/Sensitivity_tornado_0730.RDS")




## plot
xv2name<-list("disp.com"="Dispersion parameter ",
        "si_omega"="S.D. of the generation interval ",
        "detect_sen"="Sensitivity of case detection ",
        "prop.ascertain"="Prob. of contact ascertainment ",
        "prop.asym"="Prob. of asymp. infection ",
        "prop.presym"="Prop. of pre-symp. transmission ",
        "mean_incu"="Incubation period ",
        "r0community"=bquote("Counterfactual"~R[0]~""),
        "mean_delay"="Onset-to-isolation interval ",
        "relR.asym" = "Relative transmissibility ",
        "num.initial.cases" = "No. of initial introduction ") 




## plots
theme_set(theme_classic(base_size = 15) %+replace%
            theme(legend.position = c(1,0), legend.justification = c(1,0),
                  legend.background = element_blank(),
                  legend.direction = "vertical", 
                  axis.text.y =element_text(size = 14,hjust = 1, vjust = .5), 
                  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                  plot.tag = element_text(face ="bold", size=26))
)

plot.G1 <- function(df, gtag = NULL){
  
  xv <- "mu_Re"
  gdf0 <- df %>% filter(type == "baseline")
  gdf <- df %>% filter(type != "baseline") %>% 
    select(c("type","variable", xv)) %>%
    tidyr::pivot_wider(names_from = c("type"), values_from = xv) %>%
    mutate( r = abs(upper - lower)) %>%
    arrange(r)%>%
    mutate(y = 1:NROW(.),
           x0 = gdf0[,xv],
           xvlab = xv2name[variable],
           v.lower = purrr::map(para.range[variable],1),
           v.upper = purrr::map(para.range[variable],2)
            ) 
  gdf$upper[gdf$variable%in%c("detect_sen","prop.ascertain")] <- gdf0[,xv]
  
  colset <- c("Increased value"="#B22222", "Decreased value"="#4078a5")
 
  x01 <- gdf0[,xv]
  G1 <- ggplot(gdf)+
    geom_bar(aes(x=y, y=lower-x0, fill="Decreased value"), stat = "identity")+
    geom_bar(aes(x=y, y=upper-x0, fill="Increased value"), stat = "identity")+
    scale_x_continuous(breaks = 1:NROW(gdf), labels = gdf$xvlab)+
    scale_y_continuous(labels = function(x) round(x01+as.numeric(x), 2), expand = c(0,0.1)) +
    scale_fill_manual(values = colset, name=NULL)+
    labs(x = "", y="Effective reproduction number", tag = gtag)+
    geom_hline(yintercept = 0, linetype=1, color="grey40", size=1)+
    geom_text(aes(x=y,y=lower-x0,label=v.lower, hjust=if_else(lower-x0<=0.003, 1.2, -.1) ), size=3)+
    geom_text(aes(x=y,y=upper-x0,label=v.upper, hjust=if_else(lower-x0<=0.003, -.1, 1.2) ), size=3)+
    coord_flip()
  G1

  
}
plot.G2 <- function(df, gtag = NULL, move.leg = F){
  xv <- "l_extinct"
  gdf0 <- df %>% filter(type == "baseline")
  gdf2 <- df %>% filter(type != "baseline") %>% 
    select(c("type","variable", xv)) %>%
    tidyr::pivot_wider(names_from = c("type"), values_from = xv) 
  gdf2$upper[gdf2$variable%in%c("detect_sen","prop.ascertain")] <- gdf0[,xv]
  gdf2 <- gdf2 %>%
    mutate( r = abs(upper - lower)) %>%
    arrange(r)%>%
    mutate(y = 1:NROW(.),
           x0 = gdf0[,xv],
           xvlab = xv2name[variable],
           v.lower = purrr::map(para.range[variable],1),
           v.upper = purrr::map(para.range[variable],2)) 
  
  colset <- c("Increased value"="#B22222", "Decreased value"="#4078a5")
  
  x02 <- gdf0[,xv]
  G2 <- ggplot(gdf2)+
    geom_bar(aes(x=y, y=lower-x0, fill="Decreased value"), stat = "identity")+
    geom_bar(aes(x=y, y=upper-x0, fill="Increased value"), stat = "identity")+
    scale_x_continuous(breaks = 1:NROW(gdf2), labels = gdf2$xvlab)+
    scale_y_continuous(labels = function(x) round(x02+as.numeric(x), 2), expand = c(0,0.1)) +
    scale_fill_manual(values = colset, name=NULL)+
    labs(x = "", y = "Exctinction probability", tag = gtag)+
    geom_hline(yintercept = 0, linetype=1, color="grey40", size=1)+
    geom_text(aes(x=y,y=lower-x0,label=v.lower, hjust=if_else(lower-x0<=0, 1.2, -.1) ), size=3)+
    geom_text(aes(x=y,y=upper-x0,label=v.upper, hjust=if_else(lower-x0<=0, -.1, 1.2) ), size=3)+
    coord_flip() 
  if(move.leg){
    G2 + theme(legend.position = c(0,0), legend.justification = c(0,0))
  } else {
    G2
  }

}

## figure a, b
df1 <- readRDS("mid/Sensitivity_tornado_0730.RDS")

G1a <- plot.G1(df1, "a")
G1a
G1b <- plot.G2(df1, "b")

ggsave(filename = "output/Sensitivity_tornado_0730_R0.jpeg", 
       gridExtra::arrangeGrob(G1a, G1b,nrow=1, ncol=2),
       width = 24.7, height = 8, units = "cm", dpi = 300, scale = 1.8)



##1
df2 <- readRDS("mid/Sensitivity_tornado_0730_r1.RDS")

G1c <- plot.G1(df2, "c")
G1d <- plot.G2(df2, "d", move.leg = T)


ggsave(filename = "output/Sensitivity_tornado_0730_Rp.jpeg", 
       gridExtra::arrangeGrob(G1a, G1b, G1c, G1d, nrow=2, ncol=2),
       width = 24.7, height = 16, units = "cm", dpi = 300, scale = 1.7)
