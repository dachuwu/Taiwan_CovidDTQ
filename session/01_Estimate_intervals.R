rm(list=ls())
setwd("C:/Users/user/Desktop/GitDaChu/COVID19_epi_parameters")
source('./R/Functions.R')
require(purrr)
require(rstan)
require(dplyr)


### import data
load( file = "mid/Observation_DTQ158_v0729.RData" )

### read results
if(T){
  DateUpdate <<- "2020-07-29"
  # SI
  fps <- grep("_DTQ158_v0729",dir("mid/est_for_DTQ"),  value = T) 
  lapply(fps, function(fp){
    tmp <- readRDS(paste0("mid/est_for_DTQ/", fp))
    Res <- tmp$Res
    nn <- NROW(tmp$sdf)
    G <- plot.dur(Res)
    exportPlot(G, sub(".RDS", "",fp),"#gEPICc")
  })
 
  
}

### intervals estimation
if(F){
  
  # estimation function
  est.dur <- function(t0L, t0R, t1L, t1R, min_dur=0, dlbs = NULL,
                      g=NULL, glbs = NULL, opt=NULL, 
                      MOD = MOD_dur, ...){
    
    ori <- min(c(t0L,t0R,t1L,t1R),na.rm = T)
    rscl.dt <- function(x) as.numeric(difftime(x, ori, units = "days"))+1
    
    if(is.null(g)){ 
      g <- rep(1, length(t0L))
      glbs <- "Overall"
    }
    glvs <- as.character(sort(unique(g)))
    
    D <- list(
      N = length(t0L), 
      t0L = rscl.dt(t0L),
      t0R = rscl.dt(t0R),
      t1L = rscl.dt(t1L),
      t1R = rscl.dt(t1R),
      g = g,
      Ng = length(unique(g)),
      Neval=200,
      Xeval = seq(min_dur+0.001, min_dur+25, length.out = 200),
      min_dur=as.numeric(min_dur)
    )
    
    
    fit <- sampling(object = MOD, data = D, save_warmup=F, ... = ...)
    
    
    tmp <- extract(fit)
    d <- lapply(1:D$Ng, function(k){
      lapply(1:length(D$Xeval), function(j){
        xxs <- tmp$PDFeval[,k,j]
        qxxs <- quantile(xxs, c(0.025,0.5,0.975),na.rm=T)
        xx0s <- tmp$pri_PDFeval[,1,j]
        qxx0s <- quantile(xx0s, c(0.025,0.5,0.975),na.rm=T)
        data.frame("mean"=mean(xxs,na.rm=T),
                   "q2.5"= qxxs[[1]],"q50"= qxxs[[2]],"q97.5"= qxxs[[3]],
                   "pri_mean"=mean(xx0s,na.rm=T),
                   "pri_q2.5"= qxx0s[[1]],"pri_q50"= qxx0s[[2]],"pri_q97.5"= qxx0s[[3]],
                   xeval=D$Xeval[j], g=glbs[k], 
                   stringsAsFactors = F
        ) 
        
      }) %>% bind_rows()
    }) %>% bind_rows() 
    
    R1<-ggplot(d)+
      geom_ribbon(aes(x=xeval,ymin=q2.5,ymax=q97.5,group=g,fill=g),alpha=.25)+
      geom_path(aes(x=xeval,y=q50,group=g,color=g),size=1)+
      geom_ribbon(aes(x=xeval,ymin=pri_q2.5,ymax=pri_q97.5,fill="Prior"),alpha=.25)+
      geom_path(aes(x=xeval,y=pri_q50,color="Prior"),size=1)+
      facet_wrap(~g, ncol = 3)+ 
      coord_cartesian(ylim=c(0,0.2))+
      scale_fill_d3(breaks=c(glbs,"Prior"),name=NULL)+ 
      scale_color_d3(breaks=c(glbs,"Prior"),name=NULL)+
      labs(x="Duration (day)", y="Density")+
      theme(legend.position = "bottom")
    
    R2 <- summary(fit,pars=c("g_mean","g_sd","shp", "scl"))$summary %>% 
      as.data.frame() %>%
      round(2)%>%
      mutate(`95% CrI`= paste0(`50%`," (",`2.5%`,",",`97.5%`,")"))%>%
      mutate(group = rep(glbs, times=4))%>%
      mutate(parameter = rep(c("mean","sd","shape", "scale"),each=D$Ng)) %>%
      arrange(group) %>%
      select(c("group","parameter","95% CrI")) %>%
      pivot_wider(values_from = 3, names_from = 2)
    
    Res <- list(D=D, fit=fit, gdf=d, R1=R1, R2=R2, 
                glbs=glbs,dlbs=dlbs)
    if(!is.null(opt)){ 
      if(!dir.exists("mid")) dir.create("mid")
      saveRDS(Res, paste0("mid/Est_",opt,".RDS"))}
    return(Res)
  }
  
  # model specification
  MOD_dur1 <- stan_model('stan/Mod_duration_single.stan')
  
  
  
  
  
}





### 01. serial interval
df <- case_obs %>% 
  filter(!is.na(S0L) & !is.na(S1L) & !is.na(S0R) & !is.na(S1R) )  

si_obs <- as.numeric(difftime(df$S1R,df$S0R, units = "days"))
mean(si_obs)

Res <- est.dur(t0L = df$S0L, t0R = df$S0R, t1L = df$S1L, t1R = df$S1R,
               dlbs = df$to,
               g = NULL, glbs=NULL,
               MOD = MOD_dur1, opt=NULL, min_dur = -1,
               iter = 5000, warmup = 4000, chains =2, cores =1,
               control = list(adapt_delta=0.95,max_treedepth=12)
)
Res$R2
Res$R1
saveRDS(list(sdf=df, Res=Res), paste0("mid/Est_si_case_DTQ170.RDS"))



# Incubation period
df <- case_obs %>%
  filter(
    !is.na(S1L) & !is.na(S1R) &  !is.na(CL) & !is.na(CR)) 

Res <- est.dur(t0L = df$CL, t0R = df$CR, t1L = df$S1L, t1R = df$S1R,
               dlbs = df$to,
               g = NULL, glbs=NULL,
               MOD = MOD_dur1, opt=NULL, min_dur = 0,
               iter = 5000, warmup = 4000, chains = 2, cores = 1,
               control = list(adapt_delta=0.95,max_treedepth=12))
Res$R1
Res$R2
saveRDS(list(sdf=df, Res=Res), paste0("mid/Est_incu_case_DTQ170.RDS"))


# Delay ( total )
df <- df0 %>% getCaseInfo()%>%
  filter( !is.na(S1L) & !is.na(S1R) & !is.na(Test) )  # & (S0L==S0R) & (S1L==S1R)


Res <- est.dur(t0L = df$S1L, t0R = df$S1R, t1L = df$Test, t1R = df$Test,
               dlbs = df$Info,
               g = NULL, glbs=NULL,
               MOD = MOD_dur1, opt=NULL, min_dur = -1,
               iter = 3500, warmup = 2500, chains = 4, cores = 4,
               control = list(adapt_delta=0.95,max_treedepth=12)
)
Res$R2
Res$R1
saveRDS(list(sdf=df, Res=Res), paste0("mid/Est_delay_case_DTQ.RDS"))


# Delay ( local )
df <- df0 %>% getCaseInfo()%>%
  filter( !is.na(S1L) & !is.na(S1R) & !is.na(Test) & Imported !=1)  # & (S0L==S0R) & (S1L==S1R)


Res <- est.dur(t0L = df$S1L, t0R = df$S1R, t1L = df$Test, t1R = df$Test,
               dlbs = df$Info,
               g = NULL, glbs=NULL,
               MOD = MOD_dur1, opt=NULL, min_dur = -1,
               iter = 3500, warmup = 2500, chains = 4, cores = 1,
               control = list(adapt_delta=0.95,max_treedepth=12)
)
Res$R2
Res$R1
saveRDS(list(sdf=df, Res=Res), paste0("mid/Est_delay_local_case_DTQ.RDS"))


