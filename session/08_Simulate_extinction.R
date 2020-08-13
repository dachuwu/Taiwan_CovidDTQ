#### simulation part 0629
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

### Baseline setting
B <- list(
  
  cap_gen = 20,
  cap_max_days = Inf,
  cap_cases = 3000,
  
  r0isolated = 0,
  disp.iso = 1,
  relR.asym = .5, prop.asym = .15,

  
  incu_shape = 2.84, incu_scale = 1.93,
  delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158
  
  prop.presym = 0.524, si_omega = 2.689,
  disp.com = 19.196,
  
  detect_sen=.95,
  prop.ascertain = .9,
  quarant.days = 14,
  quarant.retro.days = 4
)


## Parameterise fixed paramters
sim_with_params <- purrr::partial(scenario_sim,
                                  save.case_data = F,
                                  n.sim = 100,
                                  
                                  cap_gen = 20,
                                  cap_max_days = Inf,
                                  cap_cases = 3000,
                                  
                                  r0isolated = 0,
                                  disp.iso = 1,
                                  relR.asym = .5, prop.asym = .15,
                                  
                                  
                                  incu_shape = 2.84, incu_scale = 1.93,
                                  delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158
                                  
                                  prop.presym = 0.524, si_omega = 2.689,
                                  disp.com = 19.196,
                                  
                                  detect_sen=.95,
                                  prop.ascertain = .9,
                                  quarant.days = 14,
                                  quarant.retro.days = 4
                                  )

## RUN
#
d1 = seq(.9,2.5,by = .05)
d2 = seq(5, 1000, by = 5)
length(d1)*length(d2)


ddp1 <- list(D1=d1, D2=d2[1:50])
ddp2 <- list(D1=d1, D2=d2[51:100])
ddp3 <- list(D1=d1, D2=d2[101:150])
ddp4 <- list(D1=d1, D2=d2[151:200])

#4
DDP <- list(ddp1, ddp2, ddp3, ddp4)
tag <- c("0730_ddp1","0730_ddp2","0730_ddp3","0730_ddp4")

for(i in 2:4){

  scenarios <- tidyr::expand_grid(
    Xr0community = DDP[[i]]$D1,
    Xnum.initial.cases = DDP[[i]]$D2
  ) %>% mutate(scenario = 1:dplyr::n())
  
  plan(multiprocess, workers = 10)
  RES <- future_lapply(1:NROW(scenarios), function(i){
    sim_with_params(r0community = scenarios$Xr0community[i],
                    num.initial.cases = scenarios$Xnum.initial.cases[i])
  })
  plan(sequential)
  
  df <- compile.sumdata(RES, scenarios)
  save(df, RES, scenarios, B, file = paste0("mid/ScenarioSim_Extinct_",tag[i],".RData"))

}



#######################################################
theme_set(theme_classic(base_size = 15) %+replace%
            theme(legend.position = "right", legend.justification = c(0,1),
                  legend.background = element_blank(),
                  legend.direction = "vertical", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                  plot.tag = element_text(face ="bold", size=26))
)


###
fps <- c("mid/ScenarioSim_Extinct_0730_ddp1.RData",
         "mid/ScenarioSim_Extinct_0730_ddp2.RData",
         "mid/ScenarioSim_Extinct_0730_ddp3.RData",
         "mid/ScenarioSim_Extinct_0730_ddp4.RData"
         )

c_df <- list()
c_scenario <- list()
for(fp in fps){
  load(fp)
  c_df <- rbind(c_df, df)
  c_scenario <- rbind(c_scenario, scenarios)
  rm(RES, df, scenarios)
}


#
gdf1 <- c_df %>%
  filter(est_name=="l_extinct" & Xr0community >=0.9) %>%
  arrange(Xnum.initial.cases,Xr0community)

gdfs <- gdf1 %>%
  group_by(Xnum.initial.cases) %>%
  mutate(ind1 = est >= 0.9)%>%
  filter(ind1)%>%
  mutate(ind1 = Xr0community == max(Xr0community))%>%
  filter(ind1)

G1 <- ggplot(gdf1)+
  geom_tile(aes(x=Xnum.initial.cases, y=Xr0community, fill=est),height=.05, width=10)+
  scale_fill_viridis_c(name="Probability of extinction")+
  coord_cartesian(expand = F, xlim=c(0,500))+
  guides(fill = guide_colorbar(title.position = "right", label.position = "left", barheight = unit(.8,"npc"),
                               title.theme = element_text(angle=270)))+
  labs(x="No. of initial introductions", 
       y=bquote("Effective reproduction number ( "~R[e]~")"))
G1

G2 <- readRDS( "mid/Gplot_SimEpiCurve_v0810.RDS")


require(gridExtra)
grid.arrange(G2+
               labs(tag="a")+
               theme(legend.position = c(1,.5), legend.justification = c(1,1),
                     legend.direction = "vertical"),
             G1+
               labs(tag="b"),
             nrow=1)


gunit <- 24.7
ggsave("output/NPI_Fig4_epicurv_extinction_v0810.jpeg",
       
       
       grid.arrange(G2+
                      labs(tag="a")+
                      theme(legend.position = c(1,.5), legend.justification = c(1,1),
                            legend.direction = "vertical"),
                    G1+
                      labs(tag="b"),
                    nrow=1)
       ,
       width = gunit, height = gunit*.4, units = "cm", dpi = 300, scale=1.7)

