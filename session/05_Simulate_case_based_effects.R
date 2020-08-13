#### simulation part 0629
### DTQ analysis 2020-05-20
library(DTQbp)
library(ggplot2)
library(parallel)
library(dplyr)
library(future)
library(future.apply)

R_0 <<- c(2.5, 2.5, 2.5)
R_p <<- c(0.7, 0.518, 0.903)


### simulation material
if(T){
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
    #num.initial.cases = 500,
    cap_gen = 20,
    cap_max_days = Inf,
    cap_cases = 10000,
    
    r0isolated = 0,
    disp.iso = 1, ##
    
    relR.asym = .5, prop.asym = .15,
    
    incu_shape = 2.84, incu_scale = 1.93,
    #delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158
    delay_shape = 0.3642872, delay_scale = 9.607803,  delay_off = -1,  # when delay = 2.5
    
    prop.presym = 0.524, si_omega = 2.689,
    disp.com = 19.196,
    #r0community = 0.729,
    
    #detect_sen=1,
    #prop.ascertain = .95,
    #quarant.days = 14,
    quarant.retro.days = 4
  )
  
  
  ### multi scenarios
  scenarios <- tidyr::expand_grid(
    Xr0community = c(R_p[1], R_0[1]),
    Xnum.initial.cases = c(100),
    data.frame(Xdetect_sen = c(0, .95, .95, .95, .95),
               Xprop.ascertain = c(0, 0, .9, .9, .9),
               Xquarant.days = c(0, 0, 0, 7, 14))
  ) %>% mutate(scenario = 1:dplyr::n())
  
  
  ## Parameterise fixed paramters
  sim_with_params <- purrr::partial(scenario_sim,
                                    save.case_data = F,
                                    n.sim = 100,
                                    #num.initial.cases = B$num.initial.cases,
                                    cap_gen = B$cap_gen,
                                    cap_max_days = B$cap_max_days,
                                    cap_cases = B$cap_cases,
                                    
                                    r0isolated = B$r0isolated,
                                    #r0community = B$r0community,
                                    disp.iso = B$disp.iso,
                                    disp.com = B$disp.com,
                                    
                                    prop.asym = B$prop.asym,
                                    relR.asym = B$relR.asym,
                                    prop.presym = B$prop.presym,
                                    si_omega = B$si_omega,
                                    
                                    incu_shape = B$incu_shape, incu_scale = B$incu_scale,
                                    delay_shape = B$delay_shape, delay_scale = B$delay_scale,
                                    delay_off = B$delay_off,
                                    #detect_sen = B$detect_sen,
                                    #prop.ascertain = B$prop.ascertain,
                                    #quarant.days = B$quarant.days,
                                    quarant.retro.days = B$quarant.retro.days)
  
}


## RUN
CL <- makeCluster(10)
clusterExport(CL, ls())
clusterEvalQ(CL, {library(DTQbp)})

RES <- parLapply(CL, 1:NROW(scenarios), function(i){
  sim_with_params(prop.ascertain = scenarios$Xprop.ascertain[i],
                  detect_sen = scenarios$Xdetect_sen[i],
                  quarant.days = scenarios$Xquarant.days[i],
                  r0community = scenarios$Xr0community[i],
                  num.initial.cases = scenarios$Xnum.initial.cases[i])})
stopCluster(CL)

df <- compile.sumdata(RES, scenarios[10,])
save(df, RES, scenarios, B, file = "mid/ScenSim_I100_0730.RData")

####################################################################################
# rm(list=ls())


fplot1 <- function(gdf, note = NULL){
  
  theme_set(theme_classic(base_size = 15) %+replace%
              theme(legend.position = c(1,1), legend.justification = c(1,1),
                    legend.background = element_blank(),
                    legend.direction = "horizontal", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                    plot.tag = element_text(face ="bold", size=26))
  )
  
  colset1 <- c("R0"="grey60","Re"="#8B0000")
  colset2 <- viridisLite::viridis(4)[1:3]
  
  
  
  gdf1 <- gdf %>% filter(est_name=="mu_Re")
  ymax1 <- max(gdf1$est_u)
  ymin1 <- min(min(gdf1$est_l),.8)
  gdf4 <- gdf %>% filter(est_name=="l_extinct") 
  ymax4 <- max(max(gdf4$est_u), 0.05)
  ymin4 <- min(gdf4$est_l)
  func.y2 <- function(x){
    (x-ymin4)*(ymax1-ymin1)/(ymax4-ymin4) + ymin1
  }
  func.y2inv <- function(x){
    (x - ymin1)/(ymax1-ymin1)*(ymax4-ymin4)+ ymin4
  }
  yn1 <- bquote("Effective reproduction number"~~(italic(R[e]))) #
  yn2 <- bquote("Extinction probability "~~(italic(p[ext]))) #~p^{ext}~
  
  G1 <- ggplot(gdf1,aes(x=interv, y=est, ymin=est_l, ymax=est_u))+
    geom_pointrange(color=colset1["Re"],
                    position = position_nudge(x = -0.05))+
    geom_hline(yintercept = 1, linetype=2, color="red4", alpha=.5)+
    
    geom_pointrange(data = gdf4,
                    aes(x=interv, y=func.y2(est), ymin=func.y2(est_l), ymax=func.y2(est_u)),
                    color="black",
                    position = position_nudge(x = 0.05))+
    scale_y_continuous(sec.axis = sec_axis(func.y2inv, 
                                           name = yn2))+
    
    coord_cartesian(ylim=c(ymin1,ymax1))+
    labs(y =yn1, x="Interventions\n")+
    theme(axis.title.y.left = element_text(color = "red4"), 
          axis.text.y.left= element_text(color = "red4"),
          axis.line.y.left= element_line(color = "red4"),
          axis.title.y.right = element_text(color = "black"), 
          axis.text.y.right= element_text(color = "black"),
          axis.line.y.right= element_line(color = "black")
          )
  G1
  
  
  gdf2 <- gdf %>% filter(est_name=="k_Re")
  ymax2 <- max(gdf2$est_u)
  ymin2 <- min(gdf2$est_l)
  G2 <- ggplot(gdf2,aes(x=interv, y=est, ymin=est_l, ymax=est_u))+
    geom_pointrange(color=colset1["Re"])+
    coord_cartesian(ylim=c(ymin2, ymax2))+
    labs(y ="Dispersion parameter", x="Interventions\n")
  G2
  
  
  
  
  gdf3 <- gdf %>% filter(est_name=="pgt1_re")
  ymax3 <- max(gdf3$est_u)
  ymin3 <- min(gdf3$est_l)
  G3 <- ggplot(gdf3, aes(x=interv, y=est, ymin=est_l, ymax=est_u))+
    geom_pointrange(color=colset1["Re"])+
    coord_cartesian(ylim=c(ymin3,ymax3))+
    labs(y ="Prob. R>1", x="Interventions\n")
  G3
  
  
  gdf4 <- gdf %>% filter(est_name=="l_extinct")
  ymax4 <- max(max(gdf4$est_u), 0.05)
  ymin4 <- min(gdf4$est_l)
  G4 <- ggplot(gdf4)+
    geom_pointrange(aes(x=interv, y=est, ymin=est_l, ymax=est_u), color=colset1["Re"])+
    coord_cartesian(ylim=c(0,ymax4))+
    labs(y="Probability of extinction", x="Interventions\n")
  G4
  
  
  gdf5_1 <- gdf %>% filter(est_name=="psrc_ngen_gt0")
  gdf5_2 <- gdf %>% filter(est_name=="psrc_ngen_gt1")
  gdf5_3 <- gdf %>% filter(est_name=="psrc_ngen_gt2")
  collab <- c(
    bquote("">=~"2nd-"),
    bquote("">=~"3rd-"),
    bquote("">=~"4th-wave   ")
  )
  
  ymax5 <- min(max(gdf5_1$est_u)*1.3, 1)
  G5 <- ggplot(gdf5_1, aes(x=interv, y=est, ymin=est_l, ymax=est_u))+
    geom_pointrange(aes(color="2nd-wave"))+
    geom_pointrange(data=gdf5_2, aes(color="3rd-wave"),
                    position = position_nudge(x = 0.1))+
    geom_pointrange(data=gdf5_3, aes(color="4th-wave"),
                    position = position_nudge(x = 0.2))+
    scale_color_manual(values = colset2, name=NULL, labels=collab)+
    guides(color = guide_legend(override.aes = list(linetype=0)) )+
    scale_y_continuous(breaks = seq(0, 1, .25),sec.axis = dup_axis(~.,  name = " "))+
    coord_cartesian(ylim=c(0, 1.1))+
    labs(y="Probability of subsequent \ntransmission (initial)",
         x="Interventions\n")
  G5
  
  if(T){
    gdf6_1 <- gdf %>% filter(est_name=="pasymp_ngen_gt0")
    gdf6_2 <- gdf %>% filter(est_name=="pasymp_ngen_gt1")
    gdf6_3 <- gdf %>% filter(est_name=="pasymp_ngen_gt2")
    
    ymax6 <- min(max(gdf6_1$est_u)*1.3, 1)
    G6 <-  ggplot(gdf6_1,
                  aes(x=interv, y=est, ymin=est_l, ymax=est_u))+
      geom_pointrange(aes(color="2nd-wave"))+
      geom_pointrange(data=gdf6_2, aes(color="3rd-wave"),
                      position = position_nudge(x = 0.1))+
      geom_pointrange(data=gdf6_3, aes(color="4th-wave"),
                      position = position_nudge(x = 0.2))+
      scale_color_manual(values = colset2,name=NULL, labels=collab)+
      scale_y_continuous(breaks = seq(0, 1, .25), sec.axis = dup_axis(~.,  name = " "))+
      coord_cartesian(ylim=c(0, 1.1))+
      labs(y="Probability of subsequent \ntransmission (asymptomatic)",
           x="Interventions\n")
    G6
  } else{
    gdf6_1 <- gdf %>% filter(est_name=="pasymp_ngen_gt0")
  gdf6_2 <- gdf %>% filter(est_name=="pasymp_ngen_gt1")
  gdf6_3 <- gdf %>% filter(est_name=="pasymp_ngen_gt2")
  gdf6_4 <- gdf %>% filter(est_name=="psymp_ngen_gt0")
  gdf6_5 <- gdf %>% filter(est_name=="psymp_ngen_gt1")
  gdf6_6 <- gdf %>% filter(est_name=="psymp_ngen_gt2")
  
  ymax6 <- min(max(gdf6_1$est_u, gdf6_4$est_u)*1.3, 1)
  
  G6 <-  ggplot(gdf6_1,
                aes(x=interv, y=est, ymin=est_l, ymax=est_u))+
    geom_pointrange(aes(color="2nd-wave", linetype="asym", shape="asym", group="asym"), fill="white",
                    position = position_nudge(x = -0.1))+
    geom_pointrange(data=gdf6_2, aes(color="3rd-wave", linetype="asym", shape="asym", group="asym"), fill="white",
                    position = position_nudge(x = -0.05))+
    geom_pointrange(data=gdf6_3, aes(color="4th-wave", linetype="asym", shape="asym", group="asym"), fill="white",
                    position = position_nudge(x = 0))+
    
    geom_pointrange(data=gdf6_4,aes(color="2nd-wave", linetype="sym", shape="sym", group="sym"),
                    position = position_nudge(x = 0.05))+
    geom_pointrange(data=gdf6_5, aes(color="3rd-wave", linetype="sym", shape="sym", group="sym"),
                    position = position_nudge(x = 0.1))+
    geom_pointrange(data=gdf6_6, aes(color="4th-wave", linetype="sym", shape="sym", group="sym"),
                    position = position_nudge(x = 0.15))+
    
    scale_shape_manual(values = c("asym"=21, "sym"=19), name="", labels=c("asymptomatic","symptomatic"))+
    scale_linetype_manual(values = c("asym"=2, "sym"=1), name="", labels=c("asymptomatic","symptomatic"))+
    
    scale_color_manual(values = colset2,name=NULL, labels=collab)+
    scale_y_continuous(breaks = seq(0, 1, .25), sec.axis = dup_axis(~.,  name = " "))+
    guides(color = guide_legend(override.aes = list(linetype=0)),
           linetype= guide_legend(override.aes = list(fill="white")))+
    coord_cartesian(ylim=c(0, 1.1))+
    labs(y="Probability of subsequent \ntransmission ",
         x="Interventions\n")+
    theme(legend.spacing = unit(0, "npc"), legend.margin = margin(2, 4, 0, 0))
  G6
  }
  
  
  #gridExtra::grid.arrange(G1, G2, G3, G4, G5, G6,nrow=2,ncol=3, top=note)
  # gridExtra::grid.arrange(G1,
  #                         G4,
  #                         G5,
  #                         G6,
  #                         nrow=1,ncol=4, top=note)
  list(G1, G2, G3, G4, G5, G6)
}
fplot_save <- function(GG, pt){
  
  gunit <- 24.7/length(G$widths)
  
  ggsave(filename = pt, plot = GG,
         width = gunit*length(G$widths), height = .8*gunit*length(G$heights), 
         units = "cm", dpi = 300, scale=1.8)
}


###
load(file = "mid/ScenSim_I100_0730.RData")

# R_0
gdf <- df %>%
  mutate(interv = case_when(
    Xdetect_sen==0&Xprop.ascertain==0&Xquarant.days==0 ~ "N",
    Xdetect_sen==.95&Xprop.ascertain==0&Xquarant.days==0 ~ "D",
    Xdetect_sen==.95&Xprop.ascertain==.9&Xquarant.days==0 ~ "DT",
    Xdetect_sen==.95&Xprop.ascertain==.9&Xquarant.days==7 ~ "DTQ7",
    Xdetect_sen==.95&Xprop.ascertain==.9&Xquarant.days==14 ~ "DTQ14",
    T ~ "x"
  )) %>%
  filter(interv!="x" & Xnum.initial.cases == 100 & Xr0community==R_0[1]) %>%
  mutate(interv = factor(interv, levels = c("N","D","DT","DTQ7","DTQ14"))) 


tmp <- gdf %>% filter(est_name == "mu_Re" & interv=="DTQ14")
R_c <- c(tmp$est, tmp$est_l, tmp$est_u)

Gs <- fplot1(gdf)
G <- gridExtra::grid.arrange(Gs[[1]] + labs(tag="a"),
                             Gs[[5]] + labs(tag="b"),
                             Gs[[6]] + labs(tag="c"),
                             nrow=1,ncol=3)

fplot_save(G, paste0("output/NPI_Fig2_0730_R_0_Rc.jpeg"))


# R_p
gdf <- df %>%
  mutate(interv = case_when(
    Xdetect_sen==0&Xprop.ascertain==0&Xquarant.days==0 ~ "N",
    Xdetect_sen==.95&Xprop.ascertain==0&Xquarant.days==0 ~ "D",
    Xdetect_sen==.95&Xprop.ascertain==.9&Xquarant.days==0 ~ "DT",
    Xdetect_sen==.95&Xprop.ascertain==.9&Xquarant.days==7 ~ "DTQ7",
    Xdetect_sen==.95&Xprop.ascertain==.9&Xquarant.days==14 ~ "DTQ14",
    T ~ "x"
  )) %>%
  filter(interv!="x" & Xnum.initial.cases == 100 & Xr0community==R_p[1]) %>%
  mutate(interv = factor(interv, levels = c("N","D","DT","DTQ7","DTQ14"))) 


tmp <- gdf %>% filter(est_name == "mu_Re" & interv=="DTQ14")
R_pc <- c(tmp$est, tmp$est_l, tmp$est_u)

Gs <- fplot1(gdf)
G <- gridExtra::grid.arrange(Gs[[1]] + labs(tag="a"),
                             Gs[[5]] + labs(tag="b"),
                             Gs[[6]] + labs(tag="c"),
                             nrow=1,ncol=3)

fplot_save(G, paste0("output/NPI_Fig2_0730_R_p_Rpc.jpeg"))


out <- matrix(c(R_0, R_p, R_c, R_pc), nrow = 4, ncol = 3, byrow = T) %>%
  round(2)%>%
  as.data.frame(row.names = c("R_0", "R_p", "R_c", "R_pc"))
colnames(out) <- c("mean", "CI_l", "CI_u")


