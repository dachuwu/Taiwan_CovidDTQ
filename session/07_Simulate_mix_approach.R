library(DTQbp)
library(ggplot2)
library(purrr)
library(data.table)
library(future)
library(future.apply)
library(furrr)
require(dplyr)
require(gridExtra)

theme_set(theme_classic(base_size = 15) %+replace%
            theme(legend.position = c(0,1), legend.justification = c(0,1),
                  legend.background = element_blank(),
                  legend.direction = "vertical", 
                  panel.grid.minor = element_line(linetype = 2, colour = "grey90"), 
                  panel.grid.major = element_line(linetype = 2, colour = "grey90"),
                  plot.tag = element_text(face ="bold", size=26))
)

R_0 <<- c(2.5, 2.5, 2.5)
R_p <<- c(0.7, 0.518, 0.903)


###

fixed.para <- list(
  num.initial.cases = 100,
  cap_gen = Inf,
  #cap_max_days = 120,
  #cap_cases = 30E4,
  
  r0isolated = 0,
  disp.iso = 1, 
  
  relR.asym = .5, prop.asym = .15,
  
  incu_shape = 2.84, incu_scale = 1.93,
  delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158
  
  prop.presym = 0.524, si_omega = 2.689,
  disp.com = 19.196,
  
  quarant.retro.days = 4,
  save.case_data = T)

vary.para <- list(

  "(i) Do nothing (R0=2.5)" = list(r0community = R_0[1],
                                 detect_sen =0, prop.ascertain = 0, quarant.days = 0),
  
   "(ii) Detection + tracing + 14-day quarantine (R0=2.5)" = list(r0community = R_0[1],
                                                               detect_sen=.95, prop.ascertain = .9, quarant.days = 14),
  "(iii) Population-based interventions" = list(r0community = R_p[1],
                                         detect_sen=0, prop.ascertain = 0, quarant.days = 0),
  "(iv) Combining (ii) and (iii)" = list(r0community = R_p[1],
                              detect_sen=.95, prop.ascertain = .9, quarant.days = 14)
)



plan(multiprocess, workers = 10)
wrap.iter <- function(j,  cap.para = NULL){
  wrap.f <- lift(outbreak_model)
  gdf <- future_map_dfr(1:Niter, function(i){
   
    res <- wrap.f(c(fixed.para, vary.para[[j]], cap.para))
    df <- res$case_data
    df <- df[infector>0]
    m <- 3
    df[,`:=`(
      t = (onset %/% m)*m,
      cases = 1/m,
      i = i
    )][, .(cases=sum(cases)) ,by = .(t,i)][order(t)]
  })
  gdf0 <- gdf[, .(cases=median(cases)), by = .(t)][,i:=0]
  gdf <- rbind(gdf0, gdf)
  gdf[, scen := names(vary.para)[j]]
  
  return(gdf)
}
Niter <- 300
res1 <- lapply(1:2, wrap.iter,  
              cap.para = list(cap_max_days = 300, cap_cases = 200E4))
Niter <- 500
res2 <- lapply(3:length(vary.para), wrap.iter,  
              cap.para = list(cap_max_days = 300, cap_cases = 50E4))

plan(sequential)
saveRDS(list(res1=res1, res2=res2), "mid/SimEpiCurve_0810_v1.RDS")



tmp <- readRDS("mid/SimEpiCurve_0810_v1.RDS")

gdf1 <- rbindlist(tmp$res1)[t<=60]
gdf2 <- rbindlist(tmp$res2)
gdf <- rbind(gdf1, gdf2)


if(F){
colset <- c("grey20", "#1b4965", "#fb8b24", "#9a031e")
  
names(colset) <- names(vary.para) 
G1 <- ggplot(gdf[i>0])+
  geom_path(aes(x=t, y=cases, group=interaction(i,scen), color=scen), alpha=5/max(gdf$i))+
  geom_path(data=gdf[i==0], aes(x=t, y=cases, group=interaction(i,scen), color=scen), size=1.2)+
  # geom_point(data=gdf[i==0], aes(x=t, y=cases, group=interaction(i,scen), fill=scen),
  #            color="white", size=2, shape=21)+
  scale_color_manual(values = colset, name=NULL)+
  scale_fill_manual(values = colset, name=NULL)+
  coord_trans(y="pseudo_log", xlim = c(0,60), ylim = c(0, 1500), expand = F)+
  scale_y_continuous(breaks = c(0, 10, 50, seq(100, 1500, 300)))+
  labs(x = "days", y = "Daily cases")
G1
}



colset <- c("grey20", "#4682B4","#CD853F", "#800000")
colset <- c("grey20", "#1b4965", "#fb8b24", "#9a031e")
names(colset) <- names(vary.para) 
glabeler <- function(x){
  tmp <- list(
    "(i) Do nothing (R0=2.5)" ="(i)  Do nothing",

    "(ii) Detection + tracing + 14-day quarantine (R0=2.5)" = 
      bquote("(ii) Case-based interventions (DTQ14)"),
    "(iii) Population-based interventions" = 
      bquote("(iii)  Population-based interventions"),
    "(iv) Combining (ii) and (iii)" = 
      bquote( "(iv) Combining (ii) and (iii)")
  )
  tmp[x]
}

ggdf <- as.data.frame(gdf[i>0]) %>%
  group_by(t, scen)%>%
  summarise(cases_m = median(cases), cases_l = quantile(cases,.025),cases_u = quantile(cases,.975))

G1 <- ggplot(ggdf)+
  geom_ribbon(aes(x=t,ymin=cases_l, ymax=cases_u, fill=scen, group=scen), alpha=.25)+
  geom_path(aes(x=t,y=cases_m, color=scen, group=scen), size=1)+
  scale_color_manual(values = colset, name=NULL, labels= glabeler)+
  scale_fill_manual(values = colset, name=NULL, labels= glabeler)+
  coord_cartesian(xlim = c(0,60), ylim=c(0, 4500), expand = F)+
  labs(x = "Days after introduction", y = "Daily indicent cases (log scale)")
G1


G2 <- ggplot(ggdf)+
  geom_ribbon(aes(x=t,ymin=cases_l, ymax=cases_u, fill=scen, group=scen), alpha=.25)+
  geom_path(aes(x=t,y=cases_m, color=scen, group=scen), size=1)+
  scale_color_manual(values = colset, name=NULL, labels= glabeler)+
  scale_fill_manual(values = colset, name=NULL, labels= glabeler)+
  #coord_cartesian(xlim = c(0,65), ylim=c(0, 500), expand = F)+
  coord_trans(y="pseudo_log", xlim = c(0,60), ylim = c(0, 10000), expand = F)+
  scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000),minor_breaks=NULL)+
  labs(x = "Days after introduction", y = "Daily cases (log scale)")
G2


grid.arrange(G1, G2, nrow=1)


saveRDS(G2, "mid/Gplot_SimEpiCurve_v0810.RDS")



gunit <- 24.7

ggsave("output/NPI_Fig4_epicurves_tmp.jpeg", G2,
       width = 24.7, height = 16, units = "cm", dpi = 300, scale=1.6)


#

ggdf <- as.data.frame(gdf[i>0]) %>%
  group_by(scen,i)%>%
  mutate(cases = cumsum(cases)*3)%>%
  ungroup()%>%
  group_by(t, scen)%>%
  summarise(cases_m = median(cases), cases_l = quantile(cases,.025),cases_u = quantile(cases,.975))
G2 <- ggplot(ggdf)+
  geom_ribbon(aes(x=t,ymin=cases_l, ymax=cases_u, fill=scen, group=scen), alpha=.25)+
  geom_path(aes(x=t,y=cases_m, color=scen, group=scen), size=1)+
  scale_color_manual(values = colset, name=NULL, labels= glabeler)+
  scale_fill_manual(values = colset, name=NULL, labels= glabeler)+
  #coord_cartesian(xlim = c(0,65), ylim=c(0, 500), expand = F)+
  coord_trans(y="pseudo_log", xlim = c(0,60), ylim = c(0, 100000), expand = F)+
  scale_y_continuous(breaks = c(0, 10^c(1:5)))+
  labs(x = "Days after introduction", y = "Cumulative cases (log scale)")
G2

