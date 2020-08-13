### 20200610 NPI Taiwan figures
require(ggplot2)
require(dplyr)
require(tidyr)
require(gridExtra)
require(R0)
Sys.setlocale("LC_ALL","English")
theme_set(theme_classic(base_size = 15) %+replace% 
            theme(plot.tag = element_text(face ="bold", size=22),
                  legend.position = c(1,1), legend.justification = c(1,1),
                  legend.background = element_blank(), panel.grid.minor = element_blank()))

## Figure 2. (a) epidemic curve + Rt. 
## (b) distribution of the incubation period 
## (c) distribution of the delay from onset to testing/ diagnosis, (overall vs. indigenous cases)
## (d) distribution of the generation interval distribution
df0 <- readRDS("data/Taiwan_line_list.RDS")

df_test <- readr::read_csv("data/covid19_tw_specimen.csv") 
colnames(df_test) <- c("dates", "x_noti","x_quar","x_screen", "x_total")
df_test <- df_test %>%
  filter(!is.na(x_total))%>%
  mutate(dates = as.Date(dates, format = "%Y/%m/%d"),
         xx = case_when(
           x_total<=50 ~ "<50",
           x_total<=100~ "51~100",
           x_total<=500~ "101~500",
           x_total<1000~ "501~1000",
           T ~ ">1000"
         ))%>%
  mutate(xx = factor(xx, c("<50","51~100","101~500","501~1000",">1000")))


tmp <- df0 %>%
  mutate(S1 = if_else(is.na(S1L), Confirm, S1L+0.5*difftime(S1R, S1L) )) %>%
  mutate(dates = as.Date(S1), count= 1, 
         src = if_else(!is.na(Cluster) & Cluster == "Dunmu fleet", "naval", if_else(Imported == 1, "imported", "local"))) %>%
  group_by(dates, src) %>%
  summarise(count = sum(count)) %>% 
  pivot_wider(id_cols = dates, names_from = src, names_prefix = "Import_", 
              values_from = count, values_fill = list(count=0))
colnames(tmp) <- c("dates","imported","local","naval")

gdf0 <- data.frame(dates = seq.Date(as.Date("2020-01-11"), as.Date("2020-06-01"), by = "day")) %>%
  left_join(tmp) %>%
  left_join(df_test) %>%
  replace_na(replace = list(imported=0, local=0, naval=0, x_total=0, xx="<50"))%>%
  mutate(xt = dates-min(dates)+1)

source('./R/Est_Rt_WT.R')
pdf_gi_covid <- function(dt){
  dgamma(dt, shape = 1.07, scale = 5.5)
}
res <- est.Rt(t = gdf0$xt, I_local = gdf0$local+gdf0$naval, I_import = gdf0$imported, 
              pdf_gi = pdf_gi_covid, h = 7, nsim = 1000)


gdf0$Rt <- res$Rt
gdf0$Rt_l <- res$Rt_l
gdf0$Rt_u <- res$Rt_u


gdf1 <- gdf0 %>%
  pivot_longer(cols = c("imported","local","naval"), names_to = "src", values_to = "case") %>%
  mutate(src = factor(src, levels = c("imported","naval","local")))

colset <- c("imported"="#4682B4","local"="#FF7F50","naval"="#e1b137")
ys <- 8
G1a <- ggplot(gdf1)+
  geom_bar(aes(x=dates, y=case, fill=src), stat="identity", width = 1)+
  scale_fill_manual(values = colset, name="Source", 
                    label=c("Imported","Naval outbreak","Locally acquired"))+
  labs(x="Onset date", y="Daily cases")+
  
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x=dates, y=33, fill=xx))+
  scale_fill_viridis_d(name="No. of people tested")+

  geom_ribbon(data = gdf0, 
              aes(x = dates, ymin=ys*Rt_l, ymax=ys*Rt_u), 
              fill="#8e2525", alpha=.15)+
  geom_path(data = gdf0, 
            aes(x =dates, y=ys*Rt), color="#8e2525")+
  scale_y_continuous(sec.axis = sec_axis(~./ys, 
                                         name = expression(paste("Time-varying reproduction number ( ", R[t],")"))
                                         ))+
  
  geom_hline(yintercept = 1*ys, linetype=3, color="#8e2525")+
  coord_cartesian(ylim = c(0, 32))+
  guides(fill = guide_legend(order = 2),
         fill = guide_legend(order = 1))+
  theme(legend.position = c(.95,.9), legend.justification = c(1,1),
        legend.background = element_blank(),
        axis.title.y.right = element_text(color = "#8e2525"),
        axis.text.y.right = element_text(color = "#8e2525"),
        axis.line.y.right = element_line(color = "#8e2525"))
G1a

100*(sum(gdf0$imported)+sum(gdf0$local)+sum(gdf0$naval))/sum(gdf0$x_total)

binom.test(441, 68976)

## b 
est <- readRDS("mid/est_for_DTQ/Est_incu_case_DTQ158_v0729.RDS")
gdf <- est$Res$gdf

est2 <- readRDS("mid/est_for_DTQ/Est_si_case_DTQ158_v0729.RDS")
gdf2 <- est2$Res$gdf

est3 <- readRDS("mid/est_for_DTQ/Est_delay_case_DTQ158_v0729.RDS")
gdf3 <- est3$Res$gdf #%>% filter(xeval >= 0)

colset <- c("Serial interval"="#ff7f50", 
            "Incubation period"="#2b506e",
            "Onset-to-isolation interval"="#465927")
G2b <- ggplot(gdf) +
    geom_ribbon(data=gdf2,
                aes(x= xeval, ymin = q2.5, ymax = q97.5, 
                    fill="Serial interval",group="Serial interval"), alpha=.2)+
    geom_path(data=gdf2,
              aes(x=xeval, y=q50, color="Serial interval",group="Serial interval"), size=1)+
  
    geom_ribbon(data=gdf3,
                aes(x= xeval, ymin = q2.5, ymax = q97.5, 
                    fill="Onset-to-isolation interval",group="Onset-to-isolation interval"), alpha=.2)+
    geom_path(data=gdf3,
              aes(x=xeval, y=q50, color="Onset-to-isolation interval",group="Onset-to-isolation interval"), size=1)+
  
    
  
    geom_ribbon(aes(x= xeval, ymin = q2.5, ymax = q97.5, 
                    fill="Incubation period", group="Incubation period"), alpha=.2)+
    geom_path(aes(x=xeval, y=q50, color="Incubation period", group="Incubation period"), size=1)+
  
  
    coord_cartesian( xlim = c(0, 15), ylim=c(0,.25), expand = T)+
    scale_fill_manual(values = colset, name=NULL)+
    scale_color_manual(values = colset, name=NULL)+
    labs( x= "Incubation period, days", y = "Probability density")

G2b 



gdfb <- bind_rows(gdf %>% mutate(g = "Incubation period"),
                  gdf2 %>% mutate(g = "Serial interval"),
                  gdf3 %>% mutate(g = "Onset-to-isolation interval")) %>%
  mutate(g = factor(g, c("Incubation period", "Serial interval","Onset-to-isolation interval")))
if(F){ggplot(gdfb) +
  geom_ribbon(aes(x= xeval, ymin = q2.5, ymax = q97.5,fill=g, group=g), alpha=.2)+
  geom_path(aes(x=xeval, y=q50, color=g,group=g), size=1)+
  facet_wrap(~g,nrow = 3)+
  coord_cartesian( xlim = c(0, 15), ylim=c(0,.25), expand = T)+
  scale_fill_manual(values = colset, name=NULL)+
  scale_color_manual(values = colset, name=NULL)+
  labs( x= "Days", y = "Probability density")+
  theme(strip.background = element_blank(), strip.text =  element_blank())
}

est$Res$R2
qgamma(c(0.025,.5,0.975), shape = 2.84, scale = 1.93)

est2$Res$R2
qgamma(c(0.025,.5,0.975), shape = 1.31, scale = 5.22) -1 

est3$Res$R2
qgamma(c(0.025,.5,0.975), shape = 1.08, scale = 5.58) -1



adf_annot <- data.frame(g =  c("Onset-to-isolation interval", "Serial interval","Incubation period"),
                        est = c(5.02, 5.86, 5.5), 
                        est_l = c(-0, -0, 1.063246),
                        est_u = c(20.5113168, 21.5147570, 13.450666))%>%
  mutate(g = factor(g, c("Incubation period", "Serial interval","Onset-to-isolation interval")))
G2b <- ggplot(gdfb) +
  geom_violin(aes(x= g, y=xeval, violinwidth = q50*6,fill=g, group=g), color=NA,
              alpha=.2, stat="identity")+
  geom_pointrange(data = adf_annot, aes(x = g, y=est, ymin=est_l, ymax = est_u,color=g), size=.8)+
  scale_x_discrete(labels = c("Incubation\nperiod", 
                              "Serial\ninterval","Onset-to-isolation\ninterval"))+
  scale_fill_manual(values = colset, name=NULL, guide=F)+
  scale_color_manual(values = colset, name=NULL, guide=F)+
  coord_cartesian(ylim = c(0,25))+
  labs(y = "Days", x="")
G2b


## c
if(F){
est <- readRDS("mid/est_for_DTQ/Est_delay_case_DTQ2.RDS")
gdf1 <- est$Res$gdf %>% mutate(group = "Overall") %>%
  filter(xeval>0)
est$Res$R2
qgamma(c(0.025,.5,0.975), shape = 1.26, scale = 4.85) -1

est <- readRDS("mid/est_for_DTQ/Est_delay_local_case_DTQ2.RDS")
gdf2 <- est$Res$gdf %>% mutate(group = "Indigenous") %>%
  filter(xeval>0)
est$Res$R2
qgamma(c(0.025,.5,0.975), shape = 1.05, scale = 7.8) -1
}




date2period <- function(x){
  dcuts <- c("2020-01-01", "2020-02-15", "2020-03-01", "2020-03-15", "2020-04-01", "2020-04-15", "2020-05-01") %>% as.Date()
  tmp <- format(dcuts, "%b-%d")
  lcuts <- paste0(tmp[-length(tmp)], "~", tmp[-1])
  cut.Date(x, breaks = dcuts, right = F, labels = lcuts)
}

df <- readRDS("data/Taiwan_line_list.RDS")%>%
  filter(!is.na(S1L) & !is.na(S1R) & !is.na(Test) & (Test >= S1R)) %>%
  mutate(S1 = as.Date(S1L + 0.5*difftime(S1R, S1L))) %>%
  mutate(Xperiod = date2period(as.Date(S1)),
         delay = as.integer(as.Date(Test) - S1))%>%
  filter(!is.na(Xperiod))

G2c <- ggplot(df)+
  geom_boxplot(aes(y=delay, x=Xperiod), width=.2, color="#465927", size=.8, fill="#e3e6df")+
  labs(x="Onset date", y="Onset-to-isolation interval, days")+
  scale_x_discrete(labels = function(m)sub("~","~\n",m))+
  scale_y_continuous(breaks = seq(0,60,10))+
  theme(axis.text.x = element_text(angle=0, hjust=.5))


G2c


df %>%
  group_by(Xperiod)%>%
  summarise_at("delay", median)
  

## d : Serial interval
if(F){
est <- readRDS("mid/est_for_DTQ/Est_si_case_DTQ2.RDS")
gdf <- est$Res$gdf

G2d <- ggplot(gdf) +
  geom_ribbon(aes(x= xeval, ymin = q2.5, ymax = q97.5), fill="#2b506e", alpha=.3)+
  geom_path(aes(x=xeval, y=q50), color="#2b506e", size=1)+
  coord_cartesian( xlim = c(0, 20), ylim =c(0,0.25), expand = T)+
  labs( x= "Serial interval, days", y = "Probability density")

G2d 

est$Res$R2
qgamma(c(0.025,.5,0.975), shape = 1.07, scale = 5.5)

}



Gout <- arrangeGrob(
  G1a+ labs(tag="a"),
  G2b+ labs(tag="b"),
  G2c+ labs(tag="c"), 
  #G2d+ labs(tag="d"),
  layout_matrix = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,2,2,3,3),5,4,byrow = T))


ggsave(filename = "output/NPI_Fig1_v0810.jpeg",plot = Gout, 
       width = 24.7, height = 22, units = "cm", dpi = 300, scale=1.3)
