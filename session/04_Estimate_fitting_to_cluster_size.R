# test single bp model realization
require(DTQbp)
require(ggplot2)
require(data.table)
require(magrittr)
load(file = "./mid/Observation_DTQ158_v0729.RData")
si_obs <<- si_obs_DTQ %>% sort()
cs_obs <<- csize_obsDTQ_l %>% sort()
cs_obs <<- cs_obs[cs_obs > 1]

XXprop.presym <- 0.524
XXsi_omega <- 2.689

###
pfilter_fit_chain <- function(model = NULL,
                          fitted.paras = NULL, fixed.paras = NULL,
                          measure = NULL,
                          max.iter = 200, n.sam = 1000, eps = .1){

  sams_l <- purrr::map(fitted.paras, 1)
  sams_u <- purrr::map(fitted.paras, 2)

  sams <- pfilter_lhs(fitted.paras, n = n.sam)
  sams0 <- sams

  ksDs.med_vect <- c(0, 0)
  endloop <- F
  iter <- 1
  future::plan(future::multiprocess, workers = max(future::availableCores() - 2, 1))
  while(!endloop & iter <= max.iter){

    wrap.f <- function(x){
      #xx <- purrr::map_dbl(1:10, function(i) purrr::lift_dl(model)(x) %>% measure() )
      xx <- c()
      for(i in 1:10){
        xx <- c(xx, purrr::lift_dl(model)(x) %>% measure())
      }
      # purrr::lift_dl(model)(x) %>% measure()
      median(xx)
    }

    paraLS <- lapply(1:n.sam,  function(i){
      c(fixed.paras, as.list(sams[i,]))
    })

    ksDs <- future.apply::future_sapply(paraLS, wrap.f)

    cat("\n iteration ", iter, ".")
    iter <- iter +1

    #sel.id <- which(ksDs < max(ksDs,na.rm = T)*.5) #
    #ksDs <- ksDs[sel.id]

    ksDs.med <- median(ksDs)
    ksDs.med_1 <- ksDs.med_vect[length(ksDs.med_vect)]
    ksDs.med_2 <- ksDs.med_vect[length(ksDs.med_vect)-1]
    ksDs.med_vect <- c(ksDs.med_vect, ksDs.med)

    pmar <- .1
    endloop <- ksDs.med > ksDs.med_1*(1-pmar) & ksDs.med < ksDs.med_1*(1+pmar) &
      ksDs.med > ksDs.med_2*(1-pmar) & ksDs.med < ksDs.med_2*(1+pmar) & ksDs.med <= eps
    if(endloop){
      cat("\n Converged.")
      break
    }

    #sams <- sams[sel.id, ]
    #new.sams <- sams[sample.int(NROW(sams), n.sam-NROW(sams),prob = 1/ksDs, replace = T),]
    #new.sams <- new.sams[, lapply(.SD, function(x){ x*runif(length(x), min = 0.95, 1.05) })]
    #sams <- rbind(sams, new.sams)
    iwts <- pmax(ksDs, 0.0001)^2
    new.sams <- sams[sample.int(NROW(sams), n.sam, prob = 1/iwts, replace = T),]
    sams <- new.sams[, lapply(.SD, function(x){ x*runif(length(x), min = 0.9, 1.1) })]


    for(m in colnames(sams)) sams[[m]] <- pmin(pmax(sams[[m]], sams_l[[m]]), sams_u[[m]])

  }
  future::plan(future::sequential)

  return(list(
    fitted.sams = sams,
    prior.sams = sams0,
    measure.sumvect = ksDs.med_vect, measure.stats = ksDs))
}

pfilter_fitXX <- function(model = NULL,
                          fitted.paras = NULL, fixed.paras = NULL,
                          measure = NULL,
                          max.iter = 200, n.sam = 1000, eps = .1, n.chain = 1){

  res <- lapply(1:n.chain, function(i){
    cat("\n Sampling chain ", i, " ===============")
    cfit <- pfilter_fit_chain(model = model,
                              fitted.paras = fitted.paras,
                              fixed.paras = fixed.paras,
                              measure = measure,
                              max.iter = max.iter, n.sam = n.sam, eps = eps)
    cat("\n Completed. ")
    return(cfit)
    })
  names(res) <- paste0("chain.",1:n.chain)
  return(res)
}

###

plot.pffit <- function(fit = NULL, fitted.paras = NULL, tts = NULL){

  if(is.null(tts)) tts <- names(fitted.paras)
  tts <- c("Median distance", tts)

  theme_set(theme_classic(base_size = 15) %+replace%
              theme(legend.position = c(1,1), legend.justification = c(1,1),
                    legend.background = element_blank(),
                    legend.direction = "horizontal", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                    plot.tag = element_text(face ="bold", size=26))
  )
  gdf <- lapply(1:length(fit), function(i){
    X <- fit[[i]]
    dpos <- melt.data.table(X$fitted.sams, measure.vars = colnames(X$fitted.sams))[,`:=`(
      type= "posterior",
      chain = i)]
    dpri <- melt.data.table(X$prior.sams, measure.vars = colnames(X$fitted.sams))[,`:=`(
      type= "prior",
      chain = i)]
    rbind(dpos, dpri)
  }) %>% rbindlist()

  gdf0 <- lapply(1:length(fit), function(i){
    data.table(x = fit[[i]]$measure.sumvect[-c(1,2)] , chain = i)[, iteration := 1:.N]
  })%>% rbindlist()

  colset <-  colorRampPalette(c( "khaki2","red4"))(max(gdf$chain))
  G0 <- ggplot(gdf0) +
    geom_path(aes(x=iteration, y=x, group = as.character(chain), color = as.character(chain)))+
    labs(x="Iteration", y="Median KS statistics", tag = "a")+
    scale_y_continuous(breaks =seq(0,1,.05))+
    scale_color_manual(values = colset, name="chain")

  xvs<-unique(gdf$variable)
  G <- lapply(1:length(xvs), function(m){
    xv <- xvs[m]
    sgdf <- gdf[variable == xv & type == "posterior"]
    ggplot(sgdf)+
      geom_density(aes(x = value, group=as.character(chain),
                       fill = as.character(chain)), color = "grey70", alpha = .6, adjust = 2)+
      geom_density(data = gdf[variable == xv & type == "prior" & chain == 1],
                   aes(x = value), fill = NA, color = "grey50", linetype=2,adjust = 1/5)+
      coord_cartesian(xlim = fitted.paras[[xv]], expand = T)+
      scale_y_continuous(expand = expansion(mult = c(0.05,0.15)))+
      scale_fill_manual(values = colset, name="chain")+
      labs(x = tts[m+1], y="Density", tag = letters[m+1])
  })

  #out <- gridExtra::arrangeGrob(grobs = c(list(G0), G) , ncol = 1 + length(G))

  return(c(list(G0), G))
}


sum.fit <- function(fit = NULL, nd = 3){


  sams <- lapply(fit, function(X) X$fitted.sams) %>% rbindlist()

  sdf <- lapply(1:NCOL(sams), function(i){
    xx <- sams[[i]]
    quants <- quantile(xx, c(0.025,  0.5,  0.975)) %>% round(nd)
    data.frame(
      variable = colnames(sams)[i],
      mean = mean(xx)%>% round(nd), sd = sd(xx)%>% round(nd), q2.5 = quants[1], q50 = quants[2], q97.5 = quants[3])
  }) %>% rbindlist()
  return(sdf)
}


output.wrap_cs <- function(TAG){
  load(file = paste0(TAG,".RData"))
  tmp <- sum.fit(fit)
  XXr0community <- tmp$mean[tmp$variable == "r0community"]
  XXdisp.com <- tmp$mean[tmp$variable == "disp.com"]

  Gs <- plot.pffit(fit, fitted.paras, tts = paras.lab[names(fitted.paras)])

  Niter <- 500
  sim.cs <- lapply(1:Niter, function(i){
      dat <- DTQbp::pf_outbreak_model(num.initial.cases = fixed.paras$num.initial.cases,
                                      cap_gen = fixed.paras$cap_gen,
                                      r0isolated = fixed.paras$r0isolated,
                                      disp.iso = fixed.paras$disp.iso,

                                      incu_shape = fixed.paras$incu_shape,
                                      incu_scale = fixed.paras$incu_scale,
                                      delay_shape = fixed.paras$delay_shape,
                                      delay_scale = fixed.paras$delay_scale,
                                      delay_off = fixed.paras$delay_off,  #0729
                                      prop.presym = fixed.paras$prop.presym,
                                      si_omega = fixed.paras$si_omega,

                                      detect_sen=fixed.paras$detect_sen,
                                      prop.ascertain = fixed.paras$prop.ascertain,
                                      quarant.days = fixed.paras$quarant.days,
                                      quarant.retro.days = fixed.paras$quarant.retro.days,

                                      relR.asym = fixed.paras$relR.asym,
                                      prop.asym = fixed.paras$prop.asym,
                                      r0community = XXr0community, disp.com = XXdisp.com )
      tmp <- sort(as.integer(table(dat$cluster_obs)))
      tmp <- density(tmp[tmp>1],from=0,to=20, adjust = 3)
      out <- data.frame(x= tmp$x, y=tmp$y)
      out$i <- i
      return(out)
    }) %>% rbindlist()

  ys <- 10/quantile(sim.cs$y, .95)
  G0 <- ggplot(data.frame(x=cs_obs))+
    geom_bar(aes(x=x, fill="Observed"), width = .5)+
    scale_fill_manual(values = c("Observed" = "grey40"), name = NULL)+
    scale_x_continuous(breaks = 1:15)+
    geom_path(data=sim.cs, aes(x=x, y=y*ys, group=i, color="Simulated"),alpha = 8/Niter)+
    scale_color_manual(values = c("Simulated" = "#8e2525"), name = NULL)+
    guides(color = guide_legend(override.aes = list(alpha = .7)))+
    scale_y_continuous(breaks=seq(0,15,3),sec.axis = sec_axis( ~ ./ys, name = "Density"))+
    coord_cartesian(xlim=c(2, 12))+
    labs(x = "Cluster size", y="Frequency")+
    theme(legend.box.margin = margin(0, 5, 0, 0),
          legend.box.spacing = unit(0, "npc"),
          legend.box = "horizontal",
          legend.spacing = unit(0, "npc"),
          axis.title.y.right = element_text(color = "#8e2525"),
          axis.text.y.right = element_text(color = "#8e2525"),
          axis.line.y.right = element_line(color = "#8e2525"))

  G0

  GG <- gridExtra::arrangeGrob(grobs = list(Gs[[3]] + labs(tag="a",
                                                           x=bquote(R[p])),
                               Gs[[2]] + labs(tag="b"),
                               Gs[[1]] + labs(tag="c"),
                               G0 + labs(tag="d")),
                               ncol = 2, nrow = 2)
  #plot(GG)
  gunit <- 24.7/NCOL(GG)

  ggsave(filename = paste0(TAG,".jpeg"), plot = GG,
         width = gunit*NCOL(GG), height = gunit*NROW(GG), units = "cm", dpi = 300, scale=1)

  write.csv(x = sum.fit(fit), file = paste0(TAG,".csv"))
}

###
model <- pf_outbreak_model
measure_cs <- function(dat){

  srcid <- dat$caseid[dat$infector==0]
  xx<- dat$cluster_obs[dat$cluster_obs%in%srcid]
  if(length(xx)>0){
    sxx <- as.vector(table(xx))
    sxx <- sxx[sxx > 1]
    tmp <- purrr::quietly(ks.test)(sxx, cs_obs)
    tmp$result$statistic[["D"]]
    # pf_ks.test(as.vector(table(xx)), rep(J, N))
  } else {
    1
  }
}

paras.lab <- c(
  "si_omega" = "S.D. of the generation interval",
  "prop.presym" = "Prop. of pre-symptomatic transmission",
  "disp.com" = "Dispersion parameter",
  "r0community" = "Basic reproduction number",
  "prop.asym" = "Prob. of asymptomatic infection"
)



### fit to cluster size
if(T){
  fitted.paras <- list(
    disp.com = c(0.001, 50),
    r0community = c(0.1, 3)
  )
  fixed.paras <- list(
    num.initial.cases = 200,
    cap_gen = 5,
    r0isolated = 0,
    disp.iso = 1,

    incu_shape = 2.84, incu_scale = 1.93,
    delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158

    prop.presym = XXprop.presym, si_omega = XXsi_omega, ######

    detect_sen=0.95, prop.ascertain = 0.9,
    quarant.days = 14, quarant.retro.days = 4,


    relR.asym = .5, prop.asym = 0.15
  )
  fit <- pfilter_fitXX(model, fitted.paras, fixed.paras, measure_cs,
                       max.iter = 20, n.sam = 1000, eps = 0.1, n.chain = 4)

  tmp <- sum.fit(fit)
  XXr0community <- tmp$mean[tmp$variable == "r0community"]
  XXdisp.com <- tmp$mean[tmp$variable == "disp.com"]
  save(fit, fitted.paras, fixed.paras, measure_cs, pfilter_fit_chain, pfilter_fitXX,
       file = "PFfit_DTQ158csl_gt2_v0729_3.RData")

}
if(T){
  TAG <- "PFfit_DTQ158csl_gt2_v0729_3"
  output.wrap_cs(TAG)
}
