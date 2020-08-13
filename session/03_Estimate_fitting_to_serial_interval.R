# test single bp model realization
require(DTQbp)
require(ggplot2)
require(data.table)
require(magrittr)
load(file = "./mid/Observation_DTQ158_v0729.RData")
si_obs <<- si_obs_DTQ %>% sort()


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
if(F){
  plot.pffit <- function(fit){

  suplot <- function(i){
    X <- fit[[i]]
    sams0 <- X$prior.sams
    sams <- X$fitted.sams

    theme_set(theme_classic() %+replace%
                theme(legend.position = c(1,1), legend.justification = c(1,1)))
    colset <- c("fitted"="red4","prior"="grey70")

    gs <- lapply(1:NCOL(sams), function(j){
      xv <- colnames(sams)[j]

      gdf <- data.frame(value = c(sams[[xv]], sams0[[xv]]),
                        distribution = rep(c("fitted", "prior"), times = c(NROW(sams),NROW(sams0))))


      rr0 <- quantile(sams0[[xv]], c(.025,.975))
      rr1 <- quantile(sams[[xv]], c(.025,.975))
      rr <- c(min(rr0[1], rr1[1]), max(min(rr0[2], rr1[2])))
      rr <- extendrange(r = rr, f = .1)

      ggplot(gdf)+
        geom_density(aes(x=value, group=distribution, fill=distribution, color=distribution), alpha=.7,
                     n=2000/diff(rr))+
        coord_cartesian(xlim= rr)+
        labs(title=xv)+
        scale_fill_manual(values = colset, name="")+
        scale_color_manual(values = colset, name="")

    })


    measure.vect <- X$measure.sumvect[-c(1,2)]
    gs0 <- ggplot(data.table(x = 1:length(measure.vect),y = measure.vect)) +
      geom_path(aes(x=x, y=y))+
      labs(x="iteration", y="distance measure", title= paste0("Chain ",i))


    G <- gridExtra::arrangeGrob(grobs = c(list(gs0), gs), ncol = NCOL(sams)+1)
    invisible(G)
  }

  Gs <- lapply(1:length(fit), suplot)
  out <- gridExtra::arrangeGrob(grobs = Gs , nrow = length(fit))
  plot(out, nrow=length(fit))
  invisible(out)
}

}

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

output.wrap <- function(TAG){
  load(file = paste0(TAG,".RData"))

  G <- plot.pffit2(fit, fitted.paras, tts = paras.lab[names(fitted.paras)])

  gunit <- 24.7/NCOL(G)

  ggsave(filename = paste0(TAG,".jpeg"), plot = G,
         width = gunit*NCOL(G), height = gunit, units = "cm", dpi = 300, scale=1.5)
  write.csv(x = sum.fit(fit), file = paste0(TAG,".csv"))
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

output.wrap_si <- function(TAG){

  load(file = paste0(TAG,".RData"))
  tmp <- sum.fit(fit)
  XXprop.presym <- tmp$mean[tmp$variable == "prop.presym"]
  XXsi_omega <- tmp$mean[tmp$variable == "si_omega"]

  Gs <- plot.pffit(fit, fitted.paras, tts = paras.lab[names(fitted.paras)])

  Niter<-500
  sim.si <- lapply(1:Niter, function(i){
    dat <- DTQbp::pf_outbreak_model(num.initial.cases = fixed.paras$num.initial.cases,
                                    cap_gen = fixed.paras$cap_gen,
                                    r0isolated = fixed.paras$r0isolated,
                                    disp.iso = fixed.paras$disp.iso,

                                    incu_shape = fixed.paras$incu_shape,
                                    incu_scale = fixed.paras$incu_scale,
                                    delay_shape = fixed.paras$delay_shape,
                                    delay_scale = fixed.paras$delay_scale,
                                    delay_off = fixed.paras$delay_off,  #0729
                                    prop.presym = XXprop.presym, si_omega = XXsi_omega,

                                    detect_sen=fixed.paras$detect_sen,
                                    prop.ascertain = fixed.paras$prop.ascertain,
                                    quarant.days = fixed.paras$quarant.days,
                                    quarant.retro.days = fixed.paras$quarant.retro.days,

                                    relR.asym = fixed.paras$relR.asym,
                                    prop.asym = fixed.paras$prop.asym,
                                    r0community = fixed.paras$r0community,
                                    disp.com = fixed.paras$disp.com )

    Sfrom <- dat$onset[dat$infector[dat$infector >0]]
    Sto <- dat$onset[dat$infector >0]
    xx <- Sto - Sfrom

    dat <- dat[dat$infector >0]
    asinfector <- dat[asym == T, caseid]
    incld <- which(!dat$infector %in% asinfector)
    tmp <- sort(xx[incld])
    tmp <- density(tmp, adjust = 3)
    out <- data.frame(x= tmp$x, y=tmp$y)
    out$i <- i
    return(out)
  }) %>% rbindlist()

  ys <- 7/quantile(sim.si$y, .95)
  G0 <- ggplot(data.frame(x=si_obs))+
    geom_bar(aes(x=x, fill="Observed"), width = .5)+
    scale_fill_manual(values = c("Observed" = "grey40"), name = NULL)+
    scale_x_continuous(breaks = seq(0,25,5))+
    geom_path(data=sim.si, aes(x=x, y=y*ys, group=i, color="Simulated"),alpha = 5/Niter)+
    scale_color_manual(values = c("Simulated" = "#8e2525"), name = NULL)+
    guides(color = guide_legend(override.aes = list(alpha = .7)))+
    scale_y_continuous(breaks=seq(0,10,2),sec.axis = sec_axis( ~ ./ys, name = "Density"))+
    coord_cartesian(xlim=c(-1, 25), ylim = c(0, 9))+
    labs(x = "Serial interval", y="Frequency")+
    theme(legend.box.margin = margin(0, 5, 0, 0),
          legend.box.spacing = unit(0, "npc"),
          legend.box = "horizontal",
          legend.spacing = unit(0, "npc"),
          axis.title.y.right = element_text(color = "#8e2525"),
          axis.text.y.right = element_text(color = "#8e2525"),
          axis.line.y.right = element_line(color = "#8e2525"))

  G0

  GG <- gridExtra::arrangeGrob(grobs = list(Gs[[3]] + labs(tag="a"),
                                            Gs[[2]] + labs(tag="b"),
                                            Gs[[1]] + labs(tag="c"),
                                            G0 + labs(tag="d")),
                               ncol = 2, nrow = 2)
  #plot(GG)
  gunit <- 24.7/NCOL(GG)

  ggsave(filename = paste0(TAG,".jpeg"), plot = GG,
         width = gunit*NCOL(GG), height = gunit*NROW(GG), units = "cm",
         dpi = 300, scale=1)

  write.csv(x = sum.fit(fit), file = paste0(TAG,".csv"))
}

###
model <- pf_outbreak_model
measure_si <- function(dat){

  Sfrom <- dat$onset[dat$infector[dat$infector >0]]
  Sto <- dat$onset[dat$infector >0]
  xx <- Sto - Sfrom

  dat <- dat[dat$infector >0]
  asinfector <- dat[asym == T, caseid]
  incld <- which(!dat$infector %in% asinfector)
  xx <- xx[incld]

  if(length(xx)>0){
    tmp <- purrr::quietly(ks.test)(xx, si_obs)
    tmp$result$statistic[["D"]]
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



### fit to  serial interval
if(F){
  dat <- DTQbp::pf_outbreak_model( num.initial.cases = 50,
                                   cap_gen = 3,
                                   r0isolated = 0,
                                   disp.iso = 1,

                                   incu_shape = 2.84, incu_scale = 1.93,
                                   delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # overall delay

                                   prop.ascertain = 0.9, detect_sen=0.95,
                                   quarant.days = 14, quarant.retro.days = 4,
                                   prop.asym = 0,
                                   relR.asym = 1,
                                   disp.com = 1,
                                   r0community = 2.5,
                                   prop.presym = .4, si_omega = 2
  )

  Sfrom <- dat$onset[dat$infector[dat$infector >0]]
  Sto <- dat$onset[dat$infector >0]
  xx <- Sto - Sfrom

  dat <- dat[dat$infector >0]
  asinfector <- dat[asym == T, caseid]
  incld <- which(!dat$infector %in% asinfector)
  si <- xx[incld]

  ggplot(data.frame(fitted = si))+
    geom_density(aes(x=fitted), fill="red4")+
    geom_density(data = data.frame(obs = si_obs),
                 aes(x=obs), fill="grey50", alpha=.5)


}
if(T){
  fitted.paras <- list(
   si_omega = c(0.01, 5), prop.presym = c(0.01, 0.99)

  )
  fixed.paras <- list(
    num.initial.cases = 100,
    cap_gen = 3,
    r0isolated = 0,
    disp.iso = 1,

    incu_shape = 2.84, incu_scale = 1.93,
    delay_shape = 1.08, delay_scale = 5.58,  delay_off = -1,  # 0729 DTQ158

    prop.ascertain = 0.9, detect_sen=0.95,
    quarant.days = 14, quarant.retro.days = 4,
    prop.asym = 0.15,
    relR.asym = .5,
    disp.com = 20,
    r0community = 2.5
  )
  fit <- pfilter_fitXX(model, fitted.paras, fixed.paras, measure_si,
                       max.iter = 20, n.sam = 1000, eps = 0.1, n.chain = 4)
  sum.fit(fit)
  save(fit, fitted.paras, fixed.paras, measure_si, pfilter_fit_chain, pfilter_fitXX,
       file = "PFfit_DTQ158si_0729v1.RData")
}

if(T){
  TAG <- "PFfit_DTQ158si_0729v1"
  output.wrap_si(TAG)
}



if(T){
  tmp <- sum.fit(fit)
  XXprop.presym <- tmp$mean[tmp$variable == "prop.presym"]
  XXsi_omega <- tmp$mean[tmp$variable == "si_omega"]
  source('C:/Users/user/Desktop/GitDaChu/DTQbp/tests/test_pf_realization_cs_v0729.R')
}

