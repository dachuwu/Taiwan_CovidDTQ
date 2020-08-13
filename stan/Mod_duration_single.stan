data {
  int <lower=0> N; 
  real<lower=0> t0L[N]; 
  real<lower=0> t0R[N]; 
  real<lower=0> t1L[N]; 
  real<lower=0> t1R[N]; 
  real min_dur;
  int <lower=0> Neval; 
  real Xeval[Neval]; 
}
transformed data {
  real<lower=0> Xeval_til[Neval]; 
  for(i in 1:Neval) Xeval_til[i] = Xeval[i] - min_dur;
}


parameters {
  real<lower=0.01> shp; 
  real<lower=0.01> scl; 

  real<lower=0, upper=1> t0_til[N];
  real<lower=0, upper=1> t1_til[N];
}

transformed parameters {
  real<lower=0> SI_til[N]; 
  real<lower=0> t0e[N];
  real<lower=0> t1e[N];
  real<lower=0> iscl; 

  iscl = 1/scl;
  
  for (i in 1:N) {
    
    if(t0L[i]==t0R[i]){
      t0e[i] = t0L[i];
    }else{
      t0e[i] = t0L[i] + t0_til[i]*(t0R[i]-t0L[i]);
    }
    if(t1L[i]==t1R[i]){
      t1e[i] = t1L[i];
    }else{
      t1e[i] = t1L[i] + t1_til[i]*(t1R[i]-t1L[i]);
    }

    SI_til[i] = fmax(t1e[i] - t0e[i] - min_dur, 0.1);

  }
}

model {
  t0_til ~ uniform(0, 1);
  t1_til ~ uniform(0, 1);
  shp ~ exponential(0.001);//uniform(0.1,100);//normal(5,10);//
  scl ~ exponential(0.001);//uniform(0.1,10);//normal(5,10);//
  for (i in 1:N) {
    target += gamma_lpdf(SI_til[i]| shp, iscl);
  }
  
}
generated quantities {
  real SI[N]; 
  real g_mean; 
  real<lower=0> g_sd; 
  matrix[1,Neval] PDFeval;
  real<lower=0> pri_shp; 
  real<lower=0> pri_iscl;
  matrix[1,Neval] pri_PDFeval;
  
  for(i in 1:N)SI[i] = SI_til[i] + min_dur - 0.001;
  
  pri_shp = exponential_rng(0.001);//fmax(uniform_rng(0.1,100),0);//
  pri_iscl = 1/exponential_rng(0.001);//1/fmax(uniform_rng(0.1,10),0);//
  for(j in 1:Neval){
        pri_PDFeval[1,j] = exp(gamma_lpdf(Xeval_til[j] | pri_shp, pri_iscl));
  }
  g_mean = shp/iscl + min_dur ;
  g_sd = sqrt(shp/(iscl*iscl));
  for(j in 1:Neval){
        PDFeval[1,j] = exp(gamma_lpdf(Xeval_til[j] | shp, iscl));
  }
  

}

