data{
    int N;                              // number of samples
    vector[N] conc;                     // antibiotic concentration (predictor)
    vector[N] comm;                     // dummy variable for comm/no comm (predictor)
    vector[N] fitness;                  // fitness (response variable)
    int Nnew;                           // length of predictions
    vector[Nnew] conc_new;              // new conc values to predict from
    vector[Nnew] comm_new;              // new comm values to predict from    
}
parameters{
    real<lower=0> sigma;         // standard deviation
    real b;             // steepness of increase
    real c;             // lower asymptote
    real d;             // intercept (shift curve up and down)
    real d_comm;        // effect of community (Y/N) on the intercept
}
model{
    vector[N] mu;
    b ~ normal(0,10);
    c ~ normal(0,10);
    d ~ normal(0,10);
    d_comm ~ normal(0,10);
    sigma ~ cauchy(0,2);

    mu = c*exp(b*log(conc)) + d + d_comm*comm;
    fitness ~ normal(mu, sigma);
}

