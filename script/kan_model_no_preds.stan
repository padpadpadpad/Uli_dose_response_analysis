data{
    int N;                              // number of samples
    vector[N] conc;                     // antibiotic concentration (predictor)
    vector[N] comm;                     // dummy variable for comm/no comm (predictor)
    vector[N] fitness;                  // fitness (response variable)
    int Nnew;                           // length of predictions
    vector[Nnew] conc_new;              // new conc values to predict from
    vector[Nnew] comm_new;              // new comm values to predict from    
}
transformed data{
vector[N] comm_inter;
comm_inter = conc .* comm;
}
parameters{
    real<lower=0> sigma;         // standard deviation
    real b;             // steepness of increase
    real b_comm;        // effect of community on steepness of increase
    real c;             // lower asymptote
    real c_comm;        // effect of community on the lower asymptote c
    real d;             // intercept (shift curve up and down)
}
model{
    vector[N] mu;
    b ~ normal(0,10);
    b_comm ~ normal(0,10);
    c ~ normal(0,10);
    c_comm ~ normal(0,10);
    d ~ normal(0,10);
    sigma ~ cauchy(0,2);

    for(i in 1:N){
        mu[i] = (c+c_comm*comm[i])*exp(log(conc[i])*(b + b_comm*comm[i])) + d;
    }
    fitness ~ normal(mu, sigma);
}

