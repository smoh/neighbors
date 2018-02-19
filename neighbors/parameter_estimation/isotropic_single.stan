functions {

  // constant density prior on distance
  real constdens_lpdf(real r) {
    real ln_prob;
    if (r<1)
      ln_prob = 2*log(r);
    else
      ln_prob = negative_infinity();
    return ln_prob;
  }

}

data {
  int<lower=2> N;     // number of stars

  // Units:
  //  ra [rad]
  //  dec [rad]
  //  parallax [mas]
  //  pmra [mas/yr]
  //  pmdec [mas/yr]
  real ra[N];
  real dec[N];
  vector[3] a[N];        // parallax, pmra, pmdec
  cov_matrix[3] C[N];    // covariance matrix
}

transformed data {
  matrix[2,3] M[N];      // to equitorial rectangular coordinates
  for(i in 1:N) {
    M[i,1,1] = -sin(ra[i]);
    M[i,1,2] =  cos(ra[i]);
    M[i,1,3] = 0.;
    M[i,2,1] = -sin(dec[i])*cos(ra[i]);
    M[i,2,2] = -sin(dec[i])*sin(ra[i]);
    M[i,2,3] = cos(dec[i]);
  }
}

parameters {
  vector<lower=0>[N] d;     // true distances [kpc]
  vector[3] v0;             // mean velocity  [km/s]
  real<lower=0> sigv;       // dispersion     [km/s]
}

transformed parameters {
  cov_matrix[3] D[N];       // modified covariance matrix
  vector[3] tmp[N];

  for(i in 1:N) {
    D[i] = C[i];
    D[i,2,2] = D[i,2,2] + sigv^2 / d[i]^2 / 4.74^2;
    D[i,3,3] = D[i,3,3] + sigv^2 / d[i]^2 / 4.74^2;
  }
  //print(D)

  for(i in 1:N) {
    tmp[i][1] = 1./d[i];
    tmp[i][2] = M[i,1] * v0 / d[i] / 4.74;
    tmp[i][3] = M[i,2] * v0 / d[i] / 4.74;
  }
}

model {
  // priors
  for(i in 1:N) {
    d[i] ~ constdens();
  }
  // v0 is a vector -- stan knows normal is multidimensional
  // automagically?
  v0 ~ normal(0, 30);
  sigv ~ uniform(0., 50);


  // likelihood
  for(i in 1:N) {
    a[i] ~ multi_normal(tmp[i], D[i]);
  }

}

generated quantities {
  // Generate data of star i from posterior of v0
  // vi[i] = [pmra, pmdec, rv]
  vector[3] vi[N];
  for(i in 1:N) {
    vi[i][1] = M[i,1] * v0 / d[i] / 4.74;
    vi[i][2] = M[i,2] * v0 / d[i] / 4.74;
    vi[i][3] = cos(dec[i])*cos(ra[i])*v0[1] + cos(dec[i])*sin(ra[i])*v0[2] + sin(dec[i])*v0[3];
  }
}

