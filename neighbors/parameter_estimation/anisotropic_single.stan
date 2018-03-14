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
  //  ra [deg]
  //  dec [deg]
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
  real ra_rad[N];
  real dec_rad[N];
  for (i in 1:N) {
    ra_rad[i] = ra[i] * pi() / 180.;
    dec_rad[i] = dec[i] * pi() / 180.;
  }

  for(i in 1:N) {
    M[i,1,1] = -sin(ra_rad[i]);
    M[i,1,2] =  cos(ra_rad[i]);
    M[i,1,3] = 0.;
    M[i,2,1] = -sin(dec_rad[i])*cos(ra_rad[i]);
    M[i,2,2] = -sin(dec_rad[i])*sin(ra_rad[i]);
    M[i,2,3] = cos(dec_rad[i]);
  }
}

parameters {
  vector<lower=0>[N] d;                // true distances [kpc]
  vector[3] v0;                        // mean velocity  [km/s]
  vector<lower=0>[3] sigv;             // dispersion     [km/s]
  corr_matrix[3] Omega;                // correlation
}

transformed parameters {
  vector[3] a_model[N];

  for(i in 1:N) {
    a_model[i][1] = 1./d[i];
    a_model[i][2] = M[i,1] * v0 / d[i] / 4.74;
    a_model[i][3] = M[i,2] * v0 / d[i] / 4.74;
  }
}

model {

  // priors
  v0 ~ normal(0, 30);
  sigv ~ uniform(0., 50);
  Omega ~ lkj_corr(2);
  for(i in 1:N) {
    d[i] ~ constdens();
  }

  matrix[3,3] D[N];             // modified covariance matrix
  for(i in 1:N) {
    D[i] = C[i];
  }
  for(i in 1:N) {
    D[i,2:3,2:3] = D[i,2:3,2:3]
      + M[i] * quad_form_diag(Omega, sigv) * transpose(M[i]) / d[i]^2 / 4.74^2;
  }

  // likelihood
  for(i in 1:N) {
    a[i] ~ multi_normal(a_model[i], D[i]);
  }

}

generated quantities {
  // Generate data of star i from posterior of v0
  // vi[i] = [pmra, pmdec, rv]
  vector[3] vi[N];
  for(i in 1:N) {
    vi[i][1] = M[i,1] * v0 / d[i] / 4.74;
    vi[i][2] = M[i,2] * v0 / d[i] / 4.74;
    vi[i][3] = cos(dec_rad[i])*cos(ra_rad[i])*v0[1] + cos(dec_rad[i])*sin(ra_rad[i])*v0[2] + sin(dec_rad[i])*v0[3];
  }
}
