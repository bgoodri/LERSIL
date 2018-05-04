functions { // you can use these in R following `rstan::expose_stan_functions("foo.stan")`
  /*
    Fills in the elements of a coefficient matrix containing some mix of 
    totally free, free subject to a sign constraint, and fixed elements
    
    @param free_elements vector of unconstrained elements
    @param skeleton matrix of the same dimensions as the output whose elements are
      not_a_number(): if output element is totally free
      positive_infinity(): if output element is positive but otherwise free
      negative_infinity(): if output element is negative but otherwise free
      other: if output element is fixed to that number
    @return matrix of coefficients
  */
  matrix fill_matrix(vector free_elements, matrix skeleton) {
    int R = rows(skeleton);
    int C = cols(skeleton);
    matrix[R, C] out;
    int pos = 1;
    for (c in 1:C) for (r in 1:R) {
      real rc = skeleton[r, c];
      if (is_nan(rc)) { // totally free
        out[r,c] = free_elements[pos];
        pos += 1;
      }
      else if (is_inf(rc)) { // free subject to sign constraint
        out[r,c] = rc > 0 ? exp(free_elements[pos]) : -exp(free_elements[pos]);
        pos += 1;
      }
      else out[r,c] = skeleton[r, c]; // fixed, so do not bump pos
    }
    return out;
  }
  
  /*
   * This is a bug-free version of csr_to_dense_matrix and has the same arguments
   */
  matrix to_dense_matrix(int m, int n, vector w, int[] v, int[] u) {
    matrix[m, n] out = rep_matrix(0, m, n);
    int pos = 1;
    for (i in 1:m) {
      int start = u[i];
      int nnz = u[i + 1] - start;
      for (j in 1:nnz) {
        out[i, v[pos]] = w[pos];
        pos += 1;
      }
    }
    return out;
  }
}
data {
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA2#v=onepage&q=LISREL&f=false
  int<lower=0> p; // number of manifest response variables
  int<lower=0> q; // number of manifest predictors
  int<lower=0> m; // number of latent endogenous variables
  int<lower=0> n; // number of latent exogenous variables
  cov_matrix[p + q] S;     // sample covariance matrix among all manifest variables
  int<lower= p + q + 1> N; // number of observations that S is estimated from
  int<lower=0, upper=1> has_data;     // are the raw (centered) data on y and x available?
  vector[p + q] YX[has_data ? N : 0]; // if so, include them but concatenated together

  /* these are rate hyperparameters that go into the respective exponential priors */  
  vector<lower=0>[p] epsilon_sd_rate;
  vector<lower=0>[q] delta_sd_rate;

  /* sparse matrix representations of skeletons of coefficient matrices, 
     which is not that interesting but necessary because you cannot pass
     missing values into the data block of a Stan program from R */
  int<lower=0> len_w1;        // number of not totally free elements in Lambda_y
  vector[len_w1] w1;          // values of not totally free elements in Lambda_y
  int<lower=1> v1[len_w1];    // index  of not totally free elements in Lambda_y
  int<lower=1> u1[p + 1];     // index  of not totally free elements in Lambda_y
  int<lower=0> len_small_w1;           // number of "small" elements in Lambda_y
  int<lower=1> small_w1[len_small_w1]; // index  of "small" elements in Lambda_y
  /* N.B.: "small" elements are totally free but likely to be close to zero,
           otherwise you would tend to constrain their signs */
  real<lower=0> sd1; // prior scale for "small" elements in Lambda_y

  // same things but for Lambda_x
  int<lower=0> len_w2;
  vector[len_w2] w2;
  int<lower=1> v2[len_w2];
  int<lower=1> u2[q + 1];
  int<lower=0> len_small_w2;
  int<lower=1> small_w2[len_small_w2];
  real<lower=0> sd2;
  
  // same things but for Gamma
  int<lower=0> len_w3;
  vector[len_w3] w3;
  int<lower=1> v3[len_w3];
  int<lower=1> u3[m + 1];
  int<lower=0> len_small_w3;
  int<lower=1> small_w3[len_small_w3];
  real<lower=0> sd3;

  // same things but for B
  int<lower=0> len_w4;
  vector[len_w4] w4;
  int<lower=1> v4[len_w4];
  int<lower=1> u4[m + 1];
  int<lower=0> len_small_w4;
  int<lower=1> small_w4[len_small_w4];
  real<lower=0> sd4;
}
transformed data { // (re)construct skeleton matrices in Stan (not that interesting)
  matrix[p, m] Lambda_y_skeleton = to_dense_matrix(p, m, w1, v1, u1);
  matrix[q, n] Lambda_x_skeleton = to_dense_matrix(q, n, w2, v2, u2);
  matrix[m, n] Gamma_skeleton = to_dense_matrix(m, n, w3, v3, u3);
  matrix[m, m] B_skeleton = to_dense_matrix(m, m, w4, v4, u4);

  matrix[m, m] I = diag_matrix(rep_vector(1, m));

  int len_free1 = len_small_w1;
  int len_free2 = len_small_w2;
  int len_free3 = len_small_w3;
  int len_free4 = len_small_w4;
  
  real NaN = not_a_number();
  // replace totally free elements in Lambda_y_skeleton with NaN
  for (i in 1:p) {
    int pos = 1;
    for (j in 1:m) {
      if (is_inf(Lambda_y_skeleton[i,j])) len_free1 += 1;
      if (Lambda_y_skeleton[i,j] == 0 && (len_w1 == 0 || w1[pos] != 0)) {
        Lambda_y_skeleton[i,j] = NaN;
      } else pos += 1;
    }
  }

  // same thing but for Lambda_x_skeleton
  for (i in 1:q) {
    int pos = 1;
    for (j in 1:n) {
      if (is_inf(Lambda_x_skeleton[i,j])) len_free2 += 1;
      if (Lambda_x_skeleton[i,j] == 0 && (len_w2 == 0 || w2[pos] != 0)) {
        Lambda_x_skeleton[i,j] = NaN;
      } else pos += 1;
    }
  }
  
  // same thing but for Gamma_skeleton
  for (i in 1:m) {
    int pos = 1;
    for (j in 1:n) {
      if (is_inf(Gamma_skeleton[i,j])) len_free3 += 1;
      if (Gamma_skeleton[i,j] == 0 && (len_w3 == 0 || w3[pos] != 0)) {
        Gamma_skeleton[i,j] = NaN;
      } else pos += 1;
    }
  }

  // same thing but for B_skeleton
  for (i in 1:m) {
    int pos = 1;
    for (j in 1:m) {
      if (j >= i && B_skeleton[i,j] != 0) 
        reject("B_skeleton must be strictly lower triangular for this model");
      // relaxing this triangularity constraint is quite possible in Stan
      if (is_inf(B_skeleton[i,j])) len_free4 += 1;
      if (B_skeleton[i,j] == 0 && w4[pos] != 0) B_skeleton[i,j] = NaN;
      else pos += 1;
    }
  }
}
parameters {
  // free elements (possibly with inequality constraints) for coefficient matrices
  vector[len_free1] Lambda_y_free;
  vector[len_free2] Lambda_x_free;
  vector[len_free3] Gamma_free;
  vector[len_free4] B_free;

  // Psi = L_Psi * L_Psi' is the correlation matrix of zeta where L_Psi is a Cholesky factor
  cholesky_factor_corr[m] L_Psi;
  // Phi = L_Phi * L_Phi' is the correlation matrix of xi where L_Phi is a Cholesky factor
  cholesky_factor_corr[n] L_Phi;

  // you do not need to any of the above in the output but you should save everything below

  vector<lower=0>[p] epsilon_sd; // standard deviations of measurement error in y
  vector<lower=0>[q] delta_sd;   // standard deviations of measurement error in x
}
transformed parameters {
  matrix[p, m] Lambda_y = fill_matrix(Lambda_y_free, Lambda_y_skeleton);
  matrix[q, n] Lambda_x = fill_matrix(Lambda_x_free, Lambda_x_skeleton);
  matrix[m, n] Gamma = fill_matrix(Gamma_free, Gamma_skeleton);
  matrix[m, m] B = fill_matrix(B_free, B_skeleton);
  // multiply_lower_tri_self_transpose(T) = T * T' for a triangular matrix T but more efficient
  matrix[m, m] Psi = multiply_lower_tri_self_transpose(L_Psi);
  matrix[n, n] PHI = multiply_lower_tri_self_transpose(L_Phi);
  // Phi is a built-in function in Stan so use the symbol PHI for the covariance of xi
}
model { // N.B.: things declared in the model block do not get saved in the output, which is okay here
  matrix[p, m] Lambda_y_A = mdivide_right_tri_low(Lambda_y, I - B);     // = Lambda_y * (I - B)^{-1}
  matrix[n, q] Lambda_xt = transpose(Lambda_x);                         // copies so do it just once
  
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA3#v=onepage&q=LISREL&f=false
  matrix[p + q, p + q] Sigma;                                           // model covariance matrix
  // * is overloaded for matrix-matrix, matrix-vector, etc. multiplication
  matrix[p, q] top_right = Lambda_y_A * Gamma * PHI * Lambda_xt;        // top right block of Sigma
  // quad_form(V, W) = V' * W * V but is more computationally efficient
  Sigma[1:p, 1:p] = quad_form(quad_form(PHI, transpose(Gamma)) + Psi, transpose(Lambda_y_A));
  for (i in 1:p) Sigma[i,i] += square(epsilon_sd[i]);                   // top left block w/ error
  Sigma[1:p, (p + 1):(p + q)] = top_right;
  Sigma[(p + 1):(p + q), 1:p] = transpose(top_right);
  Sigma[(p + 1):(p + q), (p + 1):(p + q)] = quad_form(PHI, Lambda_xt);
  for (i in 1:q) {
    int ip = i + p;
    Sigma[ip, ip] += square(delta_sd[i]);                               // bottom right block w/ error
  }
  
  /* log-likelihood */
  target += wishart_lpdf(S | N - 1, 0.5 * (Sigma + transpose(Sigma))); // force Sigma to be symmetric
  
  /* prior densities (can modify) in log-units */
  target += normal_lpdf(Lambda_y_free[small_w1] | 0, sd1);
  target += normal_lpdf(Lambda_x_free[small_w2] | 0, sd2);
  target += normal_lpdf(Gamma_free[small_w3]    | 0, sd3);
  target += normal_lpdf(B_free[small_w4]        | 0, sd4);
  
  target += lkj_corr_cholesky_lpdf(L_Phi | 1); // implies correlation matrix of xi is uniform
  target += lkj_corr_cholesky_lpdf(L_Psi | 1); // implies correlation matrix of zeta is uniform

  target += exponential_lpdf(epsilon_sd | epsilon_sd_rate);
  target += exponential_lpdf(delta_sd | delta_sd_rate);
}
generated quantities { // these matrices are saved in the output but do not figure into the likelihood
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA34#v=onepage&q=LISREL&f=false
  matrix[m, m] A = mdivide_left_tri_low(I - B, I); // = (I - B)^{-1}
  matrix[m, n] total_xi_eta = A * Gamma;
  matrix[m, n] indirect_xi_eta = total_xi_eta - Gamma;
  matrix[m, m] total_eta_eta = A - I;
  matrix[m, m] indirect_eta_eta = total_eta_eta - B;
  matrix[p, m] total_eta_y = Lambda_y * A;
  matrix[p, m] indirect_eta_y = total_eta_y - Lambda_y;
  matrix[p, n] total_xi_y = total_eta_y * Gamma; // = indirect_xi_y since there is no direct effect
  
  matrix[N, has_data ? m : 0] eta;
  matrix[N, has_data ? n : 0] xi;
  
  if (has_data) { // all matrices defined in this local block are not saved in the output
    matrix[m, m] Psi_star = multiply_lower_tri_self_transpose(A * L_Psi);
    matrix[n, m] Pi_t = transpose(total_xi_eta);
    matrix[m, p] Lambda_yt = transpose(Lambda_y);
    matrix[n, q] Lambda_xt = transpose(Lambda_x);
    matrix[n, m] cov_eta_xi = PHI * Pi_t;
    matrix[q, m] cov_x_eta = Lambda_x * cov_eta_xi;
    matrix[n, p] cov_y_xi = cov_eta_xi * Lambda_yt;
    matrix[q, p] cov_y_x = Lambda_x * cov_y_xi;
    matrix[n, q] cov_x_xi = PHI * Lambda_xt;
    matrix[m, m] cov_eta = quad_form_sym(PHI, Pi_t) + Psi_star;
    
    matrix[p + q, p + q] top_left = append_row(
      append_col(quad_form_sym(cov_eta, Lambda_yt) +  diag_matrix(square(epsilon_sd)), 
                 transpose(cov_y_x)), 
      append_col(cov_y_x, quad_form_sym(PHI, Lambda_xt) + diag_matrix(square(delta_sd))) );
      
    matrix[m + n, p + q] corner = transpose(append_col(
      append_row(cov_y_xi, cov_eta * Lambda_yt),
      append_row(cov_x_xi, transpose(cov_x_eta))));
    
    matrix[m + n, m + n] bottom_right = append_row(
      append_col(cov_eta, transpose(cov_eta_xi)), append_col(cov_eta_xi, PHI) );
    
    matrix[p + q, p + q] precision = inverse_spd(top_left);
    matrix[m + n, m + n] L = cholesky_decompose(bottom_right - quad_form(precision, corner));
    matrix[m + n, p + q] beta = corner * precision;
    for (i in 1:N) {
      row_vector[m + n] latents = transpose(multi_normal_cholesky_rng(beta * YX[i], L));
      eta[i, ] = head(latents, m);
      xi[i,  ] = tail(latents, n);
    }
  }
} // end a with a completely blank line (not even whitespace)
