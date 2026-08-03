// pg_apfB.cpp
// Particle Gibbs (PG) with Auxiliary Particle Filter (APF) + Backward Sampling
// for theta_t1, and exact block (Chan-method) sampling of (theta_02, theta2)
// via a tridiagonal LDL^T factorization (the extended precision matrix is
// tridiagonal, so no sparse-matrix machinery is needed).
//
// Author of the original R code: Cleiton Moya de Almeida
// C++/Rcpp port for maximum performance.
//
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------

// log-sum-exp of a std::vector<double>
static inline double logsumexp(const std::vector<double>& x) {
    double m = x[0];
    for (std::size_t i = 1; i < x.size(); ++i) if (x[i] > m) m = x[i];
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += std::exp(x[i] - m);
    return m + std::log(s);
}

// Poisson log-likelihood contribution: y*theta - exp(theta)
static inline double log_p_yt(double yt, double theta) {
    double res = yt * theta - std::exp(theta);
    if (!std::isfinite(res)) res = -std::numeric_limits<double>::infinity();
    return res;
}

// log-density of N(x; mean, sd)
static inline double log_dnorm(double x, double mean, double sd) {
    static const double LOG_SQRT_2PI = 0.9189385332046727;
    double z = (x - mean) / sd;
    return -LOG_SQRT_2PI - std::log(sd) - 0.5 * z * z;
}

// Draw one index in [0, K) from (unnormalized) log-weights, using a single
// uniform draw and the cumulative distribution (inverse-CDF sampling).
static inline int sample_index_from_logw(const std::vector<double>& logw, int K) {
    double m = logw[0];
    for (int k = 1; k < K; ++k) if (logw[k] > m) m = logw[k];
    std::vector<double> w(K);
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { w[k] = std::exp(logw[k] - m); tot += w[k]; }
    double u = unif_rand() * tot;
    double cum = 0.0;
    for (int k = 0; k < K; ++k) {
        cum += w[k];
        if (u <= cum) return k;
    }
    return K - 1;
}

// Draw K indices (with replacement) in [0, K) from (unnormalized) log-weights.
// Uses the same normalized-weight vector for all K draws (as R's sample()).
static inline void sample_indices_from_logw(const std::vector<double>& logw, int K,
                                             std::vector<int>& out) {
    double m = logw[0];
    for (int k = 1; k < K; ++k) if (logw[k] > m) m = logw[k];
    std::vector<double> w(K);
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { w[k] = std::exp(logw[k] - m); tot += w[k]; }
    // cumulative sum
    std::vector<double> cum(K);
    cum[0] = w[0];
    for (int k = 1; k < K; ++k) cum[k] = cum[k - 1] + w[k];
    double total = cum[K - 1];
    for (int j = 0; j < K; ++j) {
        double u = unif_rand() * total;
        // binary search
        int lo = 0, hi = K - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (u <= cum[mid]) hi = mid; else lo = mid + 1;
        }
        out[j] = lo;
    }
}

// ---------------------------------------------------------------------
// Tridiagonal LDL^T factorization/solve (Chan method), 0-based, size n.
// A is symmetric tridiagonal: diagonal d[0..n-1], sub/super-diagonal e[0..n-2]
// (e[i] connects rows i and i+1). Produces L (unit lower-bidiagonal,
// L[i] = multiplier at (i, i-1), i = 1..n-1) and D (diagonal, D[0..n-1]).
// ---------------------------------------------------------------------
struct TriLDLT {
    std::vector<double> L; // size n, L[0] unused
    std::vector<double> D; // size n
};

static inline TriLDLT tri_ldlt_factor(const std::vector<double>& d,
                                       const std::vector<double>& e) {
    int n = (int) d.size();
    TriLDLT f;
    f.L.assign(n, 0.0);
    f.D.assign(n, 0.0);
    f.D[0] = d[0];
    for (int i = 1; i < n; ++i) {
        f.L[i] = e[i - 1] / f.D[i - 1];
        f.D[i] = d[i] - f.L[i] * f.L[i] * f.D[i - 1];
    }
    return f;
}

// Solve L^T x = w  (backward substitution; L upper-bidiagonal transpose)
static inline void tri_solve_Lt(const TriLDLT& f, const std::vector<double>& w,
                                 std::vector<double>& x) {
    int n = (int) f.D.size();
    x[n - 1] = w[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = w[i] - f.L[i + 1] * x[i + 1];
}

// Solve A x = b  (L D L^T x = b): forward, diag, backward(=tri_solve_Lt)
static inline void tri_solve_A(const TriLDLT& f, const std::vector<double>& b,
                                std::vector<double>& x) {
    int n = (int) f.D.size();
    std::vector<double> z(n), w(n);
    z[0] = b[0];
    for (int i = 1; i < n; ++i) z[i] = b[i] - f.L[i] * z[i - 1];
    for (int i = 0; i < n; ++i) w[i] = z[i] / f.D[i];
    tri_solve_Lt(f, w, x);
}

// ---------------------------------------------------------------------
// TEST-ONLY export: solve a tridiagonal SPD system A x = b given diag d
// and off-diag e, for validation against a dense reference in R.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List test_tridiag_solve(NumericVector d, NumericVector e, NumericVector b,
                          NumericVector w_for_Lt) {
    int n = d.size();
    std::vector<double> dv(d.begin(), d.end()), ev(e.begin(), e.end()), bv(b.begin(), b.end());
    TriLDLT fac = tri_ldlt_factor(dv, ev);
    std::vector<double> x(n);
    tri_solve_A(fac, bv, x);

    std::vector<double> wv(w_for_Lt.begin(), w_for_Lt.end());
    std::vector<double> xlt(n);
    tri_solve_Lt(fac, wv, xlt);

    return List::create(
        _["x_A"] = NumericVector(x.begin(), x.end()),
        _["D"] = NumericVector(fac.D.begin(), fac.D.end()),
        _["L"] = NumericVector(fac.L.begin(), fac.L.end()),
        _["x_Lt"] = NumericVector(xlt.begin(), xlt.end())
    );
}
// [[Rcpp::export]]
List test_chan_assembly(NumericVector theta1, double phi1, double phi2,
                          double mu_02, double sigma2_02, double theta_01) {
    int Tt = theta1.size();
    int Ttp1 = Tt + 1;
    std::vector<double> d(Ttp1), e(Ttp1 - 1), b(Ttp1);

    d[0] = phi2 + (1.0 / sigma2_02 + phi1);
    for (int i = 1; i < Tt; ++i) d[i] = 2.0 * phi2 + phi1;
    d[Tt] = phi2;
    for (int i = 0; i < Ttp1 - 1; ++i) e[i] = -phi2;

    b[0] = mu_02 / sigma2_02 + phi1 * (theta1[0] - theta_01);
    for (int i = 1; i < Tt; ++i) b[i] = (theta1[i] - theta1[i - 1]) * phi1;
    b[Tt] = 0.0;

    return List::create(_["d"] = NumericVector(d.begin(), d.end()),
                         _["e"] = NumericVector(e.begin(), e.end()),
                         _["b"] = NumericVector(b.begin(), b.end()));
}

// ---------------------------------------------------------------------
// Main PG-APF sampler
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List pg_apfB_cpp(NumericVector y,
                  int K, int N,
                  double mu_01, double sigma2_01,
                  double mu_02, double sigma2_02,
                  double nu_01, double eta_01,
                  double nu_02, double eta_02,
                  NumericVector theta1_init, NumericVector theta2_init,
                  double theta_01_init, double theta_02_init,
                  double W1_init, double W2_init,
                  bool verbose = true, int print_every = 100) {

    RNGScope rngScope; // sync with R's RNG / set.seed

    const int Tt = y.size();
    const int Ttp1 = Tt + 1;

    // ---- state (valores iniciais recebidos por parametro) ----
    std::vector<double> theta1(theta1_init.begin(), theta1_init.end());
    std::vector<double> theta2(theta2_init.begin(), theta2_init.end());
    double theta_01 = theta_01_init;
    double theta_02 = theta_02_init;
    double W1 = W1_init;
    double W2 = W2_init;

    // ---- histories (row-major access via NumericMatrix, R-compatible) ----
    NumericMatrix theta1_hist(N, Tt), theta2_hist(N, Tt);
    NumericVector theta_01_hist(N), theta_02_hist(N);
    NumericVector W1_hist(N), W2_hist(N);
    NumericVector ess_smc(N);

    // particle-system buffers (reused across iterations)
    std::vector<std::vector<double> > theta_1_k(Tt, std::vector<double>(K, 0.0));
    std::vector<std::vector<double> > log_w_tilde(Tt, std::vector<double>(K, 0.0));
    std::vector<double> log_lambda_k(K), theta_hat_t1_k(K), log_aux(K), log_w_t(K);
    std::vector<int> A(K);
    std::vector<double> log_bw(K);

    // Chan-method buffers (size Ttp1)
    std::vector<double> d(Ttp1), e(Ttp1 - 1), b(Ttp1), theta_hat(Ttp1);
    std::vector<double> u(Ttp1), w(Ttp1), xsol(Ttp1), draw2(Ttp1);
    std::vector<double> z(Tt - 1);

    for (int n = 0; n < N; ++n) {

        if (verbose && ((n + 1) % print_every == 0)) {
            Rcout << "Iteration " << (n + 1) << " / " << N << "\n";
        }

        // ---- Sample W1 ----
        double ssq1 = 0.0;
        {
            double prev1 = theta_01, prev2 = theta_02;
            for (int t = 0; t < Tt; ++t) {
                double dif1 = theta1[t] - prev1;
                double diffs1 = dif1 - prev2;
                ssq1 += diffs1 * diffs1;
                prev1 = theta1[t];
                prev2 = theta2[t];
            }
        }
        double nu_01_bar = nu_01 + Tt / 2.0;
        double eta_01_bar = eta_01 + 0.5 * ssq1;
        double phi1 = R::rgamma(nu_01_bar, 1.0 / eta_01_bar); // rate -> scale=1/rate
        W1 = 1.0 / phi1;
        double sd_W1 = std::sqrt(W1);

        // ---- Sample W2 ----
        double ssq2 = 0.0;
        {
            double prev = theta_02;
            for (int t = 0; t < Tt; ++t) {
                double diffs2 = theta2[t] - prev;
                ssq2 += diffs2 * diffs2;
                prev = theta2[t];
            }
        }
        double nu_02_bar = nu_02 + Tt / 2.0;
        double eta_02_bar = eta_02 + 0.5 * ssq2;
        double phi2 = R::rgamma(nu_02_bar, 1.0 / eta_02_bar);
        W2 = 1.0 / phi2;

        // ---- Sample theta_01 (conjugate Normal) ----
        double sigma2_01_bar = 1.0 / (1.0 / sigma2_01 + 1.0 / W1);
        double mu_01_bar = sigma2_01_bar * (mu_01 / sigma2_01 + (theta1[0] - theta_02) / W1);
        theta_01 = R::norm_rand() * std::sqrt(sigma2_01_bar) + mu_01_bar;

        // ---- Conditional SMC (APF) for theta_t1 ----
        // t = 0 (R: t=1)
        for (int k = 0; k < K; ++k) {
            theta_1_k[0][k] = R::norm_rand() * sd_W1 + (theta_01 + theta_02);
        }
        theta_1_k[0][K - 1] = theta1[0];
        {
            std::vector<double> log_w_1(K);
            for (int k = 0; k < K; ++k) log_w_1[k] = log_p_yt(y[0], theta_1_k[0][k]);
            double lse = logsumexp(log_w_1);
            for (int k = 0; k < K; ++k) log_w_tilde[0][k] = log_w_1[k] - lse;
        }

        // t = 1 .. Tt-1 (R: t = 2..Tt)
        for (int t = 1; t < Tt; ++t) {
            for (int k = 0; k < K; ++k) {
                theta_hat_t1_k[k] = theta_1_k[t - 1][k] + theta2[t - 1];
                log_lambda_k[k] = y[t] * theta_hat_t1_k[k] - std::exp(theta_hat_t1_k[k]);
            }
            for (int k = 0; k < K; ++k) log_aux[k] = log_w_tilde[t - 1][k] + log_lambda_k[k];

            sample_indices_from_logw(log_aux, K, A);
            A[K - 1] = K - 1; // reference path

            for (int k = 0; k < K; ++k) {
                double mean_k = theta_1_k[t - 1][A[k]] + theta2[t - 1];
                theta_1_k[t][k] = R::norm_rand() * sd_W1 + mean_k;
            }
            theta_1_k[t][K - 1] = theta1[t];

            for (int k = 0; k < K; ++k)
                log_w_t[k] = log_p_yt(y[t], theta_1_k[t][k]) - log_lambda_k[A[k]];
            double lse = logsumexp(log_w_t);
            for (int k = 0; k < K; ++k) log_w_tilde[t][k] = log_w_t[k] - lse;
        }

        // ESS at final time
        {
            double s2 = 0.0;
            for (int k = 0; k < K; ++k) s2 += std::exp(2.0 * log_w_tilde[Tt - 1][k]);
            ess_smc[n] = 1.0 / s2;
        }

        // ---- Backward sampling for theta1 ----
        {
            int k_final = sample_index_from_logw(log_w_tilde[Tt - 1], K);
            theta1[Tt - 1] = theta_1_k[Tt - 1][k_final];
        }
        for (int t = Tt - 2; t >= 0; --t) {
            for (int k = 0; k < K; ++k) {
                double mean_k = theta_1_k[t][k] + theta2[t];
                log_bw[k] = log_w_tilde[t][k] + log_dnorm(theta1[t + 1], mean_k, sd_W1);
            }
            int b_idx = sample_index_from_logw(log_bw, K);
            theta1[t] = theta_1_k[t][b_idx];
        }

        // ---- (theta_02, theta2) jointly via Chan method (tridiagonal LDL^T) ----
        for (int i = 0; i < Tt - 1; ++i) z[i] = theta1[i + 1] - theta1[i];

        d[0] = phi2 + (1.0 / sigma2_02 + phi1);
        for (int i = 1; i < Tt; ++i) d[i] = 2.0 * phi2 + phi1;
        d[Tt] = phi2;
        for (int i = 0; i < Ttp1 - 1; ++i) e[i] = -phi2;

        b[0] = mu_02 / sigma2_02 + phi1 * (theta1[0] - theta_01);
        for (int i = 1; i < Tt; ++i) b[i] = z[i - 1] * phi1;
        b[Tt] = 0.0;

        TriLDLT fac = tri_ldlt_factor(d, e);
        tri_solve_A(fac, b, theta_hat);

        for (int i = 0; i < Ttp1; ++i) u[i] = R::norm_rand();
        for (int i = 0; i < Ttp1; ++i) w[i] = u[i] / std::sqrt(fac.D[i]);
        tri_solve_Lt(fac, w, xsol);
        for (int i = 0; i < Ttp1; ++i) draw2[i] = theta_hat[i] + xsol[i];

        theta_02 = draw2[0];
        for (int i = 0; i < Tt; ++i) theta2[i] = draw2[i + 1];

        // ---- Store results ----
        theta_01_hist[n] = theta_01;
        theta_02_hist[n] = theta_02;
        W1_hist[n] = W1;
        W2_hist[n] = W2;
        for (int t = 0; t < Tt; ++t) {
            theta1_hist(n, t) = theta1[t];
            theta2_hist(n, t) = theta2[t];
        }
    }

    return List::create(
        _["theta1_hist"] = theta1_hist,
        _["theta2_hist"] = theta2_hist,
        _["theta_01_hist"] = theta_01_hist,
        _["theta_02_hist"] = theta_02_hist,
        _["W1_hist"] = W1_hist,
        _["W2_hist"] = W2_hist,
        _["ess_smc"] = ess_smc
    );
}
