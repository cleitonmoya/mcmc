// pg_asB.cpp
// Particle Gibbs with Ancestral Sampling (PG-AS): Bootstrap Particle Filter
// with Ancestral Sampling for theta_t1, reconstructed by simple Backward
// Tracking (genealogy lookup, no Backward Simulation), and exact block
// (Chan-method) sampling of (theta_02, theta2) via a tridiagonal LDL^T
// factorization (the extended precision matrix is tridiagonal, so no
// sparse-matrix machinery is needed).
//
// PERFORMANCE NOTES (optimized version, same algorithm/results as before):
//  - theta_1_k is stored as a flat row-major buffer (index t*K + k) instead
//    of std::vector<std::vector<double>>, removing one level of pointer
//    indirection and improving cache locality.
//  - Unlike PG-APF (which needs the full T x K weight history for Backward
//    Simulation), PG-AS only ever reads log_w_tilde at t-1 (to build the
//    resampling/ancestral-sampling probabilities) and at the final time T
//    (to seed the backward tracking). So the weights are kept as two K-sized
//    buffers (log_w_prev, log_w_curr) swapped every t, instead of a T x K
//    matrix -- an O(K) memory footprint instead of O(K*T), specific to this
//    algorithm's structure.
//  - The resampling helpers (sample_index_from_logw / sample_indices_from_logw)
//    now take pre-allocated scratch buffers instead of allocating a new
//    std::vector internally on every call (they were being called T*N times).
//  - sample_indices_from_logw now performs SYSTEMATIC resampling (Kitagawa,
//    1996; Doucet, de Freitas & Gordon, 2001) instead of multinomial
//    resampling for the K-1 non-reference particles: a single unif_rand()
//    draw generates all K-1 indices via one forward sweep through the
//    cumulative weights, instead of K-1 independent unif_rand() draws + K-1
//    binary searches. This is a standard, unbiased resampling scheme with
//    LOWER variance than multinomial resampling (it is commonly preferred
//    in the SMC literature, not merely a speed shortcut). IMPORTANT: because
//    it consumes the RNG stream differently, results are NOT bit-identical
//    to the multinomial version for the same seed -- statistical validation
//    (posterior mean/ESS equivalence) was done separately, not bitwise
//    reproducibility. The single-index Ancestral Sampling draw
//    (sample_index_from_logw) is unchanged (it was already a single
//    unif_rand() draw, so there is nothing to gain there).
//  - No other change to the RNG call sequence, so aside from the resampling
//    scheme above, results are bit-identical to the original implementation
//    for the same seed.
//
// Author of the original R code: Cleiton Moya de Almeida
// C++/Rcpp port for maximum performance.
//
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <utility>
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
// logw is a raw pointer, so it can point directly into a scratch weight
// buffer without a copy. w_buf is a caller-provided scratch buffer of size
// >= K, reused across calls to avoid a heap allocation on every invocation.
static inline int sample_index_from_logw(const double* logw, int K,
                                          std::vector<double>& w_buf) {
    double m = logw[0];
    for (int k = 1; k < K; ++k) if (logw[k] > m) m = logw[k];
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { w_buf[k] = std::exp(logw[k] - m); tot += w_buf[k]; }
    double u = unif_rand() * tot;
    double cum = 0.0;
    for (int k = 0; k < K; ++k) {
        cum += w_buf[k];
        if (u <= cum) return k;
    }
    return K - 1;
}

// Draw out.size() indices (with replacement) in [0, K) via SYSTEMATIC
// resampling (Kitagawa, 1996): a single unif_rand() draw u0 ~ U(0, total/M)
// (M = out.size(), here K - 1 for the non-reference particles) determines M
// equally-spaced points u_j = u0 + total*j/M, located by ONE forward sweep
// through the cumulative weights (u_j is increasing in j, so the search
// pointer never needs to go backward). This replaces the previous
// multinomial scheme (M independent unif_rand() draws + M binary searches)
// with 1 draw + a single O(K) linear scan. Systematic resampling is
// unbiased and has lower variance than multinomial resampling (Doucet, de
// Freitas & Gordon, 2001), so this is a standard SMC choice, not merely a
// speed shortcut.
static inline void sample_indices_from_logw(const double* logw, int K,
                                             std::vector<double>& w_buf,
                                             std::vector<double>& cum_buf,
                                             std::vector<int>& out) {
    double m = logw[0];
    for (int k = 1; k < K; ++k) if (logw[k] > m) m = logw[k];
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { w_buf[k] = std::exp(logw[k] - m); tot += w_buf[k]; }
    cum_buf[0] = w_buf[0];
    for (int k = 1; k < K; ++k) cum_buf[k] = cum_buf[k - 1] + w_buf[k];
    double total = cum_buf[K - 1];
    int ndraws = (int) out.size();
    double u0 = unif_rand() * (total / ndraws);
    int k = 0;
    for (int j = 0; j < ndraws; ++j) {
        double u = u0 + total * ((double) j / ndraws);
        while (cum_buf[k] < u) ++k;
        out[j] = k;
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
// Main PG-AS sampler
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List pg_as_cpp(NumericVector y,
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
    // theta_1_k is flat, row-major (index t*K + k) to avoid the
    // pointer-chasing/cache-miss cost of std::vector<std::vector<double> >.
    std::vector<double> theta_1_k(Tt * (size_t)K, 0.0);
    // log_w_tilde only needs the previous and current time slice (see the
    // performance note above) -- kept as two K-sized buffers, swapped every
    // t, instead of a T x K matrix.
    std::vector<double> log_w_prev(K), log_w_curr(K);
    // A_hist (ancestor indices, for Backward Tracking) is flat, row-major
    // (index t*K + k), for the same cache-locality reason as theta_1_k.
    std::vector<int> A_hist(Tt * (size_t)K, 0);
    std::vector<double> log_as(K);
    std::vector<int> A(K - 1); // ancestors of the K - 1 non-reference particles
    std::vector<double> samp_w_buf(K), samp_cum_buf(K); // scratch for the resampling helpers

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

        // ---- Conditional SMC (Bootstrap PF + Ancestral Sampling) for theta_t1 ----
        // t = 0 (R: t=1)
        for (int k = 0; k < K; ++k) {
            theta_1_k[k] = R::norm_rand() * sd_W1 + (theta_01 + theta_02);
        }
        theta_1_k[K - 1] = theta1[0];
        for (int k = 0; k < K; ++k) log_w_curr[k] = log_p_yt(y[0], theta_1_k[k]);
        {
            double lse = logsumexp(log_w_curr);
            for (int k = 0; k < K; ++k) log_w_curr[k] -= lse;
        }
        std::swap(log_w_prev, log_w_curr); // log_w_prev now holds the t = 0 weights

        // t = 1 .. Tt-1 (R: t = 2..Tt)
        for (int t = 1; t < Tt; ++t) {
            const size_t base_prev = (size_t)(t - 1) * K;
            const size_t base_curr = (size_t)t * K;

            // resampling + propagation for the K - 1 non-reference particles
            // (bootstrap PF: resample directly from log_w_prev)
            sample_indices_from_logw(log_w_prev.data(), K, samp_w_buf, samp_cum_buf, A);
            for (int k = 0; k < K - 1; ++k) {
                double mean_k = theta_1_k[base_prev + A[k]] + theta2[t - 1];
                theta_1_k[base_curr + k] = R::norm_rand() * sd_W1 + mean_k;
            }

            // reference particle: fixed value, ancestor drawn via Ancestral Sampling
            theta_1_k[base_curr + (K - 1)] = theta1[t];
            for (int j = 0; j < K; ++j) {
                double mean_j = theta_1_k[base_prev + j] + theta2[t - 1];
                log_as[j] = log_w_prev[j] + log_dnorm(theta1[t], mean_j, sd_W1);
            }
            int a_ref = sample_index_from_logw(log_as.data(), K, samp_w_buf);

            for (int k = 0; k < K - 1; ++k) A_hist[base_curr + k] = A[k];
            A_hist[base_curr + (K - 1)] = a_ref;

            // updated weights (bootstrap: likelihood only)
            for (int k = 0; k < K; ++k)
                log_w_curr[k] = log_p_yt(y[t], theta_1_k[base_curr + k]);
            double lse = logsumexp(log_w_curr);
            for (int k = 0; k < K; ++k) log_w_curr[k] -= lse;

            std::swap(log_w_prev, log_w_curr); // log_w_prev now holds the weights at time t
        }

        // ESS at final time
        // (after the loop above, log_w_prev holds the weights at t = Tt - 1)
        {
            double s2 = 0.0;
            for (int k = 0; k < K; ++k) s2 += std::exp(2.0 * log_w_prev[k]);
            ess_smc[n] = 1.0 / s2;
        }

        // ---- Backward tracking for theta1 ----
        // PG-AS: the genealogy already encodes the backward path, so we
        // only need to look up the stored ancestor indices -- no Backward
        // Simulation (no density evaluations here).
        int b_idx = sample_index_from_logw(log_w_prev.data(), K, samp_w_buf);
        theta1[Tt - 1] = theta_1_k[(size_t)(Tt - 1) * K + b_idx];
        for (int t = Tt - 2; t >= 0; --t) {
            b_idx = A_hist[(size_t)(t + 1) * K + b_idx];
            theta1[t] = theta_1_k[(size_t)t * K + b_idx];
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
