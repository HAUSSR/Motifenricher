//! Statistical routines: log-Gamma, regularized incomplete beta,
//! exact binomial upper-tail probabilities and Benjamini-Hochberg FDR.

/// Lanczos coefficients (g = 7, n = 9).
const LG_COEF: [f64; 9] = [
    0.999_999_999_999_809_93,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_13,
    -176.615_029_162_140_59,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_571_6e-6,
    1.505_632_735_149_311_6e-7,
];

/// Natural log of the Gamma function (Lanczos approximation).
pub fn lgamma(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x < 0.5 {
        let sin = (std::f64::consts::PI * x).sin();
        if sin == 0.0 {
            return f64::INFINITY;
        }
        return (std::f64::consts::PI / sin).ln() - lgamma(1.0 - x);
    }
    let z = x - 1.0;
    let mut ag = LG_COEF[0];
    for i in 1..9 {
        ag += LG_COEF[i] / (z + i as f64);
    }
    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + ag.ln()
}

/// Continued fraction for the regularized incomplete beta (Numerical Recipes).
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const MAXIT: usize = 300;
    const EPS: f64 = 3e-14;
    const FPMIN: f64 = 1e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c: f64 = 1.0;
    let mut d: f64 = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAXIT {
        let mf = m as f64;
        let m2 = 2.0 * mf;
        // even step
        let mut aa = mf * (b - mf) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        // odd step
        aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h
}

/// Regularized incomplete beta function I_x(a, b).
pub fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if a <= 0.0 || b <= 0.0 || x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let ln_bt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (1.0 - x).ln();
    let bt = ln_bt.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

/// P(X >= k) for X ~ Binomial(n, p), computed exactly via the
/// regularized incomplete beta function: I_p(k, n-k+1).
pub fn binom_sf(k: i64, n: i64, p: f64) -> f64 {
    if n <= 0 {
        return if k <= 0 { 1.0 } else { 0.0 };
    }
    if k <= 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    betainc(k as f64, (n - k + 1) as f64, p)
}

/// Benjamini-Hochberg step-up FDR-adjusted q-values.
pub fn bh_fdr(pvals: &[f64]) -> Vec<f64> {
    let m = pvals.len();
    if m == 0 {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..m).collect();
    idx.sort_by(|&a, &b| pvals[a].partial_cmp(&pvals[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut q = vec![1.0f64; m];
    let mut prev = 1.0f64;
    for rank in (1..=m).rev() {
        let i = idx[rank - 1];
        let val = (pvals[i] * m as f64 / rank as f64).min(1.0);
        prev = prev.min(val);
        q[i] = prev;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn lgamma_values() {
        assert!(close(lgamma(1.0), 0.0, 1e-10));
        assert!(close(lgamma(5.0), (24.0f64).ln(), 1e-9));
        assert!(close(lgamma(0.5), 0.5 * std::f64::consts::PI.ln(), 1e-9));
        assert!(close(lgamma(10.5), lgamma(9.5) + (9.5f64).ln(), 1e-9));
    }

    #[test]
    fn binom_sf_known() {
        // P(X >= 7), n = 10, p = 0.5  =>  176/1024 = 0.171875
        assert!(close(binom_sf(7, 10, 0.5), 176.0 / 1024.0, 1e-12));
        // P(X >= 1), n = 10, p = 0.5  =>  1023/1024
        assert!(close(binom_sf(1, 10, 0.5), 1023.0 / 1024.0, 1e-12));
        // P(X >= 60), n = 100, p = 0.5 => 0.02844397 (checked against R)
        assert!(close(binom_sf(60, 100, 0.5), 0.0284439668, 1e-8));
        assert_eq!(binom_sf(0, 5, 0.3), 1.0);
        assert_eq!(binom_sf(6, 5, 0.3), 0.0);
        assert_eq!(binom_sf(3, 5, 0.0), 0.0);
        assert_eq!(binom_sf(3, 5, 1.0), 1.0);
    }

    #[test]
    fn bh_fdr_monotone() {
        let ps = vec![0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205];
        let q = bh_fdr(&ps);
        // classic textbook example: q1 = 0.008, and q is monotone in p order
        assert!(close(q[0], 0.008, 1e-12));
        let mut sorted: Vec<(f64, f64)> = ps.iter().cloned().zip(q.iter().cloned()).collect();
        sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for w in sorted.windows(2) {
            assert!(w[0].1 <= w[1].1 + 1e-15);
        }
        assert!(q.iter().all(|&v| (0.0..=1.0).contains(&v)));
    }
}
