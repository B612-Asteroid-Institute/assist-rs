//! Coordinate frame transformations.
//!
//! ASSIST works in **barycentric equatorial ICRF** (AU, AU/day).
//! THOR's propagator interface uses **heliocentric ecliptic J2000**.
//! These functions convert between the two.

/// J2000 obliquity of the ecliptic (IAU 2006, radians).
/// ε = 23°26′21.448″ = 23.4392911° = 0.40909280422 rad
pub const OBLIQUITY_J2000: f64 = 0.409_092_804_22;

/// cos(ε) — computed from OBLIQUITY_J2000 via Python's math.cos()
const COS_EPS: f64 = 0.917_482_062_070_108_2;
/// sin(ε)
const SIN_EPS: f64 = 0.397_777_155_929_776_9;

/// Rotate a 6-element state vector [x,y,z,vx,vy,vz] from equatorial to ecliptic.
///
/// Rotation about the x-axis by +ε (equatorial → ecliptic):
/// ```text
/// x_ecl =  x_eq
/// y_ecl =  cos(ε) * y_eq + sin(ε) * z_eq
/// z_ecl = -sin(ε) * y_eq + cos(ε) * z_eq
/// ```
pub fn equatorial_to_ecliptic(state: &[f64; 6]) -> [f64; 6] {
    let [x, y, z, vx, vy, vz] = *state;
    [
        x,
        COS_EPS * y + SIN_EPS * z,
        -SIN_EPS * y + COS_EPS * z,
        vx,
        COS_EPS * vy + SIN_EPS * vz,
        -SIN_EPS * vy + COS_EPS * vz,
    ]
}

/// Rotate a 6-element state vector from ecliptic to equatorial.
///
/// Inverse rotation about the x-axis by -ε:
/// ```text
/// x_eq =  x_ecl
/// y_eq =  cos(ε) * y_ecl - sin(ε) * z_ecl
/// z_eq =  sin(ε) * y_ecl + cos(ε) * z_ecl
/// ```
pub fn ecliptic_to_equatorial(state: &[f64; 6]) -> [f64; 6] {
    let [x, y, z, vx, vy, vz] = *state;
    [
        x,
        COS_EPS * y - SIN_EPS * z,
        SIN_EPS * y + COS_EPS * z,
        vx,
        COS_EPS * vy - SIN_EPS * vz,
        SIN_EPS * vy + COS_EPS * vz,
    ]
}

/// Rotate a 6×6 matrix from equatorial to ecliptic: R × M × R^T.
pub fn rotate_matrix_eq_to_ecl(m: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    // Build the 6×6 block-diagonal rotation matrix
    let r = rotation_matrix_eq_to_ecl();
    mat6x6_multiply(&mat6x6_multiply(&r, m), &transpose6x6(&r))
}

fn rotation_matrix_eq_to_ecl() -> [[f64; 6]; 6] {
    let mut r = [[0.0f64; 6]; 6];
    // Upper-left 3×3
    r[0][0] = 1.0;
    r[1][1] = COS_EPS;
    r[1][2] = SIN_EPS;
    r[2][1] = -SIN_EPS;
    r[2][2] = COS_EPS;
    // Lower-right 3×3
    r[3][3] = 1.0;
    r[4][4] = COS_EPS;
    r[4][5] = SIN_EPS;
    r[5][4] = -SIN_EPS;
    r[5][5] = COS_EPS;
    r
}

fn transpose6x6(m: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut t = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            t[i][j] = m[j][i];
        }
    }
    t
}

fn mat6x6_multiply(a: &[[f64; 6]; 6], b: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// Convert Cartesian difference vector to spherical (range, RA, Dec) and rates.
///
/// Input: topocentric Cartesian position `dx` and velocity `dv` in equatorial frame.
/// Output: `[rho, ra_rad, dec_rad, drho, dra, ddec]`
///
/// RA is in [0, 2π), Dec is in [-π/2, π/2].
pub fn cartesian_to_spherical(dx: [f64; 3], dv: [f64; 3]) -> [f64; 6] {
    let [x, y, z] = dx;
    let [vx, vy, vz] = dv;

    let rho_sq = x * x + y * y + z * z;
    let rho = rho_sq.sqrt();
    let xy_sq = x * x + y * y;
    let xy = xy_sq.sqrt();

    let mut ra = y.atan2(x);
    if ra < 0.0 {
        ra += std::f64::consts::TAU;
    }
    let dec = (z / rho).asin();

    // Rates (analytic derivatives)
    let drho = (x * vx + y * vy + z * vz) / rho;
    let dra = if xy_sq > 0.0 {
        (x * vy - y * vx) / xy_sq
    } else {
        0.0
    };
    let ddec = if xy > 0.0 {
        (vz * rho_sq - z * (x * vx + y * vy + z * vz)) / (rho_sq * xy)
    } else {
        0.0
    };

    [rho, ra, dec, drho, dra, ddec]
}

/// Compute the 6×6 Jacobian of the Cartesian-to-spherical transformation.
///
/// Maps perturbations in `(x, y, z, vx, vy, vz)` to perturbations in
/// `(ρ, α, δ, ρ̇, α̇, δ̇)`. Row index = output, column index = input.
///
/// The transformation is singular at the celestial poles (`ρ_xy = 0`), where
/// RA and its rate are not defined. Rows 1, 4 and the dδ/∂{x,y} column entries
/// of row 2 use `1/ρ_xy` factors and are left at zero when `ρ_xy = 0` rather
/// than producing NaNs.
pub fn cartesian_to_spherical_jacobian(dx: [f64; 3], dv: [f64; 3]) -> [[f64; 6]; 6] {
    let [x, y, z] = dx;
    let [vx, vy, vz] = dv;

    let rho_sq = x * x + y * y + z * z;
    let rho = rho_sq.sqrt();
    let rho3 = rho_sq * rho;
    let xy_sq = x * x + y * y;
    let xy = xy_sq.sqrt();

    let s = x * vx + y * vy + z * vz; // ρ · ρ̇
    let a = x * vy - y * vx; // ρ_xy² · α̇

    let mut jac = [[0.0f64; 6]; 6];

    // Row 0: ∂ρ/∂(x,y,z,vx,vy,vz)
    jac[0][0] = x / rho;
    jac[0][1] = y / rho;
    jac[0][2] = z / rho;

    // Row 1: ∂α/∂(x,y,z,vx,vy,vz) — α = atan2(y, x)
    if xy_sq > 0.0 {
        jac[1][0] = -y / xy_sq;
        jac[1][1] = x / xy_sq;
    }

    // Row 2: ∂δ/∂(x,y,z,vx,vy,vz) — δ = asin(z/ρ)
    if xy > 0.0 {
        jac[2][0] = -x * z / (rho_sq * xy);
        jac[2][1] = -y * z / (rho_sq * xy);
        jac[2][2] = xy / rho_sq;
    }

    // Row 3: ∂ρ̇/∂(x,y,z,vx,vy,vz) — ρ̇ = s/ρ
    jac[3][0] = (vx * rho_sq - x * s) / rho3;
    jac[3][1] = (vy * rho_sq - y * s) / rho3;
    jac[3][2] = (vz * rho_sq - z * s) / rho3;
    jac[3][3] = x / rho;
    jac[3][4] = y / rho;
    jac[3][5] = z / rho;

    // Row 4: ∂α̇/∂(x,y,z,vx,vy,vz) — α̇ = A/ρ_xy²
    if xy_sq > 0.0 {
        let xy4 = xy_sq * xy_sq;
        jac[4][0] = (vy * xy_sq - 2.0 * x * a) / xy4;
        jac[4][1] = (-vx * xy_sq - 2.0 * y * a) / xy4;
        jac[4][3] = -y / xy_sq;
        jac[4][4] = x / xy_sq;
    }

    // Row 5: ∂δ̇/∂(x,y,z,vx,vy,vz) — δ̇ = f/g, f = vz·ρ² - z·s, g = ρ²·ρ_xy
    if xy > 0.0 {
        let f = vz * rho_sq - z * s;
        let g = rho_sq * xy;
        // ∂f/∂x = 2x·vz - z·vx;     ∂g/∂x = x·(2·ρ_xy² + ρ²)/ρ_xy
        // ∂f/∂y = 2y·vz - z·vy;     ∂g/∂y = y·(2·ρ_xy² + ρ²)/ρ_xy
        // ∂f/∂z = vz·z - s;         ∂g/∂z = 2z·ρ_xy
        let dgdx = x * (2.0 * xy_sq + rho_sq) / xy;
        let dgdy = y * (2.0 * xy_sq + rho_sq) / xy;
        let dgdz = 2.0 * z * xy;
        let dfdx = 2.0 * x * vz - z * vx;
        let dfdy = 2.0 * y * vz - z * vy;
        let dfdz = vz * z - s;
        let g2 = g * g;
        jac[5][0] = (dfdx * g - f * dgdx) / g2;
        jac[5][1] = (dfdy * g - f * dgdy) / g2;
        jac[5][2] = (dfdz * g - f * dgdz) / g2;
        // Velocity partials: ∂g/∂v = 0, so ∂δ̇/∂v = (∂f/∂v) / g.
        jac[5][3] = -z * x / g;
        jac[5][4] = -z * y / g;
        jac[5][5] = xy_sq / g; // = ρ_xy / ρ²
    }

    jac
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_eq_ecl() {
        let state = [1.0, 0.5, -0.3, 0.01, 0.005, -0.003];
        let ecl = equatorial_to_ecliptic(&state);
        let back = ecliptic_to_equatorial(&ecl);
        for i in 0..6 {
            assert!((state[i] - back[i]).abs() < 1e-14, "element {i} mismatch");
        }
    }

    #[test]
    fn spherical_basic() {
        // Object along +x axis
        let dx = [2.0, 0.0, 0.0];
        let dv = [0.0, 0.0, 0.0];
        let sph = cartesian_to_spherical(dx, dv);
        assert!((sph[0] - 2.0).abs() < 1e-14, "rho");
        assert!(sph[1].abs() < 1e-14, "ra");
        assert!(sph[2].abs() < 1e-14, "dec");
    }

    /// Central-difference numerical partial of component `k` of the spherical
    /// output w.r.t. input `j` (0..=2 = x,y,z, 3..=5 = vx,vy,vz).
    fn numerical_partial(dx: [f64; 3], dv: [f64; 3], k: usize, j: usize, h: f64) -> f64 {
        let mut dx_p = dx;
        let mut dx_m = dx;
        let mut dv_p = dv;
        let mut dv_m = dv;
        if j < 3 {
            dx_p[j] += h;
            dx_m[j] -= h;
        } else {
            dv_p[j - 3] += h;
            dv_m[j - 3] -= h;
        }
        let plus = cartesian_to_spherical(dx_p, dv_p);
        let minus = cartesian_to_spherical(dx_m, dv_m);
        let mut d = plus[k] - minus[k];
        // Handle RA wrap: keep difference in (-π, π].
        if k == 1 {
            d = (d + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                - std::f64::consts::PI;
        }
        d / (2.0 * h)
    }

    fn max_row_err(
        analytic: &[[f64; 6]; 6],
        dx: [f64; 3],
        dv: [f64; 3],
        row: usize,
        h: f64,
    ) -> f64 {
        let mut worst = 0.0f64;
        for j in 0..6 {
            let num = numerical_partial(dx, dv, row, j, h);
            let ana = analytic[row][j];
            let err = (num - ana).abs();
            if err > worst {
                worst = err;
            }
        }
        worst
    }

    #[test]
    fn jacobian_matches_finite_differences_generic() {
        // A generic geometry: no coordinate axis, no special symmetry.
        let dx = [0.7, -0.4, 0.3];
        let dv = [0.02, 0.015, -0.008];
        let jac = cartesian_to_spherical_jacobian(dx, dv);
        let h = 1e-5;
        // Each row should agree with central-difference to roughly h² ≈ 1e-10
        // (times a geometry-dependent constant of order 1).
        for row in 0..6 {
            let err = max_row_err(&jac, dx, dv, row, h);
            assert!(err < 5e-8, "row {row}: max |analytic - FD| = {err:.3e}");
        }
    }

    #[test]
    fn jacobian_matches_finite_differences_near_equator() {
        // Small z — tests the δ̇ row where numerator approaches 0 but denom
        // is finite.
        let dx = [1.5, 0.8, 0.02];
        let dv = [0.01, -0.005, 0.003];
        let jac = cartesian_to_spherical_jacobian(dx, dv);
        let h = 1e-5;
        for row in 0..6 {
            let err = max_row_err(&jac, dx, dv, row, h);
            assert!(err < 5e-8, "row {row}: max |analytic - FD| = {err:.3e}");
        }
    }

    #[test]
    fn jacobian_pole_is_sanitized() {
        // Exactly on the +z axis: RA / RA rate / dδ/dx,y are undefined.
        // The function must return zeros in those slots rather than NaN.
        let dx = [0.0, 0.0, 1.0];
        let dv = [0.0, 0.0, 0.1];
        let jac = cartesian_to_spherical_jacobian(dx, dv);
        for row in 0..6 {
            for col in 0..6 {
                assert!(
                    jac[row][col].is_finite(),
                    "jac[{row}][{col}] = {} is not finite",
                    jac[row][col]
                );
            }
        }
    }
}
