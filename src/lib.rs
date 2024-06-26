pub mod constants;
pub mod utils;

const MAX_N: usize = 13;

pub struct Coefficients {
    pub(crate) g_coeff: [[f32; MAX_N]; MAX_N],
    pub(crate) h_coeff: [[f32; MAX_N]; MAX_N],
    pub(crate) delta_g: [[f32; MAX_N]; MAX_N],
    pub(crate) delta_h: [[f32; MAX_N]; MAX_N],
    pub(crate) base_time: i64,
}

fn compute_schmidt_quasi_norm_factors<const N: usize>() -> [[f32; N]; N] {
    let mut schmidt_quasi_norm: [[f32; N]; N] = [[0.0; N]; N];
    schmidt_quasi_norm[0][0] = 1.0;
    for n in 1..N {
        let mut v = schmidt_quasi_norm[n - 1][0] * (2 * n - 1) as f32 / n as f32;
        schmidt_quasi_norm[n][0] = v;
        for m in 1..=n {
            v *= (((n - m + 1) as f32 * if m == 1 { 2.0 } else { 1.0 }) / (n + m) as f32).sqrt();
            schmidt_quasi_norm[n][m] = v;
        }
    }
    schmidt_quasi_norm
}

fn compute_geocentric_coordinates(
    gd_lat_deg: f32,
    gd_lon_deg: f32,
    alt_meters: f32,
) -> (f32, f32, f32) {
    let altitude_km = alt_meters / 1000.0;
    let a2 = 4.0680636E7;
    let b2 = 4.04083E7;
    let gd_lat_rad = gd_lat_deg.to_radians();
    let clat = gd_lat_rad.cos();
    let slat = gd_lat_rad.sin();
    let tlat = slat / clat;
    let lat_rad = (a2 * clat * clat + b2 * slat * slat).sqrt();
    let gc_lat_rad = (tlat * (lat_rad * altitude_km + b2) / (lat_rad * altitude_km + a2)).atan();
    let gc_lon_rad = gd_lon_deg.to_radians();
    let rad_sq = altitude_km * altitude_km
        + 2.0 * altitude_km * (a2 * clat * clat + b2 * slat * slat).sqrt()
        + (a2 * a2 * clat * clat + b2 * b2 * slat * slat) / (a2 * clat * clat + b2 * slat * slat);
    let gc_radius_km = rad_sq.sqrt();
    (gc_lat_rad, gc_lon_rad, gc_radius_km)
}

struct LegendreTable {
    p: [[f32; MAX_N]; MAX_N],
    p_deriv: [[f32; MAX_N]; MAX_N],
}

impl LegendreTable {
    pub fn new(theta_rad: f32) -> Self {
        let cos = theta_rad.cos();
        let sin = theta_rad.sin();
        let mut t = LegendreTable {
            p: [[0.; MAX_N]; MAX_N],
            p_deriv: [[0.; MAX_N]; MAX_N],
        };
        t.p[0][0] = 1.0;
        t.p_deriv[0][0] = 0.0;
        for n in 1..=(MAX_N - 1) {
            for m in 0..=n {
                if n == m {
                    t.p[n][m] = sin * t.p[n - 1][m - 1];
                    t.p_deriv[n][m] = cos * t.p[n - 1][m - 1] + sin * t.p_deriv[n - 1][m - 1];
                } else if n == 1 || m == n - 1 {
                    t.p[n][m] = cos * t.p[n - 1][m];
                    t.p_deriv[n][m] = -sin * t.p[n - 1][m] + cos * t.p_deriv[n - 1][m];
                } else {
                    assert!(n > 1 && m < n - 1);
                    let nf = n as f32;
                    let mf = m as f32;
                    let k = ((nf - 1.) * (nf - 1.) - mf * mf) / ((2. * nf - 1.) * (2. * nf - 3.));
                    t.p[n][m] = cos * t.p[n - 1][m] - k * t.p[n - 2][m];
                    t.p_deriv[n][m] =
                        -sin * t.p[n - 1][m] + cos * t.p_deriv[n - 1][m] - k * t.p_deriv[n - 2][m];
                }
            }
        }
        t
    }
}

pub struct GeomagneticField {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    gc_lat_rad: f32,
    gc_lon_rad: f32,
    gc_radius_km: f32,
}

impl GeomagneticField {
    pub fn new(
        coeff: &Coefficients,
        gd_lat_deg: f32,
        gd_lon_deg: f32,
        alt_meters: f32,
        time_ms: i64,
    ) -> Self {
        //todo static

        let schmidt_quasi_norm_factors = compute_schmidt_quasi_norm_factors::<MAX_N>();

        let mut field = {
            let lat_cut = 89.99999f32.min(gd_lat_deg).max(-89.99999);
            let (gc_lat_rad, gc_lon_rad, gc_radius_km) =
                compute_geocentric_coordinates(lat_cut, gd_lon_deg, alt_meters);
            Self {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                gc_lat_rad,
                gc_lon_rad,
                gc_radius_km,
            }
        };

        let legendre = LegendreTable::new(std::f32::consts::FRAC_PI_2 - field.gc_lat_rad);
        let mut relative_radius_power = [0.0; MAX_N + 2];
        relative_radius_power[0] = 1.0;
        relative_radius_power[1] = 6371.2 / field.gc_radius_km;
        for n in 2..=MAX_N + 1 {
            relative_radius_power[n] = relative_radius_power[n - 1] * relative_radius_power[1];
        }

        let mut sin_lon = [0.; MAX_N];
        let mut cos_lon = [0.; MAX_N];
        sin_lon[0] = 0.;
        cos_lon[0] = 1.;
        sin_lon[1] = field.gc_lon_rad.sin();
        cos_lon[1] = field.gc_lon_rad.cos();
        for m in 2..=MAX_N - 1 {
            let x = m >> 1;
            sin_lon[m] = sin_lon[m - x] * cos_lon[x] + cos_lon[m - x] * sin_lon[x];
            cos_lon[m] = cos_lon[m - x] * cos_lon[x] - sin_lon[m - x] * sin_lon[x];
        }
        let inverse_cos_latitude = 1.0 / field.gc_lat_rad.cos();
        let time_delta = time_ms - coeff.base_time;
        let years_since_base = (time_delta as f32) / 3.1536001E10; //(1000.0 * 60.0 * 60.0 * 24.0 * 365.25);
        let mut x = 0f32;
        let mut z = 0f32;
        for n in 1..MAX_N {
            for m in 0..=n {
                let g = coeff.g_coeff[n][m] + coeff.delta_g[n][m] * years_since_base;
                let h = coeff.h_coeff[n][m] + coeff.delta_h[n][m] * years_since_base;
                x += relative_radius_power[n + 2]
                    * (g * cos_lon[m] + h * sin_lon[m])
                    * legendre.p_deriv[n][m]
                    * schmidt_quasi_norm_factors[n][m];

                field.y += relative_radius_power[n + 2]
                    * (m as f32)
                    * (g * sin_lon[m] - h * cos_lon[m])
                    * legendre.p[n][m]
                    * schmidt_quasi_norm_factors[n][m]
                    * inverse_cos_latitude;

                z -= (n as f32 + 1.)
                    * relative_radius_power[n + 2]
                    * (g * cos_lon[m] + h * sin_lon[m])
                    * legendre.p[n][m]
                    * schmidt_quasi_norm_factors[n][m];
            }
            let lat_diff_rad = gd_lat_deg.to_radians() - field.gc_lat_rad;
            field.x = x * lat_diff_rad.cos() + z * lat_diff_rad.sin();
            field.z = -x * lat_diff_rad.sin() + z * lat_diff_rad.cos();
        }

        field
    }
    pub fn get_declination(&self) -> f32 {
        (self.y).atan2(self.x).to_degrees()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]

    fn test_get_declination() {
        let x = GeomagneticField::new(
            &constants::WMM2015,
            20.200628,
            -87.500114,
            0.0,
            1609545600000,
        );
        let z = x.get_declination();
        assert_eq!(z, -1.7856733);
    }

    #[test]
    fn test_compute_schmidt_quasi_norm_factors() {
        let x = compute_schmidt_quasi_norm_factors::<6>();
        assert_eq!(
            x,
            [
                [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                [1.5, 1.7320508, 0.8660254, 0.0, 0.0, 0.0],
                [2.5, 3.0618622, 1.9364917, 0.7905695, 0.0, 0.0],
                [4.375, 5.533986, 3.913119, 2.09165, 0.73950994, 0.0],
                [7.875, 10.166581, 7.6852126, 4.7062125, 2.21853, 0.70156074],
            ]
        );
    }
}
