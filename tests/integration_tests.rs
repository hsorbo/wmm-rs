use float_eq::assert_float_eq;
use serde::Deserialize;
use std::path::Path;
use time::Duration;
use wmm_rs::utils;

#[derive(Deserialize, Debug)]
struct TestData {
    pub year: f32,
    pub height: i32,
    pub lat: i32,
    pub lon: i32,
    pub decl: f32,
    pub incl: f32,
    #[serde(rename = "H")]
    pub h: f32,
    #[serde(rename = "X")]
    pub x: f32,
    #[serde(rename = "Y")]
    pub y: f32,
    #[serde(rename = "Z")]
    pub z: f32,
    #[serde(rename = "F")]
    pub f: f32,
    pub dd_dt: f32,
    pub di_dt: f32,
    pub dh_dt: f32,
    pub dx_dt: f32,
    pub dy_dt: f32,
    pub dz_dt: f32,
    pub df_dt: f32,
}

#[test]
fn fairly_similar_to_wmm_crate() {
    let coeff = wmm_rs::constants::WMM2020;
    for lat in (-89..89).step_by(7) {
        for lon in (-179..179).step_by(13) {
            for year in (20200..20250).map(|x| (x as f32) / 10.) {
                let unix = utils::decimal_year_to_unix_timestamp(year);
                let actual =
                    wmm_rs::GeomagneticField::new(&coeff, lat as f32, lon as f32, 0.0, unix * 1000);
                let mut date =
                    time::Date::from_calendar_date(year as i32, time::Month::January, 1).unwrap();
                let days = year.fract() * 356.0;
                date = date.checked_add(Duration::days(days as i64)).unwrap();
                let other = wmm::declination(date, lat as f32, lon as f32).unwrap();
                //this is not great!
                assert_float_eq!(actual.get_declination(), other, abs <= 2.0);
            }
        }
    }
}

#[test]
fn wmm2020_testdata() {
    let coeff = wmm_rs::constants::WMM2020;
    let testdata = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/wmm2020_testdata.json");
    let contents = std::fs::read_to_string(testdata).unwrap();
    let wmm2020_test_data: Vec<TestData> = serde_json::from_str(&contents).unwrap();

    for data in wmm2020_test_data.iter() {
        let unix = utils::decimal_year_to_unix_timestamp(data.year);
        let result = wmm_rs::GeomagneticField::new(
            &coeff,
            data.lat as f32,
            data.lon as f32,
            data.height as f32,
            unix * 1000,
        );
        //this is also not great!

        assert_float_eq!(result.get_declination(), data.decl, abs <= 4.0);
    }
}
