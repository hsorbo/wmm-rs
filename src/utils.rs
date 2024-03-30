// Convert from year in decimal form (2020.5) to unix timestamp
pub fn decimal_year_to_unix_timestamp(decimal_year: f32) -> i64 {
    fn secs_in_year(year: i32) -> i64 {
        let is_leap_year = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
        if is_leap_year {
            366*24*60*60
        } else {
            365*24*60*60
        }
    }
    let mut unix = 0;
    for x in 1970 ..(decimal_year.floor() as i32) {
        unix += secs_in_year(x);
    }
    unix + (decimal_year.fract() * (secs_in_year(decimal_year.floor() as i32) as f32)) as i64
}
