const std = @import("std");

pub const First_Line = struct {
    line_number: u8,
    satellite_number: u32,
    classification: u8,
    intl_designator: u64,
    epoch_year: u16,
    epoch_day: f32,
    first_der_mean_motion: f32,
    secnd_der_mean_motion: f32,
    bstar_drag: f32,
    ephem_type: u8,
    elem_number: u32,
    checksum: u8,
};

pub const Second_Line = struct {
    line_number: u8,
    satellite_number: u32,
    inclination: f32,
    right_ascension: f32,
    eccentricity: f32,
    perigee: f32,
    m_anomaly: f32,
    m_motion: f32,
    rev_num: f32,
    checksum: u8,
};

pub const TLE = struct {
    first_line: First_Line,
    second_line: Second_Line,
};

// TODO: tle testing
// test "test tle values" {
//     const test_tle =
//         \\
//         \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
//         \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
//         \\
//     ;
//
// }
