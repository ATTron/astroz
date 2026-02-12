const std = @import("std");

pub const G = 6.6743e-11; // gravity constant
pub const c = 299792458; // sol (m/s)
pub const h = 6.62607015e-34; // planck constant
pub const au = 1.49597871e+11;
pub const j2k = 2000.0;
pub const j2000Jd: f64 = 2451545.0; // julian date of J2000 epoch
pub const julianDaysPerCentury: f64 = 36525.0;

pub const twoPi: f64 = 2.0 * std.math.pi;
pub const deg2rad: f64 = std.math.pi / 180.0;
pub const rad2deg: f64 = 180.0 / std.math.pi;

pub const secondsPerMinute: f64 = 60.0;
pub const minutesPerHour: f64 = 60.0;
pub const hoursPerDay: f64 = 24.0;
pub const secondsPerHour: f64 = 3600.0;
pub const secondsPerDay: f64 = 86400.0;
pub const minutesPerDay: f64 = 1440.0;

pub const arcminutesPerDegree: f64 = 60.0;
pub const arcsecondsPerDegree: f64 = 3600.0;
pub const degreesPerHour: f64 = 15.0; // for right ascension (360/24)

// Solar Radiation Pressure
pub const solarPressure: f64 = 4.56e-6; // N/m^2 at 1 AU
pub const auKm: f64 = 1.495978707e8; // 1 AU in km

pub const Sgp4GravityModel = struct {
    radiusEarthKm: f64, // km
    mu: f64, // km^3/s^2
    j2: f64, // J2
    j3: f64, // J3
    j4: f64, // J4
    xke: f64, // derived value
    tumin: f64, // derived value
    j3oj2: f64, // derived value
};

pub const wgs72: Sgp4GravityModel = .{
    .radiusEarthKm = 6378.135,
    .mu = 398600.8,
    .j2 = 0.001082616,
    .j3 = -0.00000253881,
    .j4 = -0.00000165597,
    .xke = 0.0743669161331734132,
    .tumin = 13.44683969695931,
    .j3oj2 = -0.00234506972242078,
};

pub const wgs84Flattening: f64 = 1.0 / 298.257223563;
pub const wgs84EccentricitySq: f64 = 2.0 * wgs84Flattening - wgs84Flattening * wgs84Flattening;

pub const wgs84: Sgp4GravityModel = .{
    .radiusEarthKm = 6378.137,
    .mu = 398600.5,
    .j2 = 0.00108262998905,
    .j3 = -0.00000253215306,
    .j4 = -0.00000161098761,
    .xke = 0.07436685316871385, // 60.0 / sqrt(re^3 / mu)
    .tumin = 13.446851082044981,
    .j3oj2 = -0.00233899967218727,
};

/// deep space threshold: orbital period > 225 minutes.
/// satellites with mean motion < 6.4 rev/day need SDP4 (not yet implemented).
pub const sgp4DeepSpaceThresholdMinutes: f64 = 225.0;

pub const CelestialBody = struct {
    name: []const u8,
    mass: f64, // kg
    mu: f64, // km^3/s^2
    mFractionSolarSystem: f64 = 0,
    mRadius: ?f64 = null, // km
    eqRadius: ?f64 = null, // km
    pRadius: ?f64 = null, // km
    semiMajorAxis: f64 = 0, // km
    perihelion: f64 = 0, // km
    aphelion: f64 = 0, // km
    period: f64 = 0, // days
    velocity: f64 = 0, // km/s
    eccentricity: f64 = 0,
    inclination: f64 = 0, // degrees
    oblateness: ?f64 = null,
    j2Perturbation: f64 = 0,
    j3Perturbation: f64 = 0,
    j4Perturbation: f64 = 0,
    seaLevelDensity: f64 = 0,
    scaleHeight: f64 = 0,
    rotationRate: f64 = 0, // rad/s
};

pub const sun: CelestialBody = .{
    .name = "sun",
    .mass = 1.8842e30,
    .mu = 1.32712e11,
    .mFractionSolarSystem = 9.98657e-1,
    .eqRadius = 695700,
    .j2Perturbation = 0.0000002,
    .seaLevelDensity = 1e-12,
    .scaleHeight = 50000.0,
    .rotationRate = 2.865e-6,
};
pub const mercury: CelestialBody = .{
    .name = "mercury",
    .mass = 3.30101e23,
    .mu = 2.20319e4,
    .mFractionSolarSystem = 1.65789e-7,
    .mRadius = 2439.4,
    .eqRadius = 2440.53,
    .pRadius = 2438.26,
    .semiMajorAxis = 5.79091e7,
    .perihelion = 4.60009e7,
    .aphelion = 6.98173e7,
    .period = 87.97,
    .velocity = 47.87,
    .eccentricity = 0.20564,
    .inclination = 7.01,
    .j2Perturbation = 0.00006,
    .seaLevelDensity = 1e-12,
    .scaleHeight = 200.0,
    .rotationRate = 1.24e-6,
};
pub const venus: CelestialBody = .{
    .name = "venus",
    .mass = 4.86732e24,
    .mu = 3.24859e5,
    .mFractionSolarSystem = 2.44455e-6,
    .mRadius = 6051.8,
    .eqRadius = 6051.8,
    .pRadius = 6051.8,
    .semiMajorAxis = 1.08209e8,
    .perihelion = 1.07477e8,
    .aphelion = 1.08940e8,
    .period = 224.70,
    .velocity = 35.02,
    .eccentricity = 0.00676,
    .inclination = 3.39,
    .j2Perturbation = 0.000027,
    .seaLevelDensity = 65.0,
    .scaleHeight = 15.9,
    .rotationRate = -2.99e-7,
};
pub const earth: CelestialBody = .{
    .name = "earth",
    .mass = 5.97219e24,
    .mu = wgs84.mu,
    .mFractionSolarSystem = 2.99946e-6,
    .mRadius = 6371.0084,
    .eqRadius = wgs84.radiusEarthKm,
    .pRadius = 6356.7519,
    .semiMajorAxis = 1.49598e8,
    .perihelion = 1.47100e8,
    .aphelion = 1.52096e8,
    .period = 365.26,
    .velocity = 29.78,
    .eccentricity = 0.01670,
    .oblateness = 0.003353,
    .j2Perturbation = wgs84.j2,
    .j3Perturbation = wgs84.j3,
    .j4Perturbation = wgs84.j4,
    .seaLevelDensity = 1.225,
    .scaleHeight = 7.249,
    .rotationRate = 7.2921150e-5,
};
pub const moon: CelestialBody = .{
    .name = "moon",
    .mass = 7.34581e22,
    .mu = 4.90280e3,
    .mFractionSolarSystem = 3.68934e-8,
    .mRadius = 1737.4,
    .semiMajorAxis = 3.83398e5,
    .perihelion = 3.62106e5,
    .aphelion = 4.04689e5,
    .period = 27.18,
    .velocity = 1.03,
    .eccentricity = 0.05555,
    .inclination = 23.71,
    .oblateness = 0.0012,
    .j2Perturbation = 0.0002027,
    .seaLevelDensity = 5e-13,
    .scaleHeight = 100.0,
    .rotationRate = 2.6617e-6,
};
pub const mars: CelestialBody = .{
    .name = "mars",
    .mass = 6.41693e23,
    .mu = 4.28284e4,
    .mFractionSolarSystem = 3.22282e-7,
    .mRadius = 3389.50,
    .eqRadius = 3396.19,
    .pRadius = 3376.20,
    .semiMajorAxis = 2.27939e8,
    .perihelion = 2.06645e8,
    .aphelion = 2.49233e8,
    .period = 686.97,
    .velocity = 24.13,
    .eccentricity = 0.09342,
    .inclination = 1.85,
    .oblateness = 0.00648,
    .j2Perturbation = 0.001964,
    .seaLevelDensity = 0.020,
    .scaleHeight = 11.1,
    .rotationRate = 7.088e-5,
};
pub const jupiter: CelestialBody = .{
    .name = "jupiter",
    .mass = 1.89852e27,
    .mu = 1.26713e8,
    .mFractionSolarSystem = 9.53510e-4,
    .mRadius = 69911,
    .eqRadius = 71492,
    .pRadius = 66854,
    .semiMajorAxis = 7.78321e8,
    .perihelion = 7.40603e8,
    .aphelion = 8.16038e8,
    .period = 4332.52,
    .velocity = 13.06,
    .eccentricity = 0.04846,
    .inclination = 1.30,
    .oblateness = 0.06487,
    .j2Perturbation = 0.014736,
    .seaLevelDensity = 0.16,
    .scaleHeight = 27.0,
    .rotationRate = 1.7585e-4,
};
pub const saturn: CelestialBody = .{
    .name = "saturn",
    .mass = 5.68460e26,
    .mu = 3.79406e7,
    .mFractionSolarSystem = 2.85502e-4,
    .mRadius = 58232,
    .eqRadius = 60268,
    .pRadius = 54364,
    .semiMajorAxis = 1.42910e9,
    .perihelion = 1.35096e9,
    .aphelion = 1.50724e9,
    .period = 10783.05,
    .velocity = 9.64,
    .eccentricity = 0.05468,
    .inclination = 2.49,
    .oblateness = 0.09796,
    .j2Perturbation = 0.016298,
    .seaLevelDensity = 0.19,
    .scaleHeight = 59.5,
    .rotationRate = 1.6378e-4,
};
pub const uranus: CelestialBody = .{
    .name = "uranus",
    .mass = 8.68192e25,
    .mu = 5.79456e6,
    .mFractionSolarSystem = 4.36039e-5,
    .mRadius = 25362,
    .eqRadius = 25559,
    .pRadius = 24973,
    .semiMajorAxis = 2.87479e9,
    .perihelion = 2.73854e9,
    .aphelion = 3.01104e9,
    .period = 30768.84,
    .velocity = 6.79,
    .eccentricity = 0.04739,
    .inclination = 0.77,
    .oblateness = 0.02293,
    .j2Perturbation = 0.003343,
    .seaLevelDensity = 0.42,
    .scaleHeight = 27.7,
    .rotationRate = -1.012e-4,
};
pub const neptune: CelestialBody = .{
    .name = "neptune",
    .mass = 1.02431e26,
    .mu = 6.83653e6,
    .mFractionSolarSystem = 5.14447e-5,
    .mRadius = 24622,
    .eqRadius = 24764,
    .pRadius = 24341,
    .semiMajorAxis = 4.50489e9,
    .perihelion = 4.46384e9,
    .aphelion = 4.54594e9,
    .period = 60357.05,
    .velocity = 5.43,
    .eccentricity = 0.00911,
    .inclination = 1.77,
    .oblateness = 0.01708,
    .j2Perturbation = 0.003411,
    .seaLevelDensity = 0.45,
    .scaleHeight = 19.7,
    .rotationRate = 1.083e-4,
};
pub const pluto: CelestialBody = .{
    .name = "pluto",
    .mass = 1.46158e22,
    .mu = 9.75500e2,
    .mFractionSolarSystem = 7.34061e-9,
    .mRadius = 1188.3,
    .semiMajorAxis = 5.91540e9,
    .perihelion = 4.44212e9,
    .aphelion = 7.38868e9,
    .period = 90821.51,
    .velocity = 4.74,
    .eccentricity = 0.24906,
    .inclination = 17.14,
    .oblateness = 0.002,
    .j2Perturbation = 0.00039,
    .seaLevelDensity = 1e-6,
    .scaleHeight = 50.0,
    .rotationRate = -1.139e-5,
};

pub const allBodies = [_]CelestialBody{
    sun,
    mercury,
    venus,
    earth,
    moon,
    mars,
    jupiter,
    saturn,
    uranus,
    neptune,
    pluto,
};
