const std = @import("std");

pub const G = 6.6743e-11; // gravity constant
pub const c = 299792458; // sol (m/s)
pub const h = 6.62607015e-34; // planck constant
pub const au = 1.49597871e+11;
pub const j2k = 2000.0;

pub const twoPi: f64 = 2.0 * std.math.pi;
pub const deg2rad: f64 = std.math.pi / 180.0;
pub const rad2deg: f64 = 180.0 / std.math.pi;
pub const minutes_per_day: f64 = 1440.0;

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

pub const wgs84: Sgp4GravityModel = .{
    .radiusEarthKm = 6378.137,
    .mu = 398600.5,
    .j2 = 0.00108262998905,
    .j3 = -0.00000253215306,
    .j4 = -0.00000161098761,
    .xke = 0.07436691613317342,
    .tumin = 13.44685108204498,
    .j3oj2 = -0.00233899967218727,
};

/// Deep space threshold: orbital period > 225 minutes.
/// Satellites with mean motion < 6.4 rev/day need SDP4 (not yet implemented).
pub const sgp4DeepSpaceThresholdMinutes: f64 = 225.0;

pub const CelestialBody = struct {
    name: []const u8,
    mass: f64, // kg
    mu: f64, // km^3/s^2
    mFractionSolarSystem: f64,
    mRadius: ?f64, // km
    eqRadius: ?f64, // km
    pRadius: ?f64, //km
    semiMajorAxis: f64, // km
    perihelion: f64, // km
    aphelion: f64, // km
    period: f64, // days
    velocity: f64, // km/s
    eccentricity: f64,
    inclination: f64, // degrees
    oblateness: ?f64,
    j2Perturbation: f64,
    seaLevelDensity: f64,
    scaleHeight: f64,
    rotationRate: f64, // rad/s

    fn init(
        name: []const u8,
        mass: f64,
        mu: f64,
        mFractionSolarSystem: f64,
        mRadius: ?f64,
        eqRadius: ?f64,
        pRadius: ?f64,
        semiMajorAxis: ?f64,
        perihelion: ?f64,
        aphelion: ?f64,
        period: ?f64,
        velocity: ?f64,
        eccentricity: ?f64,
        inclination: ?f64,
        oblateness: ?f64,
        j2Perturbation: f64,
        seaLevelDensity: f64,
        scaleHeight: f64,
        rotationRate: f64,
    ) CelestialBody {
        return .{
            .name = name,
            .mass = mass,
            .mu = mu,
            .mFractionSolarSystem = mFractionSolarSystem,
            .mRadius = mRadius,
            .eqRadius = eqRadius,
            .pRadius = pRadius,
            .semiMajorAxis = semiMajorAxis orelse 0.0,
            .perihelion = perihelion orelse 0.0,
            .aphelion = aphelion orelse 0.0,
            .period = period orelse 0.0,
            .velocity = velocity orelse 0.0,
            .eccentricity = eccentricity orelse 0.0,
            .inclination = inclination orelse 0.0,
            .oblateness = oblateness,
            .j2Perturbation = j2Perturbation,
            .seaLevelDensity = seaLevelDensity,
            .scaleHeight = scaleHeight,
            .rotationRate = rotationRate,
        };
    }
};

pub const sun = CelestialBody.init("sun", 1.8842e30, 1.32712e11, 9.98657e-1, null, 695700, null, null, null, null, null, null, null, null, null, 0.0000002, 1e-12, 50000.0, 2.865e-6);
pub const mercury = CelestialBody.init("mercury", 3.30101e23, 2.20319e4, 1.65789e-7, 2439.4, 2440.53, 2438.26, 5.79091e7, 4.60009e7, 6.98173e7, 87.97, 47.87, 0.20564, 7.01, 0.000, 0.00006, 1e-12, 200.0, 1.24e-6);
pub const venus = CelestialBody.init("venus", 4.86732e24, 3.24859e5, 2.44455e-6, 6051.8, 6051.8, 6051.8, 1.08209e8, 1.07477e8, 1.08940e8, 224.70, 35.02, 0.00676, 3.39, 0.000, 0.000027, 65.0, 15.9, -2.99e-7);
pub const earth = CelestialBody.init("earth", 5.97219e24, wgs84.mu, 2.99946e-6, 6371.0084, wgs84.radiusEarthKm, 6356.7519, 1.49598e8, 1.47100e8, 1.52096e8, 365.26, 29.78, 0.01670, 0.00, 0.003353, wgs84.j2, 1.225, 7.249, 7.2921150e-5);
pub const moon = CelestialBody.init("moon", 7.34581e22, 4.90280e3, 3.68934e-8, 1737.4, null, null, 3.83398e5, 3.62106e5, 4.04689e5, 27.18, 1.03, 0.05555, 23.71, 0.0012, 0.0002027, 5e-13, 100.0, 2.6617e-6);
pub const mars = CelestialBody.init("mars", 6.41693e23, 4.28284e4, 3.22282e-7, 3389.50, 3396.19, 3376.20, 2.27939e8, 2.06645e8, 2.49233e8, 686.97, 24.13, 0.09342, 1.85, 0.00648, 0.001964, 0.020, 11.1, 7.088e-5);
pub const jupiter = CelestialBody.init("jupiter", 1.89852e27, 1.26713e8, 9.53510e-4, 69911, 71492, 66854, 7.78321e8, 7.40603e8, 8.16038e8, 4332.52, 13.06, 0.04846, 1.30, 0.06487, 0.014736, 0.16, 27.0, 1.7585e-4);
pub const saturn = CelestialBody.init("saturn", 5.68460e26, 3.79406e7, 2.85502e-4, 58232, 60268, 54364, 1.42910e9, 1.35096e9, 1.50724e9, 10783.05, 9.64, 0.05468, 2.49, 0.09796, 0.016298, 0.19, 59.5, 1.6378e-4);
pub const uranus = CelestialBody.init("uranus", 8.68192e25, 5.79456e6, 4.36039e-5, 25362, 25559, 24973, 2.87479e9, 2.73854e9, 3.01104e9, 30768.84, 6.79, 0.04739, 0.77, 0.02293, 0.003343, 0.42, 27.7, -1.012e-4);
pub const neptune = CelestialBody.init("neptune", 1.02431e26, 6.83653e6, 5.14447e-5, 24622, 24764, 24341, 4.50489e9, 4.46384e9, 4.54594e9, 60357.05, 5.43, 0.00911, 1.77, 0.01708, 0.003411, 0.45, 19.7, 1.083e-4);
pub const pluto = CelestialBody.init("pluto", 1.46158e22, 9.75500e2, 7.34061e-9, 1188.3, null, null, 5.91540e9, 4.44212e9, 7.38868e9, 90821.51, 4.74, 0.24906, 17.14, 0.002, 0.00039, 1e-6, 50.0, -1.139e-5);

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
