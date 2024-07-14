const std = @import("std");

/// Gravity
pub const G = 6.6743e-11;
/// Speed of light
pub const c = 299792458;
/// Planck constant
pub const h = 6.62607015e-34;
/// Astronomical unit
pub const au = 1.49597871e+11;
/// Year 2000
pub const j2k = 2000.0;

pub const CelestialBody = struct {
    mass: f64, // kg
    mu: f64, // km^3/s^2
    m_fraction_solar_system: f64,
    m_radius: ?f64, //km
    eq_radius: ?f64, //km
    p_radius: ?f64, //km
    semi_major_axis: f64, // km
    perihelion: f64, // km
    aphelion: f64, // km
    period: f64, // days
    velocity: f64, // km/s
    eccentricity: f64,
    inclination: f64, // degrees
    oblateness: ?f64,
    j2_perturbation: f64,
    sea_level_density: f64,
    scale_height: f64,

    fn init(mass: f64, mu: f64, m_fraction_solar_system: f64, m_radius: ?f64, eq_radius: ?f64, p_radius: ?f64, semi_major_axis: ?f64, perihelion: ?f64, aphelion: ?f64, period: ?f64, velocity: ?f64, eccentricity: ?f64, inclination: ?f64, oblateness: ?f64, j2_perturbation: f64, sea_level_density: f64, scale_height: f64) CelestialBody {
        return .{
            .mass = mass,
            .mu = mu,
            .m_fraction_solar_system = m_fraction_solar_system,
            .m_radius = m_radius,
            .eq_radius = eq_radius,
            .p_radius = p_radius,
            .semi_major_axis = semi_major_axis orelse 0.0,
            .perihelion = perihelion orelse 0.0,
            .aphelion = aphelion orelse 0.0,
            .period = period orelse 0.0,
            .velocity = velocity orelse 0.0,
            .eccentricity = eccentricity orelse 0.0,
            .inclination = inclination orelse 0.0,
            .oblateness = oblateness,
            .j2_perturbation = j2_perturbation,
            .sea_level_density = sea_level_density,
            .scale_height = scale_height,
        };
    }
};

pub const sun = CelestialBody.init(1.8842e30, 1.32712e11, 9.98657e-1, null, 695700, null, null, null, null, null, null, null, null, null, 0.0000002, 1e-12, 50000.0);
pub const mercury = CelestialBody.init(3.30101e23, 2.20319e4, 1.65789e-7, 2439.4, 2440.53, 2438.26, 5.79091e7, 4.60009e7, 6.98173e7, 87.97, 47.87, 0.20564, 7.01, 0.000, 0.00006, 1e-12, 200.0);
pub const venus = CelestialBody.init(4.86732e24, 3.24859e5, 2.44455e-6, 6051.8, 6051.8, 6051.8, 1.08209e8, 1.07477e8, 1.08940e8, 224.70, 35.02, 0.00676, 3.39, 0.000, 0.000027, 65.0, 15.9);
pub const earth = CelestialBody.init(5.97219e24, 3.98600e5, 2.99946e-6, 6371.0084, 6378.1366, 6356.7519, 1.49598e8, 1.47100e8, 1.52096e8, 365.26, 29.78, 0.01670, 0.00, 0.003353, 0.00108263, 1.225, 7.249);
pub const moon = CelestialBody.init(7.34581e22, 4.90280e3, 3.68934e-8, 1737.4, null, null, 3.83398e5, 3.62106e5, 4.04689e5, 27.18, 1.03, 0.05555, 23.71, 0.0012, 0.0002027, 5e-13, 100.0);
pub const mars = CelestialBody.init(6.41693e23, 4.28284e4, 3.22282e-7, 3389.50, 3396.19, 3376.20, 2.27939e8, 2.06645e8, 2.49233e8, 686.97, 24.13, 0.09342, 1.85, 0.00648, 0.001964, 0.020, 11.1);
pub const jupiter = CelestialBody.init(1.89852e27, 1.26713e8, 9.53510e-4, 69911, 71492, 66854, 7.78321e8, 7.40603e8, 8.16038e8, 4332.52, 13.06, 0.04846, 1.30, 0.06487, 0.014736, 0.16, 27.0);
pub const saturn = CelestialBody.init(5.68460e26, 3.79406e7, 2.85502e-4, 58232, 60268, 54364, 1.42910e9, 1.35096e9, 1.50724e9, 10783.05, 9.64, 0.05468, 2.49, 0.09796, 0.016298, 0.19, 59.5);
pub const uranus = CelestialBody.init(8.68192e25, 5.79456e6, 4.36039e-5, 25362, 25559, 24973, 2.87479e9, 2.73854e9, 3.01104e9, 30768.84, 6.79, 0.04739, 0.77, 0.02293, 0.003343, 0.42, 27.7);
pub const neptune = CelestialBody.init(1.02431e26, 6.83653e6, 5.14447e-5, 24622, 24764, 24341, 4.50489e9, 4.46384e9, 4.54594e9, 60357.05, 5.43, 0.00911, 1.77, 0.01708, 0.003411, 0.45, 19.7);
pub const pluto = CelestialBody.init(1.46158e22, 9.75500e2, 7.34061e-9, 1188.3, null, null, 5.91540e9, 4.44212e9, 7.38868e9, 90821.51, 4.74, 0.24906, 17.14, 0.002, 0.00039, 1e-6, 50.0);

test "Test Celestial Bodies Made" {
    try std.testing.expectEqual(5.97219e24, earth.mass);
}
