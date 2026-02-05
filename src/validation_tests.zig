//! Validation tests against known good values
//!
//! These tests verify the accuracy of force models and integrators against:
//! - Analytical solutions (Kepler problem)
//! - Published perturbation formulas (Vallado)
//! - Conservation laws (energy, angular momentum)

const std = @import("std");
const testing = std.testing;

const propagators = @import("propagators/propagators.zig");
const constants = @import("constants.zig");
const calculations = @import("calculations.zig");
const Spice = @import("Spice.zig");
const Sgp4 = @import("Sgp4.zig");
const Tle = @import("Tle.zig");

const FM = propagators.ForceModel;

const muEarth = 398600.4418; // km^3/s^2 (JGM-3/EGM96, used by poliastro reference data)
const rEarth = constants.wgs84.radiusEarthKm;
const j2Earth = constants.wgs84.j2;
const j3Earth = constants.wgs84.j3;
const j4Earth = constants.wgs84.j4;

test "validate two-body against analytical kepler - circular orbit" {
    // Circular orbit: after one period, should return to start
    const r0 = 7000.0; // km
    const v0 = calculations.orbitalVelocity(muEarth, r0, null);
    const period = calculations.orbitalPeriod(muEarth, r0);

    var tb = propagators.TwoBody.init(muEarth);
    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-12, 1e-14);

    const state0 = [6]f64{ r0, 0, 0, 0, v0, 0 };
    const state1 = try dp87.integrator().step(state0, 0, period, FM.wrap(propagators.TwoBody, &tb));

    // After one period, should return to initial position
    const posError = calculations.mag(.{ state1[0] - state0[0], state1[1] - state0[1], state1[2] - state0[2] });

    // Position error should be < 1 meter after one orbit
    try testing.expect(posError < 0.001); // 1 meter in km
}

test "validate two-body against analytical kepler - elliptical orbit" {
    // Elliptical orbit (e=0.1): verify specific energy is conserved
    const a = 7500.0; // semi-major axis km
    const e = 0.1; // eccentricity
    const rp = a * (1 - e); // periapsis = 6750 km
    const vp = calculations.orbitalVelocity(muEarth, rp, a);

    var tb = propagators.TwoBody.init(muEarth);
    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-12, 1e-14);
    const tbFm = FM.wrap(propagators.TwoBody, &tb);

    const state0 = [6]f64{ rp, 0, 0, 0, vp, 0 };

    // Specific orbital energy: = v^2/2 - mu/r = -mu/(2a)
    const energyExpected = -muEarth / (2 * a);

    // Propagate for half period (to apoapsis)
    const period = calculations.orbitalPeriod(muEarth, a);
    const state1 = try dp87.integrator().step(state0, 0, period / 2.0, tbFm);

    // Check energy at apoapsis
    const r1 = calculations.posMag(state1);
    const v1 = calculations.velMag(state1);
    const energy1 = v1 * v1 / 2 - muEarth / r1;

    // Energy should match to high precision
    const energyError = @abs(energy1 - energyExpected) / @abs(energyExpected);
    try testing.expect(energyError < 1e-10);

    // At apoapsis, radius should be ra = a(1+e)
    const raExpected = a * (1 + e);
    const raError = @abs(r1 - raExpected) / raExpected;
    try testing.expect(raError < 1e-8);
}

test "validate J2 acceleration magnitude" {
    // At position (r, 0, 0), J2 acceleration should be:
    // a_J2 = (3/2) * J2 * mu * R_e^2 / r^4 * (1 - 5*(z/r)^2) for x-component
    // For z=0 (equatorial): a_J2_x = (3/2) * J2 * mu * R_e^2 / r^4

    var j2 = propagators.J2.init(muEarth, j2Earth, rEarth);
    const r = 7000.0;
    const state = [6]f64{ r, 0, 0, 0, 7.5, 0 };

    const accel = FM.wrap(propagators.J2, &j2).acceleration(state, 0);

    // Expected J2 acceleration at equator (z=0)
    const expectedMag = 1.5 * j2Earth * muEarth * rEarth * rEarth / (r * r * r * r);
    const accelMag = calculations.mag(accel);

    // Should match within 1%
    const relError = @abs(accelMag - expectedMag) / expectedMag;
    try testing.expect(relError < 0.01);
}

test "validate J2 causes nodal regression" {
    // Propagate ISS-like orbit with J2 and verify angular momentum vector precesses
    const a = 6778.0;
    const i = 51.6 * constants.deg2rad;
    const v = calculations.orbitalVelocity(muEarth, a, null);

    var tb = propagators.TwoBody.init(muEarth);
    var j2 = propagators.J2.init(muEarth, j2Earth, rEarth);
    const models = [_]FM{ FM.wrap(propagators.TwoBody, &tb), FM.wrap(propagators.J2, &j2) };
    var comp = try propagators.Composite.init(testing.allocator, &models);
    defer comp.deinit();

    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-10, 1e-12);

    // Initial state: ISS-like orbit, inclined
    const state0 = [6]f64{
        a * @cos(i), 0, a * @sin(i), // position in orbital plane
        0, v, 0, // velocity perpendicular
    };

    // Propagate for 1 day
    const state1 = try dp87.integrator().step(state0, 0, constants.secondsPerDay, FM.wrap(propagators.Composite, &comp));

    // The orbit should have precessed - angular momentum vector should have rotated
    const h0 = calculations.cross(calculations.posVec(state0), calculations.velVec(state0));
    const h1 = calculations.cross(calculations.posVec(state1), calculations.velVec(state1));

    // Angular momentum direction should have changed (RAAN drift)
    const h0Mag = calculations.mag(h0);
    const h1Mag = calculations.mag(h1);

    // Magnitude should be approximately conserved
    const hMagError = @abs(h1Mag - h0Mag) / h0Mag;
    try testing.expect(hMagError < 0.001);

    // Direction should have changed slightly
    const dot = calculations.dot(h0, h1) / (h0Mag * h1Mag);
    try testing.expect(dot < 1.0); // Not exactly parallel
    try testing.expect(dot > 0.99); // But close (allows for ~8 degrees of precession)
}

test "validate J3 and J4 are smaller than J2" {
    // At any position, |a_J3| << |a_J2| and |a_J4| << |a_J2|
    // Based on coefficient ratios: J3/J2 ~= 0.002, J4/J2 ~= 0.001

    var j2 = propagators.J2.init(muEarth, j2Earth, rEarth);
    var j3 = propagators.J3.init(muEarth, j3Earth, rEarth);
    var j4 = propagators.J4.init(muEarth, j4Earth, rEarth);

    // Test at polar position (z != 0) where all terms are active
    const state = [6]f64{ 5000, 3000, 4000, 0, 6, 2 };

    const magJ2 = calculations.mag(FM.wrap(propagators.J2, &j2).acceleration(state, 0));
    const magJ3 = calculations.mag(FM.wrap(propagators.J3, &j3).acceleration(state, 0));
    const magJ4 = calculations.mag(FM.wrap(propagators.J4, &j4).acceleration(state, 0));

    // J3 should be ~1000x smaller than J2
    try testing.expect(magJ3 < magJ2 * 0.01);
    try testing.expect(magJ3 > magJ2 * 0.0001);

    // J4 should be significantly smaller than J2
    try testing.expect(magJ4 < magJ2 * 0.01);
    try testing.expect(magJ4 > 0); // Just verify it's non-zero
}

test "validate drag behavior" {
    var drag = propagators.ImprovedDrag.init(rEarth, 2.2, 10.0, 500.0, 1500.0, 150.0);
    const dragFm = FM.wrap(propagators.ImprovedDrag, &drag);

    // ISS-like orbit, prograde velocity
    const accel = dragFm.acceleration([6]f64{ 6778, 0, 0, 0, 7.66, 0 }, 0);

    // Drag should oppose velocity (negative y-component)
    try testing.expect(accel[1] < 0);

    // Drag magnitude: ~1e-8 km/s^2 at 400 km
    const mag = calculations.mag(accel);
    try testing.expect(mag > 1e-10);
    try testing.expect(mag < 1e-6);

    // Drag at 300 km should be significantly higher than at 500 km
    const mag300 = calculations.mag(dragFm.acceleration([6]f64{ rEarth + 300, 0, 0, 0, 7.73, 0 }, 0));
    const mag500 = calculations.mag(dragFm.acceleration([6]f64{ rEarth + 500, 0, 0, 0, 7.61, 0 }, 0));
    try testing.expect(mag300 > mag500 * 10);
}

test "validate SRP magnitude and shadow" {
    var srp = propagators.SolarRadiationPressure.init(1.5, 10.0, 500.0, rEarth);
    const srpFm = FM.wrap(propagators.SolarRadiationPressure, &srp);

    // In sunlight at GEO: expected ~1e-10 km/s^2, pointing away from Sun (-x)
    const accel = srpFm.acceleration([6]f64{ 42000, 0, 0, 0, 3.07, 0 }, 0);
    const mag = calculations.mag(accel);
    try testing.expect(mag > 1e-11);
    try testing.expect(mag < 1e-9);
    try testing.expect(accel[0] < 0);

    // In shadow (behind Earth): should be zero
    const accelShadow = srpFm.acceleration([6]f64{ -7000, 0, 0, 0, 7.5, 0 }, 0);
    try testing.expectEqual(@as(f64, 0), accelShadow[0]);
    try testing.expectEqual(@as(f64, 0), accelShadow[1]);
    try testing.expectEqual(@as(f64, 0), accelShadow[2]);

    // Just outside shadow cylinder: should be non-zero
    try testing.expect(calculations.mag(srpFm.acceleration([6]f64{ -7000, rEarth + 100, 0, 0, 7.5, 0 }, 0)) > 0);
}

test "validate DP87 more accurate than RK4 for long propagation" {
    var tb = propagators.TwoBody.init(muEarth);
    const tbFm = FM.wrap(propagators.TwoBody, &tb);

    const r0 = 7000.0;
    const v0 = calculations.orbitalVelocity(muEarth, r0, null);
    const state0 = [6]f64{ r0, 0, 0, 0, v0, 0 };
    const period = calculations.orbitalPeriod(muEarth, r0);

    // Propagate 10 orbits with RK4 (10 second steps)
    var rk4 = propagators.Rk4{};
    var rk4State = state0;
    const rk4Dt = 10.0;
    var t: f64 = 0;
    const totalTime = 10 * period;
    while (t < totalTime) {
        const step = @min(rk4Dt, totalTime - t);
        rk4State = try rk4.integrator().step(rk4State, t, step, tbFm);
        t += step;
    }

    // Propagate 10 orbits with DP87
    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-10, 1e-12);
    const dp87State = try dp87.integrator().step(state0, 0, totalTime, tbFm);

    // Check position error (should return to start after 10 orbits)
    const rk4Err = calculations.mag(.{
        rk4State[0] - state0[0],
        rk4State[1] - state0[1],
        rk4State[2] - state0[2],
    });
    const dp87Err = calculations.mag(.{
        dp87State[0] - state0[0],
        dp87State[1] - state0[1],
        dp87State[2] - state0[2],
    });

    // DP87 should be significantly more accurate
    try testing.expect(dp87Err < rk4Err);
    // DP87 error should be very small (< 10 meters after 10 orbits)
    try testing.expect(dp87Err < 0.01);
}

test "validate SPICE JD/ET conversion against known epochs" {
    // J2000 epoch: 2000-01-01T12:00:00 TDB = JD 2451545.0 = ET 0
    const etJ2000 = Spice.jdToEt(constants.j2000Jd);
    try testing.expectApproxEqAbs(@as(f64, 0.0), etJ2000, 1e-10);

    // One Julian century = 36525 days
    const etCentury = Spice.jdToEt(constants.j2000Jd + constants.julianDaysPerCentury);
    const expected = constants.julianDaysPerCentury * constants.secondsPerDay;
    try testing.expectApproxEqAbs(expected, etCentury, 1e-6);
}

test "validate SPICE enabled check" {
    if (!Spice.enabled) {
        try testing.expectError(Spice.SpiceError.NotEnabled, Spice.loadKernel("test.bsp"));
    }
}

test "validate composite model combines forces correctly" {
    var tb = propagators.TwoBody.init(muEarth);
    var j2 = propagators.J2.init(muEarth, j2Earth, rEarth);
    const tbFm = FM.wrap(propagators.TwoBody, &tb);
    const j2Fm = FM.wrap(propagators.J2, &j2);

    const state = [6]f64{ 7000, 0, 0, 0, 7.5, 0 };

    const aTb = tbFm.acceleration(state, 0);
    const aJ2 = j2Fm.acceleration(state, 0);

    const models = [_]FM{ tbFm, j2Fm };
    var comp = try propagators.Composite.init(testing.allocator, &models);
    defer comp.deinit();

    const aComp = FM.wrap(propagators.Composite, &comp).acceleration(state, 0);

    // Composite should be sum of individual accelerations
    try testing.expectApproxEqAbs(aTb[0] + aJ2[0], aComp[0], 1e-15);
    try testing.expectApproxEqAbs(aTb[1] + aJ2[1], aComp[1], 1e-15);
    try testing.expectApproxEqAbs(aTb[2] + aJ2[2], aComp[2], 1e-15);
}

test "validate ISS-like orbit propagation" {
    // ISS approximate orbital parameters:
    // Altitude: ~420 km, Inclination: 51.6, Period: ~92.7 min

    const alt = 420.0;
    const r = rEarth + alt;
    const v = calculations.orbitalVelocity(muEarth, r, null);
    const incl = 51.6 * constants.deg2rad;

    // Initial state in inclined orbit
    const state0 = [6]f64{
        r, 0,              0,
        0, v * @cos(incl), v * @sin(incl),
    };

    var tb = propagators.TwoBody.init(muEarth);
    var j2 = propagators.J2.init(muEarth, j2Earth, rEarth);
    const models = [_]FM{ FM.wrap(propagators.TwoBody, &tb), FM.wrap(propagators.J2, &j2) };
    var comp = try propagators.Composite.init(testing.allocator, &models);
    defer comp.deinit();

    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-10, 1e-12);

    // Propagate for 1 orbit (~5565 seconds)
    const period = calculations.orbitalPeriod(muEarth, r);
    const state1 = try dp87.integrator().step(state0, 0, period, FM.wrap(propagators.Composite, &comp));

    // Radius should be approximately conserved (circular orbit with J2)
    const r0 = calculations.posMag(state0);
    const r1 = calculations.posMag(state1);

    // Radius should stay within ~10 km (J2 causes small oscillations)
    try testing.expect(@abs(r1 - r0) < 10.0);

    // Velocity magnitude should be approximately conserved
    const v0 = calculations.velMag(state0);
    const v1 = calculations.velMag(state1);

    try testing.expect(@abs(v1 - v0) < 0.01); // Within 10 m/s
}

// Reference values generated by python-sgp4 v2.25
// TLE: ISS (NORAD 25544)
// Line 1: 1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9993
// Line 2: 2 25544  51.6400 200.0000 0001234  90.0000 270.0000 15.50000000000017

const sgp4Reference = [_]struct {
    minutes: f64,
    pos: [3]f64,
    vel: [3]f64,
}{
    .{ .minutes = 0.0, .pos = .{ -6388.6251039457, -2316.4105418668, -10.5105608841 }, .vel = .{ 1.6184721938, -4.4688535806, 6.0095711744 } },
    .{ .minutes = 30.0, .pos = .{ 4120.4071597960, -2524.2256992174, 4767.3396263188 }, .vel = .{ 5.7416339274, 4.3213445779, -2.6670882178 } },
    .{ .minutes = 60.0, .pos = .{ 2737.1301383237, 4555.4713955647, -4241.6250638171 }, .vel = .{ -6.7027938836, 0.6446517922, -3.6368442517 } },
    .{ .minutes = 90.0, .pos = .{ -6553.5121437682, -1479.8966428166, -1028.6875809764 }, .vel = .{ 0.1847982377, -4.8867524243, 5.8962546577 } },
    .{ .minutes = 120.0, .pos = .{ 3053.5363712212, -3227.0877061800, 5131.1902821670 }, .vel = .{ 6.5442323039, 3.6613485711, -1.5875959823 } },
    .{ .minutes = 180.0, .pos = .{ -6467.5318299250, -588.4786603082, -2008.8485417963 }, .vel = .{ -1.2580863427, -5.1084756785, 5.5652923317 } },
    .{ .minutes = 360.0, .pos = .{ -5564.4087742750, 1216.0553986884, -3713.3763010911 }, .vel = .{ -3.9513152824, -4.9410186329, 4.3086655820 } },
    .{ .minutes = 720.0, .pos = .{ -1479.3844792026, 3958.1883289274, -5328.0625952509 }, .vel = .{ -7.2451444567, -2.4477674475, 0.1933342585 } },
    .{ .minutes = 1440.0, .pos = .{ 6475.4254340451, 2031.6698883792, -355.1829280487 }, .vel = .{ -1.7172177336, 4.4491247647, -5.9962431939 } },
};

test "validate SGP4 against python-sgp4 reference" {
    // TLE matching the python-sgp4 reference above
    const testTle =
        \\1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9993
        \\2 25544  51.6400 200.0000 0001234  90.0000 270.0000 15.50000000000017
    ;

    var tle = try Tle.parse(testTle, testing.allocator);
    defer tle.deinit();

    const sgp4 = try Sgp4.init(tle, constants.wgs84);

    const posTol = 0.1; // 100 meters
    const velTol = 1e-4;

    for (sgp4Reference) |ref| {
        const result = try sgp4.propagate(ref.minutes);
        for (0..3) |i| {
            try testing.expectApproxEqAbs(ref.pos[i], result[0][i], posTol);
            try testing.expectApproxEqAbs(ref.vel[i], result[1][i], velTol);
        }
    }
}

test "validate two-body energy conservation over 100 orbits" {
    // Start with elliptical orbit
    const a = 8000.0;
    const e = 0.2;
    const rp = a * (1 - e);
    const vp = calculations.orbitalVelocity(muEarth, rp, a);

    const state0 = [6]f64{ rp, 0, 0, 0, vp, 0 };

    // Expected energy (should be conserved)
    const energyExpected = -muEarth / (2 * a);

    var tb = propagators.TwoBody.init(muEarth);
    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-12, 1e-14);
    const tbFm = FM.wrap(propagators.TwoBody, &tb);

    const period = calculations.orbitalPeriod(muEarth, a);
    const totalTime = 100 * period;

    // Propagate 100 orbits
    const stateFinal = try dp87.integrator().step(state0, 0, totalTime, tbFm);

    // Calculate final energy
    const rFinal = calculations.posMag(stateFinal);
    const vFinal = calculations.velMag(stateFinal);
    const energyFinal = vFinal * vFinal / 2 - muEarth / rFinal;

    // Energy should be conserved to 1e-10 relative error
    const relErr = @abs(energyFinal - energyExpected) / @abs(energyExpected);
    try testing.expect(relErr < 1e-10);
}
