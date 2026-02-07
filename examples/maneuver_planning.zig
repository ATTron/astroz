//! Force Model Propagation Example
//!
//! Demonstrates the composable force model and adaptive integrator APIs:
//! - TwoBody and J2 force models with ForceModel.wrap()
//! - Composite force model combining multiple perturbations
//! - DormandPrince 8(7) adaptive integrator

const std = @import("std");
const astroz = @import("astroz");
const constants = astroz.constants;
const propagators = astroz.propagators;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const mu = constants.earth.mu;
    const rEarth = constants.wgs84.radiusEarthKm;
    const j2 = constants.wgs84.j2;

    // Initial state: circular orbit at 400 km, 51.6 deg inclination
    const a = rEarth + 400;
    const vCirc = @sqrt(mu / a);
    const inc = 51.6 * constants.deg2rad;
    const state0 = [6]f64{
        a, 0,                 0,
        0, vCirc * @cos(inc), vCirc * @sin(inc),
    };
    const oneDay: f64 = 86400.0;

    // Two-body propagation
    var twobody = propagators.TwoBody.init(mu);
    const twoBForceModel = propagators.ForceModel.wrap(propagators.TwoBody, &twobody);
    var dp87 = propagators.DormandPrince87.initWithTolerance(1e-10, 1e-12);
    const stateTwoBody = try dp87.integrator().step(state0, 0, oneDay, twoBForceModel);

    // Composite: two-body + J2 perturbation
    var j2Model = propagators.J2.init(mu, j2, rEarth);
    const j2ForceModel = propagators.ForceModel.wrap(propagators.J2, &j2Model);
    const models = [_]propagators.ForceModel{ twoBForceModel, j2ForceModel };
    var composite = try propagators.Composite.init(allocator, &models);
    defer composite.deinit();

    const compositeForceModel = propagators.ForceModel.wrap(propagators.Composite, &composite);
    var dp87J2 = propagators.DormandPrince87.initWithTolerance(1e-10, 1e-12);
    const stateJ2 = try dp87J2.integrator().step(state0, 0, oneDay, compositeForceModel);

    // Compare
    const rTwoBody = @sqrt(stateTwoBody[0] * stateTwoBody[0] + stateTwoBody[1] * stateTwoBody[1] + stateTwoBody[2] * stateTwoBody[2]);
    const rJ2 = @sqrt(stateJ2[0] * stateJ2[0] + stateJ2[1] * stateJ2[1] + stateJ2[2] * stateJ2[2]);
    const dx = stateJ2[0] - stateTwoBody[0];
    const dy = stateJ2[1] - stateTwoBody[1];
    const dz = stateJ2[2] - stateTwoBody[2];
    const posDiff = @sqrt(dx * dx + dy * dy + dz * dz);

    std.debug.print("=== Two-body vs J2 Composite (1 day) ===\n", .{});
    std.debug.print("Two-body final radius: {d:.3} km\n", .{rTwoBody});
    std.debug.print("J2 final radius:       {d:.3} km\n", .{rJ2});
    std.debug.print("Position difference:   {d:.1} km\n", .{posDiff});
}
