const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;

    const rawTle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var testTle = try Tle.parse(rawTle, allocator);
    defer testTle.deinit();
    var sc = Spacecraft.init("dummy_sc", testTle, 300.000, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc.deinit();

    const dt = 120.0;
    const simulationTime = 3 * 24 * 60 * 60.0;
    const orbitalPeriod = 90 * 60.0;
    var t: f64 = 0;

    std.debug.print("{s:>6} {s:>34} {s:>30}\n", .{ "hr", "quaternion", "angular velocity (rad/s)" });
    while (t < simulationTime) : (t += dt) {
        // made-up periodic torques, applied directly to the angular velocity
        sc.angularVelocity[0] += 0.001 * @sin(2 * std.math.pi * t / (orbitalPeriod * 2)) * dt;
        sc.angularVelocity[1] += 0.0005 * @cos(2 * std.math.pi * t / (orbitalPeriod * 3)) * dt;
        sc.angularVelocity[2] += 0.0002 * @sin(2 * std.math.pi * t / orbitalPeriod) * dt;

        sc.updateAttitude();
        sc.propagateAttitude(dt);

        if (@mod(t, 6 * 3600.0) == 0) {
            const q = sc.quaternion;
            const w = sc.angularVelocity;
            std.debug.print("{d:6.0} [{d:7.4} {d:7.4} {d:7.4} {d:7.4}] [{d:8.4} {d:8.4} {d:8.4}]\n", .{ t / 3600.0, q[0], q[1], q[2], q[3], w[0], w[1], w[2] });
        }
    }
}
