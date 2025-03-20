const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const rawTle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var testTle = try Tle.parse(rawTle, allocator);
    defer testTle.deinit();
    var sc = Spacecraft.init("dummy_sc", testTle, 300.000, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc.deinit();

    sc.angularVelocity = .{ 0.0, 0.0, 0.0 };

    const dt = 120.0; // 2 mins time step
    const simulationTime = 3 * 24 * 60 * 60.0; // 3 days in seconds
    const orbitalPeriod = 90 * 60.0; // 90 minutes orbital period
    var t: f64 = 0;

    while (t < simulationTime) : (t += dt) {
        // Simulate a dramatic torque effect
        const torqueX = 0.001 * @sin(2 * std.math.pi * t / (orbitalPeriod * 2));
        const torqueY = 0.0005 * @cos(2 * std.math.pi * t / (orbitalPeriod * 3));
        const torqueZ = 0.0002 * @sin(2 * std.math.pi * t / orbitalPeriod);

        // Update angular velocity based on torque (simplified)
        sc.angularVelocity[0] += torqueX * dt;
        sc.angularVelocity[1] += torqueY * dt;
        sc.angularVelocity[2] += torqueZ * dt;

        // Update attitude
        sc.updateAttitude();
        sc.propagateAttitude(dt);

        // Simulate simple circular orbit
        const orbitRadius = 7000.0;
        const x = orbitRadius * @cos(2 * std.math.pi * t / orbitalPeriod);
        const y = orbitRadius * @sin(2 * std.math.pi * t / orbitalPeriod);
        const z = 0.0;

        std.log.debug("Showing orbiting info: {d},{d},{d},{d},{d},{d},{d},{d},{d},{d},{d}\n", .{
            t,
            sc.quaternion[0],
            sc.quaternion[1],
            sc.quaternion[2],
            sc.quaternion[3],
            sc.angularVelocity[0],
            sc.angularVelocity[1],
            sc.angularVelocity[2],
            x,
            y,
            z,
        });
    }
}
