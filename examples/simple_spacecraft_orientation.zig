const std = @import("std");
const astroz = @import("astroz");
const Tle = astroz.Tle;
const constants = astroz.constants;
const Spacecraft = astroz.Spacecraft;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const raw_tle =
        \\1 55909U 23035B   24187.51050877  .00023579  00000+0  16099-2 0  9998
        \\2 55909  43.9978 311.8012 0011446 278.6226  81.3336 15.05761711 71371
    ;
    var test_tle = try Tle.parse(raw_tle, allocator);
    defer test_tle.deinit();
    var sc = Spacecraft.init("dummy_sc", test_tle, 300.000, Spacecraft.SatelliteSize.Cube, constants.earth, allocator);
    defer sc.deinit();

    sc.angular_velocity = .{ 0.0, 0.0, 0.0 };

    const dt = 120.0; // 2 mins time step
    const simulation_time = 3 * 24 * 60 * 60.0; // 3 days in seconds
    const orbital_period = 90 * 60.0; // 90 minutes orbital period
    var t: f64 = 0;

    while (t < simulation_time) : (t += dt) {
        // Simulate a dramatic torque effect
        const torque_x = 0.001 * @sin(2 * std.math.pi * t / (orbital_period * 2));
        const torque_y = 0.0005 * @cos(2 * std.math.pi * t / (orbital_period * 3));
        const torque_z = 0.0002 * @sin(2 * std.math.pi * t / orbital_period);

        // Update angular velocity based on torque (simplified)
        sc.angular_velocity[0] += torque_x * dt;
        sc.angular_velocity[1] += torque_y * dt;
        sc.angular_velocity[2] += torque_z * dt;

        // Update attitude
        sc.updateAttitude();
        sc.propagateAttitude(dt);

        // Simulate simple circular orbit
        const orbit_radius = 7000.0;
        const x = orbit_radius * @cos(2 * std.math.pi * t / orbital_period);
        const y = orbit_radius * @sin(2 * std.math.pi * t / orbital_period);
        const z = 0.0;

        std.log.debug("Showing orbiting info: {d},{d},{d},{d},{d},{d},{d},{d},{d},{d},{d}\n", .{
            t,
            sc.quaternion[0],
            sc.quaternion[1],
            sc.quaternion[2],
            sc.quaternion[3],
            sc.angular_velocity[0],
            sc.angular_velocity[1],
            sc.angular_velocity[2],
            x,
            y,
            z,
        });
    }
}
