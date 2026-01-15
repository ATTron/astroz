//! Modular orbit propagation system

// core types
pub const Propagator = @import("Propagator.zig").Propagator;
pub const StateTime = @import("Propagator.zig").StateTime;

// integrators
pub const Integrator = @import("Integrator.zig").Integrator;
pub const Rk4 = @import("Integrator.zig").Rk4;
pub const Sgp4Integrator = @import("Integrator.zig").Sgp4Integrator;

// force models
pub const ForceModel = @import("ForceModel.zig").ForceModel;
pub const TwoBody = @import("ForceModel.zig").TwoBody;
pub const J2 = @import("ForceModel.zig").J2;
pub const Drag = @import("ForceModel.zig").Drag;
pub const Composite = @import("ForceModel.zig").Composite;

test {
    const std = @import("std");
    std.testing.refAllDecls(@This());
}
