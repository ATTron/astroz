//! Error codes for C API

pub const Code = enum(i32) {
    ok = 0,
    // TLE errors
    bad_tle_length = -1,
    bad_checksum = -2,
    // SGP4 errors
    deep_space_not_supported = -10,
    invalid_eccentricity = -11,
    satellite_decayed = -12,
    // Orbital mechanics errors
    value_error = -20,
    // General errors
    alloc_failed = -100,
    null_pointer = -101,
    not_initialized = -102,
    unknown = -999,
};

/// Convert zig errors to C error code
pub fn fromError(err: anyerror) Code {
    return switch (err) {
        error.OutOfMemory => .alloc_failed,
        else => .unknown,
    };
}
