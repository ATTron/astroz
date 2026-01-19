//! Error codes for C API

pub const Code = enum(i32) {
    ok = 0,
    // TLE errors
    badTleLength = -1,
    badChecksum = -2,
    // SGP4 errors
    deepSpaceNotSupported = -10,
    invalidEccentricity = -11,
    satelliteDecayed = -12,
    // Orbital mechanics errors
    valueError = -20,
    // General errors
    allocFailed = -100,
    nullPointer = -101,
    notInitialized = -102,
    unknown = -999,
};

/// Convert zig errors to C error code
pub fn fromError(err: anyerror) Code {
    return switch (err) {
        error.OutOfMemory => .allocFailed,
        else => .unknown,
    };
}
