"""Exception classes for astroz errors."""


class AstrozError(Exception):
    """Base exception for astroz errors."""

    pass


class TleParseError(AstrozError):
    """TLE parsing failed."""

    pass


class Sgp4Error(AstrozError):
    """SGP4 propagation error."""

    pass


class OrbitalError(AstrozError):
    """Orbital mechanics calculation error."""

    pass


_ERROR_MAP = {
    -1: (TleParseError, "Bad TLE length"),
    -2: (TleParseError, "Bad checksum"),
    -10: (Sgp4Error, "Deep space not supported"),
    -11: (Sgp4Error, "Invalid eccentricity"),
    -12: (Sgp4Error, "Satellite decayed"),
    -20: (OrbitalError, "Invalid orbital parameters"),
    -100: (AstrozError, "Allocation failed"),
    -101: (AstrozError, "Null pointer"),
    -102: (AstrozError, "Not initialized"),
}


def check(code: int) -> None:
    """Raise exception if error code is non-zero."""
    if code == 0:
        return
    err = _ERROR_MAP.get(code)
    if err:
        raise err[0](err[1])
    raise AstrozError(f"Unknown error: {code}")
