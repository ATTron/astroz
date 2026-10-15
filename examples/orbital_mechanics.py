"""Orbital mechanics helpers: circular-orbit numbers, transfers and numerical propagation."""

import math

import numpy as np

from astroz import (
    EARTH_J2,
    EARTH_MU,
    EARTH_R_EQ,
    bi_elliptic_transfer,
    escape_velocity,
    hohmann_transfer,
    orbital_period,
    orbital_velocity,
    propagate_numerical,
)

def main():
    leo, geo = EARTH_R_EQ + 400.0, 42_164.0
    print(f"LEO (400 km): v = {orbital_velocity(EARTH_MU, leo):.3f} km/s, "
          f"period = {orbital_period(EARTH_MU, leo) / 60:.1f} min, "
          f"escape = {escape_velocity(EARTH_MU, leo):.3f} km/s")

    h = hohmann_transfer(EARTH_MU, leo, geo)
    print(f"LEO -> GEO Hohmann: {h['total_dv']:.3f} km/s over {h['transfer_time'] / 3600:.2f} h")

    # Bi-elliptic beats Hohmann once the radius ratio passes ~11.9
    far = 15 * leo
    hf = hohmann_transfer(EARTH_MU, leo, far)
    be = bi_elliptic_transfer(EARTH_MU, leo, far, 2 * far)
    print(f"LEO -> 15x LEO: Hohmann {hf['total_dv']:.3f} km/s, bi-elliptic {be['total_dv']:.3f} km/s")

    # One day in LEO with and without J2; J2 turns the orbit plane
    v = orbital_velocity(EARTH_MU, leo)
    inc = math.radians(51.6)
    state = (leo, 0.0, 0.0, 0.0, v * math.cos(inc), v * math.sin(inc))
    _, two_body = propagate_numerical(state, 0.0, 86400.0, 600.0, EARTH_MU)
    _, with_j2 = propagate_numerical(state, 0.0, 86400.0, 600.0, EARTH_MU, j2=EARTH_J2, r_eq=EARTH_R_EQ)
    gap = np.linalg.norm(np.asarray(two_body)[-1, :3] - np.asarray(with_j2)[-1, :3])
    print(f"After one day, J2 moves the satellite {gap:.0f} km from the two-body position")


if __name__ == "__main__":
    main()
