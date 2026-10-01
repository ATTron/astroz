import { describe, expect, test } from "bun:test";
import * as sjs from "satellite.js";
import { AstrozError, Catalog, Satellite, parseTleText } from "../src/index";

const ISS = `ISS (ZARYA)
1 25544U 98067A   24127.82853009  .00015698  00000+0  27310-3 0  9995
2 25544  51.6393 160.4574 0003580 140.6673 205.7250 15.50957674452123`;

const GPS = `GPS BIIR-2  (PRN 13)
1 24876U 97035A   24186.50000000  .00000010  00000+0  00000+0 0  9998
2 24876  55.6000 100.0000 0050000  50.0000 310.0000  2.00560000200000`;

const MOLNIYA = `MOLNIYA 1-91
1 25485U 98054A   24186.50000000  .00000100  00000+0  10000-3 0  9991
2 25485  62.8000 200.0000 7000000 270.0000  10.0000  2.00600000200000`;

const OBSERVER = { lat: 40.7, lon: -74.0, alt: 0.01 };
const deg = 180 / Math.PI;

function reference(text: string, time: Date) {
  const [, l1, l2] = text.split("\n");
  const satrec = sjs.twoline2satrec(l1, l2);
  const pv = sjs.propagate(satrec, time) as { position: sjs.EciVec3<number>; velocity: sjs.EciVec3<number> };
  const gmst = sjs.gstime(time);
  const ecf = sjs.eciToEcf(pv.position, gmst);
  const look = sjs.ecfToLookAngles(
    { latitude: OBSERVER.lat / deg, longitude: OBSERVER.lon / deg, height: OBSERVER.alt },
    ecf,
  );
  const gd = sjs.eciToGeodetic(pv.position, gmst);
  return { pv, look, gd };
}

describe("matches satellite.js", () => {
  const times = [0, 90, 720, 1440].map((m) => new Date(Date.UTC(2024, 6, 5, 0, 0, 0) + m * 60000));

  for (const [label, text] of [["ISS (SGP4)", ISS], ["GPS (SDP4)", GPS], ["Molniya (SDP4)", MOLNIYA]]) {
    test(label, () => {
      const sat = Satellite.fromTle(text);
      for (const t of times) {
        const ref = reference(text, t);
        const { position, velocity } = sat.propagate(t);
        expect(position[0]).toBeCloseTo(ref.pv.position.x, 3);
        expect(position[1]).toBeCloseTo(ref.pv.position.y, 3);
        expect(position[2]).toBeCloseTo(ref.pv.position.z, 3);
        expect(velocity[0]).toBeCloseTo(ref.pv.velocity.x, 6);

        const look = sat.lookAngles(OBSERVER, t);
        expect(look.azimuth).toBeCloseTo(ref.look.azimuth * deg, 3);
        expect(look.elevation).toBeCloseTo(ref.look.elevation * deg, 3);
        expect(look.range).toBeCloseTo(ref.look.rangeSat, 2);

        const gd = sat.geodetic(t);
        expect(gd.lat).toBeCloseTo(ref.gd.latitude * deg, 3);
        expect(gd.alt).toBeCloseTo(ref.gd.height, 2);
      }
    });
  }
});

describe("Satellite", () => {
  test("metadata", () => {
    const iss = Satellite.fromTle(ISS);
    expect(iss.name).toBe("ISS (ZARYA)");
    expect(iss.noradId).toBe(25544);
    expect(iss.isDeepSpace).toBe(false);
    expect(iss.epoch.toISOString()).toBe("2024-05-06T19:53:05.000Z");
    expect(Satellite.fromTle(GPS).isDeepSpace).toBe(true);
  });

  test("accepts separate lines and unix milliseconds", () => {
    const [, l1, l2] = ISS.split("\n");
    const sat = Satellite.fromTle(l1, l2);
    const t = new Date("2024-05-07T00:00:00Z");
    expect(sat.propagate(t.getTime())).toEqual(sat.propagate(t));
    expect(sat.name).toBeUndefined();
  });

  test("range rate matches the change in range", () => {
    const sat = Satellite.fromTle(ISS);
    const t = Date.UTC(2024, 4, 7, 12);
    const a = sat.lookAngles(OBSERVER, t);
    const b = sat.lookAngles(OBSERVER, t + 1000);
    expect(b.range - a.range).toBeCloseTo((a.rangeRate + b.rangeRate) / 2, 3);
  });

  test("rejects bad input", () => {
    expect(() => Satellite.fromTle("not a tle")).toThrow(AstrozError);
    const [, l1, l2] = ISS.split("\n");
    expect(() => Satellite.fromTle(l1, l2.replace("25544", "25545"))).toThrow("invalid TLE");
  });

  test("decayed satellite throws on propagate", () => {
    const sat = Satellite.fromTle(ISS);
    expect(() => sat.propagate(new Date("2034-01-01T00:00:00Z"))).toThrow(AstrozError);
  });
});

describe("Catalog", () => {
  const text = [ISS, GPS, "garbage line", MOLNIYA].join("\n");

  test("parses mixed input and reports rejects", () => {
    const [, l1, l2] = ISS.split("\n");
    const cat = Catalog.fromTle(text + `\n${l1}\n${l2.slice(0, 40)}`);
    expect(cat.length).toBe(3);
    expect(cat.names).toEqual(["ISS (ZARYA)", "GPS BIIR-2  (PRN 13)", "MOLNIYA 1-91"]);
    expect(Array.from(cat.noradIds)).toEqual([25544, 24876, 25485]);
    expect(cat.noradIds.indexOf(24876)).toBe(1);
    expect(parseTleText(text).length).toBe(3);
  });

  test("batch (SIMD) results match single-satellite (scalar) results", () => {
    const cat = Catalog.fromTle(text);
    for (const t of [new Date("2024-07-05T06:00:00Z"), new Date("2024-07-01T00:00:00Z"), new Date("2024-07-09T00:00:00Z")]) {
      for (const frame of ["teme", "ecef", "geodetic"] as const) {
        const { positions, velocities, status } = cat.propagate(t, frame);
        expect(Array.from(status)).toEqual([0, 0, 0]);
        for (let i = 0; i < cat.length; i++) {
          const single = cat.satellite(i).propagate(t, frame);
          for (let k = 0; k < 3; k++) {
            expect(positions[i * 3 + k]).toBeCloseTo(single.position[k], 2);
            expect(velocities[i * 3 + k]).toBeCloseTo(single.velocity[k], 5);
          }
        }
      }
      const { values } = cat.lookAngles(OBSERVER, t);
      for (let i = 0; i < cat.length; i++) {
        const single = cat.satellite(i).lookAngles(OBSERVER, t);
        expect(values[i * 4]).toBeCloseTo(single.azimuth, 3);
        expect(values[i * 4 + 1]).toBeCloseTo(single.elevation, 3);
        expect(values[i * 4 + 2]).toBeCloseTo(single.range, 2);
      }
    }
  });

  test("batch path matches satellite.js across a mixed catalog", () => {
    const sats = [ISS, GPS, MOLNIYA];
    // 20 satellites: several full SIMD batches plus partial ones for SGP4 and SDP4
    const cat = Catalog.fromTle(Array.from({ length: 20 }, (_, i) => sats[i % 3]).join("\n"));
    const t = new Date("2024-07-05T12:34:56Z");
    const { positions } = cat.propagate(t);
    for (let i = 0; i < cat.length; i++) {
      const ref = reference(sats[i % 3], t).pv.position;
      expect(positions[i * 3]).toBeCloseTo(ref.x, 2);
      expect(positions[i * 3 + 1]).toBeCloseTo(ref.y, 2);
      expect(positions[i * 3 + 2]).toBeCloseTo(ref.z, 2);
    }
  });

  test("empty catalog", () => {
    const cat = Catalog.fromTle("");
    expect(cat.length).toBe(0);
    expect(cat.propagate(new Date()).positions.length).toBe(0);
    expect(cat.visible(OBSERVER, new Date())).toEqual([]);
  });

  test("visible filters by elevation", () => {
    const cat = Catalog.fromTle(text);
    const t = new Date("2024-07-05T06:00:00Z");
    const { values } = cat.lookAngles(OBSERVER, t);
    const expected = [0, 1, 2].filter((i) => values[i * 4 + 1] >= 10);
    expect(cat.visible(OBSERVER, t, 10)).toEqual(expected);
  });

  test("per-satellite errors do not affect the rest", () => {
    const cat = Catalog.fromTle(text);
    const { positions, status } = cat.propagate(new Date("2034-01-01T00:00:00Z"));
    expect(status[0]).not.toBe(0);
    expect(Number.isNaN(positions[0])).toBe(true);
    expect(status[1]).toBe(0);
  });

  test("free releases the catalog", () => {
    const cat = Catalog.fromTle(ISS);
    cat.free();
    cat.free();
    expect(cat.length).toBe(0);
    expect(() => cat.propagate(new Date())).toThrow("freed");
  });

  test("bad frame and index", () => {
    const cat = Catalog.fromTle(ISS);
    expect(() => cat.propagate(new Date(), "itrf" as never)).toThrow(AstrozError);
    expect(() => cat.satellite(1)).toThrow(RangeError);
  });
});
