import { wasmBase64 } from "./wasm.js";

/** Calendar time as a Date or Unix milliseconds (e.g. `Date.now()`) */
export type Time = Date | number;

/** Output frame for positions. Geodetic is [lat deg, lon deg, alt km]. */
export type Frame = "teme" | "ecef" | "geodetic";

export type GravityModel = "wgs72" | "wgs84";

/** Ground site: latitude and longitude in degrees, altitude above WGS84 in km */
export interface Observer {
  lat: number;
  lon: number;
  alt?: number;
}

export interface State {
  /** km, or [lat deg, lon deg, alt km] for geodetic */
  position: [number, number, number];
  /** km/s (ECEF velocity for geodetic) */
  velocity: [number, number, number];
}

export interface LookAngles {
  /** degrees clockwise from north, [0, 360) */
  azimuth: number;
  /** degrees above the horizon */
  elevation: number;
  /** slant range, km */
  range: number;
  /** km/s, positive when moving away */
  rangeRate: number;
}

export interface CatalogState {
  /** 3 values per satellite, see Frame */
  positions: Float64Array;
  velocities: Float64Array;
  /** per satellite: 0 = ok, otherwise a propagation error (outputs are NaN) */
  status: Uint8Array;
}

export interface CatalogLookAngles {
  /** 4 values per satellite: azimuth deg, elevation deg, range km, range rate km/s */
  values: Float64Array;
  status: Uint8Array;
}

export interface TleRecord {
  name?: string;
  line1: string;
  line2: string;
}

export interface CatalogOptions {
  gravity?: GravityModel;
}

export class AstrozError extends Error {}

interface Exports {
  memory: WebAssembly.Memory;
  alloc(len: number): number;
  free(ptr: number, len: number): void;
  catalog_new(grav: number): number;
  catalog_free(cat: number): void;
  catalog_len(cat: number): number;
  catalog_add(cat: number, l1: number, len1: number, l2: number, len2: number): number;
  catalog_norad_id(cat: number, i: number): number;
  catalog_epoch_jd(cat: number, i: number): number;
  catalog_deep_space(cat: number, i: number): number;
  catalog_propagate(cat: number, jd: number, frame: number, pos: number, vel: number, status: number): number;
  catalog_propagate_one(cat: number, i: number, jd: number, frame: number, pos: number, vel: number, status: number): void;
  catalog_look_angles(cat: number, jd: number, lat: number, lon: number, alt: number, out: number, status: number): number;
  catalog_look_angles_one(cat: number, i: number, jd: number, lat: number, lon: number, alt: number, out: number, status: number): void;
}

function decodeBase64(b64: string): Uint8Array {
  const bin = atob(b64);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return bytes;
}

const { instance } = await WebAssembly.instantiate(decodeBase64(wasmBase64).buffer as ArrayBuffer);
const wasm = instance.exports as unknown as Exports;

/** Status codes, in sync with bindings/javascript/zig/root.zig */
const STATUS_MESSAGES = ["ok", "invalid TLE", "invalid eccentricity", "satellite decayed", "out of memory"];
const FRAMES: Record<Frame, number> = { teme: 0, ecef: 1, geodetic: 2 };

function checkFrame(frame: Frame): number {
  const code = FRAMES[frame];
  if (code === undefined) throw new AstrozError(`unknown frame "${frame}"`);
  return code;
}

function statusMessage(status: number): string {
  return STATUS_MESSAGES[status] ?? `error ${status}`;
}

function checkStatus(status: number): void {
  if (status !== 0) throw new AstrozError(statusMessage(status));
}

function toJd(time: Time): number {
  const ms = typeof time === "number" ? time : time.getTime();
  return ms / 86400000 + 2440587.5;
}

function jdToDate(jd: number): Date {
  return new Date(Math.round((jd - 2440587.5) * 86400000));
}

function malloc(len: number): number {
  const ptr = wasm.alloc(len);
  if (ptr === 0) throw new AstrozError("out of memory");
  return ptr;
}

const encoder = new TextEncoder();

/** Split TLE text (2-line or 3-line format, mixed) into records */
export function parseTleText(text: string): TleRecord[] {
  const lines = text.split(/\r?\n/).map((l) => l.trimEnd()).filter((l) => l.trim().length > 0);
  const records: TleRecord[] = [];
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (line.startsWith("1 ") && lines[i + 1]?.startsWith("2 ")) {
      const prev = lines[i - 1];
      const name = prev && !prev.startsWith("1 ") && !prev.startsWith("2 ") ? prev.replace(/^0 /, "").trim() : undefined;
      records.push({ name, line1: line, line2: lines[i + 1] });
      i++;
    }
  }
  return records;
}

// Buffers in WASM memory are released when a catalog is garbage collected,
// unless free() was called first.
const registry = new FinalizationRegistry<() => void>((release) => release());

/**
 * A set of satellites propagated together. SGP4 or SDP4 is chosen per
 * satellite from its orbital period.
 */
export class Catalog {
  /** names from 3-line TLEs, undefined for 2-line records */
  readonly names: (string | undefined)[];
  readonly noradIds: Uint32Array;
  /** records that failed to parse or initialize */
  readonly rejected: { record: TleRecord; reason: string }[];

  #state: { handle: number; buf: number; bufLen: number };

  private constructor(records: TleRecord[], options: CatalogOptions) {
    const grav = options.gravity === "wgs84" ? 1 : 0;
    const handle = wasm.catalog_new(grav);
    if (handle === 0) throw new AstrozError("out of memory");
    const state = { handle, buf: 0, bufLen: 0 };
    this.#state = state;
    this.names = [];
    this.rejected = [];

    for (const record of records) {
      const l1 = encoder.encode(record.line1.trim());
      const l2 = encoder.encode(record.line2.trim());
      const ptr = malloc(l1.length + l2.length);
      new Uint8Array(wasm.memory.buffer, ptr, l1.length).set(l1);
      new Uint8Array(wasm.memory.buffer, ptr + l1.length, l2.length).set(l2);
      const status = wasm.catalog_add(handle, ptr, l1.length, ptr + l1.length, l2.length);
      wasm.free(ptr, l1.length + l2.length);
      if (status === 0) this.names.push(record.name);
      else this.rejected.push({ record, reason: statusMessage(status) });
    }

    this.noradIds = new Uint32Array(this.length);
    for (let i = 0; i < this.length; i++) this.noradIds[i] = wasm.catalog_norad_id(handle, i);

    registry.register(this, () => {
      if (state.buf) wasm.free(state.buf, state.bufLen);
      wasm.catalog_free(state.handle);
    }, this);
  }

  /** Parse TLE text in 2-line or 3-line format. Invalid records end up in `rejected`. */
  static fromTle(text: string, options: CatalogOptions = {}): Catalog {
    return new Catalog(parseTleText(text), options);
  }

  static fromRecords(records: TleRecord[], options: CatalogOptions = {}): Catalog {
    return new Catalog(records, options);
  }

  /**
   * Download a CelesTrak group, e.g. "stations", "starlink", "gps-ops", "active".
   * See https://celestrak.org/NORAD/elements/
   */
  static async fetch(group: string, options: CatalogOptions = {}): Promise<Catalog> {
    const url = `https://celestrak.org/NORAD/elements/gp.php?GROUP=${encodeURIComponent(group)}&FORMAT=tle`;
    const res = await fetch(url);
    if (!res.ok) throw new AstrozError(`CelesTrak request failed: ${res.status} ${res.statusText}`);
    return Catalog.fromTle(await res.text(), options);
  }

  get length(): number {
    return this.#state.handle === 0 ? 0 : wasm.catalog_len(this.#state.handle);
  }

  /** TLE epoch of satellite i */
  epoch(i: number): Date {
    this.#check(i);
    return jdToDate(wasm.catalog_epoch_jd(this.#state.handle, i));
  }

  /** True if satellite i uses SDP4 (period >= 225 minutes) */
  isDeepSpace(i: number): boolean {
    this.#check(i);
    return wasm.catalog_deep_space(this.#state.handle, i) !== 0;
  }

  /** Single-satellite view of entry i */
  satellite(i: number): Satellite {
    this.#check(i);
    return new Satellite(this, i);
  }

  /** Propagate every satellite to `time` */
  propagate(time: Time, frame: Frame = "teme"): CatalogState {
    const n = this.#checkOpen();
    if (n === 0) return { positions: new Float64Array(0), velocities: new Float64Array(0), status: new Uint8Array(0) };
    const f = checkFrame(frame);
    const buf = this.#scratch(n * 49);
    const pos = buf, vel = buf + n * 24, status = buf + n * 48;
    checkStatus(wasm.catalog_propagate(this.#state.handle, toJd(time), f, pos, vel, status));
    const mem = wasm.memory.buffer;
    return {
      positions: new Float64Array(mem, pos, n * 3).slice(),
      velocities: new Float64Array(mem, vel, n * 3).slice(),
      status: new Uint8Array(mem, status, n).slice(),
    };
  }

  /** Look angles from a ground observer to every satellite at `time` */
  lookAngles(observer: Observer, time: Time): CatalogLookAngles {
    const n = this.#checkOpen();
    if (n === 0) return { values: new Float64Array(0), status: new Uint8Array(0) };
    const buf = this.#scratch(n * 33);
    const out = buf, status = buf + n * 32;
    checkStatus(wasm.catalog_look_angles(this.#state.handle, toJd(time), observer.lat, observer.lon, observer.alt ?? 0, out, status));
    const mem = wasm.memory.buffer;
    return {
      values: new Float64Array(mem, out, n * 4).slice(),
      status: new Uint8Array(mem, status, n).slice(),
    };
  }

  /** @internal used by Satellite; runs the scalar propagator for one entry */
  _propagateOne(i: number, time: Time, frame: Frame): State {
    this.#check(i);
    const f = checkFrame(frame);
    const buf = this.#scratch(49);
    wasm.catalog_propagate_one(this.#state.handle, i, toJd(time), f, buf, buf + 24, buf + 48);
    const mem = wasm.memory.buffer;
    checkStatus(new Uint8Array(mem, buf + 48, 1)[0]);
    const p = new Float64Array(mem, buf, 6);
    return { position: [p[0], p[1], p[2]], velocity: [p[3], p[4], p[5]] };
  }

  /** @internal */
  _lookAnglesOne(i: number, observer: Observer, time: Time): LookAngles {
    this.#check(i);
    const buf = this.#scratch(33);
    wasm.catalog_look_angles_one(this.#state.handle, i, toJd(time), observer.lat, observer.lon, observer.alt ?? 0, buf, buf + 32);
    const mem = wasm.memory.buffer;
    checkStatus(new Uint8Array(mem, buf + 32, 1)[0]);
    const v = new Float64Array(mem, buf, 4);
    return { azimuth: v[0], elevation: v[1], range: v[2], rangeRate: v[3] };
  }

  /** Indices of satellites at or above `minElevation` degrees for the observer */
  visible(observer: Observer, time: Time, minElevation = 0): number[] {
    const { values, status } = this.lookAngles(observer, time);
    const out: number[] = [];
    for (let i = 0; i < status.length; i++) {
      if (status[i] === 0 && values[i * 4 + 1] >= minElevation) out.push(i);
    }
    return out;
  }

  /** Release WASM memory now instead of waiting for garbage collection */
  free(): void {
    const s = this.#state;
    if (s.handle === 0) return;
    registry.unregister(this);
    if (s.buf) wasm.free(s.buf, s.bufLen);
    wasm.catalog_free(s.handle);
    s.handle = s.buf = s.bufLen = 0;
  }

  /** Reusable output buffer, 8-byte aligned */
  #scratch(bytes: number): number {
    const s = this.#state;
    if (s.bufLen < bytes) {
      if (s.buf) wasm.free(s.buf, s.bufLen);
      s.buf = s.bufLen = 0;
      s.buf = malloc(bytes);
      s.bufLen = bytes;
    }
    return s.buf;
  }

  #check(i: number): void {
    if (this.#state.handle === 0) throw new AstrozError("catalog has been freed");
    if (!Number.isInteger(i) || i < 0 || i >= this.length) throw new RangeError(`satellite index ${i} out of range`);
  }

  #checkOpen(): number {
    if (this.#state.handle === 0) throw new AstrozError("catalog has been freed");
    return this.length;
  }
}

/** One satellite. Create with `Satellite.fromTle` or `catalog.satellite(i)`. */
export class Satellite {
  readonly catalog: Catalog;
  readonly index: number;

  /** @internal */
  constructor(catalog: Catalog, index: number) {
    this.catalog = catalog;
    this.index = index;
  }

  /** Parse a single TLE. Accepts the two lines separately or one 2/3-line string. */
  static fromTle(lineOrText: string, line2?: string, options: CatalogOptions = {}): Satellite {
    const records = line2 === undefined ? parseTleText(lineOrText) : [{ line1: lineOrText, line2 }];
    if (records.length === 0) throw new AstrozError("no TLE found in input");
    const catalog = Catalog.fromRecords(records.slice(0, 1), options);
    if (catalog.length === 0) throw new AstrozError(catalog.rejected[0].reason);
    return catalog.satellite(0);
  }

  get name(): string | undefined {
    return this.catalog.names[this.index];
  }

  get noradId(): number {
    return this.catalog.noradIds[this.index];
  }

  get epoch(): Date {
    return this.catalog.epoch(this.index);
  }

  get isDeepSpace(): boolean {
    return this.catalog.isDeepSpace(this.index);
  }

  /** Position and velocity at `time`. Throws if propagation fails. */
  propagate(time: Time, frame: Frame = "teme"): State {
    return this.catalog._propagateOne(this.index, time, frame);
  }

  /** Sub-satellite point and altitude */
  geodetic(time: Time): { lat: number; lon: number; alt: number } {
    const [lat, lon, alt] = this.propagate(time, "geodetic").position;
    return { lat, lon, alt };
  }

  /** Azimuth, elevation, range and range rate from a ground observer */
  lookAngles(observer: Observer, time: Time): LookAngles {
    return this.catalog._lookAnglesOne(this.index, observer, time);
  }
}
