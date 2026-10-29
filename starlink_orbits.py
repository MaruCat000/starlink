#!/usr/bin/env python3
"""
Starlink 衛星の軌道データ(TLE)を取得し、地球の3Dモデル上に現在位置と軌道を描く。

2026-10-29 作成
  きっかけ: 「SpaceXの大量の衛星の軌道を地球モデルの上に表示できたらかっこいい」
  方針    : himawari_gif.py と同じく 1 ファイル・引数で挙動を変える CLI。
            軌道計算は自作せず sgp4 ライブラリ(SGP4 の標準実装)に任せる。
  やらないこと: ブラウザ上で回る地球儀(CesiumJS / three.js)は次の段階。
            ここでは matplotlib の静止画(PNG)までを一周させる。
"""
import argparse
import math
import os
import sys
import time
import urllib.request
from datetime import datetime, timezone

import numpy as np
from sgp4.api import Satrec, SatrecArray, jday

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TLE_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle"
CACHE = "starlink.tle"
CACHE_MAX_AGE_H = 6          # CelesTrak の更新間隔は約2時間。頻繁に取りに行かない
EARTH_R = 6371.0             # km(平均半径。描画用なので楕円体は使わない)


# ---------- データ取得 ----------

def load_tle(path=CACHE, max_age_h=CACHE_MAX_AGE_H, refresh=False):
    """TLE を読み込む。キャッシュが新しければネットワークに行かない。"""
    fresh = os.path.exists(path) and (time.time() - os.path.getmtime(path)) < max_age_h * 3600
    if refresh or not fresh:
        print("fetching:", TLE_URL)
        req = urllib.request.Request(TLE_URL, headers={"User-Agent": "starlink-orbits-sample/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        if not data.lstrip().startswith(b"STARLINK") and b"\n1 " not in data:
            raise SystemExit("TLE らしくない応答でした。CelesTrak 側の制限か仕様変更の可能性があります")
        with open(path, "wb") as f:
            f.write(data)
    else:
        age_min = (time.time() - os.path.getmtime(path)) / 60
        print(f"using cache: {path} ({age_min:.0f} min old)")

    with open(path, encoding="ascii", errors="replace") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]
    tles = []
    for i in range(0, len(lines) - 2, 3):
        name, l1, l2 = lines[i], lines[i + 1], lines[i + 2]
        if l1.startswith("1 ") and l2.startswith("2 "):
            tles.append((name.strip(), l1, l2))
    if not tles:
        raise SystemExit("TLE を1件も読めませんでした")
    return tles


# ---------- 軌道計算と座標変換 ----------

def gmst_rad(jd, fr):
    """グリニッジ平均恒星時(ラジアン)。TEME 座標を地球固定座標へ回すのに使う。"""
    t = ((jd - 2451545.0) + fr) / 36525.0
    sec = (67310.54841
           + (876600.0 * 3600 + 8640184.812866) * t
           + 0.093104 * t * t
           - 6.2e-6 * t * t * t)
    return math.radians((sec % 86400.0) / 240.0)


def teme_to_ecef(r, theta):
    """慣性系(TEME)の位置 [N,3] を、地球固定系(ECEF)へ z 軸まわりに -theta 回転。"""
    c, s = math.cos(theta), math.sin(theta)
    x = c * r[:, 0] + s * r[:, 1]
    y = -s * r[:, 0] + c * r[:, 1]
    return np.stack([x, y, r[:, 2]], axis=1)


def ecef_to_lla(r):
    """ECEF [N,3] → 緯度・経度(度)・高度(km)。描画用なので地心緯度で済ませる。"""
    x, y, z = r[:, 0], r[:, 1], r[:, 2]
    lon = np.degrees(np.arctan2(y, x))
    lat = np.degrees(np.arctan2(z, np.hypot(x, y)))
    alt = np.linalg.norm(r, axis=1) - EARTH_R
    return lat, lon, alt


def propagate(tles, when):
    """全衛星の位置を when(UTC) で計算。戻り値は (名前配列, ECEF位置[km], 恒星時)。"""
    sats = SatrecArray([Satrec.twoline2rv(l1, l2) for _, l1, l2 in tles])
    jd, fr = jday(when.year, when.month, when.day, when.hour, when.minute,
                  when.second + when.microsecond / 1e6)
    err, r, _v = sats.sgp4(np.array([jd]), np.array([fr]))
    err, r = err[:, 0], r[:, 0, :]
    ok = err == 0
    bad = int((~ok).sum())
    if bad:
        print(f"skipped {bad} sats (sgp4 error: decayed or bad elements)")
    theta = gmst_rad(jd, fr)
    names = np.array([n for n, _, _ in tles])[ok]
    return names, teme_to_ecef(r[ok], theta), theta


def orbit_track(l1, l2, when, theta, minutes=96, step_min=1):
    """1衛星の約1周分の軌道。恒星時は when の値で固定し、その瞬間の地球に対する軌道面を描く。"""
    sat = Satrec.twoline2rv(l1, l2)
    pts = []
    for m in range(0, minutes + 1, step_min):
        t = when.timestamp() + m * 60
        d = datetime.fromtimestamp(t, tz=timezone.utc)
        jd, fr = jday(d.year, d.month, d.day, d.hour, d.minute, d.second)
        e, r, _ = sat.sgp4(jd, fr)
        if e == 0:
            pts.append(r)
    if not pts:
        return None
    return teme_to_ecef(np.array(pts), theta)


# ---------- 描画 ----------

def draw_earth(ax):
    u = np.linspace(0, 2 * np.pi, 72)
    v = np.linspace(0, np.pi, 36)
    x = EARTH_R * np.outer(np.cos(u), np.sin(v))
    y = EARTH_R * np.outer(np.sin(u), np.sin(v))
    z = EARTH_R * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color="#1f4e79", alpha=0.35, linewidth=0, shade=True)
    # 緯度経度の格子(30度おき)
    for lat in range(-60, 61, 30):
        ph = np.radians(lat)
        ax.plot(EARTH_R * np.cos(ph) * np.cos(u), EARTH_R * np.cos(ph) * np.sin(u),
                np.full_like(u, EARTH_R * np.sin(ph)), color="#9fb8d0", lw=0.4, alpha=0.6)
    for lon in range(0, 180, 30):
        th = np.radians(lon)
        ax.plot(EARTH_R * np.sin(v) * np.cos(th), EARTH_R * np.sin(v) * np.sin(th),
                EARTH_R * np.cos(v), color="#9fb8d0", lw=0.4, alpha=0.6)
    # 本初子午線とグリニッジの目印(ECEF なので x 軸正方向が経度0度)
    ax.plot([0, EARTH_R * 1.15], [0, 0], [0, 0], color="#e0a030", lw=1)


def visible_mask(r, elev, azim):
    """視点から見て地球の裏側に隠れる点を落とす。matplotlib の 3D は前後関係を隠してくれないため。"""
    e, a = math.radians(elev), math.radians(azim)
    view = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    along = r @ view                                   # 視線方向の成分
    perp = np.linalg.norm(r - np.outer(along, view), axis=1)  # 視線からの距離
    return (along > 0) | (perp > EARTH_R)


def plot_3d(names, r, tles_for_tracks, when, theta, out, elev, azim):
    fig = plt.figure(figsize=(10, 10), facecolor="black")
    ax = fig.add_subplot(111, projection="3d", facecolor="black")
    draw_earth(ax)

    _lat, _lon, alt = ecef_to_lla(r)
    vis = visible_mask(r, elev, azim)
    sc = ax.scatter(r[vis, 0], r[vis, 1], r[vis, 2], c=alt[vis], cmap="plasma", s=1.0,
                    vmin=300, vmax=600, depthshade=False)

    for _name, l1, l2 in tles_for_tracks:
        tr = orbit_track(l1, l2, when, theta)
        if tr is not None:
            tr = tr.astype(float)
            tr[~visible_mask(tr, elev, azim)] = np.nan   # 裏側は線を切る
            ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], color="#7fe0ff", lw=0.9, alpha=0.9)

    lim = EARTH_R + 1500
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("altitude [km]", color="white")
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb.ax.get_yticklabels(), color="white")
    jst = when.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    ax.set_title(f"Starlink  {len(names):,} satellites   {jst}\n"
                 f"data: CelesTrak (GP/TLE), propagation: SGP4", color="white", fontsize=11)
    fig.savefig(out, dpi=150, facecolor="black", bbox_inches="tight")
    print("saved:", out)


def plot_map(names, r, when, out):
    lat, lon, alt = ecef_to_lla(r)
    fig, ax = plt.subplots(figsize=(14, 7), facecolor="black")
    ax.set_facecolor("black")
    sc = ax.scatter(lon, lat, c=alt, cmap="plasma", s=1.5, vmin=300, vmax=600)
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 30)); ax.set_yticks(range(-90, 91, 30))
    ax.tick_params(colors="white"); ax.grid(color="#334", lw=0.5)
    for sp in ax.spines.values():
        sp.set_color("#667")
    cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("altitude [km]", color="white")
    plt.setp(cb.ax.get_yticklabels(), color="white")
    ax.set_title(f"Starlink  {len(names):,} satellites   "
                 f"{when.strftime('%Y-%m-%d %H:%M UTC')}", color="white")
    fig.savefig(out, dpi=150, facecolor="black", bbox_inches="tight")
    print("saved:", out)


# ---------- エントリポイント ----------

def main():
    ap = argparse.ArgumentParser(description="Starlink 衛星の位置と軌道を地球モデル上に描く")
    ap.add_argument("--time", help="UTC 時刻 (例: 2026-10-29T12:00)。省略時は現在")
    ap.add_argument("--tracks", type=int, default=12, help="軌道線を描く衛星の数(均等に抜き出す)。0 で描かない")
    ap.add_argument("--map", action="store_true", help="3D の代わりに緯度経度の2D散布図を出す")
    ap.add_argument("--out", default=None, help="出力 PNG (既定: starlink_3d.png / starlink_map.png)")
    ap.add_argument("--elev", type=float, default=25, help="3D の視点仰角(度)")
    ap.add_argument("--azim", type=float, default=135, help="3D の視点方位角(度)")
    ap.add_argument("--refresh", action="store_true", help="キャッシュを無視して TLE を取り直す")
    ap.add_argument("--max-age", type=float, default=CACHE_MAX_AGE_H, help="キャッシュを使う上限(時間)")
    args = ap.parse_args()

    when = (datetime.fromisoformat(args.time).replace(tzinfo=timezone.utc)
            if args.time else datetime.now(timezone.utc))

    tles = load_tle(refresh=args.refresh, max_age_h=args.max_age)
    print(f"loaded {len(tles):,} TLEs")
    names, r, theta = propagate(tles, when)
    lat, lon, alt = ecef_to_lla(r)
    print(f"altitude: min {alt.min():.0f} km / median {np.median(alt):.0f} km / max {alt.max():.0f} km")

    if args.map:
        plot_map(names, r, when, args.out or "starlink_map.png")
    else:
        n = max(args.tracks, 0)
        picks = [tles[i] for i in np.linspace(0, len(tles) - 1, n, dtype=int)] if n else []
        plot_3d(names, r, picks, when, theta, args.out or "starlink_3d.png", args.elev, args.azim)


if __name__ == "__main__":
    sys.exit(main())
