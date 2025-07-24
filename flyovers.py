import os
from skyfield.api import load, EarthSatellite, wgs84
from datetime import datetime, timedelta, timezone, date
from zoneinfo import ZoneInfo
import requests

# ───── CONFIGURATION ──────────────────────────────────────────────────────────
OWM_API_KEY = os.environ.get("OWM_API_KEY")

CONFIG = {
    # UTC date to examine passes for
    "date_utc": date(2025, 7, 23),

    # elevation threshold = your imaging pointing accuracy (degrees)
    "min_elevation_deg": 0.5,

    # sites to evaluate: name, (lat, lon), and local timezone
    "sites": [
        {"name": "Vienna",          "lat": 48.2082, "lon": 16.3738,   "tz": "Europe/Vienna", "elevation_m": 170},
        {"name": "Munich",          "lat": 48.1351, "lon": 11.5820,   "tz": "Europe/Berlin", "elevation_m": 520},
        {"name": "Santa Barbara",   "lat": 34.4208, "lon": -119.6982, "tz": "America/Los_Angeles", "elevation_m": 15},
        {"name": "Milan",           "lat": 45.4642, "lon": 9.1900,    "tz": "Europe/Rome", "elevation_m": 120},
        {"name": "Helgoland",       "lat": 54.1833, "lon": 7.9167,    "tz": "Europe/Berlin", "elevation_m": 4},
    ],
    
    # satellite's TLE
    "tle": [
        "1 99999U 25007CV  25203.60326757  .00004178  00000-0  22378-3 0  9992",
        "2 99999  97.4471 317.6742 0001204 152.1272 151.6357 15.17801187 59108"
    ]
}
# ────────────────────────────────────────────────────────────────────────────────

def fetch_cloud_forecast(lat, lon):
    """
    Fetch hourly cloud-cover (%) for the next 48h at the given lat/lon using One Call API 3.0.
    Requires a One Call by Call subscription: https://openweathermap.org/api/one-call-3
    Returns a dict mapping UTC-epoch-hour → cloud_cover_percent.
    """
    url = "https://api.openweathermap.org/data/3.0/onecall"
    params = {
        "lat": lat,
        "lon": lon,
        "exclude": "current,minutely,daily,alerts",
        "appid": OWM_API_KEY,
        "units": "metric"
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()
    clouds = {}
    for hour in data.get("hourly", []):
        dt_hour = datetime.fromtimestamp(hour["dt"], tz=timezone.utc).replace(minute=0, second=0, microsecond=0)
        clouds[dt_hour] = hour.get("clouds", None)
    return clouds

def get_flyover_times_local(cfg, site, clouds):
    """
    Compute passes and annotate each event with cloud cover % (nearest hour).
    """
    ts = load.timescale()
    sat = EarthSatellite(cfg["tle"][0], cfg["tle"][1], 'SAT', ts)
    obs = wgs84.latlon(site["lat"], site["lon"], site.get("elevation_m", 0))

    # build UTC window
    d = cfg["date_utc"]
    start_utc = datetime(d.year, d.month, d.day, tzinfo=timezone.utc)
    end_utc   = start_utc + timedelta(days=1)
    t0, t1 = ts.from_datetime(start_utc), ts.from_datetime(end_utc)

    times, events = sat.find_events(obs, t0, t1,
                                     altitude_degrees=cfg["min_elevation_deg"])
    names = {0: 'rise', 1: 'culminate', 2: 'set'}
    tz = ZoneInfo(site["tz"])

    output = [
        f"\n=== {site['name']} on {d.isoformat()} ({site['tz']}) ===",
        f"Pointing accuracy threshold: {cfg['min_elevation_deg']}°",
        f"{'Event':10s}  {'Time':20s}  {'Clouds (%)':>10s}",
        "-"*50
    ]
    for t, e in zip(times, events):
        utc_dt   = t.utc_datetime().replace(tzinfo=timezone.utc)
        local_dt = utc_dt.astimezone(tz)
        # round down to the hour for matching forecast
        hour_key = utc_dt.replace(minute=0, second=0, microsecond=0)
        cloud_pct = clouds.get(hour_key, "N/A")
        output.append(f"{names[e]:10s}  {local_dt.strftime('%Y-%m-%d %H:%M')}  {str(cloud_pct):>10s}")
    return "\n".join(output)

if __name__ == '__main__':
    cfg = CONFIG
    for site in cfg["sites"]:
        cloud_map = fetch_cloud_forecast(site["lat"], site["lon"])
        print(get_flyover_times_local(cfg, site, cloud_map))
