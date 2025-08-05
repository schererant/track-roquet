import os
from skyfield.api import load, EarthSatellite, wgs84
from datetime import datetime, timedelta, timezone, date
from zoneinfo import ZoneInfo
import requests

# ───── CONFIGURATION ──────────────────────────────────────────────────────────
OWM_API_KEY = os.environ.get("OWM_API_KEY", "")

CONFIG = {
    # UTC date to examine passes for
    "date_utc": date(2025, 8, 7),

    # elevation threshold = your imaging pointing accuracy (degrees)
    "min_elevation_deg": 0.5,

    # sites to evaluate: name, (lat, lon), and local timezone
    "sites": [
        {"name": "Vienna",          "lat": 48.2082, "lon": 16.3738,   "tz": "Europe/Vienna", "elevation_m": 170},
        {"name": "Munich",          "lat": 48.1351, "lon": 11.5820,   "tz": "Europe/Berlin", "elevation_m": 520},
        {"name": "Santa Barbara",   "lat": 34.4208, "lon": -119.6982, "tz": "America/Los_Angeles", "elevation_m": 15},
        {"name": "Milan",           "lat": 45.4642, "lon": 9.1900,    "tz": "Europe/Rome", "elevation_m": 120},
        {"name": "Helgoland",       "lat": 54.1833, "lon": 7.9167,    "tz": "Europe/Berlin", "elevation_m": 4},
        {"name": "Istanbul",        "lat": 41.0082, "lon": 28.9784,   "tz": "Europe/Istanbul", "elevation_m": 40},
        {"name": "Tokyo",           "lat": 35.6762, "lon": 139.6503,  "tz": "Asia/Tokyo", "elevation_m": 40},
        {"name": "Dubai",           "lat": 25.2048, "lon": 55.2708,   "tz": "Asia/Dubai", "elevation_m": 5},
        {"name": "Sydney",          "lat": -33.8688, "lon": 151.2093, "tz": "Australia/Sydney", "elevation_m": 58},
        {"name": "Brasilia",        "lat": -15.8267, "lon": -47.9218, "tz": "America/Sao_Paulo", "elevation_m": 1172},
        {"name": "New York",        "lat": 40.7128, "lon": -74.0060,  "tz": "America/New_York", "elevation_m": 10},
        {"name": "Auckland",        "lat": -36.8509, "lon": 174.7645, "tz": "Pacific/Auckland", "elevation_m": 23},
        {"name": "Brisbane",        "lat": -27.4698, "lon": 153.0251, "tz": "Australia/Brisbane", "elevation_m": 28},
        {"name": "Singapore",       "lat": 1.3521,  "lon": 103.8198,  "tz": "Asia/Singapore", "elevation_m": 15},
        {"name": "Paris",           "lat": 48.8566, "lon": 2.3522,    "tz": "Europe/Paris", "elevation_m": 35},
        {"name": "Porto",           "lat": 41.1579, "lon": -8.6291,   "tz": "Europe/Lisbon", "elevation_m": 104},
        {"name": "Portland",        "lat": 45.5152, "lon": -122.6784, "tz": "America/Los_Angeles", "elevation_m": 15},
        {"name": "Seattle",         "lat": 47.6062, "lon": -122.3321, "tz": "America/Los_Angeles", "elevation_m": 53},
        {"name": "San Francisco",   "lat": 37.7749, "lon": -122.4194, "tz": "America/Los_Angeles", "elevation_m": 16},

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
    Returns empty dict if API key is missing or request fails.
    """
    if not OWM_API_KEY:
        print("Warning: No OpenWeatherMap API key found. Set OWM_API_KEY environment variable for cloud data.")
        return {}

    url = "https://api.openweathermap.org/data/3.0/onecall"
    params = {
        "lat": lat,
        "lon": lon,
        "exclude": "current,minutely,daily,alerts",
        "appid": OWM_API_KEY,
        "units": "metric"
    }

    try:
        resp = requests.get(url, params=params)
        resp.raise_for_status()
        data = resp.json()
        clouds = {}
        for hour in data.get("hourly", []):
            dt_hour = datetime.fromtimestamp(hour["dt"], tz=timezone.utc).replace(minute=0, second=0, microsecond=0)
            clouds[dt_hour] = hour.get("clouds", None)

        # Print forecast range info for debugging
        if clouds:
            start = min(clouds.keys()).strftime('%Y-%m-%d %H:%M UTC')
            end = max(clouds.keys()).strftime('%Y-%m-%d %H:%M UTC')
            print(f"Cloud forecast available from {start} to {end}")
        return clouds
    except requests.exceptions.RequestException as e:
        print(f"Warning: Cloud data fetch failed: {e}")
        return {}

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
        f"{'Event':10s}  {'Local Time':20s}  {'UTC Time':20s}  {'Clouds (%)':>10s}",
        "-"*80
    ]
    for t, e in zip(times, events):
        utc_dt   = t.utc_datetime().replace(tzinfo=timezone.utc)
        local_dt = utc_dt.astimezone(tz)

        # Try to find the closest hour for cloud data
        # First try exact hour match
        hour_key = utc_dt.replace(minute=0, second=0, microsecond=0)
        cloud_pct = clouds.get(hour_key)

        # If no exact match, look for closest hour within 1 hour (for better forecast coverage)
        if cloud_pct is None and clouds:
            for offset in [-1, 1]:
                nearby_hour = hour_key + timedelta(hours=offset)
                if nearby_hour in clouds:
                    cloud_pct = clouds[nearby_hour]
                    break

        cloud_str = str(cloud_pct) if cloud_pct is not None else "N/A"
        if cloud_str == "N/A" and clouds:
            cloud_str = "N/A*"  # Mark as outside forecast window

        output.append(f"{names[e]:10s}  {local_dt.strftime('%Y-%m-%d %H:%M'):20s}  {utc_dt.strftime('%Y-%m-%d %H:%M'):20s}  {cloud_str:>10s}")

    # Add a note about N/A values if we have some cloud data but not for all passes
    if clouds and any("N/A*" in line for line in output):
        output.append("")
        output.append("Note: N/A* indicates times outside the available weather forecast window")

    return "\n".join(output)

if __name__ == '__main__':
    cfg = CONFIG

    # Print date being analyzed
    print(f"Analyzing satellite passes for {cfg['date_utc'].isoformat()}")
    print("Cloud forecasts typically only available for 48 hours ahead")
    print("-" * 70)

    for site in cfg["sites"]:
        cloud_map = fetch_cloud_forecast(site["lat"], site["lon"])
        print(get_flyover_times_local(cfg, site, cloud_map))
        print()  # Add a blank line between sites for better readability
