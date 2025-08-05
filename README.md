# Satellite TLE Tools - Visualization, Overfly, and Sun Exposure Calculators

This repository contains tools for satellite orbit analysis using TLE data from SatNOGS:

1. `tle.py` - Visualizes satellite orbits and calculates overfly times for selected cities
2. `sun_exposure.py` - Calculates when a satellite is exposed to sunlight vs. in Earth's shadow
3. `flyovers.py` - Calculates satellite pass times for specific locations with cloud cover forecasts

## Features
- Fetches the latest TLE and satellite info from SatNOGS (by NORAD ID or SatNOGS UUID)
- Falls back to a local TLE file if SatNOGS is unavailable
- Prints all available SatNOGS satellite info
- Calculates next overfly for a selected city (with configurable radius)
- Plots the satellite's ground track and current position on a map
- Supports multiple cities and overfly radius options
- Calculates satellite flyover times with cloud cover forecasts using OpenWeatherMap API

## Requirements
Install dependencies with:

```
pip install -r requirements.txt
```

## Usage

```
python tle.py [--city CITY] [--norad-id NORAD_ID] [--satnogs-id SATNOGS_ID] [--overfly-radius RADIUS]
```

### Command-line Arguments
- `--city CITY`           : City for overfly calculation (choose from: Santa Barbara, Helgoland, Berlin, Vienna, Milan)
- `--norad-id NORAD_ID`   : NORAD Catalog ID of the satellite (default: 98581)
- `--satnogs-id SATNOGS_ID`: SatNOGS internal satellite ID/UUID (default: JLJE-0670-3801-0857-8118)
- `--overfly-radius RADIUS`: Overfly radius in km (choices: 500, 1000, 2000; default: 1000)

### Overfly Radius Options
- `500 km`   : High elevation (best for communication)
- `1000 km`  : Practical visibility (default)
- `2000 km`  : Broad overfly (includes low elevation)

### Example Usage
- Show map and SatNOGS info for default satellite:
  ```
  python tle.py
  ```
- Calculate next overfly for Berlin with default radius:
  ```
  python tle.py --city Berlin
  ```
- Use a different NORAD ID:
  ```
  python tle.py --norad-id 25544
  ```
- Use a different SatNOGS UUID:
  ```
  python tle.py --satnogs-id JLJE-0670-3801-0857-8118
  ```
- Use a custom overfly radius:
  ```
  python tle.py --city Berlin --overfly-radius 500
  ```

## Sun Exposure Calculator Usage (sun_exposure.py)

The sun_exposure.py script calculates when a satellite is in sunlight and when it's in Earth's shadow.

```
python sun_exposure.py --start-date START_DATE --end-date END_DATE [options]
```

### Command-line Arguments
- `--start-date START_DATE` : Start date in ISO format (YYYY-MM-DDTHH:MM:SS)
- `--end-date END_DATE`     : End date in ISO format (YYYY-MM-DDTHH:MM:SS)
- `--norad-id NORAD_ID`     : NORAD Catalog ID of the satellite
- `--satnogs-id SATNOGS_ID` : SatNOGS internal satellite ID/UUID
- `--tle-file TLE_FILE`     : Path to a file containing TLE data
- `--interval INTERVAL`     : Time interval in seconds for analysis (default: 60)
- `--output OUTPUT`         : Output JSON filename (default: sun_exposure_YYYYMMDD-YYYYMMDD.json)

### Example Usage
- Calculate sun exposure using a local TLE file:
  ```
  python sun_exposure.py --start-date 2025-06-24T00:00:00 --end-date 2025-06-24T12:00:00 --tle-file latest_tle.txt
  ```
- Calculate sun exposure using NORAD ID:
  ```
  python sun_exposure.py --start-date 2025-06-24T00:00:00 --end-date 2025-06-24T12:00:00 --norad-id 98581
  ```
- Use a different time interval (e.g., 30 seconds):
  ```
  python sun_exposure.py --start-date 2025-06-24T00:00:00 --end-date 2025-06-24T12:00:00 --tle-file latest_tle.txt --interval 30
  ```
- By default, the output filename includes the date range:
  ```
  python sun_exposure.py --start-date 2025-06-24T00:00:00 --end-date 2025-06-24T12:00:00 --tle-file latest_tle.txt
  # Creates: sun_exposure_20250624-20250624.json
  ```
- Save output to a custom filename:
  ```
  python sun_exposure.py --start-date 2025-06-24T00:00:00 --end-date 2025-06-24T12:00:00 --tle-file latest_tle.txt --output my_satellite_exposure.json
  ```

### Output Format
The script generates a JSON file with the following structure:

```json
{
    "sun": [
        ["2025-06-24T00:12:00+00:00", "2025-06-24T01:12:00+00:00"],
        ["2025-06-24T01:47:00+00:00", "2025-06-24T02:47:00+00:00"]
    ],
    "dark": [
        ["2025-06-24T00:00:00+00:00", "2025-06-24T00:12:00+00:00"],
        ["2025-06-24T01:12:00+00:00", "2025-06-24T01:47:00+00:00"]
    ]
}
```

The CLI output provides a more readable format with duration information:
```
Sample sun exposure periods:
  1. 2025-06-24 00:12:00 → 01:12:00 (duration: 1h 0m 0s)
  2. 2025-06-24 01:47:00 → 02:47:00 (duration: 1h 0m 0s)

Sample shadow periods:
  1. 2025-06-24 00:00:00 → 00:12:00 (duration: 0h 12m 0s)
  2. 2025-06-24 01:12:00 → 01:47:00 (duration: 0h 35m 0s)

Total duration analyzed: 12h 0m 0s
Time in sunlight: 7h 43m 0s (64.3%)
Time in shadow: 4h 17m 0s (35.7%)
```

## Flyover Calculator Usage (flyovers.py)

The `flyovers.py` script calculates satellite passes for specific locations and optionally provides cloud cover data for each pass, allowing you to plan observations when skies are clear. Cloud cover data is typically only available for approximately 48 hours into the future.

### Features
- Calculates satellite rise, culmination, and set times for configured locations
- Integrates with OpenWeatherMap API to fetch cloud cover forecasts for each site (optional)
- Displays pass times in local time zones for easy planning
- Configurable elevation threshold based on your imaging pointing accuracy

### Configuration
The script uses a configuration dictionary that can be modified directly in the script:

```
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
```

### Requirements
- Python with the skyfield library and other dependencies
- **Optional**: An OpenWeatherMap API key with access to the One Call API 3.0 (for cloud cover data)
  - You need the "One Call by Call" subscription: https://openweathermap.org/api/one-call-3
  - Set your API key as an environment variable: `export OWM_API_KEY=your_api_key_here`
  - The script will still work without an API key, but won't show cloud cover data
  - Cloud forecasts are typically only available for ~48 hours ahead; passes beyond this window will show "N/A*"

### Usage
Run the script directly:

```
python flyovers.py
```

### Sample Output

```
=== Vienna on 2025-07-23 (Europe/Vienna) ===
Pointing accuracy threshold: 0.5°
Event       Time                  Clouds (%)
--------------------------------------------------
rise        2025-07-23 03:15           10
culminate   2025-07-23 03:18           10
set         2025-07-23 03:22           10
rise        2025-07-23 04:52           15
culminate   2025-07-23 04:57           15
set         2025-07-23 05:02           15
rise        2025-07-23 14:30           N/A*
culminate   2025-07-23 14:36           N/A*
set         2025-07-23 14:41           N/A*

Note: N/A* indicates times outside the available weather forecast window
```

## Notes
- Both scripts will fall back to the local `latest_tle.txt` file if SatNOGS is unavailable.
- For best results, keep `latest_tle.txt` up to date if using fallback mode.
- The sun exposure calculations use a simplified model of Earth's shadow.
- The `flyovers.py` script will calculate flyovers without an API key, but requires an OpenWeatherMap API key for cloud cover forecasts.
- By default, `flyovers.py` calculates passes for tomorrow's date, as weather forecasts are typically only available for ~48 hours ahead.