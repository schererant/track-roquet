# Satellite TLE Tools - Visualization, Overfly, and Sun Exposure Calculators

This repository contains tools for satellite orbit analysis using TLE data from SatNOGS:

1. `tle.py` - Visualizes satellite orbits and calculates overfly times for selected cities
2. `sun_exposure.py` - Calculates when a satellite is exposed to sunlight vs. in Earth's shadow

## Features
- Fetches the latest TLE and satellite info from SatNOGS (by NORAD ID or SatNOGS UUID)
- Falls back to a local TLE file if SatNOGS is unavailable
- Prints all available SatNOGS satellite info
- Calculates next overfly for a selected city (with configurable radius)
- Plots the satellite's ground track and current position on a map
- Supports multiple cities and overfly radius options

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

## Notes
- Both scripts will fall back to the local `latest_tle.txt` file if SatNOGS is unavailable.
- For best results, keep `latest_tle.txt` up to date if using fallback mode.
- The sun exposure calculations use a simplified model of Earth's shadow.