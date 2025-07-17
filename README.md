# tle.py - Satellite TLE Visualization and Overfly Calculator

This script visualizes satellite orbits and calculates overfly times for selected cities using the latest TLE data from SatNOGS. It also displays detailed satellite information from SatNOGS and plots the satellite's ground track on a retro-style map.

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

## Notes
- The script will always show the map and satellite info, even if no city is selected.
- If SatNOGS is unavailable, it will fall back to the local `latest_tle.txt` file.
- For best results, keep `latest_tle.txt` up to date if using fallback mode. 