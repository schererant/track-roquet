# tle.py - Satellite TLE Analysis and Visualization

This project provides a Python script (`tle.py`) for analyzing satellite Two-Line Element (TLE) data, visualizing the satellite's ground track, and calculating the next overfly time for selected cities.

## Features
- Reads TLE data from a file (`latest_tle.txt`) or uses sample ISS data as fallback
- Analyzes orbital parameters (inclination, eccentricity, period, altitude, etc.)
- Calculates and displays the next time the satellite will overfly a selected city
- Plots the satellite's ground track (previous 90 min and next 90 min) on a retro-styled world map
- Supports multiple cities for overfly calculation

## Requirements
- Python 3.7+
- The following Python packages:
  - `sgp4`
  - `numpy`
  - `matplotlib`
  - `cartopy`

Install dependencies with:
```bash
pip install -r requirements.txt
```

## Usage

### 1. Prepare TLE Data
Place your TLE data in a file named `latest_tle.txt` in the project directory. The file should contain two lines:
```
<line1>
<line2>
```
If `latest_tle.txt` is missing, the script will use sample ISS (ZARYA) TLE data.

### 2. Run the Script
You can run the script from the command line:
```bash
python tle.py
```
By default, the script calculates the next overfly for Santa Barbara, USA.

#### To select a different city:
```bash
python tle.py --city Berlin
python tle.py --city Vienna
python tle.py --city Helgoland
python tle.py --city Milan
python tle.py --city "Santa Barbara"
```

### 3. Output
- **Orbit Analysis:** Prints orbital parameters and type.
- **Current Position:** Prints ECI position, velocity, speed, and geodetic coordinates.
- **City Overfly Calculation:** Prints the next time the satellite will be closest to the selected city (within 2000 km), including time, distance, and satellite position/speed at that moment.
- **Map Visualization:** Opens a retro-styled world map showing:
  - The satellite's previous and future ground track (±90 min)
  - All available cities (highlighting the selected one)
  - The satellite's current position

## Available Cities
- Santa Barbara, USA
- Helgoland, Germany
- Berlin, Germany
- Vienna, Austria
- Milan, Italy

## Example
```
$ python tle.py --city Berlin
Using TLE data from latest_tle.txt:
Line 1: ...
Line 2: ...
==================================================
ORBIT ANALYSIS
...
==================================================
SCRIPT USAGE INFORMATION
...
==================================================
CITY OVERFLY CALCULATION
Available cities for overfly calculation:
- Santa Barbara, USA (34.4208°, -119.6982°)
- Helgoland, Germany (54.1833°, 7.8833°)
- Berlin, Germany (52.52°, 13.405°)
...
Next Berlin overfly:
  Time: 2024-06-10 12:34:56 UTC
  In: 2 hours, 15 minutes
  Closest distance: 1234.5 km
  Satellite position at overfly:
    Latitude: 52.123°, Longitude: 13.456°, Altitude: 420.1 km
    Speed: 7.670 km/s (27,612.0 km/h)
==================================================
```
A map window will open showing the satellite's ground track and city locations.

## Notes
- The script uses a simple proximity check (within 2000 km) for overfly calculations.
- The map uses a retro style with green and orange highlights.
- If you encounter issues with `cartopy` installation, refer to [Cartopy installation docs](https://scitools.org.uk/cartopy/docs/latest/installing.html).

## License
MIT License 