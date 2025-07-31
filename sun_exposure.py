#!/usr/bin/env python3
"""
Script to calculate satellite sun exposure times.
Determines when a satellite is in sunlight or in Earth's shadow.
"""

import argparse
import json
import os
from datetime import datetime, timezone, timedelta
from sgp4.api import Satrec, jday
from sgp4.earth_gravity import wgs84
import numpy as np
import requests
from pprint import pprint


def read_tle_from_file(filename):
    """Read TLE data from a file."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # The TLE data already includes the line numbers, so just strip whitespace
    line1 = lines[0].strip()
    line2 = lines[1].strip()

    return line1, line2


def fetch_tle_and_info_from_satnogs(norad_id=None, satnogs_id=None):
    """Fetch the latest TLE and satellite info for a given NORAD ID or SatNOGS satellite ID from SatNOGS API."""
    tle_url = None
    info_url = None
    if norad_id:
        tle_url = f'https://db.satnogs.org/api/tle/?norad_cat_id={norad_id}'
        info_url = f'https://db.satnogs.org/api/satellites/?norad_cat_id={norad_id}'
    elif satnogs_id:
        tle_url = f'https://db.satnogs.org/api/tle/?satellite={satnogs_id}'
        info_url = f'https://db.satnogs.org/api/satellites/{satnogs_id}'
    else:
        return None, None, None

    tle1, tle2, sat_info = None, None, None
    try:
        response = requests.get(tle_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data:
                tle1, tle2 = data[0]['tle1'], data[0]['tle2']
    except Exception as e:
        print(f"Error fetching TLE from SatNOGS: {e}")

    try:
        response = requests.get(info_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            # /api/satellites/?norad_cat_id=... returns a list, /api/satellites/<uuid> returns a dict
            if isinstance(data, list) and data:
                sat_info = data[0]
            elif isinstance(data, dict) and data:
                sat_info = data
    except Exception as e:
        print(f"Error fetching satellite info from SatNOGS: {e}")

    return tle1, tle2, sat_info


def eci_to_geodetic(x, y, z, jd):
    """
    Convert Earth-centered inertial coordinates to latitude, longitude, and altitude.

    Args:
        x, y, z: ECI coordinates in km
        jd: Julian date

    Returns:
        latitude (degrees), longitude (degrees), altitude (km)
    """
    from math import sqrt, atan2, asin, pi, sin, cos

    # WGS-84 ellipsoid constants
    a = 6378.137  # equatorial radius in km
    f = 1/298.257223563  # flattening
    b = a*(1-f)  # polar radius in km
    e2 = 1 - (b*b)/(a*a)  # square of eccentricity

    # Calculate GMST (Greenwich Mean Sidereal Time)
    tut1 = (jd - 2451545.0) / 36525.0
    # Simplified GMST calculation
    gmst = 67310.54841 + (876600.0*3600 + 8640184.812866)*tut1 + 0.093104*tut1*tut1 - 6.2e-6*tut1*tut1*tut1
    # Convert to degrees and normalize to 0-360
    gmst = (gmst % 86400.0) * 360.0 / 86400.0

    # Calculate position vector magnitude
    r = sqrt(x*x + y*y + z*z)

    # Calculate longitude
    lon = atan2(y, x) * 180.0/pi - gmst
    while lon > 180.0:
        lon -= 360.0
    while lon < -180.0:
        lon += 360.0

    # Calculate latitude (approximation)
    lat = asin(z/r) * 180.0/pi

    # Iterative calculation for more accurate latitude and altitude
    phi = lat * pi/180.0
    c = 0.0
    for i in range(20):  # Usually converges in a few iterations
        sin_phi = sin(phi)
        c = 1.0 / sqrt(1.0 - e2*sin_phi*sin_phi)
        phi_new = atan2(z + a*c*e2*sin_phi, sqrt(x*x + y*y))
        if abs(phi - phi_new) < 1e-10:
            break
        phi = phi_new

    lat = phi * 180.0/pi
    alt = sqrt(x*x + y*y) / cos(phi) - a*c

    return lat, lon, alt


def is_in_sunlight(pos, jd):
    """
    Determine if a satellite is in sunlight or Earth's shadow.

    Args:
        pos: satellite position vector in ECI coordinates (km)
        jd: Julian date

    Returns:
        True if in sunlight, False if in Earth's shadow
    """
    # Earth's radius in km
    earth_radius = 6378.137

    # Calculate Sun position vector in ECI
    # Use simplified solar position algorithm
    t = jd - 2451545.0  # Julian days since J2000

    # Mean longitude of the Sun
    L = (280.460 + 0.9856474 * t) % 360

    # Mean anomaly of the Sun
    g = (357.528 + 0.9856003 * t) % 360

    # Ecliptic longitude of the Sun
    lambda_sun = L + 1.915 * np.sin(np.radians(g)) + 0.020 * np.sin(np.radians(2*g))

    # Obliquity of the ecliptic
    epsilon = 23.439 - 0.0000004 * t

    # Sun position in ECI (simplified, assuming 1 AU distance)
    sun_dist = 149597870.7  # km
    x_sun = sun_dist * np.cos(np.radians(lambda_sun))
    y_sun = sun_dist * np.sin(np.radians(lambda_sun)) * np.cos(np.radians(epsilon))
    z_sun = sun_dist * np.sin(np.radians(lambda_sun)) * np.sin(np.radians(epsilon))

    sun_vector = np.array([x_sun, y_sun, z_sun])
    sat_vector = np.array(pos)

    # Normalize the sun vector
    sun_unit = sun_vector / np.linalg.norm(sun_vector)

    # Calculate angle between satellite position and sun
    sat_to_sun = sun_vector - sat_vector
    sat_to_sun_dist = np.linalg.norm(sat_to_sun)

    # Dot product between satellite position and sun direction
    dot_product = np.dot(sat_vector, sun_unit)

    # Calculate the perpendicular distance from satellite to sun-Earth line
    perp_dist = np.linalg.norm(sat_vector - dot_product * sun_unit)

    # Check if satellite is in Earth's shadow
    if dot_product < 0 and perp_dist < earth_radius:
        return False
    else:
        return True


def calculate_sun_exposure(satellite, start_time, end_time, interval_seconds=60):
    """
    Calculate times when a satellite is in sunlight and in Earth's shadow.

    Args:
        satellite: SGP4 satellite object
        start_time: datetime object for start of analysis
        end_time: datetime object for end of analysis
        interval_seconds: time step for analysis in seconds

    Returns:
        Dictionary with "sun" and "dark" lists containing timestamp ranges
    """
    current_time = start_time
    in_sunlight = None
    sun_periods = []
    dark_periods = []
    period_start = None

    while current_time <= end_time:
        # Convert current time to Julian date for SGP4
        year, month, day = current_time.year, current_time.month, current_time.day
        hour, minute, second = current_time.hour, current_time.minute, current_time.second
        jd, fr = jday(year, month, day, hour, minute, second)

        # Get satellite position
        err, pos, vel = satellite.sgp4(jd, fr)

        if err != 0:
            print(f"Error calculating position at {current_time}")
            current_time += timedelta(seconds=interval_seconds)
            continue

        # Check if satellite is in sunlight
        sun_status = is_in_sunlight(pos, jd + fr)

        # Initialize the status if this is the first iteration
        if in_sunlight is None:
            in_sunlight = sun_status
            period_start = current_time
        # If the status has changed, record the period
        elif in_sunlight != sun_status:
            period_end = current_time
            if in_sunlight:
                sun_periods.append([period_start.isoformat(), period_end.isoformat()])
            else:
                dark_periods.append([period_start.isoformat(), period_end.isoformat()])
            period_start = current_time
            in_sunlight = sun_status

        current_time += timedelta(seconds=interval_seconds)

    # Add the final period if it exists
    if period_start is not None and period_start < current_time:
        if in_sunlight:
            sun_periods.append([period_start.isoformat(), end_time.isoformat()])
        else:
            dark_periods.append([period_start.isoformat(), end_time.isoformat()])

    return {
        "sun": sun_periods,
        "dark": dark_periods
    }


def main():
    parser = argparse.ArgumentParser(description='Calculate satellite sun exposure.')
    parser.add_argument('--start-date', type=str, required=True,
                        help='Start date in ISO format (YYYY-MM-DDTHH:MM:SS)')
    parser.add_argument('--end-date', type=str, required=True,
                        help='End date in ISO format (YYYY-MM-DDTHH:MM:SS)')
    parser.add_argument('--norad-id', type=int, help='NORAD ID of the satellite')
    parser.add_argument('--satnogs-id', type=str, help='SatNOGS ID of the satellite')
    parser.add_argument('--tle-file', type=str, help='File containing TLE data')
    parser.add_argument('--interval', type=int, default=60,
                        help='Time interval in seconds for analysis (default: 60)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output file name (default: sun_exposure_YYYYMMDD-YYYYMMDD.json)')



    args = parser.parse_args()

    # Parse dates
    try:
        start_time = datetime.fromisoformat(args.start_date.replace('Z', '+00:00')).replace(tzinfo=timezone.utc)
        end_time = datetime.fromisoformat(args.end_date.replace('Z', '+00:00')).replace(tzinfo=timezone.utc)
    except ValueError:
        print("Error: Invalid date format. Please use ISO format (YYYY-MM-DDTHH:MM:SS)")
        return

    # Get TLE data
    line1, line2 = None, None
    sat_info = None

    if args.tle_file:
        line1, line2 = read_tle_from_file(args.tle_file)
    elif args.norad_id:
        line1, line2, sat_info = fetch_tle_and_info_from_satnogs(norad_id=args.norad_id)
    elif args.satnogs_id:
        line1, line2, sat_info = fetch_tle_and_info_from_satnogs(satnogs_id=args.satnogs_id)
    else:
        print("Error: Please provide either a TLE file, NORAD ID, or SatNOGS ID")
        return

    if not line1 or not line2:
        print("Error: Failed to retrieve TLE data")
        return

    # Create satellite object
    satellite = Satrec.twoline2rv(line1, line2)

    # Calculate sun exposure
    print(f"Calculating sun exposure from {start_time} to {end_time}...")
    exposure_data = calculate_sun_exposure(satellite, start_time, end_time, args.interval)

    # Create output directory if it doesn't exist
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
    os.makedirs(out_dir, exist_ok=True)

    # Generate output filename with timestamps if not provided
    if args.output is None:
        start_date_str = start_time.strftime("%Y%m%d")
        end_date_str = end_time.strftime("%Y%m%d")
        output_filename = f"sun_exposure_{start_date_str}-{end_date_str}.json"
    else:
        output_filename = args.output

    # Save results to file
    output_path = os.path.join(out_dir, output_filename)
    with open(output_path, 'w') as f:
        json.dump(exposure_data, f, indent=4)

    print(f"Results saved to {output_path}")

    # Pretty print results
    print("\nSummary of Results:")
    print(f"Total sun exposure periods: {len(exposure_data['sun'])}")
    print(f"Total shadow periods: {len(exposure_data['dark'])}")

    # Print first few periods as examples with readable timestamps and duration
    print("\nSample sun exposure periods:")
    for i, period in enumerate(exposure_data['sun'][:3]):
        start = datetime.fromisoformat(period[0])
        end = datetime.fromisoformat(period[1])
        duration = end - start
        hours, remainder = divmod(duration.total_seconds(), 3600)
        minutes, seconds = divmod(remainder, 60)
        duration_str = f"{int(hours)}h {int(minutes)}m {int(seconds)}s"
        start_str = start.strftime("%Y-%m-%d %H:%M:%S")
        end_str = end.strftime("%H:%M:%S")  # Just time for end to make it more compact
        print(f"  {i+1}. {start_str} → {end_str} (duration: {duration_str})")

    print("\nSample shadow periods:")
    for i, period in enumerate(exposure_data['dark'][:3]):
        start = datetime.fromisoformat(period[0])
        end = datetime.fromisoformat(period[1])
        duration = end - start
        hours, remainder = divmod(duration.total_seconds(), 3600)
        minutes, seconds = divmod(remainder, 60)
        duration_str = f"{int(hours)}h {int(minutes)}m {int(seconds)}s"
        start_str = start.strftime("%Y-%m-%d %H:%M:%S")
        end_str = end.strftime("%H:%M:%S")  # Just time for end to make it more compact
        print(f"  {i+1}. {start_str} → {end_str} (duration: {duration_str})")

    # Calculate total exposure time
    sun_duration = timedelta()
    for period in exposure_data['sun']:
        start = datetime.fromisoformat(period[0])
        end = datetime.fromisoformat(period[1])
        sun_duration += end - start

    dark_duration = timedelta()
    for period in exposure_data['dark']:
        start = datetime.fromisoformat(period[0])
        end = datetime.fromisoformat(period[1])
        dark_duration += end - start

    total_duration = end_time - start_time
    sun_percentage = (sun_duration.total_seconds() / total_duration.total_seconds()) * 100

    # Format duration strings
    def format_duration(td):
        hours, remainder = divmod(td.total_seconds(), 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{int(hours)}h {int(minutes)}m {int(seconds)}s"

    total_str = format_duration(total_duration)
    sun_str = format_duration(sun_duration)
    dark_str = format_duration(dark_duration)

    print(f"\nTotal duration analyzed: {total_str}")
    print(f"Time in sunlight: {sun_str} ({sun_percentage:.1f}%)")
    print(f"Time in shadow: {dark_str} ({100-sun_percentage:.1f}%)")


if __name__ == "__main__":
    main()
