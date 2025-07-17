from sgp4.api import Satrec
from sgp4.api import jday
from datetime import datetime, timezone, timedelta
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from matplotlib.colors import LinearSegmentedColormap
import argparse

# Read current TLE data from file
def read_tle_from_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # The TLE data already includes the line numbers, so just strip whitespace
    line1 = lines[0].strip()
    line2 = lines[1].strip()

    return line1, line2

# Read TLE data from file
try:
    line1, line2 = read_tle_from_file('latest_tle.txt')
    print(f"Using TLE data from latest_tle.txt:")
    print(f"Line 1: {line1}")
    print(f"Line 2: {line2}")
except FileNotFoundError:
    print("latest_tle.txt not found, using sample ISS (ZARYA) TLE data")
    # Fallback to sample data
    line1 = "1 98580C 25135BL  25174.99035276 -.00080323  00000+0 -80437-2 0 00"
    line2 = "2 98580  97.7470 289.3930 0004228  87.9095 217.7754 14.91521210 14"

# City coordinates dictionary - available globally
CITIES = {
    "Santa Barbara": {"lat": 34.4208, "lon": -119.6982, "country": "USA"},
    "Helgoland": {"lat": 54.1833, "lon": 7.8833, "country": "Germany"},
    "Berlin": {"lat": 52.5200, "lon": 13.4050, "country": "Germany"},
    "Vienna": {"lat": 48.2082, "lon": 16.3738, "country": "Austria"},
    "Milan": {"lat": 45.4642, "lon": 9.1900, "country": "Italy"}
}

def analyze_orbit_from_tle(line1, line2):
    """Analyze TLE data to determine orbit type and characteristics"""
    # Parse TLE line 2 for orbital elements
    # Format: 2 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
    #         ^^^^^^^ ^^^^^^^^ ^^^^^^^^^^^^^^^^ ^^^^^^^^^^^ ^^^^^^^^ ^^^^^^^^ ^ ^^^^^^
    #         NORAD   Epoch    Mean Motion      B*          Incl     RAAN     Ecc  Arg Perigee  Mean Anomaly  Rev

    try:
        # Extract orbital elements
        inclination = float(line2[8:16])  # Inclination (degrees)
        raan = float(line2[17:25])        # Right Ascension of Ascending Node (degrees)
        eccentricity = float('0.' + line2[26:33])  # Eccentricity
        arg_perigee = float(line2[34:42]) # Argument of Perigee (degrees)
        mean_anomaly = float(line2[43:51]) # Mean Anomaly (degrees)
        mean_motion = float(line2[52:63])  # Mean Motion (revolutions per day)

        # Calculate orbital period and semi-major axis
        period_minutes = 1440 / mean_motion  # 1440 minutes in a day
        period_hours = period_minutes / 60

        # Calculate altitude (approximate)
        # For circular orbits: altitude ≈ (period_minutes/84.4)^(2/3) - 6378 km
        # This is a rough approximation
        altitude_km = (period_minutes/84.4)**(2/3) - 6378

        # Determine orbit type
        orbit_type = "Unknown"
        if inclination < 30:
            orbit_type = "Low Inclination (Equatorial)"
        elif inclination < 60:
            orbit_type = "Medium Inclination"
        elif inclination < 90:
            orbit_type = "High Inclination"
        elif inclination < 110:
            orbit_type = "Polar"
        else:
            orbit_type = "Retrograde"

        # Check for special orbit types
        if abs(inclination - 63.4) < 5:
            orbit_type += " (Molniya-type)"
        elif abs(inclination - 90) < 5:
            orbit_type += " (Sun-synchronous)"
        elif eccentricity > 0.1:
            orbit_type += " (Elliptical)"
        elif eccentricity < 0.01:
            orbit_type += " (Near-circular)"

        return {
            'inclination': inclination,
            'eccentricity': eccentricity,
            'period_minutes': period_minutes,
            'period_hours': period_hours,
            'altitude_km': altitude_km,
            'orbit_type': orbit_type,
            'raan': raan,
            'arg_perigee': arg_perigee,
            'mean_anomaly': mean_anomaly,
            'mean_motion': mean_motion
        }
    except:
        return None

# Create satellite record
satellite = Satrec.twoline2rv(line1, line2)

# Analyze orbit characteristics
print("\n" + "="*50)
print("ORBIT ANALYSIS")
print("="*50)
orbit_info = analyze_orbit_from_tle(line1, line2)
if orbit_info:
    print(f"Orbit Type: {orbit_info['orbit_type']}")
    print(f"Inclination: {orbit_info['inclination']:.2f}°")
    print(f"Eccentricity: {orbit_info['eccentricity']:.6f}")
    print(f"Orbital Period: {orbit_info['period_minutes']:.1f} minutes ({orbit_info['period_hours']:.2f} hours)")
    print(f"Mean Motion: {orbit_info['mean_motion']:.6f} rev/day")
    print(f"Approximate Altitude: {orbit_info['altitude_km']:.0f} km")
    print(f"Right Ascension of Ascending Node: {orbit_info['raan']:.2f}°")
    print(f"Argument of Perigee: {orbit_info['arg_perigee']:.2f}°")
    print(f"Mean Anomaly: {orbit_info['mean_anomaly']:.2f}°")
else:
    print("Could not analyze orbit from TLE data")
print("="*50 + "\n")

# Get current UTC time (timezone-aware)
now = datetime.now(timezone.utc)
yr, mo, dy, hr, mi, se = now.year, now.month, now.day, now.hour, now.minute, now.second
jd, fr = jday(yr, mo, dy, hr, mi, se)

# Propagate to current time
error_code, position, velocity = satellite.sgp4(jd, fr)

def eci_to_geodetic(x, y, z, jd):
    # WGS-84 constants
    a = 6378.137  # Equatorial radius in km
    f = 1 / 298.257223563  # Flattening
    b = a * (1 - f)
    e2 = 1 - (b**2 / a**2)

    # Calculate GMST (Greenwich Mean Sidereal Time)
    T = (jd - 2451545.0) / 36525.0
    GMST = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T**2 - (T**3) / 38710000.0
    GMST = np.deg2rad(GMST % 360)

    # Rotate ECI to ECEF
    x_ecef = x * np.cos(GMST) + y * np.sin(GMST)
    y_ecef = -x * np.sin(GMST) + y * np.cos(GMST)
    z_ecef = z

    # Longitude
    lon = np.arctan2(y_ecef, x_ecef)

    # Latitude (iterative)
    r = np.sqrt(x_ecef**2 + y_ecef**2)
    E2 = a**2 - b**2
    F = 54 * b**2 * z_ecef**2
    G = r**2 + (1 - e2) * z_ecef**2 - e2 * E2
    c = (e2**2 * F * r**2) / (G**3)
    s = np.cbrt(1 + c + np.sqrt(c**2 + 2 * c))
    P = F / (3 * (s + 1/s + 1)**2 * G**2)
    Q = np.sqrt(1 + 2 * e2**2 * P)
    r_0 = -(P * e2 * r) / (1 + Q) + np.sqrt(0.5 * a**2 * (1 + 1/Q) - P * (1 - e2) * z_ecef**2 / (Q * (1 + Q)) - 0.5 * P * r**2)
    U = np.sqrt((r - e2 * r_0)**2 + z_ecef**2)
    V = np.sqrt((r - e2 * r_0)**2 + (1 - e2) * z_ecef**2)
    Z_0 = b**2 * z_ecef / (a * V)
    lat = np.arctan((z_ecef + Z_0) / r)
    alt = U * (1 - b**2 / (a * V))

    lat_deg = np.degrees(lat)
    lon_deg = np.degrees(lon)
    return lat_deg, lon_deg, alt

def calculate_distance(lat1, lon1, lat2, lon2):
    """Calculate distance between two points on Earth's surface in km"""
    R = 6371  # Earth's radius in km
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    a = np.sin(dlat/2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = R * c

    return distance

def calculate_speed(velocity):
    """Calculate satellite speed in km/s and km/h"""
    speed_kmps = np.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2)
    speed_kmph = speed_kmps * 3600  # Convert to km/h
    return speed_kmps, speed_kmph

def find_next_city_overfly(satellite, start_time, city_name=None, max_hours=48):
    """Find the next time the satellite overflies a selected city

    Args:
        satellite: SGP4 satellite object
        start_time: datetime object representing the start time for the search
        city_name: Name of the city to find overfly for (must be in CITIES dict)
        max_hours: Maximum number of hours to search for an overfly

    Returns:
        Tuple of (overfly_time, closest_distance, city_name)

    Available cities:
    - Santa Barbara, USA (34.4208°N, 119.6982°W)
    - Helgoland, Germany (54.1833°N, 7.8833°E)
    - Berlin, Germany (52.5200°N, 13.4050°E)
    - Vienna, Austria (48.2082°N, 16.3738°E)
    - Milan, Italy (45.4642°N, 9.1900°E)
    """
    # City coordinates dictionary - this is also used elsewhere in the code
    # so it's defined at the global level
    global CITIES

    # Default to Santa Barbara if no city is specified or if city is not in the list
    if not city_name or city_name not in CITIES:
        city_name = "Santa Barbara"

    # Get the coordinates for the selected city
    city_lat = CITIES[city_name]["lat"]
    city_lon = CITIES[city_name]["lon"]

    # Search for the next 48 hours (every 5 minutes)
    closest_time = None
    closest_distance = float('inf')

    for minute in range(0, max_hours * 60, 5):  # Check every 5 minutes
        check_time = start_time + timedelta(minutes=minute)
        yr, mo, dy, hr, mi, se = check_time.year, check_time.month, check_time.day, check_time.hour, check_time.minute, check_time.second
        jd, fr = jday(yr, mo, dy, hr, mi, se)

        error_code, position, velocity = satellite.sgp4(jd, fr)
        if error_code == 0:
            sat_lat, sat_lon, alt = eci_to_geodetic(position[0], position[1], position[2], jd + fr)
            distance = calculate_distance(sat_lat, sat_lon, city_lat, city_lon)

            # Consider it an overfly if within 2000 km (much broader coverage)
            if distance < 2000 and distance < closest_distance:
                closest_distance = distance
                closest_time = check_time

    return closest_time, closest_distance, city_name

if error_code == 0:
    print(f"ECI position (km): x={position[0]:.3f}, y={position[1]:.3f}, z={position[2]:.3f}")
    print(f"ECI velocity (km/s): x={velocity[0]:.6f}, y={velocity[1]:.6f}, z={velocity[2]:.6f}")

    # Calculate current speed
    speed_kmps, speed_kmph = calculate_speed(velocity)
    print(f"Speed: {speed_kmps:.3f} km/s ({speed_kmph:.1f} km/h)")

    # Calculate current geodetic position
    lat, lon, alt = eci_to_geodetic(position[0], position[1], position[2], jd + fr)
    print(f"Latitude: {lat:.3f}°, Longitude: {lon:.3f}°, Altitude: {alt:.3f} km")

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='TLE Orbit Visualization and City Overfly Calculation')
    parser.add_argument('--city', type=str, choices=CITIES.keys(), default="Santa Barbara",
                        help='City for which to calculate satellite overfly')
    args = parser.parse_args()

    # Print script usage information
    print("\n" + "="*50)
    print("SCRIPT USAGE INFORMATION")
    print("="*50)
    print("This script analyzes satellite TLE data and calculates overflys for selected cities.")
    print("To select a different city, use the --city command-line argument:")
    print("  python tle.py --city Berlin")
    print("  python tle.py --city Vienna")
    print("  python tle.py --city Helgoland")
    print("  python tle.py --city Milan")
    print("  python tle.py --city \"Santa Barbara\"")
    print("="*50 + "\n")

    # Find next city overfly
    print("\n" + "="*50)
    print("CITY OVERFLY CALCULATION")
    print("="*50)
    # List of available cities
    print("Available cities for overfly calculation:")
    for city, data in CITIES.items():
        print(f"- {city}, {data['country']} ({data['lat']}°, {data['lon']}°)")

    # Use the city from command-line argument
    selected_city = args.city
    print(f"Using city: {selected_city}")

    next_overfly_time, closest_distance, city_name = find_next_city_overfly(satellite, now, selected_city)

    if next_overfly_time:
        time_until_overfly = next_overfly_time - now
        hours = int(time_until_overfly.total_seconds() // 3600)
        minutes = int((time_until_overfly.total_seconds() % 3600) // 60)

        print(f"Next {city_name} overfly:")
        print(f"  Time: {next_overfly_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
        print(f"  In: {hours} hours, {minutes} minutes")
        print(f"  Closest distance: {closest_distance:.1f} km")

        # Calculate satellite position and speed at overfly time
        o_yr, o_mo, o_dy, o_hr, o_mi, o_se = next_overfly_time.year, next_overfly_time.month, next_overfly_time.day, next_overfly_time.hour, next_overfly_time.minute, next_overfly_time.second
        o_jd, o_fr = jday(o_yr, o_mo, o_dy, o_hr, o_mi, o_se)
        o_err, o_pos, o_vel = satellite.sgp4(o_jd, o_fr)
        if o_err == 0:
            o_lat, o_lon, o_alt = eci_to_geodetic(o_pos[0], o_pos[1], o_pos[2], o_jd + o_fr)
            o_speed_kmps, o_speed_kmph = calculate_speed(o_vel)
            print(f"  Satellite position at overfly:")
            print(f"    Latitude: {o_lat:.3f}°, Longitude: {o_lon:.3f}°, Altitude: {o_alt:.1f} km")
            print(f"    Speed: {o_speed_kmps:.3f} km/s ({o_speed_kmph:.1f} km/h)")
    else:
        print(f"No {city_name} overfly found in the next 48 hours")
    print("="*50 + "\n")

    # Calculate trajectory for the next 90 minutes (every minute)
    lats, lons = [], []
    for minute in range(0, 91):
        future_time = now + timedelta(minutes=minute)
        f_yr, f_mo, f_dy, f_hr, f_mi, f_se = future_time.year, future_time.month, future_time.day, future_time.hour, future_time.minute, future_time.second
        f_jd, f_fr = jday(f_yr, f_mo, f_dy, f_hr, f_mi, f_se)
        err, pos, vel = satellite.sgp4(f_jd, f_fr)
        if err == 0:
            plat, plon, _ = eci_to_geodetic(pos[0], pos[1], pos[2], f_jd + f_fr)
            lats.append(plat)
            lons.append(plon)
        else:
            lats.append(np.nan)
            lons.append(np.nan)

    # Calculate previous trajectory for the last 90 minutes (every minute)
    prev_lats, prev_lons = [], []
    for minute in range(-90, 1):
        past_time = now + timedelta(minutes=minute)
        p_yr, p_mo, p_dy, p_hr, p_mi, p_se = past_time.year, past_time.month, past_time.day, past_time.hour, past_time.minute, past_time.second
        p_jd, p_fr = jday(p_yr, p_mo, p_dy, p_hr, p_mi, p_se)
        err, pos, vel = satellite.sgp4(p_jd, p_fr)
        if err == 0:
            plat, plon, _ = eci_to_geodetic(pos[0], pos[1], pos[2], p_jd + p_fr)
            prev_lats.append(plat)
            prev_lons.append(plon)
        else:
            prev_lats.append(np.nan)
            prev_lons.append(np.nan)

    # Handle longitude wrapping for plotting
    lons = np.array(lons)
    lats = np.array(lats)
    prev_lons = np.array(prev_lons)
    prev_lats = np.array(prev_lats)

    # Unwrap longitudes to handle crossing the 180/-180 meridian
    lons_unwrapped = np.unwrap(np.radians(lons))
    lons_unwrapped = np.degrees(lons_unwrapped)
    prev_lons_unwrapped = np.unwrap(np.radians(prev_lons))
    prev_lons_unwrapped = np.degrees(prev_lons_unwrapped)

    # Plot on map
    plt.figure(figsize=(16, 8))  # Wider figure for better world view
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Set map boundaries to show the full world
    ax.set_global()
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    # Set retro style with black background and green Earth
    ax.set_facecolor('black')
    plt.gcf().set_facecolor('black')

    # Create a custom colormap for Earth (green tones)
    earth_colors = ['#0a3d0a', '#1a5f1a', '#2d7a2d', '#4a9a4a', '#6bb86b']
    earth_cmap = LinearSegmentedColormap.from_list('retro_earth', earth_colors)

    # Add a simple green Earth representation
    ax.add_feature(cartopy.feature.LAND, alpha=0.8)
    ax.add_feature(cartopy.feature.OCEAN, alpha=0.9)
    ax.add_feature(cartopy.feature.COASTLINE, color='#4a9a4a', linewidth=0.8)
    ax.add_feature(cartopy.feature.BORDERS, color='#6bb86b', linewidth=0.5, alpha=0.7)

    # Add city markers
    # Plot all cities
    for city, coords in CITIES.items():
        ax.plot(coords["lon"], coords["lat"], 'o', color='#4aff4a', markersize=5, transform=ccrs.PlateCarree())
        ax.text(coords["lon"] + 1, coords["lat"], city, color='#4aff4a',
                fontsize=8, transform=ccrs.PlateCarree(), backgroundcolor='black')

    # Highlight the selected city with a larger marker
    if selected_city in CITIES:
        city_coords = CITIES[selected_city]
        ax.plot(city_coords["lon"], city_coords["lat"], 'o', color='#ffff00', markersize=8,
                transform=ccrs.PlateCarree())
        ax.text(city_coords["lon"] + 1, city_coords["lat"], selected_city, color='#ffff00',
                fontsize=10, fontweight='bold', transform=ccrs.PlateCarree(), backgroundcolor='black')

    # Set grid lines with retro styling
    ax.gridlines(color='#4a9a4a', alpha=0.3, linewidth=0.5)

    # Apply green tint to the map using a custom style
    ax.add_feature(cartopy.feature.LAND, alpha=0.8, facecolor='#2d7a2d')
    ax.add_feature(cartopy.feature.OCEAN, alpha=0.9, facecolor='#0a1a2a')

    ax.set_title(f'Satellite Ground Track (Previous 90 min + Next 90 min)\nSelected City: {selected_city}', color='#6bb86b', fontsize=14, fontweight='bold')

    # Plot previous trajectory with dotted line
    ax.plot(prev_lons_unwrapped, prev_lats, color='#ff6b35', linewidth=1.5, linestyle=':', label='Previous Path (90 min)', alpha=0.7, transform=ccrs.PlateCarree())

    # Plot current trajectory with solid line
    ax.plot(lons_unwrapped, lats, color='#ff6b35', linewidth=2, label='Future Path (90 min)', alpha=0.9, transform=ccrs.PlateCarree())

    # Plot current position
    ax.plot(lon, lat, color='#ff4757', marker='o', markersize=10, label='Current Position', alpha=0.9, transform=ccrs.PlateCarree())

    # Style the legend
    legend = ax.legend(facecolor='black', edgecolor='#4a9a4a', framealpha=0.8)
    for text in legend.get_texts():
        text.set_color('#6bb86b')

    plt.show()
else:
    print(f"SGP4 propagation error code: {error_code}")
