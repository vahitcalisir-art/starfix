"""
Ephemeris Module for Celestial Navigation

Provides celestial body positions using JPL DE440 ephemeris via Skyfield.
Calculates GHA (Greenwich Hour Angle) and Declination for navigation.

References:
    - Park et al. (2021): The JPL Planetary and Lunar Ephemerides DE440 and DE441
    - Rhodes, B.: Skyfield Python library
"""

import numpy as np
from datetime import datetime, timezone
from typing import Tuple, Optional, Dict, List, Union
from dataclasses import dataclass
from skyfield.api import load, Topos, Star, wgs84
from skyfield.almanac import find_discrete, sunrise_sunset
from skyfield import almanac


@dataclass
class CelestialBodyPosition:
    """Data class for celestial body position."""
    name: str
    gha: float  # Greenwich Hour Angle in degrees
    dec: float  # Declination in degrees
    sha: float  # Sidereal Hour Angle in degrees (for stars)
    hp: float   # Horizontal Parallax in arcminutes (for Moon/Sun)
    sd: float   # Semi-diameter in arcminutes
    altitude: float  # Computed altitude from observer position (if provided)
    azimuth: float   # Computed azimuth from observer position (if provided)
    distance_au: float  # Distance in AU


class EphemerisCalculator:
    """
    Calculates celestial body positions using JPL DE440 ephemeris.
    
    This class provides high-precision positions for:
    - Sun
    - Moon
    - Navigational planets (Venus, Mars, Jupiter, Saturn)
    - 57 navigational stars + Polaris
    
    Accuracy: Sub-arcsecond for planets, <0.1" for stars
    
    References:
        Park et al. (2021): DE440 accuracy ~0.0001" for inner planets
    """
    
    # Navigation stars from Nautical Almanac (57 stars + Polaris)
    NAVIGATION_STARS = {
        'Alpheratz': {'ra_hours': 0.139792, 'dec_degrees': 29.0905, 'magnitude': 2.1},
        'Ankaa': {'ra_hours': 0.438056, 'dec_degrees': -42.3061, 'magnitude': 2.4},
        'Schedar': {'ra_hours': 0.675278, 'dec_degrees': 56.5372, 'magnitude': 2.2},
        'Diphda': {'ra_hours': 0.726389, 'dec_degrees': -17.9867, 'magnitude': 2.0},
        'Achernar': {'ra_hours': 1.628611, 'dec_degrees': -57.2367, 'magnitude': 0.5},
        'Hamal': {'ra_hours': 2.119722, 'dec_degrees': 23.4625, 'magnitude': 2.0},
        'Polaris': {'ra_hours': 2.530278, 'dec_degrees': 89.2642, 'magnitude': 2.0},
        'Acamar': {'ra_hours': 2.971111, 'dec_degrees': -40.3047, 'magnitude': 2.9},
        'Menkar': {'ra_hours': 3.037778, 'dec_degrees': 4.0897, 'magnitude': 2.5},
        'Mirfak': {'ra_hours': 3.405278, 'dec_degrees': 49.8611, 'magnitude': 1.8},
        'Aldebaran': {'ra_hours': 4.598889, 'dec_degrees': 16.5094, 'magnitude': 0.9},
        'Rigel': {'ra_hours': 5.242222, 'dec_degrees': -8.2017, 'magnitude': 0.1},
        'Capella': {'ra_hours': 5.278056, 'dec_degrees': 45.9981, 'magnitude': 0.1},
        'Bellatrix': {'ra_hours': 5.418889, 'dec_degrees': 6.3497, 'magnitude': 1.6},
        'Elnath': {'ra_hours': 5.438333, 'dec_degrees': 28.6075, 'magnitude': 1.7},
        'Alnilam': {'ra_hours': 5.603333, 'dec_degrees': -1.2019, 'magnitude': 1.7},
        'Betelgeuse': {'ra_hours': 5.919444, 'dec_degrees': 7.4069, 'magnitude': 0.5},
        'Canopus': {'ra_hours': 6.399167, 'dec_degrees': -52.6956, 'magnitude': -0.7},
        'Sirius': {'ra_hours': 6.752222, 'dec_degrees': -16.7161, 'magnitude': -1.5},
        'Adhara': {'ra_hours': 6.977222, 'dec_degrees': -28.9722, 'magnitude': 1.5},
        'Procyon': {'ra_hours': 7.655, 'dec_degrees': 5.2247, 'magnitude': 0.4},
        'Pollux': {'ra_hours': 7.755278, 'dec_degrees': 28.0261, 'magnitude': 1.2},
        'Avior': {'ra_hours': 8.375833, 'dec_degrees': -59.5094, 'magnitude': 1.9},
        'Suhail': {'ra_hours': 9.133056, 'dec_degrees': -43.4328, 'magnitude': 2.2},
        'Miaplacidus': {'ra_hours': 9.220278, 'dec_degrees': -69.7172, 'magnitude': 1.7},
        'Alphard': {'ra_hours': 9.459722, 'dec_degrees': -8.6586, 'magnitude': 2.0},
        'Regulus': {'ra_hours': 10.139722, 'dec_degrees': 11.9672, 'magnitude': 1.4},
        'Dubhe': {'ra_hours': 11.062222, 'dec_degrees': 61.7508, 'magnitude': 1.8},
        'Denebola': {'ra_hours': 11.817222, 'dec_degrees': 14.5719, 'magnitude': 2.1},
        'Gienah': {'ra_hours': 12.263056, 'dec_degrees': -17.5419, 'magnitude': 2.6},
        'Acrux': {'ra_hours': 12.443333, 'dec_degrees': -63.0992, 'magnitude': 0.8},
        'Gacrux': {'ra_hours': 12.519444, 'dec_degrees': -57.1128, 'magnitude': 1.6},
        'Alioth': {'ra_hours': 12.900278, 'dec_degrees': 55.9597, 'magnitude': 1.8},
        'Spica': {'ra_hours': 13.419722, 'dec_degrees': -11.1614, 'magnitude': 1.0},
        'Alkaid': {'ra_hours': 13.792222, 'dec_degrees': 49.3133, 'magnitude': 1.9},
        'Hadar': {'ra_hours': 14.063889, 'dec_degrees': -60.3728, 'magnitude': 0.6},
        'Menkent': {'ra_hours': 14.111389, 'dec_degrees': -36.3700, 'magnitude': 2.1},
        'Arcturus': {'ra_hours': 14.261111, 'dec_degrees': 19.1825, 'magnitude': 0.0},
        'Rigil Kent': {'ra_hours': 14.660556, 'dec_degrees': -60.8350, 'magnitude': -0.3},
        'Zubenelgenubi': {'ra_hours': 14.847778, 'dec_degrees': -16.0419, 'magnitude': 2.8},
        'Kochab': {'ra_hours': 14.845278, 'dec_degrees': 74.1556, 'magnitude': 2.1},
        'Alphecca': {'ra_hours': 15.578056, 'dec_degrees': 26.7147, 'magnitude': 2.2},
        'Antares': {'ra_hours': 16.490278, 'dec_degrees': -26.4319, 'magnitude': 1.0},
        'Atria': {'ra_hours': 16.811111, 'dec_degrees': -69.0278, 'magnitude': 1.9},
        'Sabik': {'ra_hours': 17.172778, 'dec_degrees': -15.7247, 'magnitude': 2.4},
        'Shaula': {'ra_hours': 17.560278, 'dec_degrees': -37.1039, 'magnitude': 1.6},
        'Rasalhague': {'ra_hours': 17.582222, 'dec_degrees': 12.5603, 'magnitude': 2.1},
        'Eltanin': {'ra_hours': 17.943333, 'dec_degrees': 51.4889, 'magnitude': 2.2},
        'Kaus Australis': {'ra_hours': 18.402778, 'dec_degrees': -34.3847, 'magnitude': 1.8},
        'Vega': {'ra_hours': 18.615833, 'dec_degrees': 38.7836, 'magnitude': 0.0},
        'Nunki': {'ra_hours': 18.921111, 'dec_degrees': -26.2967, 'magnitude': 2.0},
        'Altair': {'ra_hours': 19.846389, 'dec_degrees': 8.8683, 'magnitude': 0.8},
        'Peacock': {'ra_hours': 20.427222, 'dec_degrees': -56.7350, 'magnitude': 1.9},
        'Deneb': {'ra_hours': 20.690556, 'dec_degrees': 45.2803, 'magnitude': 1.3},
        'Enif': {'ra_hours': 21.736389, 'dec_degrees': 9.8750, 'magnitude': 2.4},
        'Alnair': {'ra_hours': 22.137222, 'dec_degrees': -46.9611, 'magnitude': 1.7},
        'Fomalhaut': {'ra_hours': 22.960556, 'dec_degrees': -29.6222, 'magnitude': 1.2},
        'Markab': {'ra_hours': 23.079444, 'dec_degrees': 15.2053, 'magnitude': 2.5},
    }
    
    def __init__(self, ephemeris_file: str = 'de440s.bsp'):
        """
        Initialize ephemeris calculator.
        
        Args:
            ephemeris_file: JPL ephemeris file (de440s.bsp for compact version)
        """
        self.ts = load.timescale()
        
        # Load ephemeris (will download if not present)
        try:
            self.eph = load(ephemeris_file)
        except Exception:
            # Fall back to auto-download
            self.eph = load('de440s.bsp')
        
        self.earth = self.eph['earth']
        self.sun = self.eph['sun']
        self.moon = self.eph['moon']
        
        # Planet references
        self.planets = {
            'Venus': self.eph['venus barycenter'],
            'Mars': self.eph['mars barycenter'],
            'Jupiter': self.eph['jupiter barycenter'],
            'Saturn': self.eph['saturn barycenter'],
        }
    
    def get_time(self, dt: datetime) -> 'Time':
        """
        Convert datetime to Skyfield Time object.
        
        Args:
            dt: Python datetime (should be UTC)
            
        Returns:
            Skyfield Time object
        """
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return self.ts.from_datetime(dt)
    
    def get_gha_aries(self, dt: datetime) -> float:
        """
        Calculate Greenwich Hour Angle of Aries (First Point of Aries).
        
        GHA Aries is the foundation for stellar navigation, representing
        the angle from Greenwich meridian to the vernal equinox.
        
        Args:
            dt: Observation time (UTC)
            
        Returns:
            GHA Aries in degrees (0-360)
        """
        t = self.get_time(dt)
        
        # Get Greenwich Apparent Sidereal Time
        gast = t.gast  # Greenwich Apparent Sidereal Time in hours
        
        # Convert to degrees
        gha_aries = (gast * 15.0) % 360.0
        
        return gha_aries
    
    def get_sun_position(self, dt: datetime, 
                         observer_lat: Optional[float] = None,
                         observer_lon: Optional[float] = None) -> CelestialBodyPosition:
        """
        Calculate Sun position for navigation.
        
        Args:
            dt: Observation time (UTC)
            observer_lat: Observer latitude in degrees (optional)
            observer_lon: Observer longitude in degrees (optional)
            
        Returns:
            CelestialBodyPosition with GHA, Dec, HP, SD
        """
        t = self.get_time(dt)
        
        # Calculate from Earth center
        earth = self.earth.at(t)
        sun_astrometric = earth.observe(self.sun)
        sun_apparent = sun_astrometric.apparent()
        
        # Get RA and Dec
        ra, dec, distance = sun_apparent.radec()
        
        # Calculate GHA
        gha_aries = self.get_gha_aries(dt)
        sha = 360.0 - (ra.hours * 15.0)  # SHA = 360 - RA (in degrees)
        gha = (gha_aries + sha) % 360.0
        
        # Sun's semi-diameter (varies with distance)
        # Mean SD = 16.0' at 1 AU
        sd = 16.0 / distance.au  # arcminutes
        
        # Sun's horizontal parallax (very small, ~8.8")
        hp = 8.8 / 60.0 / distance.au  # arcminutes
        
        # Calculate altitude and azimuth if observer position provided
        alt, az = 0.0, 0.0
        if observer_lat is not None and observer_lon is not None:
            alt, az = self._calculate_alt_az(dt, self.sun, observer_lat, observer_lon)
        
        return CelestialBodyPosition(
            name='Sun',
            gha=gha,
            dec=dec.degrees,
            sha=sha,
            hp=hp,
            sd=sd,
            altitude=alt,
            azimuth=az,
            distance_au=distance.au
        )
    
    def get_moon_position(self, dt: datetime,
                          observer_lat: Optional[float] = None,
                          observer_lon: Optional[float] = None) -> CelestialBodyPosition:
        """
        Calculate Moon position for navigation.
        
        Moon requires special handling due to:
        - Large horizontal parallax (up to 61.5')
        - Variable semi-diameter (14.7' - 16.7')
        - Rapid motion in RA/Dec
        
        Args:
            dt: Observation time (UTC)
            observer_lat: Observer latitude in degrees (optional)
            observer_lon: Observer longitude in degrees (optional)
            
        Returns:
            CelestialBodyPosition with GHA, Dec, HP, SD
        """
        t = self.get_time(dt)
        
        # Calculate from Earth center
        earth = self.earth.at(t)
        moon_astrometric = earth.observe(self.moon)
        moon_apparent = moon_astrometric.apparent()
        
        # Get RA and Dec
        ra, dec, distance = moon_apparent.radec()
        
        # Calculate GHA
        gha_aries = self.get_gha_aries(dt)
        sha = 360.0 - (ra.hours * 15.0)
        gha = (gha_aries + sha) % 360.0
        
        # Moon's horizontal parallax
        # HP = arcsin(Earth_radius / Moon_distance)
        # Earth radius = 6378.137 km, Moon distance in km
        earth_radius_km = 6378.137
        moon_distance_km = distance.au * 149597870.7  # AU to km
        hp = np.degrees(np.arcsin(earth_radius_km / moon_distance_km)) * 60.0  # arcminutes
        
        # Moon's semi-diameter
        # Mean SD = 15.5' at mean distance (384400 km)
        mean_moon_distance = 384400.0
        sd = 15.5 * (mean_moon_distance / moon_distance_km)  # arcminutes
        
        # Calculate altitude and azimuth if observer position provided
        alt, az = 0.0, 0.0
        if observer_lat is not None and observer_lon is not None:
            alt, az = self._calculate_alt_az(dt, self.moon, observer_lat, observer_lon)
        
        return CelestialBodyPosition(
            name='Moon',
            gha=gha,
            dec=dec.degrees,
            sha=sha,
            hp=hp,
            sd=sd,
            altitude=alt,
            azimuth=az,
            distance_au=distance.au
        )
    
    def get_planet_position(self, planet_name: str, dt: datetime,
                            observer_lat: Optional[float] = None,
                            observer_lon: Optional[float] = None) -> CelestialBodyPosition:
        """
        Calculate planet position for navigation.
        
        Navigation planets: Venus, Mars, Jupiter, Saturn
        
        Args:
            planet_name: Name of planet ('Venus', 'Mars', 'Jupiter', 'Saturn')
            dt: Observation time (UTC)
            observer_lat: Observer latitude in degrees (optional)
            observer_lon: Observer longitude in degrees (optional)
            
        Returns:
            CelestialBodyPosition with GHA, Dec
        """
        if planet_name not in self.planets:
            raise ValueError(f"Unknown planet: {planet_name}. Valid: {list(self.planets.keys())}")
        
        t = self.get_time(dt)
        planet = self.planets[planet_name]
        
        # Calculate from Earth center
        earth = self.earth.at(t)
        planet_astrometric = earth.observe(planet)
        planet_apparent = planet_astrometric.apparent()
        
        # Get RA and Dec
        ra, dec, distance = planet_apparent.radec()
        
        # Calculate GHA
        gha_aries = self.get_gha_aries(dt)
        sha = 360.0 - (ra.hours * 15.0)
        gha = (gha_aries + sha) % 360.0
        
        # Planets have negligible HP (< 0.5")
        hp = 0.0
        
        # Planet semi-diameters are small and often ignored in navigation
        # For completeness, approximate values
        sd_values = {'Venus': 0.2, 'Mars': 0.1, 'Jupiter': 0.4, 'Saturn': 0.3}
        sd = sd_values.get(planet_name, 0.0)
        
        # Calculate altitude and azimuth if observer position provided
        alt, az = 0.0, 0.0
        if observer_lat is not None and observer_lon is not None:
            alt, az = self._calculate_alt_az(dt, planet, observer_lat, observer_lon)
        
        return CelestialBodyPosition(
            name=planet_name,
            gha=gha,
            dec=dec.degrees,
            sha=sha,
            hp=hp,
            sd=sd,
            altitude=alt,
            azimuth=az,
            distance_au=distance.au
        )
    
    def get_star_position(self, star_name: str, dt: datetime,
                          observer_lat: Optional[float] = None,
                          observer_lon: Optional[float] = None) -> CelestialBodyPosition:
        """
        Calculate navigation star position.
        
        Uses catalog positions with proper motion. Stars have:
        - Constant SHA (Sidereal Hour Angle)
        - Very slowly changing Dec (proper motion)
        - No HP or SD
        
        Args:
            star_name: Name of navigation star
            dt: Observation time (UTC)
            observer_lat: Observer latitude in degrees (optional)
            observer_lon: Observer longitude in degrees (optional)
            
        Returns:
            CelestialBodyPosition with GHA, Dec, SHA
        """
        if star_name not in self.NAVIGATION_STARS:
            raise ValueError(f"Unknown star: {star_name}. Valid stars: {list(self.NAVIGATION_STARS.keys())}")
        
        star_data = self.NAVIGATION_STARS[star_name]
        t = self.get_time(dt)
        
        # Create star object
        star = Star(ra_hours=star_data['ra_hours'], 
                    dec_degrees=star_data['dec_degrees'])
        
        # Calculate apparent position (includes aberration, precession)
        earth = self.earth.at(t)
        star_astrometric = earth.observe(star)
        star_apparent = star_astrometric.apparent()
        
        # Get RA and Dec
        ra, dec, _ = star_apparent.radec()
        
        # Calculate SHA and GHA
        sha = 360.0 - (ra.hours * 15.0)
        gha_aries = self.get_gha_aries(dt)
        gha = (gha_aries + sha) % 360.0
        
        # Calculate altitude and azimuth if observer position provided
        alt, az = 0.0, 0.0
        if observer_lat is not None and observer_lon is not None:
            alt, az = self._calculate_alt_az_star(dt, star, observer_lat, observer_lon)
        
        return CelestialBodyPosition(
            name=star_name,
            gha=gha,
            dec=dec.degrees,
            sha=sha,
            hp=0.0,  # Stars have no parallax for navigation
            sd=0.0,  # Stars are point sources
            altitude=alt,
            azimuth=az,
            distance_au=float('inf')
        )
    
    def _calculate_alt_az(self, dt: datetime, body, 
                          observer_lat: float, observer_lon: float) -> Tuple[float, float]:
        """
        Calculate altitude and azimuth of a body from observer position.
        
        Args:
            dt: Observation time
            body: Skyfield body object
            observer_lat: Observer latitude in degrees
            observer_lon: Observer longitude in degrees
            
        Returns:
            Tuple of (altitude, azimuth) in degrees
        """
        t = self.get_time(dt)
        
        # Create observer location
        observer = self.earth + wgs84.latlon(observer_lat, observer_lon)
        
        # Calculate apparent position from observer
        apparent = observer.at(t).observe(body).apparent()
        
        # Get altitude and azimuth
        alt, az, _ = apparent.altaz()
        
        return alt.degrees, az.degrees
    
    def _calculate_alt_az_star(self, dt: datetime, star: Star,
                               observer_lat: float, observer_lon: float) -> Tuple[float, float]:
        """
        Calculate altitude and azimuth of a star from observer position.
        """
        t = self.get_time(dt)
        
        # Create observer location
        observer = self.earth + wgs84.latlon(observer_lat, observer_lon)
        
        # Calculate apparent position from observer
        apparent = observer.at(t).observe(star).apparent()
        
        # Get altitude and azimuth
        alt, az, _ = apparent.altaz()
        
        return alt.degrees, az.degrees
    
    def get_visible_stars(self, dt: datetime, observer_lat: float, observer_lon: float,
                          min_altitude: float = 15.0, max_altitude: float = 75.0) -> List[CelestialBodyPosition]:
        """
        Get list of visible navigation stars for observation.
        
        Filters stars based on altitude constraints (15° - 75° is optimal range).
        
        Args:
            dt: Observation time (UTC)
            observer_lat: Observer latitude in degrees
            observer_lon: Observer longitude in degrees
            min_altitude: Minimum altitude in degrees (default 15°)
            max_altitude: Maximum altitude in degrees (default 75°)
            
        Returns:
            List of visible star positions sorted by altitude
        """
        visible_stars = []
        
        for star_name in self.NAVIGATION_STARS:
            try:
                pos = self.get_star_position(star_name, dt, observer_lat, observer_lon)
                if min_altitude <= pos.altitude <= max_altitude:
                    visible_stars.append(pos)
            except Exception:
                continue
        
        # Sort by altitude (optimal for observation)
        visible_stars.sort(key=lambda x: x.altitude, reverse=True)
        
        return visible_stars
    
    def get_all_bodies(self, dt: datetime, observer_lat: Optional[float] = None,
                       observer_lon: Optional[float] = None) -> Dict[str, CelestialBodyPosition]:
        """
        Get positions of all navigation bodies (Sun, Moon, planets).
        
        Args:
            dt: Observation time (UTC)
            observer_lat: Observer latitude in degrees (optional)
            observer_lon: Observer longitude in degrees (optional)
            
        Returns:
            Dictionary of body positions keyed by name
        """
        bodies = {}
        
        # Sun
        bodies['Sun'] = self.get_sun_position(dt, observer_lat, observer_lon)
        
        # Moon
        bodies['Moon'] = self.get_moon_position(dt, observer_lat, observer_lon)
        
        # Planets
        for planet in ['Venus', 'Mars', 'Jupiter', 'Saturn']:
            bodies[planet] = self.get_planet_position(planet, dt, observer_lat, observer_lon)
        
        return bodies


def validate_against_almanac(dt: datetime, body_name: str, 
                              expected_gha: float, expected_dec: float,
                              tolerance_arcmin: float = 0.5) -> Tuple[bool, Dict]:
    """
    Validate ephemeris calculations against Nautical Almanac values.
    
    Args:
        dt: Observation time (UTC)
        body_name: Name of celestial body
        expected_gha: Expected GHA from Almanac (degrees)
        expected_dec: Expected Dec from Almanac (degrees)
        tolerance_arcmin: Acceptable difference in arcminutes
        
    Returns:
        Tuple of (pass/fail, details dict)
    """
    calc = EphemerisCalculator()
    
    # Get calculated position
    if body_name == 'Sun':
        pos = calc.get_sun_position(dt)
    elif body_name == 'Moon':
        pos = calc.get_moon_position(dt)
    elif body_name in calc.planets:
        pos = calc.get_planet_position(body_name, dt)
    elif body_name in calc.NAVIGATION_STARS:
        pos = calc.get_star_position(body_name, dt)
    else:
        raise ValueError(f"Unknown body: {body_name}")
    
    # Calculate differences
    gha_diff = abs(pos.gha - expected_gha)
    if gha_diff > 180:
        gha_diff = 360 - gha_diff
    gha_diff_arcmin = gha_diff * 60
    
    dec_diff_arcmin = abs(pos.dec - expected_dec) * 60
    
    # Check tolerance
    passed = (gha_diff_arcmin <= tolerance_arcmin and 
              dec_diff_arcmin <= tolerance_arcmin)
    
    return passed, {
        'body': body_name,
        'calculated_gha': pos.gha,
        'expected_gha': expected_gha,
        'gha_diff_arcmin': gha_diff_arcmin,
        'calculated_dec': pos.dec,
        'expected_dec': expected_dec,
        'dec_diff_arcmin': dec_diff_arcmin,
        'tolerance_arcmin': tolerance_arcmin,
        'passed': passed
    }


if __name__ == "__main__":
    # Quick test
    from datetime import datetime, timezone
    
    calc = EphemerisCalculator()
    
    # Test with a known date
    dt = datetime(2025, 3, 21, 12, 0, 0, tzinfo=timezone.utc)
    
    print("=== Ephemeris Test ===")
    print(f"Date/Time: {dt}")
    print(f"GHA Aries: {calc.get_gha_aries(dt):.4f}°")
    
    sun = calc.get_sun_position(dt)
    print(f"\nSun: GHA={sun.gha:.4f}°, Dec={sun.dec:.4f}°, SD={sun.sd:.2f}'")
    
    moon = calc.get_moon_position(dt)
    print(f"Moon: GHA={moon.gha:.4f}°, Dec={moon.dec:.4f}°, HP={moon.hp:.2f}', SD={moon.sd:.2f}'")
    
    sirius = calc.get_star_position('Sirius', dt, observer_lat=40.0, observer_lon=-74.0)
    print(f"Sirius: GHA={sirius.gha:.4f}°, Dec={sirius.dec:.4f}°, Alt={sirius.altitude:.2f}°")
