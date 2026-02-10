"""
Sight Reduction Module for Celestial Navigation

Core calculations for converting sextant observations to navigation data.
Includes altitude corrections and the fundamental navigation triangle solution.

References:
    - Bowditch (2019): American Practical Navigator
    - Karl (2011): Celestial Navigation in the GPS Age
    - Hohenkerk et al. (2012): Astro Navigation Remembered
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional
from enum import Enum


class BodyType(Enum):
    """Type of celestial body for correction purposes."""
    SUN_LOWER = "sun_lower"
    SUN_UPPER = "sun_upper"
    MOON_LOWER = "moon_lower"
    MOON_UPPER = "moon_upper"
    STAR = "star"
    PLANET = "planet"


@dataclass
class AltitudeCorrections:
    """Container for all altitude corrections applied."""
    dip: float           # Dip correction (arcminutes)
    refraction: float    # Atmospheric refraction (arcminutes)
    parallax: float      # Parallax in altitude (arcminutes)
    semidiameter: float  # Semi-diameter correction (arcminutes)
    total: float         # Total correction (arcminutes)
    

@dataclass
class SightReductionResult:
    """Result of sight reduction calculation."""
    Hc: float           # Computed altitude (degrees)
    Zn: float           # Azimuth (degrees, 0-360 from North)
    altitude_intercept: float  # a = Ho - Hc (nautical miles)
    toward_away: str    # 'toward' or 'away' from GP
    

def calculate_dip(height_of_eye_m: float) -> float:
    """
    Calculate dip correction for observer height above sea level.
    
    The dip (height of eye) correction accounts for the observer
    seeing a horizon that is below the true horizontal.
    
    Formula: Dip = 1.76 * sqrt(height_m)  [in arcminutes]
    
    Reference: Bowditch (2019), Table A2
    
    Args:
        height_of_eye_m: Height of eye above sea level in meters
        
    Returns:
        Dip correction in arcminutes (always negative)
    """
    if height_of_eye_m < 0:
        raise ValueError("Height of eye cannot be negative")
    
    # Standard formula: Dip = 1.76 * sqrt(h) arcminutes
    # where h is in meters
    dip = 1.76 * np.sqrt(height_of_eye_m)
    
    return -dip  # Always negative (horizon is depressed)


def calculate_refraction(apparent_altitude_deg: float, 
                         temperature_c: float = 10.0,
                         pressure_mb: float = 1010.0) -> float:
    """
    Calculate atmospheric refraction correction.
    
    Refraction causes celestial bodies to appear higher than they actually are.
    This function uses the Bennett formula with meteorological corrections.
    
    Standard conditions: T=10°C, P=1010 mb
    
    Reference: 
        - Bennett (1982): Calculation of Astronomical Refraction
        - Meeus (1998): Astronomical Algorithms
    
    Args:
        apparent_altitude_deg: Apparent (observed) altitude in degrees
        temperature_c: Air temperature in Celsius (default 10°C)
        pressure_mb: Atmospheric pressure in millibars (default 1010 mb)
        
    Returns:
        Refraction correction in arcminutes (always negative for correction)
    """
    if apparent_altitude_deg < -1:
        return 0.0  # Below horizon, no meaningful refraction
    
    # Ensure altitude is at least 0 for formula stability
    h = max(apparent_altitude_deg, 0.0)
    
    # Bennett formula for refraction in arcminutes
    # R = 1 / tan(h + 7.31/(h + 4.4)) [in arcminutes]
    # This is accurate to about 0.07' for h > 5°
    
    h_rad = np.radians(h)
    
    if h > 0:
        cot_term = h + 7.31 / (h + 4.4)
        R0 = 1.0 / np.tan(np.radians(cot_term))
    else:
        R0 = 34.5  # Refraction at horizon (approximately)
    
    # Apply meteorological corrections
    # Correction factor = (P/1010) * (283/(273+T))
    correction_factor = (pressure_mb / 1010.0) * (283.0 / (273.0 + temperature_c))
    
    refraction = R0 * correction_factor
    
    return -refraction  # Negative because we subtract from observed altitude


def calculate_parallax(altitude_deg: float, horizontal_parallax_arcmin: float) -> float:
    """
    Calculate parallax in altitude correction.
    
    Parallax is the difference between the topocentric (observer) and
    geocentric (Earth center) positions. Significant for Moon, minor for Sun.
    
    Formula: P = HP * cos(altitude)
    
    Reference: Bowditch (2019)
    
    Args:
        altitude_deg: Altitude of body in degrees
        horizontal_parallax_arcmin: Horizontal parallax in arcminutes
            - Moon: ~57-61'
            - Sun: ~0.15'
            - Stars/Planets: negligible
            
    Returns:
        Parallax correction in arcminutes (always positive)
    """
    # Parallax in altitude = HP * cos(altitude)
    parallax = horizontal_parallax_arcmin * np.cos(np.radians(altitude_deg))
    
    return parallax  # Positive (body appears lower than actual)


def calculate_semidiameter_correction(sd_arcmin: float, limb: str = 'lower') -> float:
    """
    Calculate semi-diameter correction for Sun and Moon observations.
    
    When observing the limb (edge) of the Sun or Moon, a correction
    is needed to find the center altitude.
    
    Args:
        sd_arcmin: Semi-diameter in arcminutes
            - Sun: ~15.9-16.3'
            - Moon: ~14.7-16.7'
        limb: 'lower' or 'upper' limb observed
        
    Returns:
        Semi-diameter correction in arcminutes
    """
    if limb.lower() == 'lower':
        return sd_arcmin  # Add SD for lower limb
    elif limb.lower() == 'upper':
        return -sd_arcmin  # Subtract SD for upper limb
    else:
        raise ValueError("Limb must be 'lower' or 'upper'")


def apply_altitude_corrections(
    sextant_altitude_deg: float,
    height_of_eye_m: float,
    body_type: BodyType,
    horizontal_parallax_arcmin: float = 0.0,
    semidiameter_arcmin: float = 0.0,
    index_error_arcmin: float = 0.0,
    temperature_c: float = 10.0,
    pressure_mb: float = 1010.0
) -> Tuple[float, AltitudeCorrections]:
    """
    Apply all altitude corrections to convert sextant altitude to observed altitude.
    
    Correction sequence:
    1. Index error (instrument)
    2. Dip (height of eye)
    3. Refraction (atmosphere)
    4. Semi-diameter (Sun/Moon limb)
    5. Parallax (topocentric correction)
    
    The result is Ho (Observed Altitude), ready for sight reduction.
    
    Reference: Bowditch (2019), Karl (2011)
    
    Args:
        sextant_altitude_deg: Raw sextant reading in degrees
        height_of_eye_m: Observer height above sea level in meters
        body_type: Type of celestial body observed
        horizontal_parallax_arcmin: HP for Moon/Sun (arcminutes)
        semidiameter_arcmin: Semi-diameter for Sun/Moon (arcminutes)
        index_error_arcmin: Sextant index error (arcminutes)
        temperature_c: Air temperature (Celsius)
        pressure_mb: Atmospheric pressure (millibars)
        
    Returns:
        Tuple of (Ho in degrees, AltitudeCorrections details)
    """
    # Start with sextant altitude in arcminutes for precision
    hs_arcmin = sextant_altitude_deg * 60.0
    
    # 1. Apply index error
    ha_arcmin = hs_arcmin - index_error_arcmin  # Apparent altitude
    
    # 2. Calculate dip
    dip = calculate_dip(height_of_eye_m)
    ha_arcmin += dip
    
    # 3. Calculate refraction (using current apparent altitude)
    apparent_alt_deg = ha_arcmin / 60.0
    refraction = calculate_refraction(apparent_alt_deg, temperature_c, pressure_mb)
    ha_arcmin += refraction
    
    # 4. Semi-diameter correction (for Sun and Moon)
    sd_correction = 0.0
    if body_type in [BodyType.SUN_LOWER, BodyType.MOON_LOWER]:
        sd_correction = calculate_semidiameter_correction(semidiameter_arcmin, 'lower')
    elif body_type in [BodyType.SUN_UPPER, BodyType.MOON_UPPER]:
        sd_correction = calculate_semidiameter_correction(semidiameter_arcmin, 'upper')
    ha_arcmin += sd_correction
    
    # 5. Parallax correction (mainly for Moon)
    parallax = 0.0
    if horizontal_parallax_arcmin > 0:
        current_alt = ha_arcmin / 60.0
        parallax = calculate_parallax(current_alt, horizontal_parallax_arcmin)
        ha_arcmin += parallax
    
    # Total correction
    total_correction = dip + refraction + sd_correction + parallax - index_error_arcmin
    
    # Convert to degrees
    Ho = ha_arcmin / 60.0
    
    corrections = AltitudeCorrections(
        dip=dip,
        refraction=refraction,
        parallax=parallax,
        semidiameter=sd_correction,
        total=total_correction
    )
    
    return Ho, corrections


def compute_altitude_azimuth(
    observer_lat_deg: float,
    observer_lon_deg: float,
    gha_deg: float,
    dec_deg: float
) -> Tuple[float, float]:
    """
    Compute computed altitude (Hc) and azimuth (Zn) using navigation triangle.
    
    This is the core sight reduction calculation, solving the spherical
    navigation triangle for altitude and azimuth.
    
    Navigation Triangle:
        - Vertex P: Elevated pole (North or South)
        - Vertex Z: Zenith (observer's position)
        - Vertex X: Celestial body (GP - geographical position)
    
    Formulas:
        sin(Hc) = sin(Lat) * sin(Dec) + cos(Lat) * cos(Dec) * cos(LHA)
        
        tan(Z) = sin(LHA) / (cos(Lat) * tan(Dec) - sin(Lat) * cos(LHA))
        
    Reference: 
        - Bowditch (2019): Chapter 20
        - Karl (2011): Celestial Navigation in the GPS Age
    
    Args:
        observer_lat_deg: Assumed/DR latitude in degrees (+ North, - South)
        observer_lon_deg: Assumed/DR longitude in degrees (+ East, - West)
        gha_deg: Greenwich Hour Angle of body in degrees
        dec_deg: Declination of body in degrees (+ North, - South)
        
    Returns:
        Tuple of (Hc in degrees, Zn in degrees 0-360 from North)
    """
    # Convert to radians
    lat = np.radians(observer_lat_deg)
    dec = np.radians(dec_deg)
    
    # Calculate Local Hour Angle
    # LHA = GHA + Longitude (East positive)
    lha_deg = (gha_deg + observer_lon_deg) % 360.0
    lha = np.radians(lha_deg)
    
    # Calculate computed altitude (Hc)
    sin_Hc = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(lha)
    
    # Clamp to valid range for arcsin
    sin_Hc = np.clip(sin_Hc, -1.0, 1.0)
    Hc = np.degrees(np.arcsin(sin_Hc))
    
    # Calculate azimuth angle (Z)
    # Using: tan(Z) = sin(LHA) / (cos(Lat) * tan(Dec) - sin(Lat) * cos(LHA))
    # Or: cos(Z) = (sin(Dec) - sin(Lat) * sin(Hc)) / (cos(Lat) * cos(Hc))
    
    cos_Hc = np.cos(np.radians(Hc))
    
    if abs(cos_Hc) < 1e-10:
        # Body at zenith or nadir - azimuth undefined
        Zn = 0.0
    else:
        # Calculate azimuth using sine and cosine for quadrant determination
        sin_Z = np.sin(lha) * np.cos(dec) / cos_Hc
        cos_Z = (np.sin(dec) - np.sin(lat) * np.sin(np.radians(Hc))) / (np.cos(lat) * cos_Hc)
        
        # Clamp values
        sin_Z = np.clip(sin_Z, -1.0, 1.0)
        cos_Z = np.clip(cos_Z, -1.0, 1.0)
        
        # Use atan2 for proper quadrant
        Z = np.degrees(np.arctan2(sin_Z, cos_Z))
        
        # Convert to true azimuth (0-360 from North)
        if observer_lat_deg >= 0:
            # Northern hemisphere
            if Z >= 0:
                Zn = Z
            else:
                Zn = 360.0 + Z
        else:
            # Southern hemisphere
            Zn = (180.0 - Z) % 360.0
    
    return Hc, Zn


def sight_reduction(
    observer_lat_deg: float,
    observer_lon_deg: float,
    gha_deg: float,
    dec_deg: float,
    observed_altitude_deg: float
) -> SightReductionResult:
    """
    Perform complete sight reduction to get intercept and azimuth.
    
    This function combines altitude/azimuth computation with intercept
    calculation, producing all data needed for plotting a line of position.
    
    The intercept (a) is the difference between observed and computed altitudes:
        a = Ho - Hc
        
    A positive intercept means the observer is closer to the GP (toward).
    A negative intercept means the observer is farther from the GP (away).
    
    Each arcminute of intercept equals one nautical mile.
    
    Reference: Marcq Saint-Hilaire intercept method (1875)
    
    Args:
        observer_lat_deg: Assumed position latitude (degrees)
        observer_lon_deg: Assumed position longitude (degrees)
        gha_deg: Greenwich Hour Angle (degrees)
        dec_deg: Declination (degrees)
        observed_altitude_deg: Observed altitude Ho after corrections (degrees)
        
    Returns:
        SightReductionResult containing Hc, Zn, intercept, toward/away
    """
    # Compute Hc and Zn
    Hc, Zn = compute_altitude_azimuth(observer_lat_deg, observer_lon_deg, gha_deg, dec_deg)
    
    # Calculate altitude intercept in arcminutes (= nautical miles)
    intercept_arcmin = (observed_altitude_deg - Hc) * 60.0
    
    # Determine toward or away
    if intercept_arcmin >= 0:
        toward_away = "toward"
    else:
        toward_away = "away"
    
    return SightReductionResult(
        Hc=Hc,
        Zn=Zn,
        altitude_intercept=intercept_arcmin,
        toward_away=toward_away
    )


def calculate_line_of_position(
    assumed_lat_deg: float,
    assumed_lon_deg: float,
    azimuth_deg: float,
    intercept_nm: float
) -> Tuple[float, float, float]:
    """
    Calculate the Line of Position (LOP) from intercept and azimuth.
    
    The LOP is perpendicular to the azimuth, passing through a point
    that is 'intercept' nautical miles from the assumed position
    in the direction of (or opposite to) the azimuth.
    
    Args:
        assumed_lat_deg: Assumed position latitude
        assumed_lon_deg: Assumed position longitude
        azimuth_deg: Azimuth to celestial body GP
        intercept_nm: Altitude intercept in nautical miles
        
    Returns:
        Tuple of (LOP_lat, LOP_lon, LOP_bearing)
        where LOP passes through (LOP_lat, LOP_lon) with bearing LOP_bearing
    """
    # Calculate point on LOP
    # Move from AP in direction of azimuth by intercept distance
    
    azimuth_rad = np.radians(azimuth_deg)
    lat_rad = np.radians(assumed_lat_deg)
    
    # Distance in degrees (1 nm = 1/60 degree)
    d = intercept_nm / 60.0
    
    # Calculate new position
    delta_lat = d * np.cos(azimuth_rad)
    delta_lon = d * np.sin(azimuth_rad) / np.cos(lat_rad)
    
    lop_lat = assumed_lat_deg + delta_lat
    lop_lon = assumed_lon_deg + delta_lon
    
    # LOP bearing is perpendicular to azimuth
    lop_bearing = (azimuth_deg + 90.0) % 360.0
    
    return lop_lat, lop_lon, lop_bearing


def validate_sight_reduction(
    test_lat: float = 34.0,
    test_lon: float = -120.0,
    test_gha: float = 45.0,
    test_dec: float = 20.0,
    expected_Hc: Optional[float] = None,
    expected_Zn: Optional[float] = None
) -> dict:
    """
    Validate sight reduction calculation with known values.
    
    This can be used to verify the implementation against
    published tables (H.O. 229) or other references.
    
    Args:
        test_lat: Test latitude
        test_lon: Test longitude  
        test_gha: Test GHA
        test_dec: Test declination
        expected_Hc: Expected computed altitude (optional)
        expected_Zn: Expected azimuth (optional)
        
    Returns:
        Dictionary with calculated values and differences
    """
    Hc, Zn = compute_altitude_azimuth(test_lat, test_lon, test_gha, test_dec)
    
    result = {
        'test_lat': test_lat,
        'test_lon': test_lon,
        'test_gha': test_gha,
        'test_dec': test_dec,
        'LHA': (test_gha + test_lon) % 360.0,
        'computed_Hc': Hc,
        'computed_Zn': Zn
    }
    
    if expected_Hc is not None:
        result['expected_Hc'] = expected_Hc
        result['Hc_diff_arcmin'] = (Hc - expected_Hc) * 60
        
    if expected_Zn is not None:
        result['expected_Zn'] = expected_Zn
        result['Zn_diff_deg'] = abs(Zn - expected_Zn)
    
    return result


if __name__ == "__main__":
    # Quick validation test
    print("=== Sight Reduction Module Test ===\n")
    
    # Test altitude corrections
    print("1. Altitude Corrections Test:")
    print("-" * 40)
    
    Ho, corr = apply_altitude_corrections(
        sextant_altitude_deg=35.5,
        height_of_eye_m=3.0,
        body_type=BodyType.SUN_LOWER,
        horizontal_parallax_arcmin=0.15,
        semidiameter_arcmin=16.0,
        index_error_arcmin=0.0
    )
    
    print(f"  Sextant Altitude: 35° 30.0'")
    print(f"  Height of Eye: 3.0 m")
    print(f"  Body: Sun (Lower Limb)")
    print(f"  Dip: {corr.dip:.2f}'")
    print(f"  Refraction: {corr.refraction:.2f}'")
    print(f"  Semi-diameter: +{corr.semidiameter:.2f}'")
    print(f"  Parallax: +{corr.parallax:.2f}'")
    print(f"  Total Correction: {corr.total:.2f}'")
    print(f"  Observed Altitude (Ho): {Ho:.4f}°")
    
    # Test sight reduction
    print("\n2. Sight Reduction Test:")
    print("-" * 40)
    
    result = sight_reduction(
        observer_lat_deg=34.0,
        observer_lon_deg=-120.0,
        gha_deg=315.0,
        dec_deg=20.0,
        observed_altitude_deg=45.5
    )
    
    print(f"  Assumed Position: 34°N, 120°W")
    print(f"  GHA: 315°, Dec: N20°")
    print(f"  LHA: {(315.0 - 120.0) % 360}°")
    print(f"  Computed Altitude (Hc): {result.Hc:.4f}°")
    print(f"  Azimuth (Zn): {result.Zn:.1f}°")
    print(f"  Intercept: {abs(result.altitude_intercept):.1f}' {result.toward_away}")
    
    # Validate against known values
    print("\n3. Validation Test:")
    print("-" * 40)
    validation = validate_sight_reduction(
        test_lat=40.0,
        test_lon=-74.0,
        test_gha=45.0,
        test_dec=23.0
    )
    print(f"  Test Position: 40°N, 74°W")
    print(f"  GHA: 45°, Dec: N23°")
    print(f"  LHA: {validation['LHA']}°")
    print(f"  Hc: {validation['computed_Hc']:.4f}°")
    print(f"  Zn: {validation['computed_Zn']:.1f}°")
