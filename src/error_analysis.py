"""
Error Analysis Module for Celestial Navigation

Implements error analysis and uncertainty quantification:
1. HDOP/VDOP/GDOP calculations
2. Confidence ellipse computation
3. Error propagation
4. Monte Carlo simulation support

References:
    - Hoover (1984): Algorithms for Confidence Circles and Ellipses
    - Swaszek et al. (2019): Rethinking Star Selection
    - Kaplan (1995): Error ellipse methodology
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional
from scipy import stats
from scipy.linalg import svd, inv


@dataclass
class ConfidenceEllipse:
    """Parameters defining a confidence ellipse for position uncertainty."""
    center_lat: float
    center_lon: float
    semi_major_nm: float     # Semi-major axis in nautical miles
    semi_minor_nm: float     # Semi-minor axis in nautical miles
    orientation_deg: float   # Orientation of major axis (degrees from North)
    confidence_level: float  # 0.95 for 95% confidence
    area_sq_nm: float        # Area of ellipse in square nautical miles


@dataclass
class ErrorMetrics:
    """Complete error metrics for a position fix."""
    hdop: float              # Horizontal Dilution of Precision
    vdop: float              # Vertical DOP (for 3D solutions)
    pdop: float              # Position DOP
    gdop: float              # Geometric DOP
    sigma_lat_nm: float      # Standard error in latitude (nm)
    sigma_lon_nm: float      # Standard error in longitude (nm)
    cep_nm: float            # Circular Error Probable (50%)
    r95_nm: float            # 95% radius
    drms_nm: float          # Distance RMS (67%)


@dataclass
class Observation:
    """Observation data for error calculations."""
    azimuth_deg: float      # Azimuth to celestial body
    altitude_error_arcmin: float = 1.0  # Assumed observation error


def calculate_hdop_from_azimuths(azimuths_deg: List[float]) -> float:
    """
    Calculate HDOP from observation azimuths.
    
    HDOP indicates how observation geometry affects position accuracy.
    Lower HDOP = better geometry = more accurate fix.
    
    Optimal geometry: observations evenly distributed in azimuth
    Poor geometry: observations clustered in similar azimuths
    
    Reference: Swaszek et al. (2019)
    
    Args:
        azimuths_deg: List of azimuths to celestial bodies (degrees)
        
    Returns:
        HDOP value (dimensionless, typical range 1-10)
    """
    if len(azimuths_deg) < 2:
        return float('inf')
    
    n = len(azimuths_deg)
    
    # Build geometry matrix G
    # Each row represents the unit vector toward a celestial body
    G = np.zeros((n, 2))
    
    for i, az in enumerate(azimuths_deg):
        az_rad = np.radians(az)
        G[i, 0] = np.cos(az_rad)  # North component
        G[i, 1] = np.sin(az_rad)  # East component
    
    try:
        # DOP matrix = (G^T G)^-1
        GtG = G.T @ G
        GtG_inv = inv(GtG)
        
        # HDOP = sqrt(σ_N² + σ_E²)
        hdop = np.sqrt(GtG_inv[0, 0] + GtG_inv[1, 1])
        
        return hdop
    
    except np.linalg.LinAlgError:
        return float('inf')


def calculate_dop_full(azimuths_deg: List[float], altitudes_deg: Optional[List[float]] = None) -> ErrorMetrics:
    """
    Calculate complete DOP values (HDOP, VDOP, PDOP, GDOP).
    
    For 2D celestial navigation, HDOP is most relevant.
    VDOP/PDOP/GDOP are included for completeness and 3D applications.
    
    Args:
        azimuths_deg: Azimuths to observed bodies
        altitudes_deg: Altitudes of observed bodies (optional, for 3D DOP)
        
    Returns:
        ErrorMetrics with all DOP values
    """
    n = len(azimuths_deg)
    
    if n < 2:
        return ErrorMetrics(
            hdop=float('inf'), vdop=float('inf'), pdop=float('inf'), 
            gdop=float('inf'), sigma_lat_nm=float('inf'), sigma_lon_nm=float('inf'),
            cep_nm=float('inf'), r95_nm=float('inf'), drms_nm=float('inf')
        )
    
    # Build 2D geometry matrix
    G = np.zeros((n, 2))
    for i, az in enumerate(azimuths_deg):
        az_rad = np.radians(az)
        G[i, 0] = np.cos(az_rad)
        G[i, 1] = np.sin(az_rad)
    
    try:
        GtG = G.T @ G
        GtG_inv = inv(GtG)
        
        hdop = np.sqrt(GtG_inv[0, 0] + GtG_inv[1, 1])
        
        # For celestial navigation, assume 1' observation error
        sigma_obs_nm = 1.0  # 1 arcminute = 1 nm
        
        sigma_lat = np.sqrt(GtG_inv[0, 0]) * sigma_obs_nm
        sigma_lon = np.sqrt(GtG_inv[1, 1]) * sigma_obs_nm
        
        # CEP (Circular Error Probable) - 50% confidence
        # For bivariate normal: CEP ≈ 0.589 * (σ_x + σ_y)
        cep = 0.589 * (sigma_lat + sigma_lon)
        
        # DRMS (Distance RMS) - ~63-68% confidence
        drms = np.sqrt(sigma_lat**2 + sigma_lon**2)
        
        # R95 (95% confidence circle)
        # R95 ≈ 2.45 * sigma for circular error
        r95 = 2.45 * (sigma_lat + sigma_lon) / 2
        
        return ErrorMetrics(
            hdop=hdop,
            vdop=0.0,  # Not applicable for 2D
            pdop=hdop,  # For 2D, PDOP = HDOP
            gdop=hdop,
            sigma_lat_nm=sigma_lat,
            sigma_lon_nm=sigma_lon,
            cep_nm=cep,
            r95_nm=r95,
            drms_nm=drms
        )
        
    except np.linalg.LinAlgError:
        return ErrorMetrics(
            hdop=float('inf'), vdop=float('inf'), pdop=float('inf'),
            gdop=float('inf'), sigma_lat_nm=float('inf'), sigma_lon_nm=float('inf'),
            cep_nm=float('inf'), r95_nm=float('inf'), drms_nm=float('inf')
        )


def compute_confidence_ellipse(
    azimuths_deg: List[float],
    observation_error_arcmin: float = 1.0,
    center_lat: float = 0.0,
    center_lon: float = 0.0,
    confidence_level: float = 0.95
) -> ConfidenceEllipse:
    """
    Compute confidence ellipse for position uncertainty.
    
    The confidence ellipse represents the region within which the true
    position lies with specified probability, accounting for the
    geometry of observations.
    
    Reference: Hoover (1984): Algorithms for Confidence Circles and Ellipses
    
    Args:
        azimuths_deg: Azimuths to observed bodies
        observation_error_arcmin: Standard error per observation (arcminutes)
        center_lat: Fix latitude (for output)
        center_lon: Fix longitude (for output)
        confidence_level: Confidence probability (default 0.95)
        
    Returns:
        ConfidenceEllipse with ellipse parameters
    """
    n = len(azimuths_deg)
    
    if n < 2:
        return ConfidenceEllipse(
            center_lat=center_lat,
            center_lon=center_lon,
            semi_major_nm=float('inf'),
            semi_minor_nm=float('inf'),
            orientation_deg=0.0,
            confidence_level=confidence_level,
            area_sq_nm=float('inf')
        )
    
    # Build geometry matrix
    G = np.zeros((n, 2))
    for i, az in enumerate(azimuths_deg):
        az_rad = np.radians(az)
        G[i, 0] = np.cos(az_rad)  # North
        G[i, 1] = np.sin(az_rad)  # East
    
    try:
        # Covariance matrix of position errors
        # Cov = σ² (G^T G)^-1
        GtG = G.T @ G
        GtG_inv = inv(GtG)
        
        sigma = observation_error_arcmin  # 1' = 1 nm
        Cov = (sigma ** 2) * GtG_inv
        
        # Eigenvalue decomposition of covariance matrix
        eigenvalues, eigenvectors = np.linalg.eig(Cov)
        
        # Sort by eigenvalue (largest first)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        # Chi-squared scaling for confidence level
        # For 2D, use chi2 distribution with 2 degrees of freedom
        chi2_scale = stats.chi2.ppf(confidence_level, df=2)
        
        # Semi-axes of confidence ellipse
        semi_major = np.sqrt(eigenvalues[0] * chi2_scale)
        semi_minor = np.sqrt(eigenvalues[1] * chi2_scale)
        
        # Orientation angle (from North, clockwise)
        orientation = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
        
        # Normalize to 0-180 (ellipse is symmetric)
        if orientation < 0:
            orientation += 180
        
        # Area of ellipse
        area = np.pi * semi_major * semi_minor
        
        return ConfidenceEllipse(
            center_lat=center_lat,
            center_lon=center_lon,
            semi_major_nm=semi_major,
            semi_minor_nm=semi_minor,
            orientation_deg=orientation,
            confidence_level=confidence_level,
            area_sq_nm=area
        )
        
    except np.linalg.LinAlgError:
        return ConfidenceEllipse(
            center_lat=center_lat,
            center_lon=center_lon,
            semi_major_nm=float('inf'),
            semi_minor_nm=float('inf'),
            orientation_deg=0.0,
            confidence_level=confidence_level,
            area_sq_nm=float('inf')
        )


def error_propagation_altitude(
    sextant_error_arcmin: float = 0.3,
    index_error_arcmin: float = 0.1,
    timing_error_sec: float = 1.0,
    dip_error_arcmin: float = 0.2,
    refraction_error_arcmin: float = 0.2
) -> float:
    """
    Calculate total observation error from component errors.
    
    Uses root-sum-square (RSS) combination of independent error sources.
    
    Reference: 
        - Ross (1994): Minimizing Errors in Celestial Positioning
        - Gordon (1964): The Attainment of Precision
    
    Args:
        sextant_error_arcmin: Random sextant reading error
        index_error_arcmin: Residual index error after correction
        timing_error_sec: Time error (converted to altitude error)
        dip_error_arcmin: Dip calculation error
        refraction_error_arcmin: Refraction model error
        
    Returns:
        Total altitude observation error in arcminutes
    """
    # Timing error contribution
    # Earth rotates 15°/hour = 0.25'/second at equator
    # Effect on altitude depends on azimuth, use average factor
    timing_contrib = timing_error_sec * 0.25 * 0.7  # Approximate factor
    
    # RSS combination
    total_error = np.sqrt(
        sextant_error_arcmin**2 +
        index_error_arcmin**2 +
        timing_contrib**2 +
        dip_error_arcmin**2 +
        refraction_error_arcmin**2
    )
    
    return total_error


def optimal_observation_angles(n_observations: int) -> List[float]:
    """
    Calculate optimal azimuth angles for n observations.
    
    Optimal geometry has observations evenly distributed in azimuth.
    This minimizes HDOP and produces a circular error distribution.
    
    Reference: Swaszek et al. (2019)
    
    Args:
        n_observations: Number of observations planned
        
    Returns:
        List of optimal azimuths in degrees
    """
    return [i * 360.0 / n_observations for i in range(n_observations)]


def analyze_geometry_quality(azimuths_deg: List[float]) -> dict:
    """
    Analyze the quality of observation geometry.
    
    Provides diagnostic information about the observation set.
    
    Args:
        azimuths_deg: Actual observation azimuths
        
    Returns:
        Dictionary with geometry analysis
    """
    n = len(azimuths_deg)
    
    if n < 2:
        return {
            'n_observations': n,
            'geometry_quality': 'insufficient',
            'hdop': float('inf'),
            'recommendation': 'Need at least 2 observations'
        }
    
    # Calculate HDOP
    hdop = calculate_hdop_from_azimuths(azimuths_deg)
    
    # Calculate azimuth spread
    azimuths_sorted = sorted(azimuths_deg)
    gaps = []
    for i in range(len(azimuths_sorted)):
        gap = (azimuths_sorted[(i+1) % n] - azimuths_sorted[i]) % 360
        gaps.append(gap)
    max_gap = max(gaps)
    
    # Optimal HDOP for n observations (evenly spaced)
    optimal_hdop = calculate_hdop_from_azimuths(optimal_observation_angles(n))
    hdop_ratio = hdop / optimal_hdop if optimal_hdop > 0 else float('inf')
    
    # Quality assessment
    if hdop < 1.5:
        quality = 'excellent'
    elif hdop < 2.5:
        quality = 'good'
    elif hdop < 4.0:
        quality = 'acceptable'
    elif hdop < 6.0:
        quality = 'marginal'
    else:
        quality = 'poor'
    
    # Recommendations
    recommendations = []
    if max_gap > 180:
        recommendations.append(f"Large gap of {max_gap:.0f}° in azimuth coverage")
    if hdop_ratio > 1.5:
        recommendations.append("Consider observations with better azimuth distribution")
    if n < 3:
        recommendations.append("Additional observations would improve reliability")
    
    return {
        'n_observations': n,
        'geometry_quality': quality,
        'hdop': hdop,
        'optimal_hdop': optimal_hdop,
        'hdop_ratio': hdop_ratio,
        'max_azimuth_gap_deg': max_gap,
        'azimuths': azimuths_deg,
        'recommendations': recommendations if recommendations else ['Geometry is acceptable']
    }


def monte_carlo_position_error(
    azimuths_deg: List[float],
    observation_error_arcmin: float = 1.0,
    n_simulations: int = 10000
) -> dict:
    """
    Estimate position error distribution using Monte Carlo simulation.
    
    Simulates random observation errors and computes resulting
    position fix errors to characterize the error distribution.
    
    Args:
        azimuths_deg: Observation azimuths
        observation_error_arcmin: Standard deviation of altitude errors
        n_simulations: Number of Monte Carlo trials
        
    Returns:
        Dictionary with error statistics
    """
    n_obs = len(azimuths_deg)
    
    if n_obs < 2:
        return {'error': 'Insufficient observations'}
    
    # Build geometry matrix
    G = np.zeros((n_obs, 2))
    for i, az in enumerate(azimuths_deg):
        az_rad = np.radians(az)
        G[i, 0] = np.cos(az_rad)
        G[i, 1] = np.sin(az_rad)
    
    try:
        # Pseudo-inverse for least squares
        G_pinv = np.linalg.pinv(G)
    except:
        return {'error': 'Cannot compute pseudo-inverse'}
    
    # Monte Carlo simulation
    position_errors = np.zeros((n_simulations, 2))
    
    for sim in range(n_simulations):
        # Generate random altitude errors
        altitude_errors = np.random.normal(0, observation_error_arcmin, n_obs)
        
        # Compute position error (in nautical miles)
        pos_error = G_pinv @ altitude_errors
        position_errors[sim] = pos_error
    
    # Compute statistics
    lat_errors = position_errors[:, 0]
    lon_errors = position_errors[:, 1]
    total_errors = np.sqrt(lat_errors**2 + lon_errors**2)
    
    return {
        'n_simulations': n_simulations,
        'observation_error_arcmin': observation_error_arcmin,
        'mean_lat_error_nm': np.mean(lat_errors),
        'mean_lon_error_nm': np.mean(lon_errors),
        'std_lat_nm': np.std(lat_errors),
        'std_lon_nm': np.std(lon_errors),
        'mean_total_error_nm': np.mean(total_errors),
        'median_error_nm': np.median(total_errors),
        'percentile_50_nm': np.percentile(total_errors, 50),
        'percentile_95_nm': np.percentile(total_errors, 95),
        'percentile_99_nm': np.percentile(total_errors, 99),
        'max_error_nm': np.max(total_errors)
    }


if __name__ == "__main__":
    print("=== Error Analysis Module Test ===\n")
    
    # Test HDOP calculation
    print("1. HDOP Calculation:")
    print("-" * 40)
    
    # Optimal geometry (4 evenly spaced observations)
    optimal_azimuths = [0, 90, 180, 270]
    hdop_optimal = calculate_hdop_from_azimuths(optimal_azimuths)
    print(f"  Optimal (4 at 90° spacing): HDOP = {hdop_optimal:.3f}")
    
    # Poor geometry (clustered observations)
    poor_azimuths = [30, 45, 60, 75]
    hdop_poor = calculate_hdop_from_azimuths(poor_azimuths)
    print(f"  Poor (4 clustered): HDOP = {hdop_poor:.3f}")
    
    # Test confidence ellipse
    print("\n2. Confidence Ellipse:")
    print("-" * 40)
    
    ellipse = compute_confidence_ellipse(
        azimuths_deg=[45, 135, 225, 315],
        observation_error_arcmin=1.0,
        confidence_level=0.95
    )
    print(f"  Semi-major axis: {ellipse.semi_major_nm:.2f} nm")
    print(f"  Semi-minor axis: {ellipse.semi_minor_nm:.2f} nm")
    print(f"  Orientation: {ellipse.orientation_deg:.1f}°")
    print(f"  Area: {ellipse.area_sq_nm:.2f} sq nm")
    
    # Test error propagation
    print("\n3. Error Propagation:")
    print("-" * 40)
    
    total_error = error_propagation_altitude(
        sextant_error_arcmin=0.3,
        index_error_arcmin=0.1,
        timing_error_sec=1.0,
        dip_error_arcmin=0.2,
        refraction_error_arcmin=0.2
    )
    print(f"  Component errors:")
    print(f"    Sextant: 0.3'")
    print(f"    Index: 0.1'")
    print(f"    Timing (1s): ~0.2'")
    print(f"    Dip: 0.2'")
    print(f"    Refraction: 0.2'")
    print(f"  Total (RSS): {total_error:.2f}'")
    
    # Test geometry analysis
    print("\n4. Geometry Analysis:")
    print("-" * 40)
    
    analysis = analyze_geometry_quality([30, 120, 210, 300])
    print(f"  Quality: {analysis['geometry_quality']}")
    print(f"  HDOP: {analysis['hdop']:.2f}")
    print(f"  Optimal HDOP: {analysis['optimal_hdop']:.2f}")
    print(f"  Max gap: {analysis['max_azimuth_gap_deg']:.0f}°")
    
    # Test Monte Carlo
    print("\n5. Monte Carlo Simulation:")
    print("-" * 40)
    
    mc_results = monte_carlo_position_error(
        azimuths_deg=[0, 90, 180, 270],
        observation_error_arcmin=1.0,
        n_simulations=10000
    )
    print(f"  Mean error: {mc_results['mean_total_error_nm']:.2f} nm")
    print(f"  50th percentile: {mc_results['percentile_50_nm']:.2f} nm")
    print(f"  95th percentile: {mc_results['percentile_95_nm']:.2f} nm")
