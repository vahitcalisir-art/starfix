"""
Position Fix Module for Celestial Navigation

Implements multiple position fixing methods:
1. Two-body direct fix (Chiesa & Chiesa, Gery)
2. Multi-body least squares fix (SVD decomposition)
3. Running fix with motion correction
4. Iterative refinement

References:
    - Chiesa & Chiesa (1990): A Mathematical Method of obtaining an Astronomical Vessel Position
    - Gery (1997): The Direct Fix of Latitude and Longitude from Two Observed Altitudes
    - Nguyen & Im (2014): SVD-Least Square Algorithm
    - Kaplan (1995): Determining the Position of a Vessel from Celestial Observations
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
from scipy import linalg


@dataclass
class Observation:
    """A single celestial observation."""
    body_name: str
    gha: float          # Greenwich Hour Angle (degrees)
    dec: float          # Declination (degrees)
    Ho: float           # Observed altitude after corrections (degrees)
    time_offset: float = 0.0  # Time offset from reference in hours (for running fix)
    weight: float = 1.0       # Observation weight for least squares


@dataclass
class PositionFix:
    """Result of position fix calculation."""
    latitude: float
    longitude: float
    method: str
    iterations: int = 1
    residual_rms: float = 0.0  # RMS of altitude residuals (arcminutes)
    converged: bool = True
    observations_used: int = 0
    hdop: float = 0.0  # Horizontal Dilution of Precision


@dataclass
class RunningFixResult:
    """Result of running fix with motion correction."""
    latitude: float
    longitude: float
    course: float       # Course in degrees (if estimated)
    speed: float        # Speed in knots (if estimated)
    residual_rms: float
    converged: bool


def two_body_fix_direct(
    obs1: Observation, 
    obs2: Observation,
    dr_lat: float = None
) -> PositionFix:
    """
    Calculate position from two celestial observations using direct method.
    
    This implements the analytical solution for the intersection of two
    circles of equal altitude (position circles).
    
    The method solves the system of two spherical equations:
        sin(Ho1) = sin(lat)*sin(dec1) + cos(lat)*cos(dec1)*cos(LHA1)
        sin(Ho2) = sin(lat)*sin(dec2) + cos(lat)*cos(dec2)*cos(LHA2)
    
    Reference: Chiesa & Chiesa (1990), Gery (1997)
    
    Args:
        obs1: First observation
        obs2: Second observation
        
    Returns:
        PositionFix with calculated position
        
    Note:
        Returns the solution closer to the assumed DR position if available,
        otherwise returns the fix in the expected hemisphere.
    """
    # Convert to radians
    dec1 = np.radians(obs1.dec)
    dec2 = np.radians(obs2.dec)
    Ho1 = np.radians(obs1.Ho)
    Ho2 = np.radians(obs2.Ho)
    
    # Geographical positions (substellar points) of bodies
    gp1_lat = obs1.dec  # GP latitude = declination
    gp1_lon = -obs1.gha  # GP longitude = -GHA (West positive convention)
    gp2_lat = obs2.dec
    gp2_lon = -obs2.gha
    
    # Convert GPs to ECEF unit vectors
    def latlon_to_ecef(lat_deg, lon_deg):
        lat = np.radians(lat_deg)
        lon = np.radians(lon_deg)
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        return np.array([x, y, z])
    
    # Unit vectors to GPs
    n1 = latlon_to_ecef(gp1_lat, gp1_lon)
    n2 = latlon_to_ecef(gp2_lat, gp2_lon)
    
    # Zenith distances (arc from zenith to body)
    zd1 = np.pi/2 - Ho1  # zenith distance = 90° - altitude
    zd2 = np.pi/2 - Ho2
    
    # The observer lies on both circles:
    # n1 · P = cos(zd1)
    # n2 · P = cos(zd2)
    # |P| = 1
    
    # Cross product gives direction perpendicular to both
    n_cross = np.cross(n1, n2)
    n_cross_norm = np.linalg.norm(n_cross)
    
    if n_cross_norm < 1e-10:
        # Bodies are at same or opposite positions - no unique solution
        raise ValueError("Bodies are collinear - cannot compute fix")
    
    n_cross = n_cross / n_cross_norm
    
    # Find the two intersection points
    # Using the formula for intersection of two small circles on a sphere
    
    cos_zd1 = np.cos(zd1)
    cos_zd2 = np.cos(zd2)
    
    n1_dot_n2 = np.dot(n1, n2)
    
    # Coefficients for linear combination
    denom = 1 - n1_dot_n2**2
    
    if abs(denom) < 1e-10:
        raise ValueError("Cannot solve - bodies too close together")
    
    alpha = (cos_zd1 - cos_zd2 * n1_dot_n2) / denom
    beta = (cos_zd2 - cos_zd1 * n1_dot_n2) / denom
    
    # In-plane component
    p_plane = alpha * n1 + beta * n2
    
    # Check if solution exists
    p_plane_norm_sq = np.dot(p_plane, p_plane)
    
    if p_plane_norm_sq > 1:
        # No intersection (circles don't intersect)
        raise ValueError("Circles do not intersect - check observations")
    
    # Out-of-plane component magnitude
    gamma = np.sqrt(1 - p_plane_norm_sq)
    
    # Two possible solutions
    P1 = p_plane + gamma * n_cross
    P2 = p_plane - gamma * n_cross
    
    # Convert back to lat/lon
    def ecef_to_latlon(P):
        lat = np.degrees(np.arcsin(P[2]))
        lon = np.degrees(np.arctan2(P[1], P[0]))
        return lat, lon
    
    lat1, lon1 = ecef_to_latlon(P1)
    lat2, lon2 = ecef_to_latlon(P2)
    
    # Choose the more reasonable solution based on DR position
    # If DR is provided, choose closest solution; otherwise use hemisphere
    
    if dr_lat is not None:
        # Choose solution closest to DR latitude
        if abs(lat1 - dr_lat) <= abs(lat2 - dr_lat):
            final_lat, final_lon = lat1, lon1
        else:
            final_lat, final_lon = lat2, lon2
    else:
        # Default: choose solution in same hemisphere as average declination
        avg_dec = (obs1.dec + obs2.dec) / 2.0
        if (lat1 * avg_dec) >= 0:  # Same sign as avg declination
            final_lat, final_lon = lat1, lon1
        else:
            final_lat, final_lon = lat2, lon2
    
    # Normalize longitude to -180 to 180
    if final_lon > 180:
        final_lon -= 360
    elif final_lon < -180:
        final_lon += 360
    
    return PositionFix(
        latitude=final_lat,
        longitude=final_lon,
        method="two_body_direct",
        iterations=1,
        residual_rms=0.0,  # Direct method has no residual
        converged=True,
        observations_used=2,
        hdop=0.0  # Not calculated for two-body
    )


def multi_body_fix_lsq(
    observations: List[Observation],
    initial_lat: float = 0.0,
    initial_lon: float = 0.0,
    max_iterations: int = 20,
    tolerance_nm: float = 0.01
) -> PositionFix:
    """
    Calculate position from multiple observations using least squares.
    
    This method iteratively refines the position estimate by minimizing
    the sum of squared altitude residuals. Uses SVD for numerical stability.
    
    The linearized observation equation:
        ΔHo = (∂Hc/∂lat) Δlat + (∂Hc/∂lon) Δlon
        
    Partial derivatives (from navigation triangle):
        ∂Hc/∂lat = cos(Zn)  
        ∂Hc/∂lon = -sin(Zn) * cos(lat)  # Note: negative sign verified numerically
    
    Reference: 
        - Nguyen & Im (2014): SVD-Least Square Algorithm
        - Kaplan (1995): Least squares with differential correction
    
    Args:
        observations: List of celestial observations
        initial_lat: Initial latitude estimate (degrees)
        initial_lon: Initial longitude estimate (degrees)
        max_iterations: Maximum iterations
        tolerance_nm: Convergence tolerance in nautical miles
        
    Returns:
        PositionFix with calculated position and statistics
    """
    if len(observations) < 2:
        raise ValueError("At least 2 observations required")
    
    # Initialize position in degrees
    lat = initial_lat
    lon = initial_lon
    
    for iteration in range(max_iterations):
        n_obs = len(observations)
        
        # Build design matrix A and residual vector b
        # A is in units of [degrees Ho per degree position]
        # b is in units of [degrees]
        A = np.zeros((n_obs, 2))
        b = np.zeros(n_obs)
        
        for i, obs in enumerate(observations):
            # Calculate Hc and Zn for current position estimate
            Hc, Zn = _compute_hc_zn(lat, lon, obs.gha, obs.dec)
            
            # Residual in degrees (Ho - Hc)
            b[i] = (obs.Ho - Hc) * obs.weight
            
            # Partial derivatives of Hc with respect to lat/lon
            # From spherical trig: sin(Hc) = sin(lat)sin(dec) + cos(lat)cos(dec)cos(LHA)
            # where LHA = GHA + lon (with lon negative for west)
            # Taking partial derivatives:
            # ∂Hc/∂lat = cos(Zn)  [degrees per degree]
            # ∂Hc/∂lon = -sin(Zn) * cos(lat)  [degrees per degree]
            # Note: Verified numerically - the negative sign is required
            
            Zn_rad = np.radians(Zn)
            lat_rad = np.radians(lat)
            
            # Design matrix coefficients
            A[i, 0] = np.cos(Zn_rad) * obs.weight  # ∂Hc/∂lat
            A[i, 1] = -np.sin(Zn_rad) * np.cos(lat_rad) * obs.weight  # ∂Hc/∂lon
        
        # Solve A * [dlat, dlon]^T = b using pseudo-inverse
        try:
            # Use numpy's lstsq for stability
            solution, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
            
            dlat = solution[0]  # degrees
            dlon = solution[1]  # degrees
            
        except Exception:
            # Fall back to SVD
            try:
                U, s, Vt = linalg.svd(A, full_matrices=False)
                s_inv = np.where(s > 1e-10, 1.0/s, 0.0)
                solution = Vt.T @ np.diag(s_inv) @ (U.T @ b)
                dlat = solution[0]
                dlon = solution[1]
            except:
                return PositionFix(
                    latitude=lat,
                    longitude=lon,
                    method="multi_body_lsq",
                    iterations=iteration + 1,
                    residual_rms=999.0,
                    converged=False,
                    observations_used=n_obs,
                    hdop=0.0
                )
        
        # Apply correction
        lat += dlat
        lon += dlon
        
        # Normalize longitude
        if lon > 180:
            lon -= 360
        elif lon < -180:
            lon += 360
        
        # Check convergence (correction < tolerance)
        correction_nm = np.sqrt((dlat * 60)**2 + (dlon * 60 * np.cos(np.radians(lat)))**2)
        
        if correction_nm < tolerance_nm:
            # Calculate final RMS residual
            residuals = []
            for obs in observations:
                Hc, _ = _compute_hc_zn(lat, lon, obs.gha, obs.dec)
                residuals.append((obs.Ho - Hc) * 60.0)  # arcminutes
            
            rms = np.sqrt(np.mean(np.array(residuals)**2))
            
            # Calculate HDOP from geometry
            hdop = _calculate_hdop(observations, lat, lon)
            
            return PositionFix(
                latitude=lat,
                longitude=lon,
                method="multi_body_lsq",
                iterations=iteration + 1,
                residual_rms=rms,
                converged=True,
                observations_used=n_obs,
                hdop=hdop
            )
    
    # Max iterations reached
    residuals = []
    for obs in observations:
        Hc, _ = _compute_hc_zn(lat, lon, obs.gha, obs.dec)
        residuals.append((obs.Ho - Hc) * 60.0)
    rms = np.sqrt(np.mean(np.array(residuals)**2))
    hdop = _calculate_hdop(observations, lat, lon)
    
    return PositionFix(
        latitude=lat,
        longitude=lon,
        method="multi_body_lsq",
        iterations=max_iterations,
        residual_rms=rms,
        converged=False,
        observations_used=len(observations),
        hdop=hdop
    )


def running_fix(
    observations: List[Observation],
    course_deg: float,
    speed_knots: float,
    initial_lat: float = 0.0,
    initial_lon: float = 0.0,
    reference_time: float = 0.0
) -> PositionFix:
    """
    Calculate position from observations taken at different times.
    
    This advances/retards each observation to a common reference time
    based on the vessel's course and speed, then solves for position.
    
    Each observation's circle of position is rotated to account for
    vessel motion between observation time and reference time.
    
    Reference: 
        - Metcalf (1991): Advancing Celestial Circles of Position
        - Kaplan (1995): Running fix methodology
    
    Args:
        observations: List of observations with time_offset field set
        course_deg: Vessel course in degrees (0-360)
        speed_knots: Vessel speed in knots
        initial_lat: Initial position estimate
        initial_lon: Initial position estimate
        reference_time: Reference time (hours), observations adjusted to this
        
    Returns:
        PositionFix adjusted for vessel motion
    """
    if len(observations) < 2:
        raise ValueError("At least 2 observations required")
    
    # Adjust observations for vessel motion
    adjusted_obs = []
    
    for obs in observations:
        # Time difference from reference
        dt_hours = obs.time_offset - reference_time
        
        # Distance traveled
        distance_nm = speed_knots * dt_hours
        
        # Adjust GHA for vessel motion
        # The GP moves 15°/hour due to Earth rotation, already in GHA
        # We need to advance the observer's position
        
        # Calculate position offset in lat/lon
        course_rad = np.radians(course_deg)
        lat_rad = np.radians(initial_lat)
        
        dlat_deg = (distance_nm / 60.0) * np.cos(course_rad)
        dlon_deg = (distance_nm / 60.0) * np.sin(course_rad) / np.cos(lat_rad)
        
        # Create adjusted observation
        # Effectively, we adjust the GHA to "move" the GP relative to moving observer
        adjusted_gha = obs.gha + dlon_deg  # Simplified adjustment
        
        adjusted_obs.append(Observation(
            body_name=obs.body_name,
            gha=adjusted_gha,
            dec=obs.dec,
            Ho=obs.Ho,
            time_offset=reference_time,
            weight=obs.weight
        ))
    
    # Solve using least squares
    return multi_body_fix_lsq(adjusted_obs, initial_lat, initial_lon)


def _compute_hc_zn(lat_deg: float, lon_deg: float, gha_deg: float, dec_deg: float) -> Tuple[float, float]:
    """
    Compute altitude and azimuth (internal helper).
    
    Args:
        lat_deg: Observer latitude
        lon_deg: Observer longitude
        gha_deg: Greenwich Hour Angle
        dec_deg: Declination
        
    Returns:
        (Hc, Zn) in degrees
    """
    lat = np.radians(lat_deg)
    dec = np.radians(dec_deg)
    lha = np.radians((gha_deg + lon_deg) % 360.0)
    
    # Altitude
    sin_Hc = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(lha)
    sin_Hc = np.clip(sin_Hc, -1.0, 1.0)
    Hc = np.degrees(np.arcsin(sin_Hc))
    
    # Azimuth
    cos_Hc = np.cos(np.radians(Hc))
    if abs(cos_Hc) < 1e-10:
        Zn = 0.0
    else:
        sin_Zn = np.sin(lha) * np.cos(dec) / cos_Hc
        cos_Zn = (np.sin(dec) - np.sin(lat) * np.sin(np.radians(Hc))) / (np.cos(lat) * cos_Hc)
        
        sin_Zn = np.clip(sin_Zn, -1.0, 1.0)
        cos_Zn = np.clip(cos_Zn, -1.0, 1.0)
        
        Zn = np.degrees(np.arctan2(sin_Zn, cos_Zn))
        
        if Zn < 0:
            Zn += 360.0
    
    return Hc, Zn


def _calculate_hdop(observations: List[Observation], lat_deg: float, lon_deg: float) -> float:
    """
    Calculate Horizontal Dilution of Precision (HDOP).
    
    HDOP indicates the geometric quality of the fix. Lower values
    indicate better geometry (observations well-distributed in azimuth).
    
    Reference: Swaszek et al. (2019): Rethinking Star Selection
    
    Args:
        observations: List of observations
        lat_deg: Position latitude
        lon_deg: Position longitude
        
    Returns:
        HDOP value (dimensionless)
    """
    if len(observations) < 2:
        return float('inf')
    
    # Build geometry matrix
    n_obs = len(observations)
    G = np.zeros((n_obs, 2))
    
    for i, obs in enumerate(observations):
        _, Zn = _compute_hc_zn(lat_deg, lon_deg, obs.gha, obs.dec)
        Zn_rad = np.radians(Zn)
        lat_rad = np.radians(lat_deg)
        
        # Unit vector pointing toward GP
        G[i, 0] = np.cos(Zn_rad)  # North component
        G[i, 1] = np.sin(Zn_rad)  # East component (scaled by cos(lat) for distance)
    
    try:
        # (G^T G)^-1 gives the covariance factor
        GtG = G.T @ G
        GtG_inv = np.linalg.inv(GtG)
        
        # HDOP = sqrt(trace of horizontal components)
        hdop = np.sqrt(GtG_inv[0, 0] + GtG_inv[1, 1])
        
        return hdop
    except:
        return float('inf')


def iterative_refinement(
    observations: List[Observation],
    initial_fix: PositionFix,
    refinement_iterations: int = 3
) -> PositionFix:
    """
    Refine a position fix using iterative least squares.
    
    Starting from an initial fix (e.g., two-body), refines the
    position using all available observations.
    
    Args:
        observations: All available observations
        initial_fix: Initial position estimate
        refinement_iterations: Number of refinement passes
        
    Returns:
        Refined PositionFix
    """
    return multi_body_fix_lsq(
        observations,
        initial_lat=initial_fix.latitude,
        initial_lon=initial_fix.longitude,
        max_iterations=refinement_iterations * 5
    )


def select_optimal_observations(
    observations: List[Observation],
    lat_deg: float,
    lon_deg: float,
    n_select: int = 4
) -> List[Observation]:
    """
    Select optimal subset of observations for best fix geometry.
    
    Uses azimuth distribution criterion: observations should be
    well-distributed around the horizon for minimum HDOP.
    
    Reference: Swaszek et al. (2019)
    
    Args:
        observations: Available observations
        lat_deg: Approximate position
        lon_deg: Approximate position
        n_select: Number of observations to select
        
    Returns:
        List of selected observations
    """
    if len(observations) <= n_select:
        return observations
    
    # Calculate azimuths for all observations
    azimuths = []
    for obs in observations:
        _, Zn = _compute_hc_zn(lat_deg, lon_deg, obs.gha, obs.dec)
        azimuths.append(Zn)
    
    # Greedy selection: start with first, add most distant in azimuth
    selected = [0]
    selected_azimuths = [azimuths[0]]
    
    while len(selected) < n_select:
        best_idx = -1
        best_min_dist = -1
        
        for i, az in enumerate(azimuths):
            if i in selected:
                continue
            
            # Minimum angular distance to any selected observation
            min_dist = min(
                min(abs(az - sel_az), 360 - abs(az - sel_az))
                for sel_az in selected_azimuths
            )
            
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_idx = i
        
        if best_idx >= 0:
            selected.append(best_idx)
            selected_azimuths.append(azimuths[best_idx])
    
    return [observations[i] for i in selected]


if __name__ == "__main__":
    print("=== Position Fix Module Test ===\n")
    
    # Test two-body fix
    print("1. Two-Body Direct Fix:")
    print("-" * 40)
    
    obs1 = Observation(
        body_name="Sun",
        gha=45.0,
        dec=23.0,
        Ho=45.0
    )
    
    obs2 = Observation(
        body_name="Sirius",
        gha=120.0,
        dec=-16.7,
        Ho=55.0
    )
    
    try:
        fix = two_body_fix_direct(obs1, obs2)
        print(f"  Position: {fix.latitude:.4f}°N, {abs(fix.longitude):.4f}°W")
        print(f"  Method: {fix.method}")
    except Exception as e:
        print(f"  Error: {e}")
    
    # Test multi-body fix
    print("\n2. Multi-Body Least Squares Fix:")
    print("-" * 40)
    
    # Create test observations (simulated from known position)
    true_lat = 34.0
    true_lon = -120.0
    
    observations = [
        Observation("Star1", gha=30.0, dec=20.0, Ho=0.0),
        Observation("Star2", gha=120.0, dec=-10.0, Ho=0.0),
        Observation("Star3", gha=210.0, dec=45.0, Ho=0.0),
        Observation("Star4", gha=300.0, dec=5.0, Ho=0.0),
    ]
    
    # Calculate what the observed altitudes should be
    for obs in observations:
        Hc, _ = _compute_hc_zn(true_lat, true_lon, obs.gha, obs.dec)
        # Add small "observation error"
        obs.Ho = Hc + np.random.normal(0, 0.02)  # 0.02° ≈ 1.2' error
    
    fix = multi_body_fix_lsq(
        observations,
        initial_lat=35.0,  # Start with offset
        initial_lon=-122.0
    )
    
    print(f"  True Position: {true_lat}°N, {abs(true_lon)}°W")
    print(f"  Calculated: {fix.latitude:.4f}°N, {abs(fix.longitude):.4f}°W")
    print(f"  Error: {np.sqrt((fix.latitude-true_lat)**2 + ((fix.longitude-true_lon)*np.cos(np.radians(true_lat)))**2)*60:.2f} nm")
    print(f"  Iterations: {fix.iterations}")
    print(f"  RMS Residual: {fix.residual_rms:.2f}'")
    print(f"  HDOP: {fix.hdop:.2f}")
    print(f"  Converged: {fix.converged}")
