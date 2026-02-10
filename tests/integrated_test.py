"""
Integrated End-to-End Validation with Real Ephemeris

This test uses the actual ephemeris module to create realistic observations
and validates the complete pipeline from star selection to position fix.
"""

import numpy as np
import pandas as pd
from datetime import datetime, timezone
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ephemeris import EphemerisCalculator
from sight_reduction import apply_altitude_corrections, BodyType
from position_fix import Observation, multi_body_fix_lsq, _compute_hc_zn
from error_analysis import calculate_hdop_from_azimuths, monte_carlo_position_error


def run_integrated_test():
    """
    Run integrated end-to-end test using real ephemeris data.
    """
    print("="*60)
    print("INTEGRATED END-TO-END VALIDATION")
    print("="*60)
    
    # Initialize ephemeris
    calc = EphemerisCalculator()
    
    # Test scenario: Evening star observations at sea
    test_time = datetime(2025, 6, 15, 3, 0, 0, tzinfo=timezone.utc)  # Evening in Pacific
    true_lat = 34.0
    true_lon = -135.0
    height_of_eye = 3.0  # meters
    
    print(f"\nTest Scenario:")
    print(f"  Time: {test_time}")
    print(f"  True Position: {true_lat}°N, {abs(true_lon)}°W")
    print(f"  Height of Eye: {height_of_eye} m")
    
    # Get visible stars
    visible_stars = calc.get_visible_stars(
        test_time, true_lat, true_lon,
        min_altitude=15.0, max_altitude=70.0
    )
    
    print(f"\n  Visible Navigation Stars: {len(visible_stars)}")
    
    if len(visible_stars) < 3:
        print("  Not enough visible stars for this test")
        return None
    
    # Select best 4 stars for observation
    selected_stars = visible_stars[:4]
    
    print("\n  Selected Stars for Observation:")
    observations = []
    azimuths = []
    
    for i, star in enumerate(selected_stars):
        # Calculate true Hc from navigation triangle (this is what algorithm computes)
        true_Hc, _ = _compute_hc_zn(true_lat, true_lon, star.gha, star.dec)
        
        # Simulate sextant observation: Ho = true_Hc + observation_error
        # Use different seed for each star
        np.random.seed(42 + i)
        observation_error = np.random.normal(0, 0.5/60)  # 0.5' error in degrees
        Ho = true_Hc + observation_error
        
        print(f"    {star.name}: GHA={star.gha:.2f}°, Dec={star.dec:.2f}°, "
              f"Alt={true_Hc:.1f}°, Az={star.azimuth:.1f}°")
        
        observations.append(Observation(
            body_name=star.name,
            gha=star.gha,
            dec=star.dec,
            Ho=Ho
        ))
        azimuths.append(star.azimuth)
    
    
    # Calculate geometry quality
    hdop = calculate_hdop_from_azimuths(azimuths)
    print(f"\n  Observation Geometry:")
    print(f"    HDOP: {hdop:.2f}")
    print(f"    Azimuths: {[f'{a:.0f}°' for a in azimuths]}")
    
    # Perform multi-body fix
    print("\n  Multi-Body Fix Results:")
    
    # Start with DR position offset by ~5 nm
    dr_lat = true_lat + 0.08  # ~5 nm north
    dr_lon = true_lon - 0.05  # ~3 nm west
    
    fix = multi_body_fix_lsq(
        observations,
        initial_lat=dr_lat,
        initial_lon=dr_lon,
        max_iterations=20,
        tolerance_nm=0.001
    )
    
    # Calculate error
    lat_error_nm = (fix.latitude - true_lat) * 60
    lon_error_nm = (fix.longitude - true_lon) * 60 * np.cos(np.radians(true_lat))
    total_error = np.sqrt(lat_error_nm**2 + lon_error_nm**2)
    
    print(f"    DR Position: {dr_lat:.4f}°N, {abs(dr_lon):.4f}°W")
    print(f"    Fix Position: {fix.latitude:.4f}°N, {abs(fix.longitude):.4f}°W")
    print(f"    True Position: {true_lat:.4f}°N, {abs(true_lon):.4f}°W")
    print(f"    Position Error: {total_error:.2f} nm")
    print(f"    RMS Residual: {fix.residual_rms:.2f}'")
    print(f"    Iterations: {fix.iterations}")
    print(f"    Converged: {fix.converged}")
    
    # Compare with Monte Carlo prediction
    print("\n  Monte Carlo Comparison (10,000 trials):")
    mc = monte_carlo_position_error(azimuths, observation_error_arcmin=0.5, n_simulations=10000)
    print(f"    Predicted Mean Error: {mc['mean_total_error_nm']:.2f} nm")
    print(f"    Predicted 95%: {mc['percentile_95_nm']:.2f} nm")
    print(f"    Actual Error: {total_error:.2f} nm")
    
    # Test with multiple noise realizations
    print("\n  Repeated Trials (20 different noise realizations):")
    errors = []
    
    for trial in range(20):
        trial_obs = []
        for i, star in enumerate(selected_stars):
            np.random.seed(trial * 10 + i + 100)
            # Calculate true Hc and add observation error
            true_Hc, _ = _compute_hc_zn(true_lat, true_lon, star.gha, star.dec)
            observation_error = np.random.normal(0, 0.5/60)  # 0.5' error
            Ho = true_Hc + observation_error
            trial_obs.append(Observation(star.name, star.gha, star.dec, Ho))
        
        trial_fix = multi_body_fix_lsq(trial_obs, dr_lat, dr_lon)
        
        lat_err = (trial_fix.latitude - true_lat) * 60
        lon_err = (trial_fix.longitude - true_lon) * 60 * np.cos(np.radians(true_lat))
        errors.append(np.sqrt(lat_err**2 + lon_err**2))
    
    errors = np.array(errors)
    print(f"    Mean Error: {np.mean(errors):.2f} nm")
    print(f"    Std Dev: {np.std(errors):.2f} nm")
    print(f"    Min/Max: {np.min(errors):.2f} / {np.max(errors):.2f} nm")
    print(f"    Expected (MC): {mc['mean_total_error_nm']:.2f} ± {mc['std_lat_nm']:.2f} nm")
    
    # Summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    
    if np.mean(errors) < 2.0:  # Expected ~0.5-1.0 nm with 0.5' observation error
        print("  STATUS: PASS")
        print("  Algorithm produces positions within expected error bounds.")
    else:
        print("  STATUS: NEEDS REVIEW")
        print("  Errors higher than expected - check algorithm.")
    
    # Save results
    results = {
        'test_time': test_time.isoformat(),
        'true_lat': true_lat,
        'true_lon': true_lon,
        'n_stars': len(selected_stars),
        'hdop': hdop,
        'single_trial_error_nm': total_error,
        'mean_error_20_trials_nm': np.mean(errors),
        'std_error_nm': np.std(errors),
        'mc_predicted_mean_nm': mc['mean_total_error_nm'],
        'mc_predicted_95pct_nm': mc['percentile_95_nm']
    }
    
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.DataFrame([results])
    df.to_csv(os.path.join(output_dir, '09_integrated_validation.csv'), index=False)
    print(f"\n  Results saved to: results/09_integrated_validation.csv")
    
    return results


if __name__ == "__main__":
    run_integrated_test()
