# StarFix - Open-Source Celestial Navigation Algorithm

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An open-source Python-based sight reduction algorithm for celestial navigation applications. The algorithm integrates high-precision ephemeris calculations using JPL DE440 data, altitude corrections for all celestial body types, and multi-body position fixing using iterative least squares optimization with SVD for numerical stability.

## Features

- **High-Precision Ephemeris**: Sun, Moon, planets (Venus, Mars, Jupiter, Saturn), and 57 navigation stars using JPL DE440
- **Complete Altitude Corrections**: Dip, refraction, parallax, semidiameter, and Moon augmentation
- **Multi-Body Position Fixing**: Iterative least squares with SVD for numerical stability
- **HDOP Calculation**: Horizontal Dilution of Precision for observation quality assessment
- **Monte Carlo Simulation**: Error analysis and uncertainty quantification

## Performance

| Metric | Value |
|--------|-------|
| Position accuracy (1.0' obs error, optimal geometry) | 0.89 nm |
| Ephemeris accuracy (stellar) | < 0.6' |
| Convergence | 2-4 iterations |
| Execution time | < 2 ms |

## Installation

```bash
git clone https://github.com/vahitcalisir-art/starfix.git
cd starfix
pip install -r requirements.txt
```

### Dependencies

- Python 3.10+
- NumPy 1.26+
- SciPy 1.12+
- Skyfield 1.48+
- Astropy 6.0+
- Pandas 2.1+

## Quick Start

```python
from src.ephemeris import CelestialBody
from src.sight_reduction import SightReduction
from src.position_fix import PositionFix

# Calculate ephemeris for a star
body = CelestialBody("Sirius")
gha, dec = body.get_gha_dec("2025-06-15 03:00:00")

# Perform sight reduction
sr = SightReduction(lat=34.0, lon=-120.0)
hc, zn = sr.calculate(gha, dec)

# Multi-body position fix
fix = PositionFix()
fix.add_observation(ho=45.5, gha=120.0, dec=23.4)
fix.add_observation(ho=35.2, gha=210.0, dec=-15.6)
fix.add_observation(ho=55.8, gha=330.0, dec=45.2)
lat, lon, hdop = fix.solve(dr_lat=34.0, dr_lon=-120.0)
```

## Project Structure

```
starfix/
├── src/
│   ├── ephemeris.py        # Celestial body position calculations
│   ├── sight_reduction.py  # Altitude and azimuth calculations
│   ├── position_fix.py     # Multi-body least squares solver
│   └── error_analysis.py   # HDOP and Monte Carlo simulation
├── tests/                   # Comprehensive test suite
├── data/                    # Star catalogs and test data
└── results/                 # Validation results
```

## Validation

The algorithm has been validated through:

1. **Ephemeris Accuracy**: Computed positions vs. Nautical Almanac values
2. **Sight Reduction**: Analytical verification at special geometries
3. **Two-Body Fix**: Five globally distributed locations (both hemispheres)
4. **Multi-Body Fix**: 3-6 observations with noise levels 0.0-1.0'
5. **Monte Carlo**: 10,000 simulations per configuration

## Citation

If you use this software in your research, please cite:

```bibtex
@software{starfix2025,
  author = {Çalışır, Vahit},
  title = {StarFix: Open-Source Celestial Navigation Algorithm},
  year = {2025},
  url = {https://github.com/vahitcalisir-art/starfix}
}
```

## Author

**Vahit Çalışır**  
Iskenderun Technical University  
Maritime Transportation Engineering Department  
📧 vahit.calisir@iste.edu.tr  
🆔 [ORCID: 0000-0001-6575-8988](https://orcid.org/0000-0001-6575-8988)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work builds upon the theoretical foundations established by Chiesa (1990), Kaplan (1995), and Nguyen (2014). Ephemeris data provided by JPL DE440 via the Skyfield library.
