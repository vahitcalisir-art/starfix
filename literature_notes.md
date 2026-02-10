# Literature Notes - Phase 2

## Research Topic
**Development and Validation of an Open-Source Python-Based Sight Reduction Algorithm for Celestial Navigation**

---

## Accepted Articles

### Article #1
**Citation:** Vanin, G. (2022). *The beginnings of celestial navigation: early techniques and instruments.* Preprint.

**Relevance:** Historical and conceptual foundations of celestial navigation

#### Content for Introduction:
- Historical development of celestial navigation from antiquity to Modern Age
- Transition from qualitative to quantitative navigation methods
- Role of astronomical observations in maritime exploration

#### Content for Literature Review:
- **Latitude determination methods:**
  - Polaris altitude measurement
  - Sun's meridian altitude with declination tables
  - Regimento do Norte (North Star correction rules based on Guards position - 8 positions)
  
- **Fundamental relationship:** Observed altitude + declination → latitude calculation

- **Historical instruments:**
  - Nautical quadrant (~25 cm radius, wooden/brass)
  - Nautical astrolabe (evolved from planispheric astrolabe)
  - Cross-staff (adapted from Levi Ben Gerson's 1342 design)
  - Kamal (Arab navigation instrument with knotted string)

#### Historical Accuracy Benchmarks:
| Source | Period | Average Error |
|--------|--------|---------------|
| Magellan voyage (Pigafetta) | 1519-1522 | 16' (~0.27°) |
| Magellan voyage (Albo) | 1519-1522 | 18' (~0.30°) |
| Regiomontanus astronomical | ~1480 | 19' (~0.32°) |
| Cross-staff at sea (modern test) | - | ~20' (~0.33°) |

#### Key Quotes:
- "before even the nautical cross-staff came into use, the first phase of celestial navigation aiming at measuring latitude at sea with sufficient precision, came to an end 500 years ago"
- Columbus fourth voyage: "error of less than half a degree, the result though of the average of several measurements"

---

### Article #2
**Citation:** Li, C., Chen, Z., Liu, X., Chen, B., Zheng, Y., Tong, S., & Wang, R. (2022). Adaptively robust filtering algorithm for maritime celestial navigation. *The Journal of Navigation, 75*(1), 200–212.

**Relevance:** Modern computational algorithms for maritime celestial navigation

#### Content for Introduction:
- "Celestial navigation is an autonomous navigation technology based on the inertial frame of indestructible natural celestial bodies"
- Advantages: no electromagnetic interference, no dependence on external facilities, no accumulation of navigation errors over time
- Still important despite GNSS due to electromagnetic interference vulnerability

#### Content for Literature Review:
- **Historical Development:**
  - 1875: Altitude difference method proposed by Saint-Hilaire (marks maturity of modern celestial navigation)
  - Traditional method: observe celestial bodies with sextant, calculate AVP from measured vs calculated altitude
  - Evolution: manual observation → automatic star sensors

- **Previous Computational Methods:**
  - Position circles intersection method (Chiesa & Chiesa, 1990)
  - Astronomical running fix algorithm (Williams, 1990; Matti, 1990)
  - SEEM - Simultaneous Equal-Altitude Equations Method (Hsu et al., 2003, 2005)
  - Least Squares method for multi-leg voyages (Kaplan, 1995, 1996)
  - NAVSSI - Navigation Sensor System Interface for Navy use (Kaplan, 1999)

#### Mathematical Models for Methodology:
- **Celestial Observation Equation:**
  - Zenith distance as observation: $z_i$ for each star
  - Design matrix elements:
    - $a_i = \frac{-\cos\varphi_0 \sin\delta_i - \sin\varphi_0 \cos\delta_i \cos H_i}{\sqrt{1-P_{i0}^2}}$
    - $b_i = \frac{\cos\varphi_0 \cos\delta_i \sin H_i}{\sqrt{1-P_{i0}^2}}$
    - $P_{i0} = \sin\varphi_0 \sin\delta_i + \cos\varphi_0 \cos\delta_i \cos H_i$
  - Where: $\varphi_0$ = latitude, $\delta_i$ = star declination, $H_i$ = star hour angle

- **Ship Motion Model:**
  - State vector: position $(\varphi, \lambda)$ and velocity $(V_N, V_E)$
  - Kinematic model for uniform motion

#### Accuracy Benchmarks:
| Method | Position RMSE | Speed RMSE | Course RMSE |
|--------|---------------|------------|-------------|
| Least Squares | 0.42 n.miles | - | - |
| Extended Kalman Filter | 0.27 n.miles | 0.19 knots | 0.64° |
| Robust Adaptive EKF | 0.11 n.miles | 0.10 knots | 0.25° |

- Historical benchmark: "standard deviation of celestial AVP observed by sextant was about 1525 m" (Bovens, 1994)

#### Key Algorithms Referenced:
- Least Squares (LS) estimation
- Extended Kalman Filter (EKF)
- Robust Adaptive Extended Kalman Filter (RAEKF)

---

### Article #3
**Citation:** Critchley-Marrows, J.J.R., & Mortari, D. (2023). A Return to the Sextant—Maritime Navigation Using Celestial Bodies and the Horizon. *Sensors, 23*, 4869.

**Relevance:** Modern computational approaches to maritime celestial navigation using horizon and star observations

#### Content for Introduction:
- "The sextant could be argued to be one of the most influential inventions of the second millennium"
- Invented by John Bird in 1759
- GPS spoofing/jamming attacks motivate return to celestial navigation:
  - Black Sea cargo ships reported position outside Moscow
  - Port of Shanghai incidents
  - Red Sea GPS disruptions
- IALA requirements: <100m accuracy for ocean/coastal phases

#### Content for Literature Review:
- **Historical context:**
  - Sextant: optical telescope + measuring tools for angle from star/Sun to horizon
  - Traditional accuracy: "few kilometres"
  - Modern optical systems can achieve 100m level

- **Previous computational methods cited:**
  - Genetic algorithm for celestial navigation (Tsou, 2012)
  - New ideas for celestial navigation (Vulfovich & Fogilev, 2010)
  - Adaptively robust filtering (Li et al., 2022) - cross-reference!

#### Mathematical Models for Methodology:
- **Ellipsoid representation:**
  - WGS-84 model
  - $\frac{X^2}{a^2} + \frac{Y^2}{b^2} + \frac{Z^2}{c^2} = 1$
  - Matrix form: $\mathbf{p}^T \mathbf{A} \mathbf{p} = 1$

- **Four position estimation methods:**
  1. **Nonlinear Least Squares (NLLS):** Iterative solution using $x_{k+1} = x_k - (H^T H)^{-1} H^T h$
  2. **Christian-Robinson Method:** Cholesky factorization to transform ellipsoid to unit sphere
  3. **Near Linear Horizon:** Assumes $s^T g = 0$ (horizon perpendicular to gravity)
  4. **Geodetic Horizon:** Accounts for ellipsoidal curvature - BEST PERFORMANCE

- **Gravity vector (geocentric model):**
  - $g = [\cos\lambda\cos\phi, \sin\lambda\cos\phi, \sin\phi]^T$

#### Accuracy Benchmarks:
| Method | Horizontal Error | East Error | North Error |
|--------|------------------|------------|-------------|
| NLLS | 41.19 km | 34.19 km | 22.97 km |
| Christian-Robinson | 11.47 km | 9.52 km | 6.39 km |
| Near Linear Horizon | 10.97 km | 4.09 km | 10.18 km |
| **Geodetic Horizon** | **0.27 km** | **0.19 km** | **0.13 km** |

- Star tracker accuracy assumed: σ = 10 arc-sec
- Geoid variation ±100m causes only 0.0041 arc-sec error (negligible)
- Best method achieves **sub-100m accuracy** in horizontal plane

#### Key Quotes:
- "performances are demonstrated that exceed the traditional approaches of a sextant, delivering an accuracy and precision under 100 m"
- "The best approach considers a geodetic-based horizon"

#### Technical Notes:
- Atmospheric deflection: max 34 arc-min at horizon
- Astronomical seeing noise: ~1 arc-sec
- Field of view: 40°, Sensor resolution: 2048×2048 px

---

### Article #4
**Citation:** Tsai, K.-C., Tseng, W.-K., Chen, C.-L., & Sun, Y.-J. (2022). A Novel Analytical Solution Method for Celestial Positioning. *Journal of Marine Science and Engineering, 10*, 771.

**Relevance:** Computational algorithms for two-body celestial fix using vector methods - directly applicable to Python implementation

#### Content for Introduction:
- STCW 2010 Manila amendments require seafarers to use "Electronic Nautical Almanac and celestial navigation calculation software"
- Traditional intercept method (Marcq St. Hilaire, 1875) has shortcomings
- GPS spoofing/jamming motivates celestial navigation alternatives
- Two-body fix is "generally a more favorable option for seafarers"

#### Content for Literature Review:
- **Previous computational methods:**
  - Closed analytic solution (Van Allen, 1981)
  - Simultaneous equations (Chen et al., 2003)
  - Vector method + spherical triangle (Gonzalez, 2008)
  - Genetic Algorithm (Tsou, 2012) - cross-reference!
  - Particle Swarm Optimization (Tsou, 2015)
  - Sumner line methods (Chen et al., 2015; Hsu et al., 2017)

#### Mathematical Models for Methodology:
- **ECEF Cartesian Coordinates:**
  - Position vector: $\mathbf{P} = r[\cos(lat)\cos(lng), \cos(lat)\sin(lng), \sin(lat)]^T$
  - GP vector: $\mathbf{GP}_i = [\cos(dec_i)\cos(gha_i), -\cos(dec_i)\sin(gha_i), \sin(dec_i)]$

- **Circle of Position (CoP):**
  - Small circle on unit sphere
  - Zenith distance: $Z_d = 90° - H_o$
  - Two CoPs intersect at two points (one is vessel position)

- **Scenario 1: Newton's Method (Iterative)**
  - Nonlinear system: $F(X) = [\mathbf{GP}_1 \cdot X, \mathbf{GP}_2 \cdot X, X \cdot X]^T - H = 0$
  - Where $H = [\sin(h_1), \sin(h_2), 1]^T$
  - Jacobian: $J_F(X_n) = [\mathbf{GP}_1; \mathbf{GP}_2; X_n/|X_n|]$
  - Iteration: $X_{n+1} = J_F(X_n)^{-1} \cdot H$
  - Convergence criterion: $|X_1 - X_0| < \varepsilon$ (e.g., $10^{-8}$)

- **Scenario 2: Closed Analytical Solution (Direct)**
  - Direction vector: $N_{12} = \mathbf{GP}_1 \times \mathbf{GP}_2$
  - Parametric equation: $X = X_1 + N_{12} \cdot t$
  - Solve quadratic: $|N_{12}|^2 t^2 + 2(X_1 \cdot N_{12})t + |X_1|^2 = 1$
  - Solution: $t = \frac{-(X_1 \cdot N_{12}) \pm \sqrt{(X_1 \cdot N_{12})^2 - |N_{12}|^2(|X_1|^2 - 1)}}{|N_{12}|^2}$

- **Solution conditions:**
  - No solution: $\mathbf{GP}_1 \cdot \mathbf{GP}_2 > \cos(\pi - h_1 - h_2)$
  - One solution: $\mathbf{GP}_1 \cdot \mathbf{GP}_2 = \cos(\pi - h_1 - h_2)$
  - Two solutions: $\mathbf{GP}_1 \cdot \mathbf{GP}_2 < \cos(\pi - h_1 - h_2)$

#### Accuracy Benchmarks:
| Method | Example 1 Result | Execution Time |
|--------|------------------|----------------|
| Intercept Method | 41°38.6'N, 017°08.1'W | Manual |
| Newton's Method | 41°39.135'N, 017°07.313'W | 40-50 μs |
| Closed Analytical | 41°39.135'N, 017°07.313'W | 20-30 μs |

- Both methods achieve same accuracy
- Closed analytical ~10 μs faster (no iteration)
- Newton's method typically converges in 3+ iterations

#### Key Implementation Notes:
- Uses atan2 for quadrant identification and reduced round-off error
- MATLAB and JavaScript implementations demonstrated
- "Coordinate transformation and vector algebra instead of spherical trigonometry"
- Suitable for integration with ECDIS

#### Key Quotes:
- "The two methods' capability to avoid complicated manual calculations and chart work processes makes them highly suitable as an alternative to the use of global navigation satellite systems"
- "The derivation of the equation is straightforward and intuitive"

---

### Article #5
**Citation:** Yang, S., Feng, W., Wang, S., & Li, J. (2022). A SINS/CNS integrated navigation scheme with improved mathematical horizon reference. *Measurement, 195*, 111028.

**Relevance:** Integrated navigation using intercept method (Saint-Hilaire) with mathematical models for celestial positioning

#### Content for Introduction:
- "Since Saint-Hilaire founded the astronomical positioning theory based on the intercept method, the development of celestial positioning technology mainly pursues higher measurement accuracy"
- CNS measures elevation angle and azimuth of celestial bodies in local geographic coordinate system
- Applicable to ships, aircraft, and spacecraft
- Horizon reference accuracy is the bottleneck for high-precision celestial navigation

#### Content for Literature Review:
- **Horizon Reference Methods:**
  1. Direct horizon-sensing (sextant, infrared horizon sensor, inclinometer)
  2. Indirect horizon-sensing (starlight refraction)
  3. SINS-aided method
- Infrared horizon sensors: 0.02° precision → ~2 km positioning error
- Space-sextant CNS: 0.9 km positioning accuracy
- **Integration modes:**
  - Tightly integrated (uses elevation/azimuth directly)
  - Deep integrated (uses Mathematical Horizon Reference)

#### Mathematical Models for Methodology:
- **Coordinate Systems:**
  - Body frame (b)
  - Earth-Centered Inertial (ECI) frame (i)
  - Earth-Centered Earth-Fixed (ECEF) frame (e)
  - Local East-North-Up (ENU) frame (t)

- **Position Vector Conversion:**
  - Spherical coordinates: $(r, \alpha_d, \delta_d)$ = (radius, right ascension, declination)
  - $\bar{r} = \sqrt{\bar{x}^2 + \bar{y}^2 + \bar{z}^2}$
  - $\bar{\alpha}_d = \arctan(\bar{y}/\bar{x})$
  - $\bar{\delta}_d = \arcsin(\bar{z}/\bar{r})$

- **Geographic Position from Inertial:**
  - $\lambda = \alpha_d - \text{GHA}_{Aries}$
  - $L = \delta_d$

- **Direction Cosine Matrix (ENU to ECEF):**
$$C^t_e = \begin{bmatrix} -\sin\lambda & \cos\lambda & 0 \\ -\sin L \cos\lambda & -\sin L \sin\lambda & \cos L \\ \cos L \cos\lambda & \cos L \sin\lambda & \sin L \end{bmatrix}$$

- **Mathematical Horizon Reference:**
  - $C^t_b = C^t_e C^e_i C^i_b$
  - Error correction: $\breve{C}^i_b = (I + \hat{\phi}\times)\bar{C}^i_b$

#### Accuracy Benchmarks:
| Method | Latitude RMSE | Longitude RMSE | Attitude RMSE |
|--------|---------------|----------------|---------------|
| SINS only | 670.95 m | 112.41 m | 0.18' |
| INM-1 (SINS horizon) | 241.02 m | 151.28 m | 0.29' |
| INM-2 (MHR geographic) | 65.4 m | 67.3 m | 0.04' |
| **Proposed MHR** | **21.2 m** | **24.9 m** | **0.005'** |

- Uses Kalman filter for data fusion
- Star sensor noise: 3 arcsec
- Barometric altimeter noise: 5 m

#### Key Implementation Notes:
- GHA_Aries (Greenwich Hour Angle of Aries) needed for coordinate conversion
- Decoupling of position and attitude errors improves accuracy
- Intercept method for astronomical positioning

#### Key Quotes:
- "CNS measures the elevation angle and azimuth of the celestial body in the local geographic coordinate system, and then determines the latitude and longitude information of the observation point"
- "The proposed approach has the best converging speed as well as the highest navigation precision"

---

### Article #6
**Citation:** Villmoare, B. (2022). Determination of latitude by two fixed-altitude sightings. *The Journal of Navigation, 75*(4), 805–812. DOI: 10.1017/S0373463322000406

**Relevance:** Double-altitude latitude determination method with spherical triangle mathematics and computational implementation

#### Content for Introduction:
- Traditional double-altitude / two ex-meridian methods have centuries of history (Lax 1799, Bowditch 1802, Moore 1807)
- Modern navigators still use noon sights with sextant to confirm GPS position
- Computational approaches enable complex calculations that were previously impractical

#### Content for Literature Review:
- **Double-altitude method (two ex-meridians):**
  - Two sightings at equal altitudes, equidistant from meridian passage
  - Times generate longitudinal fix; combined with altitude data → latitude
  
- **Spherical triangle solution:**
  - Uses Equatorial and Horizontal (Azimuthal) coordinate system intersection
  - Three elements: sides b (90° ± declination) and c (90° - altitude), angle C (half angular time between observations)
  - Solves for side a → latitude = 90° - a

#### Mathematical Models (Napier's spherical triangle rules):
```
Angle B = 180° - arcsin[(sin b × sin C) / sin c]

Side a = 2 × arctan{tan[½(b - c)] × [sin(½(B + C)) / sin(½(B - C))]}

Latitude = 90° - a
```
Where:
- b = 90° ± declination of celestial body
- c = 90° - altitude of observation
- C = half the angular time between two sightings (time difference × 15°/hr ÷ 2)

#### Worked Example from Paper:
- Observations: 10:30:36 AM and 1:16:59 PM
- Altitude: 40°, Declination: -10.005°, Time zone: +8 GMT
- Side b: 100.005° (90 - (-10.005))
- Side c: 50° (90 - 40)
- Angle C: 20.7979° (2h 46m 23s = 2.773056 hr × 15°/hr ÷ 2)
- **Result:** Latitude 35°59.4' N, Longitude 118°27' E

#### Computational Implementation:
- Author developing smartphone/tablet app for latitude/longitude calculation
- Excel spreadsheet tables for various sextant angles and declinations
- Solution space analysis performed for 20°, 40°, 60° altitude angles
- Method applicable to traditional marine sextant, bubble sextant, artificial horizon

#### Key Quotes:
- "the need for traditional plotting tools and workspace, typically required for identifying the location of the 'cocked hat', is eliminated"
- "This method produces a result that is a specific position, using latitude and longitude"
- "Because this method is computationally cumbersome, it is most convenient when used in a computer or tablet application, or with tables"

#### Relevance to Research:
- Provides classical spherical triangle formulas adaptable for Python implementation
- Demonstrates computational approach to sight reduction
- References historical methods (Bowditch, Nautical Almanac principles)
- Shows how algorithmic solutions can replace manual plotting methods

---

### Article #7
**Citation:** Barbot, L., Ferrari, M., Montel, J., Roehlli, Y., Gach, J.-L., Thuillot, W., & Dohlen, K. (2022). Towards a Daytime and Low-altitude Stellar Positioning System: Challenges and First Results. *2022 International Technical Meeting of The Institute of Navigation*, 1371–1379. DOI: 10.33012/2022.18263

**Relevance:** Modern sensor-based approach to celestial navigation for maritime applications (MARIS STELLA project)

#### Content for Introduction:
- MARIS STELLA project aims to provide stellar receiver for maritime navigation
- Bridge between traditional sextant navigation and modern sensor-based systems
- Position, attitude, and time determination from star observations

#### Content for Literature Review:
- **Celestial fix principles (computational):**
  - Line of Position is a circle on Earth where same star observed at same altitude
  - Three celestial bodies required for position fix (intersection of three circles)
  - Modern implementation: camera-based star detection replaces human eye + sextant

- **Accuracy analysis:**
  - Angular resolution of 1 arcsec = 31 m position error
  - Multi-star observation enhances accuracy by factor √n (n = number of stars)
  - Star tracker example: 35 arcsec resolution → 0.5 arcsec pointing with 24 stars

- **Star identification algorithms:**
  | Algorithm | Minimum Stars Required |
  |-----------|------------------------|
  | Angle | 2 |
  | Combine triangle | 3 |
  | Pyramid | 4 |
  | Quadrilateral | 4 |
  | Pentagon | 5 |

#### Technical Implementation Notes:
- **Star catalogs:** HIPPARCOS (visible), 2MASS (J, H, K infrared), GAIA
- **Photometric bands:** B, V, R, I (visible), J, H, K (infrared)
- **Nautical Almanac:** Contains 173 navigation stars (81 stars with mag ≤ 2.8)
- **Sensor comparison:**
  | Parameter | Silicon | InGaAs | HgCdTe |
  |-----------|---------|--------|--------|
  | Pixel pitch | ~3 μm | ~15 μm | ~15 μm |
  | Angular resolution | 2 arcsec | 10 arcsec | 10 arcsec |
  | Definition | 9602×6498 | 640×512 | 1280×1024 |
  | Band | 0.3-1.0 μm | 0.9-1.7 μm | 1.0-2.5 μm |

#### Horizon Reference & Corrections:
- Sextant altitude (Hs) requires atmospheric refraction correction
- Additional corrections: temperature, pressure (accuracy 0.1 arcmin)
- Bubble sextant adds Coriolis correction (accuracy 1 arcmin)
- Deflection of Vertical (DoV): ±6 arcsec at sea, ±10 arcsec on land
- EGM2008 gravity model: ~1.3 arcsec accuracy

#### Key Quotes:
- "Celestial navigation on a ship is considered as a backup system if satellite positioning is not available"
- "the line of positions, where the same star can be observed, at the same time, with the same zenithal distance, is a circle on the Earth"
- "Sailors navigating with sextant used to shoot three celestial bodies to fix their position on the sole point of intersection of three circles"

#### Relevance to Research:
- Provides context for modernizing celestial navigation algorithms
- Star catalog integration relevant for Python implementation (HIPPARCOS, 2MASS, GAIA)
- Accuracy benchmarks for validation targets
- Circle of Position / multi-body fix principles applicable to algorithm design
- Horizon reference and correction methods for altitude processing

---

### Article #8
**Citation:** Kotlarić, S. (1975). Simplification in observation and computation of a two-star fix without use of the altitude difference method. *International Hydrographic Review, 52*(1), 157–170.

**Relevance:** Direct Method for two-star fix computation as alternative to intercept method

#### Content for Introduction:
- Direct Method eliminates need for computing altitudes and azimuths for intercepts
- Fix coordinates (Latitude and Longitude) obtained directly from observed star altitudes
- Computational approach designed for practical maritime navigation

#### Content for Literature Review:
- **Direct Method vs. Intercept Method (Marcq St. Hilaire):**
  - Traditional: compute Hc, Az → intercept → plot position lines → find intersection
  - Direct: observed altitudes → Tables K₁₁ → Latitude and Longitude directly
  - Eliminates graphical plotting of "cocked hat"

- **Two-star fix principles:**
  - Two celestial bodies observed at known times
  - Non-simultaneous observations require corrections for elapsed time
  - Ship's run correction between sights (speed × time × heading angle)

#### Mathematical Models:
**Error propagation formulas:**
```
dB = -dHo₁ · sec Ho₂ · cosec(AZ₁ - AZ₂)
dLat = -dm · cos Ho₂ · sin AZ₂
dLHA₂ = -dm · cos Ho₂ · cos AZ₂ · sec Lat
```
Where:
- dHo₁ = error in first star's observed altitude
- dm = error in parallactic angle
- AZ₁, AZ₂ = azimuths of first and second stars

**Altitude correction for elapsed time:**
```
diff.Alt = diff.LHA × cos Lat × sin AZ
```

#### Computational Workflow:
1. Compute approximate LHA Aries from GMT and DR longitude
2. Select star pair from Tables K₁₁ based on Lat and LHA
3. Observe both stars (up to 4 minutes apart)
4. Correct sextant altitudes (Index, dip, refraction)
5. Enter Main Tables with observed altitudes → tabulated Lat, LHA Aries
6. Apply corrections using multiplication table (indices × differences)
7. Add SHA₂ to LHA Aries → LHA₂
8. Longitude = LHA₂ - GHA₂

#### Key Implementation Notes:
- Uses Nautical Almanac for GHA Aries, SHA, declinations
- Correction indices for: altitude (IV₁, IV₂), declination (ID₁, ID₂), SHA difference (ISU)
- Non-simultaneous observation correction via Table IV C
- Discusses future computer/calculator applications for Direct Method

#### Worked Example from Paper:
- Date: June 12, 1968, Zone Time 19h50m
- Stars: SPICA and REGULUS
- Observations: hs₁ = 37°19.1', hs₂ = 36°30.7'
- **Result:** Lat 41°46.1' N, Long 17°07.4' E

#### Key Quotes:
- "Tables K₁₁ enable calculation of the fix coordinates direct from the Tables themselves"
- "eliminate the precise graphical work on the plotting sheet or an approximate sketch"
- "the Direct Method of finding fix coordinates will be of more interest than the St. Hilaire Method for computer applications"

#### Relevance to Research:
- Alternative algorithmic approach to sight reduction
- Mathematical framework directly adaptable for Python implementation
- Error analysis formulas useful for accuracy validation
- Historical bridge between table-based and computational methods
- Explicitly discusses electronic calculator/computer applications

---

### Article #9
**Citation:** Gery, S. W. (1997). The Direct Fix of Latitude and Longitude from Two Observed Altitudes. *NAVIGATION: Journal of The Institute of Navigation, 44*(1), 15–23. DOI: 10.1002/j.2161-4296.1997.tb01935.x

**Relevance:** Complete algorithm for direct two-body fix computation with calculator/computer implementation

#### Content for Introduction:
- Direct method computes latitude and longitude without assumed position or plotting
- Observation of two celestial bodies places observer at intersection of two constant altitude circles
- Modern pocket calculators make direct computation practical (<30 seconds)

#### Content for Literature Review:
- **Constant altitude circle principle:**
  - Observation of altitude h places observer on small circle centered at GP
  - Circle radius = zenith distance = 90° - h (complement of altitude)
  - Two circles intersect at two points; one is fix, other is extraneous

- **Direct fix geometry (three triangle approach):**
  1. **Polar Triangle:** Formed by North Pole and two GPs; yields pseudoaltitude h12
  2. **Zenith Triangles:** Formed by two GPs and each intersection point; yields angle B
  3. **Navigational Triangle:** "Side-angle-side" known; yields latitude and meridian angle

#### Mathematical Models (Law of Cosines based):
**Altitude subroutine:**
```
alt(a,b,c) = arcsin[cos(a)·cos(b)·cos(c) + sin(b)·sin(c)]
```
h = alt(a,b,c) is complement of side opposite angle a

**Azimuth subroutine:**
```
azi(x,y,z) = arccos[(sin(x) - sin(y)·sin(z)) / (cos(y)·cos(z))]
```
Z = azi(x,y,z) is angle opposite side that is complement of x

**Core algorithm:**
```
1. t12 = GHA2 - GHA1                    # meridian angle in polar triangle
2. h12 = alt(t12, dec1, dec2)           # pseudoaltitude
3. A = azi(dec2, dec1, h12)             # angle at GP1 in polar triangle
4. B = azi(h2, h1, h12)                 # angle at GP1 in zenith triangle
5. P1 = A - B                           # parallactic angle (upper)
6. L1 = alt(P1, dec1, h1)               # latitude (upper intersection)
7. t1 = azi(h1, dec1, L1)               # meridian angle
8. Lo1 = GHA1 ± t1                       # longitude (sign from bearing)
```

**Latitude formula (law of cosines for altitude):**
```
sin(L) = cos(P)·cos(dec)·cos(h) - sin(dec)·sin(h)
```

**Meridian angle formula (time-sight form):**
```
cos(t) = [sin(h) - sin(L)·sin(dec)] / [cos(L)·cos(dec)]
```

#### Accuracy Analysis:
- 8 decimal places in trigonometric computations provide acceptable results
- Error at 7 decimal places: ~0.06 arcmin (~0.06 nmi)
- Error at 8 decimal places: ~0.002 arcmin
- Typical pocket calculator (HP-15C): 12 decimal places, "quite adequate"

#### Worked Examples from Paper:
**Example 1 (Artificial):**
- Bodies: dec1=75°, GHA1=30°, h1=60°; dec2=30°, GHA2=320°, h2=45°
- Result: L1=68.527°N, Lo1=80.291°E

**Example 2 (Real observations, Feb 4, 1995):**
- Bodies: Spica (h=47°33.8', dec=-11°08.2') and Venus (h=28°54.8', dec=-20°47.7')
- Result: 24°35.6'N, 81°46.4'W
- Matches Marcq St-Hilaire method to nearest 0.1'

#### Calculator Implementation:
- Hewlett-Packard 15C and Radio Shack EC-4026 programmed
- Running time: <30 seconds (exclusive of data entry)
- Complete code provided in appendix (easily adaptable to Python)

#### Key Quotes:
- "The intersections of these two constant altitude circles are the two possible points from which the pair of observed altitudes could have been measured"
- "The compact, low-power, lightweight, low-cost, fast, accurate... electronic calculator makes this direct fix approach practical"
- "eliminating the need for tables; and eliminating the need for plotting"

#### Relevance to Research:
- **PRIMARY SOURCE** for Python algorithm implementation
- Complete step-by-step algorithm with subroutines
- Law of cosines formulas directly translatable to code
- Precision analysis guides numerical implementation
- Multiple worked examples for validation testing
- Explicit discussion of calculator/computer implementation

---

### Article #10
**Citation:** Metcalf, T. R. (1991). Advancing Celestial Circles of Position. *NAVIGATION: Journal of The Institute of Navigation, 38*(3), 285–288. DOI: 10.1002/j.2161-4296.1991.tb01861.x

**Relevance:** Rigorous equations for running fix computation with vessel motion correction

#### Content for Introduction:
- Running fix requires advancing earlier lines of position for vessel motion
- Graphical methods approximate; computer methods need rigorous formulation
- Correcting GHA and declination maintains exact circle of position geometry

#### Content for Literature Review:
- **Problem with approximate methods:**
  - Simple altitude correction may alter LOP direction (especially near zenith)
  - Errors increase with distance traveled and high-altitude observations
  - Example: 600 nmi travel → 21.2 nmi GP error, 27.2° azimuth error

- **Rigorous approach:**
  - Rotate coordinate system rather than move GP along vessel course
  - Preserves relationship between GP and vessel position
  - "The circle of position is 'dragged' along with the vessel"

#### Mathematical Models:
**Vector rotation formula (Rodrigues' rotation):**
```
P' = P·cos(θ) + n̂(n̂·P)[1 - cos(θ)] + (P × n̂)·sin(θ)
```
Where:
- P = initial GP vector, P' = rotated GP vector
- n̂ = unit rotation axis vector
- θ = rotation angle (latitude or longitude change)

**Rhumb line motion (spherical Earth):**
```
θ_lat = D · cos(α)                                    # latitude change
θ_long = -ln[tan(45° + (L + D·cos(α))/2) / tan(45° + L/2)] × tan(α)   # longitude change
```
Where:
- D = distance traveled (degrees = nmi/60)
- α = true course (degrees clockwise from north)
- L = initial latitude

**Special case (course near E-W):**
```
θ_long = -D·sin(α) / cos(L)    when |D·cos(α)| << |cot(L)|
```

**Coordinate vectors:**
```
Rotation axis: n̂ = [cos(90° + λ + θ_long), sin(90° + λ + θ_long), 0]
GP vector:     P  = [cos(d)·cos(GHA + θ_long), cos(d)·sin(GHA + θ_long), sin(d)]
```

**Extract new GHA' and declination d':**
```
d' = arcsin(P'_3)
GHA' = arctan(P'_2 / P'_1) - θ_long
```

#### Error Analysis from Paper:
| Distance Traveled | GP Error (approx method) | Azimuth Error (altitude correction) |
|-------------------|--------------------------|--------------------------------------|
| 60 nmi | 0.06 nmi | -7.1° |
| 120 nmi | 0.53 nmi | -12.3° |
| 300 nmi | 4.60 nmi | -21.6° |
| 600 nmi | 21.20 nmi | -27.2° |

*Test case: Body at 85°17' altitude (near zenith), Course 045°*

#### Key Quotes:
- "By correcting the GHA and declination rather than the altitude, the circle of position is 'dragged' along with the vessel, making the correction exact"
- "a rigorous formulation is not much more demanding computationally and, with the advent of the handheld computer, may be of interest"

#### Relevance to Research:
- Essential for multi-sight running fix implementation in Python
- Vector rotation approach ideal for NumPy/SciPy implementation
- Error analysis provides validation benchmarks
- Addresses practical navigation scenario (non-simultaneous observations)
- Complements Gery (1997) two-body fix algorithm

---

### Article #11
**Citation:** Hoover, W. E. (1984). *Algorithms for Confidence Circles and Ellipses.* NOAA Technical Report NOS 107 C&GS 3. National Oceanic and Atmospheric Administration.

**Relevance:** Error quantification and uncertainty analysis for position fixing from intersecting Lines of Position

#### Content for Introduction:
- Position location traditionally determined by intersection of two Lines of Position (LOPs)
- LOPs may be derived from celestial observations, trilateration, LORAN signals, satellite signals
- Two fundamental questions: (1) probability of true position within given radius, (2) radius for given probability

#### Content for Literature Review:
- **Confidence circles vs. confidence ellipses:**
  - Confidence ellipse provides same probability over significantly smaller area
  - Error ellipse: center at observed position, oriented by angle θ from L₁
  - 50% ellipse: scale factor k ≈ 1.1774
  - 95% ellipse: scale factor k ≈ 2.4477

- **Error model assumptions (Section 2.2):**
  1. Earth flat in local region, LOPs are straight lines
  2. Measurement errors normally distributed with correlation ρ₁₂, zero means, standard deviations σ₁, σ₂
  3. Bivariate error distribution constant throughout region

- **Nonorthogonal to orthogonal transformation:**
  - u₁ perpendicular to L₁, u₂ perpendicular to L₂
  - Crossing angle α (0 < α < π)
  - Transform to x-y Cartesian system with orientation θ

#### Mathematical Formulas for Python Implementation:

**Error Ellipse Parameters (Algorithm 1, Steps 1-2):**
```
a₁ = σ₁²sin(2α) + 2ρ₁₂σ₁σ₂sin(α)
a₂ = σ₁² - σ₂² + 2ρ₁₂σ₁σ₂cos(α)
a₃ = σ₁² + σ₂² - 2ρ₁₂σ₁σ₂cos(α)
a₄ = 2sin²(α)

σₓ = [(a₃ + √(a₁² + a₂²)) / a₄]^(1/2)
σᵧ = [(a₃ - √(a₁² + a₂²)) / a₄]^(1/2)
θ = ½ arctan(a₁/a₂)    [use atan2 for proper quadrant]
```

**Confidence Ellipse Scale Factor:**
```
p = 1 - e^(-½k²)           [probability for given k]
k = √[-2·ln(1-p)]          [scale factor for given probability]
Semimajor axis: k·σₓ
Semiminor axis: k·σᵧ
```

**Confidence Circle Probability (Algorithm 1, Steps 3-8):**
```
c = σᵧ/σₓ
K = R/σₓ
δ = 2c/π
γ = (K/2c)²

w(φ) = (c² - 1)cos(φ) - (c² + 1)
f(φ) = [e^(γ·w(φ)) - 1] / w(φ)

p(R) = δ·∫₀^π f(φ)dφ    [numerical integration]
```

**Trapezoidal Rule with Simpson Correction:**
```
h = π/n
T₁ = f(0) + f(π)
T₂ = Σᵢ₌₁ⁿ⁻¹ f(i·h)
T₃ = Σᵢ₌₁ⁿ f((i-½)h)

p = ¼[T₁ + T₂ + T₃]·δ/n     [derivative-corrected Simpson's rule]
```

**Secant Method for R(p) (Algorithm 2):**
```
K_{i+1} = K_i - g_i · (K_i - K_{i-1}) / (g_i - g_{i-1})
where g_i = p_i - p (target probability)

R = σₓ · K    [final radius]
```

#### Special Cases (Section 5.1):
**When σ₁ = σ₂ = σ and ρ₁₂ = 0:**
```
σₓ = σ·csc(α)·[1 + |cos(α)|]^(1/2)
σᵧ = σ·csc(α)·[1 - |cos(α)|]^(1/2)
θ = α/2
```

**Circular normal (σ₁ = σ₂ = σ, ρ₁₂ = 0, α = π/2):**
```
σₓ = σᵧ = σ
p(R) = 1 - e^(-(R/σ)²/2)
R(p) = σ·√[-2·ln(1-p)]
```

#### Numerical Examples from Paper:
| Parameters | σₓ (nm) | σᵧ (nm) | θ | 95% Circle R | 95% Ellipse Area |
|------------|---------|---------|-------|--------------|------------------|
| σ₁=2, σ₂=1, α=30° | 4.3778 | 0.9137 | 24.56° | 8.6302 nm | 75.3 nm² |
| σ₁=σ₂=1, α=30° | 2.7321 | 0.7321 | 15° | 5.4069 | 37.6 |
| σ₁=σ₂=1, α=90° | 1.0000 | 1.0000 | 45° | 2.4477 | 18.8 |

**Key Finding:** 95% confidence circle area can be 200%+ larger than 95% ellipse for same probability.

#### Key Terms:
- **CEP (Circular Error Probable):** Radius of 50% confidence circle
- **1dRMS (Radial error):** √(σₓ² + σᵧ²)
- **2dRMS:** Upper bound for 95% confidence circle radius

#### Key Quotes:
- "The LOPs may be derived from celestial observations, trilateration, LORAN signals, satellite signals, etc."
- "Confidence ellipses should be considered as a superior alternative...since much less computation is required"
- "The area of a confidence ellipse is generally substantially less than the area of a confidence circle having the same associated probability"

#### Relevance to Research:
- Essential for SRQ4 (validation framework) and H4 (accuracy benchmarks)
- Provides complete algorithms for quantifying position fix uncertainty
- Implements error propagation from LOP measurements to position uncertainty
- Circular error probability table (Appendix) useful for validation
- Numerical quadrature methods directly implementable in Python (NumPy/SciPy)
- Error model applicable to celestial navigation two-body fix

---

### Article #12
**Citation:** Nguyen, V. S., Im, N. K., & Dao, Q. D. (2017). Azimuth method for ship position in celestial navigation. *International Journal of e-Navigation and Maritime Economy, 7*, 55–62.

**Relevance:** Alternative celestial navigation method using azimuth instead of altitude - works when horizon invisible

#### Content for Introduction:
- Traditional celestial methods require both celestial bodies AND visible horizon (twilight only)
- Azimuth method overcomes limitation: can determine position at night when horizon invisible
- Does not require sextant equipment
- Uses great-circle azimuth of observed celestial bodies

#### Content for Literature Review:
- **Problem with altitude-based methods:**
  - Only applicable during twilight when stars and horizon visible simultaneously
  - Cannot be used at night even when stars visible
  - Requires sextant for altitude measurement

- **Azimuth method concept:**
  - Uses compass bearing to celestial bodies instead of altitude
  - Great circle from observer through star defines Circle of Position
  - Light travels shortest arc (great circle) from star to observer
  - Azimuth angle formed by meridian arc and great circle through star

#### Mathematical Formulas for Python Implementation:

**Cartesian Coordinate System:**
```
Ship position vector: OP = Xi + Yj + Zk
Celestial body vector: OC = x_C·i + y_C·j + z_C·k

Spherical to Cartesian:
[X, Y, Z] = [cos(φ)cos(λ), cos(φ)sin(λ), sin(φ)]
[x_C, y_C, z_C] = [cos(δ)cos(GHA), cos(δ)sin(GHA), sin(δ)]

Cartesian to Spherical:
φ = atan(Z / √(X² + Y²))
λ = atan(Y / X)
```

**Great-Circle Azimuth Equation:**
```
cosA = [(z_C·Y - y_C·Z)Y - (x_C·Z - z_C·X)X] / 
       [√((z_C·Y - y_C·Z)² + (x_C·Z - z_C·X)²) · √((1-Z²) + (y_C·X - x_C·Y)²)]
```

**Nonlinear System for Ship Position:**
```
f_i(X,Y,Z) = [(z_Ci·Y - y_Ci·Z)Y - (x_Ci·Z - z_Ci·X)X] /
             [√((z_Ci·Y - y_Ci·Z)² + (x_Ci·Z - z_Ci·X)²) · √((1-Z²) + (y_Ci·X - x_Ci·Y)²)]
             - cos(A_i) = 0

For 3 stars: F(P) = [f₁(X,Y,Z), f₂(X,Y,Z), f₃(X,Y,Z)]ᵀ = 0
```

**Newton-Raphson Iteration:**
```
P(k+1) = P(k) - J⁻¹(P(k)) · F(P(k))

where J is the Jacobian matrix:
J_ij = ∂f_i/∂X_j  for X_j ∈ {X, Y, Z}

Convergence criterion: ||F(P(k))|| ≈ 0 or ||P(k+1) - P(k)||² ≤ ε
```

#### Numerical Experiment Results:
| Parameter | Value |
|-----------|-------|
| Test location | Gulf of Tonkin, Vietnam |
| Stars observed | Antares, Regulus, Vega |
| GMT | 13h 35m 28s (June 9, 2016) |
| Estimated position | 20°28.7'N, 107°13'E |
| True position | 20°37.4'N, 107°20.1'E |
| Calculated position | 20°35.6'N, 107°17.8'E |
| **Position error** | **2.92 NM** |

#### Key Quotes:
- "The methods of celestial navigation to fix the ship position in line with the stars are only applied in the twilight time interval when both the celestial bodies and the horizon appear simultaneously"
- "The key advantage which differentiates this method from previous ones is its ability to determine the ship position during the night when the horizon is invisible"
- "the proposed method does not demand the horizon and sextant equipment as with the previous methods"

#### Relevance to Research:
- Provides alternative method for position fixing without altitude measurement
- Newton-Raphson solver directly implementable in Python (SciPy)
- Vector calculus approach aligns with modern computational methods
- Complements altitude-based methods for complete celestial navigation system
- Real-world validation with 2.92 NM accuracy benchmark
- Addresses limitation of traditional methods (twilight-only observation)

---

### Article #13
**Citation:** Tsou, M. C. (2012). Genetic algorithm for solving celestial navigation fix problems. *Polish Maritime Research, 19*(3), 53–59.

**Relevance:** Alternative optimization approach using evolutionary computation for celestial fix

#### Content for Introduction:
- Celestial navigation remains important backup when GPS unavailable
- Traditional intercept method limited: ~20 minutes per fix, 30 NM assumed position constraint
- IMO Manila 2010 amendment encourages electronic nautical almanacs and calculation software
- e-Navigation development direction for celestial navigation

#### Content for Literature Review:
- **Limitations of intercept method:**
  1. Assumed position must be within 30 NM of actual position
  2. Celestial body altitude should not exceed 70° (curvature error)
  3. Approximate method with graphical procedures
  4. Takes ~20 minutes for professional seafarers

- **Limitations of vector-matrix methods:**
  - Require prior knowledge of observer's position
  - Gradient descent can converge to local minimum
  - Least squares requires >3 celestial bodies

- **Genetic Algorithm advantages:**
  1. Simple calculation without complicated mathematical procedures
  2. Does not converge to local optima (unlike gradient methods)
  3. More flexibility in practice
  4. Rapid convergence toward final result

#### Mathematical Formulas for Python Implementation:

**Calculated Altitude (Hc):**
```
Hc = sin⁻¹(sin(L)·sin(Dec) + cos(L)·cos(Dec)·cos(GHA - λ))

where:
  L = observer's latitude
  λ = observer's longitude
  Dec = celestial body declination
  GHA = Greenwich Hour Angle
```

**Fitness Function (RMSE):**
```
RMSE = √[Σ(Ho_i - Hc_i)² / n]

where:
  Ho_i = observed altitude (corrected for dip, refraction, etc.)
  Hc_i = calculated altitude from equation above
  n = number of observations
```

**Initial Population Generation:**
```
L_individual = IniL + Tol × RandNum
λ_individual = Iniλ + Tol × RandNum
P₁ = {(L₁,λ₁), (L₂,λ₂), ..., (L_n,λ_n)}

where:
  IniL, Iniλ = reference position (latitude, longitude)
  Tol = maximum search degree (default: 2°)
  RandNum = random number in [-1, +1]
  n = population size (default: 50)
```

**Crossover Operator:**
```
L_new = L₁ + (L₂ - L₁) × RandNum
λ_new = λ₁ + (λ₂ - λ₁) × RandNum

Crossover rate P_c = 0.6
```

**Mutation Operator:**
```
L_mutation = L_old + Noise × RandNum
λ_mutation = λ_old + Noise × RandNum

Noise = perturbation degree (default: 1°)
Mutation rate P_m = 0.05
```

#### GA Parameters (Recommended):
| Parameter | Value |
|-----------|-------|
| Population size | 50 individuals |
| Crossover rate (P_c) | 0.6 |
| Mutation rate (P_m) | 0.05 |
| Maximum search space | 4° × 4° (lat × lon) |
| Perturbation range | 1° |
| Convergence | ~50 generations, ~5 seconds |

#### Validation Results:

**Case 1 - Two-Body Fix:**
| Body | Ho | GHA | Dec |
|------|-----|-----|-----|
| Capella | 15°19.3' | 131°24.8'E | 45°58.4'N |
| Alkaid | 77°34.9' | 003°14.2'E | 49°25.7'N |

| Method | Latitude | Longitude |
|--------|----------|----------|
| Intercept | 41°38.6'N | 017°08.1'W |
| SEEM | 41°39.1'N | 017°07.3'W |
| **GA** | **41°39.1'N** | **017°07.3'W** |

*Note: Alkaid altitude 77°34.9' exceeds 70° limit of intercept method*

**Case 2 - Multi-Body Fix (4 stars):**
| Method | Latitude | Longitude |
|--------|----------|----------|
| Intercept | 35°19.0'S | 005°26.5'E |
| Vector-Matrix | 35°18.6'S | 005°26.9'E |
| **GA** | **35°18.6'S** | **005°27.0'E** |

#### Key Quotes:
- "celestial fix still serves as an important backup measure"
- "traditional methods for computing a celestial navigation fix can no longer meet the requirements of modern vessels in terms of calculation speed and precision"
- "GA shows excellent performance in optimum fix processing problems"
- "does not suffer from possible sensitivity to the initial position chosen to start the iterative process"

#### Relevance to Research:
- Provides alternative optimization method to Newton-Raphson/least squares
- Avoids local minima convergence issues of gradient-based methods
- Removes 30 NM assumed position constraint
- Complete parameter set for Python implementation (DEAP library)
- Handles high-altitude observations (>70°) that cause intercept method errors
- ~5 second convergence suitable for real-time navigation
- Directly addresses research objectives of accurate, computational celestial fix

---

### Article #14
**Citation:** Silverberg, J. (2007). Circles of Illumination, Parallels of Equal Altitude, and le Calcul du Point Observé: Nineteenth Century Advances in Celestial Navigation. *Proceedings of the Canadian Society for the History and Philosophy of Mathematics*, presented at Concordia University, Montréal.

**Relevance:** Historical foundations and mathematical derivation of Sumner line and intercept method

#### Content for Introduction:
- By early 19th century, latitude/longitude determination at sea was possible but flawed
- Problem: Longitude calculation depends on latitude, which changes during voyage
- Dead reckoning updates latitude but inherently inaccurate
- Two key innovations: Sumner's Line of Position (1837) and Saint-Hilaire's Intercept Method (1873)

#### Content for Literature Review:
- **Limitations of pre-Sumner methods:**
  - Latitude determined at noon (meridian observation)
  - Longitude determined when sun bears E/W (time sight)
  - Ship may move hundreds of miles between observations
  - Dead reckoning latitude used for longitude calculation introduces errors
  - Errors of 1°-2° longitude possible (dangerous near coast)

- **Sumner's Discovery (December 1837):**
  - Approaching Welsh coast in gale without adequate position
  - Calculated longitude for three different assumed latitudes
  - Found all three positions lay on a line pointing toward Small's Light
  - Realized: "the observed altitude must have happened at all the three points, and at Small's light, and at the ship, at the same instant of time"
  - Published 1843: "New and Accurate Method of Finding a Ship's Position at Sea"

- **Sumner's concepts:**
  - Circle of Illumination = boundary between illuminated/dark hemispheres
  - Parallels of Equal Altitude = small circles where sun has same altitude
  - Pole of Illumination = geographic position (GP) where sun is at zenith
  - Line of Position = tangent to circle of equal altitude (appears straight on chart)

- **Saint-Hilaire's Intercept Method (1873-1875):**
  - Published in Revue Maritime et Coloniale
  - Uses estimated latitude AND longitude (not just latitude)
  - Calculates altitude and zenith angle for assumed position
  - Compares calculated altitude with observed altitude
  - Moves position toward/away from sun's GP by altitude difference
  - Constructs line perpendicular to azimuth through corrected position

#### Mathematical Formulas for Python Implementation:

**Sumner's Calculation (Time Sight Formula):**
```
cos(t) = [sin(h) - sin(L)·sin(δ)] / [cos(L)·cos(δ)]

Alternative form for logarithmic calculation:
1 - cos(t) = [cos(L-δ) - sin(h)] · sec(L) · sec(δ)

Bowditch/Saint-Hilaire form:
1 - cos(t) = 2·cos(s-δ)·cos(s-L)·sec(L)·sec(δ)
where s = (L + δ + h) / 2

t = hour angle (difference in longitude between GP and observer)
h = observed altitude
L = latitude
δ = declination
```

**Saint-Hilaire's Intercept Method:**
```
Calculated Altitude (Hc):
sin(Hc) = sin(δ)·sin(L) + cos(δ)·cos(L)·cos(t)

Zenith Angle (Zn):
sin(Zn) = cos(δ)·sec(Hc)·sin(t)

Intercept:
a = Ho - Hc  (in minutes of arc = nautical miles)

If a > 0: move TOWARD the sun's GP
If a < 0: move AWAY from sun's GP

Line of Position: perpendicular to azimuth at intercept distance
```

**French Practice (Numerical Position):**
```
ΔLat = a · cos(Zn)
Departure = a · sin(Zn)
ΔLon = Departure / cos(Lat_new)

Lat_new = Lat_assumed + ΔLat
Lon_new = Lon_assumed + ΔLon
```

**Sight Reduction Tables (20th century development):**
- Assumed latitude: rounded to whole degrees
- Assumed longitude: adjusted so LHA is whole degrees
- LHA = GHA - Assumed_Longitude
- Tables provide Hc and Zn for given (LHA, Dec, Lat)

#### Historical Timeline:
| Year | Development | Contributor |
|------|-------------|-------------|
| 1837 | Line of Position discovered | Thomas Sumner |
| 1843 | Sumner's method published | Thomas Sumner |
| 1873 | "Note sur la Détermination du Point" | Marcq Saint-Hilaire |
| 1875 | "Calcul du Point Observé" (2 parts) | Marcq Saint-Hilaire |
| 1876 | Method published by Thomson | Sir William Thomson |
| 1920-1940 | Sight Reduction Tables developed | Various |

#### Key Quotes:
- Sumner: "the observed altitude must have happened at all the three points, and at Small's light, and at the ship, at the same instant of time"
- "Sir W. Thomson...a déclaré que ce serait la plus grande bénédiction pour les navigateurs jeunes et vieux, si l'on supprimait toutes les méthods autres que celles de Sumner" [would be greatest blessing to navigators if all methods except Sumner's were suppressed]
- "Modern calculators would not make this compromise, and in fact make use of tables in the first place completely unnecessary"

#### Relevance to Research:
- Essential historical context for Introduction section
- Provides original mathematical derivations of intercept method
- Explains evolution from graphical to computational methods
- Formulas directly implementable in Python
- Shows how modern computers eliminate need for sight reduction tables
- Documents the transition from Sumner's method to Saint-Hilaire's intercept method
- Supports SRQ1 (appropriate algorithms) with historical algorithm development

---

### Article #15
**Citation:** Nguyen, V.-S., & Im, N.-K. (2014). Development of Computer Program for Solving Astronomical Ship Position Based on Circle of Equal Altitude Equation and SVD-Least Square Algorithm. *Journal of Navigation and Port Research*, 38(2), 89-96. DOI: 10.5394/KINPR.2014.38.2.89

**Relevance:** Complete matrix-based algorithm for celestial fix using SVD decomposition - directly implementable in Python/NumPy

#### Content for Introduction:
- Traditional intercept method (LOP) is approximation of circle of equal altitude (COP)
- LOP introduces error because tangent line approximates curved COP
- Modern computational approach can solve COP equations directly
- SVD provides numerical stability superior to normal decomposition
- STCW 2010 Manila amendments encourage software-based celestial navigation

#### Content for Literature Review:
- **Limitations of intercept method identified:**
  - Co-altitude constraint: must not exceed 25 minutes
  - Distance constraint: ship must be within 30 NM of assumed position
  - LOP is tangent line approximation → inherent accuracy loss
  - Graphical plotting introduces additional errors

- **Circle of Equal Altitude (COP) advantages:**
  - COP represents exact locus of ship position at observation time
  - Direct solution without graphical procedure
  - No assumed position distance constraint
  - Higher accuracy than tangent line approximation

- **SVD-Least Squares advantages over normal decomposition:**
  - More numerically stable
  - Handles rank-deficient matrices
  - Avoids problems when celestial matrix is poorly conditioned
  - Better solution for overdetermined systems

#### Mathematical Formulas for Python Implementation:

**COP Equation (Spherical Coordinates):**
```
sin(δ)·sin(φ) + cos(δ)·cos(φ)·cos(t) = sin(h)

Where:
δ = declination of celestial body
φ = observer's latitude
t = local hour angle (LHA)
h = observed altitude
```

**COP Equation (Cartesian Coordinates) - Key Innovation:**
```
Xx + Yy + Zz = sin(h)      [Equation 8]

Where celestial body coordinates:
X = cos(δ)·cos(GHA)
Y = cos(δ)·sin(GHA)
Z = sin(δ)

Ship position coordinates (on unit sphere):
x = cos(φ)·cos(λ)
y = cos(φ)·sin(λ)
z = sin(φ)

Note: Celestial sphere radius R = 1 (normalized)
```

**Coordinate System Rules:**
```
Longitude (λ): East positive, West negative
                -180° ≤ λ ≤ 180°

GHA: 0° ≤ GHA ≤ 360° (measured westward from Greenwich)

Latitude (φ): -90° ≤ φ ≤ 90°

Declination (δ): -90° ≤ δ ≤ 90°
```

**Equation System (Matrix Form):**
```
AX = b + ε      [Equation 13]

Where:
A = | X₁  Y₁  Z₁ |    (celestial body matrix, n×3)
    | X₂  Y₂  Z₂ |
    | ..  ..  .. |
    | Xₙ  Yₙ  Zₙ |

X = | x |    (ship position vector, 3×1)
    | y |
    | z |

b = | sin(h₁) |    (observed altitudes, n×1)
    | sin(h₂) |
    | ..      |
    | sin(hₙ) |

ε = error vector
```

**Least Squares Solution:**
```
Minimize: ||AX - b||²

Normal equation solution:
X = (AᵀA)⁻¹Aᵀb      [Equation 16]
```

**SVD Decomposition - More Stable:**
```
A = U × D × Vᵀ      [Equation 17]

Where:
U = m×m column orthonormal matrix (UᵀU = I)
    columns are eigenvectors of AAᵀ

D = m×n diagonal matrix of singular values
    σ₁ ≥ σ₂ ≥ ... ≥ σₙ ≥ 0
    σᵢ = √(eigenvalues of AᵀA)

V = n×n orthonormal matrix (VᵀV = I)
    columns are eigenvectors of AᵀA
```

**SVD Solution for Ship Position:**
```
X = V × D⁻¹ × Uᵀ × b      [Equation 24]

Python/NumPy implementation:
import numpy as np
from numpy.linalg import svd

U, s, Vt = svd(A, full_matrices=False)
X = Vt.T @ np.diag(1/s) @ U.T @ b

# Or use numpy.linalg.lstsq which uses SVD internally
X, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
```

**Virtual Celestial Body (for 2-body observations):**
```
Problem: Need ≥3 equations for 3 unknowns (x, y, z)

Solution: Use assumed position as "virtual body"

Virtual body altitude:
sin(hᵥ) = (S·t) / (111.11·cos(φ))·sin(φ)      [Equation 10]

Where:
S = ship speed (knots)
t = time between observations
φ = assumed latitude
111.11 = conversion factor (1° latitude ≈ 60 NM)

Additional equation:
Xᵥx + Yᵥy + Zᵥz = sin(hᵥ)      [Equation 11]

Where (Xᵥ, Yᵥ, Zᵥ) are Cartesian coordinates of assumed position
```

**Conversion to Geodetic Coordinates (WGS-84):**
```
From Cartesian (x, y, z) to Lat/Lon:

φ = arcsin(z)           [latitude]
λ = arctan2(y, x)       [longitude]

WGS-84 parameters [Equation 25]:
a = 6,378,137 m  (equatorial radius)
b = 6,356,752.3 m  (polar radius)

Note: For celestial navigation, unit sphere conversion:
φ = arcsin(z / √(x² + y² + z²))
λ = arctan2(y, x)
```

#### Experimental Results (Accuracy Comparison):

| Observation | Intercept Method | SVD-Least Squares | Improvement |
|-------------|------------------|-------------------|-------------|
| Sun (2 obs) | 4.45 NM | 1.94 NM | 56% |
| 2 Stars | 3.42 NM | 2.46 NM | 28% |
| 3 Stars | 1.649 NM | 1.146 NM | 31% |

**Test Conditions:**
- True position verified on chart
- Stars used: Rasalhague, Vega, Deneb
- Location: ~20°39'N, 107°E (Gulf of Tonkin region)
- Ship speed: 12 knots, Course: 98°-112°

#### Algorithm Flow (for Python implementation):
```
1. Input: Sextant altitudes at times t₁, t₂, ... tₙ
         Altitude correction parameters (IE, dip, refraction)
         Assumed position (for virtual body if n < 3)

2. Calculate observed altitudes:
   Ho = Hs + corrections

3. Calculate celestial body Cartesian coordinates:
   For each body at observation time:
   - Get GHA, declination from almanac
   - X = cos(δ)·cos(GHA)
   - Y = cos(δ)·sin(GHA)
   - Z = sin(δ)

4. Build equation system:
   A = [X₁ Y₁ Z₁; X₂ Y₂ Z₂; ...]
   b = [sin(Ho₁); sin(Ho₂); ...]

5. If n < 3: Add virtual body equation

6. Solve using SVD:
   X = V × D⁻¹ × Uᵀ × b

7. Convert to spherical coordinates:
   φ = arcsin(z)
   λ = arctan2(y, x)

8. Output: Latitude, Longitude
```

#### Key Quotes:
- "the ship position given LOP has low accuracy" (comparing LOP tangent approximation to true COP)
- "SVD was proved more numerically stable than normal decomposition"
- "the accuracy of ship position based on new method is better than the intercept method"
- "if the altitude of celestial bodies are determined continuously by electrical equipments, this method has the potential to be integrated to electrical chart and the e-navigation system"

#### Relevance to Research:
- **Directly addresses SRQ1:** Complete algorithm for celestial position fix
- **Matrix formulation ready for NumPy:** A, b matrices directly implementable
- **SVD via scipy.linalg.svd or numpy.linalg.lstsq**
- **Quantified accuracy improvements:** 28-56% better than intercept method
- **Virtual body technique:** Solves 2-observation problem elegantly
- **WGS-84 conversion included**
- **Validates hypothesis H1:** Computational approach matches/exceeds traditional accuracy
- Comparison baseline for evaluating proposed algorithm
- Flow diagram directly translatable to Python code

---

### Article #16
**Citation:** Chiesa, A., & Chiesa, R. (1990). A Mathematical Method of obtaining an Astronomical Vessel Position. *The Journal of Navigation*, 43(1), 125-129. DOI: 10.1017/S0373463300013862

**Relevance:** Direct analytical solution for COP intersection - no estimated position required; complete spherical trigonometry formulas for Python implementation

#### Content for Introduction:
- Traditional intercept method uses tangent lines (LOPs) to approximate position circles
- Tangent line intersection does NOT exactly coincide with true COP intersection
- With programmable calculators, no reason to use approximation instead of exact solution
- Key advantage: No estimated/assumed position required to start calculation
- Method proven on Atlantic sailing yacht passages

#### Content for Literature Review:
- **Limitations of LOP method identified:**
  - LOPs are tangent lines to position circles
  - Intersection of tangent lines ≠ intersection of circles on Earth's surface
  - Requires estimated position to start
  - Uses tabular methods (HO 214, HO 249, HO 229) designed for manual computation

- **Direct COP intersection advantages:**
  - Exact vessel position by mathematical definition
  - No estimated position needed
  - Simply enter sextant readings → get position coordinates
  - Works with any celestial body type
  - Supports running fix procedure

- **Two-body observation reality:**
  - Two circles always intersect at two points (physical reality)
  - Both points mathematically valid
  - Points typically thousands of miles apart → easy to choose correct one
  - Example: Atlantic fix gave P₁ near Newfoundland, P₂ near Madagascar (6500+ miles apart)

#### Mathematical Formulas for Python Implementation:

**Input Data (for each observed body):**
```
GHA = Greenwich Hour Angle (from almanac)
Dec = Declination (from almanac)
H = True altitude (corrected sextant altitude)
ZD = 90° - H  (zenith distance = radius of position circle)

Geographical Position (GP) on Earth:
φ_GP = Dec      (latitude of GP)
λ_GP = GHA      (longitude of GP, measured westward)
```

**Step 1: Orthodromic Distance and Course Between Two GPs:**
```
Given two GPs: S₁(φ₁, λ₁) and S₂(φ₂, λ₂)

Initial course angle from S₁ to S₂:
tan(R) = sin(λ₂ - λ₁) / [tan(φ₂)·cos(φ₁) - sin(φ₁)·cos(λ₂ - λ₁)]

Orthodromic (great circle) distance:
cos(D) = sin(φ₁)·sin(φ₂) + cos(φ₁)·cos(φ₂)·cos(λ₂ - λ₁)

Where:
φ₁ = Dec₁, φ₂ = Dec₂  (declinations)
λ₁ = GHA₁, λ₂ = GHA₂  (hour angles)
```

**Step 2: Spherical Triangle Solution (Half-Angle Formula):**
```
Known sides of spherical triangle S₁P₁S₂:
- Side S₁S₂ = D (distance between GPs)
- Side S₁P = ZD₁ = 90° - H₁
- Side S₂P = ZD₂ = 90° - H₂

Angle α at vertex S₁:
sin(α/2) = √[(cos(H₁ + H₂ - D) - cos(H₁ - H₂ + D)) / (2·sin(D)·cos(H₁))]

Simplified form using half-angle identities:
sin(α/2) = √[sin((ZD₁ + ZD₂ - D)/2)·sin((ZD₂ - ZD₁ + D)/2)] / (sin(D)·cos(H₁))
```

**Step 3: Course Angles to Two Intersection Points:**
```
R₁ = R - α    (course from S₁ to first intersection P₁)
R₂ = R + α    (course from S₁ to second intersection P₂)
```

**Step 4: Position Coordinates of Intersection Points:**
```
Given:
- Starting point S₁ with coordinates (φ₁, λ₁)
- Initial course angle R₁ or R₂
- Distance ZD₁ = 90° - H₁

Latitude of intersection point P:
φ_P = λ₁ + arctan[tan(R/2) · cos((ZD₁ - 90 + φ₁)/2) / cos((ZD₁ + 90 - φ₁)/2)]
      + arctan[tan(R/2) · sin((ZD₁ - 90 + φ₁)/2) / sin((ZD₁ + 90 - φ₁)/2)]

Longitude of intersection point P:
λ_P = arccos[sin(R)·sin(ZD₁) / sin(φ_P - λ₁)]

(Apply for both R₁ and R₂ to get P₁ and P₂)
```

**Python Implementation Skeleton:**
```python
import numpy as np

def two_body_fix(gha1, dec1, h1, gha2, dec2, h2):
    """
    Calculate vessel position from two celestial body observations.
    Returns two possible positions (choose based on DR position).
    
    Parameters:
    - gha1, dec1, h1: GHA, Dec, true altitude of body 1 (degrees)
    - gha2, dec2, h2: GHA, Dec, true altitude of body 2 (degrees)
    
    Returns:
    - (lat1, lon1), (lat2, lon2): Two possible vessel positions
    """
    # Convert to radians
    phi1, lam1 = np.radians(dec1), np.radians(gha1)
    phi2, lam2 = np.radians(dec2), np.radians(gha2)
    zd1, zd2 = np.radians(90 - h1), np.radians(90 - h2)
    
    # Step 1: Distance and course between GPs
    delta_lam = lam2 - lam1
    R = np.arctan2(np.sin(delta_lam),
                   np.tan(phi2)*np.cos(phi1) - np.sin(phi1)*np.cos(delta_lam))
    D = np.arccos(np.sin(phi1)*np.sin(phi2) + 
                  np.cos(phi1)*np.cos(phi2)*np.cos(delta_lam))
    
    # Step 2: Half-angle formula for alpha
    h1_rad, h2_rad = np.radians(h1), np.radians(h2)
    num = np.cos(h1_rad + h2_rad - D) - np.cos(h1_rad - h2_rad + D)
    sin_alpha_half = np.sqrt(num / (2 * np.sin(D) * np.cos(h1_rad)))
    alpha = 2 * np.arcsin(sin_alpha_half)
    
    # Step 3: Two course angles
    R1, R2 = R - alpha, R + alpha
    
    # Step 4: Calculate both positions
    pos1 = calculate_position(phi1, lam1, R1, zd1)
    pos2 = calculate_position(phi1, lam1, R2, zd1)
    
    return pos1, pos2
```

#### Three-Body Observation Algorithm:
```
With 3 celestial bodies:
- 3 circles of position
- 6 intersection points (each pair of circles intersects at 2 points)
- Only ONE point of "multiple coincidence" where all 3 circles meet

Due to observation errors:
- No exact coincidence point exists
- Instead, find triplet of points that are closest to each other
- Return average of triplet coordinates as vessel position

Error Detection:
- If no triplet within preset tolerance → display "ERROR"
- Indicates one or more observations have excessive error
- When position appears → observations valid within tolerance
```

#### Running Fix Procedure:
```
For observations at different times:
1. Observe celestial body at time t₁ → get GP₁
2. Observe same/different body at time t₂ → get GP₂
3. Transfer GP₁ according to course and distance sailed between t₁ and t₂
4. Calculate intersection using transferred GP₁ and GP₂

Advantage: Simply transfer GP coordinates before intersection calculation
(Equivalent to advancing the first LOP in traditional method)
```

#### Worked Example from Paper:
```
Date: 15 September 1988, 08:58:00 GMT
Location: Atlantic passage, near Newfoundland

Observations:
- Venus: H = 34°54.5'
- Sirius: H = 22°05.0'

Results (two mathematically valid positions):
- P₁: 46°33.6'N, 55°18.8'W  (Atlantic, ~60 miles S of Newfoundland) ← CORRECT
- P₂: 18°58.7'S, 43°56.7'E  (Indian Ocean, near Madagascar)

Distance between solutions: ~6500 nautical miles
→ No ambiguity in choosing correct position
```

#### Key Quotes:
- "a vessel's position is obtained according to its exact definition, as the intersection point of the position circles of the observed celestial bodies"
- "No introduction of an estimated position is required to start the calculation"
- "Simple geometrical considerations demonstrate that the intersection of such straight lines [LOPs] does not exactly coincide with the intersection point of the position circles"
- "the navigators of three sailing yachts have made extensive use of the method described, programmed on small minicalculators, during their Atlantic passages"

#### Relevance to Research:
- **Directly addresses SRQ1:** Alternative algorithm to intercept method
- **No estimated position requirement:** Major advantage for autonomous systems
- **Exact solution:** Position circles, not tangent line approximation
- **Spherical trigonometry formulas:** Complete and ready for Python
- **Three-body algorithm:** Handles practical observation errors
- **Running fix included:** Essential for single-body (Sun) observations
- **BASIC implementation proven:** Direct port to Python feasible
- **Validates H1:** Computational method matches/exceeds table-based accuracy
- Comparison baseline for evaluating SVD and other matrix methods
- Real-world validation on Atlantic yacht passages

---

### Article #17
**Citation:** Nguyen, V. S. (2019). A Theoretical Approach of Astronomical Ship Positioning Using a Single Celestial Body and Secant Technique. *International Journal of Civil Engineering and Technology*, 10(10), [pages not specified].

**Relevance:** Novel method requiring only ONE celestial body observation (altitude + azimuth) - solved via Secant numerical method; directly implementable in Python

#### Content for Introduction:
- Traditional celestial methods require either:
  - At least 2 celestial bodies, OR
  - 1 body observed at 2 different times
- This method: 1 body, 1 observation time, using BOTH altitude and azimuth simultaneously
- Useful when only one celestial body is identifiable (cloudy conditions, limited visibility)
- Secant method provides rapid convergence (3-4 iterations)

#### Content for Literature Review:
- **Problem addressed:**
  - Existing methods require 2+ observations or 2+ bodies
  - Time-consuming when many steps needed before fixing
  - What if only ONE celestial body is identifiable?

- **Novel contribution:**
  - Uses altitude (H) AND azimuth (A) from single body at single time
  - Both measurements provide independent constraints on position
  - Circle of equal altitude + great circle of observed azimuth intersection

- **Numerical technique:**
  - Secant method (root-finding for nonlinear equations)
  - Faster than bisection, doesn't require derivative like Newton-Raphson
  - Converges in 3-4 iterations for navigation problems

#### Mathematical Formulas for Python Implementation:

**Spherical Triangle Setup:**
```
Spherical triangle: P_N (North Pole) - P (Observer) - C (Celestial Body)

Known:
- GHA, δ (declination) of celestial body from almanac
- H (true altitude) from sextant
- A (azimuth) from compass (corrected for variation/deviation)

Unknown:
- φ (latitude of observer)
- λ (longitude of observer)
```

**Fundamental Relations:**
```
Local Hour Angle:
LHA = GHA - λ     (if observer West of Greenwich)
LHA = GHA + λ     (if observer East of Greenwich)
General: LHA = GHA ± λ

Rules:
- If LHA > 180° → use (LHA - 180°)
- If LHA > 360° → use (LHA - 360°)
```

**Equation 1 - Altitude Equation (Cosine Formula):**
```
sin(H) = sin(φ)·sin(δ) + cos(φ)·cos(δ)·cos(LHA)      [Eq. 5]

Solving for longitude:
λ = GHA - arccos[(sin(H) - sin(φ)·sin(δ)) / (cos(φ)·cos(δ))]     [Eq. 8]
```

**Equation 2 - Azimuth Equation (Four-Part Formula):**
```
cot(A)·sin(LHA) = cos(φ)·tan(δ) - sin(φ)·cos(LHA)      [Eq. 7]

Rearranged:
cos(φ)·tan(δ) - sin(φ)·cos(LHA) - cot(A)·sin(LHA) = 0
```

**Combined Function F(φ) to Solve:**
```
Substitute λ from Eq. 8 into Eq. 7:

F(φ) = cos(φ)·tan(δ) - sin(φ)·cos(LHA(φ)) - cot(A)·sin(LHA(φ)) = 0

where LHA(φ) = GHA - arccos[(sin(H) - sin(φ)·sin(δ)) / (cos(φ)·cos(δ))]

Objective: Find φ such that F(φ) = 0
```

**Secant Method Iteration:**
```
φ_{n+1} = φ_n - F(φ_n) · (φ_n - φ_{n-1}) / (F(φ_n) - F(φ_{n-1}))     [Eq. 11]

Convergence criterion: |φ_{n+1} - φ_n| < ε (typically ε = 1e-6)

Initial values:
φ_0 = DR_latitude
φ_1 = DR_latitude ± 1° (depending on course direction)
```

**Python Implementation:**
```python
import numpy as np
from scipy.optimize import brentq  # Alternative to manual Secant

def single_body_fix(gha, dec, H, A, dr_lat):
    """
    Calculate ship position from single celestial body observation.
    Uses both altitude (H) and azimuth (A) simultaneously.
    
    Parameters:
    - gha: Greenwich Hour Angle (degrees)
    - dec: Declination (degrees)
    - H: True altitude (degrees)
    - A: Azimuth (degrees, 0-360)
    - dr_lat: Dead reckoning latitude for initial guess (degrees)
    
    Returns:
    - (latitude, longitude) in degrees
    """
    gha_r = np.radians(gha)
    dec_r = np.radians(dec)
    H_r = np.radians(H)
    A_r = np.radians(A)
    
    def longitude_from_lat(phi):
        """Calculate longitude from latitude using altitude equation."""
        cos_lha = (np.sin(H_r) - np.sin(phi)*np.sin(dec_r)) / \
                  (np.cos(phi)*np.cos(dec_r))
        cos_lha = np.clip(cos_lha, -1, 1)  # Numerical stability
        lha = np.arccos(cos_lha)
        return gha_r - lha  # Assuming East longitude
    
    def F(phi):
        """Function to find root of (azimuth equation residual)."""
        lam = longitude_from_lat(phi)
        lha = gha_r - lam
        residual = np.cos(phi)*np.tan(dec_r) - \
                   np.sin(phi)*np.cos(lha) - \
                   (1/np.tan(A_r))*np.sin(lha)
        return residual
    
    # Secant method
    phi0 = np.radians(dr_lat)
    phi1 = np.radians(dr_lat + 1)
    
    for _ in range(20):  # Max iterations
        f0, f1 = F(phi0), F(phi1)
        if abs(f1 - f0) < 1e-12:
            break
        phi_new = phi1 - f1 * (phi1 - phi0) / (f1 - f0)
        if abs(phi_new - phi1) < 1e-8:  # Convergence
            break
        phi0, phi1 = phi1, phi_new
    
    latitude = np.degrees(phi1)
    longitude = np.degrees(longitude_from_lat(phi1))
    
    return latitude, longitude
```

#### Worked Examples from Paper:

**Example 1: MIRFAK Star**
```
Date: 01 April 2019
Time: UTC 11:40:07
DR Position: 20°44.6'N, 107°06.4'E
True Position: 20°37.5'N, 107°16.3'E

Observations:
- Star: MIRFAK
- Altitude H = 34°31.9'
- Azimuth A = 317°18'

Celestial Body Data:
- GHA = 313°09.5'
- δ = 49°55.7'N

Secant Iteration Results:
| Iter | φ input    | F value      | φ solution  | Error      |
|------|------------|--------------|-------------|------------|
| 1    | 21°, 20°   | 0.00741...   | 20.6299°    | 0.6299°    |
| 2    | 20.6299°   | 9.39e-05     | 20.6252°    | 0.00466°   |
| 3    | 20.6252°   | 1.18e-06     | 20.6252°    | 5.91e-05°  |

Result: φ = 20°37.5'N, λ = 107°16.3'E
Error: 0.09 NM from true position
```

**Example 2: DUBHE Star**
```
Date: 30 April 2019
Time: UTC 10:51:22
DR Position: 20°40.9'N, 107°12.8'E
True Position: 20°45.2'N, 107°20.6'E

Observations:
- Star: DUBHE
- Altitude H = 41°27.1'
- Azimuth A = 22°93'

Celestial Body Data:
- GHA = 214°42.5'
- δ = 61°39.1'N

Result: φ = 20°45.2'N, λ = 107°20.4'E
Error: 0.17 NM from true position
Iterations: 4
```

#### Algorithm Flowchart:
```
1. Input: DR position P₁(φ₁, λ₁)
         Celestial body identified at time UTC
         
2. Determine celestial body coordinates:
   GHA, δ from nautical almanac at UTC
   
3. Measure from ship:
   H = true altitude (sextant + corrections)
   A = azimuth (compass + corrections)
   
4. Initialize Secant method:
   φ₀ = φ₁ (DR latitude)
   φ₁ = φ₁ ± 1°
   
5. Iterate using Eq. 11 until |φ_{n+1} - φ_n| < ε

6. Calculate longitude using Eq. 8

7. Output: Ship position (φ, λ)
```

#### Limitations Noted:
- Requires accurate azimuth measurement (compass error affects result)
- Large error in observed azimuth significantly impacts position
- Future work: investigate effect of azimuth error on position accuracy

#### Key Quotes:
- "this proposed technique used both parameters which are the altitude and the azimuth of a single celestial body"
- "it can be useful for some cases where there is only an identified single celestial body"
- "applying secant technique to solve the system is seen as one of new approaches of advanced mathematics in celestial navigation"

#### Relevance to Research:
- **Addresses SRQ1:** Novel algorithm for celestial position fix
- **Single-body capability:** Unique approach when only one body visible
- **Secant method:** Standard numerical technique (scipy.optimize available)
- **Accuracy:** 0.09-0.17 NM demonstrated in examples
- **Practical limitation identified:** Azimuth measurement error sensitivity
- **Complements other methods:** Use when conditions limit body visibility
- **Same author as SVD and azimuth methods:** Consistent notation/approach
- **Python implementation straightforward:** Direct translation from equations

---

## Article #18: Venkat (2022) - Autonomous Celestial Navigation Challenges

**Citation:** Venkat, R. (2022). Challenges to overcome limitations of human centric practices in celestial navigation at sea. *International Journal of Health Sciences*, 6(S6), 11080–11093. https://doi.org/10.53730/ijhs.v6nS6.13034

**Relevance:** Review paper providing contextual framework for why autonomous celestial navigation systems are needed. Useful for Introduction and Literature Review sections.

**Primary Value:** Background/motivation rather than new algorithms.

### Content for Introduction Section:

#### GPS Vulnerabilities and Need for Backup:
- "GPS is operated by the US Navy and allows merchant vessels of all countries to use it. US agencies can tweak GPS accuracy or introduce positional errors during conflicts/wars, international sanctions, or political imbroglio."
- "Satellite systems are prone to security threats like computer malware, electromagnetic pulse attacks and jamming"
- "satellites can be damaged or disabled during Sunspot activity"
- "This vulnerability of GPS necessitates the backup system for marine navigation, which is standalone, fool-proof and highly accurate"

#### INS Limitations:
- "Inertial Navigation System (INS) is an alternative to GPS, capable of accurately maintaining the ship's dead reckoning. Still, the position's accuracy is not guaranteed for extended periods."
- "INS needs regular alignment with an external reference system"
- "Celestial navigation may not be feasible during a cloudy day; however, INS is regarded as an 'excellent bad-weather flywheel'"

### Content for Literature Review Section:

#### Accuracy Comparison Table (Table 1):
| Factor | Celestial Navigation | Satellite Navigation |
|--------|---------------------|---------------------|
| Average Fix Accuracy | 1 to 1.5 NM (1850-2750 m) | 100-300 m |
| Reliability | 60% | 95% |
| Time Lag | Yes | Instantaneous |
| External Dependence | Autonomous | Depends on satellite signals |
| Weather Dependence | Restricted window | No |

#### Error Sources in Celestial Navigation:
1. **Instrument error** (sextant)
2. **Manual observation error** (improper observation)
3. **Unclear horizon/limb** (affects altitude reading)
4. **Refraction error** (atmospheric)
5. **Height of eye / dip error**
6. **Parallax error**
7. **Manual sight reduction of azimuth**

### Content for Methodology Section:

#### Azimuth Formulas (from Nguyen 2017):

**Equation 1 - When altitude known:**
```
Z = arccos[(sin(δ) - sin(φ) × sin(Hc)) / (cos(φ) × cos(Hc))]
```

**Equation 2 - When Hour Angle and Declination known:**
```
Z = arctan[sin(LHA) / (cos(φ)·tan(δ) - sin(φ)·cos(LHA))]

where:
LHA = GHA - λ (Local Hour Angle)
φ = observer latitude
δ = declination of celestial body
```

**Equation 3 - With assumed latitude:**
```
Z = arccos[(sin(φ₂) - sin(φ₁)·cos(δ)) / (cos(φ₁)·sin(δ))]
```

### Technology Discussion (for context):

#### Automated Observation Systems:
- **Star identification apps:** Google Sky, Sky Safari, Star Tracker, Skyview
- **Camera sensors:** Active Pixel Sensor (APS) vs Charge-Coupled Device (CCD)
- **Wide Field of View (WFOV)** cameras for star pattern recognition
- **Digital Zenith Camera System (DZCS)** for plumb line direction measurement
- **AURIGA software** for geodetic astronomy image processing

#### Projective Camera Geometry:
- Pinhole camera model for star field imaging
- Invariant asterism descriptors for pattern matching
- Star catalogue matching algorithms (e.g., Jiangbo 2019)

### Key Quotes:

- "The main ingredient of accurate celestial position fixing is the identification of celestial bodies to be observed"
- "Prevention of human observational error refines the positional accuracy of celestial navigation"
- "if a celestial body's azimuth is observed with reliable accuracy using a technologically advanced optical device instead of the current practice of human observation, the position fixing by celestial navigation at sea becomes highly reliable and accurate"
- "Many advanced computer-based algorithms/software are currently available to compute a ship's position coordinates based on celestial bodies' altitude and azimuth"

### References to Follow Up:
- Li (2022) - Adaptive robust filtering algorithm *(Already accepted as Article #10)*
- Nguyen (2017) - Azimuth method *(Already accepted as Article #12)*
- Parish et al. (2010) - Stellar Positioning System (SPS)
- Benjamin B. et al. (2009) - Star identification algorithms survey
- Jiangbo (2019) - Star identification based on oriented singular value

### Relevance to Research:
- **Introduction context:** GPS vulnerabilities, need for autonomous backup systems
- **Literature review:** Error source taxonomy, accuracy benchmarks
- **Motivation:** Why algorithmic automation matters
- **NOT for Methodology:** Formulas are from other sources (Nguyen 2017)
- **Supports hypothesis:** Computational methods can match/exceed traditional accuracy if observation errors eliminated

### Limitations of This Paper:
- Published in unusual venue (health sciences journal)
- Review/conceptual paper, not empirical research
- No new algorithms - relies on citations to other work
- No validation or simulation results
- Proposes future research directions rather than presenting solutions

---

## Rejected Articles

1. **PySilSub** (Martin et al., 2023) - Vision science/photoreceptor research, not navigation
2. **Asteroid (2059) Baboquivari** - Asteroid characterization, not navigation
3. **Space Debris GPU Detection** - Space surveillance, not position fixing
4. **Autonomous satellite navigation using starlight refraction** (Ning et al., 2013) - Satellite orbital navigation, not traditional celestial navigation
5. **The Mariner's Sextant and the Royal Society** (Cotter, 1978) - History of physical instrument, no sight reduction algorithms
6. **A New Navigation Computer** (Richey, 1978) - Product review without implementable algorithms
7. **Parallelization of Bresenham's Line and Circle Algorithms** (Wright, 1990) - Computer graphics, not navigation
8. **Intersection algorithms for lines and circles** (Middleditch et al., 1989) - CAD/computational geometry, not navigation
9. **The high precision positioning algorithm of circular landmark center** (Cui et al., 2014) - Computer vision/camera calibration, not navigation
10. **On a Circle Placement Problem** (Chazelle & Lee, 1986) - Operations research/facility siting, not navigation
11. **Searching for the Centre of a Circle** (Biedl et al.) - Robotics/avalanche rescue algorithms, not navigation

---

## Article #19: Kaplan (1995) - Moving Observer Celestial Fix

**Citation:** Kaplan, G. H. (1995). Determining the Position of a Vessel from Celestial Observations and Motion. *NAVIGATION: Journal of The Institute of Navigation*, 42(4), 631-648. https://doi.org/10.1002/j.2161-4296.1995.tb01911.x

**Relevance:** Foundational paper from U.S. Naval Observatory presenting complete algorithm for celestial fix with moving observer. Treats fix as "orbit correction problem" solving for 4 parameters simultaneously: latitude, longitude, course, speed.

**Primary Value:** Complete mathematical framework with all partial derivatives for Python implementation.

### Core Concept:

Celestial navigation as **orbit correction problem:**
- Moving object = observer's ship
- Observations = one-dimensional sextant altitudes
- "Orbital elements" = latitude, longitude, course, speed (4 parameters)
- Laws of motion = rhumb-line sailing formulas
- Solution method = differential correction with least squares

### Content for Methodology Section:

#### Altitude Equation:
```
Hc = arcsin[sin(δ)·sin(φ) + cos(δ)·cos(φ)·cos(GHA + λ)]     [Eq. A1]

where:
  Hc = computed altitude
  δ = declination of celestial body
  φ = observer latitude
  GHA = Greenwich Hour Angle
  λ = observer longitude
```

#### Partial Derivatives of Altitude:
```
∂Hc/∂φ = sec(Hc)·[sin(δ)·cos(φ) - cos(δ)·sin(φ)·cos(GHA + λ)]     [Eq. A2a]

∂Hc/∂λ = -sec(Hc)·[cos(δ)·cos(φ)·sin(GHA + λ)]                    [Eq. A2b]
```

#### Conditional Equation (Core Algorithm):
```
a = (∂Hc/∂φ · ∂f/∂φ₀ + ∂Hc/∂λ · ∂g/∂φ₀)·Δφ₀
  + (∂Hc/∂φ · ∂f/∂λ₀ + ∂Hc/∂λ · ∂g/∂λ₀)·Δλ₀
  + (∂Hc/∂φ · ∂f/∂C  + ∂Hc/∂λ · ∂g/∂C)·ΔC
  + (∂Hc/∂φ · ∂f/∂S  + ∂Hc/∂λ · ∂g/∂S)·ΔS                        [Eq. 4]

where:
  a = Ho - Hc (altitude intercept)
  f, g = sailing formulas for latitude and longitude
  Δφ₀, Δλ₀, ΔC, ΔS = corrections to initial estimates
```

#### Earth Curvature (WGS-84):
```
M = a(1 - e²) / (1 - e²·sin²φ)^(3/2)    (meridian radius of curvature)
N = a / (1 - e²·sin²φ)^(1/2)             (prime vertical radius)

where:
  a = 6378.137 km (equatorial radius)
  f = 1/298.257223563 (flattening)
  e² = 2f - f² (eccentricity squared)
```

#### Simplified Sailing Formulas (short tracks):
```
φ' = φ₀ + S(t - t₀)·cos(C) / M₀
λ' = λ₀ + S(t - t₀)·sin(C) / (N₀·cos(φ_avg))                      [Eq. A5]
```

#### Full Sailing Formulas (rhumb line, oblate earth):
```
φ = φ₀ + M₀/(a(1-e²))·[(1 - 3e²/4)(φ' - φ₀) + (3e²/8)(sin(2φ') - sin(2φ₀))]

λ = λ₀ + tan(C)·{ln[tan(π/4 + φ/2)] + (e/2)·ln[(1 - e·sin(φ))/(1 + e·sin(φ))]
       - ln[tan(π/4 + φ₀/2)] - (e/2)·ln[(1 - e·sin(φ₀))/(1 + e·sin(φ₀))]}    [Eq. A7]
```

#### Partial Derivatives for Sailing Formulas:
```
∂f/∂φ₀, ∂f/∂λ₀, ∂f/∂C, ∂f/∂S  (4 derivatives for latitude function)
∂g/∂φ₀, ∂g/∂λ₀, ∂g/∂C, ∂g/∂S  (4 derivatives for longitude function)

Total: 10 partial derivatives needed per observation
(Full expressions given in Appendix A6-A9 of paper)
```

### Algorithm Procedure:

```
1. Determine reference time t₀ for fix
   Set initial estimates: φ₀, λ₀, C, S

2. Decide whether to solve for course/speed
   (Requires ≥8 observations, good azimuth distribution, several hours span)

3. Correct each observation: hs → Ho
   (index error, dip, refraction)

4. For each observation at time t:
   a) Compute estimated position (φ, λ) using sailing formulas
   b) Compute Hc from ephemeris data (GHA, Dec)
   c) Compute altitude intercept: a = Ho - Hc
   d) Evaluate all 10 partial derivatives
   e) Form conditional equation (Eq. 4)
   f) Add row to design matrix for least squares

5. Solve least squares for Δφ₀, Δλ₀, ΔC, ΔS
   Compute variance-covariance matrix

6. Update estimates:
   φ₀ ← φ₀ + Δφ₀
   λ₀ ← λ₀ + Δλ₀
   C ← C + ΔC
   S ← S + ΔS

7. Iterate if needed (usually 2 iterations sufficient)

8. Output: Fix position (φ₀, λ₀), track (C, S), error ellipse
```

### Sample Calculation Results:

```
Test Case:
- 8 observations over 8 hours (Sun, Moon, stars)
- Ship: 45°N, 50°W, Course 330°, Speed 20 kn
- Observation noise: σ = 0.7 arcmin

Results after 2 iterations:
- Position error: 0.41 NM
- CEP: 0.63 NM
- Course recovery: within 0.1°
- Speed recovery: within 0.03 kn

Comparison (without solving for C, S):
- Standard fix (twilight only): 0.78 NM error
- Standard fix (all 8 obs): 3.67 NM error (much worse!)
```

#### Algorithm Accuracy vs Observation Error:
| Obs Error (arcmin) | New Alg CEP | New Alg Actual | Std Fix CEP | Std Fix Actual |
|-------------------|-------------|----------------|-------------|----------------|
| 0.7               | 0.63 NM     | 0.41 NM        | 0.80 NM     | 0.78 NM        |
| 0.3               | 0.27 NM     | 0.17 NM        | 0.51 NM     | 0.53 NM        |
| 0.1               | 0.090 NM    | 0.058 NM       | 0.36 NM     | 0.45 NM        |
| 0.02              | 0.018 NM    | 0.013 NM       | 0.31 NM     | 0.42 NM        |
| 0                 | 0.001 NM    | 0.001 NM       | 0.29 NM     | 0.42 NM        |

### Python Implementation Skeleton:

```python
import numpy as np
from scipy.linalg import lstsq

# WGS-84 constants
A_EARTH = 6378.137  # km
F_EARTH = 1/298.257223563
E2 = 2*F_EARTH - F_EARTH**2

def curvature_radii(phi):
    """Compute M (meridian) and N (prime vertical) radii."""
    sin_phi = np.sin(phi)
    denom = np.sqrt(1 - E2 * sin_phi**2)
    M = A_EARTH * (1 - E2) / denom**3
    N = A_EARTH / denom
    return M, N

def computed_altitude(phi, lam, gha, dec):
    """Compute Hc from observer position and body coordinates."""
    lha = gha + lam  # Local Hour Angle
    sin_Hc = (np.sin(dec) * np.sin(phi) + 
              np.cos(dec) * np.cos(phi) * np.cos(lha))
    return np.arcsin(sin_Hc)

def altitude_partials(phi, lam, gha, dec, Hc):
    """Partial derivatives of Hc with respect to phi and lambda."""
    lha = gha + lam
    sec_Hc = 1 / np.cos(Hc)
    
    dHc_dphi = sec_Hc * (np.sin(dec) * np.cos(phi) - 
                         np.cos(dec) * np.sin(phi) * np.cos(lha))
    dHc_dlam = -sec_Hc * np.cos(dec) * np.cos(phi) * np.sin(lha)
    
    return dHc_dphi, dHc_dlam

def sailing_position(phi0, lam0, C, S, dt, M0, N0):
    """Compute position at time t = t0 + dt using sailing formulas."""
    # Simplified formulas for short tracks
    phi_avg = phi0 + 0.5 * S * dt * np.cos(C) / M0
    phi = phi0 + S * dt * np.cos(C) / M0
    lam = lam0 + S * dt * np.sin(C) / (N0 * np.cos(phi_avg))
    return phi, lam

def sailing_partials(phi0, lam0, C, S, dt, phi, M0, N0):
    """
    Compute 8 partial derivatives of sailing formulas.
    Returns: df_dphi0, df_dlam0, df_dC, df_dS,
             dg_dphi0, dg_dlam0, dg_dC, dg_dS
    """
    # Simplified version - see paper Appendix for full oblate earth formulas
    cos_C = np.cos(C)
    sin_C = np.sin(C)
    cos_phi = np.cos(phi)
    
    # Latitude partials
    df_dphi0 = 1.0  # + small correction terms
    df_dlam0 = 0.0
    df_dC = -S * dt * sin_C / M0
    df_dS = dt * cos_C / M0
    
    # Longitude partials
    dg_dphi0 = S * dt * sin_C * np.tan(phi) / (N0 * cos_phi)
    dg_dlam0 = 1.0
    dg_dC = S * dt * cos_C / (N0 * cos_phi)
    dg_dS = dt * sin_C / (N0 * cos_phi)
    
    return (df_dphi0, df_dlam0, df_dC, df_dS,
            dg_dphi0, dg_dlam0, dg_dC, dg_dS)

def build_design_matrix(observations, phi0, lam0, C, S, t0):
    """
    Build design matrix for least squares solution.
    
    observations: list of (t, Ho, gha, dec) tuples
    Returns: design matrix A, observation vector b
    """
    M0, N0 = curvature_radii(phi0)
    
    A = []
    b = []
    
    for t, Ho, gha, dec in observations:
        dt = t - t0  # hours
        
        # Estimated position at observation time
        phi, lam = sailing_position(phi0, lam0, C, S, dt, M0, N0)
        
        # Computed altitude and partials
        Hc = computed_altitude(phi, lam, gha, dec)
        dHc_dphi, dHc_dlam = altitude_partials(phi, lam, gha, dec, Hc)
        
        # Sailing partials
        (df_dphi0, df_dlam0, df_dC, df_dS,
         dg_dphi0, dg_dlam0, dg_dC, dg_dS) = sailing_partials(
             phi0, lam0, C, S, dt, phi, M0, N0)
        
        # Coefficients for conditional equation (Eq. 4)
        coef_phi0 = dHc_dphi * df_dphi0 + dHc_dlam * dg_dphi0
        coef_lam0 = dHc_dphi * df_dlam0 + dHc_dlam * dg_dlam0
        coef_C = dHc_dphi * df_dC + dHc_dlam * dg_dC
        coef_S = dHc_dphi * df_dS + dHc_dlam * dg_dS
        
        A.append([coef_phi0, coef_lam0, coef_C, coef_S])
        b.append(Ho - Hc)  # altitude intercept
    
    return np.array(A), np.array(b)

def solve_fix(observations, phi0_init, lam0_init, C_init, S_init, t0,
              solve_motion=True, max_iter=3):
    """
    Iterative least squares solution for celestial fix.
    
    Parameters:
    - observations: list of (t, Ho, gha, dec) tuples
    - phi0_init, lam0_init: initial position estimate (radians)
    - C_init, S_init: initial course (radians) and speed (km/h)
    - t0: reference time for fix (hours)
    - solve_motion: if True, solve for C and S as well
    - max_iter: maximum iterations
    
    Returns:
    - phi0, lam0, C, S: solution
    - cov: covariance matrix
    """
    phi0, lam0, C, S = phi0_init, lam0_init, C_init, S_init
    
    for iteration in range(max_iter):
        A, b = build_design_matrix(observations, phi0, lam0, C, S, t0)
        
        if not solve_motion:
            A = A[:, :2]  # Only position columns
        
        # Least squares solution
        result, residuals, rank, s = lstsq(A, b)
        
        # Update parameters
        phi0 += result[0]
        lam0 += result[1]
        if solve_motion and len(result) >= 4:
            C += result[2]
            S += result[3]
        
        # Check convergence
        if np.all(np.abs(result) < 1e-8):
            break
    
    # Compute covariance matrix
    residuals = b - A @ result
    sigma2 = np.sum(residuals**2) / (len(observations) - len(result))
    cov = sigma2 * np.linalg.inv(A.T @ A)
    
    return phi0, lam0, C, S, cov
```

### Key Quotes:

- "The object of celestial navigation, as traditionally practiced, is the determination of the latitude and longitude of a vessel at a specific time"
- "Our problem is quite similar to 'orbit correction' problems faced by astronomers who deal with the dynamics of solar system bodies"
- "Given enough observations, suitably distributed in time and azimuth, we should be able to obtain an estimate of the average over-bottom track of the vessel as part of the solution for the fix"
- "The most important weakness of such procedures is well known: they require data on the motion of the observer's ship over bottom, and the course and speed values used may not be accurate"
- "For hand-held sextant observations, difficulties may arise for observations spread over more than about 0.5 h"

### Important References Cited:
- De Wit (1974) - Optimal estimation of multi-star fix
- Metcalf & Metcalf (1991) - Overdetermined celestial fix *(Already accepted)*
- Severance (1989) - Overdetermined fix by iteration
- A'Hearn & Rossano (1977) - Two body fix by calculator
- Van Allen (1981) - Two star sight analytical solution

### Relevance to Research:
- **Addresses SRQ1:** Complete algorithm for position fix
- **Addresses SRQ2:** Handles multiple observations over time
- **Moving observer:** Solves running fix problem directly
- **Least squares:** Standard scipy.linalg.lstsq applicable
- **Error ellipse:** Covariance matrix output for confidence region
- **Oblate earth:** WGS-84 compatible (Skyfield default)
- **Rhumb line:** Standard vessel track model
- **4-parameter solution:** Position + motion in one solve
- **Iterative convergence:** 2 iterations typically sufficient
- **U.S. Naval Observatory:** Authoritative source, Navy-tested algorithm

### Algorithm Advantages:
1. Self-correcting for course/speed errors
2. Exploits information content in extended observation sets
3. Properly propagates positional error along track
4. Works with automated star trackers (arcsecond accuracy)
5. Handles day+twilight combined observations

### Limitations:
- Requires good azimuth distribution
- Needs ≥8 observations for motion solution
- Assumes constant course/speed (rhumb line)
- East-west courses need special handling (tan C singularity)

---

## Article #20: Mederos (2024) - Computational Celestial Navigation

**Citation:** Mederos, L. (2024). Computational celestial navigation at sea. *JCAAC*, 1, 5-14.

**Relevance:** Modern computational approach to celestial navigation with clean vector-matrix formulation for circle of position intersection. Provides explicit closed-form solution for two-body fix.

**Primary Value:** Clean mathematical formulation directly implementable in Python.

### Core Concept:

Replace traditional graphical methods with computational solution:
- Two circles of equal altitude → find intersection points
- Vector-matrix approach using Earth-centered reference frame
- Closed-form solution (no iteration needed for basic case)
- Handle non-simultaneous observations with radius correction

### Content for Methodology Section:

#### Reference Frame:
```
Earth-centered reference system:
- X axis: pointing towards Greenwich meridian
- Z axis: pointing to North Celestial Pole
- Y axis: completing right-handed system
```

#### Unit Position Vectors for Stars:
```
r̂₁ = (x₁, y₁, z₁) = (cos(δ₁)·cos(α₁ - θ₁), cos(δ₁)·sin(α₁ - θ₁), sin(δ₁))
r̂₂ = (x₂, y₂, z₂) = (cos(δ₂)·cos(α₂ - θ₂), cos(δ₂)·sin(α₂ - θ₂), sin(δ₂))    [Eq. 3.2]

where:
  δᵢ = declination of star i
  αᵢ = right ascension of star i
  θᵢ = sidereal time at observation i
```

#### Circle Intersection Equations:
```
r̂ · r̂₁ = sin(a₁)     (observer on circle 1)
r̂ · r̂₂ = sin(a₂)     (observer on circle 2)
|r̂|² = 1             (on unit sphere)                                         [Eq. 3.3]

where:
  aᵢ = true altitude of star i
  r̂ = unit position vector of observer
```

#### Parametric Solution (using z as parameter):
```
Matrix form:
┌       ┐   ┌   ┐     ┌              ┐
│ x₁ y₁ │   │ x │     │ sin(a₁) - z·z₁│
│ x₂ y₂ │ × │ y │  =  │ sin(a₂) - z·z₂│                                        [Eq. 3.4]
└       ┘   └   ┘     └              ┘

Solution:
x = (bₓ + cₓ·z) / cᵤ
y = (bᵧ + cᵧ·z) / cᵤ                                                           [Eq. 3.5]

where:
  c = r̂₁ × r̂₂                    (cross product)
  b = r₂⊥·sin(a₁) - r₁⊥·sin(a₂)
  rᵢ⊥ = (yᵢ, -xᵢ, 0)
```

#### Quadratic Equation for z:
```
z² + 2Bz + C = 0                                                               [Eq. 3.7]

where:
  B = (b · c) / c² = (b/c)·cos(γ)
  C = (b/c)² - (cᵤ/c)²                                                         [Eq. 3.8]
  γ = angle between vectors b and c
```

#### Final Solution (two intersection points):
```
z± = (b/c)·(-cos(γ) ± √[(cᵤ/b)² - sin²(γ)])                                    [Eq. 3.9]

x± = (bₓ + cₓ·z±) / cᵤ
y± = (bᵧ + cᵧ·z±) / cᵤ

Convert to spherical coordinates:
φ± = arcsin(z±)                    (latitude)
λ± = arctan2(y±, x±)               (longitude)
```

### Non-Simultaneous Observations:

#### Running Fix Correction:
```
If boat sailed distance D at course Rᵥ between observations:

Modified radius for first circle:
  90° - a₁ ± x

where:
  x = D · cos(α)                                                               [Eq. 3.10]
  α = angle between course Rᵥ and star azimuth Z
  
+ if sailing away from circle center
- if sailing toward circle center
```

#### Azimuth Calculation:
When azimuth measurement not available:
1. Use estimated position to calculate azimuth (error negligible due to large circle radius)
2. Or: First solve assuming stationary observer, then iterate with corrected radius

### Python Implementation:

```python
import numpy as np

def star_position_vector(dec, ra, sidereal_time):
    """
    Compute unit position vector for star in Earth-centered frame.
    All angles in radians.
    """
    hour_angle = ra - sidereal_time
    x = np.cos(dec) * np.cos(hour_angle)
    y = np.cos(dec) * np.sin(hour_angle)
    z = np.sin(dec)
    return np.array([x, y, z])

def two_body_fix(dec1, ra1, theta1, alt1, dec2, ra2, theta2, alt2):
    """
    Calculate position from two star observations.
    
    Parameters:
    - dec1, ra1: declination and right ascension of star 1 (radians)
    - theta1: sidereal time at observation 1 (radians)
    - alt1: true altitude of star 1 (radians)
    - dec2, ra2, theta2, alt2: same for star 2
    
    Returns:
    - (lat1, lon1), (lat2, lon2): two possible positions (radians)
    """
    # Unit position vectors for stars (circle centers)
    r1 = star_position_vector(dec1, ra1, theta1)
    r2 = star_position_vector(dec2, ra2, theta2)
    
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    
    # Cross product c = r1 × r2
    c = np.cross(r1, r2)
    cx, cy, cz = c
    c_mag = np.linalg.norm(c)
    
    # Perpendicular vectors
    r1_perp = np.array([y1, -x1, 0])
    r2_perp = np.array([y2, -x2, 0])
    
    # Vector b
    sin_a1 = np.sin(alt1)
    sin_a2 = np.sin(alt2)
    b = r2_perp * sin_a1 - r1_perp * sin_a2
    bx, by, bz = b
    b_mag = np.linalg.norm(b)
    
    # Coefficients for quadratic equation
    cos_gamma = np.dot(b, c) / (b_mag * c_mag)
    sin_gamma_sq = 1 - cos_gamma**2
    
    B = (b_mag / c_mag) * cos_gamma
    C = (b_mag / c_mag)**2 - (cz / c_mag)**2
    
    # Solve quadratic: z² + 2Bz + C = 0
    discriminant = B**2 - C
    if discriminant < 0:
        raise ValueError("No intersection - circles don't intersect")
    
    z_plus = -B + np.sqrt(discriminant)
    z_minus = -B - np.sqrt(discriminant)
    
    # Compute x, y for each solution
    def compute_xy(z):
        x = (bx + cx * z) / cz
        y = (by + cy * z) / cz
        return x, y
    
    x_plus, y_plus = compute_xy(z_plus)
    x_minus, y_minus = compute_xy(z_minus)
    
    # Convert to spherical coordinates (lat, lon)
    def to_spherical(x, y, z):
        lat = np.arcsin(z)
        lon = np.arctan2(y, x)
        return lat, lon
    
    pos1 = to_spherical(x_plus, y_plus, z_plus)
    pos2 = to_spherical(x_minus, y_minus, z_minus)
    
    return pos1, pos2

def apply_running_fix_correction(alt, azimuth, course, distance):
    """
    Modify altitude for non-simultaneous observation.
    
    Parameters:
    - alt: true altitude (radians)
    - azimuth: star azimuth at observation (radians)
    - course: vessel course (radians)
    - distance: distance sailed since observation (same units as circle radius)
    
    Returns:
    - corrected altitude (radians)
    """
    alpha = azimuth - course
    x = distance * np.cos(alpha)
    
    # x positive = sailing away from GP, radius increases, alt decreases
    # Convert distance to angular measure (assuming nautical miles and degrees)
    # 1 NM ≈ 1 arcmin = 1/60 degree
    x_rad = np.radians(x / 60)  # if distance in NM
    
    return alt - x_rad  # Modified altitude

def select_correct_position(pos1, pos2, azimuth1, azimuth2):
    """
    Select correct position based on observed azimuths.
    The correct position should have similar azimuths to observed values.
    """
    # Calculate expected azimuths from each position
    # ... (implementation depends on star GP coordinates)
    # Choose position where calculated azimuths match observed
    pass
```

### Algorithm Summary:

```
Input:
  Star 1: (UT₁, altitude a₁, azimuth Z₁, RA α₁, Dec δ₁)
  Star 2: (UT₂, altitude a₂, azimuth Z₂, RA α₂, Dec δ₂)

Step 1: Convert UT to sidereal time θ

Step 2: Compute unit vectors r̂₁, r̂₂ for star GPs

Step 3: If non-simultaneous:
        - Apply running fix correction to a₁

Step 4: Compute vectors c = r̂₁ × r̂₂ and b

Step 5: Solve quadratic for z±

Step 6: Compute (x±, y±) from z±

Step 7: Convert to (φ±, λ±) spherical coordinates

Step 8: Select correct position using azimuths

Output: Position (φ, λ)
```

### Key Quotes:

- "Today we can use a computer to avoid these two requirements" [estimated position and graphical work]
- "The intersection points of these circles can be easily obtained by calculating the intersection (a straight line) of the two planes perpendicular to the directions r̂₁ and r̂₂"
- "the radius of the circle is enormous compared to the distance from our estimated position to the real position"
- "Using the estimated position to calculate the star azimuth will produce a virtually exact result, because the centre of the circle of altitude is so far away"

### Comparison with Other Methods:

| Method | Observations | Iteration | Output |
|--------|--------------|-----------|--------|
| Mederos (this) | 2 | No | 2 positions, select by azimuth |
| Kaplan (1995) | ≥4 | Yes | Position + motion + covariance |
| Chiesa (1990) | ≥2 | No | Direct COP intersection |
| Intercept method | ≥2 | No | LOP intersection (graphical) |

### Relevance to Research:
- **Addresses SRQ1:** Clean algorithm for two-body fix
- **Direct solution:** No iteration needed (closed-form)
- **Vector formulation:** NumPy-native implementation
- **Running fix:** Explicit radius correction formula
- **Modern source:** 2024 publication, current best practices
- **Educational value:** Clear derivation from first principles
- **Historical context:** Good introduction section on longitude problem

### Limitations:
- Only handles 2-body case (not overdetermined)
- Assumes nearly simultaneous observations (or requires correction)
- Does not provide error estimates
- Spherical earth assumption (no WGS-84 corrections)

---

## Article #21: Metcalf & Metcalf (1994) - Piloting with Celestial Algorithms

**Citation:** Metcalf, T. R. & Metcalf, F. T. (1994). Piloting with Celestial Algorithms. *NAVIGATION: Journal of The Institute of Navigation*, 41(2), 207-214. https://doi.org/10.1002/j.2161-4296.1994.tb02572.x

**Relevance:** Shows how coastal piloting observations can be converted to equivalent celestial observations, enabling a unified algorithm for both celestial and terrestrial navigation.

**Primary Value:** Extends celestial fix algorithms to handle combined celestial/piloting fixes.

### Core Concept:

Piloting observations yield circles of position, just like celestial observations:
- **Range to object** → circle centered on object
- **Bearing to object** → great circle (altitude = 0°)
- **Horizontal angle** → circle of equal angle

Convert to equivalent (GHA, Dec, H₀) → use same celestial algorithm.

### Content for Methodology Section:

#### 1. Range Observation:
```
If object at position (λ, l) and range = d nautical miles:

GHA = λ                     (add 360° if λ < 0)
Declination = l
H₀ = 90° - d/60             (d in NM, 1 NM = 1 arcmin)                      [Eq. 1]

Note: No corrections for dip/refraction needed - these are equivalent observations.
```

#### 2. Bearing Observation:
```
Bearing to known object yields great circle of position.
Great circle = celestial observation with H₀ = 0°

Objects: P_obs = (λ_obs, l_obs), bearing C_obs
DR position: P_dr = (λ_dr, l_dr)

Initial great circle course from P₁ to P₂:
C_i(P₁, P₂) = arctan[sin(λ₂ - λ₁) / (cos(l₁)·tan(l₂) - sin(l₁)·cos(λ₂ - λ₁))]   [Eq. 2]

Corrected bearing (for radio bearings over long distances):
C = C_obs + C_i(P_obs, P_dr) - C_i(P_dr, P_obs)                              [Eq. 3]

Equivalent celestial observation:
GHA = λ_obs + arcsin[tan(l_obs) / tan(l_v)]        (± based on C > 180°)
Declination = l_v + 90°                             (± based on hemisphere)
H₀ = 0°                                                                      [Eq. 4]

where vertex latitude:
l_v = ±arccos(|cos(l_obs)·sin(C)|)                                           [Eq. 5]
```

#### 3. Horizontal Angle Observation:
```
Horizontal angle A between two known objects yields circle of equal angle.

Objects at (λ₁, l₁) and (λ₂, l₂), observed angle A:

Bearing from Object 1 to Object 2:
C₁₂ = arctan[(λ₁ - λ₂) / ln(tan(π/4 + l₂/2) / tan(π/4 + l₁/2))]           [Eq. 7]

Half-distance between objects:
a = |l₂ - l₁| × 60 / |cos(C₁₂)|    (if |cos(C₁₂)| > 0.0001)
a = (λ₂ - λ₁) × cos((l₁ + l₂)/2) × 60 / sin(C₁₂)    (otherwise)          [Eq. 8]

Radius of circle: r = a / sin(A)
Distance to center: a × cot(A)

Equivalent celestial observation:
GHA = λ₁ - r × sin(C₁₂ ± (90° - A)) / cos(l₁)
Declination = l₁ + r × cos(C₁₂ ± (90° - A)) / 60
H₀ = 90° - r/60                                                             [Eq. 6]

Select + or - based on:
- Acute angle (A < 90°): use center closest to DR
- Obtuse angle (A > 90°): use center farthest from DR
```

### Circle of Equal Angle Derivation:

```
Law of cosines in plane geometry:
(2a)² = r₁² + r₂² - 2r₁r₂·cos(A)

where r₁, r₂ = distances to objects 1 and 2

Simplifies to:
x² + (y ± a·cot(A))² = (a/sin(A))²                                          [Eq. A-2]

This is a circle with:
- Radius = a/sin(A)
- Center at (0, ±a·cot(A))

Both circles pass through the two observed objects.
```

### Python Implementation:

```python
import numpy as np

def range_to_celestial(obj_lon, obj_lat, range_nm):
    """
    Convert range observation to equivalent celestial observation.
    
    Parameters:
    - obj_lon, obj_lat: object position (degrees, W positive)
    - range_nm: range to object in nautical miles
    
    Returns:
    - (GHA, Dec, H0) equivalent celestial observation
    """
    gha = obj_lon if obj_lon >= 0 else obj_lon + 360
    dec = obj_lat
    h0 = 90 - range_nm / 60  # Convert NM to degrees
    
    return gha, dec, h0

def bearing_to_celestial(obj_lon, obj_lat, bearing, dr_lon=None, dr_lat=None):
    """
    Convert bearing observation to equivalent celestial observation.
    
    Parameters:
    - obj_lon, obj_lat: object position (degrees)
    - bearing: true bearing from observer to object (degrees)
    - dr_lon, dr_lat: dead reckoning position (for correction)
    
    Returns:
    - (GHA, Dec, H0) equivalent celestial observation
    """
    # Apply bearing correction for long-distance radio bearings
    if dr_lon is not None and dr_lat is not None:
        c_obs_to_dr = great_circle_course(obj_lon, obj_lat, dr_lon, dr_lat)
        c_dr_to_obs = great_circle_course(dr_lon, dr_lat, obj_lon, obj_lat)
        bearing = bearing + c_obs_to_dr - c_dr_to_obs
        bearing = bearing % 360
    
    # Calculate vertex of great circle
    obj_lat_rad = np.radians(obj_lat)
    bearing_rad = np.radians(bearing)
    
    l_v = np.degrees(np.arccos(np.abs(np.cos(obj_lat_rad) * np.sin(bearing_rad))))
    if obj_lat < 0:
        l_v = -l_v
    
    # Calculate GHA from vertex
    # (simplified - full formula in paper Eq. 5)
    
    dec = l_v + 90 if l_v <= 0 else l_v - 90
    h0 = 0.0  # Great circle
    
    return gha, dec, h0

def great_circle_course(lon1, lat1, lon2, lat2):
    """
    Calculate initial great circle course from point 1 to point 2.
    """
    lat1_r = np.radians(lat1)
    lat2_r = np.radians(lat2)
    dlon = np.radians(lon2 - lon1)
    
    num = np.sin(dlon)
    den = np.cos(lat1_r) * np.tan(lat2_r) - np.sin(lat1_r) * np.cos(dlon)
    
    course = np.degrees(np.arctan2(num, den))
    return course % 360

def horizontal_angle_to_celestial(obj1_lon, obj1_lat, obj2_lon, obj2_lat, 
                                   angle, dr_lon, dr_lat):
    """
    Convert horizontal angle observation to equivalent celestial observation.
    
    Parameters:
    - obj1_lon, obj1_lat: first object position
    - obj2_lon, obj2_lat: second object position
    - angle: horizontal sextant angle between objects (degrees)
    - dr_lon, dr_lat: dead reckoning position
    
    Returns:
    - (GHA, Dec, H0) equivalent celestial observation
    """
    # Calculate bearing and distance between objects
    c12 = rhumb_line_bearing(obj1_lon, obj1_lat, obj2_lon, obj2_lat)
    a = half_distance(obj1_lon, obj1_lat, obj2_lon, obj2_lat, c12)
    
    angle_rad = np.radians(angle)
    
    # Radius of circle (in NM)
    radius = a / np.sin(angle_rad)
    
    # Two possible centers
    alpha_plus = c12 + (90 - angle)
    alpha_minus = c12 - (90 - angle)
    
    # Calculate center positions
    r_deg = radius / 60  # Convert NM to degrees
    
    center1_lon = obj1_lon - r_deg * np.sin(np.radians(alpha_plus)) / np.cos(np.radians(obj1_lat))
    center1_lat = obj1_lat + r_deg * np.cos(np.radians(alpha_plus))
    
    center2_lon = obj1_lon - r_deg * np.sin(np.radians(alpha_minus)) / np.cos(np.radians(obj1_lat))
    center2_lat = obj1_lat + r_deg * np.cos(np.radians(alpha_minus))
    
    # Select correct center based on angle type and DR position
    dist1 = np.sqrt((center1_lon - dr_lon)**2 + (center1_lat - dr_lat)**2)
    dist2 = np.sqrt((center2_lon - dr_lon)**2 + (center2_lat - dr_lat)**2)
    
    if angle < 90:  # Acute - use closer center
        if dist1 < dist2:
            gha, dec = center1_lon, center1_lat
        else:
            gha, dec = center2_lon, center2_lat
    else:  # Obtuse - use farther center
        if dist1 > dist2:
            gha, dec = center1_lon, center1_lat
        else:
            gha, dec = center2_lon, center2_lat
    
    if gha < 0:
        gha += 360
    
    h0 = 90 - radius / 60
    
    return gha, dec, h0

def rhumb_line_bearing(lon1, lat1, lon2, lat2):
    """Calculate rhumb line bearing from point 1 to point 2."""
    dlon = lon1 - lon2
    if dlon > 180:
        dlon -= 360
    elif dlon < -180:
        dlon += 360
    
    lat1_r = np.radians(lat1)
    lat2_r = np.radians(lat2)
    
    dphi = np.log(np.tan(np.pi/4 + lat2_r/2) / np.tan(np.pi/4 + lat1_r/2))
    
    bearing = np.degrees(np.arctan2(np.radians(dlon), dphi))
    return bearing % 360

def half_distance(lon1, lat1, lon2, lat2, bearing):
    """Calculate half the distance between two objects in NM."""
    cos_bearing = np.cos(np.radians(bearing))
    
    if abs(cos_bearing) > 0.0001:
        a = abs(lat2 - lat1) * 60 / abs(cos_bearing)
    else:
        sin_bearing = np.sin(np.radians(bearing))
        avg_lat = (lat1 + lat2) / 2
        a = abs(lon2 - lon1) * np.cos(np.radians(avg_lat)) * 60 / abs(sin_bearing)
    
    return a / 2
```

### Example from Paper:

```
DR Position: 117°41'W, 33°27'N

Observation 1 - Range:
  Object: East End Santa Catalina Island (118°20.0'W, 33°18.5'N)
  Range: 31.6 NM
  Equivalent: GHA=118°20.0', Dec=33°18.5'N, H₀=89°28.4'

Observation 2 - Bearing:
  Object: Santiago Peak (117°31.9'W, 33°42.5'N)
  Bearing: 28.5° true
  Equivalent: GHA=224°21.3', Dec=23°27.2'N, H₀=0°00.0'

Observation 3 - Horizontal Angle:
  Objects: San Onofre (117°33.5'W, 33°22.5'N) and Santiago Peak
  Angle: 102°
  Equivalent: GHA=117°30.2', Dec=33°32.3'N, H₀=89°49.8'

Computed Position: 117°41.8'W, 33°26.4'N
```

### Key Quotes:

- "Since these observations yield circles of position, just like the celestial observations, the same algorithm used to solve for a celestial fix can be used to solve for a coastal fix, or for a combination of both."
- "Once the coastal observation has been input, it can be treated exactly as a celestial observation would be: it can be advanced, retarded, and combined with any other data to yield a fix or running fix."
- "It is not necessary (and not appropriate) to correct the altitude in any way for such effects as dip and refraction."

### Accuracy Considerations:

- Plane geometry valid for distances < 100 NM
- Error in LOP < 50 m for this range
- Near poles: use great circle translation instead of rhumb line
- Radio bearings: require bearing correction for long distances

### Relevance to Research:
- **Algorithm unification:** Single fix algorithm handles all observation types
- **Practical value:** Coastal/celestial combined fixes
- **Same authors:** Continues Metcalf (1991) overdetermined fix work
- **Code reuse:** One COP intersection algorithm for everything
- **Running fix compatible:** Equivalent observations can be advanced/retarded
- **Python implementation:** Direct translation of formulas

### Limitations:
- Plane geometry limits: ~100 NM max object distance
- Horizontal angle yields partial circle (ambiguity possible)
- Near-pole special handling needed

---

## Article #22: Swaszek et al. (2019) - HDOP Optimization for Celestial Navigation

### Citation
Swaszek, P. F., Hartnett, R. J., & Seals, K. C. (2019). Rethinking star selection in celestial navigation. *Proceedings of the 2019 International Technical Meeting of The Institute of Navigation*, Reston, Virginia, January 2019, pp. 522-535. https://doi.org/10.33012/2019.16678

### Relevance to Research
- **Geometry optimization:** HDOP minimization for celestial body selection
- **Performance bounds:** Theoretical lower limits on position accuracy
- **Algorithm design:** Real-time star selection for optimal geometry
- **GNSS connection:** Shows altitude-intercept = first iteration of linearized least squares
- **Python implementation:** Complete HDOP calculation and star selection algorithms

### Core Mathematical Framework

#### Direction Cosine Matrix
For m celestial bodies at azimuths θ₁, θ₂, ..., θₘ:

$$\mathbf{G} = \begin{bmatrix} \sin\theta_1 & \cos\theta_1 \\ \sin\theta_2 & \cos\theta_2 \\ \vdots & \vdots \\ \sin\theta_m & \cos\theta_m \end{bmatrix} = \begin{bmatrix} e_1 & n_1 \\ e_2 & n_2 \\ \vdots & \vdots \\ e_m & n_m \end{bmatrix}$$

where $e_k = \sin\theta_k$ and $n_k = \cos\theta_k$ are East and North direction components.

#### Least Squares Position Solution
$$\begin{bmatrix} \delta_e \\ \delta_n \end{bmatrix} = (\mathbf{G}^T \mathbf{G})^{-1} \mathbf{G}^T \boldsymbol{\delta}$$

where $\boldsymbol{\delta}$ is the vector of altitude differences (measured - computed).

#### Position Error Covariance
$$\text{Cov}\begin{bmatrix} \delta_e \\ \delta_n \end{bmatrix} = \sigma^2 (\mathbf{G}^T \mathbf{G})^{-1}$$

### HDOP Formula

$$\text{HDOP} = \sqrt{\text{trace}\{(\mathbf{G}^T \mathbf{G})^{-1}\}}$$

Expanded form:

$$\text{HDOP} = \sqrt{\frac{m}{\left(\sum_{k=1}^{m} e_k^2\right)\left(\sum_{k=1}^{m} n_k^2\right) - \left(\sum_{k=1}^{m} e_k n_k\right)^2}}$$

### Theorem 1: Lower Bound to HDOP

$$\boxed{\text{HDOP} \geq \sqrt{\frac{4}{m}}}$$

This minimum is achieved when the **balance conditions** are satisfied:

$$\sum_{k=1}^{m} e_k n_k = 0 \quad \text{and} \quad \sum_{k=1}^{m} e_k^2 = \sum_{k=1}^{m} n_k^2 = \frac{m}{2}$$

### Balanced Constellation Properties

**Theorem 2:** Rotating a balanced constellation by angle φ remains balanced.

**Theorem 3:** For m ≥ 3, a regular m-gon (azimuths evenly distributed over 360°) is balanced.

**Theorem 4:** The union of two balanced constellations yields another balanced constellation.

**Implication:** Many configurations achieve minimum HDOP, not just evenly-spaced azimuths.

Examples of balanced 8-star configurations:
- 8 stars at 45° intervals
- 4 + 4: Two groups of 4, each at 90° intervals, arbitrary relative rotation
- 5 + 3: Pentagon + triangle, arbitrary rotations

### HDOP Values for Different Star Counts

| Stars (m) | Lower Bound HDOP |
|-----------|------------------|
| 2         | 1.414           |
| 3         | 1.155           |
| 4         | 1.000           |
| 5         | 0.894           |
| 6         | 0.816           |
| 7         | 0.756           |
| 8         | 0.707           |
| 10        | 0.632           |

### Python Implementation

```python
import numpy as np
from itertools import combinations

def compute_hdop(azimuths_deg):
    """
    Compute HDOP for celestial navigation given star azimuths.
    
    Parameters:
    - azimuths_deg: array of azimuth angles in degrees
    
    Returns:
    - HDOP value
    """
    azimuths = np.radians(azimuths_deg)
    m = len(azimuths)
    
    # Direction cosine components
    e = np.sin(azimuths)  # East components
    n = np.cos(azimuths)  # North components
    
    # Build G matrix
    G = np.column_stack([e, n])
    
    # Compute (G^T G)^-1
    GTG = G.T @ G
    GTG_inv = np.linalg.inv(GTG)
    
    # HDOP is square root of trace
    hdop = np.sqrt(np.trace(GTG_inv))
    
    return hdop

def compute_hdop_direct(azimuths_deg):
    """
    Compute HDOP using closed-form formula (more efficient).
    
    Parameters:
    - azimuths_deg: array of azimuth angles in degrees
    
    Returns:
    - HDOP value
    """
    azimuths = np.radians(azimuths_deg)
    m = len(azimuths)
    
    e = np.sin(azimuths)
    n = np.cos(azimuths)
    
    sum_e2 = np.sum(e**2)
    sum_n2 = np.sum(n**2)
    sum_en = np.sum(e * n)
    
    denominator = sum_e2 * sum_n2 - sum_en**2
    
    if denominator <= 0:
        return np.inf  # Degenerate geometry
    
    hdop = np.sqrt(m / denominator)
    
    return hdop

def hdop_lower_bound(m):
    """
    Theoretical HDOP lower bound for m celestial bodies.
    
    Parameters:
    - m: number of celestial bodies
    
    Returns:
    - Lower bound HDOP value
    """
    return np.sqrt(4.0 / m)

def check_balance(azimuths_deg, tolerance=0.01):
    """
    Check if constellation is balanced (achieves minimum HDOP).
    
    Parameters:
    - azimuths_deg: array of azimuth angles in degrees
    - tolerance: relative tolerance for balance check
    
    Returns:
    - dict with balance metrics and status
    """
    azimuths = np.radians(azimuths_deg)
    m = len(azimuths)
    
    e = np.sin(azimuths)
    n = np.cos(azimuths)
    
    # Balance conditions
    sum_en = np.sum(e * n)           # Should be 0
    sum_e2 = np.sum(e**2)            # Should be m/2
    sum_n2 = np.sum(n**2)            # Should be m/2
    
    is_balanced = (
        abs(sum_en) < tolerance * m and
        abs(sum_e2 - m/2) < tolerance * m and
        abs(sum_n2 - m/2) < tolerance * m
    )
    
    return {
        'sum_en': sum_en,
        'sum_e2': sum_e2,
        'sum_n2': sum_n2,
        'target_e2_n2': m / 2,
        'is_balanced': is_balanced,
        'hdop': compute_hdop_direct(azimuths_deg),
        'hdop_bound': hdop_lower_bound(m)
    }

def find_best_subset_exhaustive(azimuths_deg, m):
    """
    Find best m-star subset by exhaustive search (brute force).
    
    Parameters:
    - azimuths_deg: array of all available star azimuths
    - m: number of stars to select
    
    Returns:
    - best_azimuths: optimal subset
    - best_hdop: HDOP value achieved
    """
    best_hdop = np.inf
    best_subset = None
    
    for subset_indices in combinations(range(len(azimuths_deg)), m):
        subset = np.array([azimuths_deg[i] for i in subset_indices])
        hdop = compute_hdop_direct(subset)
        
        if hdop < best_hdop:
            best_hdop = hdop
            best_subset = subset
    
    return best_subset, best_hdop

def find_balanced_subset_fast(azimuths_deg, m):
    """
    Fast star selection algorithm using balance decomposition.
    Implements Swaszek et al. overlap-and-add method.
    
    Parameters:
    - azimuths_deg: array of all available star azimuths (integers 0-359)
    - m: number of stars to select
    
    Returns:
    - best_subset: near-optimal star azimuths
    - hdop: achieved HDOP
    """
    # Decomposition strategies for common m values
    decompositions = {
        3: [(3,)],
        4: [(4,)],
        5: [(5,)],
        6: [(6,), (3, 3)],
        7: [(7,), (4, 3)],
        8: [(8,), (4, 4), (5, 3)],
        9: [(9,), (6, 3), (5, 4)],
        10: [(10,), (5, 5), (7, 3), (6, 4)]
    }
    
    if m not in decompositions:
        # Default: try single polygon + complementary
        decompositions[m] = [(m,), (m-3, 3) if m > 3 else (m,)]
    
    best_hdop = np.inf
    best_subset = None
    
    for decomp in decompositions.get(m, [(m,)]):
        subset = find_decomposed_subset(azimuths_deg, decomp)
        if subset is not None and len(subset) == m:
            hdop = compute_hdop_direct(subset)
            if hdop < best_hdop:
                best_hdop = hdop
                best_subset = subset
    
    return best_subset, best_hdop

def find_decomposed_subset(azimuths_deg, decomposition):
    """
    Find stars matching a decomposition pattern (e.g., 4+4 or 5+3).
    
    Parameters:
    - azimuths_deg: available star azimuths
    - decomposition: tuple of subset sizes (e.g., (4, 4) or (5, 3))
    
    Returns:
    - selected star azimuths or None if not found
    """
    available = set(int(a) % 360 for a in azimuths_deg)
    selected = []
    
    for n_stars in decomposition:
        polygon_stars = find_near_regular_polygon(available, n_stars)
        if polygon_stars is None:
            return None
        selected.extend(polygon_stars)
        available -= set(polygon_stars)
    
    return np.array(selected)

def find_near_regular_polygon(available_azimuths, n_stars):
    """
    Find n_stars from available that most closely form a regular polygon.
    Uses overlap-and-add approach from Swaszek et al.
    
    Parameters:
    - available_azimuths: set of available azimuth values
    - n_stars: number of stars needed
    
    Returns:
    - list of azimuth values or None
    """
    sector_size = 360 // n_stars
    
    # Create indicator vector
    x = np.zeros(360, dtype=int)
    for a in available_azimuths:
        x[a % 360] = 1
    
    # Iteratively widen until n_stars found in balance
    for width in range(1, 180 // n_stars):
        if width > 1:
            # Widen by one degree alternating left/right
            x_new = np.zeros_like(x)
            for i in range(360):
                x_new[i] = 1 if any(x[(i + d) % 360] for d in range(-width//2, width//2 + 1)) else 0
            x = x_new
        
        # Subdivide into n_stars sectors and sum
        y = np.zeros(sector_size, dtype=int)
        for k in range(n_stars):
            y += x[k * sector_size:(k + 1) * sector_size]
        
        # Check if any offset achieves n_stars
        for offset in range(sector_size):
            if y[offset] == n_stars:
                # Found! Extract the actual azimuths
                stars = []
                for k in range(n_stars):
                    target = (k * sector_size + offset) % 360
                    # Find closest available azimuth
                    for a in available_azimuths:
                        if abs(a - target) <= width or abs(a - target - 360) <= width:
                            if a not in stars:
                                stars.append(a)
                                break
                if len(stars) == n_stars:
                    return stars
    
    return None

def position_error_ellipse(azimuths_deg, sigma_nm=1.0):
    """
    Compute 95% confidence ellipse parameters for position fix.
    
    Parameters:
    - azimuths_deg: star azimuths used in fix
    - sigma_nm: measurement standard deviation in nautical miles
    
    Returns:
    - dict with ellipse parameters
    """
    azimuths = np.radians(azimuths_deg)
    e = np.sin(azimuths)
    n = np.cos(azimuths)
    
    G = np.column_stack([e, n])
    GTG_inv = np.linalg.inv(G.T @ G)
    
    # Covariance matrix
    cov = sigma_nm**2 * GTG_inv
    
    # Eigenvalue decomposition for ellipse
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    
    # 95% confidence (chi-square with 2 DOF, p=0.05: 5.991)
    chi2_95 = 5.991
    
    # Semi-axes
    a = np.sqrt(chi2_95 * eigenvalues[1])  # Major
    b = np.sqrt(chi2_95 * eigenvalues[0])  # Minor
    
    # Orientation (angle of major axis from North)
    angle = np.degrees(np.arctan2(eigenvectors[0, 1], eigenvectors[1, 1]))
    
    return {
        'semi_major_nm': a,
        'semi_minor_nm': b,
        'orientation_deg': angle,
        'cov_ee': cov[0, 0],
        'cov_nn': cov[1, 1],
        'cov_en': cov[0, 1],
        'hdop': np.sqrt(np.trace(GTG_inv))
    }
```

### Example from Paper (2019 Nautical Almanac Data)

```
Location: 39°N, 74°W (off coast of Wildwood, NJ)
Time: Twilight, January 30, 2019
Available stars: M = 29

Star azimuths:
[1, 16, 25, 31, 64, 76, 95, 97, 114, 122,
 122, 125, 125, 135, 172, 178, 205, 209, 212, 219,
 227, 260, 265, 269, 309, 319, 322, 332, 360]

Goal: Select m = 8 stars for optimal HDOP

Results:
- Lower bound: HDOP = √(4/8) = 0.7071
- Best 4+4 selection: [1, 31, 95, 125, 178, 219, 269, 309]
  HDOP = 0.7075 (within 0.06% of bound)
- Best 3+5 selection: [31, 95, 97, 172, 212, 219, 309, 332]
  HDOP = 0.7087
- Best evenly-spaced attempt: [31, 76, 125, 172, 219, 265, 309, 360]
  HDOP = 0.7079

Worst possible: [122, 122, 125, 125, 135, 309, 319, 322]
  HDOP = 2.802 (collinear bodies)
```

### Key Quotes

- "The altitude-intercept method combined with the 'keen eye' of the observer are equivalent to the first iteration in a typical least squares position solution."
- "A lower bound to the HDOP is HDOP ≥ √(4/m) and this value is achieved by m celestial objects whose azimuths satisfy two balance conditions."
- "These conditions do not, in fact, require that the m celestial observations be separated by 360/m degrees of azimuth in order to minimize HDOP."
- "No longer requiring an even spread around the azimuth circle provides added flexibility to a navigator, especially if there is lack of a visible horizon in some direction."
- "Fully 50% have HDOP below 0.7195 and 90% with HDOP below 0.752; hence, even choosing stars at random is likely to produce a pretty good result!"

### Critical Insight: Bowditch Guidance Reconsidered

Bowditch [15] recommends: "Choose the stars and planets that give the best bearing spread... select bodies with predicted altitude between 30° and 70°."

Swaszek et al. show:
1. Even spread is **sufficient** but **not necessary** for minimum HDOP
2. Many non-uniform configurations achieve the same optimal HDOP
3. Limited sky visibility (land, fog, lighting) does not preclude optimal results
4. The "best" configuration from their example was NOT evenly spaced

### Relevance to Research Questions

- **SRQ-1 (Accuracy):** HDOP provides theoretical accuracy limits; √(4/m) × σ for m bodies
- **SRQ-2 (Multi-body methods):** Balance conditions for optimal geometry
- **SRQ-4 (Confidence intervals):** Error ellipse computation from covariance
- **H2 (3+ bodies improve accuracy):** Quantified via 1/√m improvement factor
- **H3 (Error modeling):** Position covariance matrix formulation

### Comparison with Other Articles

| Aspect | Swaszek 2019 | Nguyen 2014 | Kaplan 1995 |
|--------|--------------|-------------|-------------|
| Focus | Star selection | Position algorithm | Moving observer |
| Optimization | HDOP minimization | Least squares | Least squares |
| Geometry | Balance conditions | Not addressed | Not addressed |
| Error metric | HDOP, ellipse | Residual | Covariance |
| Pre-observation | Yes (planning) | No (post-processing) | No |

### Limitations

- Assumes equal measurement variances (σ² common to all)
- Does not account for altitude-dependent refraction errors
- No bias term included (could add column of ones to G)
- 2-D horizontal position only (no altitude estimation)

### Future Extensions Noted by Authors

- Weighted balance for varying measurement accuracy
- Bias term estimation (systematic sextant errors)
- Altitude estimation capability

---

## Article #23: Nguyen Thai Duong & Luong Tu Nam (2021) - Altitude Difference Least Squares

### Citation
Nguyen, T. D., & Luong, T. N. (2021). Determining ship's position by the celestial altitude difference using the least squares method. *Journal of Hunan University (Natural Sciences)*, 48(5), 199-207.

### Relevance to Research
- **Grid-search least squares:** Practical brute-force position search
- **Search domain construction:** Systematic approach based on DR error bounds
- **Shipboard validation:** Real experiment on M/V Thai Binh 05
- **Accuracy comparison:** 1.7 NM (proposed) vs 12.2 NM (Saint-Hilaire) vs GPS
- **Excel implementation:** Demonstrates practical applicability

### Core Mathematical Framework

#### Altitude Equation (Celestial Triangle)

$$\sin h = \sin\varphi \sin\delta + \cos\varphi \cos\delta \cos(t_G^\gamma + \lambda)$$

where:
- $h$: altitude of celestial body
- $\varphi$: observer's latitude
- $\lambda$: observer's longitude (West positive)
- $\delta$: declination of celestial body
- $t_G^\gamma$: Greenwich Hour Angle of Aries

#### System of Equations for Multi-Body Observation

$$\begin{cases}
\sin h_1 = \sin\varphi \sin\delta_1 + \cos\varphi \cos\delta_1 \cos(t_{G1} + \lambda) \\
\sin h_2 = \sin\varphi \sin\delta_2 + \cos\varphi \cos\delta_2 \cos(t_{G2} + \lambda) \\
\vdots \\
\sin h_n = \sin\varphi \sin\delta_n + \cos\varphi \cos\delta_n \cos(t_{Gn} + \lambda)
\end{cases}$$

### Algorithm: Grid Search Least Squares

#### Step 1: Define Search Domain

Based on estimated position $M_C(\varphi_C, \lambda_C)$ and estimated error:

$$\varphi_{min} \leq \varphi \leq \varphi_{max}$$
$$\lambda_{min} \leq \lambda \leq \lambda_{max}$$

where bounds are set to ensure >95% probability of containing true position.

**Radius of Mean Square Error:**
$$R = \sqrt{(S \cdot \Delta_L)^2 + (S \cdot \Delta_{TK})^2}$$

where:
- $S$: distance traveled
- $\Delta_L$: compass error
- $\Delta_{TK}$: tachometer (speed log) error

#### Step 2: Establish Grid of Assumed Positions

Create set $A = \{F_{xy}\}$ with grid spacing:
$$\varphi_{i+1} - \varphi_i = \lambda_{i+1} - \lambda_i = 0.000001° \approx 0.11 \text{ m}$$

#### Step 3: Calculate Altitude Variation for Each Grid Point

For each grid point $F_{xy}(\varphi_x, \lambda_y)$:

**Calculated altitude:**
$$h_{Cxy} = \sin^{-1}[\sin\varphi_x \sin\delta + \cos\varphi_x \cos\delta \cos(t_G + \lambda_y)]$$

**Altitude variation:**
$$\Delta h_{xy} = h_{Cxy} - h_O$$

where $h_O$ is the observed (measured) altitude.

#### Step 4: Find Single Probable Position per Body

For celestial body $C_1$, find position minimizing squared altitude variation:

$$F_{m_1n_1}: \quad (\Delta h_1)^2 = \min_{m,n}(\Delta h_{xy})^2 = \min_{m,n}(h_{Cxy} - h_{O1})^2$$

#### Step 5: Compute Most Probable Position (Multi-Body)

For $k$ celestial bodies with individual probable positions $F_{m_1n_1}, F_{m_2n_2}, \ldots, F_{m_kn_k}$:

$$\boxed{\varphi = \frac{\sum_{i=1}^{k} \varphi_{m_i}}{k} \quad \text{and} \quad \lambda = \frac{\sum_{i=1}^{k} \lambda_{n_i}}{k}}$$

### Python Implementation

```python
import numpy as np

def altitude_equation(lat, lon, dec, gha):
    """
    Calculate computed altitude using celestial triangle formula.
    
    Parameters:
    - lat: observer latitude (degrees, N positive)
    - lon: observer longitude (degrees, W positive)
    - dec: celestial body declination (degrees)
    - gha: Greenwich Hour Angle (degrees)
    
    Returns:
    - altitude in degrees
    """
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    dec_r = np.radians(dec)
    gha_r = np.radians(gha)
    
    # Local Hour Angle
    lha = gha_r + lon_r  # For West longitude
    
    sin_h = (np.sin(lat_r) * np.sin(dec_r) + 
             np.cos(lat_r) * np.cos(dec_r) * np.cos(lha))
    
    # Clamp to valid range
    sin_h = np.clip(sin_h, -1.0, 1.0)
    
    return np.degrees(np.arcsin(sin_h))

def compute_search_domain(dr_lat, dr_lon, compass_error_deg, 
                          speed_error_pct, distance_nm, confidence=0.95):
    """
    Compute search domain boundaries based on DR error analysis.
    
    Parameters:
    - dr_lat, dr_lon: dead reckoning position
    - compass_error_deg: compass adjustment error in degrees
    - speed_error_pct: speed log error as percentage (e.g., 0.02 for 2%)
    - distance_nm: distance traveled since last fix in NM
    - confidence: probability of containing true position (default 95%)
    
    Returns:
    - (lat_min, lat_max, lon_min, lon_max) in degrees
    """
    # Radius of mean square error (in NM)
    compass_error_nm = distance_nm * np.sin(np.radians(compass_error_deg))
    speed_error_nm = distance_nm * speed_error_pct
    R = np.sqrt(compass_error_nm**2 + speed_error_nm**2)
    
    # Scale for confidence level (95% -> ~2 sigma)
    if confidence >= 0.95:
        R *= 2.0
    
    # Convert NM to degrees
    R_lat = R / 60.0  # 1 NM = 1/60 degree latitude
    R_lon = R / (60.0 * np.cos(np.radians(dr_lat)))  # Adjust for latitude
    
    return (dr_lat - R_lat, dr_lat + R_lat, 
            dr_lon - R_lon, dr_lon + R_lon)

def grid_search_position(observations, dr_lat, dr_lon, 
                         search_radius_nm=10.0, grid_spacing_nm=0.1):
    """
    Find ship position using grid search least squares method.
    
    Parameters:
    - observations: list of dicts with keys 'ho', 'dec', 'gha'
      - ho: observed altitude (degrees)
      - dec: declination (degrees)
      - gha: Greenwich Hour Angle (degrees)
    - dr_lat, dr_lon: dead reckoning position
    - search_radius_nm: radius of search domain in NM
    - grid_spacing_nm: grid point spacing in NM
    
    Returns:
    - (lat, lon): most probable position
    - residuals: altitude residuals at solution
    """
    # Create search grid
    lat_range = search_radius_nm / 60.0  # degrees
    lon_range = search_radius_nm / (60.0 * np.cos(np.radians(dr_lat)))
    
    lat_spacing = grid_spacing_nm / 60.0
    lon_spacing = grid_spacing_nm / (60.0 * np.cos(np.radians(dr_lat)))
    
    lats = np.arange(dr_lat - lat_range, dr_lat + lat_range, lat_spacing)
    lons = np.arange(dr_lon - lon_range, dr_lon + lon_range, lon_spacing)
    
    # Find best position for each observation
    individual_positions = []
    
    for obs in observations:
        min_residual_sq = np.inf
        best_lat, best_lon = dr_lat, dr_lon
        
        for lat in lats:
            for lon in lons:
                h_calc = altitude_equation(lat, lon, obs['dec'], obs['gha'])
                residual_sq = (h_calc - obs['ho'])**2
                
                if residual_sq < min_residual_sq:
                    min_residual_sq = residual_sq
                    best_lat, best_lon = lat, lon
        
        individual_positions.append((best_lat, best_lon))
    
    # Average individual positions (least squares combination)
    final_lat = np.mean([p[0] for p in individual_positions])
    final_lon = np.mean([p[1] for p in individual_positions])
    
    # Compute residuals at final position
    residuals = []
    for obs in observations:
        h_calc = altitude_equation(final_lat, final_lon, obs['dec'], obs['gha'])
        residuals.append(h_calc - obs['ho'])
    
    return (final_lat, final_lon), residuals

def grid_search_combined(observations, dr_lat, dr_lon,
                         search_radius_nm=10.0, grid_spacing_nm=0.1):
    """
    Alternative: Find position minimizing total squared residuals.
    
    Parameters: (same as grid_search_position)
    
    Returns:
    - (lat, lon): position minimizing sum of squared residuals
    - total_residual_sq: sum of squared altitude differences
    """
    lat_range = search_radius_nm / 60.0
    lon_range = search_radius_nm / (60.0 * np.cos(np.radians(dr_lat)))
    
    lat_spacing = grid_spacing_nm / 60.0
    lon_spacing = grid_spacing_nm / (60.0 * np.cos(np.radians(dr_lat)))
    
    lats = np.arange(dr_lat - lat_range, dr_lat + lat_range, lat_spacing)
    lons = np.arange(dr_lon - lon_range, dr_lon + lon_range, lon_spacing)
    
    min_total_sq = np.inf
    best_lat, best_lon = dr_lat, dr_lon
    
    for lat in lats:
        for lon in lons:
            total_sq = 0.0
            for obs in observations:
                h_calc = altitude_equation(lat, lon, obs['dec'], obs['gha'])
                total_sq += (h_calc - obs['ho'])**2
            
            if total_sq < min_total_sq:
                min_total_sq = total_sq
                best_lat, best_lon = lat, lon
    
    return (best_lat, best_lon), min_total_sq

def grid_search_vectorized(observations, dr_lat, dr_lon,
                           search_radius_nm=10.0, grid_spacing_nm=0.1):
    """
    Vectorized grid search for improved performance.
    
    Uses NumPy broadcasting to evaluate all grid points simultaneously.
    """
    lat_range = search_radius_nm / 60.0
    lon_range = search_radius_nm / (60.0 * np.cos(np.radians(dr_lat)))
    
    lat_spacing = grid_spacing_nm / 60.0
    lon_spacing = grid_spacing_nm / (60.0 * np.cos(np.radians(dr_lat)))
    
    lats = np.arange(dr_lat - lat_range, dr_lat + lat_range, lat_spacing)
    lons = np.arange(dr_lon - lon_range, dr_lon + lon_range, lon_spacing)
    
    # Create 2D grid
    LAT, LON = np.meshgrid(lats, lons, indexing='ij')
    
    # Compute total squared residual for each grid point
    total_sq = np.zeros_like(LAT)
    
    for obs in observations:
        LAT_r = np.radians(LAT)
        LON_r = np.radians(LON)
        dec_r = np.radians(obs['dec'])
        gha_r = np.radians(obs['gha'])
        
        lha = gha_r + LON_r
        
        sin_h = (np.sin(LAT_r) * np.sin(dec_r) + 
                 np.cos(LAT_r) * np.cos(dec_r) * np.cos(lha))
        sin_h = np.clip(sin_h, -1.0, 1.0)
        h_calc = np.degrees(np.arcsin(sin_h))
        
        total_sq += (h_calc - obs['ho'])**2
    
    # Find minimum
    min_idx = np.unravel_index(np.argmin(total_sq), total_sq.shape)
    best_lat = LAT[min_idx]
    best_lon = LON[min_idx]
    
    return (best_lat, best_lon), total_sq[min_idx]
```

### Shipboard Experiment Results (M/V Thai Binh 05)

**Test Conditions:**
- Date: June 16, 2020
- Location: Gulf of Tonkin
- Speed: 9 knots, Course: 155°
- Time: 21:15:00 GMT (nautical twilight)

**Observations:**

| Star | GMT | H₀ (observed) | GHA | Dec |
|------|-----|---------------|-----|-----|
| Vega | 21:15:00 | 36°12.8' | 313°35.1' | 38°47.7'N |
| Mirfak | 21:18:00 | 27°29.2' | 181°08.2' | 50°05.0'N |
| Fomalhaut | 21:21:00 | 35°01.5' | 225°05.5' | 39°55.0'S |

**Results Comparison:**

| Method | Latitude | Longitude | Error vs GPS |
|--------|----------|-----------|---------------|
| GPS (reference) | 20°39.6'N | 106°59.9'E | — |
| Least Squares (proposed) | 20°38.9'N | 107°01.5'E | **1.7 NM** |
| Saint-Hilaire (traditional) | 20°50.5'N | 107°05.7'E | 12.2 NM |
| Dead Reckoning | 20°36.5'N | 106°55.5'E | 5.2 NM |

### Comparison with Other Methods

**4-Star Test (Analytical):**

| Method Comparison | Difference |
|-------------------|------------|
| Saint-Hilaire vs Genetic Algorithm | 0.9 NM / 064° |
| Saint-Hilaire vs Least Squares | 0.8 NM / 255° |
| Genetic Algorithm vs Least Squares | 0.2 NM / 202° |

### Key Advantages Claimed

1. **Single-body capable:** Can work with one celestial body observation
2. **DR-independent:** Less sensitive to dead reckoning errors than traditional method
3. **Multi-body improvement:** Accuracy improves with additional observations
4. **Simple implementation:** Excel-based, suitable for shipboard use
5. **Fast computation:** Search domain limits reduce computation time

### Key Quotes

- "The ship's position is the position with the least squares of the altitude variations from the simultaneous observation of celestial bodies."
- "This method satisfies the requirements for a backup method of ship positioning for ECDIS according to the amendments of the STCW Convention 78/2010."
- "The difference between the position by the Least Squares method and GPS position is 1.7 nm/course 255°."
- "The program of determining ship's position by the celestial altitude difference using the Least Squares Method is constructed using Excel with an easy-to-understand algorithm and simple user-oriented interface."

### Comparison with Other Articles

| Aspect | Nguyen 2021 | Nguyen 2014 | Kaplan 1995 |
|--------|-------------|-------------|-------------|
| Method | Grid search | Matrix least squares | Moving observer LS |
| Approach | Brute force | Iterative | Iterative |
| Validated | Yes (shipboard) | Simulation | Simulation |
| Accuracy | 1.7 NM vs GPS | Sub-NM claimed | 0.41 NM simulation |
| Implementation | Excel | Not specified | Not specified |

### Relevance to Research Questions

- **SRQ-1 (Accuracy):** Real-world validation: 1.7 NM vs GPS
- **SRQ-2 (Multi-body):** Grid search with averaging of individual body solutions
- **SRQ-3 (Uncertainty):** Search domain based on DR error bounds (>95% confidence)
- **H1 (Intercept equivalence):** Provides alternative to intercept method
- **H3 (Error modeling):** Compass and speed log error propagation

### Limitations

- Computationally intensive for fine grids
- Averaging individual positions may not be optimal (vs. combined minimization)
- Single observation case requires additional assumptions
- No formal covariance analysis

### Extensions for Python Implementation

1. **Hierarchical search:** Coarse grid → fine grid refinement
2. **Parallel processing:** Grid evaluation is embarrassingly parallel
3. **Combined objective:** Minimize total Σ(Δh)² instead of averaging
4. **Weighted average:** Use individual residuals as inverse weights

---

## Article #24: Barazzetti (2025) - Astronomical Azimuth with Python/Skyfield

### Citation
Barazzetti, L. (2025). Revitalizing astronomical azimuth determination: Integrating modern computing with traditional techniques. *Sensors*, 25(6), 1871. https://doi.org/10.3390/s25061871

### Relevance to Research
- **Python + Skyfield implementation:** Directly applicable to project's target stack
- **Azimuth calculation:** Complete formulas with error propagation
- **Field validation:** ±1-2 arcseconds accuracy demonstrated
- **Multiple celestial bodies:** Stars, planets (Venus, Mars, Jupiter, Saturn), Moon
- **Open-source approach:** Matches project philosophy
- **Error simulation:** Pre-observation planning capability

### Core Mathematical Framework

#### Azimuth Calculation Formula

$$A = \arctan\left(\frac{\sin H}{\cos H \sin\Phi - \tan\delta \cos\Phi}\right)$$

where:
- $A$: azimuth of celestial body
- $H$: hour angle
- $\Phi$: observer's astronomical latitude
- $\delta$: declination of celestial body

#### Hour Angle Calculation

$$H = \text{LST} - \text{RA}$$

where:
- $\text{LST}$: Local Sidereal Time
- $\text{RA}$: Right Ascension of celestial body

#### Astronomical Azimuth from Observations

Using total station measurements (Figure 3 in paper):

$$A_{SP} = A + \beta_P - \beta_{CB}$$

where:
- $A_{SP}$: astronomical azimuth from station S to target P
- $A$: computed azimuth of celestial body
- $\beta_P$: horizontal circle reading to terrestrial target
- $\beta_{CB}$: horizontal circle reading to celestial body

#### Laplace Equation (Geodetic to Astronomical Azimuth)

$$A_G = \alpha + \left(\eta \tan\varphi + (\xi\sin\alpha - \eta\cos\alpha)\cot\zeta\right)$$

where:
- $A_G$: astronomical azimuth
- $\alpha$: geodetic azimuth
- $(\xi, \eta)$: deflection of the vertical components
- $\varphi$: geodetic latitude
- $\zeta$: zenith angle

### Error Propagation Analysis

Azimuth error from latitude uncertainty:
$$\sigma_{A\Phi} = \frac{\partial A}{\partial\Phi} \cdot \sigma_\Phi$$

Azimuth error from hour angle uncertainty:
$$\sigma_{AH} = \frac{\partial A}{\partial H} \cdot \sigma_H$$

**Key Finding:** Polaris shows minimal hour angle sensitivity (~1" error max), while stars like Dubhe show up to 60-80" error due to faster apparent motion.

### Python Implementation (Skyfield-based)

```python
import numpy as np
from skyfield.api import load, Topos
from skyfield.data import hipparcos

def compute_azimuth_altitude_star(hip_id, lat_deg, lon_deg, 
                                   utc_datetime, temperature_C=10, 
                                   pressure_mbar=1013.25):
    """
    Compute azimuth and altitude of a star using Skyfield.
    
    Parameters:
    - hip_id: Hipparcos catalog number
    - lat_deg, lon_deg: Observer position (degrees)
    - utc_datetime: datetime object in UTC
    - temperature_C: Temperature for refraction correction
    - pressure_mbar: Pressure for refraction correction
    
    Returns:
    - azimuth_deg: Azimuth from north (degrees)
    - altitude_deg: Altitude above horizon (degrees)
    """
    # Load ephemeris and star catalog
    eph = load('de421.bsp')
    earth = eph['earth']
    
    # Load Hipparcos catalog
    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)
    
    star = df.loc[hip_id]
    
    # Create observer location
    observer = earth + Topos(latitude_degrees=lat_deg,
                             longitude_degrees=lon_deg)
    
    # Create timescale and time
    ts = load.timescale()
    t = ts.utc(utc_datetime.year, utc_datetime.month, utc_datetime.day,
               utc_datetime.hour, utc_datetime.minute, utc_datetime.second)
    
    # Compute apparent position
    from skyfield.api import Star
    star_obj = Star.from_dataframe(star)
    astrometric = observer.at(t).observe(star_obj)
    apparent = astrometric.apparent()
    
    # Get altitude and azimuth
    alt, az, distance = apparent.altaz(temperature_C=temperature_C,
                                        pressure_mbar=pressure_mbar)
    
    return az.degrees, alt.degrees

def compute_azimuth_planet(planet_name, lat_deg, lon_deg, utc_datetime):
    """
    Compute azimuth and altitude of a planet.
    
    Parameters:
    - planet_name: 'venus', 'mars', 'jupiter barycenter', 'saturn barycenter'
    - lat_deg, lon_deg: Observer position (degrees)
    - utc_datetime: datetime object in UTC
    
    Returns:
    - azimuth_deg, altitude_deg
    """
    eph = load('de421.bsp')
    earth = eph['earth']
    planet = eph[planet_name]
    
    observer = earth + Topos(latitude_degrees=lat_deg,
                             longitude_degrees=lon_deg)
    
    ts = load.timescale()
    t = ts.utc(utc_datetime.year, utc_datetime.month, utc_datetime.day,
               utc_datetime.hour, utc_datetime.minute, utc_datetime.second)
    
    astrometric = observer.at(t).observe(planet)
    apparent = astrometric.apparent()
    alt, az, distance = apparent.altaz()
    
    return az.degrees, alt.degrees

def compute_azimuth_moon(lat_deg, lon_deg, utc_datetime):
    """
    Compute azimuth of Moon with limb calculations.
    
    Returns:
    - center_az: Azimuth to Moon center
    - east_limb_az: Azimuth to east limb
    - west_limb_az: Azimuth to west limb
    - angular_diameter_arcmin: Moon's angular diameter
    """
    eph = load('de421.bsp')
    earth = eph['earth']
    moon = eph['moon']
    
    observer = earth + Topos(latitude_degrees=lat_deg,
                             longitude_degrees=lon_deg)
    
    ts = load.timescale()
    t = ts.utc(utc_datetime.year, utc_datetime.month, utc_datetime.day,
               utc_datetime.hour, utc_datetime.minute, utc_datetime.second)
    
    astrometric = observer.at(t).observe(moon)
    apparent = astrometric.apparent()
    alt, az, distance = apparent.altaz()
    
    # Moon's mean radius: 1737.4 km
    moon_radius_km = 1737.4
    distance_km = distance.km
    angular_radius_rad = np.arctan(moon_radius_km / distance_km)
    angular_radius_deg = np.degrees(angular_radius_rad)
    angular_diameter_arcmin = 2 * angular_radius_deg * 60
    
    # East and west limb azimuths
    east_limb_az = az.degrees + angular_radius_deg
    west_limb_az = az.degrees - angular_radius_deg
    
    return az.degrees, east_limb_az, west_limb_az, angular_diameter_arcmin

def azimuth_error_propagation(star_hip_id, lat_deg, lon_deg, utc_datetime,
                               sigma_phi_arcsec=20, sigma_H_arcsec=20):
    """
    Compute azimuth error due to uncertainties in latitude and hour angle.
    
    Parameters:
    - star_hip_id: Hipparcos ID
    - lat_deg, lon_deg: Observer position
    - utc_datetime: Observation time
    - sigma_phi_arcsec: Latitude uncertainty (arcseconds)
    - sigma_H_arcsec: Hour angle uncertainty (arcseconds)
    
    Returns:
    - sigma_A_phi: Azimuth error from latitude uncertainty (arcsec)
    - sigma_A_H: Azimuth error from hour angle uncertainty (arcsec)
    """
    # Numerical differentiation step
    delta = 1.0 / 3600.0  # 1 arcsecond in degrees
    
    # Get nominal azimuth
    az_nom, _ = compute_azimuth_altitude_star(star_hip_id, lat_deg, lon_deg, 
                                               utc_datetime)
    
    # Partial derivative w.r.t. latitude
    az_phi_plus, _ = compute_azimuth_altitude_star(star_hip_id, lat_deg + delta, 
                                                    lon_deg, utc_datetime)
    dA_dPhi = (az_phi_plus - az_nom) / delta  # degrees/degree
    
    # Partial derivative w.r.t. hour angle (via longitude)
    az_lon_plus, _ = compute_azimuth_altitude_star(star_hip_id, lat_deg, 
                                                    lon_deg + delta, utc_datetime)
    dA_dH = (az_lon_plus - az_nom) / delta  # degrees/degree
    
    # Convert uncertainties
    sigma_phi_deg = sigma_phi_arcsec / 3600.0
    sigma_H_deg = sigma_H_arcsec / 3600.0
    
    # Propagate errors
    sigma_A_phi = abs(dA_dPhi * sigma_phi_deg) * 3600  # arcseconds
    sigma_A_H = abs(dA_dH * sigma_H_deg) * 3600  # arcseconds
    
    return sigma_A_phi, sigma_A_H

def compute_astronomical_azimuth(celestial_az, beta_target, beta_celestial):
    """
    Calculate astronomical azimuth from celestial body observation.
    
    Parameters:
    - celestial_az: Computed azimuth of celestial body (degrees)
    - beta_target: Horizontal circle reading to target (degrees)
    - beta_celestial: Horizontal circle reading to celestial body (degrees)
    
    Returns:
    - azimuth_SP: Astronomical azimuth from station to target (degrees)
    """
    azimuth_SP = celestial_az + beta_target - beta_celestial
    return azimuth_SP % 360

def simulate_star_precision(star_hip_id, latitude_deg, 
                            hour_angles_deg=None,
                            sigma_phi=20, sigma_H=20):
    """
    Simulate azimuth precision for a star across hour angles.
    
    Useful for planning observations.
    """
    if hour_angles_deg is None:
        hour_angles_deg = np.linspace(-90, 90, 181)
    
    # This would require more implementation for full simulation
    # Returns framework for planning
    pass
```

### Field Validation Results (Sondrio, Italy)

**Equipment:** Leica TS30 total station (0.5" angular precision)

**Reference:** Geodetic azimuth from static GPS corrected via Laplace equation:
$$A_G = 334°59'12.5''$$

**Results (9 August 2024):**

| Celestial Body | Observations | Computed Azimuth | Error vs Reference |
|----------------|--------------|------------------|--------------------|
| Polaris | 16 | 334°59'11.1" | **-1.4"** |
| Altair | 8 | 334°59'06.8" | -5.7" |
| Dubhe | 8 | 334°59'15.4" | +2.9" |
| Moon | 12 | 334°58'44.1" | -28.4" |

**Results (30 August 2024):**

| Celestial Body | Observations | Computed Azimuth | Error vs Reference |
|----------------|--------------|------------------|--------------------|
| Polaris | 8 | 334°59'11.7" | **-0.8"** |
| Altair | 8 | 334°59'10.7" | -1.8" |

**Polaris Consistency:** Maximum difference between campaigns = 0.5"

### Apparent Motion Rates

Critical for timing precision requirements:

| Celestial Body | Azimuth Rate |
|----------------|---------------|
| Polaris | ~0.2"/s |
| Dubhe | ~4.5"/s |
| Moon | ~11.1"/s |
| Altair | ~16.4"/s |

**Implication:** Polaris requires less precise timing; 1s error ≈ 0.2" azimuth error.

### Error Analysis: Polaris vs Other Stars

**Polaris (near celestial pole):**
- Hour angle error: <1" maximum across all latitudes
- Latitude error: ~2" maximum at extreme hour angles
- Minimal sensitivity to timing errors

**Dubhe (far from pole):**
- Hour angle error: Up to 60" near 0° hour angle at high latitudes
- Latitude error: Up to 80" maximum
- High sensitivity to timing and position errors

### Key Implementation Notes

1. **Skyfield Library:** Uses JPL DE421 ephemeris for high-precision positions
2. **Hipparcos Catalog:** Milliarcsecond precision for star positions
3. **Atmospheric Refraction:** Applied to altitude (not azimuth)
4. **Topocentric Correction:** Important for Moon and planets (parallax)
5. **UTC Timing:** Critical - use NTP servers or GPS time

### Deflection of Vertical Correction

When using GPS for position instead of astronomical determination:

$$\Phi = \varphi + \xi$$
$$\Lambda = \lambda + \eta / \cos\varphi$$

where $(\xi, \eta)$ are deflection components from geoid model.

**Resources:**
- International Service for the Geoid (ISG): https://www.isgeoid.polimi.it/
- ICGEM: https://icgem.gfz-potsdam.de/

### Key Quotes

- "Astronomical methods for azimuth determination continue to rank among the most accurate techniques for establishing orientation."
- "Integration of astronomical measurements with modern computing codes allows surveyors to achieve azimuths with an accuracy of ±1–2 arcseconds."
- "Python's open-source nature fosters collaboration and innovation, making advanced astronomical tools accessible to a wider audience."
- "Accurate timing is essential for these measurements, as it significantly impacts the hour angle."
- "Polaris would, of course, always be the first choice under favorable weather conditions."

### Comparison with Other Articles

| Aspect | Barazzetti 2025 | Swaszek 2019 | Nguyen 2014 |
|--------|-----------------|--------------|-------------|
| Focus | Azimuth determination | Star selection (HDOP) | Position fixing |
| Implementation | Python/Skyfield | Theoretical | Theoretical |
| Validation | Field (total station) | Analysis | Simulation |
| Accuracy | ±1-2" azimuth | HDOP optimization | <1 NM position |
| Celestial bodies | Stars, planets, Moon | Stars | Stars |

### Relevance to Research Questions

- **SRQ-1 (Accuracy):** Field-validated ±1-2" azimuth accuracy
- **SRQ-3 (Uncertainty):** Error propagation methodology
- **H3 (Error modeling):** Complete error analysis framework
- **Technical:** Skyfield implementation directly applicable to project

### Limitations

- Azimuth only (not position fixing)
- Requires clear sky conditions
- Moon accuracy limited (~28" vs stars ~1-6")
- Sun not implemented (requires solar filter)

### Direct Applicability to Project

1. **Skyfield usage patterns:** Loading ephemeris, star catalog, planets, Moon
2. **Topocentric calculations:** Observer position handling
3. **Error propagation:** Numerical differentiation approach
4. **Refraction correction:** Built-in Skyfield methods
5. **Validation methodology:** Comparison with GPS-derived reference

---

## Pending Content Needs

### For Methodology Section:
- [x] Sight reduction mathematical models - Kotlarić 1976
- [x] Intercept method (Marcq St Hilaire) formulas - Kotlarić 1976
- [ ] GHA, LHA, declination computation algorithms
- [ ] Altitude correction procedures (refraction, dip, semi-diameter)

### For Technical Implementation:
- [x] Python libraries for ephemeris (Skyfield) - Barazzetti 2025
- [ ] Nautical Almanac data sources
- [ ] Validation datasets

### For Results/Discussion:
- [x] Modern accuracy standards for celestial navigation - Kotlarić 1976
- [x] Comparison benchmarks from contemporary studies - Kotlarić 1976

---

## Article 25: Kotlarić (1976) - Sight Reduction Tables vs Celestial Navigation Computer

### Citation
Kotlarić, S. (1976). Sight reduction — by tables or by celestial navigation computer? *International Hydrographic Review*, LIII(2), 123-156.

### Relevance to Research
**HIGH RELEVANCE** - Directly addresses computational efficiency comparison (SRQ4/H4) between different sight reduction methods. Provides historical baseline for comparing modern Python implementation efficiency.

### Key Content

#### 1. Computational Efficiency Benchmarks (SRQ4/H4)
**Critical finding for efficiency gap:**
- Traditional tabular methods (K21 tables): **"less than two minutes"** for complete sight reduction
- Early celestial navigation computers (Galaxy 1, Intercepter): **"two minutes or less"**
- Step count comparison:
  - Intercepter computer: 11 data entries required
  - K21 tables: 2 table lookups (main table + multiplication table)

**Quote:** "From the analysis of the procedures for computation of Altitude (Intercept) and Azimuth by these two celestial navigation computers it can be realised that the procedure by these mini-computers is neither shorter nor simpler than the procedure by Tables K21."

#### 2. Accuracy Standards
| Method | Altitude Accuracy | Azimuth Accuracy |
|--------|------------------|------------------|
| K21 Tables | **0.2'** (arc-minutes) | **0.2°** |
| H.O.249 | 1' | 1° |
| H.O.214 | 0.1' (with interpolation) | 0.1° |
| Galaxy 1 | 0.1' | 0.1° |
| Intercepter | 0.1' | 0.1° |

**Quote:** "Tables K21 are not limited to tabulated stars... but enable observers in all latitudes to use all celestial bodies given in the Nautical Almanac."

#### 3. Sight Reduction Procedure (Intercept Method)
Complete worked examples with step-by-step procedures:

**Sun observation example:**
```
Input: LAT 35°N, DEC -10°34.7', MA 40°W
From K21: HC=31°13.5', ID=-47.5, Z=132°W
Altitude correction: -27.5' (from multiplication table)
Final HC = 30°46.0'
Ho = 30°51.6'
Intercept = +5.6' (Toward)
AZ = N132.6°W (or Zn=227.4°)
```

**Moon observation example:**
```
Input: LAT 35°N, DEC +11°38.7', MA 9°E
From K21: HC=64°38.5', ID=+57, Z=N159°E
Altitude correction: +36.8'
Final HC = 65°15.3'
Ho = 65°27.3'
Intercept = +12.0' (Toward)
AZ = N158.3°E
```

#### 4. Method Comparison Table (from paper)
| Tables | Volumes | SUN He | SUN AZ | MOON He | MOON AZ |
|--------|---------|--------|--------|---------|---------|
| K21 | 3 | 30°46.0' | N132.6°W | 65°15.3' | N158.3°E |
| H.O.249 | 3 | 30°46' | 227.4° | 65°16' | 158.4° |
| K1 | 1 | 30°46.0' | N132.7°W | 65°15.3' | N158.4°E |
| H.O.214 | 9 | 30°45.9' | N132.7°W | 65°15.4' | N158.5°E |
| H.O.229 | 6 | 30°45.9' | 227.5° | 65°15.4' | 158.5° |

#### 5. Early Computer Specifications (1976)
**Galaxy 1:**
- Size: 134 × 315 × 143 mm
- Weight: 3.65 kg
- Power: 12-24V DC or 115-230V AC
- Price: $1,295 USD (1976)
- Features: Sight reduction, great circle, star identification

**Intercepter:**
- Size: 270 × 180 × 165 mm
- Weight: 2.7 kg
- Power: 10-24V DC or 115/230V AC
- Price: $1,250 CAD (1976)
- Accuracy: 0.1' displayed, internal accuracy higher

#### 6. Key Formulas Referenced
The paper references standard celestial navigation formulas:
- Computed altitude (He) from spherical trigonometry
- Azimuth angle (Z) conversion rules
- Intercept calculation: **a = Ho - He** (+ toward, - away)
- GHA/LHA/MA relationships

#### 7. Advantages of Tabular Methods Noted
1. No power supply required
2. No repairs needed ("Tables need no repairs")
3. Lower cost (130× cheaper than computers)
4. Portable between vessels
5. Familiar to all navigators
6. No training for special procedures

### Relevance to Python Implementation

**For SRQ4/H4 (Computational Efficiency):**
This paper establishes historical baseline:
- 1976 computers: ~2 minutes
- Manual tables: ~2 minutes
- **Modern Python with Skyfield**: 20-50 μs (Tsai 2022)
- **Improvement factor**: ~2,400,000× faster than 1976 methods

**For Algorithm Design:**
The step-by-step procedures inform Python function structure:
1. Extract ephemeris data (GHA, Dec) - now via Skyfield
2. Calculate meridian angle from assumed position
3. Compute altitude using spherical formula
4. Apply corrections for declination interpolation
5. Calculate intercept (Ho - He)
6. Determine azimuth with proper quadrant

### Key Quotes

> "Sight reduction... in two minutes or less" (p. 131)

> "Tabulation of altitude and azimuth is made for the value of meridian angle (MA) from 0° to 180°" (p. 126)

> "The Tables are much cheaper; the price is incomparably lower than that of a celestial navigation computer" (p. 152)

> "Tables need no repairs. A defect on the computer can be expected sometime" (p. 152)

### Mathematical Content for Implementation
- Altitude correction index (ID) system for interpolation
- Multiplication table for declination corrections
- Azimuth angle to true azimuth conversion rules
- Position line plotting geometry

---

## Article 26: Feldman et al. (1972) - Sight Reduction Using the Portable Sextant Computer System

### Citation
Feldman, S., Seidelmann, P. K., Stephenson, E. D., & Kells, H. C. (1972). Sight reduction using the portable sextant computer system. *Navigation: Journal of The Institute of Navigation*, 19(4), 317-321. https://doi.org/10.1002/j.2161-4296.1972.tb01701.x

### Relevance to Research
**HIGH RELEVANCE** - Provides core sight reduction formulas and historical computational efficiency baseline. Authors from U.S. Naval Observatory and Defense Mapping Agency. Directly addresses SRQ4/H4 with timing data.

### Key Content

#### 1. Computational Efficiency Historical Baseline (SRQ4/H4)
**Critical efficiency data:**
- Traditional manual sight reduction: **"hours each day on sight reduction calculations"**
- 1972 computer goal: **"time required to complete a sight reduction is reduced to seconds"**
- Per-sight estimate: 4-5 sights × 3-5 stars = 12-25 sights daily requiring "hours"

**Quote:** "The conscientious navigator, who takes 4 to 5 sights on each of 3 to 5 stars to obtain an accurate LOP intersection, or a small cocked hat, may spend hours each day on sight reduction calculations."

#### 2. Core Sight Reduction Formulas
**Azimuth formula:**
$$\tan Z = \frac{-\cos\delta \sin \text{LHA}}{\cos L \sin\delta - \sin L \cos\delta \cos \text{LHA}}$$

**Altitude formula:**
$$\sin h = \sin L \sin\delta + \cos L \cos\delta \cos \text{LHA}$$

**Hour angle relationships:**
$$\text{LHA} = \text{GHA} - \text{Long}$$
$$\text{GHA}_{\text{star}} = \text{GHA}_\gamma + \text{SHA}_{\text{star}}$$

Where:
- $Z$ = true azimuth
- $h$ = computed altitude
- $\delta$ = declination
- $L$ = latitude (+ North, - South)
- LHA = Local Hour Angle
- GHA = Greenwich Hour Angle
- SHA = Sidereal Hour Angle
- $\gamma$ = Aries (First Point of Aries)

#### 3. Precession Correction Formulas
$$\Delta\alpha = (3^s.073 + 1^s.336 \sin\alpha \tan\delta) \times \text{fraction of year}$$
$$\Delta\delta = (20''.04 \cos\alpha) \times \text{fraction of year}$$

Where $\alpha$ = right ascension = 360° - SHA

#### 4. Aberration Correction Formulas
$$\Delta\alpha = \frac{-k \sin\lambda \sin\alpha - k \cos\lambda \cos\epsilon \cos\alpha}{\cos\delta}$$
$$\Delta\delta = -k \sin\lambda \cos\alpha \sin\delta + k \cos\lambda \cos\epsilon \sin\alpha \sin\delta - k \cos\lambda \sin\epsilon \cos\delta$$

Where:
- $k = 20''.496$ (constant of aberration)
- $\lambda$ = Sun's true longitude
- $\epsilon$ = obliquity of ecliptic

#### 5. Refraction Correction Formula
$$R = \frac{460 + T}{(460 + 50)(29.83)} \times 132.272727B \times \left[z - \sin^{-1}(0.998115 \sin z)\right]$$

Where:
- $R$ = refraction in arc-seconds
- $T$ = temperature in °F
- $B$ = barometric pressure in inches Hg
- $z$ = observed zenith distance (90° - h)
- Nominal values: T = 50°F, B = 29.83 in

#### 6. Computer Program Logic Flow
```
(1) Input: DR LAT (±XX XX.X) - degrees, minutes, tenths
(2) Input: DR LONG (±XXX XX.X) - West+, East-
(3) Input: DIP (X.X) - from Nautical Almanac, always negative
(4) Input: TIME (XXX XX XX XX) - day, hours, minutes, seconds GMT
(5) Input: STAR (XX) - star number 1-57
(6) Input: ALTITUDE (XX.XXX) - uncorrected measured altitude
(7) Output: TRUE AZIMUTH (XXX.X) - degrees from North through East
(8) Output: ALTITUDE INTERCEPT (±XXX.X) - minutes of arc (Ho - Hc)
→ Return to (4) for next observation
```

#### 7. Accuracy Specifications
- Target accuracy: **0.1' of arc** (one-tenth of a minute)
- Nutation: included in mean position to required accuracy
- 57 navigation stars with SHA and declination stored in memory

#### 8. System Concept
**Components proposed:**
- Digital readout day/night marine sextant
- Portable specialized computer
- Target cost: ~$1,000 (1972) for volume production
- Size: "attache case" for brass-board model

**Quote:** "A low-cost, specialized, small portable computer that operates with the sextant can be developed to relieve the navigator of his calculation tedium and blunders and improve his navigation accuracy."

### Relevance to Python Implementation

**For SRQ4/H4 (Computational Efficiency):**
Historical progression now documented:
| Era | Method | Time per Sight |
|-----|--------|----------------|
| Pre-1972 | Manual tables | ~10-30 minutes |
| 1972 | Proposed computer | "seconds" |
| 1976 | Galaxy 1/Intercepter | ~2 minutes |
| 2022 | Python/Skyfield | 20-50 μs |

**For Algorithm Implementation:**
The formulas provided are directly implementable in Python:
```python
# Azimuth (conceptual)
tan_Z = -cos(dec) * sin(LHA) / (cos(lat) * sin(dec) - sin(lat) * cos(dec) * cos(LHA))

# Altitude (conceptual)
sin_h = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(LHA)
```

**For Corrections:**
- Precession corrections inform long-term ephemeris validation
- Aberration constant $k = 20''.496$ is standard
- Refraction formula provides alternative to Skyfield's built-in

### Key Quotes

> "The reduction of each sextant sight to true azimuth and altitude intercept requires the use of standard forms, The Nautical Almanac, and several volumes of sight reduction tables. The procedure is tedious and time consuming, and errors in data selection and arithmetic are common." (p. 317)

> "To the accuracy of interest, a tenth of a minute of arc, the effect of nutation may be included in the mean position of the star." (p. 320)

> "One of the major advantages of the portable sextant computer system in marine and air navigation is that the time required to complete a sight reduction is reduced to seconds." (p. 320)

---

## Article 27: Vulfovich & Fogilev (2010) - New Ideas for Celestial Navigation

### Citation
Vulfovich, B., & Fogilev, V. (2010). New ideas for celestial navigation in the third millennium. *The Journal of Navigation*, 63(2), 373-378. https://doi.org/10.1017/S0373463309990348

### Relevance to Research
**MODERATE-HIGH RELEVANCE** - Provides alternative algorithmic approaches for position fixing, including iteration method and weighted observation processing. Relevant to multi-body fixes (SRQ3).

### Key Content

#### 1. Weighted Observation Estimation (SRQ3 - Multi-body Fixes)
Novel approach assigning unequal weights to observations based on clustering:

**Weight calculation:**
$$p_i = \frac{K}{d_i}$$

Where $d_i$ is the sum of squared deviations from other observations:
$$d_i = \sum_{j=1}^{n-1}(x_j - x_i)^2$$

**Normalization coefficient:**
$$K = \frac{1}{\sum_{i=1}^{n}\frac{1}{d_i}}$$

**Weighted mean:**
$$\tilde{M}(x) = \sum_{i=1}^{n} p_i x_i$$

**Example from paper:**
- Traditional mean: $\bar{X} = \frac{1}{3}(2+5+9) = 5.33$
- Weighted mean: $\tilde{M} = 0.24 \times 2 + 0.55 \times 5 + 0.21 \times 9 = 5.12$

#### 2. Standard Deviation Estimation
Combined a posteriori and a priori estimate:
$$\tilde{\sigma}(x) = \alpha(n) \cdot \sigma_{apost}(x) + (1-\alpha(n)) \cdot \sigma_{apr}(x)$$

**A posteriori calculation:**
$$\sigma_{apost}(x) = \sqrt{\sum_{i=1}^{n} p_i (x_i - \tilde{M}(x))^2}$$

**Weight coefficient:**
$$\alpha(n) = 0.012 \times n^{1.893}$$

Giving $\alpha(3) = 0.1$ for 3 observations and $\alpha(9) = 0.8$ for 9 observations.

#### 3. Iteration Method for Position Fixing
Alternative to intercept method using direct iteration:

**Core altitude equation:**
$$\sin h = \sin\varphi \sin\delta + \cos\varphi \cos\delta \cos(t_{Gr} + \lambda)$$

**Longitude iteration (from first body):**
$$\cos(t_{Gr1} + \lambda_1) = \frac{\sin h_1 - \sin\varphi_c \sin\delta_1}{\cos\varphi_c \cos\delta_1}$$
$$\lambda_1 = (t_{Gr1} + \lambda_1) - t_{Gr1}$$

**Latitude iteration (from second body):**
$$\tan Q = \cot\delta_2 \cos(t_{Gr2} + \lambda_1)$$
$$\sin(\varphi_1 + Q) = \cos\delta_2 \sin h_2 \cos Q$$
$$\varphi_1 = (\varphi_1 + Q) - Q$$

**Convergence criteria:**
$$|\varphi_{i+1} - \varphi_i| < \epsilon$$
$$|\lambda_{i+1} - \lambda_i| < \epsilon$$

#### 4. Accuracy Context
**Quote:** "Because observations with hand-held sextants have typical uncertainties of about one arc minute, celestial fixes are rarely more accurate than several nautical miles."

#### 5. Algorithm Characteristics
- Alternates between latitude and longitude corrections
- Starts from DR position $M_c(\varphi_c, \lambda_c)$
- Converges to observed position $M_o(\varphi_o, \lambda_o)$
- Suitable for computer calculation
- Convergence is proven mathematically

### Relevance to Python Implementation

**For SRQ3 (Multi-body Fixes):**
The weighted averaging approach provides alternative to simple mean for combining observations:
- Outlier-resistant through inverse-distance weighting
- Applicable to altitude/intercept combinations
- Confidence estimation with sample-size-dependent weights

**For Algorithm Design:**
The iteration method offers alternative to intercept/LOP plotting:
```python
# Conceptual iteration loop
while not converged:
    lambda_new = solve_longitude(h1, phi_current, delta1, t_gr1)
    phi_new = solve_latitude(h2, lambda_new, delta2, t_gr2)
    if abs(phi_new - phi_current) < epsilon and abs(lambda_new - lambda_current) < epsilon:
        converged = True
    phi_current, lambda_current = phi_new, lambda_new
```

### Key Quotes

> "The calculations that are required for the reduction of a celestial sight, if performed by hand, are slow and error-prone and discourage the human navigator from taking sights because of the tedious work involved." (p. 377)

> "Any reasonably accurate algorithm, implemented in a user-friendly program, would encourage navigators to broaden their observational habits and obtain more sights." (p. 377)

> "Independent alternatives to GPS are needed and are required by official policy." (p. 377-378)

---

## Pending Content Needs

### For Literature Review:
- [x] GHA, LHA, declination computation algorithms - Feldman et al. 1972
- [ ] Altitude correction procedures (refraction, dip, semi-diameter)

### For Technical Implementation:
- [x] Python libraries for ephemeris (Skyfield) - Barazzetti 2025
- [x] Nautical Almanac data sources - Hohenkerk et al. 2012
- [ ] Validation datasets

### For Results/Discussion:
- [x] Modern accuracy standards for celestial navigation - Kotlarić 1976
- [x] Comparison benchmarks from contemporary studies - Kotlarić 1976, Tsai 2022, Feldman 1972
- [x] Reference software implementation - NavPac (Hohenkerk et al. 2012)

---

## Article 28: Hohenkerk, Kemp & Nibbs (2012) - Astro Navigation Remembered

### Citation
Hohenkerk, C., Kemp, J., & Nibbs, B. (2012). Astro navigation remembered. *The Journal of Navigation*, 65(3), 381-395. https://doi.org/10.1017/S0373463312000033

### Relevance to Research
**MODERATE-HIGH RELEVANCE** - Provides historical context, describes NavPac (Royal Navy's official software), practical accuracy benchmarks, and ephemeris approximation methods. Author Hohenkerk works at HM Nautical Almanac Office.

### Key Content

#### 1. NavPac Software (Reference Implementation)
Official Royal Navy sight reduction software:

**Quote:** "It has been adopted by the Royal Navy as their primary method of sight-reduction to calculate astro navigation positions."

**Core algorithm:** Least-squares solution of Marcq St Hilaire method
- First introduced in Nautical Almanac in 1989
- Uses economised polynomial coefficients for ephemerides
- Valid 1986-2015 in current version
- Handles all corrections automatically (parallax, HP, semi-diameter, refraction)

**Software features:**
- Rise/Set times including azimuth checks
- Altitude and azimuth calculations
- Great circle/rhumb line routes
- Complete sight reduction with all corrections
- 95% confidence ellipse computation

#### 2. Ephemeris Approximation Methods
**Chebyshev polynomials for planetary positions:**
- Used for Venus, Mars, Jupiter, Saturn
- Represents coordinates: RA, Dec, GHA, Semi-Diameter, Horizontal Parallax
- Target accuracy: **0.01 minutes of arc**
- Coefficients converted to standard polynomials for user convenience

**Moon ephemeris:**
- Monthly Chebyshev coefficients for ecliptic coordinates
- Transform algorithm to GHA and Dec
- Ecliptic coordinates reduce number of terms needed

**Star approximations:**
- Least-squares mixed function: linear in time + sine/cosine terms
- Represents proper motion, aberration, precession, nutation
- Different term counts for 1-year vs 5-year validity

#### 3. Practical Accuracy Benchmarks
From field observations:

| Skill Level | Accuracy Achieved |
|-------------|-------------------|
| RN top navigator | **0.4 NM** |
| Trainee | 4 NM |
| Ideal conditions | ~1 NM (estimated) |

**Quote:** "I once witnessed people of varying experience taking Sun sights in the Bay of Biscay on a lovely calm day. Needless to say, one of the Royal Navy's top navigators was within 0·4 nautical miles, while the trainee was 4 nautical miles from the GPS position!"

#### 4. Correction Procedures Documented
Complete list of altitude corrections needed:
- **Index error correction**
- **Height of eye (dip)**
- **Refraction** (temperature and pressure dependent)
- **Parallax in Altitude** (significant for Moon)
- **Augmentation of Moon's semi-diameter**
- **Upper/Lower limb** (Sun and Moon)

**Quote on Moon corrections:** "These include 'Parallax in Altitude' and 'Augmentation of the Moon's semi-diameter'. They were straightforward to apply but... the subject of distrust and suspicion."

#### 5. Marcq St Hilaire vs Longitude by Chronometer
Comparison of two sight reduction methods:

**Marcq St Hilaire (Intercept method):**
- Assumes latitude and longitude
- Calculates expected altitude at assumed position
- Intercept = Observed - Calculated altitude
- Universal method, works for all observations

**Longitude by Chronometer:**
- Assumes latitude only
- Calculates longitude where position line crosses that latitude
- Breaks down when body is near meridian
- Simpler when latitude known (e.g., from Polaris)

**Quote:** "The Marcq St Hilaire method... is universal, in that it can be used for all observations, even those where the chosen body is close to the observer's meridian."

#### 6. Historical Context
Timeline of practice:
- 1950s-1960s: Daily routine of morning/noon sights
- Pre-1989: Manual calculations with NA and trigonometric tables
- 1989: Least-squares algorithm added to Nautical Almanac
- 1990s: NavPac software first published
- Present: Astro navigation as GPS backup

#### 7. Sight Reduction Tables Reference
**AP3270/NP303 (Rapid Sight Reduction Tables):**
- Volume 1: Selected stars (7 best stars for each lat/LHA combination)
- Published every 5 years
- Pre-computed altitude and azimuth
- Author Hohenkerk rewrote explanations for marine navigators

### Relevance to Python Implementation

**For Validation:**
- NavPac provides authoritative reference for comparison
- 0.4 NM benchmark for expert-level accuracy
- 0.01' target for ephemeris accuracy

**For Algorithm Design:**
- Least-squares Marcq St Hilaire method is standard approach
- Chebyshev polynomial approximation for efficiency (vs full ephemeris)
- Automatic correction handling simplifies user interface

**For Methodology:**
- Confidence ellipse computation for uncertainty estimation
- All correction procedures documented

### Key Quotes

> "The goal of the NavPac software is not only to take the tedium out of astro navigation, but also to provide a learning tool through the 'sight-reduction' steps." (p. 392)

> "However, no-one can find their position using astro navigation unless they can take accurate sights." (p. 394)

> "A prudent mariner keeps a sextant under his bunk, an almanac in his bookcase and the knowledge of how to use them in his head." (p. 395)

> "The object of the NA is, and has been since the very first edition in 1767, 'To provide, in a convenient form, the data required for the practice of astronomical navigation'." (p. 389)

---

### Article #32
**Citation:** Li, C., Zheng, Y., Zhang, C., Yuan, Y., Lian, Y., & Zhou, P. (2014). Astronomical Vessel Position Determination Utilizing the Optical Super Wide Angle Lens Camera. *The Journal of Navigation, 67*(4), 633–649. DOI: 10.1017/S0373463314000058

**Relevance:** Automated celestial navigation using fisheye camera; robust estimation algorithm

#### Novel Approach:
- Images celestial bodies AND horizon simultaneously with fisheye camera
- No manual sextant operation required
- 24-hour capability (not limited to twilight)
- ~100 stars imaged per frame

#### Fisheye Equisolid Projection with Distortion Model:
$$\theta = 2\arcsin\left(\frac{r}{2f}\right) + k_1\left[\arcsin\left(\frac{r}{2f}\right)\right]^2 + k_2\left[\arcsin\left(\frac{r}{2f}\right)\right]^3 + k_3\left[\arcsin\left(\frac{r}{2f}\right)\right]^4$$

Where:
- $r$ = distance from projection point to principal point
- $f$ = focal length
- $k_1, k_2, k_3$ = radial distortion coefficients

#### Horizon Line Fitting (Least Squares):
Theoretical cosine of semiangular field:
$$\cos\theta = \frac{R\sin i\cos\phi - \sqrt{h^2 + 2Rh}\cos i}{R + h}$$

Where:
- $R$ = Earth radius
- $h$ = height of imaging center
- $i$ = obliquity of camera principal optical axis

#### Coordinate Transformations:
1. Image → Projection (translation + rotation by $\phi_0$)
2. Projection → Camera (Cartesian → polar → spherical)
3. Camera → Horizontal (rotation by $-i$)

#### Robust Estimation Weight Function:
$$\bar{p}(|\tilde{v}_i|) = \begin{cases} 1 & |\tilde{v}_i| < k_0 \\ \frac{k_0}{|\tilde{v}_i|}\left(\frac{k_1 - |\tilde{v}_i|}{k_1 - k_0}\right)^2 & k_0 < |\tilde{v}_i| < k_1 \\ 0 & |\tilde{v}_i| > k_1 \end{cases}$$

With $k_0 = 1.5$, $k_1 = 3.0$ (outlier rejection threshold)

#### Terrestrial Refraction Correction:
$$\rho_0 = a\tan\left(\frac{\pi}{2} - h_H\right) + b\tan^3\left(\frac{\pi}{2} - h_H\right)$$

Where $a = 60.27''$, $b = -0.0669''$

Dip correction: $\text{dip} = 0.0294\sqrt{H}$ (degrees)

#### Experimental Results (Rizhao, China, Oct 2012):
| Metric | Value |
|--------|-------|
| Stars per image | ~100 (avg 107.7) |
| Single image accuracy | 2.98 ± 1.75 NM |
| 10-image average accuracy | 0.40 ± 1.32 NM |
| Exposure time | 0.2 seconds |
| Total processing time | <13 seconds |
| Ground truth | GPS |

#### Advantages Over Sextant:
- Automated observation and processing
- 24-hour capability (not just twilight)
- High efficiency (<13 sec vs 15+ min for sextant)
- System miniaturization (no bubble level needed)
- ~100 stars per observation vs 2-3 with sextant

#### Key Quote:
> "Even a skilled navigator utilising a sextant to determine the AVP by observing 2–3 stars will generally need 15 minutes. However, when utilising the fisheye camera, it only needs 0.2 seconds of exposure time, about 11 seconds of image transmission time, and approximately 1 second of processing time."

---

### Article #31 (Non-Peer-Reviewed)
**Citation:** Garvin, M. J. (2010). *Future of Celestial Navigation and the Ocean-Going Military Navigator.* Master's Project, Old Dominion University. https://digitalcommons.odu.edu/ots_masters_projects/41

**Relevance:** Survey data on actual celestial navigation usage by military navigators; GPS vulnerability context

#### Survey Results (n=78 U.S. Army Marine Deck Warrant Officers):

**Celestial Navigation Usage:**
| Finding | Percentage |
|---------|------------|
| Believe celestial should continue | 57.1% |
| Prefer electronic navigation despite GPS errors | 68% |
| Do NOT use celestial every time at sea | 52.6% |
| Never sailed upon open ocean | 2.6% |
| More than 20 ocean voyages | 66.7% |

**Key Survey Question Results:**
- Q5: "Every time at sea, performed celestial navigation" → Mean 2.53 (Disagree)
- Q6: "Used celestial to check compass" → Mean 2.74 (Neutral)
- Q12: "Prefer electronic over celestial" → Mean 3.78 (Agree)
- Q14: "Only do celestial for school" → Mean 2.78 (Neutral)

**Would use celestial more if:**
- Easier, faster, less weather-dependent (27.4%)
- Taught with calculators/computer programs (21.9%)
- Redundant electronics unavailable (13.7%)
- Only option left (11.0%)

#### Regulatory Context:

**USCG Requirements:**
- 2002: Passing grade for celestial reduced from 90% to 80%
- Still required for "Upon Oceans" license endorsement
- 46 CFR 10.215(c) and 10.401(d)

**STCW (International Convention):**
- Requires deck officers to show celestial proficiency
- Under comprehensive review (2008-2010)
- U.S. position: reduce but not eliminate celestial requirements

**SOLAS:**
- Does NOT require ships to carry sextant (Chapter V, Regulation 19)
- Requires GPS or electronic navigation system

#### GPS Vulnerabilities (citing GAO 2009):

> "It is uncertain whether the Air Force will be able to acquire new satellites in time to maintain current GPS service without interruption."

- Aging satellites not replaced quickly enough
- 80% probability of maintaining 24 satellites (2010-2014)
- 50-80% probability (2018-2020)
- Susceptible to jamming (6 GPS jammers captured in Iraq)
- Air Force Chief Gen. Schwartz: "military needs to wean itself off GPS dependency"

#### Interview with Senior Instructor (A. Lipson, March 2010):

> "The need for continued celestial instruction is a must especially in this ever changing electronic age! As prudent mariners, the need to remain proficient in all means of navigation is only professional."

> "Instead of eliminating celestial instruction, the need to alter celestial instruction to bring together celestial and electronic means of navigation is a viable means of enticing our young officers to stay on course."

#### Historical Evolution Summary:
- Quadrant, backstaff, kamal → sextant (1732)
- Longitude solved by chronometer (late 18th century)
- Bowditch simplified calculations (1802)
- Electronic systems: LORAN-C (terminated 2010), GPS (operational 1996)

#### Recommendations from Study:
1. Integrate electronic navigation with celestial module
2. Allow navigational calculators and computer programs
3. Implement distributed learning for continuing education
4. Adopt MERPAC recommendations for revised exam requirements

---

### Article #30
**Citation:** Dunlap, G. D. (1993). Captain P. V. H. Weems and the Transition from Marine to Air Navigation. *NAVIGATION: Journal of The Institute of Navigation, 40*(1), 1–8. DOI: 10.1002/j.2161-4296.1993.tb02290.x

**Relevance:** Historical context for evolution of sight reduction methods and navigation instruments

#### Historical Evolution of Sight Reduction Methods:

| Era | Method | Characteristics |
|-----|--------|----------------|
| Pre-1920s | Logarithms | "Long and tedious mathematical procedure" |
| 1924-27 | Ogura tables (Japanese) | Basis for Weems line-of-position book |
| 1920s-30s | Short methods | Weems → Ageton → Dreisonstok |
| 1930s+ | Inspection tables | HO 214, HO 249, HO 229 |
| 1930s | Star Altitude Curves | Precomputed 3-star altitude circles |
| Modern | Pocket calculators | Same basic formulae, preprogrammed |

#### Four Equipment Items for Celestial Navigation:
1. **Observing instrument**: Marine sextant (1732), bubble octant for aviation
2. **Timekeeper**: Second-setting watch (Weems patent 1935)
3. **Almanac**: GHA tabulation (1941) replacing right ascension
4. **Sight reduction method**: Tables or computation

#### Almanac Evolution:
- **Pre-1933**: Ephemeral data with right ascension (for astronomers)
- **1933**: First Air Almanac (Naval Observatory, Weems design) - dropped after one year
- **1930s**: British publication of Air Almanac
- **1941**: U.S. Naval Observatory adopted direct GHA vs GMT tabulation

#### Star Altitude Curves (Weems' proudest achievement):
- Precomputed altitude circles of three selected stars
- Different color for each star
- Entering arguments: latitude (vertical), local sidereal time (horizontal)
- Required only yearly GHA Aries sheet
- Adopted by Army Air Corps at start of WWII
- 10 years development, privately funded
- Both patented and copyrighted

#### Key Instruments Developed:
- **Weems Aircraft Plotter (Mk II)**: 50+ years unchanged design, over 1 million manufactured
- **E6B Computer**: Dalton design, Weems marketing rights, copied worldwide
- **Second-setting watch**: Rotating seconds dial, patent 1935
- **Lindbergh watch**: Longines-Weems with hour angle graduation

#### Institute of Navigation:
- Co-founded 1945 by Captain Weems and Dr. Samuel Herrick

#### Key Quote:
> "The various methods are still explained in HO 9 Bowditch—The American Practical Navigator. Today both the almanac data and the sight reduction computation can be handled by a small pocket navigation computer."

---

### Article #29
**Citation:** Pepperday, M. (1994). The Nautical Almanac's Faulty Calculator Instructions. *The Journal of Navigation, 47*, 89. (Forum section, peer-reviewed)

**Relevance:** Critical analysis of Nautical Almanac calculator instructions with alternative formulas for sight reduction

#### Key Contributions:

**1. Tangent-Based Azimuth Formula (Single Expression, No IF-Testing):**
$$Z = \tan^{-1}\left[\frac{\sin \text{LHA}}{\cos \text{LHA} \cdot \sin \text{lat} - \tan \text{Dec} \cdot \cos \text{lat}}\right] + 180°$$
- Uses two-argument arctan (rectangular-polar key)
- Result automatically falls within 0°–360° range
- Avoids IF-testing which is "intellectually sloppy" and error-prone

**2. Complete Altitude Correction Expression:**
$$H_o = H_s + I/60 - 0.03\sqrt{h} - \frac{0.0167}{\tan(H_s + 7.31/(H_s + 4.4))} + \text{passage} + HP \cdot \cos H_s$$

Where:
- $H_s$ = sextant altitude
- $I$ = index correction (minutes)
- $h$ = height of eye (metres)
- Refraction term: $-0.0167/\tan(H_s + 7.31/(H_s + 4.4))$ works 0°–90°
- $HP$ = horizontal parallax (for Moon)

**3. Passage Correction Formula (Running Fix):**
$$\text{passage} = \cos(\text{course} - \text{azimuth}) \times (\text{fix time} - \text{sight time}) \times \text{Speed} / 60$$
- Time in hours, speed in knots
- Result in degrees
- Critical for moving vessel—omitted from RGO instructions

**4. Standard Error of Fix:**
$$S = 60 \times \sqrt{\frac{F - D \cdot d\text{lat} - E \cdot dL}{n - 2}}$$
- $S$ in nautical miles
- Provides quality measure for fix
- "If larger than a mile or two, look for an explanation"

**5. Least Squares Fix Summations:**
- $n$ = count of sights
- $B'$ = sum of $\sin 2Z$
- $C'$ = sum of $\cos 2Z$
- $F$ = sum of $p^2$ (intercepts squared)
- Derived: $B = B'/2$, $C = (n - C')/2$, $A = n - C$

#### Critical Points for Algorithm Design:

1. **Semi-diameter application:** Should apply to computed altitude, not observed
   - Makes lower/upper limb administratively different "bodies"
   - Avoids user prompts about which limb

2. **Single DR position:** All intercepts from single DR lat/long
   - Easier plotting
   - Erroneous readings obvious when multiple sights of same body
   - RGO's multiple DR positions "astonishing departure from standard procedure"

3. **Iteration unnecessary:** "After a century of use at sea... if the intercept method were deficient someone might have noticed"

4. **Formula choice guidance:**
   | Platform | Azimuth | Altitude |
   |----------|---------|----------|
   | Most computers/calculators | tan formula | cos rule |
   | Small calculators | tan formula | sine rule by-product |
   | BASIC without rec-pol | cos rule | five parts formula |

#### Commercial Navigation Computers (1994):
- Celesticomp V, CN2000, Merlin II, NC99, Petrel
- All compute least squares fix
- All apply passage correction
- Standard operation: enter clock time + sextant altitude → displays azimuth + intercept

#### Key Quotes:

> "IF-testing is intellectually sloppy, introduces the risk of overlooking some contingency, eats computer space, makes a program hard to read and is practically always unnecessary."

> "It is the great virtue of the intercept method that there is no need to be fussy about DR position. DR uncertainty is of no consequence. Not even a hundred mile error would matter."

> "Error S, without which the navigator has no measure of the quality of the fix, would also bear some discussing."

> "Even land surveyors, who are at home with iterative procedures and who look for an accuracy of metres rather than miles, do not iterate a celestial fix."

---

### Article #36
**Citation:** Seidelmann, P. K., Janiczek, P. M., & Haupt, R. F. (1976). The Almanacs—Yesterday, Today and Tomorrow. *NAVIGATION: Journal of The Institute of Navigation, 24*(4), 303–312. doi:10.1002/j.2161-4296.1976.tb00755.x

**Relevance:** Authoritative U.S. Naval Observatory paper documenting the evolution from tabular almanacs to computational methods, including Chebyshev polynomials for ephemeris calculations. Essential historical context for modern Python-based implementations.

#### Historical Timeline of Almanacs:
| Year | Development |
|------|-------------|
| 1679 | Connaissance des Temps (first official almanac, France) |
| 1767 | British Nautical Almanac (for lunar distances) |
| 1856 | American Ephemeris and Nautical Almanac |
| 1929 | GHA tabulation in arc introduced (Weems suggestion) |
| 1933 | First American Air Almanac experiment |
| 1937 | British Air Almanac first published |
| 1941 | American Air Almanac permanent publication |
| 1950 | American Nautical Almanac completely revised |
| 1976 | "Almanac for Computers" introduced |

#### Key Innovation: GHA Tabulation
Traditional method using Right Ascension (R.A.):
- Subtract R.A. from sidereal time
- Convert result from time to arc
- Error-prone and time-consuming

Modern method using GHA:
- GHA = GHA Aries + S.H.A. (simple addition)
- S.H.A. = Sidereal Hour Angle = 360° - R.A.

#### Air Almanac Features (1941):
- Single sheets for single day
- 10-minute interval GHA tabulation (Sun, Aries, planets, Moon)
- Hourly declination for Sun and planets
- 10-minute declination for Moon
- S.H.A. of stars for 4-month period
- Convenient interpolation tables

> "The American Air Almanac... marks the most important step in the simplification of the art of navigation since the work of Marcq-Saint-Hilaire."

#### Almanac for Computers (1976):
Contained polynomial expressions for:
- Sidereal hour angle and declination of Sun, Moon, planets (0'.1 accuracy)
- Greenwich Hour Angle of Aries
- Mean positions of navigational stars
- Nutation and precession expressions
- Refraction, dip, augmentation, parallax corrections
- Sunrise, sunset, twilight calculations

**Chebyshev Polynomials Section:**
- Truncatable to any desired accuracy
- Apparent right ascension and declination
- True distance, semi-diameter, ephemeris transit
- Available for Sun, Moon, and all planets
- Accuracy matching American Ephemeris

#### Computer Classification for Navigation (1976):
1. **Small general-purpose calculators:** Add/subtract/multiply/divide; assist with math
2. **Medium-range calculators:** Trigonometric functions; solve celestial triangle
3. **Programmable calculators:** Store star positions, precession/nutation; compute observable stars, azimuths, intercepts; determine position
4. **Special-purpose celestial navigation computers:** Less flexible but easier operation; permanent program storage
5. **Multipurpose navigational computers:** Combine data from inertial, Loran, Omega, satellite, celestial systems

#### Critical Insight on Calculation Speed:
> "The medium range device can be used to calculate the altitude and azimuth of a star by solution of the celestial triangle, but an experienced navigator can use the sight reduction tables, such as HO 249, faster than he can use the calculator, and with less likelihood of error."

*Note: This was true in 1976 but reversed with modern computers.*

#### Vision for Future Navigation:
- Image intensification sextants for 24-hour operation
- Combined electronic and celestial methods with computer
- Celestial as backup when electronic systems fail
- "Celestial navigation... should remain a part of any navigational system, and possibly a life saver."

#### Three Ephemeris Processes:
1. **Theory construction:** Mathematical definition, solution of equations of motion
2. **Table construction:** Reduce theory to arithmetic operations
3. **Application:** Use tables, convert and arrange numerical results

#### Key Quotes:

> "By judicious use of computers, or calculators, the navigator can make more observations and use more data to provide more accurate positions, while spending less of his time on tedious calculations."

> "The printed almanacs are indispensable for cases when equipment has failed or alternate methods are unavailable."

> "With the information published in this 'Almanac for Computers'... the navigator will be able to select the celestial data in the form best suited for his computer, or calculator."

---

### Article #35
**Citation:** McConnell, P. (2008). Slide Rules, Sextants, and the Sky: Instrumentation Milestones In Aerial Celestial Navigation. *Journal of the Oughtred Society, 17*(2), 40–50.

**Relevance:** Comprehensive historical review of celestial navigation instrumentation evolution from marine to aerial applications. Documents key formulas, computational methods, and the transition from manual to automated systems.

#### Celestial Navigation Prerequisites:
1. Means of predicting celestial body positions (almanac)
2. Accurate timekeeper (chronometer)
3. Method of computing relative positions (spherical trigonometry)
4. Observation instrument (sextant)

#### Cosine Formula (Fundamental):
For spherical triangle with sides $a$, $b$, $c$ and angle $A$:
$$\cos a = \cos b \cos c + \sin b \sin c \cos A$$

#### Cosine-Haversine Formula:
$$\text{hav } z = \text{hav}(L \mp d) + \cos L \cos d \cdot \text{hav } t$$

Where:
- $z$ = zenith distance to body
- $L$ = latitude of assumed position
- $d$ = declination of body
- $t$ = meridian angle
- $(L \mp d)$: subtract if same name, add if opposite

**Computed Altitude:** $H_c = 90° - z$

#### Bygrave Position Line Slide Rule Formulas:
Tangent-Cosine relationships (divides triangle into two right triangles):

$$\tan y = \frac{\tan d}{\cos t}$$

If $L$ and $d$ same name: $Y = (90° - L) + y$
If $L$ and $d$ opposite name: $Y = (90° - L) - y$

$$\tan A = \frac{\cos y \cdot \tan t}{\cos Y}$$

$$\tan h = \cos A \cdot \tan Y$$

Where $h$ = altitude, $A$ = azimuth, $y$ and $Y$ = intermediate values.

**Advantages:** Uses only tangent and cosine; helical scales provide 1 arcmin accuracy; solution in 2-3 minutes.

#### Historical Timeline:
| Year | Development |
|------|-------------|
| 1730 | Marine sextant (Hadley/Godfrey) |
| 1735 | Marine chronometer (Harrison) |
| 1837 | Celestial Line of Position (Sumner) |
| 1874 | Intercept Method (Marcq Saint Hilaire) |
| 1897 | Bubble sextant (Andrée) |
| 1901 | First aerial bubble sextant (Marcuse) |
| 1913 | Schwarzschild bubble sextant |
| ~1920 | Bygrave Position Line Slide Rule |
| 1931 | H.O. 208 & H.O. 211 tables |
| 1936 | H.O. 214 Tables of Computed Altitude |
| 1941 | Air Almanac resumed publication |
| 1943 | H.O. 218 for selected stars |
| 1951 | H.O. 249 Sight Reduction Tables |

#### Saint Hilaire Intercept Method Summary:
1. Select celestial body, assumed position (AP), observation time
2. Compute altitude ($H_c$) and azimuth from navigational triangle
3. Measure observed altitude ($H_o$) with sextant
4. Calculate intercept: $a = H_o - H_c$
5. LOP perpendicular to azimuth, distance $a$ from AP
6. Toward body if $H_o > H_c$; away if $H_o < H_c$

#### Worked Example Comparison (Sun sight, 41°17'N, 72°26'W):
| Method | Computed Altitude | Azimuth |
|--------|-------------------|----------|
| Cosine-Haversine | 23°25' | (not computed) |
| Bygrave Slide Rule | 22°24' | 218°02' |
| H.O. 214 Tables | 22°30' | 217°42' |
| USNO Web Service | 22.5° | 218.0° |

#### Bubble Sextant Averaging Methods:
1. **Mechanical averagers:** Sum fractions of readings (e.g., Mark IXA: 60 readings over 2 minutes)
2. **Median recorders:** Graphical record, navigator estimates median (A-7, A-10, A-12 sextants)
3. **Mechanical integrators:** Continuous running average

#### Key Quotes:

> "The principal problem of celestial navigation is that of converting the coordinates of one system to those of another."

> "By carefully selecting the optimum formula and scale layout, Bygrave reduced the solution of the Navigational Triangle to a dozen moves. An experienced navigator could use the instrument to obtain a solution in two or three minutes."

> P.V.H. Weems on Bygrave: "...is probably the most convenient mechanical computer for obtaining position lines from sextant observations."

---

### Article #34
**Citation:** Kaplan, G. H. (1996). A Navigation Solution Involving Changes to Course and Speed. *NAVIGATION: Journal of The Institute of Navigation, 43*(4), 470–482. doi:10.1002/j.2161-4296.1996.tb01933.x

**Relevance:** Extends Kaplan 1995 single-leg algorithm to handle multiple voyage legs with course/speed changes. Essential for real-world implementation where vessels maneuver during observation periods.

#### Problem Statement:
Traditional methods assume constant positional error during running fix, which is "unlikely to be true, even approximately, except over very short periods." Course/speed changes further complicate the assumption.

#### Single-Leg Conditional Equation (from Kaplan 1995):
$$a = \frac{\partial H_c}{\partial \phi} \left( \Delta\phi_0 + \frac{\partial f}{\partial \phi_0} \Delta\phi_0 + \frac{\partial f}{\partial \lambda_0} \Delta\lambda_0 + \frac{\partial f}{\partial C} \Delta C + \frac{\partial f}{\partial S} \Delta S \right) + \frac{\partial H_c}{\partial \lambda} \left( \Delta\lambda_0 + \frac{\partial g}{\partial \phi_0} \Delta\phi_0 + \frac{\partial g}{\partial \lambda_0} \Delta\lambda_0 + \frac{\partial g}{\partial C} \Delta C + \frac{\partial g}{\partial S} \Delta S \right)$$

Where:
- $a$ = altitude intercept ($H_o - H_c$)
- $\phi, \lambda$ = latitude, longitude at observation time
- $\phi_0, \lambda_0$ = latitude, longitude at fix time
- $C, S$ = course, speed
- $f, g$ = rhumb-line sailing formulas

#### Multiple-Leg Extension:
For $N$ legs, solve for $2N + 2$ free parameters:
- $\Delta\phi_1, \Delta\lambda_1$ (fix coordinates on first leg)
- $\Delta C_j, \Delta S_j$ for $j = 1$ to $N$ (course/speed each leg)

**Continuity constraint:** Start of leg $j$ must equal end of leg $j-1$.

#### Multiple-Leg Conditional Equation:
$$a_{ik} = \frac{\partial H_c}{\partial \phi} \left( \frac{\partial f_k}{\partial \phi_1} \Delta\phi_1 + \frac{\partial f_k}{\partial \lambda_1} \Delta\lambda_1 + \sum_{j=1}^{N} \left( \frac{\partial f_k}{\partial C_j} \Delta C_j + \frac{\partial f_k}{\partial S_j} \Delta S_j \right) \right) + \frac{\partial H_c}{\partial \lambda} \left( \frac{\partial g_k}{\partial \phi_1} \Delta\phi_1 + \frac{\partial g_k}{\partial \lambda_1} \Delta\lambda_1 + \sum_{j=1}^{N} \left( \frac{\partial g_k}{\partial C_j} \Delta C_j + \frac{\partial g_k}{\partial S_j} \Delta S_j \right) \right)$$

Where $k$ = leg of current observation, $i$ = observation index.

#### Partial Derivative Cases:
| Case | Condition | Evaluation |
|------|-----------|------------|
| $j > k$ | Future leg | Partials = 0 (future doesn't affect current position) |
| $j = k$ | Current leg | Use single-leg formulas from Kaplan 1995 |
| $j < k$ | Past leg | Recurrence relations through continuity constraints |

#### Time-Reversal Strategy:
For practical use, process observations backwards:
- Leg C (final) → Leg 1 (contains fix coordinates)
- Leg B → Leg 2
- Leg A (first) → Leg 3

This provides fix position on the current/final leg.

#### Sample Calculation Results (12 observations, 3 legs):

**Initial errors:** Position ≈ 100 nmi, Course ≈ 2-3°, Speed ≈ 3 kn

**After solution:**
| Parameter | Truth | Solution | Error |
|-----------|-------|----------|-------|
| Latitude | +47.932° | +47.911° | -0.021° |
| Longitude | -37.745° | -37.797° | -0.052° |
| Course (Leg 1) | 70.00° | 70.26° | +0.26° |
| Speed (Leg 1) | 15.00 kn | 14.77 kn | -0.23 kn |

**Fix error:** 2.5 nmi (with 0.7 arcmin observation scatter)

#### Scaling Results:
| Obs. Uncertainty | Fix Lat Error | Fix Lon Error |
|------------------|---------------|---------------|
| ±0.05' | +0.001° | +0.001° |
| ±0.2' | -0.001° | -0.003° |
| ±0.7' | -0.021° | -0.052° |
| ±1.5' | -0.047° | -0.113° |

"Errors in the parameters solved for are directly proportional to the scatter in the observations."

#### Implementation Notes:

1. **Minimum observations:** At least 6 per leg recommended
2. **Correlation check:** Parameter correlation matrix reveals solution reliability
3. **Empty legs:** Valid solutions possible with some empty legs
4. **Iteration:** Usually 2 iterations sufficient for convergence
5. **STELLA software:** Implemented for U.S. Navy operational use

#### Key Quotes:

> "The error in position is assumed to change, and its rate of change is assumed to be different along each leg."

> "The algorithm is conceptually simple, uses ordinary least-squares procedures, makes efficient use of the available observations, and is relatively easy to implement in computer code."

> "For altitude sights made with a hand-held sextant (typical accuracy ±1 arcmin), at least eight observations, spread over several hours, are generally needed for reliable course and speed estimates."

> "An automated shipboard star tracker might make the necessary observations and obtain a good solution in a matter of minutes."

---

### Article #33
**Citation:** Morrison, G. D. (1981). Most Probable Fix Position Reduction. *NAVIGATION: Journal of The Institute of Navigation, 28*(1), 1–8. doi:10.1002/j.2161-4296.1981.tb01435.x

**Relevance:** Foundational paper for implementing least-squares position fixing with minimum variance technique, fix bias estimation, and error analysis for outlier detection.

#### LOP Normal Form Equation:
$$X_1 \sin(Zn) + X_2 \cos(Zn) = a$$

Where:
- $X_1$ = east-west coordinate (positive east)
- $X_2$ = north-south coordinate (positive north)
- $Zn$ = azimuth (bearing from true north)
- $a$ = altitude intercept (distance from assumed position)

#### Error Model:
Total error in each observation:
$$\xi_i = r + e_i$$

Where:
- $r$ = common repeatable error (sextant miscalibration, observer bias)
- $e_i$ = stochastic variable (random error) for each observation

**Assumptions:**
1. $E(\xi_i) = r$ (expected value biased by constant)
2. $E(e_i) = 0$ (mean of unbiased errors is zero)
3. $E(e_i e_j) = 0$ for all $i \neq j$ (errors uncorrelated)
4. All observations given equal weighting

#### Measurement Error Equation:
$$X_1 \sin(Zn_i) + X_2 \cos(Zn_i) + \xi_i = a_i$$

#### Individual LOP Error:
$$e_i = a_i - r - X_1 \sin(Zn_i) - X_2 \cos(Zn_i)$$

This is the perpendicular distance from adjusted LOP to Most Probable Position.

#### Minimum Variance System (n > 3 LOPs):
Minimizes $\sum e_i^2$ by solving:

$$\begin{bmatrix} n & \sum x & \sum y \\ \sum x & \sum x^2 & \sum xy \\ \sum y & \sum xy & \sum y^2 \end{bmatrix} \begin{bmatrix} r \\ X_1 \\ X_2 \end{bmatrix} = \begin{bmatrix} \sum a \\ \sum ax \\ \sum ay \end{bmatrix}$$

Where:
- $n$ = number of LOPs
- $x = \cos(Zn_i)$
- $y = \sin(Zn_i)$
- $a = a_i$ (altitude intercept)

#### 3-LOP Solution:
For exactly 3 LOPs, uses traditional inside/outside fix approach:

$$\begin{bmatrix} \sin(Zn_1) & \cos(Zn_1) & 1 \\ \sin(Zn_2) & \cos(Zn_2) & 1 \\ \sin(Zn_3) & \cos(Zn_3) & 1 \end{bmatrix} \begin{bmatrix} X_1 \\ X_2 \\ r \end{bmatrix} = \begin{bmatrix} a_1 \\ a_2 \\ a_3 \end{bmatrix}$$

**Fix Type Determination:**
| Condition | Fix Type |
|-----------|----------|
| Azimuth spread > 180° | Inside fix |
| Azimuth spread < 180° | Outside fix |
| Determinant < 0 | Outside fix |
| Determinant > 0 | Inside fix |

#### Standard Deviation (Fix Quality):
$$\sigma = \sqrt{\frac{\sum e_i^2}{n-3}}$$

#### Error Analysis Interpretation:
| Fix Bias | Std Dev | Conclusion |
|----------|---------|------------|
| Low | Low | Good sextant calibration, accurate observer |
| High | Low | Poor sextant calibration or consistent observer bias |
| Low | High | Inconsistent; random errors present |
| High | High | No conclusions possible |

#### Outlier Detection:
- If $|e_i| \geq 2\sigma$, consider deletion from MPP computation
- Recompute fix after excluding erroneous LOP

#### Worked Example Results:
**8 LOPs with outlier:**
| Metric | Before | After LOP #3 removed |
|--------|--------|----------------------|
| Standard Deviation | 0.6 nmi | 0.2 nmi |
| Fix Bias | -1.1 nmi | -1.3 nmi |
| Fix Distance | 8.6 nmi | 8.6 nmi |
| Fix Bearing | 047.8° | 045.4° |

#### Projection Error Analysis:
LOP assumption valid when:
- Altitudes < 75°
- Assumed position within 15 nmi of actual
- Combined azimuth/projection error < 0.15 nmi RMS

**Iteration capability:** "Assumed position errors in excess of 500 nmi can be reduced to fix position accuracy less than 0.1 nmi with only a single iteration."

#### Solution Algorithm (Appendix II):
```
Compute determinant:
Δ = n·Σx²·Σy² + 2·Σx·Σy·Σxy - Σy²·(Σx)² - Σx²·(Σy)² - n·(Σxy)²

Compute intermediate values:
A₁₁ = Σy²·Σx² - (Σxy)²
A₂₂ = n·Σy² - (Σy)²
A₁₂ = Σxy·Σy - Σx·Σy²
A₁₃ = Σx·Σxy - Σy·Σx²
A₂₃ = Σx·Σy - n·Σxy

Solve:
r = (A₁₁·Σa + A₁₂·Σax + A₁₃·Σay) / Δ
X₁ = [A₂₂·(Σax - r·Σx) - A₁₂·(Σay - r·Σy)] / (A₁₁·A₂₂ - A₁₂²)
X₂ = (Σay - r·Σy - X₁·A₁₂) / A₁₁
```

#### Key Quotes:

> "A quick, accurate method for calculating lines of position from a single assumed position has already been developed. The procedure to reduce these lines of position to a most probable fix position will be discussed within this paper."

> "With more than three observations it became apparent that these techniques could provide a ship's navigator with error analyses of statistical parameters concerning errors of the sextant/observer, identification of erroneous LOP's and the reliability or credibility of a fix position."

> "A navigator always attempts to make each observation as accurately as possible, hence would not weight any given observation any higher than another."

> "Fix bias coupled with the standard deviation will define the distribution."

> "Whenever an LOP's error differs from the standard deviation by a factor of two or greater, consideration should be given to deletion from the computation of fix position."

---

### Article #37
**Citation:** Dunlap, G. D. (1971). Major Developments in Marine Navigation during the Last 25 Years. *NAVIGATION: Journal of The Institute of Navigation, 18*(1), 63–76. doi:10.1002/j.2161-4296.1971.tb00075.x

**Relevance:** Authoritative historical survey from NAVIGATION's 25th Anniversary Issue. Covers post-WWII navigation evolution including celestial navigation, electronic systems, inertial navigation, and submarine navigation developments. Essential context for understanding the transition from traditional to computational methods.

#### Celestial Navigation Status (1971):
> "Celestial navigation is still the most widely used method at sea today. In many ships operating in areas of electronic navigational coverage it has become the secondary method, but for redundancy there is a continued need for the art in these vessels."

#### Sight Reduction Tables Evolution:
- **H.O. Pub. No. 214:** "World War II vintage" - most widely used at time of writing
- **H.O. 229 Tables:** "just appearing on the scene" - noted as "merely a rearrangement of the basic mathematical data"

#### GHA Tabulation Innovation:
> "The Nautical Almanac has undergone a significant change in that the GHA of the body is tabulated directly rather than being stated as right ascension as was previously the case. The tabulation of GHA in the Almanac was introduced by P. V. H. Weems, Captain, USN (Retired)."

#### Chronometer Advances:
**Quartz Crystal Chronometer:**
- "Considerable improvement over the conventional spring-operated type"
- "Electronically, rather than mechanically operated"
- "Does not require gimbal mounting"
- "Highly resistant to acceleration"
- "Maintains an extremely steady rate"
- "Available at prices comparable to those charged for the older type"

#### Submarine Celestial Navigation:

**Problem:** Periscope moves with vessel; deck plane unstable due to roll/pitch

**Solution 1 - Submarine Celestial Altitude Recorder:**
- Used roll and pitch data from gyrocompass
- Required submarine to bring body directly ahead/astern or on beam
- Automatically applied roll/pitch corrections to observed altitude

**Solution 2 - Periscopic Sextant:**
- Self-contained gyroscopic reference platform
- Mounted on lower extremity of general purpose periscope
- Turns with periscope for observations on any bearing
- Records altitude as moving star image coincides with reticle
- Ancillary equipment continuously computes true elevation angle
- Produces printout of altitude and time

**Solution 3 - Type II-A and II-B Systems (Polaris submarines):**
- Integrated with NAVDAC computer
- Computer stores star ephemerides in memory bank
- Continuously computes altitude, azimuth, altitude rate, azimuth rate
- Keeps star within field of view for tracking
- Reference platform from ship's inertial system (SINS)
- Celestial data used to monitor SINS errors and compute optimum resets

#### Radiometric Sextant:
- Radio telescope concept for celestial bodies emitting radio waves
- Receives signals from sun and moon
- Tracks radio source and reads out altitude and azimuth
- Combined with inertial platform for "self-contained passive all-weather navigation system"

#### Star Tracker Units:
- Photoelectric telescope or TV vidicon tube as sensitive element
- Gyroscopic unit for vertical reference
- "Sensitivity is so great in some units that the brighter stars can be tracked even in daylight"

#### Light Amplification Telescope:
- Completely passive (unlike infrared)
- Amplifies ambient light tremendously
- Small enough to fit on standard marine sextant
- "Permits very accurate star observations to be obtained"
- "The horizon, as seen through the scope, is sharply defined"

#### Electronic Navigation Systems Timeline:
| System | Characteristics | Range |
|--------|-----------------|-------|
| Loran A | Pulse system, 1750-1950 kHz | 600-900 mi (day), 1400 mi (night) |
| Loran C | Multipulse, time + phase measurement | 1200-1300 mi (day), 3000 mi (night) |
| Omega | VLF 10.2 kHz, phase comparison | Global (8 stations) |
| Decca | CW 70-135 kHz, phase comparison | 400 mi (day), 250 mi (night) |
| NAVSAT | Doppler shift from satellites | Global |

#### Satellite Navigation Key Insight:
> "The comput­er compares calculated range differences from the known position of the satellite to the estimated position of the ship, with the range rates measured by the doppler shift, and the navigation fixes are obtained by continuous iteration to find the position where the calculated doppler curve of differences agrees best with the measured differences."

#### Inertial Navigation (SINS):
- First ship installation: USS COMPASS ISLAND, fall 1956
- Practical use: N6A systems in NAUTILUS and SKATE for trans-polar voyages (1958)
- "The inertial navigator is simply a highly refined dead reckoning instrumentation"
- Cumulative error with time due to uncompensated gyro drift
- Requires periodic position updates from other fixing methods

#### Speed Measurement Limitations:
**Log Types (all measure speed through water, not over ground):**
1. Force/resistance log (strut movement)
2. Pitometer log (differential pressure)
3. Impeller log (propeller rotation)
4. Electromagnetic (EM) log (Faraday's electromagnetic induction)
5. Doppler sonar (measures over-ground in shallow water only)

**Critical Quote from Griswold (1967):**
> "For purposes of position determination, a 0.5 degree compass error for 100 miles is equivalent to a 0.2 knot log system error during five hours."

#### Merchant vs Naval Navigator:
> "The merchant navigator by and large has maintained his skill in celestial navigation, while his naval counterpart has become rusty, due to overreliance on electronics."

#### Precision Requirements Note:
**Repeatability:** Reliability to return to a point using electronic LOPs
**Predictability:** Reliability to define position in geographic coordinates (often more important)

---

## Article #38: Park et al. 2021 - JPL DE440/DE441 Ephemerides

### Bibliographic Information
- **Authors:** Ryan S. Park, William M. Folkner, James G. Williams, Dale H. Boggs
- **Title:** The JPL Planetary and Lunar Ephemerides DE440 and DE441
- **Journal:** The Astronomical Journal
- **Year:** 2021, Volume 161, Number 3, Article 105
- **DOI:** 10.3847/1538-3881/abd414
- **Institution:** Jet Propulsion Laboratory, California Institute of Technology
- **Classification:** PEER-REVIEWED (Open access, 625+ citations)

### Relevance Assessment
| Criterion | Rating | Notes |
|-----------|--------|-------|
| SRQ1 (Sight Reduction Accuracy) | ★★★★★ | **CRITICAL** - Authoritative ephemeris accuracy for Skyfield validation |
| SRQ2 (Celestial Body Performance) | ★★★★★ | Per-body accuracy data for Sun, Moon, planets |
| SRQ3 (Multi-body Fix Error) | ★★★★☆ | Provides ephemeris uncertainty for error propagation |
| SRQ4 (Computational Efficiency) | ★★★☆☆ | Confirms Skyfield's DE ephemeris source |
| H1 (Accuracy Hypothesis) | ★★★★★ | **ESSENTIAL** - Defines accuracy baseline for celestial body positions |

### Abstract Summary
DE440 and DE441 are JPL planetary/lunar ephemerides computed by fitting numerically integrated orbits to ground/space observations. Compared to DE430:
- 7 years of new data added
- Improved dynamical models
- Jupiter orbit substantially improved (Juno data)
- Saturn orbit improved (Cassini data)
- Pluto orbit improved (stellar occultation + Gaia catalog)

**Key Difference:**
- **DE440:** Higher accuracy for current century, covers 1550-2650
- **DE441:** Less accurate now but no divergence in past, covers −13,200 to +17,191

### Critical Data for Celestial Navigation

#### Ephemeris Accuracy by Body (RMS Residuals):

| Celestial Body | Data Type | RMS Accuracy | Notes |
|----------------|-----------|--------------|-------|
| **Moon** | LLR (recent) | **1.3 cm** | Laser ranging to retroreflectors |
| **Moon** | LLR (1970-1980) | ~20 cm | Early data |
| **Mercury** | MESSENGER range | **0.7 m** | Radio range data |
| **Venus** | Venus Express range | **8 m** | Radio range data |
| **Mars** | MGS/ODY/MRO range | **0.7 m** | Radio range data |
| **Mars** | MEX range | **2 m** | Mars Express |
| **Jupiter** | Juno range | **13 m** | New for DE440 |
| **Saturn** | Cassini range | **3 m** | VLBA + range data |
| **Pluto** | Stellar occultation | **8-11 mas** | R.A./decl. |

#### VLBI Angular Accuracy:

| Measurement | RMS Accuracy | Notes |
|-------------|--------------|-------|
| Mars orbiter Goldstone-Madrid | 0.25 mas | R.A. orientation |
| Mars orbiter Goldstone-Canberra | 0.18 mas | Mixed orientation |
| Juno VLBA R.A. | 0.4 mas | Jupiter orientation |
| Juno VLBA decl. | 0.6 mas | Jupiter orientation |
| Cassini VLBA R.A. | 0.35-0.6 mas | Saturn orientation |
| Cassini VLBA decl. | 0.36-0.8 mas | Saturn orientation |

#### Reference Frame:
- **ICRS/ICRF3:** International Celestial Reference System, 3rd realization
- Inner planet orientations aligned to ICRF3 with ~0.2 mas average accuracy
- Based on VLBI measurements of extragalactic radio sources (quasars)

### Ephemeris Time Coverage

| Ephemeris | Time Span | Recommended Use |
|-----------|-----------|-----------------|
| **DE440** | 1550-2650 | Modern data analysis, current navigation |
| **DE441** | −13,200 to +17,191 | Historical data analysis |

**Accuracy Difference:** < 1 m for planets, < 2 m for Moon during 1970-2020 overlap.

### Bodies Integrated

The ephemeris includes:
- Sun
- 8 planetary system barycenters
- Moon
- Pluto system barycenter
- Lunar libration angles
- 343 asteroids (~90% of asteroid belt mass)
- 30 Kuiper Belt Objects (KBOs)
- KBO ring (modeled as 36 point masses at 44 AU)

### Mass Parameters (GM values):

| Body | GM (km³/s²) |
|------|-------------|
| Sun | 132,712,440,041.279 |
| Mercury | 22,031.869 |
| Venus | 324,858.592 |
| Earth | 398,600.436 |
| Moon | 4,902.800 |
| Mars System | 42,828.376 |
| Jupiter System | 126,712,764.100 |
| Saturn System | 37,940,584.842 |
| Ceres | 62.629 |
| Vesta | 17.288 |

### Coordinate Systems

#### Solar System Barycenter (SSBC):
```
r_ssbc = (Σ mᵢ* rᵢ) / (Σ mᵢ*)

where:
mᵢ* = GMᵢ/c² [1 + (1/2)vᵢ²/c² - (1/2)Σⱼ≠ᵢ GMⱼ/(rᵢⱼc²)]
```

**Important Note:** Earth orbits the Sun, not the SSBC. Solar inertial motion does not significantly affect Earth-Sun distance.

#### Time Systems:
- **TDB:** Barycentric Dynamical Time (integration time)
- **UTC:** Coordinated Universal Time (observation time tags)
- **TAI:** International Atomic Time = UTC + leap seconds
- **TT:** Terrestrial Time = TAI + 32.184 s
- Complex TT→TDB conversion provided (Equation 3)

### Lunar Orientation Model

#### Retroreflector Positions (DE440 Principal Axes frame):

| Retroreflector | X (m) | Y (m) | Z (m) |
|----------------|-------|-------|-------|
| Apollo 11 | 1,591,967 | 690,699 | 21,004 |
| Apollo 14 | 1,652,689 | -520,998 | -109,730 |
| Apollo 15 | 1,554,678 | 98,094 | 765,006 |
| Lunokhod 2 | 1,339,364 | 801,871 | 756,359 |
| Lunokhod 1 | 1,114,291 | -781,299 | 1,076,059 |

#### Libration Angles:
- **φₘ:** Angle from inertial X-axis to mantle equator intersection
- **θₘ:** Inclination of mantle equator from inertial XY plane
- **ψₘ:** Longitude from intersection to prime meridian

### Dynamical Model Updates (vs DE430)

1. **Kuiper Belt Objects:** 30 individual KBOs + ring at 44 AU
2. **Lense-Thirring Effect:** Gravito-magnetic GTR effect from Sun's angular momentum
3. **Vondrak Precession:** More accurate for ±1000 year integrations
4. **Geodetic Precession:** Effect on lunar librations added
5. **Solar Radiation Pressure:** Force on Earth-Moon system

### Earth Tidal Effects on Moon

- Second-degree gravitational Love numbers k₂ⱼ (j = 0,1,2)
- Time-delay tidal model for dissipation
- Different delays for long-period/diurnal/semi-diurnal components
- Causes Moon to recede from Earth
- Transfers angular momentum from Earth rotation to lunar orbit

### Implications for Skyfield/Celestial Navigation

#### Accuracy Guarantee for Navigation:
```
For celestial bodies used in navigation:
- Sun: ~0.7-2 m positional accuracy (via Mars/Mercury range)
- Moon: 1.3 cm to surface, ~2 m orbital
- Planets: 0.7-13 m depending on body

At typical sextant distances (~1 AU):
- 1 m positional error ≈ 0.001 arcsec angular error
- Negligible compared to sextant precision (±0.1-0.3')
```

#### Skyfield Implementation:
- Skyfield uses DE ephemeris files directly
- Default: DE421 or DE440 depending on version
- Can specify DE440/DE441 explicitly
- Provides sub-meter planetary positions for navigation

#### Critical Quote:
> "The high-precision orbits and lunar rotations around the three axes have a wide range of practical and fundamental applications" including "U.S. Nautical Almanac Office & Her Majesty's Nautical Almanac Office."

### Observation Data Summary

#### Moon (LLR):
- 50 years of data (1970-2020)
- 5 retroreflectors: Apollo 11/14/15, Lunokhod 1/2
- 5 stations: McDonald, OCA, Apache Point, Matera, Wettzell

#### Inner Planets:
- Mercury: MESSENGER 2011-2016
- Venus: Venus Express 2006-2014
- Mars: Viking, Pathfinder, MGS, ODY, MRO, MEX (1976-2020)

#### Outer Planets:
- Jupiter: Juno 2016-2020 (major improvement)
- Saturn: Cassini 2004-2018 (major improvement)
- Uranus/Neptune: Voyager flybys + astrometry
- Pluto: Stellar occultations + Gaia catalog

### Key Equations for Implementation

#### Point-Mass Acceleration (PPN formulation - Eq. 27):
Full post-Newtonian gravitational interaction with β, γ parameters constrained to unity per GTR.

#### Extended Body Interaction (Eq. 28):
Spherical harmonic expansion for non-spherical gravity:
- Earth: Zonal harmonics to degree 5
- Moon: Degree/order 6 gravity field
- Sun: J₂ (solar oblateness)

#### Lunar Moment of Inertia (Eq. 54):
Time-varying due to tidal/spin distortion with k₂,M Love number.

### Validation for Celestial Navigation Algorithm

**For SRQ1/H1:**
- Ephemeris accuracy: Sub-meter for navigational bodies
- Angular accuracy: Sub-milliarcsecond for planets
- Moon accuracy: 1.3 cm (surface), negligible for navigation

**For SRQ2/H2:**
- Per-body accuracy documented
- Sun/Moon data most relevant for navigation
- All navigational stars: Fixed positions from star catalogs

**Error Budget Contribution:**
```
Ephemeris error contribution to celestial fix:
- At 1 AU: 1 m ≈ 0.0007 arcsec
- At Moon distance: 2 m ≈ 0.001 arcsec
- Sextant precision: 0.1-0.3 arcmin = 6-18 arcsec

Conclusion: Ephemeris error is NEGLIGIBLE (~0.01% of total)
- Dominant errors: Sextant measurement, horizon dip, refraction
```

---

## Article #39: Standish & Fienga 2002 - Ephemeris Accuracy Limits from Asteroid Masses

### Bibliographic Information
- **Authors:** E. Myles Standish, Agnès Fienga
- **Title:** Accuracy limit of modern ephemerides imposed by the uncertainties in asteroid masses
- **Journal:** Astronomy & Astrophysics
- **Year:** 2002, Volume 384, Pages 322-328
- **DOI:** 10.1051/0004-6361:20011821
- **Institutions:** JPL, California Institute of Technology; IMCCE, Paris
- **Classification:** PEER-REVIEWED

### Relevance Assessment
| Criterion | Rating | Notes |
|-----------|--------|-------|
| SRQ1 (Sight Reduction Accuracy) | ★★★★☆ | Quantifies ephemeris uncertainty limits |
| SRQ2 (Celestial Body Performance) | ★★★★☆ | Mars-specific accuracy analysis |
| SRQ3 (Multi-body Fix Error) | ★★★☆☆ | Error propagation context |
| SRQ4 (Computational Efficiency) | ★★☆☆☆ | Integration methods discussed |
| H1 (Accuracy Hypothesis) | ★★★★★ | **CRITICAL** - Defines fundamental accuracy limits |

### Abstract Summary
Accuracy limits in inner planet ephemerides from asteroid mass uncertainties are investigated using Monte Carlo simulations. Key findings:
- Asteroid perturbations can reach several kilometers for Mars
- Present-day ranging measurements accurate to ~10 meters
- Asteroid mass modeling would need 1% accuracy to match observations
- With ranging data fully weighted: orbit distortion up to 5 km in RA/Dec
- With ranging de-weighted to VLBI level: uncertainties of 2-3 km
- Sub-kilometer accuracy possible only for short-term extrapolation

### Three Types of Ephemeris Uncertainties

1. **Relative angles and distances** between bodies
   - Determined by ranging measurements over sufficient geometry change

2. **Orientation onto external reference frame** (ICRF)
   - Determined by VLBI measurements of spacecraft

3. **Mean motions with respect to inertial space**
   - Most affected by asteroid perturbations
   - **PRIMARY SOURCE OF EPHEMERIS ERROR**

### Asteroid Mass Sources

| Method | Application | Accuracy |
|--------|-------------|----------|
| Direct dynamical | Big 3 (Ceres, Pallas, Vesta) | Best available |
| Star occultations | Selected asteroids | Direct diameter |
| Radar echo | Selected asteroids | Direct diameter |
| IRAS survey | Most large asteroids | σ/D typically 5-30% |
| Ground-based photometry | Remaining asteroids | Requires albedo assumption |

### Taxonomic Classification System

| Class | Type | Density (g/cm³) | σ |
|-------|------|-----------------|---|
| C | Carbonaceous | Variable | ±20% |
| S | Silicate | Variable | ±20% |
| M | Metallic | Variable | ±20% |
| V | (Vesta-like) | 3.44 | ±0.12 |
| B | (B-type) | 2.71 | ±0.11 |
| G | (G-type) | 2.12 | ±0.04 |

### Monte Carlo Methodology

#### Differential Equation of Motion (Eq. 1):
```
δr̈ = -G(m + m₀)[δr/r³ - 3r(r·δr)/r⁵] + Σᵢ Gmᵢ[ρᵢ'/|ρᵢ'|³ - rᵢ/rᵢ³]

where:
- δr = difference between perturbed and base ephemeris
- m = central body mass
- m₀ = perturbed body mass
- mᵢ = mass of ith perturbing asteroid
- ρᵢ' = rᵢ - r - δr
```

#### Mass Uncertainty Simulation:
```
For each asteroid:
1. Generate Gaussian distribution N(μᵣ, σᵣ) for diameter
2. Generate Gaussian distribution N(ρ, 0.2ρ) for density
3. Compute mass: M = (4/3)π(D/2)³ × ρ
4. Create δM set of mass perturbations
```

### Observational Data Used (Table 1)

| Period | Type | Planet | n | σ (Fit 1) | σ (Fit 2) |
|--------|------|--------|---|-----------|-----------|
| 1990-94 | Magellan VLBI | Venus | 18 | 3 km | 3 km |
| 2001-03 | MGS+Odyssey VLBI | Mars | 34 | 1 km | 1 km |
| 1976-81 | Viking Lander range | Mars | 1282 | 10 m | 1 km |
| 1997 | Pathfinder range | Mars | 90 | 10 m | 1 km |
| 1998-2001 | MGS range | Mars | 200 | 10 m | 1 km |

### Key Results: 100 Monte Carlo Runs

#### Pre-fit Residuals:
- Initial value: 0 at epoch JED 2440400.5 (1969 Jun 28)
- Growth: 2+ km over few decades
- Spikes at Mars oppositions (geometric)
- Declination/range reflect RA perturbations

#### Post-fit 1 (Full ranging weight):
- **Range tightly fit** during Viking (1976-82) and MGS (1997-2003)
- **RA and Dec WORSE than pre-fit**
- Ephemeris distorts to accommodate highly accurate ranging
- Result: **Up to 5 km errors in RA/Dec**

#### Post-fit 2 (Ranging de-weighted to VLBI level):
- Residuals small during VLBI interval (1990-2003)
- Larger errors for earlier decades
- Result: **2-3 km uncertainty overall**

### Critical Finding: Accuracy vs Time Trade-off

```
SHORT-TERM (years): Sub-kilometer accuracy achievable
- Use only recent data with full weight
- VLBI + ranging establish current ephemeris

LONG-TERM (decades): 2-5 km uncertainty unavoidable
- Asteroid mass uncertainties accumulate
- Cannot connect ephemeris over 20 years at measurement accuracy
- Mis-modeled forces cause orbit distortion
```

### IRAS-FBM Comparison

Comparison between IRAS and Free Beaming Model (Hasegawa 1999):
- IRAS diameters appear under-estimated for small objects
- FBM generally gives larger diameters
- Demonstrates observational/model bias in diameter estimation
- Positive diameter bias → positive RA drift
- Negative bias → negative RA drift

### Sample Asteroids from Monte Carlo

| Asteroid | Diameter (km) | 1σ (km) | σ/D |
|----------|---------------|---------|-----|
| (5) Astraea | 119 | 7 | 6% |
| (127) Johanna | 42 | 2 | 5% |
| (875) Nymphe | 14 | 1 | 7% |
| (683) Lanzia | 82 | 22 | 27% |
| (393) Lampetia | 97 | 31 | 32% |

### Implications for Ephemeris Improvement

1. **Best known asteroid masses** should be carefully selected
2. **Solving for individual masses** not yet realistic (need longer accurate observations)
3. **Independent mass determinations** difficult (spacecraft rare, dynamical poor success)
4. **Frequency grouping** possible alternative to taxonomic grouping
5. **Massive ring** for small asteroids (Krasinsky et al. 2001) - mass comparable to Ceres

### Implications for Celestial Navigation

#### Ephemeris Accuracy Context:
```
Mars ephemeris uncertainty: 2-5 km (over decades)
At Mars distance (~1.5 AU): 
  5 km ≈ 0.02 arcsec angular error

For navigation using Sun/Moon (much closer):
  Ephemeris errors MUCH smaller relative to user
  
Sextant precision: 0.1-0.3 arcmin = 6-18 arcsec
Ephemeris contribution: <0.1 arcsec

CONCLUSION: Asteroid-induced ephemeris errors 
are negligible for celestial navigation
```

#### Why This Matters for SRQ1/H1:
1. Confirms ephemeris is NOT the limiting factor in celestial fix accuracy
2. For Sun/Moon (primary navigation bodies): errors sub-arcsecond
3. Dominant error sources remain: sextant, refraction, chronometer
4. Validates using Skyfield/DE440 without ephemeris uncertainty concern

### Relationship to Park et al. 2021

| Aspect | Standish & Fienga 2002 | Park et al. 2021 |
|--------|------------------------|------------------|
| Focus | Uncertainty sources | Achieved accuracy |
| Method | Monte Carlo simulation | Observational residuals |
| Result | 2-5 km Mars uncertainty | Sub-meter to 13m by body |
| Era | DE ephemerides ~2000 | DE440/DE441 (2021) |
| Value | Explains limitations | Documents improvements |

**Together these papers establish:**
- Ephemeris errors are well-characterized
- Primary limitation is asteroid masses (for Mars)
- For Sun/Moon: negligible errors for navigation
- 20 years of improvement from DE~2000 to DE440

---

## Article #40: Ross 1994 - Minimizing Errors in Celestial Positioning

### Bibliographic Information
- **Author:** Paul F. Ross
- **Title:** Minimizing Errors in Celestial Positioning
- **Journal:** NAVIGATION: Journal of The Institute of Navigation
- **Year:** 1994, Volume 41, Number 3, Pages 251-264
- **DOI:** 10.1002/j.2161-4296.1994.tb01879.x
- **Institution:** Charles River Power Squadron, Lexington, Massachusetts
- **Classification:** PEER-REVIEWED

### Relevance Assessment
| Criterion | Rating | Notes |
|-----------|--------|-------|
| SRQ1 (Sight Reduction Accuracy) | ★★★★★ | **CRITICAL** - Quantifies all sextant error sources |
| SRQ2 (Celestial Body Performance) | ★★★☆☆ | General to all bodies |
| SRQ3 (Multi-body Fix Error) | ★★★★☆ | Error propagation to fix |
| SRQ4 (Computational Efficiency) | ★★☆☆☆ | Focus on observation not computation |
| H1 (Accuracy Hypothesis) | ★★★★★ | **ESSENTIAL** - Defines achievable accuracy limits |

### Abstract Summary
Statistical procedures for reducing systematic and random errors in sextant altitude observations. Key claims:
- Error reduction by factor of **3 to 8** compared to standard practices
- Probable error reducible from ~10 nmi to 1-3 nmi
- Statistical approach missing from Bowditch and Dutton
- Practical value for yachtsmen and professional mariners

### Circle of Probable Error Concept

**Definition:** If 100 fixes taken from stationary position:
- Circle of probable error (centered at true position) encloses ~50 fixes
- Other 50 fixes lie outside the circle
- May be elliptical in specialized situations

**Note:** Chapter XXIX on navigational errors from 1964 Bowditch was **dropped** from 1981/1984 editions.

### Error Sources in Line of Position (Table 1)

| Source of Error | Estimated Typical σ (arcmin) |
|-----------------|------------------------------|
| 1. Height of eye (dip) | 0.1 |
| 2. **Index error** | **2.0** |
| 3. Refraction (>15° altitude) | 0.1 |
| 4. Parallax | 0.1 |
| 5. Semidiameter | 0.1 |
| 6. **Personal random error** | **0.5** |
| 7. Timing error | 0.3 |
| **Total (RSS)** | **2.1** |

**Critical Finding:** Index error is:
- 4× larger than personal random error
- 20× larger than other error sources
- **PRIMARY TARGET FOR ERROR REDUCTION**

### Conversion Factors
```
1 arcmin error in observed altitude = 1 nmi error in LOP
σ = 2.1 arcmin → Probable error = 2.1 × 0.6745 = 1.4 arcmin = 1.4 nmi

At equator: 4 seconds timing error = 1 nmi position error
```

### Index Error and Temperature Effects

#### Key Discovery:
Index error changes **rapidly** with temperature changes.

#### Experimental Data (Davis Mark 15 sextant):

**October Series (Fig. 1):**
- Sextant at 60°F, air temp 80°F
- First reading: +2.8 arcmin
- After 9 min: -1.3 arcmin
- After 30 min: -0.3 arcmin
- **Change of 4 arcmin in first 5 minutes**

**March Series (Fig. 2):**
- Sextant at 70°F, air temp 50°F
- Similar rapid change pattern

**September Series (Fig. 3):**
- Even 5°F cooling (75→70°F) over 2 hours produced time trend

#### Three-Material Hypothesis:
Three slopes in temperature response may relate to:
- Plastic (frame)
- Glass (mirrors)
- Steel (arc/mechanism)

### Personal Random Error Analysis

**Method:** Fit straight line to observations 10-30 in temperature series, measure deviation from line.

**Result:** σ = 0.5 arcmin for personal random error

**Comparison to Bowditch claim:**
- Bowditch: "well-constructed sextant capable of 0.1' instrument error"
- Actual personal random error: 0.5 arcmin
- Personal error **5× larger** than instrument precision

### Prismatic Error of Filters

#### Discovery:
Most navigators assume filter prismatic error is zero. **This is wrong.**

#### Measurement Procedure:
1. Observe index error 3-7 times WITHOUT filter
2. Insert Filter A
3. Observe index error 3-7 times WITH filter
4. Remove Filter A
5. Observe index error 3-7 times WITHOUT filter
6. Prismatic error = Median(step 3) - Mean(Median step 1, Median step 5)

#### Results for Davis Mark 15 (Table 2):

| Filter | Color | Path | Prismatic Error (arcmin) |
|--------|-------|------|--------------------------|
| a | gray | horizon glass | -1.4 |
| b | orange | horizon glass | -2.1 |
| c | blue | horizon glass | -0.5 |
| d | gray (1st) | index mirror | -1.1 |
| e | orange (2nd) | index mirror | -0.4 |
| f | blue (3rd) | index mirror | +0.2 |
| g | blue (4th) | index mirror | +0.3 |

**Validation:** Multiple filters combined show cumulative prismatic error matching predictions (Table 3).

### Recommended Procedures to Minimize Error

#### 1. Mirror Adjustment:
- Check and adjust mirrors at beginning of each season
- Leave unchanged for rest of season

#### 2. Drum Rotation:
- **Always rotate in same direction** (bringing body down to horizon)
- Reduces gear backlash error
- Note: Dutton incorrectly advises "up and down" for index error

#### 3. Temperature Awareness:
- Keep sextant shaded from direct sunlight
- Learn your sextant's temperature sensitivity
- Test with 20-40°F temperature change

#### 4. Sight-Taking Procedure:
```
a) Observe index error 3-5 times (with all filters in place)
   → Choose median value
   
b) Immediately observe sextant altitude 3-5 times
   → Record time for each shot
   → Choose median altitude with its time
   OR: Plot 5-7 altitudes vs time, fit line, choose best point
   
c) Immediately observe index error 3-5 times again
   → Choose median value
   
d) Calculate index error for shot time by interpolating
   between medians from (a) and (c)
```

#### 5. Filter Correction:
- If possible, observe index error WITH filters in place
- If not possible, prepare filter prismatic error table
- Apply corrections based on which filters used

#### 6. Timing:
- Observe time to nearest second
- Cannot take repeated measures of time
- 4 seconds error at equator = 1 nmi error

### Error Reduction Achievable

| Error Source | Before | After | Reduction |
|--------------|--------|-------|-----------|
| Index error variation | 2.0 arcmin | 0.1 arcmin | 20× |
| Personal random error | 0.5 arcmin | 0.3 arcmin | 1.7× |
| **Total error (RSS)** | **2.1 arcmin** | **0.5 arcmin** | **4×** |

**Positioning accuracy improvement:** From ~10 nmi to 1-3 nmi probable error

### Critical Quote on Education

> "The education of navigators needs to be revised to teach the underlying principles of statistical errors in measurement, the need to measure each error component, and the benefit of working to minimize the errors contributing most to errors in positioning."

### Implications for Algorithm Development

#### For SRQ1/H1 (Accuracy Validation):
```
Error Budget for Celestial Fix:

Algorithm errors (ephemeris, computation): < 0.1 arcmin
Sextant/observation errors: 0.5-2.1 arcmin

DOMINANT ERROR SOURCE: Observation, not algorithm

Implication: Algorithm accuracy of 0.1 arcmin or better
is sufficient - further improvement has diminishing returns
```

#### For SRQ3 (Multi-body Fix):
```
With n observations:
- Personal random error reduces by √n
- Systematic errors (index, filter) must be corrected individually
- Least squares fitting appropriate for random errors
- Systematic errors require proper calibration
```

#### Error Model for Simulation:
```python
# Based on Ross 1994 Table 1
error_sources = {
    'height_of_eye': 0.1,      # arcmin σ
    'index_error': 2.0,         # arcmin σ (uncorrected)
    'index_error_corrected': 0.1,  # arcmin σ (with procedures)
    'refraction': 0.1,          # arcmin σ
    'parallax': 0.1,           # arcmin σ
    'semidiameter': 0.1,       # arcmin σ
    'personal_random': 0.5,    # arcmin σ (single shot)
    'personal_median_of_5': 0.3,  # arcmin σ (median of 5)
    'timing': 0.3              # arcmin σ (1 second error)
}

# Total error (RSS)
import numpy as np
total_uncorrected = np.sqrt(sum(v**2 for k, v in error_sources.items() 
                                 if 'corrected' not in k and 'median' not in k))
# ≈ 2.1 arcmin

# With corrections
total_corrected = np.sqrt(0.1**2 + 0.1**2 + 0.1**2 + 0.1**2 + 0.3**2 + 0.3**2)
# ≈ 0.5 arcmin
```

### Relationship to Other Error Sources

| This Paper | vs Ephemeris (Park 2021, Standish 2002) |
|------------|----------------------------------------|
| **Observation errors** dominate | Ephemeris errors negligible |
| 0.5-2.1 arcmin σ | <0.001 arcmin contribution |
| Reducible by procedures | Already optimized |
| Navigator responsibility | Algorithm handles automatically |

**Conclusion:** Ross 1994 quantifies the **dominant error sources** in celestial navigation. The algorithm's ephemeris/computation accuracy is far superior to observation accuracy. Error reduction comes from **proper observation procedures**, not algorithm improvement.

---

### Article #41
**Citation:** Gordon, R. B. (1964). The Attainment of Precision in Celestial Navigation. *The Journal of Navigation, 17*(2), 125–147.

**Relevance:** Controlled experimental study separating error sources through 500 observations - foundational quantitative error analysis

#### Study Design:
- **Location:** Long Island Sound, Connecticut (shore and sea observations)
- **Position:** 41°16'06"N, 72°40'04"W
- **Sample size:** ~500 observations
- **Instruments:** Hughes 6-inch sextant (10" graduation), Mark II Navy sextant
- **Telescopes tested:** 2.5×, 3×, 8× magnification
- **Methodology:** Compare observed vs computed altitudes under controlled conditions

#### Index Error Determination:
| Horizon Quality | σ (arcmin) |
|-----------------|------------|
| Excellent/Very Good | 0.28 |
| Good/Fair | 0.35 |
| Poor | 0.50 |

**Key finding:** "Some or all of what is called 'personal error' in earlier studies is, in fact, variable index error"

**Methods comparison:**
1. Star superposition: Poor (5' range due to diffraction)
2. Horizon method: σ = 0.28-0.50' depending on horizon quality
3. Sun tangency: Most precise, difficult on moving vessel

#### Random Errors - Daylight Observations (Standard Telescope):
| Conditions | σ (arcmin) | R/n (rejected) |
|------------|------------|----------------|
| **Ideal conditions** | **0.33** | 2/39 |
| Cloud cover only | 0.32 | 1/15 |
| Motion - yawing at anchor | 0.32 | 0/9 |
| Motion - reaching Force 3 | 0.31 | 1/16 |
| Motion - beating Force 4 | 1.6 | 0/6 |
| Horizon - good | 0.52 | 0/6 |
| Horizon - poor | 0.70 | 0/38 |
| Poor horizon + motion | 0.43-0.44 | variable |
| **Moon (ideal)** | **0.40** | 0/19 |

**Critical insight:** "Horizon quality is the most important single factor influencing the standard deviation of daylight observations"

#### Random Errors - Twilight Observations:
| Condition | σ (arcmin) |
|-----------|------------|
| Star sights (average) | 0.58 |
| Range | 0.24 - 1.6 |

**Twilight parameter t*:**
- t* = 0: sunset
- t* = 1.0: end of nautical twilight
- Practical limit with 2.5× telescope: t* ≈ 0.55
- With 8×30 monocular: useful observations to t* > 1.0 (after nautical twilight)

#### Telescope Comparison:
| Telescope | σ Ideal | σ Poor Horizon | Notes |
|-----------|---------|----------------|-------|
| 2.5× (standard) | 0.33 | 0.70 | Best overall |
| 3× Galilean | 0.52 | 0.22 | Inferior optics |
| 8× Monocular | 0.37 | 0.92 | High quality optics |

**Key finding:** "The use of high-power telescopes does not offer any substantial advantage in daylight observations but ... it is desirable to have a high degree of optical quality"

#### Systematic Errors (Δh):
| Conditions | Telescope | Δh (arcmin) | Range |
|------------|-----------|-------------|-------|
| **Ideal** | 2.5× or 8× | **0.11 ± 0.12** | +0.08 to -0.15 |
| Cloud cover only | 2.5× | +0.5 ± 0.1 | - |
| Motion (reaching) | Various | +0.2 to +0.9 | - |
| Full visibility, poor horizon | 2.5× | 0.7 ± 0.12 | +0.6 to -0.9 |
| Restricted visibility | 2.5× | **5.7 ± 1** | +6.1 to -6.9 |
| Restricted visibility | 8× | **7.4** | +9.0 to -5.9 |
| Abnormal refraction | 2.5× | **1.0 ± 0.1** | +2.2 to -0.1 |

**Critical insight:** Under restricted visibility, systematic errors dominate (up to 9 arcmin) due to:
1. Uncertainty in distance to apparent horizon
2. Difficulty judging tangency
3. Abnormal dip corrections

#### Atmospheric Refraction Effects:
- **Normal conditions:** Air temperature > water temperature → negligible error
- **Abnormal conditions:** Air temperature < water temperature → significant errors
- **Observed dip departure:** Up to 6.6' from tabulated values
- **Mechanism:** Heated air layer near surface creates "boiling" effect visible through telescope

**Refraction formula cited:**
```
γ = (16"17/T) × tan ζ - 0.07 tan ζ sec² ζ × B
```
where ζ = apparent zenith distance, B = pressure (millibars), T = absolute temp (°C)

#### Twilight Systematic Errors:
| Star | t* | Δh (arcmin) |
|------|-----|-------------|
| Arcturus | 0.21 | +5.7 |
| Arcturus | 0.52 | +5.4 |
| Altair | 0.55 | +6.8 |
| Altair | 1.01 | +8.13 |
| Altair | 1.17 | +10.3 |

**Trend:** Δh increases as twilight deepens - residual error ~3' due to visual difficulty judging dimly illuminated horizon

#### Moonlight Effects:
| Body | Azimuth rel. to Moon | Δh - ΔD (arcmin) |
|------|---------------------|------------------|
| Saturn | 18° | +3.1 |
| Altair | 101° | +3.0 |
| Moon | 0° (on moonpath) | +4.0 |

**Insight:** Moonlight on horizon lowers apparent horizon by ~1 arcmin

#### Practical Recommendations:
1. **Use computed altitude from sine-cosine formula** rather than H.O. 214 triple interpolation
   - Direct formula: sin h = sin φ sin d + cos φ cos d cos t
   - More accurate than tabular methods

2. **Balanced observations on reciprocal azimuths** eliminate systematic dip errors

3. **Night sextant design proposed:**
   - Small mirrors (only bright stars needed)
   - Remove horizon glass frame
   - Large objective night monocular
   - Filtered red light for arc reading

4. **Index error checks:** 25 readings on poor horizon needed for ±0.26' accuracy at 99% confidence

#### Comparison with Other Studies:
| Source | Method | σ Reported |
|--------|--------|------------|
| Smiley 1951 | 3× telescope | 0.37' (probable error 0.25') |
| Shufeldt 1962 | 20× telescope | 0.05-0.26' |
| Shufeldt 1962 | 6× telescope | 0.09-0.89' |
| **Gordon 1964** | 2.5× (ideal) | **0.33'** |
| **Gordon 1964** | 8× (ideal) | **0.37'** |

**Discrepancy explained:** Superior optical quality, not magnification, accounts for Shufeldt's better results with high-power telescopes

#### Minimum Achievable Error Analysis:
| Source | Contribution |
|--------|-------------|
| Sextant graduation (10") | ±0.06' |
| Timing error (±4s) | ±0.1' |
| Index error σ | 0.28-0.50' |
| **Theoretical minimum** | **~0.15'** |

#### Content for Paper Sections:

**For Error Analysis:**
- Random error under ideal conditions: σ = 0.33' (Sun), 0.40' (Moon), 0.58' (stars)
- Horizon quality is dominant factor for daylight observations
- Systematic errors can reach 7-9' under restricted visibility

**For Algorithm Validation:**
- Target accuracy of 0.1-0.2' computation error is well below observation limit
- Algorithm errors negligible compared to observation errors (0.33-0.70')

**For Discussion:**
- Telescope optical quality more important than magnification
- Balanced observations on reciprocal azimuths eliminate systematic dip errors
- Night observations possible with proper instrument design

#### Key Quotes:
- "Under very good or ideal conditions of observation during daylight, altitudes can be measured to an accuracy comparable with that of which the instrument itself is capable"
- "It is doubtful if the improvement attained with a better instrument would be worth the additional trouble and expense involved"
- "The standard deviation of a sextant altitude made in daylight is 0'4 and increases only slightly under poor observing conditions"
- "In careful navigation, the average of several (preferably five or more) observations should always be used"

#### Relationship to Ross 1994:
| Gordon 1964 | Ross 1994 |
|-------------|----------|
| Controlled experiments | Statistical analysis |
| σ = 0.33' ideal daylight | σ = 0.5' personal random |
| Index error σ = 0.28-0.50' | Index error σ = 2.0' (uncorrected) |
| Horizon quality dominant | Index error dominant |
| 500 observations | Different methodology |

**Synthesis:** Gordon's controlled conditions achieve lower errors than Ross's practical conditions. Combined, they bracket the expected error range: 0.3-2.1 arcmin depending on procedures and conditions.

---

### Article #42
**Citation:** Shufeldt, H. H. (1961). *Report on Precision Celestial Navigation Experiments.* Technical Report, Contract N Onr-2449(00). Office of Naval Research.

**Relevance:** Foundational precision sextant study (~5,000 observations) comparing telescope magnifications - the study Gordon 1964 referenced

#### Study Design:
- **Period:** January 1958 - May 1961
- **Locations:** Key West, FL; USS VALLEY FORGE (CVS-45); USS TUTUILA (ARG-4); Yachts PARAKOUU and DOUBLOON
- **Sample size:** ~5,000 celestial observations
- **Principal investigators:** H. H. Shufeldt, Capt. P.V.H. Weems USN (Ret.), CDR William S. Brown USN
- **Sponsor:** Office of Naval Research

#### Equipment Tested:
**Sextants:**
- U.S. Navy Mark II Type A (Aluminum) with 3× telescope
- Hughes "Gothic" (brass frame, ~4 lb 10 oz)
- 3× Plath sextants (custom-built for study)
  - Large micrometer drums reading to 0.1'
  - Vernier readable to 0.05'
  - Arc graduation error <10" (vs 30" acceptable in US manufacture)
  - Large index/horizon mirrors matching 50mm objectives
  - Two alloy frames (~3 lb 12 oz), one brass (~4 lb 10 oz)

**Telescopes:**
| Power | Type | Objective | Field | Application |
|-------|------|-----------|-------|-------------|
| 3× | Erect | Standard | Wide | Navy standard (inferior) |
| 6× | Prismatic | 50mm | Medium | Comparison baseline |
| 7× | Prismatic | 50mm | 7° | Star observations at dusk |
| 16× | Prismatic | 50mm | ~4° | High power alternative |
| 20× | Prismatic | 50mm | ~3° | Optimal for daylight |

**Dip Meters:**
- Gavrisheff dip meters with 6× and 16× telescopes
- Reading to 0.2' with interpolation to 0.1'
- Index error cancelled by 180° rotation about optical axis

#### Key Findings - Telescope Comparison:

**Summary of Findings (quoted):**
1. "The 20 x 50 mm prismatic monocular of good manufacture, is superior to any sextant telescope of lesser power during daylight hours"
2. "The 7 x 50 mm prismatic monocular is superior to any other telescope available for test when the horizon is poorly defined, due to darkness"
3. "A dip meter, fitted with a high power telescope, and capable of determining empirically the value of the dip to one-tenth of a minute of arc is an essential adjunct to the sextant, where accuracy is of primary importance"

#### Quantitative Results - Random Errors:

**Observations at Sea:**
| 20× Telescope | 6× Telescope | Body |
|---------------|--------------|------|
| 0.042' | 0.073' | Sun |
| 0.046' | 0.113' | Sun |
| 0.050' | 0.138' | Sun |
| 0.054' | 0.192' | Sun |
| 0.061' | 0.108' | Sun |
| 0.073' | 0.222' | Sun |
| 0.157' | 0.364' | Sun |

**Observations from Beach:**
| 20× Telescope | 6× Telescope | Body |
|---------------|--------------|------|
| **0.027'** | 0.060' | Sun |
| 0.031' | 0.089' | Venus |
| 0.036' | 0.078' | Sun |
| 0.058' | 0.054' | Moon |
| 0.063' | 0.127' | Sun |
| 0.077' | 0.227' | Sun |
| 0.083' | 0.123' | Sun |
| 0.085' | 0.323' | Sun |
| 0.092' | 0.273' | Sun |
| 0.119' | **0.716'** | Sun |
| 0.127' | 0.196' | Sun |
| 0.130' | 0.200' | Sun |
| 0.150' | 0.231' | Sun |
| 0.162' | 0.265' | Sun |
| 0.173' | 0.327' | Sun |
| 0.178' | 0.112' | Sun |
| 0.204' | 0.373' | Sun |

**Summary Statistics:**
- 20× best result: **0.027'** mean aberration
- 20× typical range: 0.042' - 0.204'
- 6× best result: 0.054'
- 6× typical range: 0.073' - 0.716'
- **20× consistently 2-5× better than 6×**

#### Dip Meter Results:

**USS VALLEY FORGE observations (January 1959):**
| Date | Conditions | Dip by Meter | Dip by Table | Difference |
|------|------------|--------------|--------------|------------|
| 3 Jan | O-5 level, Air 48°, Water 60° | 7.25' | 8.9' | -1.65' |
| 5 Jan | O-6 level, Air 53°, Water 68° | 9.5' | 8.9' | +0.6' |
| 11 Jan | O-6 level | 9.625' | 8.9' | +0.725' |
| 9 Jan | O-6 level, Air 74°, Water 70° | 8.7' | 8.9' | -0.2' |

**Position accuracy comparison (7 miles from Bermuda):**
- With dip meter correction: **0.115 nmi** mean error
- With tabulated dip: **0.85 nmi** mean error
- **Improvement factor: 7.4×**

**Diurnal variation discovered:**
- USS TUTUILA, August 1960: 1' diurnal variation in dip
- Dip 1' greater in afternoon/evening than morning
- Persisted over several days without temperature changes

#### Star Observation Timing:

**Bright twilight advantage with 20× telescope:**
- 5 stars observed within 10 minutes of sunset (30 Aug 1960):
  - Altair (mag 0.9), Alt 46°, Az 110°
  - Antares (mag 1.2), Alt 34°, Az 194°
  - Arcturus (mag 0.2), Alt 50°, Az 265°
  - Deneb (mag 1.3), Alt 45°, Az 054°
  - Vega (mag 0.1), Alt 69°, Az 160°

- "With a 20x telescope, the observations for the evening star fix were completed at about the time a navigator using the ordinary sextant and telescope would start making his observations"

**7×50 superiority at dim horizon:**
- "Under such conditions of poor horizon illumination, the 7 x 50 monocular was found very considerably superior to any other telescope tested"

#### Daytime Star Observations:

Attempted but unsuccessful with conventional sextant design due to mirror light loss:
- Sirius observed 15 min before sunset
- Arcturus observed 25 min after sunrise
- Recommended sextant redesign: direct star view with reflected horizon

#### Reduction Method Recommendations:

**For precision work:**
- Use Ephemeris, not Nautical Almanac (accuracy limitation noted)
- Direct computation: `sin h = sin L sin d + cos L cos d cos t`
- Six-place natural functions preferred over logarithms
- Reduce from estimated position, not assumed position
- Smart (1919) showed 30' position error can cause 1.0 nmi fix error at 75° altitude

**Nautical Almanac accuracy limits (quoted):**
- GHA/Dec error: up to 0.2' (0.25' for Sun GHA, 0.3' for Moon)
- "Only one-third of values will have errors >0.05'"
- Insufficient for <0.2 nmi position accuracy

#### Position Accuracy Conclusions:

| Equipment | Achievable Accuracy |
|-----------|--------------------|
| 20× telescope + dip meter | **~0.25 nmi** |
| 6× and 3× telescopes (conventional) | ~0.4 nmi |
| **Improvement** | **~40%** |

#### Observation Technique Recommendations:

1. **String size:** 10-13 observations per string (accuracy deteriorates after due to eye fatigue)
2. **String duration:** 2.5-3.5 minutes
3. **Rest period:** Brief rest between strings restores accuracy
4. **Plotting scale:** 1" = 10 seconds time, 1" = 1' altitude
5. **Observer selection:** "The observer should be selected for his visual acuity; he need have no knowledge of navigation"
6. **Expert marksmen** often make excellent observers

#### Equipment Recommendations:

**Sextant improvements:**
1. Maximum constant error: 5" (not 30" as US standard)
2. Front-coated mirrors to reduce light loss
3. Remove clear glass portion of horizon mirror
4. Rubber eye cups essential for high-power telescopes
5. Remote altitude/time readout system recommended

**Recommended telescope configuration:**
- Daylight: 20×50 prismatic monocular
- Twilight (dim horizon): 7×50 prismatic monocular
- Interchangeable eyepieces on single 50mm objective

#### Content for Paper Sections:

**For Error Analysis:**
- Best achievable random error: 0.027' with 20× telescope
- Typical improvement: 2-5× better than conventional 6× telescope
- Dip meter reduces systematic error by factor of 7

**For Algorithm Validation:**
- Algorithm accuracy target of 0.1' is well within the observation capability
- Best observers achieve 0.027-0.073' mean aberration

**For Discussion:**
- Precision navigation to 0.25 nmi demonstrated in 1961
- Technology limited by Nautical Almanac accuracy, not observation capability
- Observer skill selection critical for precision work

#### Key Quotes:
- "The accuracy obtained with the 3x telescope were so obviously inferior to those obtained from a 6x telescope, that early in the study it was decided to use only the 6x glass for purposes of comparison"
- "Without exception, once the observers became accustomed to the new glasses, they preferred them to those of 6 power, when accuracy was particularly required"
- "The instruments are available for refined celestial navigation; data for the accurate and rapid reduction of sights should also be made available"
- "A ship's position may be fixed at sea, under good observational conditions, to about 0.25 miles, by a round of multiple observations of stars"

#### Relationship to Gordon 1964:
| Shufeldt 1961 | Gordon 1964 |
|---------------|-------------|
| ~5,000 observations | 500 observations |
| Multiple telescopes tested | Standard + 8× tested |
| 20× best: 0.027-0.073' | Standard σ = 0.33' |
| Dip meter essential | Dip meter used |
| Navy/commercial vessels | Shore + small boat |
| Operational conditions | Controlled conditions |

**Synthesis:** Shufeldt demonstrates that with specialized equipment (20× telescope + Gavrisheff dip meter), precision 2-10× better than conventional equipment is achievable. Gordon's controlled experiments with standard equipment establish the baseline that Shufeldt's specialized equipment exceeds.

---

### Article #43
**Citation:** Pinčetić, T., Lušić, Z., & Krančević, P. E. (2024). The Role of Celestial Navigation in Modern Day and Future Navigation. In *MT'24: 10th International Conference on Maritime Transport* (pp. 1–11). Barcelona, Spain.

**Relevance:** Current regulatory framework (STCW/IMO) for celestial navigation education and GPS backup justification

#### Regulatory Framework:

**STCW Requirements (Tables A-II/1 and A-II/2):**
- "User must be able to obtain and determine a position by means of celestial observation"
- Criteria: "Fix obtained by celestial observations within the accepted accuracy level"
- Required competencies:
  - Nautical astronomy
  - Sextant adjustment and reading
  - Sight reduction computation
  - Plotting lines of position
  - Meridian altitude of sun
  - Latitude by Polaris
  - Sunrise/sunset determination
  - Compass error by azimuth/amplitude
  - Star identification and selection
  - UT determination

**Critical Regulatory Inconsistency:**
- STCW: Requires deck officers to demonstrate celestial navigation proficiency
- SOLAS Chapter V, Regulation 19: Does **NOT** require vessels to carry a sextant
- "Inconsistency between competencies students are required to learn during schooling and those required in professional practice"

#### IMO Member Positions (2008 STCW Review):

| Country | Position |
|---------|----------|
| **Norway** | Proposed **canceling** celestial navigation as mandatory |
| **China** | Simplify: sun/stars only + electronic almanac + calculation software |
| **United States** | Maintain as GPS backup + compass error; pare down to "bare essentials" |
| **UK (2023)** | Modernize methodology; keep sextant handling; remove obsolete parts |

**Manila STCW Amendments (2010) Outcome:**
- Celestial navigation retained as essential competency
- May be excluded for restricted certificates (coastal voyages only)
- Primary role: backup system for satellites in open sea navigation

#### US Naval Academy History:
- 1998: Announced students would no longer be taught celestial navigation
- Proposed adding computer-based navigation instead
- **Idea was finally abandoned**

#### UK Proposed Modernization (2023):
1. **Keep essential:** Sextant handling techniques
2. **Modernize:** Twilight prediction, meridian passage, sunrise/sunset using celestial navigation software
3. **Remove:** Parts deemed obsolete

#### Method Classification (Dachev and Panov 2017):
1. **Traditional manual:** Handheld sextant + printed almanac + manual forms
2. **Traditional computer-based:** Sextant + software for sight planning/reduction
3. **Fully automated:** Automatic electronic sextant or star-tracking systems

#### Current Methods Referenced:
- **Most common:** Intercept method (Marcq St. Hilaire)
- **Other:** Sumner line, Dozier/direct method, longitude methods, ex-meridian method
- **Azimuth-only methods:** Based on isoazimuthal curves (imprecise with standard diopter)
- **Single-body methods:** Vertical angle + azimuth (rarely used in practice)

#### Modern Solutions Discussed:

**ECDIS Integration:**
- Most current ECDIS do not feature ephemerides or direct celestial positioning
- Drawing capabilities allow LOP plotting via intercept method
- "Most desired modern solution would be to upgrade ECDIS to solve problems in celestial navigation"
- GIS-based celestial fix demonstrated (Tsou 2020)

**Other Technologies:**
- Digital sextant (PSM Instrumentation)
- Artificial horizon aids
- Star tracker systems (military: Polaris, Poseidon, Trident, MX missiles)
- Smartphone apps (Celestial Navigator, CamSextant)
- Pulsar-based navigation (Adamson 2022)

#### SWOT Analysis:

| **Strengths** | **Weaknesses** |
|---------------|----------------|
| Independent and reliable positioning | Low position accuracy |
| Only real alternative to GNSS | Dependence on weather, sky cover |
| Low price and easy maintenance | Time of day, visible horizon |
| Simplified use with electronic aids | Lack of celestial bodies during day |
| | Instrument imperfection |
| | Complicated/demanding methods |
| | High price of new alternative systems |

| **Opportunities** | **Threats** |
|-------------------|-------------|
| Improvement of system integration | Development of alternative positioning systems |
| Improvement of ECDIS | Removal from STCW as mandatory competency |
| Development of military/commercial technologies | Fully autonomous vessels |
| Application of modern electronic aids/software | |
| Automatic star tracking systems | |

#### Key Arguments for Maintaining Celestial Navigation:

1. **GPS backup:** "GNSS liability implies a need for alternative means of fixing a position"
2. **War navigation:** Independent of electrical supply; no external signals required
3. **Cost-effective:** Low price and easy maintenance
4. **Proven reliability:** "Long-lasting proven success in use"
5. **Philosophy:** "Always have a reliable backup to anything, let alone means of navigation"

#### Educational Challenges:

- "Traditional methods thought by students as somewhat challenging in terms of nomenclature"
- "Mathematics involved is what is usually pointed out as challenging"
- Students perceive sextant as "old fashioned version of GPS"
- Teachers face challenge of "teaching a subject thought of as difficult and even as outdated"
- "Once students finish education, likely they will no longer use sight reduction practice"

#### Proposed Solutions:

1. **Visual methods:** Google Maps for demonstrating GP and COP concepts (Bensky 2017)
2. **Electronic aids:** ECDIS integration with celestial navigation software
3. **Simplified teaching:** "Reducing the methods and problems to the necessary minimum"
4. **Policy enforcement:** "Conventions, inspections and company policies require officers to perform some form of celestial navigation daily"

#### Content for Paper Sections:

**For Introduction:**
- STCW/IMO regulatory framework justifies continued relevance
- GPS backup motivation: "only real alternative to existing GNSS"
- Research addresses gap between required competencies and practical use

**For Literature Review:**
- Current regulatory status (STCW Tables A-II/1, A-II/2)
- Method classification (traditional manual → computer-based → automated)
- ECDIS integration as future direction

**For Discussion:**
- Algorithm contributes to "computer-based" method category
- Supports modernization trend in celestial navigation education
- Addresses gap between educational requirements and practical application

#### Key Quotes:
- "The only real alternative to existing GNSS systems is still celestial navigation"
- "The philosophy of any good seafarer is to always have a reliable backup to anything"
- "Teaching celestial navigation must adapt to the future"
- "Knowledge of the basics of celestial navigation... belong to classic nautical skills and should remain as such"

#### Relationship to Research:
| Paper Contribution | Research Application |
|-------------------|----------------------|
| STCW requires celestial competency | Justifies algorithm development |
| GPS backup motivation | Addresses research problem |
| ECDIS integration trend | Algorithm could integrate with ECDIS |
| Simplification needed | Python algorithm simplifies computation |

---

### Article #44
**Citation:** Zalewski, P., Bąk, A., & Bergmann, M. (2022). Evolution of Maritime GNSS and RNSS Performance Standards. *Remote Sensing*, 14(21), 5291. https://doi.org/10.3390/rs14215291

**Relevance:** Establishes GNSS performance baselines and limitations that justify celestial navigation as backup system; provides accuracy/integrity standards for comparison

#### Key IMO GNSS Performance Requirements (Res. A.915(22)):

| Parameter | Requirement | Notes |
|-----------|-------------|-------|
| **Accuracy** | 10 m (95%) horizontal | Most applications |
| **Alert Limit (AL)** | 25 m | Port approach/restricted waters |
| **Time to Alarm (TTA)** | 10 s | Maximum delay before warning |
| **Integrity Risk (IR)** | 10⁻⁵ per 3 h | Probability of undetected error |
| **Continuity** | 99.97% over 3 h | Service availability |
| **Availability** | 99.8% per 30 days | System uptime |

#### GNSS Accuracy Standards by System:

| System | Horizontal Accuracy (95%) | Standard |
|--------|---------------------------|----------|
| GPS | 100 m (HDOP=4) → 8 m global avg | MSC.112(73), SPS 2020 |
| GLONASS | 45 m (HDOP=4) | MSC.113(73) |
| DGPS/DGLONASS | 10 m | MSC.114(73) |
| GPS+GLONASS | 35 m non-diff, 10 m differential | MSC.115(73) |
| Galileo (L1) | 15 m horizontal, 35 m vertical | MSC.233(82) |
| Galileo (dual-freq) | 10 m horizontal, 10 m vertical | MSC.233(82) |
| Beidou | 25 m horizontal, 30 m vertical | MSC.379(93) |
| IRNSS | 25 m horizontal, 30 m vertical | MSC.449(99) |
| QZSS | 50.4 m (HDOP≤6.7) | MSC.480(102) |

#### GPS SPS Performance Standards (2020):

| Metric | Value | Condition |
|--------|-------|----------|
| Global average | ≤8 m (95%) horizontal | Representative user |
| Worst site | ≤15 m (95%) horizontal | Position solution available |
| Global average | ≤13 m (95%) vertical | Representative user |
| Worst site | ≤33 m (95%) vertical | Position solution available |

#### Key Definitions (Res. A.915(22)):

- **Accuracy:** Degree of conformance between estimated/measured position and true position
- **Integrity:** Ability to provide warnings within specified time when system should not be used
- **RAIM:** Receiver Autonomous Integrity Monitoring - uses redundant GNSS info to monitor signal integrity
- **CAIM:** Craft Autonomous Integrity Monitoring - cross-checks multiple navigation sensors
- **Alert Limit (AL):** Maximum allowable error before alarm triggered
- **Protection Level (PL):** Upper confidence bound on position error
- **Integrity Risk (IR):** Probability of undetected position error larger than AL

#### Error Model Components (TSE):

```
TSE² = NSE² + VTE²
where:
  NSE = Navigation System Error = CE + PE
  CE = Chart Error (surveying/chart generation inaccuracies)
  PE = Position Error (GNSS estimate error)
  VTE = Vessel Technical Error (difference from commanded position)
  TSE = Total System Error
```

#### Protection Level Calculations:

**Horizontal Protection Level (HPL):**
```
HPL = k_HPL × √[(σ_E² + σ_N²)/2 + √((σ_E² - σ_N²)/2)² + σ_EN²]
```

**Coverage factors for different integrity levels:**
| Probability | 1-p | k factor |
|-------------|-----|----------|
| 95% | 5×10⁻² | 2.45 |
| 99.99992% | 8.33×10⁻⁷ | 5.29 |

**Critical finding:** For 8m 95% accuracy, HPL at 99.99992% = 24.47m
- IF 95% error slightly exceeds 8.18m → HPL exceeds 25m AL
- Protection levels can exceed alert limits under normal conditions

#### Measured Data Example (Baltic Sea, 15-min observation):

| Parameter | Value |
|-----------|-------|
| Number of positions | 900 (1 Hz) |
| 1σ x coordinates | 4.17 m |
| 1σ y coordinates | 2.67 m |
| 1.73×DRMS (95%) | 8.55 m |
| HPL (99.99992%) | 26.16 m |
| Alert Limit | 25 m |
| **Result** | HPL exceeds AL |

#### GNSS Vulnerabilities Identified:

1. **Jamming:** Intentional RF interference
2. **Spoofing:** False GNSS signals to mislead receivers
3. **Multipath:** Signal reflections causing errors
4. **Solar activity:** Ionospheric effects on signals
5. **Satellite geometry:** High DOP conditions
6. **System faults:** Satellite, ground segment, receiver failures

#### SWOT Analysis of Maritime GNSS Standardization:

| **Strengths** | **Weaknesses** |
|---------------|----------------|
| 20 years of IMO GNSS policy research | Some requirements impossible to justify |
| Established standardization guidelines | Accuracy for ocean nav >> chart accuracy |
| SOLAS vessel navigation critical for safety | Lack of operation-specific requirements |

| **Opportunities** | **Threats** |
|-------------------|-------------|
| EU Space Programme support | Maritime community reluctance to change |
| e-Navigation and MASS initiatives | Lengthy IMO procedures |
| Multi-system receiver standards | Environmental/physical constraints |

#### Implications for Celestial Navigation Backup:

**Why GNSS backup is needed:**
1. Protection levels can exceed alert limits → integrity violations possible
2. Jamming/spoofing risks → system denial possible
3. Single points of failure → no redundancy without backup
4. 99.8% availability = 1.44 h/month downtime expected

**Celestial navigation advantages as backup:**
| GNSS | Celestial |
|------|----------|
| 10 m accuracy (95%) | ~1 nm typical, 0.25 nm precision |
| Subject to jamming/spoofing | Passive, no external signals |
| Requires satellites | Requires clear sky + horizon |
| Electronic, power-dependent | Can be manual/mechanical |
| Global coverage | Global coverage |

#### Accuracy Comparison Context:

| Source | Accuracy | Notes |
|--------|----------|-------|
| GNSS (GPS) | 8-15 m (95%) | Best case, no interference |
| DGPS | 10 m (95%) | With differential corrections |
| Celestial (conventional) | 0.5-2.1' obs → 1-2 nm fix | Ross 1994 |
| Celestial (precision) | 0.3' obs → 0.25 nm fix | Shufeldt 1961 |
| Algorithm target | 0.1' (<200 m LOP error) | Research goal |

**Key insight:** GNSS provides 8-15m accuracy in ideal conditions, but integrity monitoring shows real position uncertainty (HPL) can reach 25-50m. Celestial navigation provides 0.5-2 nm accuracy - complementary rather than competitive - serving as independent backup when GNSS unavailable or untrusted.

#### Content for Paper Sections:

**For Introduction:**
- GNSS accuracy standards (10m at 95%) establish electronic navigation baseline
- Integrity vulnerabilities (jamming, spoofing) motivate backup systems
- IMO recognizes need for resilient PNT (Position, Navigation, Timing)

**For Literature Review:**
- IMO performance standards evolution (MSC resolutions)
- Protection level vs alert limit framework
- RAIM integrity monitoring limitations

**For Discussion:**
- Celestial navigation fills GNSS backup gap
- Algorithm accuracy (0.1') → observation accuracy (0.3-2.1') dominates
- GNSS integrity uncertainties support maintaining traditional skills

#### Key Quotes:
- "The quality of shipborne GNSS... data forming the backbone of AIS and other global and regional remote monitoring and traffic service systems remains an open issue"
- "It is not possible to meet the requirements for any of the applications specified in the IMO resolution A.915(22) without SBAS support"
- "Neither GPS nor other GNSS systems... meet the IMO future requirements for satellite navigation without any problems"
- "The sample measurements taken in the Baltic Sea prove that HPL calculated according to (5) or (3) can easily exceed the AL of 25 m"

#### Relationship to Other Sources:

| This Paper | Related Source | Connection |
|------------|----------------|------------|
| GNSS 10m accuracy | Ross 1994 2.1' obs | GNSS more accurate but vulnerable |
| Protection level theory | Hoover 1984 | Confidence ellipse mathematics |
| IMO requirements | Lušić 2024 | Regulatory framework context |
| Backup system need | Garvin 2010 | Military navigator GPS concerns |

---

### Article #45
**Citation:** Dachev, Y., & Panov, A. (2017). 21st Century Celestial Navigation Systems. *Proceedings of the International Scientific Conference*, Nikola Vaptsarov Naval Academy, Varna, Bulgaria.

**Relevance:** Classification framework for celestial navigation methods; GPS vulnerability analysis; positions research algorithm in method taxonomy

#### Method Classification Framework:

| Category | Observation | Sight Reduction | Examples |
|----------|-------------|-----------------|----------|
| **Traditional, manual** | Handheld sextant | Printed almanac + manual forms | Classical navigation |
| **Traditional, computer-based** | Handheld sextant | Software for planning/reduction | **This research** |
| **Fully automated** | Electronic star tracker | Software + INS integration | AST-201, CIPP |

**Key insight:** This research falls into the "traditional, computer-based" category - maintaining human sextant observation while automating the computational sight reduction.

#### GPS Dependencies in Maritime Transportation:

| System | Positioning | Navigation | Timing |
|--------|-------------|------------|--------|
| AIS | ✓ | ✓ | ✓ |
| ECDIS | ✓ | ✓ | ✓ |
| EPIRB | ✓ | | |
| GMDSS | ✓ | ✓ | |
| Voyage Data Recorder | ✓ | ✓ | ✓ |
| Vessel Traffic Services | ✓ | ✓ | ✓ |
| Hydrographic survey | ✓ | ✓ | ✓ |

#### GPS Disruption Characteristics:

| Disruption Type | Intentional? | Predictable? | Local/Widespread |
|-----------------|--------------|--------------|------------------|
| **Spectrum encroachment** | Unintentional | Unpredictable | Local |
| **Solar weather** | Unintentional | Predictable | Widespread |
| **GPS infrastructure** | Either | Either | Widespread |
| **Jamming** | Intentional | Unpredictable | Local |
| **Spoofing** | Intentional | Unpredictable | Local |

#### Accuracy Comparison:

| Method | Angular Accuracy | Position Accuracy |
|--------|------------------|-------------------|
| Handheld sextant | 1-2 arcminutes | 1-2 nm |
| Star tracker (modern) | Sub-arcsecond (<1") | <30 m |
| GPS (standard) | N/A | 8-15 m |

#### Star Tracker Technology Examples:

**AST-201 (Lockheed Martin, 1998):**
- Size: 15 × 15 × 30 cm
- Weight: ~4 kg
- Field of view: 8.8°
- Detection limit: magnitude 7 (fainter than naked eye)
- Accuracy: several arcseconds
- MTBF: >700,000 hours
- Space qualified, no moving parts

**CIPP System (Rockwell Collins, 2015):**
- Size: 6.6 × 4.8 × 1.5 cm
- Weight: <100 g
- Accuracy: 0.1° pointing
- Power: <2 W peak, <0.3 W standby
- Update rate: 40-50 Hz
- Combines cameras + MEMS inertial sensors

#### INS/Celestial Complementary Characteristics:

| Characteristic | INS | Celestial |
|----------------|-----|----------|
| Self-contained after init | ✓ | |
| Link to external reference | | ✓ (fundamental) |
| Requires initial alignment | ✓ | |
| Completely autonomous | | ✓ |
| Accuracy degrades with time | ✓ | |
| Time-independent accuracy | | ✓ |
| Weather independent | ✓ | |
| Sensitive to clouds | | ✓ |
| Passive (no emissions) | ✓ | ✓ |
| Jam-proof | ✓ | ✓ |
| No shore/space dependency | ✓ | ✓ |

**Synergy:** INS provides continuous navigation but drifts; celestial provides periodic absolute fixes. Combined system is GPS-independent and jam-proof.

#### SOLAS Requirements Cited:

- "All ships irrespective of size to have a receiver for a global navigation satellite system or a terrestrial radionavigation system, or other means"
- ECDIS can satisfy chart carriage requirements
- Back-up arrangements required for electronic navigation

#### Key Arguments for Celestial as GPS Alternative:

1. **Independence:** "Celestial navigation remains one of the few independent alternatives to GPS"
2. **Single-point failure:** "GPS dependence... single-point-failure risk for safe navigation"
3. **No external dependency:** "Both INS and celestial are passive, jam-proof, and... not dependent on shore or space components"
4. **Historical proof:** "As much of a role to play in the future as it has in the past"

#### Sextant Accuracy Context:

- Observations accurate to "about 1-2 arcminutes (0.017-0.033 degrees)"
- "Most methods of sight reduction... incorporate approximations and non-rigorous assumptions as a means to simplify the computations"
- Implication: Algorithm can be rigorous since computational cost is now negligible

#### Content for Paper Sections:

**For Introduction:**
- Method classification positions research as "traditional, computer-based"
- GPS vulnerability analysis motivates backup system need
- SOLAS requirements establish regulatory context

**For Literature Review:**
- Three-category classification framework for celestial navigation
- GPS disruption taxonomy (jamming, spoofing, solar weather)
- Star tracker technology evolution

**For Discussion:**
- Algorithm addresses "computer-based" method category
- Rigorous computation now feasible (vs historical approximations)
- Complements potential future INS integration

#### Key Quotes:
- "Celestial navigation is often overlooked as an alternative to GPS because of the drawbacks of its traditional practice"
- "With reduced costs and enhanced reliability, such systems may be practical on many platforms not previously considered, including commercial and naval ships"
- "Independent alternatives to GPS are needed and are required"
- "The question what to do if GPS is not available is still unanswered firmly"

#### Relationship to Other Sources:

| This Paper | Related Source | Connection |
|------------|----------------|------------|
| Method classification | Lušić 2024 | Same 3-category framework |
| GPS vulnerabilities | Zalewski 2022 | Detailed GNSS standards |
| Sextant 1-2' accuracy | Ross 1994 | Error analysis confirms |
| Star tracker technology | Barbot 2022 | MARIS STELLA daytime system |
| INS integration | Yang 2022 | SINS/CNS integration |

---

### Article #46
**Citation:** Reid, T. G. R., Chan, B., Goel, A., Gunning, K., Manning, B., Martin, J., Neish, A., Perkins, A., & Tarantino, P. (2020). Satellite Navigation for the Age of Autonomy. *2020 IEEE/ION Position, Location and Navigation Symposium (PLANS)*. Xona Space Systems.

**Relevance:** Historical navigation accuracy evolution showing 10× improvement per 30 years; positions celestial navigation in technology lineage; autonomous vehicle requirements and GPS vulnerabilities

#### Historical Navigation Accuracy Trend:

| Era | Technology | RMS Accuracy | Driver |
|-----|------------|--------------|--------|
| Early 1900s | Celestial navigation | km-level | Maritime trade |
| 1940s | Ground-based radio (GEE, LORAN-C, DME) | ~100s m | WWII aviation |
| 1964 | Transit satellite | ~10s m | Cold War submarines |
| 1995 | GPS (SA removed 2000) | ~meters | Military precision |
| 2020s | LEO/enhanced GNSS | ~decimeters | Autonomous vehicles |

**Key finding:** "A tenfold improvement in location accuracy every thirty years"

#### Celestial Navigation Historical Context:

- "At the turn of the 20th century, celestial navigation was the state-of-the-art in location technology"
- "Developed in-part to help mariners find their way, it has been used at sea since at least the 16th century"
- "In the early 20th century, this technology was capable of delivering kilometer-level accuracy, sufficient for ships to get within sight of land"
- Weems 1951 cited for early 20th century accuracy assessment

#### Autonomous Vehicle Localization Requirements:

| Application | Lateral Accuracy | Longitudinal Accuracy | Reliability |
|-------------|------------------|----------------------|-------------|
| Highway driving | 50 cm | 50 cm | 99.999999% |
| City driving | 30 cm | 30 cm | 99.999999% |
| Full autonomy | 10 cm (95%) | 10 cm (95%) | ASIL D |

**ASIL D requirement:** One failure per billion miles driven

#### GPS/GNSS Vulnerabilities Documented:

**Jamming:**
- "More than 50,000 disruptions were recorded in Europe alone between 2016-2018"
- Motivated by privacy concerns (fleet tracking avoidance)
- Low-cost jammers readily available

**Spoofing:**
- "Just less than ten years ago, GNSS spoofing required specialized expertise and equipment costs of upwards of $50,000"
- "Now, with open source software and more accessible hardware, spoofing attacks can be accomplished for as little as $100"
- 2019 Geneva Motor Show: Multiple automakers (Audi, Peugeot, Renault, Rolls-Royce, VW, Daimler-Benz, BMW) reported GPS showing Buckingham, England in 2036
- 2019: Tesla Model S and Model 3 spoofing demonstrated by Regulus Cyber, creating "unsafe behavior of the autopilot"

#### LiDAR Limitations (Justifies Multi-Sensor Approach):

- Fog, rain, snow: 25% reduction in detection range
- 70% of US roads in snowy regions (FHWA)
- HD maps require frequent updates due to construction/seasonal changes
- Occlusions from large vehicles
- Sparse environments (open highways) yield poor localization

#### GNSS in Autonomous Driving Context:

**Strengths:**
- Unaffected by rain, snow, fog (microwave signals)
- Performs best in open sparse environments (highways)
- Complementary to LiDAR

**Baidu Apollo Framework Results:**
| Configuration | Alert Limit Achievement |
|---------------|------------------------|
| LiDAR only | 95% |
| LiDAR + IMU | 99.99% |
| LiDAR + IMU + precision GNSS | 100% (in test data) |

**Current GNSS Correction Services:**
- Continent-scale deployment of monitoring stations
- Delivery via cellular connectivity
- Accuracy approaching lane-determination requirements

#### LEO Navigation Concept (Future Context):

- 20-40× closer to Earth than MEO GNSS
- ~30 dB (1000×) less path loss → stronger signals
- Faster geometry change → rapid carrier phase resolution
- Encryption/authentication possible (no legacy constraints)
- ~300 satellites needed for GPS-equivalent coverage

#### Implications for Celestial Navigation:

1. **Historical baseline:** Celestial represents starting point of navigation evolution
2. **Accuracy context:** km-level (1900s) → 1-2 nm (modern sextant) shows improvement within celestial methods
3. **Backup justification:** GNSS vulnerabilities documented extensively
4. **Technology progression:** Each accuracy step required new infrastructure investment

#### Accuracy Comparison Table:

| Technology | Accuracy | Era |
|------------|----------|-----|
| Celestial (early) | ~km | 1900s |
| Celestial (modern) | 1-2 nm | Current |
| GPS (standard) | 8-15 m | 1995-present |
| DGPS/corrections | ~1 m | 2000s |
| Autonomous needs | 10 cm | 2020s |
| Algorithm target | ~0.1' (185 m LOP) | Research |

#### Content for Paper Sections:

**For Introduction:**
- Historical trend positions celestial navigation as foundation of modern navigation
- "10× improvement every 30 years" provides context for algorithm development
- GPS vulnerabilities motivate maintaining backup capabilities

**For Literature Review:**
- Navigation technology evolution timeline
- Autonomous vehicle requirements (10 cm 95%)
- GNSS vulnerability documentation

**For Discussion:**
- Celestial navigation's role in technology lineage
- Algorithm accuracy (0.1') vs autonomous needs (10 cm)
- Complementary role: celestial for GPS backup, not replacement

#### Key Quotes:
- "At the turn of the 20th century, celestial navigation was the state-of-the-art in location technology"
- "With every new order of magnitude in position accuracy, a new investment in infrastructure was required"
- "Decimeters are needed for humans and robotic systems to coexist and to share the same physical spaces"
- "The challenge facing auto makers is meeting the required level of reliability at 99.999999%"

#### Relationship to Other Sources:

| This Paper | Related Source | Connection |
|------------|----------------|------------|
| Historical trend | Dunlap 1971 | Navigation evolution |
| GPS vulnerabilities | Zalewski 2022 | GNSS standards/integrity |
| GPS vulnerabilities | Dachev 2017 | Disruption taxonomy |
| Celestial as backup | Lušić 2024 | Regulatory framework |
| Autonomous needs | Algorithm validation context |

---

### Article #47
**Citation:** Critchley-Marrows, J. J. R., Wu, X., & Cairns, I. H. (2023). An architecture for a visual-based PNT alternative. *Acta Astronautica*, 210, 601–609. doi:10.1016/j.actaastro.2023.05.022

**Relevance:** Modern automated celestial navigation using star trackers and horizon sensing; documents GNSS vulnerabilities and PNT requirements for maritime, aviation, and lunar applications; achieves <100m performance for LEO satellites

#### Core Technology - CROSS Star Tracker:

| Parameter | Value |
|-----------|-------|
| Cross-Axis Accuracy | 10″ |
| Boresight Accuracy | 25″ |
| Field of View | 20° (star sense), 40° (horizon sense) |
| Star Brightness Cut-off | 5.0 mag |
| Sky Availability | 99% |
| Sun Exclusion Angle | 50° |
| Update Rate | 2 Hz |
| Power | 2.5W nominal, 4.5W peak |

**Architecture:** Dual-sensor optical assembly
- "Star sense": Attitude determination from star field observations
- "Horizon sense": Position reference from planetary limb detection

#### PNT Performance Requirements by Domain:

**Maritime (from IALA R-129):**
| Phase | Accuracy | Time-to-Alert |
|-------|----------|---------------|
| Oceanic | 1 km | 1 min |
| Coastal | 100 m | 30 s |
| Port Approach | 10 m | 10 s |
| Port | 1 m | 10 s |
| Inland Waterways | 10 m | 10 s |

**Aviation (from ICAO GNSS Manual):**
| Phase | Accuracy | Time-to-Alert |
|-------|----------|---------------|
| Oceanic | 7.4 km | 5 min |
| Continental | 3.7 km | 5 min |
| Terminal | 1.85 km | 15 s |
| Approach | 556 m | 10 s |

**Lunar (from NASA LunaNet):**
| Phase | Accuracy |
|-------|----------|
| Lunar Transfer | 13 km |
| Orbit | 1 km |
| Initial Descent | 300 m |
| Final Descent | 100 m |
| Surface Position | 10 m |

#### Visual-Based Navigation Performance:

**Simulation Results (LEO orbit):**
| Altitude | Cross-Track | Along-Track | Radial |
|----------|-------------|-------------|--------|
| 300 km | 0.62 m (σ=39.1 m) | 0.14 m (σ=65.9 m) | -12.6 m (σ=25.0 m) |
| 600 km | -19.18 m (σ=41.3 m) | -32.54 m (σ=43.5 m) | -16.2 m (σ=31.0 m) |
| 1200 km | 2.83 m (σ=31.0 m) | 68.4 m (σ=99.9 m) | -51.4 m (σ=61.2 m) |

**Key finding:** Steady state achieved within one orbit (~90 min), margin below 100m

#### Atmospheric Offset Problem:

- Earth atmosphere creates ~25 km offset in observed horizon vs true ellipsoid
- Problem not well-documented in literature (Christian 2016 cited)
- Solution: Estimate ρ (atmospheric offset) as state parameter in EKF
- Offset reduced to ~10 m average within 10 minutes of filter convergence

**Modified ellipsoid matrix:**
```
Ã = diag{1/(a+ρ)², 1/(b+ρ)², 1/(c+ρ)²}
```

#### Algorithm Implementation:

**Position Localization:**
1. Star tracker provides inertial attitude (TPB)
2. Horizon edge identified via strip search + Sobel gradient + Zernike moments
3. Hyperbola fitted to horizon points (Fitzgibbon conic fitting)
4. RANSAC outlier rejection
5. Christian-Robinson method for position (Cholesky parameterization)
6. Extended Kalman Filter for state estimation

**EKF State Vector:**
```
x = {r, v, ρ}ᵀ  (position, velocity, atmospheric offset)
```

**Horizon Edge Precision:** ~2 px RMS (from ISS imagery analysis)

#### GNSS Vulnerabilities Context:

**Jamming/Spoofing concerns:**
- "These mechanisms [GNSS] are also easily jammed and spoofed from malicious sources"
- References C4ADS study on GPS spoofing in Russia/Syria
- Black Sea cargo ships reporting position outside Moscow
- Shanghai port spoofing incidents
- Red Sea disruptions

**Visual-based advantages:**
- "Disrupting a camera would require a very close or extremely powerful light emitting device"
- Natural references not dependent on man-made infrastructure
- Harder to manipulate or deceive

#### LEO Ranging Service Application:

**User Range Error (σURE):**
```
σ²URE = σ²eph + σ²clock + σ²atm + σ²rcv + σ²sat + σ²other
```

- Most terms on order 0.1–1 m for LEO systems
- Visual-based ephemeris error: ~50 m
- GDOP for LEO constellations: ~1.0
- **User position accuracy: ~50 m** (meets 100 m backup requirement)

#### Technology Positioning:

**Evolution of autonomous celestial navigation:**
| Era | Technology | Accuracy |
|-----|------------|----------|
| Traditional | Sextant + manual computation | 1-2 nm |
| Computer-based | Sextant + software | sub-nm possible |
| Star tracker | CCD/CMOS + autonomous ID | sub-arcsecond attitude |
| Integrated visual PNT | Star + horizon sensing | <100 m position |

**Key quote:** "The current state-of-art visual navigation systems are two orders of magnitude worse than their RF-based equivalents."

#### Comparison: Sextant vs Star Tracker:

| Aspect | Traditional Sextant | CROSS Star Tracker |
|--------|--------------------|--------------------|  
| Observation | Manual altitude measurement | Automatic star centroid |
| Accuracy (attitude) | 0.33-0.70' (Gordon 1964) | 10-25" |
| Position | 1-2 nm (multi-body fix) | <100 m (with horizon) |
| Horizon | Visual (dip limited) | Limb detection algorithm |
| Atmospheric correction | Tables | EKF state estimation |
| Update rate | Minutes | 2 Hz |

#### Content for Paper Sections:

**For Introduction:**
- Modern evolution of celestial navigation to automated systems
- Visual-based PNT as GNSS alternative architecture
- Performance requirements for backup navigation systems

**For Literature Review:**
- Star tracker technology in celestial navigation context
- Maritime/aviation PNT requirements (IALA, ICAO standards)
- Horizon-based position determination algorithms

**For Discussion:**
- Traditional sextant vs automated star tracker comparison
- Algorithm accuracy (0.1') fits between sextant (0.33') and star tracker (0.003')
- Research positions in celestial navigation technology evolution

#### Relationship to Other Sources:

| This Paper | Related Source | Connection |
|------------|----------------|------------|
| Star tracker tech | Dachev 2017 | Method classification (automated category) |
| GNSS vulnerabilities | Reid 2020 | Complementary vulnerability documentation |
| GNSS vulnerabilities | Zalewski 2022 | IMO standards context |
| Horizon sensing | Critchley-Marrows 2023 (Sensors) | Same authors, maritime sextant focus |
| Maritime requirements | Lušić 2024 | STCW/IMO regulatory framework |
| Position algorithms | Christian (cited) | Horizon-based optical navigation |

---

## Pending Content Needs

### For Literature Review:
- [x] GHA, LHA, declination computation algorithms - Feldman et al. 1972
- [x] Altitude correction procedures - Hohenkerk et al. 2012, Pepperday 1994
- [x] Sextant error sources and magnitudes - Ross 1994

### For Technical Implementation:
- [x] Python libraries for ephemeris (Skyfield) - Barazzetti 2025
- [x] Nautical Almanac data sources - Hohenkerk et al. 2012
- [x] Validation datasets - Park et al. 2021 (DE440/DE441 accuracy)
- [x] Ephemeris uncertainty analysis - Standish & Fienga 2002

### For Results/Discussion:
- [x] Modern accuracy standards - Kotlarić 1976, Hohenkerk 2012
- [x] Comparison benchmarks - Kotlarić 1976, Tsai 2022, Feldman 1972
- [x] Reference software implementation - NavPac (Hohenkerk et al. 2012)
- [x] Ephemeris accuracy baseline - Park et al. 2021, Standish & Fienga 2002
- [x] Error budget and propagation - Ross 1994

---

*Last updated: Phase 2 in progress - 47 peer-reviewed + 3 non-peer-reviewed sources accepted*
