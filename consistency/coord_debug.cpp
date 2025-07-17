/*
 * Coordinate Transformation Debug Tool
 * 
 * This program tests various coordinate transformation algorithms against
 * known good test cases from NASA JPL Horizons to identify the correct
 * implementation for your mount system.
 * 
 * Compile with:
 * g++ -std=c++11 -o coord_debug coord_debug.cpp -lm
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

struct TestCase {
    double lat, lon;           // Observer location (degrees)
    double jd;                 // Julian Date
    double ra, dec;            // RA/Dec (degrees)
    double expected_az, expected_alt;  // Expected Alt/Az (degrees)
    double expected_lst;       // Expected LST (hours)
    double expected_ha;        // Expected Hour Angle (hours)
};

class CoordinateDebugger {
private:
    std::vector<TestCase> test_cases;
    
    // Different coordinate transformation algorithms to test
    void meeus_radec_to_altaz(double ra, double dec, double lat, double lon, double jd,
                             double& alt, double& az, double& lst, double& ha) const;
    void expsup_radec_to_altaz(double ra, double dec, double lat, double lon, double jd,
                              double& alt, double& az, double& lst, double& ha) const;
    void green_radec_to_altaz(double ra, double dec, double lat, double lon, double jd,
                             double& alt, double& az, double& lst, double& ha) const;
    void taff_radec_to_altaz(double ra, double dec, double lat, double lon, double jd,
                            double& alt, double& az, double& lst, double& ha) const;
    
    // Different Alt/Az to RA/Dec algorithms
    void altaz_to_radec_simple(double alt, double az, double lat, double lst,
                              double& ra, double& dec, double& ha) const;
    void altaz_to_radec_correct(double alt, double az, double lat, double lst,
                               double& ra, double& dec, double& ha) const;
    
    // Utility functions
    double greenwichMeanSiderealTime(double jd) const;
    double earthRotationAngle(double jd) const;
    double normalizeAngle(double angle) const;
    double normalizeLST(double lst) const;
    double normalizeHA(double ha) const;
    
public:
    void loadTestCases();
    void runAllTests();
    void testSpecificCase(double alt, double az, double lat, double lst, 
                         double expected_ra, double expected_dec);
};

double CoordinateDebugger::earthRotationAngle(double jd) const {
    // IERS Technical Note No. 32
    const double t = jd - 2451545.0;
    const double f = fmod(jd, 1.0);
    
    double theta = 2.0 * M_PI * (f + 0.7790572732640 + 0.00273781191135448 * t);
    return normalizeAngle(theta);
}

double CoordinateDebugger::greenwichMeanSiderealTime(double jd) const {
    // "Expressions for IAU 2000 precession quantities"
    const double t = (jd - 2451545.0) / 36525.0;
    
    double gmst = earthRotationAngle(jd) + 
                  (0.014506 + 4612.156534*t + 1.3915817*t*t - 0.00000044*t*t*t - 
                   0.000029956*t*t*t*t - 0.0000000368*t*t*t*t*t) / 60.0 / 60.0 * M_PI / 180.0;
    
    return normalizeAngle(gmst);
}

double CoordinateDebugger::normalizeAngle(double angle) const {
    angle = fmod(angle, 2.0 * M_PI);
    if (angle < 0) angle += 2.0 * M_PI;
    return angle;
}

double CoordinateDebugger::normalizeLST(double lst) const {
    while (lst < 0) lst += 24.0;
    while (lst >= 24.0) lst -= 24.0;
    return lst;
}

double CoordinateDebugger::normalizeHA(double ha) const {
    while (ha > 12.0) ha -= 24.0;
    while (ha < -12.0) ha += 24.0;
    return ha;
}

void CoordinateDebugger::meeus_radec_to_altaz(double ra, double dec, double lat, double lon, double jd,
                                             double& alt, double& az, double& lst, double& ha) const {
    // Meeus "Astronomical Algorithms" 13.5 and 13.6
    const double gmst = greenwichMeanSiderealTime(jd);
    lst = normalizeLST((gmst * 12.0 / M_PI) + (lon / 15.0));
    
    // Convert to radians
    const double ra_rad = ra * M_PI / 180.0;
    const double dec_rad = dec * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    const double lst_rad = lst * M_PI / 12.0;
    
    // Calculate hour angle
    double H = lst_rad - ra_rad;
    if (H < 0) H += 2.0 * M_PI;
    if (H > M_PI) H = H - 2.0 * M_PI;
    
    ha = normalizeHA(H * 12.0 / M_PI);
    
    // Calculate altitude and azimuth
    double az_rad = atan2(sin(H), cos(H) * sin(lat_rad) - tan(dec_rad) * cos(lat_rad));
    double alt_rad = asin(sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(H));
    
    az_rad -= M_PI;  // Meeus adjustment
    
    alt = alt_rad * 180.0 / M_PI;
    az = normalizeAngle(az_rad) * 180.0 / M_PI;
}

void CoordinateDebugger::expsup_radec_to_altaz(double ra, double dec, double lat, double lon, double jd,
                                              double& alt, double& az, double& lst, double& ha) const {
    // Explanatory Supplement equation 7.16
    const double gmst = greenwichMeanSiderealTime(jd);
    lst = normalizeLST((gmst * 12.0 / M_PI) + (lon / 15.0));
    
    const double ra_rad = ra * M_PI / 180.0;
    const double dec_rad = dec * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    const double lst_rad = lst * M_PI / 12.0;
    
    double H = lst_rad - ra_rad;
    if (H < 0) H += 2.0 * M_PI;
    if (H > M_PI) H = H - 2.0 * M_PI;
    
    ha = normalizeHA(H * 12.0 / M_PI);
    
    double alt_rad = asin(sin(dec_rad) * sin(lat_rad) + cos(dec_rad) * cos(H) * cos(lat_rad));
    
    double cosaz = (sin(dec_rad) * cos(lat_rad) - cos(dec_rad) * cos(H) * sin(lat_rad)) / cos(alt_rad);
    double sinaz = (-cos(dec_rad) * sin(H)) / cos(alt_rad);
    
    double az_rad;
    if (cosaz < 0 && sinaz < 0) {
        az_rad = 2.0 * M_PI - acos(cosaz);
    } else if (cosaz < 0 && sinaz >= 0) {
        az_rad = acos(cosaz);
    } else {
        az_rad = asin(sinaz);
    }
    
    alt = alt_rad * 180.0 / M_PI;
    az = normalizeAngle(az_rad) * 180.0 / M_PI;
}

void CoordinateDebugger::altaz_to_radec_simple(double alt, double az, double lat, double lst,
                                              double& ra, double& dec, double& ha) const {
    // Your current algorithm - let's test this
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    // Calculate declination
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    // Calculate hour angle
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    // Determine sign based on azimuth
    if (sin(az_rad) > 0) {  // Your current logic
        ha_rad = ha_rad;    // Positive HA
    } else {
        ha_rad = -ha_rad;   // Negative HA
    }
    
    ha = ha_rad * 12.0 / M_PI;
    
    // Calculate RA
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void CoordinateDebugger::altaz_to_radec_correct(double alt, double az, double lat, double lst,
                                               double& ra, double& dec, double& ha) const {
    // Corrected algorithm based on standard astronomy
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    // Calculate declination
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    // Calculate hour angle using correct quadrant determination
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    // Correct quadrant determination for hour angle
    if (sin(az_rad) > 0) {  // Azimuth in eastern half (object setting)
        ha_rad = -ha_rad;   // Negative HA (setting)
    } else {                // Azimuth in western half (object rising)
        ha_rad = ha_rad;    // Positive HA (rising)
    }
    
    ha = ha_rad * 12.0 / M_PI;
    
    // Calculate RA
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void CoordinateDebugger::loadTestCases() {
    // Load a subset of test cases from the HTML file
    test_cases = {
        // lat, lon, jd, ra, dec, az, alt, lst, ha
        {77.2069702520977, 118.639627806683, 2258936.199571759, 156.3041080, 9.9129372, 10.1558, -2.6785, 23.1069400896, -11.313333776},
        {-45.0396581387729, 52.6267552330626, 2279088.433946759, 198.2076359, -5.3974508, 124.4669, -22.8271, 4.5307070225, -8.683135370},
        {75.8072259286161, 308.733126034821, 2652064.354988426, 303.4519425, -21.2960037, 231.4286, -12.7578, 23.8918905112, 3.661761012},
        {81.9121742352522, 248.896246076003, 2435346.210578704, 152.0646154, 12.7006088, 145.5122, 19.4223, 7.9250404113, -2.212600613},
        {68.4688217628369, 60.2632226629549, 2326837.745578703, 116.6004466, 21.5438150, 149.0091, 40.6459, 6.1176176006, -1.655745503},
        {51.507400, -0.127800, 2460319.42603009, 10.912770, 41.212212, 107.241425, 41.754635, 5.47957209, 4.752054}  // Your case
    };
}

void CoordinateDebugger::testSpecificCase(double alt, double az, double lat, double lst, 
                                         double expected_ra, double expected_dec) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nTesting Alt/Az → RA/Dec conversion:" << std::endl;
    std::cout << "Input: Alt=" << alt << "° Az=" << az << "° Lat=" << lat << "° LST=" << lst << "h" << std::endl;
    std::cout << "Expected: RA=" << expected_ra << "° Dec=" << expected_dec << "°" << std::endl;
    
    double ra1, dec1, ha1, ra2, dec2, ha2;
    
    // Test your current algorithm
    altaz_to_radec_simple(alt, az, lat, lst, ra1, dec1, ha1);
    std::cout << "\nYour algorithm:" << std::endl;
    std::cout << "  RA=" << ra1 << "° Dec=" << dec1 << "° HA=" << ha1 << "h" << std::endl;
    std::cout << "  RA error=" << fabs(ra1 - expected_ra) << "° Dec error=" << fabs(dec1 - expected_dec) << "°" << std::endl;
    
    // Test corrected algorithm
    altaz_to_radec_correct(alt, az, lat, lst, ra2, dec2, ha2);
    std::cout << "\nCorrected algorithm:" << std::endl;
    std::cout << "  RA=" << ra2 << "° Dec=" << dec2 << "° HA=" << ha2 << "h" << std::endl;
    std::cout << "  RA error=" << fabs(ra2 - expected_ra) << "° Dec error=" << fabs(dec2 - expected_dec) << "°" << std::endl;
}

void CoordinateDebugger::runAllTests() {
    std::cout << "Testing coordinate transformation algorithms against NASA JPL Horizons data:" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    double total_error_meeus = 0, total_error_expsup = 0;
    int test_count = 0;
    
    for (const auto& tc : test_cases) {
        std::cout << "\nTest case " << (test_count + 1) << ":" << std::endl;
        std::cout << "Lat=" << tc.lat << "° Lon=" << tc.lon << "° JD=" << std::fixed << std::setprecision(6) << tc.jd << std::endl;
        std::cout << "RA=" << tc.ra << "° Dec=" << tc.dec << "°" << std::endl;
        std::cout << "Expected: Alt=" << tc.expected_alt << "° Az=" << tc.expected_az << "°" << std::endl;
        
        double alt1, az1, lst1, ha1, alt2, az2, lst2, ha2;
        
        // Test Meeus algorithm
        meeus_radec_to_altaz(tc.ra, tc.dec, tc.lat, tc.lon, tc.jd, alt1, az1, lst1, ha1);
        double error1 = sqrt(pow(alt1 - tc.expected_alt, 2) + pow(az1 - tc.expected_az, 2));
        
        // Test Explanatory Supplement algorithm
        expsup_radec_to_altaz(tc.ra, tc.dec, tc.lat, tc.lon, tc.jd, alt2, az2, lst2, ha2);
        double error2 = sqrt(pow(alt2 - tc.expected_alt, 2) + pow(az2 - tc.expected_az, 2));
        
        std::cout << "Meeus:   Alt=" << alt1 << "° Az=" << az1 << "° Error=" << error1 << "°" << std::endl;
        std::cout << "ExpSup:  Alt=" << alt2 << "° Az=" << az2 << "° Error=" << error2 << "°" << std::endl;
        
        total_error_meeus += error1;
        total_error_expsup += error2;
        test_count++;
        
        // Test round-trip for the most accurate algorithm
        if (error2 < error1) {  // ExpSup is better
            std::cout << "\nTesting Alt/Az → RA/Dec round-trip (using ExpSup reference):" << std::endl;
            testSpecificCase(alt2, az2, tc.lat, lst2, tc.ra, tc.dec);
        }
    }
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "Summary:" << std::endl;
    std::cout << "Meeus average error: " << (total_error_meeus / test_count) << "°" << std::endl;
    std::cout << "ExpSup average error: " << (total_error_expsup / test_count) << "°" << std::endl;
    
    // Test your specific case
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "Testing your specific FITS case:" << std::endl;
    testSpecificCase(42.041003, 286.852567, 51.507400, 5.47957059, 10.912770, 41.212212);
}

int main() {
    CoordinateDebugger debugger;
    debugger.loadTestCases();
    debugger.runAllTests();
    
    return 0;
}
