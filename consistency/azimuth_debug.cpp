/*
 * Azimuth Convention Debug Tool
 * 
 * This program tests different azimuth conventions and coordinate transformations
 * to identify the correct one for your mount system.
 * 
 * Compile with:
 * g++ -std=c++11 -o azimuth_debug azimuth_debug.cpp -lm
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

struct AzimuthConvention {
    std::string name;
    std::function<double(double)> convert;
    std::string description;
};

class AzimuthDebugger {
private:
    std::vector<AzimuthConvention> conventions;
    
    // Utility functions
    double normalizeAngle(double angle) const;
    double greenwichMeanSiderealTime(double jd) const;
    double earthRotationAngle(double jd) const;
    
    // Standard coordinate transformations
    void radecToAltAz(double ra, double dec, double lat, double lst,
                     double& alt, double& az) const;
    void altazToRaDec(double alt, double az, double lat, double lst,
                     double& ra, double& dec, double& ha) const;
    
public:
    AzimuthDebugger();
    void testAllConventions(double your_alt, double your_az, double lat, double lst,
                           double expected_ra, double expected_dec);
    void testRoundTrip(double ra, double dec, double lat, double lst);
};

AzimuthDebugger::AzimuthDebugger() {
    // Define different azimuth conventions to test
    conventions = {
        {"Original", [](double az) { return az; }, "Your current azimuth (no change)"},
        {"North+180", [](double az) { return fmod(az + 180.0, 360.0); }, "Add 180° (South=0° instead of North=0°)"},
        {"Complement", [](double az) { return 360.0 - az; }, "360° - azimuth (reverse direction)"},
        {"North-180", [](double az) { return fmod(az - 180.0 + 360.0, 360.0); }, "Subtract 180°"},
        {"East=0", [](double az) { return fmod(az + 90.0, 360.0); }, "Add 90° (East=0° instead of North=0°)"},
        {"West=0", [](double az) { return fmod(az - 90.0 + 360.0, 360.0); }, "Subtract 90° (West=0°)"},
        {"Flip+180", [](double az) { return fmod(360.0 - az + 180.0, 360.0); }, "360° - azimuth + 180°"},
        {"Astronomical", [](double az) { return fmod(az + 180.0, 360.0); }, "Astronomical convention (South=0°, East=90°)"}
    };
}

double AzimuthDebugger::normalizeAngle(double angle) const {
    angle = fmod(angle, 360.0);
    if (angle < 0) angle += 360.0;
    return angle;
}

double AzimuthDebugger::earthRotationAngle(double jd) const {
    const double t = jd - 2451545.0;
    const double f = fmod(jd, 1.0);
    double theta = 2.0 * M_PI * (f + 0.7790572732640 + 0.00273781191135448 * t);
    return fmod(theta, 2.0 * M_PI);
}

double AzimuthDebugger::greenwichMeanSiderealTime(double jd) const {
    const double t = (jd - 2451545.0) / 36525.0;
    double gmst = earthRotationAngle(jd) + 
                  (0.014506 + 4612.156534*t + 1.3915817*t*t - 0.00000044*t*t*t - 
                   0.000029956*t*t*t*t - 0.0000000368*t*t*t*t*t) / 60.0 / 60.0 * M_PI / 180.0;
    return fmod(gmst, 2.0 * M_PI);
}

void AzimuthDebugger::radecToAltAz(double ra, double dec, double lat, double lst,
                                  double& alt, double& az) const {
    // Standard conversion using Explanatory Supplement method
    const double ra_rad = ra * M_PI / 180.0;
    const double dec_rad = dec * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    const double lst_rad = lst * M_PI / 12.0;
    
    double H = lst_rad - ra_rad;
    if (H < 0) H += 2.0 * M_PI;
    if (H > M_PI) H = H - 2.0 * M_PI;
    
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
    az = normalizeAngle(az_rad * 180.0 / M_PI);
}

void AzimuthDebugger::altazToRaDec(double alt, double az, double lat, double lst,
                                  double& ra, double& dec, double& ha) const {
    // Standard Alt/Az to RA/Dec conversion
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    // Calculate declination
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    // Calculate hour angle
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    // Determine sign of hour angle
    if (sin(az_rad) > 0) {
        ha_rad = -ha_rad;  // Standard: sin(az) > 0 means setting (negative HA)
    }
    
    ha = ha_rad * 12.0 / M_PI;
    
    // Calculate RA
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void AzimuthDebugger::testAllConventions(double your_alt, double your_az, double lat, double lst,
                                        double expected_ra, double expected_dec) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Testing azimuth conventions for Alt/Az → RA/Dec conversion:" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "Input: Alt=" << your_alt << "° Az=" << your_az << "° Lat=" << lat << "° LST=" << lst << "h" << std::endl;
    std::cout << "Expected: RA=" << expected_ra << "° Dec=" << expected_dec << "°" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    double best_error = 1000.0;
    std::string best_convention;
    
    for (const auto& conv : conventions) {
        double converted_az = conv.convert(your_az);
        double ra, dec, ha;
        
        altazToRaDec(your_alt, converted_az, lat, lst, ra, dec, ha);
        
        double ra_error = fabs(ra - expected_ra);
        double dec_error = fabs(dec - expected_dec);
        
        // Handle RA wraparound
        if (ra_error > 180.0) ra_error = 360.0 - ra_error;
        
        double total_error = sqrt(ra_error * ra_error + dec_error * dec_error);
        
        std::cout << std::setw(12) << conv.name << ": ";
        std::cout << "Az=" << std::setw(8) << converted_az << "° → ";
        std::cout << "RA=" << std::setw(9) << ra << "° ";
        std::cout << "Dec=" << std::setw(9) << dec << "° ";
        std::cout << "Error=" << std::setw(8) << total_error << "° ";
        
        if (total_error < best_error) {
            best_error = total_error;
            best_convention = conv.name;
            std::cout << " ← BEST";
        }
        
        std::cout << std::endl;
        std::cout << std::setw(15) << " " << conv.description << std::endl;
    }
    
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "RESULT: Best convention is '" << best_convention << "' with error = " << best_error << "°" << std::endl;
    
    // Test the best convention with round-trip
    std::cout << "\nTesting round-trip with best convention (" << best_convention << "):" << std::endl;
    for (const auto& conv : conventions) {
        if (conv.name == best_convention) {
            double converted_az = conv.convert(your_az);
            double ra, dec, ha;
            altazToRaDec(your_alt, converted_az, lat, lst, ra, dec, ha);
            
            // Now convert back
            double back_alt, back_az;
            radecToAltAz(ra, dec, lat, lst, back_alt, back_az);
            
            std::cout << "Forward:  Alt=" << your_alt << "° Az=" << your_az << "° (original)" << std::endl;
            std::cout << "         Alt=" << your_alt << "° Az=" << converted_az << "° (converted)" << std::endl;
            std::cout << "         → RA=" << ra << "° Dec=" << dec << "° HA=" << ha << "h" << std::endl;
            std::cout << "Reverse: RA=" << ra << "° Dec=" << dec << "°" << std::endl;
            std::cout << "         → Alt=" << back_alt << "° Az=" << back_az << "°" << std::endl;
            
            double alt_error = fabs(your_alt - back_alt);
            double az_error = fabs(converted_az - back_az);
            if (az_error > 180.0) az_error = 360.0 - az_error;
            
            std::cout << "Round-trip errors: Alt=" << alt_error << "° Az=" << az_error << "°" << std::endl;
            break;
        }
    }
}

void AzimuthDebugger::testRoundTrip(double ra, double dec, double lat, double lst) {
    std::cout << "\nTesting RA/Dec → Alt/Az → RA/Dec round-trip:" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Input: RA=" << ra << "° Dec=" << dec << "° Lat=" << lat << "° LST=" << lst << "h" << std::endl;
    
    // Forward conversion
    double alt, az;
    radecToAltAz(ra, dec, lat, lst, alt, az);
    std::cout << "→ Alt=" << alt << "° Az=" << az << "°" << std::endl;
    
    // Test all azimuth conventions for reverse conversion
    std::cout << "\nTesting reverse conversions with different azimuth conventions:" << std::endl;
    
    double best_error = 1000.0;
    std::string best_convention;
    
    for (const auto& conv : conventions) {
        double converted_az = conv.convert(az);
        double back_ra, back_dec, back_ha;
        
        altazToRaDec(alt, converted_az, lat, lst, back_ra, back_dec, back_ha);
        
        double ra_error = fabs(ra - back_ra);
        double dec_error = fabs(dec - back_dec);
        
        // Handle RA wraparound
        if (ra_error > 180.0) ra_error = 360.0 - ra_error;
        
        double total_error = sqrt(ra_error * ra_error + dec_error * dec_error);
        
        std::cout << std::setw(12) << conv.name << ": ";
        std::cout << "Az=" << std::setw(8) << converted_az << "° → ";
        std::cout << "RA=" << std::setw(9) << back_ra << "° ";
        std::cout << "Dec=" << std::setw(9) << back_dec << "° ";
        std::cout << "Error=" << std::setw(8) << total_error << "°";
        
        if (total_error < best_error) {
            best_error = total_error;
            best_convention = conv.name;
            std::cout << " ← BEST";
        }
        
        std::cout << std::endl;
    }
    
    std::cout << "\nBest round-trip convention: " << best_convention << " (error = " << best_error << "°)" << std::endl;
}

int main() {
    AzimuthDebugger debugger;
    
    std::cout << "Azimuth Convention Debug Tool" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    // Test with your specific FITS data
    std::cout << "\n1. Testing your FITS file data:" << std::endl;
    debugger.testAllConventions(
        42.041003,   // STELLALT
        286.852567,  // STELLAZ (already corrected to South-based)
        51.507400,   // OBSLAT
        5.47957059,  // LST
        10.912770,   // CRVAL1 (expected RA)
        41.212212    // CRVAL2 (expected Dec)
    );
    
    // Test round-trip with known good coordinates
    std::cout << "\n\n2. Testing round-trip with your WCS coordinates:" << std::endl;
    debugger.testRoundTrip(
        10.912770,   // CRVAL1
        41.212212,   // CRVAL2
        51.507400,   // OBSLAT
        5.47957059   // LST
    );
    
    // Test with a few reference cases from the HTML test data
    std::cout << "\n\n3. Testing with reference cases:" << std::endl;
    
    // Case 1: High latitude, different quadrant
    std::cout << "\nReference case 1 (High latitude):" << std::endl;
    debugger.testAllConventions(19.4223, 145.5122, 81.9121742352522, 7.9250404113, 152.0646154, 12.7006088);
    
    // Case 2: Southern hemisphere
    std::cout << "\nReference case 2 (Southern hemisphere):" << std::endl;
    debugger.testAllConventions(-22.8271, 124.4669, -45.0396581387729, 4.5307070225, 198.2076359, -5.3974508);
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "Analysis complete. Check which convention gives the smallest errors." << std::endl;
    std::cout << "The correct azimuth convention should give errors < 1° for all test cases." << std::endl;
    
    return 0;
}
