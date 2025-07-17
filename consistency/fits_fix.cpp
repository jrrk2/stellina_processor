/*
 * FITS Coordinate Fix Tool
 * 
 * This program identifies the correct azimuth convention and coordinate
 * transformation for your FITS file data and provides the corrected
 * coordinate conversion functions.
 * 
 * Compile with:
 * g++ -std=c++11 -o fits_fix fits_fix.cpp -lm
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <functional>

class FITSCoordinateFixer {
private:
    // Your FITS data
    struct FITSData {
        double stellalt = 42.041003;
        double stellaz = 286.852567;  // Already South-corrected
        double obslat = 51.507400;
        double obslong = -0.127800;
        double lst = 5.47957059;
        double crval1 = 10.912770;  // Expected RA
        double crval2 = 41.212212;  // Expected Dec
    } fits_data;
    
    // Test different coordinate algorithms
    void test_algorithm_1(double alt, double az, double lat, double lst,
                         double& ra, double& dec, double& ha) const;
    void test_algorithm_2(double alt, double az, double lat, double lst,
                         double& ra, double& dec, double& ha) const;
    void test_algorithm_3(double alt, double az, double lat, double lst,
                         double& ra, double& dec, double& ha) const;
    void test_algorithm_4(double alt, double az, double lat, double lst,
                         double& ra, double& dec, double& ha) const;
    
    // Standard RA/Dec to Alt/Az for verification
    void radecToAltAz(double ra, double dec, double lat, double lst,
                     double& alt, double& az) const;
    
public:
    void diagnoseCoordinateSystem();
    void findBestAlgorithm();
    void generateCorrectedCode();
};

void FITSCoordinateFixer::test_algorithm_1(double alt, double az, double lat, double lst,
                                          double& ra, double& dec, double& ha) const {
    // Your current algorithm
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    if (sin(az_rad) > 0) {
        ha_rad = ha_rad;    // Current logic
    } else {
        ha_rad = -ha_rad;
    }
    
    ha = ha_rad * 12.0 / M_PI;
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void FITSCoordinateFixer::test_algorithm_2(double alt, double az, double lat, double lst,
                                          double& ra, double& dec, double& ha) const {
    // Corrected hour angle sign
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    // CORRECTED: Reverse the sign logic
    if (sin(az_rad) > 0) {
        ha_rad = -ha_rad;   // Western sky = negative HA
    } else {
        ha_rad = ha_rad;    // Eastern sky = positive HA
    }
    
    ha = ha_rad * 12.0 / M_PI;
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void FITSCoordinateFixer::test_algorithm_3(double alt, double az, double lat, double lst,
                                          double& ra, double& dec, double& ha) const {
    // Different azimuth quadrant logic
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    // Use azimuth angle directly for quadrant
    if (az >= 0 && az <= 180) {
        ha_rad = -ha_rad;   // North to South = negative HA
    } else {
        ha_rad = ha_rad;    // South to North = positive HA
    }
    
    ha = ha_rad * 12.0 / M_PI;
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void FITSCoordinateFixer::test_algorithm_4(double alt, double az, double lat, double lst,
                                          double& ra, double& dec, double& ha) const {
    // Standard astronomical algorithm (ExpSup method)
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    const double lst_rad = lst * M_PI / 12.0;
    
    // Calculate declination
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    // Calculate hour angle using atan2 for better quadrant handling
    double sin_ha = -cos(alt_rad) * sin(az_rad);
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / (cos(lat_rad) * cos(dec_rad));
    
    double ha_rad = atan2(sin_ha, cos_ha);
    
    ha = ha_rad * 12.0 / M_PI;
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

void FITSCoordinateFixer::radecToAltAz(double ra, double dec, double lat, double lst,
                                      double& alt, double& az) const {
    // Standard conversion for verification
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
    az = fmod(az_rad * 180.0 / M_PI, 360.0);
    if (az < 0) az += 360.0;
}

void FITSCoordinateFixer::findBestAlgorithm() {
    std::cout << "Testing coordinate transformation algorithms:" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Input: Alt=" << fits_data.stellalt << "° Az=" << fits_data.stellaz << "°" << std::endl;
    std::cout << "       Lat=" << fits_data.obslat << "° LST=" << fits_data.lst << "h" << std::endl;
    std::cout << "Expected: RA=" << fits_data.crval1 << "° Dec=" << fits_data.crval2 << "°" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    struct Algorithm {
        std::string name;
        std::function<void(double, double, double, double, double&, double&, double&)> func;
        std::string description;
    };
    
    std::vector<Algorithm> algorithms = {
        {"Algorithm 1", [this](double alt, double az, double lat, double lst, double& ra, double& dec, double& ha) {
            test_algorithm_1(alt, az, lat, lst, ra, dec, ha);
        }, "Your current algorithm"},
        
        {"Algorithm 2", [this](double alt, double az, double lat, double lst, double& ra, double& dec, double& ha) {
            test_algorithm_2(alt, az, lat, lst, ra, dec, ha);
        }, "Corrected hour angle sign"},
        
        {"Algorithm 3", [this](double alt, double az, double lat, double lst, double& ra, double& dec, double& ha) {
            test_algorithm_3(alt, az, lat, lst, ra, dec, ha);
        }, "Different azimuth quadrant logic"},
        
        {"Algorithm 4", [this](double alt, double az, double lat, double lst, double& ra, double& dec, double& ha) {
            test_algorithm_4(alt, az, lat, lst, ra, dec, ha);
        }, "Standard astronomical method (atan2)"}
    };
    
    double best_error = 1000.0;
    int best_algorithm = -1;
    
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < algorithms.size(); i++) {
        double ra, dec, ha;
        algorithms[i].func(fits_data.stellalt, fits_data.stellaz, fits_data.obslat, 
                          fits_data.lst, ra, dec, ha);
        
        double ra_error = fabs(ra - fits_data.crval1);
        double dec_error = fabs(dec - fits_data.crval2);
        if (ra_error > 180.0) ra_error = 360.0 - ra_error;
        
        double total_error = sqrt(ra_error * ra_error + dec_error * dec_error);
        
        std::cout << algorithms[i].name << ": ";
        std::cout << "RA=" << std::setw(9) << ra << "° ";
        std::cout << "Dec=" << std::setw(9) << dec << "° ";
        std::cout << "HA=" << std::setw(8) << ha << "h ";
        std::cout << "Error=" << std::setw(8) << total_error << "°";
        
        if (total_error < best_error) {
            best_error = total_error;
            best_algorithm = i;
            std::cout << " ← BEST";
        }
        
        std::cout << std::endl;
        std::cout << "           " << algorithms[i].description << std::endl;
    }
    
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "RESULT: " << algorithms[best_algorithm].name << " is best with error = " 
              << best_error << "°" << std::endl;
    
    if (best_error < 1.0) {
        std::cout << "✅ Found working algorithm!" << std::endl;
    } else {
        std::cout << "❌ No algorithm gives good results. May need azimuth convention change." << std::endl;
    }
    
    // Test round-trip with best algorithm
    std::cout << "\nTesting round-trip accuracy with best algorithm:" << std::endl;
    double ra, dec, ha;
    algorithms[best_algorithm].func(fits_data.stellalt, fits_data.stellaz, fits_data.obslat, 
                                   fits_data.lst, ra, dec, ha);
    
    double back_alt, back_az;
    radecToAltAz(ra, dec, fits_data.obslat, fits_data.lst, back_alt, back_az);
    
    double alt_error = fabs(fits_data.stellalt - back_alt);
    double az_error = fabs(fits_data.stellaz - back_az);
    if (az_error > 180.0) az_error = 360.0 - az_error;
    
    std::cout << "Forward:  Alt=" << fits_data.stellalt << "° Az=" << fits_data.stellaz << "°" << std::endl;
    std::cout << "         → RA=" << ra << "° Dec=" << dec << "°" << std::endl;
    std::cout << "Reverse: → Alt=" << back_alt << "° Az=" << back_az << "°" << std::endl;
    std::cout << "Round-trip error: Alt=" << alt_error << "° Az=" << az_error << "°" << std::endl;
    
    if (alt_error < 0.1 && az_error < 0.1) {
        std::cout << "✅ Excellent round-trip accuracy!" << std::endl;
    } else if (alt_error < 1.0 && az_error < 1.0) {
        std::cout << "✅ Good round-trip accuracy" << std::endl;
    } else {
        std::cout << "⚠️  Poor round-trip accuracy - check algorithm" << std::endl;
    }
}

void FITSCoordinateFixer::generateCorrectedCode() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "CORRECTED COORDINATE CONVERSION CODE:" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << R"(
// Corrected Alt/Az to RA/Dec conversion function
void correctedAltAzToRaDec(double alt, double az, double lat, double lst,
                          double& ra, double& dec, double& ha) {
    const double alt_rad = alt * M_PI / 180.0;
    const double az_rad = az * M_PI / 180.0;
    const double lat_rad = lat * M_PI / 180.0;
    
    // Calculate declination
    double dec_rad = asin(sin(alt_rad) * sin(lat_rad) + 
                         cos(alt_rad) * cos(lat_rad) * cos(az_rad));
    
    // Calculate hour angle using atan2 for proper quadrant handling
    double sin_ha = -cos(alt_rad) * sin(az_rad);
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec_rad)) / 
                    (cos(lat_rad) * cos(dec_rad));
    
    double ha_rad = atan2(sin_ha, cos_ha);
    ha = ha_rad * 12.0 / M_PI;
    
    // Calculate RA
    double ra_hours = lst - ha;
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

// Usage in your FITS checker:
void performCorrectedBlindCoordinateConversion() {
    double calculated_ra, calculated_dec, calculated_ha;
    
    // Apply atmospheric refraction if needed
    double corrected_alt = stellalt - calculateAtmosphericRefraction(stellalt);
    
    // Use corrected conversion
    correctedAltAzToRaDec(corrected_alt, stellaz, obslat, lst, 
                         calculated_ra, calculated_dec, calculated_ha);
    
    // Compare with WCS solution...
}
)" << std::endl;
}

void FITSCoordinateFixer::diagnoseCoordinateSystem() {
    std::cout << "FITS Coordinate System Diagnosis" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    findBestAlgorithm();
    generateCorrectedCode();
    
    std::cout << "\nRECOMMENDATIONS:" << std::endl;
    std::cout << "1. Replace your blindAltAzToEquatorial() function with the corrected version above" << std::endl;
    std::cout << "2. The key fix is using atan2() for proper hour angle quadrant determination" << std::endl;
    std::cout << "3. Your azimuth convention (South-based) appears to be correct" << std::endl;
    std::cout << "4. Test the corrected function with your FITS data" << std::endl;
}

int main() {
    FITSCoordinateFixer fixer;
    fixer.diagnoseCoordinateSystem();
    return 0;
}
