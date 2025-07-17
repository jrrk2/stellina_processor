/*
 * FITS WCS Consistency Checker - C++ Version
 * 
 * This program reads a plate-solved FITS file and checks the consistency
 * of World Coordinate System (WCS) keywords, particularly focusing on
 * comparing Alt/Az mount coordinates with RA/Dec WCS solution.
 * 
 * Dependencies:
 * - cfitsio
 * - libnova
 * 
 * Compile with:
 * g++ -std=c++11 -o fits_wcs_checker fits_wcs_checker.cpp -lcfitsio -lnova -lm
 * 
 * Usage:
 * ./fits_wcs_checker <fits_file_path>
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <fitsio.h>
#include <libnova/libnova.h>

struct FITSKeywords {
    // Mount coordinates
    double stellalt = 0.0;     // STELLALT - altitude (degrees)
    double stellaz = 0.0;      // STELLAZ - azimuth (degrees)
    
    // Observer location
    double obslat = 0.0;       // OBSLAT - latitude (degrees)
    double obslong = 0.0;      // OBSLONG - longitude (degrees)
    
    // Time information
    double julian = 0.0;       // JULIAN - Julian date
    double lst = 0.0;          // LST - Local sidereal time (hours)
    double hourang = 0.0;      // HOURANG - Hour angle (hours)
    
    // WCS solution
    double crval1 = 0.0;       // CRVAL1 - RA (degrees)
    double crval2 = 0.0;       // CRVAL2 - Dec (degrees)
    
    // Image properties
    int naxis1 = 0;            // Image width
    int naxis2 = 0;            // Image height
    
    // Instrument info
    std::string instrume;      // INSTRUME
    std::string date_obs;      // DATE-OBS
    
    // Flags for data availability
    bool has_mount_data = false;
    bool has_wcs_data = false;
    bool has_time_data = false;
};

class FITSWCSChecker {
private:
    FITSKeywords keywords;
    
    bool readFITSKeywords(const std::string& filename);
    void printBasicInfo() const;
    void checkJulianDateConsistency() const;
    void checkLocalSiderealTime() const;
    void verifyLSTFromJulianDate() const;
    void calculateExpectedAltAz() const;
    void performBlindCoordinateConversion() const;
    void checkCoordinateConsistency() const;
    void altazToRADec(double alt, double az, double lat, double lst, 
                      double& ra, double& dec) const;
    void radecToAltAz(double ra, double dec, double lat, double lst, 
                      double& alt, double& az) const;
    void blindAltAzToEquatorial(double alt, double az, double latitude, double lst,
                               double& ra, double& dec, double& ha) const;
    double calculateAtmosphericRefraction(double apparent_alt_degrees) const;
    double dateObsToJulian(const std::string& date_obs) const;
    double calculateLSTFromJulian(double julian, double longitude) const;
    double calculateAngularSeparation(double ra1, double dec1, 
                                    double ra2, double dec2) const;
    double normalizeRA(double ra) const;
    double normalizeHA(double ha) const;
    
public:
    bool analyze(const std::string& filename);
};

bool FITSWCSChecker::readFITSKeywords(const std::string& filename) {
    fitsfile* fptr;
    int status = 0;
    char comment[FLEN_COMMENT];
    
    // Open FITS file
    if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)) {
        std::cerr << "Error: Cannot open FITS file " << filename << std::endl;
        fits_report_error(stderr, status);
        return false;
    }
    
    std::cout << "Reading FITS file: " << filename << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Read image dimensions
    if (fits_read_key(fptr, TINT, "NAXIS1", &keywords.naxis1, comment, &status) == 0) {
        // NAXIS1 found
    }
    if (fits_read_key(fptr, TINT, "NAXIS2", &keywords.naxis2, comment, &status) == 0) {
        // NAXIS2 found
    }
    status = 0; // Reset status
    
    // Read instrument info
    char instrume[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "INSTRUME", instrume, comment, &status) == 0) {
        keywords.instrume = instrume;
    }
    status = 0;
    
    char date_obs[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "DATE-OBS", date_obs, comment, &status) == 0) {
        keywords.date_obs = date_obs;
    }
    status = 0;
    
    // Read mount coordinates
    if (fits_read_key(fptr, TDOUBLE, "STELLALT", &keywords.stellalt, comment, &status) == 0) {
        keywords.has_mount_data = true;
    }
    status = 0;
    
    if (fits_read_key(fptr, TDOUBLE, "STELLAZ", &keywords.stellaz, comment, &status) == 0) {
        // STELLAZ found
    }
    status = 0;
    
    // Read observer location
    if (fits_read_key(fptr, TDOUBLE, "OBSLAT", &keywords.obslat, comment, &status) == 0) {
        // OBSLAT found
    }
    status = 0;
    
    if (fits_read_key(fptr, TDOUBLE, "OBSLONG", &keywords.obslong, comment, &status) == 0) {
        // OBSLONG found
    }
    status = 0;
    
    // Read time information
    if (fits_read_key(fptr, TDOUBLE, "JULIAN", &keywords.julian, comment, &status) == 0) {
        keywords.has_time_data = true;
    }
    status = 0;
    
    if (fits_read_key(fptr, TDOUBLE, "LST", &keywords.lst, comment, &status) == 0) {
        // LST found
    }
    status = 0;
    
    if (fits_read_key(fptr, TDOUBLE, "HOURANG", &keywords.hourang, comment, &status) == 0) {
        // HOURANG found
    }
    status = 0;
    
    // Read WCS solution
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &keywords.crval1, comment, &status) == 0) {
        keywords.has_wcs_data = true;
    }
    status = 0;
    
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &keywords.crval2, comment, &status) == 0) {
        // CRVAL2 found
    }
    status = 0;
    
    // Close FITS file
    fits_close_file(fptr, &status);
    
    return true;
}

void FITSWCSChecker::printBasicInfo() const {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Basic FITS Information:" << std::endl;
    std::cout << std::string(25, '-') << std::endl;
    std::cout << "Image dimensions: " << keywords.naxis1 << " x " << keywords.naxis2 << std::endl;
    std::cout << "Instrument: " << keywords.instrume << std::endl;
    std::cout << "Date-Obs: " << keywords.date_obs << std::endl;
    std::cout << std::endl;
    
    if (keywords.has_mount_data) {
        std::cout << "✅ Mount coordinate data found" << std::endl;
    } else {
        std::cout << "❌ Mount coordinate data missing" << std::endl;
    }
    
    if (keywords.has_wcs_data) {
        std::cout << "✅ WCS solution data found" << std::endl;
    } else {
        std::cout << "❌ WCS solution data missing" << std::endl;
    }
    
    if (keywords.has_time_data) {
        std::cout << "✅ Time information found" << std::endl;
    } else {
        std::cout << "❌ Time information missing" << std::endl;
    }
    
    std::cout << std::endl;
}

double FITSWCSChecker::dateObsToJulian(const std::string& date_obs) const {
    // Parse DATE-OBS format: '2024-01-09T22:13:29'
    if (date_obs.length() < 19) {
        return 0.0;  // Invalid format
    }
    
    // Extract components
    int year = std::stoi(date_obs.substr(0, 4));
    int month = std::stoi(date_obs.substr(5, 2));
    int day = std::stoi(date_obs.substr(8, 2));
    int hour = std::stoi(date_obs.substr(11, 2));
    int minute = std::stoi(date_obs.substr(14, 2));
    int second = std::stoi(date_obs.substr(17, 2));
    
    // Use libnova to calculate Julian date
    struct ln_date date;
    date.years = year;
    date.months = month;
    date.days = day;
    date.hours = hour;
    date.minutes = minute;
    date.seconds = second;
    
    return ln_get_julian_day(&date);
}

void FITSWCSChecker::checkJulianDateConsistency() const {
    std::cout << "Julian Date Consistency Analysis:" << std::endl;
    std::cout << std::string(35, '-') << std::endl;
    
    if (keywords.date_obs.empty()) {
        std::cout << "❌ DATE-OBS not available" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    if (keywords.julian == 0.0) {
        std::cout << "❌ JULIAN keyword not available" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << std::fixed << std::setprecision(8);
    
    // Parse DATE-OBS and convert to Julian date
    double calculated_julian = dateObsToJulian(keywords.date_obs);
    
    if (calculated_julian == 0.0) {
        std::cout << "❌ Failed to parse DATE-OBS: " << keywords.date_obs << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << "Time Information:" << std::endl;
    std::cout << "  DATE-OBS: " << keywords.date_obs << std::endl;
    std::cout << "  Calculated Julian Date: " << calculated_julian << std::endl;
    std::cout << "  FITS JULIAN keyword: " << keywords.julian << std::endl;
    std::cout << std::endl;
    
    // Calculate difference
    double julian_diff = fabs(calculated_julian - keywords.julian);
    double julian_diff_seconds = julian_diff * 86400.0;  // Convert to seconds
    
    std::cout << "Julian Date Comparison:" << std::endl;
    std::cout << "  Difference: " << std::setprecision(8) << julian_diff << " days" << std::endl;
    std::cout << "  Difference: " << std::setprecision(3) << julian_diff_seconds << " seconds" << std::endl;
    std::cout << std::endl;
    
    // Evaluate consistency
    std::cout << "Time Consistency Assessment:" << std::endl;
    if (julian_diff_seconds < 1.0) {
        std::cout << "✅ Excellent time consistency (<1 second)" << std::endl;
    } else if (julian_diff_seconds < 60.0) {
        std::cout << "✅ Good time consistency (<1 minute)" << std::endl;
    } else if (julian_diff_seconds < 3600.0) {
        std::cout << "⚠️  Moderate time consistency (<1 hour)" << std::endl;
    } else {
        std::cout << "❌ Poor time consistency (>1 hour)" << std::endl;
    }
    
    // Additional validation using libnova date conversion
    struct ln_date converted_date;
    ln_get_date(keywords.julian, &converted_date);
    
    std::cout << std::endl;
    std::cout << "Julian Date Verification:" << std::endl;
    std::cout << "  JULIAN → Date: " << converted_date.years << "-" 
              << std::setfill('0') << std::setw(2) << converted_date.months << "-"
              << std::setw(2) << converted_date.days << "T"
              << std::setw(2) << converted_date.hours << ":"
              << std::setw(2) << converted_date.minutes << ":"
              << std::setw(2) << (int)converted_date.seconds << std::endl;
    
    // Check if DATE-OBS is UTC (should be for astronomical observations)
    if (keywords.date_obs.find("T") != std::string::npos) {
        std::cout << "  DATE-OBS appears to be in ISO 8601 format (likely UTC)" << std::endl;
        std::cout << "✅ Proper timestamp format detected" << std::endl;
    } else {
        std::cout << "⚠️  DATE-OBS format may not be standard ISO 8601" << std::endl;
    }
    
    std::cout << std::endl;
}

double FITSWCSChecker::calculateLSTFromJulian(double julian, double longitude) const {
    // Calculate Greenwich Mean Sidereal Time (GMST) from Julian date
    double gmst = ln_get_mean_sidereal_time(julian);
    
    // Convert to Local Sidereal Time (LST)
    // LST = GMST + longitude (in hours, positive east)
    double lst = gmst + (longitude / 15.0);
    
    // Normalize to [0, 24) hours
    while (lst < 0.0) lst += 24.0;
    while (lst >= 24.0) lst -= 24.0;
    
    return lst;
}

void FITSWCSChecker::checkLocalSiderealTime() const {
    std::cout << "Local Sidereal Time Verification:" << std::endl;
    std::cout << std::string(35, '-') << std::endl;
    
    if (keywords.julian == 0.0 || keywords.obslong == 0.0) {
        std::cout << "❌ Insufficient data for LST calculation (need JULIAN and OBSLONG)" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    if (keywords.lst == 0.0) {
        std::cout << "❌ LST keyword not available" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << std::fixed << std::setprecision(8);
    
    // Calculate LST from Julian date and longitude
    double calculated_lst = calculateLSTFromJulian(keywords.julian, keywords.obslong);
    
    std::cout << "Sidereal Time Information:" << std::endl;
    std::cout << "  Observer Longitude: " << std::setprecision(6) << keywords.obslong << "°" << std::endl;
    std::cout << "  Julian Date: " << std::setprecision(8) << keywords.julian << std::endl;
    std::cout << "  Calculated LST: " << calculated_lst << " hours" << std::endl;
    std::cout << "  FITS LST keyword: " << keywords.lst << " hours" << std::endl;
    std::cout << std::endl;
    
    // Calculate difference
    double lst_diff = fabs(calculated_lst - keywords.lst);
    
    // Handle 24-hour wraparound
    if (lst_diff > 12.0) {
        lst_diff = 24.0 - lst_diff;
    }
    
    double lst_diff_seconds = lst_diff * 3600.0;  // Convert to seconds
    double lst_diff_arcsec = lst_diff * 15.0 * 3600.0;  // Convert to arcseconds on sky
    
    std::cout << "LST Comparison:" << std::endl;
    std::cout << "  Difference: " << std::setprecision(6) << lst_diff << " hours" << std::endl;
    std::cout << "  Difference: " << std::setprecision(3) << lst_diff_seconds << " seconds" << std::endl;
    std::cout << "  Sky equivalent: " << std::setprecision(1) << lst_diff_arcsec << " arcseconds" << std::endl;
    std::cout << std::endl;
    
    // Evaluate consistency
    std::cout << "LST Consistency Assessment:" << std::endl;
    if (lst_diff_seconds < 1.0) {
        std::cout << "✅ Excellent LST consistency (<1 second)" << std::endl;
    } else if (lst_diff_seconds < 60.0) {
        std::cout << "✅ Good LST consistency (<1 minute)" << std::endl;
    } else if (lst_diff_seconds < 3600.0) {
        std::cout << "⚠️  Moderate LST consistency (<1 hour)" << std::endl;
    } else {
        std::cout << "❌ Poor LST consistency (>1 hour)" << std::endl;
    }
    
    // Additional GMST information
    double gmst = ln_get_mean_sidereal_time(keywords.julian);
    std::cout << std::endl;
    std::cout << "Additional Time Information:" << std::endl;
    std::cout << "  Greenwich Mean Sidereal Time: " << std::setprecision(8) << gmst << " hours" << std::endl;
    std::cout << "  Longitude correction: " << std::setprecision(6) << (keywords.obslong / 15.0) << " hours" << std::endl;
    
    std::cout << std::endl;
}

void FITSWCSChecker::verifyLSTFromJulianDate() const {
    std::cout << "LST vs Julian Date Cross-Verification:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    if (keywords.julian == 0.0 || keywords.lst == 0.0 || keywords.obslong == 0.0) {
        std::cout << "❌ Insufficient data for LST-JD verification (need JULIAN, LST, and OBSLONG)" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << std::fixed << std::setprecision(8);
    
    // Method 1: Calculate LST from Julian Date
    double lst_from_jd = calculateLSTFromJulian(keywords.julian, keywords.obslong);
    
    // Method 2: Use libnova functions for independent verification
    double gmst = ln_get_mean_sidereal_time(keywords.julian);
    double lst_from_gmst = gmst + (keywords.obslong / 15.0);
    while (lst_from_gmst < 0.0) lst_from_gmst += 24.0;
    while (lst_from_gmst >= 24.0) lst_from_gmst -= 24.0;
    
    // Method 3: Calculate apparent sidereal time (includes nutation)
    double apparent_gmst = ln_get_apparent_sidereal_time(keywords.julian);
    double lst_apparent = apparent_gmst + (keywords.obslong / 15.0);
    while (lst_apparent < 0.0) lst_apparent += 24.0;
    while (lst_apparent >= 24.0) lst_apparent -= 24.0;
    
    std::cout << "LST Calculation Methods:" << std::endl;
    std::cout << "  FITS LST keyword: " << keywords.lst << " hours" << std::endl;
    std::cout << "  LST from JD (mean): " << lst_from_jd << " hours" << std::endl;
    std::cout << "  LST from GMST: " << lst_from_gmst << " hours" << std::endl;
    std::cout << "  LST apparent (with nutation): " << lst_apparent << " hours" << std::endl;
    std::cout << std::endl;
    
    // Calculate differences
    double diff_mean = fabs(keywords.lst - lst_from_jd);
    double diff_gmst = fabs(keywords.lst - lst_from_gmst);
    double diff_apparent = fabs(keywords.lst - lst_apparent);
    
    // Handle 24-hour wraparound
    if (diff_mean > 12.0) diff_mean = 24.0 - diff_mean;
    if (diff_gmst > 12.0) diff_gmst = 24.0 - diff_gmst;
    if (diff_apparent > 12.0) diff_apparent = 24.0 - diff_apparent;
    
    std::cout << "LST Differences:" << std::endl;
    std::cout << "  vs Mean LST: " << std::setprecision(6) << diff_mean << " hours (" 
              << std::setprecision(3) << (diff_mean * 3600.0) << " seconds)" << std::endl;
    std::cout << "  vs GMST method: " << std::setprecision(6) << diff_gmst << " hours (" 
              << std::setprecision(3) << (diff_gmst * 3600.0) << " seconds)" << std::endl;
    std::cout << "  vs Apparent LST: " << std::setprecision(6) << diff_apparent << " hours (" 
              << std::setprecision(3) << (diff_apparent * 3600.0) << " seconds)" << std::endl;
    std::cout << std::endl;
    
    // Determine which method the FITS file likely used
    double min_diff = std::min({diff_mean, diff_gmst, diff_apparent});
    std::string method_used;
    
    if (min_diff == diff_mean) {
        method_used = "Mean Sidereal Time";
    } else if (min_diff == diff_gmst) {
        method_used = "GMST conversion";
    } else {
        method_used = "Apparent Sidereal Time (with nutation)";
    }
    
    std::cout << "LST Method Assessment:" << std::endl;
    std::cout << "  Best match: " << method_used << std::endl;
    std::cout << "  Minimum difference: " << std::setprecision(3) << (min_diff * 3600.0) << " seconds" << std::endl;
    std::cout << std::endl;
    
    // Evaluate overall consistency
    if (min_diff * 3600.0 < 1.0) {
        std::cout << "✅ Excellent LST-JD consistency (<1 second)" << std::endl;
    } else if (min_diff * 3600.0 < 10.0) {
        std::cout << "✅ Good LST-JD consistency (<10 seconds)" << std::endl;
    } else if (min_diff * 3600.0 < 60.0) {
        std::cout << "⚠️  Moderate LST-JD consistency (<1 minute)" << std::endl;
    } else {
        std::cout << "❌ Poor LST-JD consistency (>1 minute)" << std::endl;
    }
    
    // Additional astronomical information
    std::cout << std::endl;
    std::cout << "Astronomical Details:" << std::endl;
    std::cout << "  Julian Date: " << std::setprecision(8) << keywords.julian << std::endl;
    std::cout << "  Days since J2000.0: " << std::setprecision(3) << (keywords.julian - 2451545.0) << std::endl;
    std::cout << "  Greenwich Mean Sidereal Time: " << std::setprecision(8) << gmst << " hours" << std::endl;
    std::cout << "  Greenwich Apparent Sidereal Time: " << std::setprecision(8) << apparent_gmst << " hours" << std::endl;
    std::cout << "  Equation of Equinoxes: " << std::setprecision(6) << (apparent_gmst - gmst) * 3600.0 << " seconds" << std::endl;
    std::cout << "  Longitude correction: " << std::setprecision(6) << (keywords.obslong / 15.0) << " hours" << std::endl;
    
    // Check if the difference between mean and apparent is significant
    double nutation_effect = fabs(apparent_gmst - gmst) * 3600.0;
    if (nutation_effect > 1.0) {
        std::cout << "  ⚠️  Significant nutation effect (" << std::setprecision(3) << nutation_effect << " seconds)" << std::endl;
    } else {
        std::cout << "  ✅ Nutation effect is minimal (" << std::setprecision(3) << nutation_effect << " seconds)" << std::endl;
    }
    
    std::cout << std::endl;
}

void FITSWCSChecker::radecToAltAz(double ra, double dec, double lat, double lst, 
                                  double& alt, double& az) const {
    // Convert degrees to radians
    double ra_rad = ra * M_PI / 180.0;
    double dec_rad = dec * M_PI / 180.0;
    double lat_rad = lat * M_PI / 180.0;
    
    // Calculate hour angle
    double ha_hours = lst - (ra / 15.0);
    
    // Normalize hour angle to [-12, 12] range
    while (ha_hours > 12.0) ha_hours -= 24.0;
    while (ha_hours < -12.0) ha_hours += 24.0;
    
    double ha_rad = ha_hours * M_PI / 12.0;
    
    // Calculate altitude
    double alt_rad = asin(sin(dec_rad) * sin(lat_rad) + 
                         cos(dec_rad) * cos(lat_rad) * cos(ha_rad));
    
    // Calculate azimuth using standard astronomical method
    double sin_az = -cos(dec_rad) * sin(ha_rad);
    double cos_az = (sin(dec_rad) * cos(lat_rad) - 
                     cos(dec_rad) * sin(lat_rad) * cos(ha_rad)) / cos(alt_rad);
    
    double az_rad = atan2(sin_az, cos_az);
    
    // Convert to degrees
    alt = alt_rad * 180.0 / M_PI;
    az = az_rad * 180.0 / M_PI;
    
    // Normalize azimuth to [0, 360) range and apply mount convention
    while (az < 0.0) az += 360.0;
    while (az >= 360.0) az -= 360.0;
}

void FITSWCSChecker::calculateExpectedAltAz() const {
    std::cout << "Expected Alt/Az from DATE-OBS and WCS:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    if (keywords.date_obs.empty() || keywords.crval1 == 0.0 || keywords.crval2 == 0.0 ||
        keywords.obslat == 0.0 || keywords.obslong == 0.0) {
        std::cout << "❌ Insufficient data for Alt/Az calculation" << std::endl;
        std::cout << "   Need: DATE-OBS, CRVAL1, CRVAL2, OBSLAT, OBSLONG" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << std::fixed << std::setprecision(6);
    
    // Calculate Julian Date from DATE-OBS
    double julian_from_dateobs = dateObsToJulian(keywords.date_obs);
    if (julian_from_dateobs == 0.0) {
        std::cout << "❌ Failed to parse DATE-OBS" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    // Calculate LST from DATE-OBS
    double lst_from_dateobs = calculateLSTFromJulian(julian_from_dateobs, keywords.obslong);
    
    // Calculate expected Alt/Az from WCS RA/Dec
    double expected_alt, expected_az;
    radecToAltAz(keywords.crval1, keywords.crval2, keywords.obslat, lst_from_dateobs,
                 expected_alt, expected_az);
    
    std::cout << "Calculation Inputs:" << std::endl;
    std::cout << "  DATE-OBS: " << keywords.date_obs << std::endl;
    std::cout << "  Julian Date from DATE-OBS: " << std::setprecision(8) << julian_from_dateobs << std::endl;
    std::cout << "  LST from DATE-OBS: " << lst_from_dateobs << " hours" << std::endl;
    std::cout << "  WCS RA (CRVAL1): " << std::setprecision(6) << keywords.crval1 << "°" << std::endl;
    std::cout << "  WCS Dec (CRVAL2): " << keywords.crval2 << "°" << std::endl;
    std::cout << "  Observer Latitude: " << keywords.obslat << "°" << std::endl;
    std::cout << "  Observer Longitude: " << keywords.obslong << "°" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Expected vs Reported Alt/Az:" << std::endl;
    std::cout << "  Expected Altitude: " << expected_alt << "°" << std::endl;
    std::cout << "  Reported Altitude (STELLALT): " << keywords.stellalt << "°" << std::endl;
    std::cout << "  Expected Azimuth: " << expected_az << "°" << std::endl;
    std::cout << "  Reported Azimuth (STELLAZ): " << keywords.stellaz << "°" << std::endl;
    std::cout << std::endl;
    
    // Calculate differences
    double alt_diff = fabs(expected_alt - keywords.stellalt);
    double az_diff = fabs(expected_az - keywords.stellaz);
    
    // Handle azimuth wraparound (0°/360°)
    if (az_diff > 180.0) {
        az_diff = 360.0 - az_diff;
    }
    
    std::cout << "Alt/Az Differences:" << std::endl;
    std::cout << "  Altitude difference: " << alt_diff << "° (" << (alt_diff * 3600.0) << " arcsec)" << std::endl;
    std::cout << "  Azimuth difference: " << az_diff << "° (" << (az_diff * 3600.0) << " arcsec)" << std::endl;
    
    // Calculate total angular separation in Alt/Az coordinates
    double altaz_separation = sqrt(alt_diff * alt_diff + 
                                  (az_diff * cos(expected_alt * M_PI / 180.0)) * 
                                  (az_diff * cos(expected_alt * M_PI / 180.0)));
    
    std::cout << "  Total Alt/Az separation: " << altaz_separation << "° (" 
              << (altaz_separation * 3600.0) << " arcsec)" << std::endl;
    std::cout << std::endl;
    
    // Evaluate consistency
    std::cout << "Alt/Az Consistency Assessment:" << std::endl;
    if (altaz_separation * 3600.0 < 60.0) {
        std::cout << "✅ Excellent Alt/Az consistency (<1 arcmin)" << std::endl;
    } else if (altaz_separation * 3600.0 < 300.0) {
        std::cout << "✅ Good Alt/Az consistency (<5 arcmin)" << std::endl;
    } else if (altaz_separation * 3600.0 < 1800.0) {
        std::cout << "⚠️  Moderate Alt/Az consistency (<30 arcmin)" << std::endl;
    } else {
        std::cout << "❌ Poor Alt/Az consistency (>30 arcmin)" << std::endl;
    }
    
    // Additional information about the observation
    std::cout << std::endl;
    std::cout << "Observation Details:" << std::endl;
    
    // Calculate hour angle
    double ha_hours = lst_from_dateobs - (keywords.crval1 / 15.0);
    while (ha_hours > 12.0) ha_hours -= 24.0;
    while (ha_hours < -12.0) ha_hours += 24.0;
    
    std::cout << "  Hour Angle: " << std::setprecision(6) << ha_hours << " hours" << std::endl;
    
    // Determine if object is rising or setting
    if (ha_hours < 0) {
        std::cout << "  Object status: Rising (east of meridian)" << std::endl;
    } else if (ha_hours > 0) {
        std::cout << "  Object status: Setting (west of meridian)" << std::endl;
    } else {
        std::cout << "  Object status: At meridian" << std::endl;
    }
    
    // Calculate airmass (sec(zenith angle))
    double zenith_angle = 90.0 - expected_alt;
    double airmass = 1.0 / cos(zenith_angle * M_PI / 180.0);
    
    std::cout << "  Altitude: " << expected_alt << "° (zenith angle: " << zenith_angle << "°)" << std::endl;
    std::cout << "  Airmass: " << std::setprecision(2) << airmass << std::endl;
    
    if (expected_alt < 30.0) {
        std::cout << "  ⚠️  Low altitude observation (high airmass)" << std::endl;
    } else if (expected_alt > 70.0) {
        std::cout << "  ✅ High altitude observation (low airmass)" << std::endl;
    } else {
        std::cout << "  ✅ Moderate altitude observation" << std::endl;
    }
    
    // Check for potential refraction effects
    if (expected_alt < 20.0) {
        std::cout << "  ⚠️  Significant atmospheric refraction expected" << std::endl;
    }
    
    std::cout << std::endl;
}

static double FITSWCSChecker::calculateAtmosphericRefraction(double apparent_alt_degrees) const {
    // Bennett's refraction formula (1982) - more accurate than simple formulas
    // Returns refraction in degrees
    
    if (apparent_alt_degrees < -1.0) {
        return 0.0;  // Below horizon
    }
    
    double refraction_arcsec;
    
    if (apparent_alt_degrees > 15.0) {
        // For elevations above 15°
        double alt_rad = apparent_alt_degrees * M_PI / 180.0;
        double tan_alt = tan(alt_rad);
        refraction_arcsec = 58.1 / tan_alt - 0.07 / (tan_alt * tan_alt * tan_alt) + 
                           0.000086 / (tan_alt * tan_alt * tan_alt * tan_alt * tan_alt);
    } else {
        // For low elevations (more complex formula needed)
        double h = apparent_alt_degrees;
        refraction_arcsec = 1735.0 + h * (-518.2 + h * (103.4 + h * (-12.79 + h * 0.711)));
        refraction_arcsec = refraction_arcsec / 3600.0;  // Convert to degrees
        return refraction_arcsec;
    }
    
    return refraction_arcsec / 3600.0;  // Convert arcseconds to degrees
}

void FITSWCSChecker::blindAltAzToEquatorial(double alt, double az, double latitude, double lst,
                                           double& ra, double& dec, double& ha) const {
    // Apply atmospheric refraction correction first
    double true_alt = alt - calculateAtmosphericRefraction(alt);
    
    // Convert to radians (NO azimuth conversion needed - mount reports correct values)
    double alt_rad = true_alt * M_PI / 180.0;
    double az_rad = az * M_PI / 180.0;
    double lat_rad = latitude * M_PI / 180.0;
    
    // Calculate declination using spherical law of cosines
    double sin_dec = sin(alt_rad) * sin(lat_rad) + 
                     cos(alt_rad) * cos(lat_rad) * cos(az_rad);
    
    // Clamp to avoid numerical errors
    sin_dec = std::max(-1.0, std::min(1.0, sin_dec));
    dec = asin(sin_dec);
    
    // Calculate hour angle
    double cos_ha = (sin(alt_rad) - sin(lat_rad) * sin(dec)) / 
                    (cos(lat_rad) * cos(dec));
    
    // Clamp to avoid numerical errors
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    double ha_rad = acos(cos_ha);
    
    // FIXED: Reverse the hour angle sign logic
    if (sin(az_rad) > 0) {  // Azimuth > 180° (western sky)
        ha_rad = -ha_rad;   // NEGATIVE HA (setting) - CORRECTED
    } else {                // Azimuth < 180° (eastern sky)
        ha_rad = ha_rad;    // POSITIVE HA (rising) - CORRECTED
    }
    
    // Convert Hour Angle to hours
    ha = ha_rad * 12.0 / M_PI;
    
    // Calculate Right Ascension using fundamental relation: RA = LST - HA
    double ra_hours = lst - ha;
    
    // Normalize RA to [0, 24) hours
    while (ra_hours < 0) ra_hours += 24.0;
    while (ra_hours >= 24.0) ra_hours -= 24.0;
    
    // Convert to degrees
    ra = ra_hours * 15.0;
    dec = dec * 180.0 / M_PI;
}

void FITSWCSChecker::performBlindCoordinateConversion() const {
    std::cout << "Blind Coordinate Conversion (Alt/Az → RA/Dec):" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    if (!keywords.has_mount_data || keywords.lst == 0.0 || keywords.obslat == 0.0) {
        std::cout << "❌ Insufficient data for blind coordinate conversion" << std::endl;
        std::cout << "   Need: STELLALT, STELLAZ, LST, OBSLAT" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Input Coordinates (from mount):" << std::endl;
    std::cout << "  Altitude: " << keywords.stellalt << "°" << std::endl;
    std::cout << "  Azimuth: " << keywords.stellaz << "°" << std::endl;
    std::cout << "  Observer Latitude: " << keywords.obslat << "°" << std::endl;
    std::cout << "  Local Sidereal Time: " << std::setprecision(8) << keywords.lst << " hours" << std::endl;
    std::cout << std::endl;
    
    // Perform blind coordinate conversion
    double calculated_ra, calculated_dec, calculated_ha;
    blindAltAzToEquatorial(keywords.stellalt, keywords.stellaz, keywords.obslat, 
                          keywords.lst, calculated_ra, calculated_dec, calculated_ha);
    
    std::cout << "Calculated Equatorial Coordinates:" << std::endl;
    std::cout << "  RA: " << std::setprecision(6) << calculated_ra << "°" << std::endl;
    std::cout << "  Dec: " << calculated_dec << "°" << std::endl;
    std::cout << "  Hour Angle: " << std::setprecision(8) << calculated_ha << " hours" << std::endl;
    std::cout << std::endl;
    
    // Compare with WCS solution if available
    if (keywords.has_wcs_data) {
        std::cout << "Comparison with WCS Solution:" << std::endl;
        std::cout << "  WCS RA (CRVAL1): " << std::setprecision(6) << keywords.crval1 << "°" << std::endl;
        std::cout << "  WCS Dec (CRVAL2): " << keywords.crval2 << "°" << std::endl;
        
        // Calculate differences
        double ra_diff = fabs(calculated_ra - keywords.crval1);
        double dec_diff = fabs(calculated_dec - keywords.crval2);
        
        // Handle RA wraparound
        if (ra_diff > 180.0) {
            ra_diff = 360.0 - ra_diff;
        }
        
        std::cout << "  ΔRA: " << ra_diff << "° (" << (ra_diff * 3600.0) << " arcsec)" << std::endl;
        std::cout << "  ΔDec: " << dec_diff << "° (" << (dec_diff * 3600.0) << " arcsec)" << std::endl;
        
        // Calculate total angular separation
        double angular_sep = sqrt(
            (ra_diff * cos(keywords.crval2 * M_PI / 180.0)) * (ra_diff * cos(keywords.crval2 * M_PI / 180.0)) +
            dec_diff * dec_diff
        ) * 3600.0;  // Convert to arcseconds
        
        std::cout << "  Angular separation: " << std::setprecision(1) << angular_sep << " arcsec" << std::endl;
        std::cout << std::endl;
        
        // Assessment
        if (angular_sep < 60.0) {
            std::cout << "✅ Excellent blind conversion accuracy (<1 arcmin)" << std::endl;
        } else if (angular_sep < 300.0) {
            std::cout << "✅ Good blind conversion accuracy (<5 arcmin)" << std::endl;
        } else if (angular_sep < 1800.0) {
            std::cout << "⚠️  Moderate blind conversion accuracy (<30 arcmin)" << std::endl;
        } else {
            std::cout << "❌ Poor blind conversion accuracy (>30 arcmin)" << std::endl;
        }
    }
    
    // Round-trip validation: RA/Dec → Alt/Az
    std::cout << std::endl;
    std::cout << "Round-Trip Validation (RA/Dec → Alt/Az):" << std::endl;
    double roundtrip_alt, roundtrip_az;
    radecToAltAz(calculated_ra, calculated_dec, keywords.obslat, keywords.lst,
                 roundtrip_alt, roundtrip_az);
    
    std::cout << "  Back-calculated Altitude: " << std::setprecision(6) << roundtrip_alt << "°" << std::endl;
    std::cout << "  Back-calculated Azimuth: " << roundtrip_az << "°" << std::endl;
    std::cout << std::endl;
    
    // Round-trip errors
    double alt_roundtrip_error = fabs(keywords.stellalt - roundtrip_alt);
    double az_roundtrip_error = fabs(keywords.stellaz - roundtrip_az);
    
    // Handle azimuth wraparound
    if (az_roundtrip_error > 180.0) {
        az_roundtrip_error = 360.0 - az_roundtrip_error;
    }
    
    std::cout << "Round-Trip Errors:" << std::endl;
    std::cout << "  Altitude error: " << alt_roundtrip_error << "° (" 
              << (alt_roundtrip_error * 3600.0) << " arcsec)" << std::endl;
    std::cout << "  Azimuth error: " << az_roundtrip_error << "° (" 
              << (az_roundtrip_error * 3600.0) << " arcsec)" << std::endl;
    
    double total_roundtrip_error = sqrt(alt_roundtrip_error * alt_roundtrip_error + 
                                       az_roundtrip_error * az_roundtrip_error) * 3600.0;
    
    std::cout << "  Total round-trip error: " << std::setprecision(3) << total_roundtrip_error << " arcsec" << std::endl;
    std::cout << std::endl;
    
    // Round-trip assessment
    std::cout << "Round-Trip Accuracy Assessment:" << std::endl;
    if (total_roundtrip_error < 1.0) {
        std::cout << "✅ Excellent numerical precision (<1 arcsec)" << std::endl;
    } else if (total_roundtrip_error < 10.0) {
        std::cout << "✅ Very good numerical precision (<10 arcsec)" << std::endl;
    } else if (total_roundtrip_error < 60.0) {
        std::cout << "✅ Good numerical precision (<1 arcmin)" << std::endl;
    } else {
        std::cout << "⚠️  Moderate numerical precision (>1 arcmin)" << std::endl;
    }
    
    // Additional information
    std::cout << std::endl;
    std::cout << "Atmospheric Refraction:" << std::endl;
    double refraction = calculateAtmosphericRefraction(keywords.stellalt);
    std::cout << "  Refraction correction: " << std::setprecision(4) << refraction << "° (" 
              << (refraction * 3600.0) << " arcsec)" << std::endl;
    
    if (keywords.stellalt < 20.0) {
        std::cout << "  ⚠️  Low altitude - significant refraction effects" << std::endl;
    } else if (keywords.stellalt > 70.0) {
        std::cout << "  ✅ High altitude - minimal refraction effects" << std::endl;
    } else {
        std::cout << "  ✅ Moderate altitude - small refraction effects" << std::endl;
    }
    
    std::cout << std::endl;
}

void FITSWCSChecker::altazToRADec(double alt, double az, double lat, double lst, 
                                  double& ra, double& dec) const {
    // Convert degrees to radians
    double alt_rad = alt * M_PI / 180.0;
    double az_rad = az * M_PI / 180.0;
    double lat_rad = lat * M_PI / 180.0;
    
    // Calculate declination
    double dec_rad = asin(
        sin(alt_rad) * sin(lat_rad) + 
        cos(alt_rad) * cos(lat_rad) * cos(az_rad)
    );
    
    // Calculate hour angle
    double cos_ha = (sin(alt_rad) - sin(dec_rad) * sin(lat_rad)) / (
        cos(dec_rad) * cos(lat_rad)
    );
    
    // Clamp to valid range to avoid numerical errors
    cos_ha = std::max(-1.0, std::min(1.0, cos_ha));
    
    double ha_rad = acos(cos_ha);
    
    // Determine sign of hour angle based on azimuth
    if (sin(az_rad) > 0) {  // Azimuth > 180° (west)
        ha_rad = -ha_rad;
    }
    
    // Convert hour angle to hours
    double ha_hours = ha_rad * 12.0 / M_PI;
    
    // Calculate RA from LST and hour angle
    double ra_hours = lst - ha_hours;
    
    // Normalize RA to [0, 24) range
    while (ra_hours < 0) {
        ra_hours += 24;
    }
    while (ra_hours >= 24) {
        ra_hours -= 24;
    }
    
    // Convert to degrees
    ra = ra_hours * 15.0;
    dec = dec_rad * 180.0 / M_PI;
}

double FITSWCSChecker::calculateAngularSeparation(double ra1, double dec1, 
                                                 double ra2, double dec2) const {
    // Convert to radians
    double ra1_rad = ra1 * M_PI / 180.0;
    double dec1_rad = dec1 * M_PI / 180.0;
    double ra2_rad = ra2 * M_PI / 180.0;
    double dec2_rad = dec2 * M_PI / 180.0;
    
    // Calculate angular separation using spherical law of cosines
    double cos_sep = sin(dec1_rad) * sin(dec2_rad) + 
                     cos(dec1_rad) * cos(dec2_rad) * cos(ra1_rad - ra2_rad);
    
    // Clamp to valid range to avoid numerical errors
    cos_sep = std::max(-1.0, std::min(1.0, cos_sep));
    
    double sep_rad = acos(cos_sep);
    return sep_rad * 180.0 / M_PI;  // Convert to degrees
}

double FITSWCSChecker::normalizeRA(double ra) const {
    while (ra < 0) ra += 24.0;
    while (ra >= 24.0) ra -= 24.0;
    return ra;
}

double FITSWCSChecker::normalizeHA(double ha) const {
    while (ha > 12.0) ha -= 24.0;
    while (ha < -12.0) ha += 24.0;
    return ha;
}

void FITSWCSChecker::checkCoordinateConsistency() const {
    std::cout << "Coordinate System Consistency Analysis:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    if (!keywords.has_mount_data || !keywords.has_wcs_data || !keywords.has_time_data) {
        std::cout << "❌ Insufficient data for coordinate consistency check" << std::endl;
        std::cout << std::endl;
        return;
    }
    
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Mount and Observatory Data:" << std::endl;
    std::cout << "  Altitude: " << keywords.stellalt << "°" << std::endl;
    std::cout << "  Azimuth: " << keywords.stellaz << "°" << std::endl;
    std::cout << "  Observer Latitude: " << keywords.obslat << "°" << std::endl;
    std::cout << "  Observer Longitude: " << keywords.obslong << "°" << std::endl;
    std::cout << "  Julian Date: " << std::setprecision(8) << keywords.julian << std::endl;
    std::cout << "  Local Sidereal Time: " << std::setprecision(8) << keywords.lst << " hours" << std::endl;
    std::cout << "  Hour Angle: " << std::setprecision(8) << keywords.hourang << " hours" << std::endl;
    std::cout << std::endl;
    
    std::cout << std::setprecision(6);
    std::cout << "WCS Solution:" << std::endl;
    std::cout << "  RA (CRVAL1): " << keywords.crval1 << "°" << std::endl;
    std::cout << "  Dec (CRVAL2): " << keywords.crval2 << "°" << std::endl;
    std::cout << std::endl;
    
    // Convert Alt/Az to RA/Dec
    double calculated_ra, calculated_dec;
    altazToRADec(keywords.stellalt, keywords.stellaz, keywords.obslat, 
                 keywords.lst, calculated_ra, calculated_dec);
    
    std::cout << "Calculated RA/Dec from Alt/Az:" << std::endl;
    std::cout << "  RA: " << calculated_ra << "°" << std::endl;
    std::cout << "  Dec: " << calculated_dec << "°" << std::endl;
    std::cout << std::endl;
    
    // Calculate differences
    double ra_diff = fabs(calculated_ra - keywords.crval1);
    double dec_diff = fabs(calculated_dec - keywords.crval2);
    
    // Handle RA wraparound (0°/360°)
    if (ra_diff > 180.0) {
        ra_diff = 360.0 - ra_diff;
    }
    
    std::cout << "Coordinate Differences:" << std::endl;
    std::cout << "  ΔRA: " << ra_diff << "° (" << (ra_diff * 3600.0) << " arcsec)" << std::endl;
    std::cout << "  ΔDec: " << dec_diff << "° (" << (dec_diff * 3600.0) << " arcsec)" << std::endl;
    
    // Convert to angular separation
    double angular_sep = sqrt(
        (ra_diff * cos(keywords.crval2 * M_PI / 180.0)) * (ra_diff * cos(keywords.crval2 * M_PI / 180.0)) + 
        dec_diff * dec_diff
    ) * 3600.0;  // arcsec
    
    std::cout << "  Angular separation: " << std::setprecision(1) << angular_sep << " arcsec" << std::endl;
    std::cout << std::endl;
    
    // Evaluate consistency
    if (angular_sep < 60) {  // 1 arcminute
        std::cout << "✅ Excellent coordinate consistency (<1 arcmin)" << std::endl;
    } else if (angular_sep < 300) {  // 5 arcminutes
        std::cout << "✅ Good coordinate consistency (<5 arcmin)" << std::endl;
    } else if (angular_sep < 1800) {  // 30 arcminutes
        std::cout << "⚠️  Moderate coordinate consistency (<30 arcmin)" << std::endl;
    } else {
        std::cout << "❌ Poor coordinate consistency (>30 arcmin)" << std::endl;
    }
    
    // Check hour angle consistency
    double ra_hours = keywords.crval1 / 15.0;  // Convert RA to hours
    double calculated_ha = keywords.lst - ra_hours;
    
    // Normalize hour angle to [-12, 12] range
    while (calculated_ha > 12) {
        calculated_ha -= 24;
    }
    while (calculated_ha < -12) {
        calculated_ha += 24;
    }
    
    double ha_diff = fabs(calculated_ha - keywords.hourang);
    if (ha_diff > 12) {
        ha_diff = 24 - ha_diff;
    }
    
    std::cout << std::endl;
    std::cout << "Hour Angle Verification:" << std::endl;
    std::cout << "  Recorded hour angle: " << std::setprecision(6) << keywords.hourang << " hours" << std::endl;
    std::cout << "  Calculated hour angle: " << calculated_ha << " hours" << std::endl;
    std::cout << "  Difference: " << ha_diff << " hours (" << (ha_diff * 3600.0) << " arcsec)" << std::endl;
    
    if (ha_diff < 0.01) {  // ~36 arcsec
        std::cout << "✅ Hour angle is consistent" << std::endl;
    } else {
        std::cout << "⚠️  Hour angle discrepancy detected" << std::endl;
    }
    
    std::cout << std::endl;
}

bool FITSWCSChecker::analyze(const std::string& filename) {
    if (!readFITSKeywords(filename)) {
        return false;
    }
    
    printBasicInfo();
    checkJulianDateConsistency();
    checkLocalSiderealTime();
    verifyLSTFromJulianDate();
    calculateExpectedAltAz();
    performBlindCoordinateConversion();
    checkCoordinateConsistency();
    
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <fits_file_path>" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    
    std::cout << "FITS WCS Consistency Checker (C++ Version)" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << std::endl;
    
    FITSWCSChecker checker;
    
    if (checker.analyze(filename)) {
        std::cout << "Analysis completed successfully!" << std::endl;
        return 0;
    } else {
        std::cout << "Analysis failed - see errors above" << std::endl;
        return 1;
    }
}
