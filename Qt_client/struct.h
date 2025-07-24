#ifndef STELLINASTRUCT_H
#define STELLINASTRUCT_H

// Processing modes
enum ProcessingMode {
    MODE_BASIC_PLATESOLVE = 0,
    MODE_DARK_CALIBRATION = 1,
    MODE_ASTROMETRIC_STACKING = 2,
    MODE_FULL_PIPELINE = 3
};

// Stacking parameters - unified structure for both WCS and traditional stacking
struct StackingParams {
    enum CombinationMethod {
        MEAN,
        MEDIAN, 
        WEIGHTED_MEAN,
        SIGMA_CLIPPED_MEAN,
        MINIMUM,
        MAXIMUM
    };
    
    enum RejectionMethod {
        NO_REJECTION,
        SIGMA_CLIPPING,
        PERCENTILE_CLIPPING,
        LINEAR_FIT_CLIPPING
    };
    
    CombinationMethod combination = WEIGHTED_MEAN;
    RejectionMethod rejection = SIGMA_CLIPPING;
    double sigma_low = 3.0;           // Low sigma clipping threshold
    double sigma_high = 3.0;          // High sigma clipping threshold
    double percentile_low = 5.0;      // Low percentile (%)
    double percentile_high = 95.0;    // High percentile (%)
    bool normalize_exposure = true;   // Normalize by exposure time
    bool apply_flat_correction = false; // Apply flat field correction
    bool apply_brightness_normalization = false; // Apply brightness normalisation
    double output_pixel_scale = 0.0;  // Override pixel scale (0 = auto)
    int output_width = 0;             // Override width (0 = auto)
    int output_height = 0;            // Override height (0 = auto)
    bool create_weight_map = true;    // Generate output weight map
    bool save_intermediate = false;   // Save reprojected images
    QString output_format = "fits";   // Output format
};

// Dark frame information
struct DarkFrame {
    QString filepath;
    int exposure;     // exposure time in seconds
    int temperature;  // sensor temperature in degrees C
    QString binning;  // binning mode (e.g., "1x1", "2x2")
    QString bayerPattern;      // NEW: bayer pattern for this dark frame
    
    DarkFrame() : exposure(0), temperature(0), bayerPattern("RGGB") {}

};

// Update the StellinaImageData structure
struct StellinaImageData {
    QString originalFitsPath;     // Original raw FITS file path
    QString originalJsonPath;     // Original JSON metadata file path  
    QString currentFitsPath;      // Current FITS file path (updated through pipeline)
    QJsonObject metadata;         // Complete JSON metadata
    double altitude;              // Stellina altitude (degrees)
    double azimuth;               // Stellina azimuth (degrees)
    QString dateObs;              // DATE-OBS from FITS header
    bool hasValidCoordinates;     // Whether Alt/Az coordinates are valid
    int exposureSeconds;          // Exposure time in seconds
    int temperatureKelvin;        // Sensor temperature in Kelvin
    QString binning;              // Binning mode (e.g., "1x1", "2x2")
    
    // NEW: Pre-calculated RA/DEC coordinates from coordinate conversion
    double calculatedRA;          // Calculated RA from Alt/Az conversion (degrees)
    double calculatedDec;         // Calculated Dec from Alt/Az conversion (degrees)
    bool hasCalculatedCoords;     // Whether calculated coordinates are available
    
    // NEW: Reversed image support
    bool isReversedImage;         // True if this is a reversed stellina image (img-0001r.fits)
    QString bayerPattern;         // Detected bayer pattern (e.g., "RGGB", "GRBG", "GBRG", "BGGR")
    QString baseName;             // Base name without 'r' suffix (e.g., "img-0001")
    bool hasValidWCS;    
    StellinaImageData() : altitude(0), azimuth(0), hasValidCoordinates(false), 
                         exposureSeconds(0), temperatureKelvin(284), binning("1x1"),
			  calculatedRA(0), calculatedDec(0), hasCalculatedCoords(false), isReversedImage(false), bayerPattern("RGGB")  {}
    
    // Convenience function to check if pre-calculated coordinates are available
    bool hasPreCalculatedCoords() const {
        return hasCalculatedCoords && calculatedRA != 0.0 && calculatedDec != 0.0;
    }
};

struct StackingCorrectionData {
    QString imageFilename;
    int imageNumber;
    qint64 acqTime;
    double minutesFromStart;
    
    // Mount position
    double stellinaAlt, stellinaAz;
    
    // Stacking corrections (pixels)
    double correctionX, correctionY, correctionRot;
    
    // Registration quality metrics
    int starsUsed;
    QString statusMessage;
    double distanceToCenter;
    
    bool isValid;
    
    StackingCorrectionData() : imageNumber(0), minutesFromStart(0.0), 
                              stellinaAlt(0), stellinaAz(0),
                              correctionX(0), correctionY(0), correctionRot(0),
                              starsUsed(0), distanceToCenter(0), isValid(false) {}
};

struct ProcessedImageData {
    QString filename;
    int imageNumber;
    double stellinaAlt, stellinaAz;
    double predictedRA, predictedDec;
    double solvedRA, solvedDec;
    QString dateObs;
    QDateTime obsTime;
    double minutesFromStart;
    bool isValid;
    // ADDED: Missing members that were causing compilation errors
    double pixelScale;     // Pixel scale in arcsec/pixel
    bool hasValidWCS;      // Whether WCS data is valid
    double correctionX;
    double correctionY;
    double correctionRot;
    int starsUsed;
    double distanceToCenter;
    
    ProcessedImageData() : stellinaAlt(0), stellinaAz(0), predictedRA(0), predictedDec(0),
                          solvedRA(0), solvedDec(0), correctionX(0), correctionY(0), 
                          correctionRot(0), starsUsed(0), distanceToCenter(0), 
                          isValid(false), pixelScale(0), hasValidWCS(false) {}
};

struct CoordinateTestCase {
    QString name;
    double lat;           // Observer latitude (degrees)
    double lon;           // Observer longitude (degrees) 
    QString date;         // UTC date/time string
    double julianDay;     // Julian Day
    double raOfDate;      // RA of date (degrees)
    double decOfDate;     // Dec of date (degrees)
    double expectedAz;    // Expected azimuth (degrees)
    double expectedAlt;   // Expected altitude (degrees)
    double siderealTime;  // Local sidereal time (hours)
    double hourAngle;     // Hour angle (hours)
    QString description;  // Test description
};

// New structure to hold either JSON-derived or FITS-derived coordinate data
struct CoordinateSource {
    enum Type {
        FROM_JSON_ALTAZ,    // Original: Alt/Az from JSON, converted to RA/Dec
        FROM_FITS_WCS,      // New: Direct RA/Dec from FITS WCS headers
        INVALID
    };
    
    Type type;
    double ra, dec;         // Final RA/Dec coordinates (degrees)
    double alt, az;         // Alt/Az if available (degrees)
    QString source_info;    // Description of coordinate source
    
    CoordinateSource() : type(INVALID), ra(0), dec(0), alt(0), az(0) {}
};

#endif STELLINASTRUCT_H
