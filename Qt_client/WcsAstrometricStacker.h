// Refactored WcsAstrometricStacker.h with modular stacking support
// Updated to support subframe processing and simplified pixel structures

#ifndef WCS_ASTROMETRIC_STACKER_H
#define WCS_ASTROMETRIC_STACKER_H

#include <QObject>
#include <QString>
#include <QList>
#include <QDateTime>
#include <QProgressBar>
#include <QLabel>
#include <QTimer>
#include <QStringList>
#include <opencv2/opencv.hpp>
#include <fitsio.h>
#include <vector>
#include <memory>
#include <array>
#include <cmath>

// Forward declarations
struct StellinaImageData;
struct PixelAccumulator;

// Stacking parameters - unified structure for both WCS and traditional stacking
typedef struct StackingParams {
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
} StackingParams;

// Simplified pixel contribution structure - only stores essential data
struct PixelContribution {
    float value;
    float weight;
    uint8_t bayer_color;  // 0=R, 1=G1, 2=G2, 3=B (for RGGB pattern)
    
    PixelContribution(float v, float w, uint8_t color) 
        : value(v), weight(w), bayer_color(color) {}
};

// Simple TAN projection WCS implementation
struct SimpleTANWCS {
    double crval1, crval2;  // Reference RA, Dec (degrees)
    double crpix1, crpix2;  // Reference pixels (1-indexed)
    double cd11, cd12, cd21, cd22;  // CD matrix (degrees/pixel)
    bool valid;
    
    SimpleTANWCS() : crval1(0), crval2(0), crpix1(0), crpix2(0),
                     cd11(0), cd12(0), cd21(0), cd22(0), valid(false) {}
    
    bool pixelToWorld(double px, double py, double& ra, double& dec) const;
    bool worldToPixel(double ra, double dec, double& px, double& py) const;
    void printDiagnostics() const;
    double getPixelScale() const; // Returns arcsec/pixel
};

// Enhanced image data structure for WCS stacking
struct WCSImageData {
    cv::Mat image;                    // Image pixel data (32-bit float)
    cv::Mat weight_map;               // Quality weight map
    SimpleTANWCS wcs;                 // Simple TAN WCS instead of wcsprm
    QString filename;                 // Original filename
    QString solved_filename;          // Plate-solved FITS file
    double quality_score;             // Overall quality (0-1)
    double exposure_time;             // Exposure in seconds
    int star_count;                   // Number of detected stars
    double background_level;          // Background ADU level
    double noise_level;               // Image noise estimate
    QDateTime obs_time;               // Observation time
    bool wcs_valid;                   // WCS successfully loaded
    
    // Stellina-specific quality metrics
    double stellina_correction_magnitude; // From stacking JSON
    int stellina_stars_used;             // Stars used in Stellina registration
    
    WCSImageData() : quality_score(1.0), exposure_time(10.0), star_count(0),
                    background_level(0.0), noise_level(0.0), wcs_valid(false),
                    stellina_correction_magnitude(0.0), stellina_stars_used(0) {}
};

class WCSAstrometricStacker : public QObject {
    Q_OBJECT

public:
    explicit WCSAstrometricStacker(QObject *parent = nullptr);
    ~WCSAstrometricStacker();
    
    // Main interface
    bool addImage(const QString &fits_file, const QString &solved_fits_file = "");
    bool addImageWithMetadata(const QString &fits_file, const struct StellinaImageData &stellina_data);
    void setStackingParameters(const StackingParams &params);
    void setProgressWidgets(QProgressBar *progress, QLabel *status);
    
    // Core stacking process - now modular
    bool computeOptimalWCS();
    bool stackImages();
    bool saveResult(const QString &output_path);
    
    // Quality analysis
    void analyzeImageQuality();
    QString getQualityReport() const;
    
    // Access results
    cv::Mat getStackedImage() const { return m_stacked_image; }
    cv::Mat getWeightMap() const { return m_weight_map; }
    cv::Mat getOverlapMap() const { return m_overlap_map; }
    cv::Mat getNoiseMap() const { return m_noise_map; }
    SimpleTANWCS getOutputWCS() const { return m_output_wcs; }
    
    // Statistics
    int getImageCount() const { return m_images.size(); }
    double getTotalExposureTime() const;
    double getAverageQuality() const;
    
    // Processing control
    void startStacking();
    void cancelStacking();
    
    // Modular stacking functions
    bool processImageSubframe(
        const std::unique_ptr<WCSImageData>& img_data,
        std::vector<PixelAccumulator>& pixel_accumulators,
        int start_row, int end_row,
        const QString& bayer_pattern,
        size_t img_idx);
    
    bool finalizeStackedImage(const std::vector<PixelAccumulator>& pixel_accumulators);
    // Post-processing functions
    void generateQualityMaps();
    void applyBrightnessNormalization();
    cv::Size getOutputSize() { return m_output_size; }
    bool loadWCSFromFITS(const QString &fits_file, WCSImageData &img_data);

signals:
    void progressUpdated(int percentage);
    void statusUpdated(const QString &message);
    void imageProcessed(const QString &filename);
    void stackingComplete(bool success);
    void errorOccurred(const QString &error);
    void qualityAnalysisComplete();
    
private:
    // Core processing functions - now modular
    bool extractImageStatistics(WCSImageData &img_data);
    bool computeImageQualityScore(WCSImageData &img_data);
    
    // Pixel processing utilities
    float calculatePixelWeight(
        const std::unique_ptr<WCSImageData>& img_data,
        int target_x, int target_y, float pixel_value);
    
    float applySigmaClipping(const PixelAccumulator& acc, int x, int y);
    
    // Bayer pattern utilities
    QString determineBayerPattern(const QString& filename);
    void logBayerStatistics(const std::vector<PixelAccumulator>& pixel_accumulators);
        
    // Legacy compatibility functions
    bool addPlatesolveDFITSFile(const QString &solved_fits_file);
    bool addImageFromStellinaData(const QString &fits_file, const StellinaImageData &stellina_data);
    void saveOverlapMap(const cv::Mat& overlap_count, const QString& output_path);
    void analyzeOverlapDistribution(const cv::Mat& overlap_count);
    void logOverlapStatistics();
    void applyGlobalBrightnessNormalization();
    void applySigmaClipping(std::vector<PixelContribution>& contributions, int x, int y);
  
    // Utility functions
    void updateProgress(int percentage, const QString &message);
    void logProcessing(const QString &message);
    void finishStacking();
    
    // Member variables
    std::vector<std::unique_ptr<WCSImageData>> m_images;
    SimpleTANWCS m_output_wcs;        // Target WCS for output
    cv::Size m_output_size;           // Output image dimensions
    double m_output_pixel_scale;      // Arcseconds per pixel
    
    // Results
    cv::Mat m_stacked_image;          // Final stacked result
    cv::Mat m_weight_map;             // Combined weight map
    cv::Mat m_overlap_map;            // Number of overlapping images per pixel
    cv::Mat m_noise_map;              // Noise estimate map
    
    // Processing parameters
    StackingParams m_params;
    QStringList m_processing_log;
    
    // UI progress tracking
    QProgressBar *m_progress_bar;
    QLabel *m_status_label;
    
    // Processing state
    bool m_stacking_active;
    int m_current_image_index;
    QTimer *m_processing_timer;
    
    // Statistics
    double m_total_processing_time;
    int m_pixels_processed;
    int m_pixels_rejected;
    
    // Modular processing constants
    static const int DEFAULT_SUBFRAME_HEIGHT = 256;
    static const int MIN_CONTRIBUTIONS_FOR_SIGMA_CLIPPING = 3;
};

// Bayer pattern utilities (can be used by other classes)
namespace BayerUtils {
    enum BayerColor { RED = 0, GREEN1 = 1, GREEN2 = 2, BLUE = 3 };
    
    uint8_t getBayerColor(int x, int y, const QString& pattern);
    QString normalizeBayerPattern(const QString& pattern);
    bool isValidBayerPattern(const QString& pattern);
}

#endif // WCS_ASTROMETRIC_STACKER_H
