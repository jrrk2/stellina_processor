// Updated WcsAstrometricStacker.cpp with TAN projection implementation
// Replace the WCS-related functions with these implementations

#include "WcsAstrometricStacker.h"
#include <QFileInfo>
#include <QTimer>
#include <QDateTime>
#include <QApplication>
#include <QTextStream>
#include <QRegularExpression>
#include <QDebug>
#include <cmath>
#include <algorithm>

// SOLUTION: List-based brightness compensation
// Each target pixel accumulates a list of contributing source pixels
// Final brightness = sum of all contributions / number of contributions

#include <vector>
#include <cmath>

void WCSAstrometricStacker::imageAccumWCS(size_t i)
      {
        updateProgress(10 + (i * 70) / m_images.size(), 
                       QString("Processing image %1 of %2: %3")
                       .arg(i+1)
                       .arg(m_images.size())
                       .arg(QFileInfo(m_images[i]->filename).fileName()));
        
        const auto& img = m_images[i];
        if (!img->wcs_valid || img->image.empty()) {
            logProcessing(QString("Skipping invalid image: %1")
                         .arg(QFileInfo(img->filename).fileName()));
            return;
        }
        
        // Calculate base image weight
        float imageWeight = img->quality_score;
        if (m_params.normalize_exposure && img->exposure_time > 0) {
            imageWeight *= img->exposure_time;
        }
        
        int pixelsProcessedThisImage = 0;
        
        // Project each source image onto the target grid
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                // Convert target pixel to world coordinates
                double ra, dec;
                if (m_output_wcs.pixelToWorld(x + 1.0, y + 1.0, ra, dec)) {
                    
                    // Convert world coordinates to source image coordinates
                    double imgpx, imgpy;
                    if (img->wcs.worldToPixel(ra, dec, imgpx, imgpy)) {
                        
                        float srcX = imgpx - 1.0f;
                        float srcY = imgpy - 1.0f;
                        
                        // Check bounds
                        if (srcX >= 0.0f && srcX < img->image.cols - 1 &&
                            srcY >= 0.0f && srcY < img->image.rows - 1) {
                            
                            // Bilinear interpolation
                            int x0 = int(floor(srcX));
                            int y0 = int(floor(srcY));
                            int x1 = x0 + 1;
                            int y1 = y0 + 1;
                            
                            if (x0 >= 0 && x1 < img->image.cols && 
                                y0 >= 0 && y1 < img->image.rows) {
                                
                                float dx = srcX - x0;
                                float dy = srcY - y0;
                                
                                float v00 = img->image.at<float>(y0, x0);
                                float v01 = img->image.at<float>(y0, x1);
                                float v10 = img->image.at<float>(y1, x0);
                                float v11 = img->image.at<float>(y1, x1);
                                
                                float v0 = v00 * (1 - dx) + v01 * dx;
                                float v1 = v10 * (1 - dx) + v11 * dx;
                                float pixelValue = v0 * (1 - dy) + v1 * dy;

				// Determine Bayer color for this pixel
                                uint8_t bayer_color = getBayerColor(x, y, img->bayer_pattern);

				
				
                                if (std::isfinite(pixelValue)) {
                                    // Add contribution to target pixel's list
                                    int targetIndex = y * m_output_size.width + x;
                                    m_pixel_lists[targetIndex].emplace_back(
                                        pixelValue, imageWeight, i, bayer_color
                                    );
                                    
                                    pixelsProcessedThisImage++;
                                    m_pixels_processed++;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        logProcessing(QString("Processed image %1: %2 pixels contributed")
                     .arg(i+1)
                     .arg(pixelsProcessedThisImage));
    }

// Replace the main stacking function with this list-based approach
bool WCSAstrometricStacker::beginStackWCSImages() {
    updateProgress(0, "Starting list-based astrometric stacking...");
    
    if (m_images.empty()) {
        emit errorOccurred("No images loaded for stacking");
        return false;
    }
    
    // Step 1: Compute optimal output WCS and dimensions
    updateProgress(5, "Computing optimal output coordinate system...");
    if (!computeOptimalWCS()) {
        emit errorOccurred("Failed to compute optimal WCS");
        return false;
    }
    
    // Initialize output images
    m_stacked_image_red = cv::Mat::zeros(m_output_size, CV_32F);
    m_stacked_image_green = cv::Mat::zeros(m_output_size, CV_32F);
    m_stacked_image_blue = cv::Mat::zeros(m_output_size, CV_32F);
    m_weight_map = cv::Mat::zeros(m_output_size, CV_32F);
    m_overlap_map = cv::Mat::zeros(m_output_size, CV_8U);
    
    // Step 2: Create pixel contribution lists for each target pixel
    updateProgress(10, "Initializing pixel contribution lists...");
    
    m_pixel_lists = std::vector<std::vector<PixelContribution>> (
        m_output_size.height * m_output_size.width
    );

 // Step 3: Accumulate contributions from all images (continues in pixelAccum)
    m_pixels_processed = 0;
    m_pixels_rejected = 0;
    
    return true;
}

bool WCSAstrometricStacker::pixelAccumWCS(int y)
    {
          bool useSigmaClipping = (m_params.rejection == StackingParams::SIGMA_CLIPPING);

        for (int x = 0; x < m_output_size.width; ++x) {
            int targetIndex = y * m_output_size.width + x;
            auto& contributions = m_pixel_lists[targetIndex];
            
            if (contributions.empty()) {
                // No contributions for this pixel
                continue;
            }
            
            // Store overlap count (number of contributing images)
            m_overlap_map.at<uchar>(y, x) = static_cast<uchar>(
                std::min(contributions.size(), size_t(255))
            );
            
            if (useSigmaClipping && contributions.size() > 2) {
                // Apply sigma clipping to the contribution list
                applySigmaClipping(contributions, x, y);
            }
            
            if (!contributions.empty()) {
                // Calculate final pixel value
                float redSum = 0.0f;
                float blueSum = 0.0f;
                float greenSum = 0.0f;

                for (const auto& contrib : contributions) {
                    switch(contrib.colour)
                    {
                        case RED:
                            redSum += contrib.value;
                            break;
                        case GREEN1:
                            greenSum += contrib.value;
                            break;
                        case GREEN2:
                            greenSum += contrib.value;
                            break;
                        case BLUE:
                            blueSum += contrib.value;
                            break;
                        default:
                            break;
                    }
                 }
                
                if (contributions.size() > 0) {
                    // KEY INSIGHT: Divide by contribution count for brightness compensation
                    m_stacked_image_red.at<float>(y, x) = redSum / contributions.size();
                    m_stacked_image_green.at<float>(y, x) = greenSum / (2 * contributions.size());
                    m_stacked_image_blue.at<float>(y, x) = blueSum / contributions.size();
                    m_weight_map.at<float>(y, x) = contributions.size();
                } else {
                    // Should not happen, but handle gracefully
                    m_stacked_image_red.at<float>(y, x) = 0.0f;
                    m_stacked_image_green.at<float>(y, x) = 0.0f;
                    m_stacked_image_blue.at<float>(y, x) = 0.0f;
                    m_weight_map.at<float>(y, x) = 0.0f;
                }
            }
        }
    return true;
    }

bool WCSAstrometricStacker::endStackWCSImages() {
    // Step 4: Process each target pixel's contribution list
    updateProgress(80, "Computing final pixel values...");
    
    // Step 5: Apply optional additional brightness normalization
    updateProgress(90, "Applying final brightness normalization...");
    
    if (m_params.apply_brightness_normalization) {
        applyGlobalBrightnessNormalization();
    }
    
    // Step 6: Generate statistics
    updateProgress(95, "Computing final statistics...");
    
    cv::Mat mask = (m_weight_map > 0);
    double minVal, maxVal;
    cv::Scalar meanVal = cv::mean(m_stacked_image_green, mask);
    cv::minMaxLoc(m_stacked_image_green, &minVal, &maxVal, nullptr, nullptr, mask);
    int validPixels = cv::countNonZero(mask);
    
    logProcessing(QString("List-based stacking complete:"));
    logProcessing(QString("  Valid pixels: %1 of %2 (%3%)")
                 .arg(validPixels)
                 .arg(m_output_size.area())
                 .arg(validPixels * 100.0 / m_output_size.area(), 0, 'f', 1));
    logProcessing(QString("  Pixel values: min=%1, max=%2, mean=%3")
                 .arg(minVal, 0, 'f', 3)
                 .arg(maxVal, 0, 'f', 3)
                 .arg(meanVal[0], 0, 'f', 3));
    
    if (m_pixels_rejected > 0) {
        logProcessing(QString("  Pixels rejected: %1 of %2 (%3%)")
                     .arg(m_pixels_rejected)
                     .arg(m_pixels_processed)
                     .arg(m_pixels_rejected * 100.0 / m_pixels_processed, 0, 'f', 1));
    }
    
    // Log overlap statistics
    logOverlapStatistics();
    
    updateProgress(100, "List-based stacking completed successfully");
    emit stackingComplete(true);
    
    return true;
}

// Helper function to apply sigma clipping to a contribution list
void WCSAstrometricStacker::applySigmaClipping(
    std::vector<PixelContribution>& contributions, 
    int x, int y) {
    
    if (contributions.size() <= 2) return;
    
    // Calculate mean and standard deviation
    double sum = 0.0, sumSq = 0.0;
    for (const auto& contrib : contributions) {
        sum += contrib.value;
        sumSq += contrib.value * contrib.value;
    }
    
    double mean = sum / contributions.size();
    double variance = (sumSq / contributions.size()) - (mean * mean);
    double stddev = sqrt(std::max(0.0, variance));
    
    // Apply sigma clipping thresholds
    float lowThreshold = mean - m_params.sigma_low * stddev;
    float highThreshold = mean + m_params.sigma_high * stddev;
    
    // Remove outliers
    auto it = std::remove_if(contributions.begin(), contributions.end(),
        [lowThreshold, highThreshold, this](const PixelContribution& contrib) {
            bool reject = (contrib.value < lowThreshold || contrib.value > highThreshold);
            if (reject) {
                m_pixels_rejected++;
            }
            return reject;
        }
    );
    
    contributions.erase(it, contributions.end());
}

// Helper function to apply global brightness normalization
void WCSAstrometricStacker::applyGlobalBrightnessNormalization() {
    // Optional: Apply gentle global correction based on overlap patterns
    // This is usually not needed with the list-based approach
    
    // Find regions with different overlap counts
    std::vector<float> meansByOverlap(256, 0.0f);
    std::vector<int> countsByOverlap(256, 0);
    
    for (int y = 0; y < m_output_size.height; ++y) {
        for (int x = 0; x < m_output_size.width; ++x) {
            uchar overlap = m_overlap_map.at<uchar>(y, x);
            if (overlap > 0) {
                meansByOverlap[overlap] += m_stacked_image_green.at<float>(y, x);
                countsByOverlap[overlap]++;
            }
        }
    }
    
    // Calculate average brightness for each overlap level
    for (int i = 1; i < 256; ++i) {
        if (countsByOverlap[i] > 0) {
            meansByOverlap[i] /= countsByOverlap[i];
        }
    }
    
    // Find the most common overlap level as reference
    int refOverlap = 1;
    int maxCount = 0;
    for (int i = 1; i < 256; ++i) {
        if (countsByOverlap[i] > maxCount) {
            maxCount = countsByOverlap[i];
            refOverlap = i;
        }
    }
    
    if (maxCount > 100) { // Only apply if we have enough reference pixels
        float refBrightness = meansByOverlap[refOverlap];
        
        logProcessing(QString("Global brightness reference: overlap=%1, brightness=%2")
                     .arg(refOverlap)
                     .arg(refBrightness, 0, 'f', 3));
        
        // Apply very gentle corrections only for significant differences
        for (int y = 0; y < m_output_size.height; ++y) {
            for (int x = 0; x < m_output_size.width; ++x) {
                uchar overlap = m_overlap_map.at<uchar>(y, x);
                if (overlap > 0 && countsByOverlap[overlap] > 50) {
                    float expectedBrightness = meansByOverlap[overlap];
                    if (expectedBrightness > 0.0f) {
                        float ratio = refBrightness / expectedBrightness;
                        
                        // Only apply gentle corrections for significant differences
                        if (ratio < 0.9f || ratio > 1.1f) {
                            float correction = 1.0f + 0.1f * (ratio - 1.0f); // 10% of the difference
                            m_stacked_image_red.at<float>(y, x) *= correction;
                            m_stacked_image_green.at<float>(y, x) *= correction;
                            m_stacked_image_blue.at<float>(y, x) *= correction;
                        }
                    }
                }
            }
        }
    }
}

// Helper function to log overlap statistics
void WCSAstrometricStacker::logOverlapStatistics() {
    std::vector<int> overlapCounts(256, 0);
    
    for (int y = 0; y < m_output_size.height; ++y) {
        for (int x = 0; x < m_output_size.width; ++x) {
            uchar overlap = m_overlap_map.at<uchar>(y, x);
            overlapCounts[overlap]++;
        }
    }
    
    logProcessing("=== OVERLAP STATISTICS ===");
    int totalPixels = m_output_size.area();
    
    for (int i = 1; i < 256; ++i) {
        if (overlapCounts[i] > 0) {
            logProcessing(QString("  %1 image(s): %2 pixels (%3%)")
                         .arg(i)
                         .arg(overlapCounts[i])
                         .arg(overlapCounts[i] * 100.0 / totalPixels, 0, 'f', 1));
        }
    }
    
    logProcessing(QString("  No coverage: %1 pixels (%2%)")
                 .arg(overlapCounts[0])
                 .arg(overlapCounts[0] * 100.0 / totalPixels, 0, 'f', 1));
}

// SimpleTANWCS implementation
bool SimpleTANWCS::pixelToWorld(double px, double py, double& ra, double& dec) const {
    if (!valid) return false;
    
    // Convert to 0-indexed and apply CD matrix
    double dx = px - crpix1;
    double dy = py - crpix2;
    
    double xi = cd11 * dx + cd12 * dy;
    double eta = cd21 * dx + cd22 * dy;
    
    // TAN projection inverse (standard formulae)
    double ra0_rad = crval1 * M_PI / 180.0;
    double dec0_rad = crval2 * M_PI / 180.0;
    double xi_rad = xi * M_PI / 180.0;
    double eta_rad = eta * M_PI / 180.0;
    
    double cos_dec0 = cos(dec0_rad);
    double sin_dec0 = sin(dec0_rad);
    
    double denom = cos_dec0 - eta_rad * sin_dec0;
    if (std::abs(denom) < 1e-12) return false;
    
    double ra_rad = ra0_rad + atan2(xi_rad, denom);
    double dec_rad = atan((sin_dec0 + eta_rad * cos_dec0) / (sqrt(xi_rad*xi_rad + denom*denom)));
    
    ra = ra_rad * 180.0 / M_PI;
    dec = dec_rad * 180.0 / M_PI;
    
    // Normalize RA
    while (ra < 0) ra += 360.0;
    while (ra >= 360.0) ra -= 360.0;
    
    return true;
}

bool SimpleTANWCS::worldToPixel(double ra, double dec, double& px, double& py) const {
    if (!valid) return false;
    
    // TAN projection forward
    double ra_rad = ra * M_PI / 180.0;
    double dec_rad = dec * M_PI / 180.0;
    double ra0_rad = crval1 * M_PI / 180.0;
    double dec0_rad = crval2 * M_PI / 180.0;
    
    double cos_dec = cos(dec_rad);
    double sin_dec = sin(dec_rad);
    double cos_dec0 = cos(dec0_rad);
    double sin_dec0 = sin(dec0_rad);
    double cos_dra = cos(ra_rad - ra0_rad);
    double sin_dra = sin(ra_rad - ra0_rad);
    
    double denom = sin_dec * sin_dec0 + cos_dec * cos_dec0 * cos_dra;
    if (std::abs(denom) < 1e-12) return false;
    
    double xi = cos_dec * sin_dra / denom;
    double eta = (sin_dec * cos_dec0 - cos_dec * sin_dec0 * cos_dra) / denom;
    
    // Convert to degrees
    xi *= 180.0 / M_PI;
    eta *= 180.0 / M_PI;
    
    // Apply inverse CD matrix
    double det = cd11 * cd22 - cd12 * cd21;
    if (std::abs(det) < 1e-12) return false;
    
    double dx = (cd22 * xi - cd12 * eta) / det;
    double dy = (-cd21 * xi + cd11 * eta) / det;
    
    px = dx + crpix1;
    py = dy + crpix2;
    
    return true;
}

void SimpleTANWCS::printDiagnostics() const {
    if (!valid) {
        qDebug() << "WCS is invalid";
        return;
    }
    
    qDebug() << "=== SimpleTAN WCS Diagnostics ===";
    qDebug() << "Reference point: RA =" << crval1 << "°, Dec =" << crval2 << "°";
    qDebug() << "Reference pixel: X =" << crpix1 << ", Y =" << crpix2;
    qDebug() << "CD matrix:";
    qDebug() << "  " << cd11 << cd12;
    qDebug() << "  " << cd21 << cd22;
    
    double det = cd11 * cd22 - cd12 * cd21;
    qDebug() << "Matrix determinant:" << det;
    qDebug() << "Pixel scale:" << getPixelScale() << "arcsec/pixel";
    
    // Test reference pixel
    double ra, dec;
    if (pixelToWorld(crpix1, crpix2, ra, dec)) {
        qDebug() << "Reference pixel maps to: RA =" << ra << "°, Dec =" << dec << "°";
        qDebug() << "Should equal CRVAL1/CRVAL2 - Error: RA =" 
                 << (ra - crval1) << "°, Dec =" << (dec - crval2) << "°";
    }
}

double SimpleTANWCS::getPixelScale() const {
    if (!valid) return 0.0;
    
    // Calculate pixel scale from CD matrix determinant
    double det = std::abs(cd11 * cd22 - cd12 * cd21);
    double scale_deg = sqrt(det);
    return scale_deg * 3600.0; // Convert to arcsec/pixel
}

// Updated WCSAstrometricStacker constructor
WCSAstrometricStacker::WCSAstrometricStacker(QObject *parent)
    : QObject(parent)
    , m_progress_bar(nullptr)
    , m_status_label(nullptr)
    , m_stacking_active(false)
    , m_current_image_index(0)
    , m_processing_timer(new QTimer(this))
    , m_total_processing_time(0.0)
    , m_pixels_processed(0)
    , m_stacked_image_red()
    , m_stacked_image_green()
    , m_stacked_image_blue()
    , m_pixels_rejected(0)
{
    // No WCSLIB initialization needed
    
    // Set up processing timer for non-blocking operation
    m_processing_timer->setSingleShot(true);
    //    connect(m_processing_timer, &QTimer::timeout, this, &WCSAstrometricStacker::processNextImage);
}

WCSAstrometricStacker::~WCSAstrometricStacker() {
    // No WCSLIB cleanup needed
}

// Updated loadWCSFromFITS function
bool WCSAstrometricStacker::loadWCSFromFITS(const QString &fits_file, WCSImageData &img_data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fits_file.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        emit errorOccurred(QString("Failed to open FITS file: %1 (status: %2)").arg(fits_file).arg(status));
        return false;
    }
    
    img_data.filename = fits_file;
    img_data.solved_filename = fits_file;
    
    // Read essential WCS keywords with error checking
    auto readKey = [&](const char* keyword, double& value) -> bool {
        int local_status = 0;
        return fits_read_key(fptr, TDOUBLE, keyword, &value, nullptr, &local_status) == 0;
    };
    
    char bayer[FLEN_VALUE];
    memset(bayer, 0, sizeof(bayer));
    bool success = !fits_read_key(fptr, TSTRING, "BAYERPAT", bayer, nullptr, &status);
    img_data.bayer_pattern = QString(bayer);
    
    success &= readKey("CRVAL1", img_data.wcs.crval1);
    success &= readKey("CRVAL2", img_data.wcs.crval2);
    success &= readKey("CRPIX1", img_data.wcs.crpix1);
    success &= readKey("CRPIX2", img_data.wcs.crpix2);
    
    // Try CD matrix first, fall back to CDELT
    if (readKey("CD1_1", img_data.wcs.cd11) && readKey("CD1_2", img_data.wcs.cd12) &&
        readKey("CD2_1", img_data.wcs.cd21) && readKey("CD2_2", img_data.wcs.cd22)) {
        // CD matrix available
        logProcessing("Using CD matrix for WCS");
    } else {
        // Use CDELT and assume no rotation
        double cdelt1, cdelt2;
        if (readKey("CDELT1", cdelt1) && readKey("CDELT2", cdelt2)) {
            img_data.wcs.cd11 = cdelt1; 
            img_data.wcs.cd12 = 0.0;
            img_data.wcs.cd21 = 0.0; 
            img_data.wcs.cd22 = cdelt2;
            logProcessing("Using CDELT for WCS (no rotation)");
        } else {
            success = false;
            logProcessing("No CD matrix or CDELT found");
        }
    }
    
    // Verify coordinate types are TAN projection
    char ctype1[FLEN_VALUE], ctype2[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "CTYPE1", ctype1, nullptr, &status) == 0 &&
        fits_read_key(fptr, TSTRING, "CTYPE2", ctype2, nullptr, &status) == 0) {
        
        QString ctype1_str = QString::fromLatin1(ctype1).trimmed().remove('\'').remove('"');
        QString ctype2_str = QString::fromLatin1(ctype2).trimmed().remove('\'').remove('"');
        
        if (!ctype1_str.contains("TAN") || !ctype2_str.contains("TAN")) {
            logProcessing(QString("Warning: Non-TAN projection detected: %1, %2").arg(ctype1_str).arg(ctype2_str));
            logProcessing("Proceeding anyway - results may be inaccurate for non-TAN projections");
        }
    }
    
    img_data.wcs.valid = success;
    img_data.wcs_valid = success;
    
    if (success) {
        logProcessing(QString("WCS loaded successfully from: %1").arg(QFileInfo(fits_file).fileName()));
        logProcessing(QString("  CRVAL: %1, %2").arg(img_data.wcs.crval1, 0, 'f', 6).arg(img_data.wcs.crval2, 0, 'f', 6));
        logProcessing(QString("  CRPIX: %1, %2").arg(img_data.wcs.crpix1, 0, 'f', 2).arg(img_data.wcs.crpix2, 0, 'f', 2));
        logProcessing(QString("  Pixel scale: %1 arcsec/pixel").arg(img_data.wcs.getPixelScale(), 0, 'f', 3));
        
        // Print full diagnostics in debug mode
        if (qEnvironmentVariableIsSet("DEBUG_WCS")) {
            img_data.wcs.printDiagnostics();
        }
    }
    
    // Read image dimensions and data
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, &status)) {
        emit errorOccurred(QString("Failed to get image size: %1 (status: %2)").arg(fits_file).arg(status));
        fits_close_file(fptr, &status);
        return false;
    }
    
    long totalPixels = naxes[0] * naxes[1];
    std::vector<float> pixels(totalPixels);
    
    if (fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, pixels.data(), nullptr, &status)) {
        emit errorOccurred(QString("Failed to read image data: %1 (status: %2)").arg(fits_file).arg(status));
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Convert to OpenCV Mat (FITS is row-major, OpenCV expects this)
    img_data.image = cv::Mat(naxes[1], naxes[0], CV_32F, pixels.data()).clone();
    
    fits_close_file(fptr, &status);
    
    logProcessing(QString("Successfully loaded WCS and image data: %1x%2 pixels")
                 .arg(naxes[0]).arg(naxes[1]));
    
    return success;
}

// Updated computeOptimalWCS function
bool WCSAstrometricStacker::computeOptimalWCS() {
    if (m_images.empty()) {
        emit errorOccurred("No images available for WCS calculation");
        return false;
    }
    
    logProcessing("Computing optimal output WCS using TAN projection...");
    
    // Calculate bounding box in world coordinates
    double ra_min = 360.0, ra_max = 0.0;
    double dec_min = 90.0, dec_max = -90.0;
    
    // Use the first valid WCS as reference
    const SimpleTANWCS* ref_wcs = nullptr;
    for (const auto& img : m_images) {
        if (img->wcs_valid) {
            ref_wcs = &img->wcs;
            break;
        }
    }
    
    if (!ref_wcs) {
        emit errorOccurred("No valid WCS found in any image");
        return false;
    }
    
    logProcessing(QString("Using reference WCS from first image: center RA=%1°, Dec=%2°")
                 .arg(ref_wcs->crval1, 0, 'f', 6).arg(ref_wcs->crval2, 0, 'f', 6));
    
    // For each image, find the corner points in world coordinates
    for (const auto& img : m_images) {
        if (!img->wcs_valid) continue;
        
        // Get image corners in pixel coordinates (1-indexed for FITS convention)
        double corners_pix[4][2] = {
            {1.0, 1.0},  // Bottom left
            {double(img->image.cols), 1.0},  // Bottom right
            {1.0, double(img->image.rows)},  // Top left
            {double(img->image.cols), double(img->image.rows)}  // Top right
        };
        
        // Convert each corner to world coordinates
        for (int i = 0; i < 4; ++i) {
            double ra, dec;
            if (img->wcs.pixelToWorld(corners_pix[i][0], corners_pix[i][1], ra, dec)) {
                // Handle RA wrapping at 0/360 degrees
                if (i == 0) {
                    ra_min = ra_max = ra;
                } else {
                    // Check if we cross the RA=0/360 boundary
                    if (std::abs(ra - ra_min) > 180.0) {
                        if (ra < ra_min) ra += 360.0;
                        else ra_min += 360.0;
                    }
                    if (std::abs(ra - ra_max) > 180.0) {
                        if (ra < ra_max) ra += 360.0;
                        else ra_max += 360.0;
                    }
                    
                    ra_min = std::min(ra_min, ra);
                    ra_max = std::max(ra_max, ra);
                }
                
                dec_min = std::min(dec_min, dec);
                dec_max = std::max(dec_max, dec);
            }
        }
    }
    
    // Normalize RA range
    while (ra_min >= 360.0) ra_min -= 360.0;
    while (ra_max >= 360.0) ra_max -= 360.0;
    
    // Add margin
    double ra_margin = (ra_max - ra_min) * 0.05;
    double dec_margin = (dec_max - dec_min) * 0.05;
    
    ra_min = std::max(0.0, ra_min - ra_margin);
    ra_max = std::min(360.0, ra_max + ra_margin);
    dec_min = std::max(-90.0, dec_min - dec_margin);
    dec_max = std::min(90.0, dec_max + dec_margin);
    
    logProcessing(QString("Field bounds: RA=%1° to %2°, Dec=%3° to %4°")
                 .arg(ra_min, 0, 'f', 4).arg(ra_max, 0, 'f', 4)
                 .arg(dec_min, 0, 'f', 4).arg(dec_max, 0, 'f', 4));
    
    // Determine pixel scale - use the first image as reference if not overridden
    if (m_params.output_pixel_scale <= 0.0) {
        m_output_pixel_scale = ref_wcs->getPixelScale();
    } else {
        m_output_pixel_scale = m_params.output_pixel_scale;
    }
    
    logProcessing(QString("Using output pixel scale: %1 arcsec/pixel").arg(m_output_pixel_scale, 0, 'f', 3));
    
    // Calculate output dimensions based on sky coverage and pixel scale
    double ra_span_deg = ra_max - ra_min;
    if (ra_span_deg < 0) ra_span_deg += 360.0;
    
    double dec_span_deg = dec_max - dec_min;
    
    // Convert to pixels using pixel scale
    int width, height;
    
    // Adjust RA span for cos(Dec) factor at center declination
    double dec_center = (dec_min + dec_max) / 2.0;
    double cos_dec_factor = std::cos(dec_center * M_PI / 180.0);
    if (cos_dec_factor < 0.01) cos_dec_factor = 0.01; // Avoid division by very small values
    
    width = int(ra_span_deg * 3600.0 / (m_output_pixel_scale / cos_dec_factor));
    height = int(dec_span_deg * 3600.0 / m_output_pixel_scale);
    
    // Override dimensions if specified
    if (m_params.output_width > 0) width = m_params.output_width;
    if (m_params.output_height > 0) height = m_params.output_height;
    
    // Limit to reasonable size
    width = std::min(std::max(width, 100), 10000);
    height = std::min(std::max(height, 100), 10000);
    
    m_output_size = cv::Size(width, height);
    
    logProcessing(QString("Output dimensions: %1 x %2 pixels").arg(width).arg(height));
    
    // Initialize output WCS as TAN projection centered on field
    double ra_center = (ra_min + ra_max) / 2.0;
    if (ra_max < ra_min) ra_center = fmod(ra_center + 180.0, 360.0); // Handle RA wrap
    
    m_output_wcs.crpix1 = width / 2.0 + 0.5;   // Center pixel X (1-indexed)
    m_output_wcs.crpix2 = height / 2.0 + 0.5;  // Center pixel Y (1-indexed)
    m_output_wcs.crval1 = ra_center;           // RA at reference pixel
    m_output_wcs.crval2 = dec_center;          // Dec at reference pixel
    m_output_wcs.cd11 = -m_output_pixel_scale / 3600.0; // RA step (degrees, negative for normal orientation)
    m_output_wcs.cd12 = 0.0;                            // No rotation
    m_output_wcs.cd21 = 0.0;                            // No rotation
    m_output_wcs.cd22 = m_output_pixel_scale / 3600.0;  // Dec step (degrees)
    m_output_wcs.valid = true;
    
    logProcessing(QString("Output WCS computed: center RA=%1°, Dec=%2°")
                 .arg(ra_center, 0, 'f', 6)
                 .arg(dec_center, 0, 'f', 6));
    
    // Test the output WCS
    double test_ra, test_dec;
    if (m_output_wcs.pixelToWorld(m_output_wcs.crpix1, m_output_wcs.crpix2, test_ra, test_dec)) {
        logProcessing(QString("WCS test: center pixel maps to RA=%1°, Dec=%2° (should match center)")
                     .arg(test_ra, 0, 'f', 6).arg(test_dec, 0, 'f', 6));
    }
    
    return true;
}

// Updated saveResult function

// Updated saveResult function to save RGB color planes
bool WCSAstrometricStacker::saveResult(const QString &output_path) {
    if (m_stacked_image_red.empty() || m_stacked_image_green.empty() || m_stacked_image_blue.empty()) {
        return false;
    }
    
    updateProgress(0, "Saving RGB TAN projection result...");
    
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = QString("!%1").arg(output_path).toLocal8Bit();
    if (fits_create_file(&fptr, pathBytes.data(), &status)) {
        return false;
    }
    
    // Create 3D FITS cube for RGB data: [width, height, 3]
    long naxes[3] = {m_stacked_image_green.cols, m_stacked_image_green.rows, 3};
    if (fits_create_img(fptr, FLOAT_IMG, 3, naxes, &status)) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Convert OpenCV Mats to FITS format
    long totalPixelsPerPlane = m_stacked_image_green.rows * m_stacked_image_green.cols;
    long totalPixels = totalPixelsPerPlane * 3;
    std::vector<float> pixels(totalPixels);
    
    // Interleave RGB data: R plane, G plane, B plane
    for (int y = 0; y < m_stacked_image_green.rows; ++y) {
        for (int x = 0; x < m_stacked_image_green.cols; ++x) {
            long baseIndex = y * m_stacked_image_green.cols + x;
            
            // Red plane (plane 0)
            pixels[baseIndex] = m_stacked_image_red.at<float>(y, x);
            
            // Green plane (plane 1)
            pixels[baseIndex + totalPixelsPerPlane] = m_stacked_image_green.at<float>(y, x);
            
            // Blue plane (plane 2)
            pixels[baseIndex + 2 * totalPixelsPerPlane] = m_stacked_image_blue.at<float>(y, x);
        }
    }
    
    if (fits_write_img(fptr, TFLOAT, 1, totalPixels, pixels.data(), &status)) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Write complete WCS keywords for TAN projection
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &m_output_wcs.crval1, "Reference RA (degrees)", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &m_output_wcs.crval2, "Reference Dec (degrees)", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &m_output_wcs.crpix1, "Reference pixel X", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &m_output_wcs.crpix2, "Reference pixel Y", &status);
    
    // Write CD matrix
    fits_write_key(fptr, TDOUBLE, "CD1_1", &m_output_wcs.cd11, "Coordinate matrix element", &status);
    fits_write_key(fptr, TDOUBLE, "CD1_2", &m_output_wcs.cd12, "Coordinate matrix element", &status);
    fits_write_key(fptr, TDOUBLE, "CD2_1", &m_output_wcs.cd21, "Coordinate matrix element", &status);
    fits_write_key(fptr, TDOUBLE, "CD2_2", &m_output_wcs.cd22, "Coordinate matrix element", &status);
    
    // Write coordinate types
    char ctype1[] = "RA---TAN";
    char ctype2[] = "DEC--TAN";
    char* ctype1_ptr = ctype1;
    char* ctype2_ptr = ctype2;
    fits_write_key(fptr, TSTRING, "CTYPE1", &ctype1_ptr, "Coordinate type", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", &ctype2_ptr, "Coordinate type", &status);
    
    // Write coordinate units
    char cunit1[] = "deg";
    char cunit2[] = "deg";
    char* cunit1_ptr = cunit1;
    char* cunit2_ptr = cunit2;
    fits_write_key(fptr, TSTRING, "CUNIT1", &cunit1_ptr, "Coordinate unit", &status);
    fits_write_key(fptr, TSTRING, "CUNIT2", &cunit2_ptr, "Coordinate unit", &status);
    
    // Add RGB-specific keywords
    char colorspace[] = "RGB";
    char* colorspace_ptr = colorspace;
    fits_write_key(fptr, TSTRING, "COLORSP", &colorspace_ptr, "Color space", &status);
    
    // Add plane descriptions
    char plane1[] = "Red";
    char plane2[] = "Green"; 
    char plane3[] = "Blue";
    char* plane1_ptr = plane1;
    char* plane2_ptr = plane2;
    char* plane3_ptr = plane3;
    fits_write_key(fptr, TSTRING, "PLANE1", &plane1_ptr, "First plane content", &status);
    fits_write_key(fptr, TSTRING, "PLANE2", &plane2_ptr, "Second plane content", &status);
    fits_write_key(fptr, TSTRING, "PLANE3", &plane3_ptr, "Third plane content", &status);
    
    // Add processing information
    int nimages = m_images.size();
    fits_write_key(fptr, TINT, "NSTACKED", &nimages, "Number of stacked images", &status);
    
    double pixelScale = m_output_wcs.getPixelScale();
    fits_write_key(fptr, TDOUBLE, "PIXSCALE", &pixelScale, "Pixel scale (arcsec/pixel)", &status);
    
    // Add processing method
    QString method = QString("TAN_PROJECTION_RGB_%1").arg(static_cast<int>(m_params.combination));
    QByteArray methodBytes = method.toLocal8Bit();
    char* methodPtr = methodBytes.data();
    fits_write_key(fptr, TSTRING, "STACKMET", &methodPtr, "Stacking method", &status);
    
    // Add history
    QString history = QString("Stacked %1 images using TAN projection astrometric alignment with RGB output").arg(nimages);
    QByteArray historyBytes = history.toLocal8Bit();
    fits_write_history(fptr, historyBytes.data(), &status);
    
    fits_close_file(fptr, &status);
    
    updateProgress(100, "RGB save complete");
    
    return (status == 0);
}

// Missing member function implementations

void WCSAstrometricStacker::setProgressWidgets(QProgressBar *progress, QLabel *status) {
    m_progress_bar = progress;
    m_status_label = status;
}

bool WCSAstrometricStacker::addImage(const QString &fits_file, const QString &solved_fits_file) {
    // Use the solved_fits_file if provided, otherwise use the fits_file
    QString fileToUse = solved_fits_file.isEmpty() ? fits_file : solved_fits_file;
    return addPlatesolveDFITSFile(fileToUse);
}

bool WCSAstrometricStacker::addImageWithMetadata(const QString &fits_file, const StellinaImageData &stellina_data) {
    auto img_data = std::make_unique<WCSImageData>();
    
    // Load FITS image and WCS from plate-solved file
    if (!loadWCSFromFITS(fits_file, *img_data)) {
        emit errorOccurred(QString("Failed to load WCS from: %1").arg(fits_file));
        return false;
    }
    
    // Integrate Stellina-specific data
    img_data->exposure_time = stellina_data.exposureSeconds;
    img_data->obs_time = QDateTime::fromString(stellina_data.dateObs, "yyyy-MM-ddThh:mm:ss");
    if (!img_data->obs_time.isValid()) {
        img_data->obs_time = QDateTime::fromString(stellina_data.dateObs, "yyyy-MM-ddThh:mm:ss.zzz");
    }
    
    // Use Stellina quality metrics if available
    if (stellina_data.hasValidCoordinates) {
        img_data->stellina_correction_magnitude = sqrt(pow(stellina_data.altitude, 2) +
                                                      pow(stellina_data.azimuth, 2));
    }
    
    // Extract image statistics and compute quality score
    extractImageStatistics(*img_data);
    computeImageQualityScore(*img_data);
    
    m_images.push_back(std::move(img_data));
    
    emit imageProcessed(fits_file);
    logProcessing(QString("Added Stellina image: %1 (quality: %2)")
                 .arg(QFileInfo(fits_file).fileName())
                 .arg(m_images.back()->quality_score, 0, 'f', 3));
    
    return true;
}

void WCSAstrometricStacker::setStackingParameters(const StackingParams &params) {
    m_params = params;
    logProcessing(QString("Updated stacking parameters: method=%1, rejection=%2")
                 .arg(static_cast<int>(params.combination))
                 .arg(static_cast<int>(params.rejection)));
}

bool WCSAstrometricStacker::addPlatesolveDFITSFile(const QString &solved_fits_file) {
    auto img_data = std::make_unique<WCSImageData>();
    
    // Load FITS image and WCS from plate-solved file
    if (!loadWCSFromFITS(solved_fits_file, *img_data)) {
        emit errorOccurred(QString("Failed to load WCS from: %1").arg(solved_fits_file));
        return false;
    }
    
    // Extract basic metadata from FITS headers
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = solved_fits_file.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status) == 0) {
        // Try to read exposure time
        double exptime = 10.0;
        if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &exptime, nullptr, &status) == 0) {
            img_data->exposure_time = exptime;
        }
        status = 0;
        
        // Try to read observation time
        char dateobs[FLEN_VALUE];
        if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) == 0) {
            QString dateStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
            img_data->obs_time = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss");
            if (!img_data->obs_time.isValid()) {
                // Try alternative formats
                img_data->obs_time = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss.zzz");
            }
        }
        
        fits_close_file(fptr, &status);
    }
    
    // Extract image statistics and compute quality score
    extractImageStatistics(*img_data);
    computeImageQualityScore(*img_data);
    
    m_images.push_back(std::move(img_data));
    
    emit imageProcessed(solved_fits_file);
    logProcessing(QString("Added FITS file: %1 (quality: %2)")
                 .arg(QFileInfo(solved_fits_file).fileName())
                 .arg(m_images.back()->quality_score, 0, 'f', 3));
    
    return true;
}

bool WCSAstrometricStacker::extractImageStatistics(WCSImageData &img_data) {
    if (img_data.image.empty()) return false;
    
    // Calculate basic statistics
    cv::Scalar mean, stddev;
    cv::meanStdDev(img_data.image, mean, stddev);
    
    img_data.background_level = mean[0];
    img_data.noise_level = stddev[0];
    
    // Estimate star count by counting bright pixels
    cv::Mat mask;
    double threshold = img_data.background_level + 3 * img_data.noise_level;
    cv::threshold(img_data.image, mask, threshold, 255, cv::THRESH_BINARY);
    
    // Convert to 8-bit for contour detection
    cv::Mat mask8;
    mask.convertTo(mask8, CV_8U);
    
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(mask8, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    
    img_data.star_count = 0;
    for (const auto &contour : contours) {
        double area = cv::contourArea(contour);
        if (area > 10 && area < 1000) {  // Reasonable star size range
            img_data.star_count++;
        }
    }
    
    return true;
}

bool WCSAstrometricStacker::computeImageQualityScore(WCSImageData &img_data) {
    double quality = 1.0;
    
    // Factor in noise level (lower noise = higher quality)
    if (img_data.noise_level > 0) {
        double snr = img_data.background_level / img_data.noise_level;
        quality *= std::min(1.0, snr / 50.0);  // Normalize to reasonable SNR range
    }
    
    // Factor in star count (more stars = better registration)
    if (img_data.star_count > 0) {
        quality *= std::min(1.0, img_data.star_count / 200.0);
    }
    
    // Factor in Stellina-specific metrics
    if (img_data.stellina_stars_used > 0) {
        quality *= std::min(1.0, img_data.stellina_stars_used / 100.0);
    }
    
    // Reduce quality for images with large correction magnitudes
    if (img_data.stellina_correction_magnitude > 100.0) {
        quality *= 0.5;  // Penalize images that needed large corrections
    }
    
    img_data.quality_score = std::max(0.1, std::min(1.0, quality));  // Clamp to [0.1, 1.0]
    
    return true;
}

/*
void WCSAstrometricStacker::startStacking() {
    if (m_images.empty()) {
        emit errorOccurred("No images loaded for stacking");
        return;
    }
    
    m_stacking_active = true;
    m_current_image_index = 0;
    
    updateProgress(0, "Starting TAN projection stacking...");
    
    if (!beginStackImages()) {
        emit errorOccurred("Stacking process failed");
        return;
    }
    
    analyzeImageQuality();
    emit qualityAnalysisComplete();
}
*/

void WCSAstrometricStacker::cancelStacking() {
    m_stacking_active = false;
    m_processing_timer->stop();
    emit statusUpdated("Stacking cancelled by user");
}

void WCSAstrometricStacker::analyzeImageQuality() {
    // Quality analysis is already done during image loading
    // This could be expanded to do more detailed analysis
    logProcessing("Image quality analysis complete");
    
    if (m_images.empty()) return;
    
    // Calculate aggregate statistics
    double totalQuality = 0.0;
    double totalExposure = 0.0;
    int totalStars = 0;
    
    for (const auto& img : m_images) {
        totalQuality += img->quality_score;
        totalExposure += img->exposure_time;
        totalStars += img->star_count;
    }
    
    double avgQuality = totalQuality / m_images.size();
    double avgStars = double(totalStars) / m_images.size();
    
    logProcessing(QString("Quality analysis: %1 images, avg quality: %2, avg stars: %3, total exposure: %4s")
                 .arg(m_images.size())
                 .arg(avgQuality, 0, 'f', 3)
                 .arg(avgStars, 0, 'f', 1)
                 .arg(totalExposure, 0, 'f', 1));
}

QString WCSAstrometricStacker::getQualityReport() const {
    QStringList report;
    
    report << "=== TAN Projection Astrometric Stacking Quality Report ===";
    report << "";
    report << QString("Total images processed: %1").arg(m_images.size());
    
    if (!m_stacked_image_green.empty()) {
        report << QString("Output dimensions: %1 x %2 pixels").arg(m_output_size.width).arg(m_output_size.height);
    }
    
    if (m_output_wcs.valid) {
        report << QString("Output pixel scale: %1 arcsec/pixel").arg(m_output_wcs.getPixelScale(), 0, 'f', 2);
        report << QString("Field center: RA=%1°, Dec=%2°").arg(m_output_wcs.crval1, 0, 'f', 6).arg(m_output_wcs.crval2, 0, 'f', 6);
    }
    
    report << QString("Total exposure time: %1 seconds").arg(getTotalExposureTime(), 0, 'f', 1);
    report << QString("Average image quality: %1").arg(getAverageQuality(), 0, 'f', 3);
    
    if (m_pixels_rejected > 0) {
        report << QString("Pixels rejected by sigma clipping: %1").arg(m_pixels_rejected);
    }
    
    report << "";
    report << "Individual Image Quality:";
    
    for (size_t i = 0; i < m_images.size(); ++i) {
        const auto &img = m_images[i];
        report << QString("  %1: Quality=%2, Stars=%3, Exposure=%4s")
                     .arg(QFileInfo(img->filename).fileName())
                     .arg(img->quality_score, 0, 'f', 3)
                     .arg(img->star_count)
                     .arg(img->exposure_time, 0, 'f', 1);
    }
    
    report << "";
    report << "Processing Log:";
    for (const QString &entry : m_processing_log) {
        report << QString("  %1").arg(entry);
    }
    
    return report.join("\n");
}

double WCSAstrometricStacker::getTotalExposureTime() const {
    double total = 0.0;
    for (const auto &img : m_images) {
        total += img->exposure_time;
    }
    return total;
}

double WCSAstrometricStacker::getAverageQuality() const {
    if (m_images.empty()) return 0.0;
    
    double total = 0.0;
    for (const auto &img : m_images) {
        total += img->quality_score;
    }
    return total / m_images.size();
}

void WCSAstrometricStacker::updateProgress(int percentage, const QString &message) {
    if (m_progress_bar) {
        m_progress_bar->setValue(percentage);
    }
    if (m_status_label) {
        m_status_label->setText(message);
    }
    
    emit progressUpdated(percentage);
    emit statusUpdated(message);
    logProcessing(message);
}

void WCSAstrometricStacker::logProcessing(const QString &message) {
    QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
    m_processing_log.append(QString("[%1] %2").arg(timestamp).arg(message));
}

void WCSAstrometricStacker::finishStacking() {
    m_stacking_active = false;
    emit stackingComplete(true);
    emit statusUpdated("TAN projection stacking completed successfully");
}

// Add this diagnostic function to WCSAstrometricStacker class
// Call it after the first pass to understand overlap patterns

void WCSAstrometricStacker::analyzeOverlapDistribution(const cv::Mat& overlap_count) {
    // Calculate overlap statistics
    double minOverlap, maxOverlap;
    cv::minMaxLoc(overlap_count, &minOverlap, &maxOverlap);
    
    cv::Scalar meanOverlap = cv::mean(overlap_count, overlap_count > 0);
    
    // Count pixels by overlap level
    int noOverlap = cv::countNonZero(overlap_count == 0);
    int singleOverlap = cv::countNonZero(overlap_count == 1);
    int doubleOverlap = cv::countNonZero(overlap_count == 2);
    int tripleOverlap = cv::countNonZero(overlap_count >= 3);
    int totalPixels = m_output_size.area();
    
    logProcessing("=== OVERLAP ANALYSIS ===");
    logProcessing(QString("Overlap range: %1 to %2 images per pixel")
                 .arg(minOverlap, 0, 'f', 1)
                 .arg(maxOverlap, 0, 'f', 1));
    logProcessing(QString("Mean overlap: %1 images per pixel")
                 .arg(meanOverlap[0], 0, 'f', 2));
    
    logProcessing(QString("Pixel distribution:"));
    logProcessing(QString("  No coverage: %1 pixels (%2%)")
                 .arg(noOverlap)
                 .arg(noOverlap * 100.0 / totalPixels, 0, 'f', 1));
    logProcessing(QString("  Single image: %1 pixels (%2%)")
                 .arg(singleOverlap)
                 .arg(singleOverlap * 100.0 / totalPixels, 0, 'f', 1));
    logProcessing(QString("  Two images: %1 pixels (%2%)")
                 .arg(doubleOverlap)
                 .arg(doubleOverlap * 100.0 / totalPixels, 0, 'f', 1));
    logProcessing(QString("  Three+ images: %1 pixels (%2%)")
                 .arg(tripleOverlap)
                 .arg(tripleOverlap * 100.0 / totalPixels, 0, 'f', 1));
    
    // This helps identify if your brightness problem correlates with overlap
    if (singleOverlap > totalPixels * 0.3) {
        logProcessing("WARNING: High percentage of single-overlap pixels may cause brightness variations");
    }
    
    if (tripleOverlap > totalPixels * 0.5) {
        logProcessing("INFO: Good overlap coverage - applying compensation for center over-brightening");
    }
}

// Also add this method to save overlap map for visual inspection
void WCSAstrometricStacker::saveOverlapMap(const cv::Mat& overlap_count, const QString& output_path) {
    // Convert overlap count to 8-bit for visualization
    cv::Mat overlap_vis;
    overlap_count.convertTo(overlap_vis, CV_8U, 255.0 / 5.0); // Scale assuming max 5 overlaps
    
    // Save as FITS for quantitative analysis
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = QString("!%1").arg("_overlap.fits").toLocal8Bit();
    if (fits_create_file(&fptr, pathBytes.data(), &status) == 0) {
        long naxes[2] = {overlap_count.cols, overlap_count.rows};
        if (fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status) == 0) {
            
            std::vector<float> pixels(overlap_count.rows * overlap_count.cols);
            for (int y = 0; y < overlap_count.rows; ++y) {
                for (int x = 0; x < overlap_count.cols; ++x) {
                    pixels[y * overlap_count.cols + x] = overlap_count.at<float>(y, x);
                }
            }
            
            if (fits_write_img(fptr, TFLOAT, 1, pixels.size(), pixels.data(), &status) == 0) {
                logProcessing(QString("Saved overlap map"));
            }
        }
        fits_close_file(fptr, &status);
    }
}

 // Utility function to determine Bayer color from pixel coordinates
uint8_t WCSAstrometricStacker::getBayerColor(int x, int y, const QString& pattern) {
    if (pattern == "RGGB") {
        if (y % 2 == 0) {
            return (x % 2 == 0) ? RED : GREEN1;  // R or G1
        } else {
            return (x % 2 == 0) ? GREEN2 : BLUE;  // G2 or B
        }
    } else if (pattern == "BGGR") {
        if (y % 2 == 0) {
            return (x % 2 == 0) ? BLUE : GREEN1;  // B or G1
        } else {
            return (x % 2 == 0) ? GREEN2 : RED;  // G2 or R
        }
    }
    return UNKNOWN; // Default to Green
}

