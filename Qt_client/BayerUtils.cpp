// BayerUtils.cpp - Utility functions for Bayer pattern handling
// Separated utilities for clean code organization

#include "WcsAstrometricStacker.h"

namespace BayerUtils {

uint8_t getBayerColor(int x, int y, const QString& pattern) {
    QString normalized = normalizeBayerPattern(pattern);
    
    if (normalized == "RGGB") {
        if (y % 2 == 0) {
            return (x % 2 == 0) ? RED : GREEN1;
        } else {
            return (x % 2 == 0) ? GREEN2 : BLUE;
        }
    } else if (normalized == "BGGR") {
        if (y % 2 == 0) {
            return (x % 2 == 0) ? BLUE : GREEN1;
        } else {
            return (x % 2 == 0) ? GREEN2 : RED;
        }
    } else if (normalized == "GRBG") {
        if (y % 2 == 0) {
            return (x % 2 == 0) ? GREEN1 : RED;
        } else {
            return (x % 2 == 0) ? BLUE : GREEN2;
        }
    } else if (normalized == "GBRG") {
        if (y % 2 == 0) {
            return (x % 2 == 0) ? GREEN1 : BLUE;
        } else {
            return (x % 2 == 0) ? RED : GREEN2;
        }
    }
    
    return GREEN1; // Default fallback
}

QString normalizeBayerPattern(const QString& pattern) {
    QString upper = pattern.toUpper().trimmed();
    
    // Handle common variations
    if (upper == "RGGB" || upper == "RG/GB") return "RGGB";
    if (upper == "BGGR" || upper == "BG/GR") return "BGGR";
    if (upper == "GRBG" || upper == "GR/BG") return "GRBG";
    if (upper == "GBRG" || upper == "GB/RG") return "GBRG";
    
    return "RGGB"; // Default
}

bool isValidBayerPattern(const QString& pattern) {
    QString normalized = normalizeBayerPattern(pattern);
    return (normalized == "RGGB" || normalized == "BGGR" || 
            normalized == "GRBG" || normalized == "GBRG");
}

} // namespace BayerUtils

// Updated integration with StellinaProcessor
// Add this method to StellinaProcessor class to use the new modular stacker

bool StellinaProcessor::performModularAstrometricStacking() {
    logMessage("Starting modular astrometric stacking with progress updates...", "blue");
    
    // Initialize the new modular stacker
    WCSAstrometricStacker stacker(this);
    
    // Set up progress tracking
    stacker.setProgressWidgets(m_progressBar, m_progressLabel);
    
    // Configure stacking parameters
    StackingParams params;
    params.combination = StackingParams::WEIGHTED_MEAN;
    params.rejection = StackingParams::SIGMA_CLIPPING;
    params.sigma_low = 3.0;
    params.sigma_high = 3.0;
    params.normalize_exposure = true;
    params.apply_brightness_normalization = true;
    params.create_weight_map = true;
    
    stacker.setStackingParameters(params);
    
    // Connect progress signals
    connect(&stacker, &WCSAstrometricStacker::progressUpdated,
            this, [this](int percentage) {
                m_progressBar->setValue(percentage);
            });
    
    connect(&stacker, &WCSAstrometricStacker::statusUpdated,
            this, [this](const QString& status) {
                logMessage(status, "blue");
            });
    
    // Add all plate-solved images to the stacker
    int images_added = 0;
    for (const QString& imagePath : m_imagesToProcess) {
        StellinaImageData* imageData = findImageDataByPath(imagePath);
        if (imageData && imageData->hasValidWCS) {
            QString solvedPath = imageData->currentFitsPath;
            
            if (stacker.addImageWithMetadata(solvedPath, *imageData)) {
                images_added++;
                logMessage(QString("Added image %1 to stacking queue").arg(
                    QFileInfo(imagePath).fileName()), "gray");
            } else {
                logMessage(QString("Failed to add image %1 to stacking").arg(
                    QFileInfo(imagePath).fileName()), "orange");
            }
        }
    }
    
    if (images_added < 2) {
        logMessage("Error: Need at least 2 images for stacking", "red");
        return false;
    }
    
    logMessage(QString("Added %1 images to modular stacker").arg(images_added), "blue");
    
    // Perform the stacking with automatic subframe processing
    if (!stacker.stackImages()) {
        logMessage("Modular stacking failed", "red");
        return false;
    }
    
    // Save the result
    QString outputDir = QFileInfo(m_imagesToProcess.first()).absolutePath();
    QString outputPath = QDir(outputDir).absoluteFilePath("stellina_stacked_modular.fits");
    
    if (!stacker.saveResult(outputPath)) {
        logMessage("Failed to save stacked result", "red");
        return false;
    }
    
    // Generate quality report
    QString qualityReport = stacker.getQualityReport();
    logMessage("=== MODULAR STACKING QUALITY REPORT ===", "blue");
    QStringList reportLines = qualityReport.split('\n');
    for (const QString& line : reportLines) {
        if (!line.trimmed().isEmpty()) {
            logMessage(line.trimmed(), "gray");
        }
    }
    
    logMessage(QString("Modular stacking completed successfully: %1").arg(outputPath), "green");
    
    // Update processing statistics
    double totalExposure = stacker.getTotalExposureTime();
    double avgQuality = stacker.getAverageQuality();
    
    logMessage(QString("Total exposure time: %1 seconds").arg(totalExposure, 0, 'f', 1), "blue");
    logMessage(QString("Average image quality: %1").arg(avgQuality, 0, 'f', 3), "blue");
    
    return true;
}

// Memory usage comparison logging
void WCSAstrometricStacker::logMemoryUsage() {
    size_t estimated_old_memory = m_output_size.area() * m_images.size() * sizeof(float) * 4; // Old list approach
    size_t estimated_new_memory = m_output_size.area() * sizeof(PixelAccumulator); // New accumulator approach
    
    double old_mb = estimated_old_memory / (1024.0 * 1024.0);
    double new_mb = estimated_new_memory / (1024.0 * 1024.0);
    double savings = ((old_mb - new_mb) / old_mb) * 100.0;
    
    logProcessing("=== MEMORY USAGE ANALYSIS ===");
    logProcessing(QString("Old list-based approach: ~%1 MB").arg(old_mb, 0, 'f', 1));
    logProcessing(QString("New accumulator approach: ~%1 MB").arg(new_mb, 0, 'f', 1));
    logProcessing(QString("Memory savings: %1%").arg(savings, 0, 'f', 1));
}

// Enhanced progress tracking for subframes
void WCSAstrometricStacker::updateSubframeProgress(
    int image_idx, int subframe_idx, int total_subframes) {
    
    int images_completed = image_idx;
    int total_images = m_images.size();
    
    // Calculate overall progress considering both images and subframes
    double image_progress = static_cast<double>(images_completed) / total_images;
    double subframe_progress = static_cast<double>(subframe_idx) / total_subframes / total_images;
    double overall_progress = (image_progress + subframe_progress) * 70.0; // 70% for stacking phase
    
    int percentage = 15 + static_cast<int>(overall_progress); // Start at 15% (after setup)
    
    QString message = QString("Processing image %1/%2, subframe %3/%4")
                     .arg(image_idx + 1)
                     .arg(total_images)
                     .arg(subframe_idx + 1)
                     .arg(total_subframes);
    
    updateProgress(percentage, message);
}