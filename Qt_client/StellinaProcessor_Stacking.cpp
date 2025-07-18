// StellinaProcessor_Stacking.cpp - Stacking operations module (symmetric)
// Extracted from StellinaProcessor_Core.cpp with symmetric processing pattern

#include "StellinaProcessor.h"
#include <QDir>
#include <QFileInfo>
#include <QApplication>

// ============================================================================
// Stacking Stage Functions (Symmetric with other stages)
// ============================================================================

bool StellinaProcessor::setupStackingStage() {
    m_currentStage = STAGE_STACKING;
    
    // Stacking uses plate-solved images
    QDir plateSolvedDir(m_plateSolvedDirectory);
    if (m_plateSolvedDirectory.isEmpty() || !plateSolvedDir.exists()) {
        logMessage("Stacking requires plate-solved images directory to be set", "red");
        return false;
    }
    
    QStringList solvedFiles = plateSolvedDir.entryList(
        QStringList() << "plate_solved_*.fits", QDir::Files);
    
    if (solvedFiles.isEmpty()) {
        logMessage("No plate-solved images found for stacking", "red");
        logMessage("Run plate solving first or select a directory with solved images", "red");
        return false;
    }
    
    // Build full paths and validate WCS
    m_imagesToProcess.clear();
    int validFiles = 0;
    
    for (const QString &solvedFile : solvedFiles) {
        QString fullPath = plateSolvedDir.absoluteFilePath(solvedFile);
        
        // Validate WCS and metadata
        if (validateSolveFieldResult(fullPath) && hasStellinaMetadata(fullPath)) {
            m_imagesToProcess.append(fullPath);
            validFiles++;
        } else {
            logMessage(QString("Skipping %1: invalid WCS or missing metadata").arg(solvedFile), "orange");
        }
    }
    
    if (validFiles < 2) {
        logMessage("Need at least 2 valid images for stacking", "red");
        return false;
    }
    
    logMessage(QString("Stacking stage: found %1 valid plate-solved images").arg(validFiles), "blue");
    
    // Initialize the WCS stacker
    m_stacker = std::make_unique<WCSAstrometricStacker>(this);
    m_stacker->setProgressWidgets(m_progressBar, m_progressLabel);
    
    // Configure stacking parameters
    StackingParams params;
    params.combination = StackingParams::WEIGHTED_MEAN;
    params.rejection = StackingParams::SIGMA_CLIPPING;
    params.sigma_low = 3.0;
    params.sigma_high = 3.0;
    params.normalize_exposure = true;
    params.apply_brightness_normalization = true;
    params.create_weight_map = true;
    
    m_stacker->setStackingParameters(params);
    
    // Add all valid images to the stacker
    int images_added = 0;
    for (const QString &imagePath : m_imagesToProcess) {
        StellinaImageData* imageData = findImageDataByPath(imagePath);
        if (imageData && imageData->hasValidWCS) {
            if (m_stacker->addImageWithMetadata(imagePath, *imageData)) {
                images_added++;
            }
        }
    }
    
    if (images_added < 2) {
        logMessage("Failed to add sufficient images to stacker", "red");
        return false;
    }
    
    // Compute optimal WCS and initialize data structures
    if (!m_stacker->computeOptimalWCS()) {
        logMessage("Failed to compute optimal WCS for stacking", "red");
        return false;
    }
    
    // Initialize pixel accumulators for symmetric processing
    cv::Size output_size = m_stacker->getOutputSize();
    m_pixel_accumulators.resize(output_size.area());
    
    m_stacking_initialized = true;
    m_stacking_subframe_row = 0;
    
    logMessage(QString("Stacking setup complete: %1 images, %2x%3 output")
              .arg(images_added)
              .arg(output_size.width)
              .arg(output_size.height), "blue");
    
    // Ensure stacked directory exists
    QDir stackedDir(m_stackedDirectory);
    if (!stackedDir.exists()) {
        if (!stackedDir.mkpath(".")) {
            logMessage("Failed to create stacked directory", "red");
            return false;
        }
    }
    
    return true;
}

bool StellinaProcessor::processImageStacking(const QString &currentFile) {
    if (!m_stacking_initialized || !m_stacker) {
        logMessage("Stacking not properly initialized", "red");
        return false;
    }
    
    QString filename = QFileInfo(currentFile).fileName();
    
    // Find the image data for this file
    const std::unique_ptr<WCSImageData> imageData;
    WCSAstrometricStacker *stack = new WCSAstrometricStacker;
    bool rslt = stack->loadWCSFromFITS(currentFile, *imageData);
    if (!rslt) {
        logMessage(QString("Skipping %1: No valid WCS data").arg(filename), "orange");
        return true; // Not an error, just skip
    }
    
    // Process one subframe of this image
    cv::Size output_size = m_stacker->getOutputSize();
    int total_subframes = (output_size.height + STACKING_SUBFRAME_HEIGHT - 1) / STACKING_SUBFRAME_HEIGHT;
    
    int start_row = m_stacking_subframe_row;
    int end_row = std::min(start_row + STACKING_SUBFRAME_HEIGHT, output_size.height);
    
    logMessage(QString("Processing %1, subframe rows %2-%3 (subframe %4/%5)")
              .arg(filename)
              .arg(start_row)
              .arg(end_row - 1)
              .arg((start_row / STACKING_SUBFRAME_HEIGHT) + 1)
              .arg(total_subframes), "gray");
    
    // Determine Bayer pattern for this image
    QString bayer_pattern = determineBayerPattern(currentFile);
    size_t img_idx = findImageIndex(currentFile);
    
    // Process this subframe for the current image
    bool success = m_stacker->processImageSubframe(
        imageData,
        m_pixel_accumulators,
        start_row,
        end_row,
        bayer_pattern,
        img_idx);
    
    if (!success) {
        logMessage(QString("Failed to process subframe for %1").arg(filename), "orange");
        return false;
    }
    
    // Move to next subframe row
    m_stacking_subframe_row += STACKING_SUBFRAME_HEIGHT;
    
    // If we've completed all subframes, reset for next image
    if (m_stacking_subframe_row >= output_size.height) {
        m_stacking_subframe_row = 0; // Reset for next image
        
        // Log completion of this image
        logMessage(QString("Completed processing all subframes for %1").arg(filename), "blue");
    }
    
    return true;
}

bool StellinaProcessor::finalizeStackingStage() {
    if (!m_stacker) {
        logMessage("No stacker to finalize", "red");
        return false;
    }
    
    logMessage("Finalizing stacking stage...", "blue");
    
    // Finalize the stacked image from accumulated pixel data
    if (!m_stacker->finalizeStackedImage(m_pixel_accumulators)) {
        logMessage("Failed to finalize stacked image", "red");
        return false;
    }
    
    // Apply post-processing
    m_stacker->applyBrightnessNormalization();
    m_stacker->generateQualityMaps();
    
    // Generate output filename
    QString outputName = QString("stellina_stacked_%1.fits")
                        .arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
    QString outputPath = QDir(m_stackedDirectory).absoluteFilePath(outputName);
    
    // Save the result
    if (!m_stacker->saveResult(outputPath)) {
        logMessage("Failed to save stacked result", "red");
        return false;
    }
    
    // Generate and log quality report
    QString qualityReport = m_stacker->getQualityReport();
    logMessage("=== STACKING QUALITY REPORT ===", "blue");
    QStringList reportLines = qualityReport.split('\n');
    for (const QString &line : reportLines) {
        if (!line.trimmed().isEmpty()) {
            logMessage(line.trimmed(), "gray");
        }
    }
    
    // Update final result tracking
    m_finalStackedImage = outputPath;
    
    logMessage(QString("Stacking completed successfully: %1").arg(outputPath), "green");
    
    // Update processing statistics
    double totalExposure = m_stacker->getTotalExposureTime();
    double avgQuality = m_stacker->getAverageQuality();
    
    logMessage(QString("Total exposure time: %1 seconds").arg(totalExposure, 0, 'f', 1), "blue");
    logMessage(QString("Average image quality: %1").arg(avgQuality, 0, 'f', 3), "blue");
    
    // Clean up
    m_pixel_accumulators.clear();
    m_stacker.reset();
    m_stacking_initialized = false;
    
    return true;
}

// ============================================================================
// Stacking Utilities (Symmetric Processing Support)
// ============================================================================

QString StellinaProcessor::determineBayerPattern(const QString &filename) {
    // First check FITS header
    QString pattern = getBayerPatternFromFits(filename);
    if (pattern != "RGGB") {
        return pattern; // Found non-default pattern in FITS
    }
    
    // Fallback to filename heuristics for Stellina images
    QString baseName = QFileInfo(filename).baseName().toLower();
    
    if (baseName.contains("r.fits") || baseName.endsWith("r")) {
        return "BGGR";  // Reversed images typically have rotated Bayer pattern
    }
    
    return "RGGB";  // Default for normal Stellina images
}

size_t StellinaProcessor::findImageIndex(const QString &filePath) {
    for (int i = 0; i < m_imagesToProcess.size(); ++i) {
        if (m_imagesToProcess[i] == filePath) {
            return static_cast<size_t>(i);
        }
    }
    return 0; // Default fallback
}

bool StellinaProcessor::initializePixelAccumulators() {
    if (!m_stacker) {
        return false;
    }
    
    cv::Size output_size = m_stacker->getOutputSize();
    m_pixel_accumulators.clear();
    m_pixel_accumulators.resize(output_size.area());
    
    logMessage(QString("Initialized %1 pixel accumulators").arg(output_size.area()), "gray");
    return true;
}

// ============================================================================
// Legacy Stacking Support (for backward compatibility)
// ============================================================================

bool StellinaProcessor::performLegacyStacking() {
    logMessage("Starting legacy stacking mode...", "blue");
    
    // This function provides backward compatibility with the old stacking method
    // Use the new modular stacker but run it in a single pass
    
    if (!setupStackingStage()) {
        return false;
    }
    
    // Process all images in one go (non-symmetric)
    if (!m_stacker->stackImages()) {
        logMessage("Legacy stacking failed", "red");
        return false;
    }
    
    // Save results
    QString outputPath = QDir(m_stackedDirectory).absoluteFilePath("stellina_stacked_legacy.fits");
    if (!m_stacker->saveResult(outputPath)) {
        logMessage("Failed to save legacy stacking result", "red");
        return false;
    }
    
    logMessage(QString("Legacy stacking completed: %1").arg(outputPath), "green");
    return true;
}

bool StellinaProcessor::performModularAstrometricStacking() {
    logMessage("Starting modular astrometric stacking with progress updates...", "blue");
    
    // This is the new recommended method that follows the symmetric pattern
    // It will be called from the main processing loop via setupStackingStage(),
    // processImageStacking(), and finalizeStackingStage()
    
    return setupStackingStage(); // Setup only - processing happens in main loop
}

// ============================================================================
// Stacking Analysis and Quality Control
// ============================================================================

void StellinaProcessor::analyzeStackingQuality() {
    if (!m_stacker) {
        logMessage("No stacker available for quality analysis", "orange");
        return;
    }
    
    // Get quality metrics from the stacker
    QString qualityReport = m_stacker->getQualityReport();
    double avgQuality = m_stacker->getAverageQuality();
    int imageCount = m_stacker->getImageCount();
    
    logMessage("=== STACKING QUALITY ANALYSIS ===", "blue");
    logMessage(QString("Images processed: %1").arg(imageCount), "blue");
    logMessage(QString("Average quality: %1").arg(avgQuality, 0, 'f', 3), "blue");
    
    // Analyze overlap patterns
    cv::Mat overlapMap = m_stacker->getOverlapMap();
    if (!overlapMap.empty()) {
        cv::Scalar meanOverlap = cv::mean(overlapMap);
        double minOverlap, maxOverlap;
        cv::minMaxLoc(overlapMap, &minOverlap, &maxOverlap);
        
        logMessage(QString("Overlap statistics: min=%1, mean=%2, max=%3")
                  .arg(minOverlap, 0, 'f', 1)
                  .arg(meanOverlap[0], 0, 'f', 1)
                  .arg(maxOverlap, 0, 'f', 1), "blue");
    }
    
    // Check for quality issues
    if (avgQuality < 0.5) {
        logMessage("Warning: Low average image quality detected", "orange");
    }
    
    if (imageCount < 10) {
        logMessage("Warning: Small number of images may affect final quality", "orange");
    }
}

bool StellinaProcessor::validateStackingInputs() {
    // Validate that we have sufficient images for stacking
    if (m_imagesToProcess.size() < 2) {
        logMessage("Need at least 2 images for stacking", "red");
        return false;
    }
    
    // Validate that all images have WCS information
    int validWCSCount = 0;
    for (const QString &imagePath : m_imagesToProcess) {
        if (validateSolveFieldResult(imagePath)) {
            validWCSCount++;
        }
    }
    
    if (validWCSCount < 2) {
        logMessage(QString("Only %1 images have valid WCS data").arg(validWCSCount), "red");
        return false;
    }
    
    // Validate output directory
    QDir stackedDir(m_stackedDirectory);
    if (!stackedDir.exists() && !stackedDir.mkpath(".")) {
        logMessage("Cannot create stacked output directory", "red");
        return false;
    }
    
    return true;
}

// ============================================================================
// Stacking Progress and Status
// ============================================================================

void StellinaProcessor::updateStackingProgress(int percentage, const QString &message) {
    // This function can be connected to the stacker's progress signals
    if (m_progressBar) {
        m_progressBar->setValue(percentage);
    }
    
    logMessage(message, "blue");
    
    // Allow UI updates
    QApplication::processEvents();
}

QString StellinaProcessor::getStackingStatusDescription() const {
    if (!m_stacker) {
        return "Stacking not initialized";
    }
    
    cv::Size outputSize = m_stacker->getOutputSize();
    int totalSubframes = (outputSize.height + STACKING_SUBFRAME_HEIGHT - 1) / STACKING_SUBFRAME_HEIGHT;
    int currentSubframe = (m_stacking_subframe_row / STACKING_SUBFRAME_HEIGHT) + 1;
    
    return QString("Stacking subframe %1/%2 (row %3/%4)")
           .arg(currentSubframe)
           .arg(totalSubframes)
           .arg(m_stacking_subframe_row)
           .arg(outputSize.height);
}
