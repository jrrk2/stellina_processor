// StellinaProcessor_DarkCalibration.cpp - Dark frame processing module
// Extracted from StellinaProcessor_Core.cpp

#include "StellinaProcessor.h"
#include <QDir>
#include <QFileInfo>
#include <QApplication>
#include <fitsio.h>
#include <opencv2/opencv.hpp>

// ============================================================================
// Dark Calibration Stage Functions (Symmetric with other stages)
// ============================================================================

bool StellinaProcessor::setupDarkCalibrationStage() {
    m_currentStage = STAGE_DARK_CALIBRATION;
    
    // Dark calibration always uses raw light frames from source directory
    if (!findStellinaImages()) {
        logMessage("No valid Stellina images found in source directory", "red");
        return false;
    }
    
    logMessage(QString("Dark calibration stage: found %1 raw images to process")
                  .arg(m_imagesToProcess.size()), "blue");
    
    if (m_darkFrames.isEmpty()) {
        logMessage("Warning: No dark frames available for calibration", "orange");
    } else {
        logMessage(QString("Using %1 dark frame groups for calibration")
                      .arg(m_darkFrames.size()), "blue");
    }
    
    return true;
}

bool StellinaProcessor::processImageDarkCalibration(const QString &currentFile) {
    QString filename = QFileInfo(currentFile).fileName();
    logMessage(QString("Dark calibrating: %1").arg(filename), "blue");
    
    // Read Stellina metadata from raw FITS file
    StellinaImageData imageData;
    imageData.originalFitsPath = currentFile;
    
    if (!readStellinaMetadataFromFits(currentFile, imageData)) {
        // Try to find corresponding JSON metadata
        QString jsonPath = currentFile;
        jsonPath.replace(".fits", ".json");
        
        if (QFile::exists(jsonPath)) {
            if (!readStellinaDataFromJSON(jsonPath, imageData)) {
                logMessage(QString("Could not read metadata for %1").arg(filename), "orange");
                return false;
            }
            imageData.originalJsonPath = jsonPath;
        } else {
            logMessage(QString("No metadata found for %1").arg(filename), "orange");
            return false;
        }
    }
    
    // Calculate equatorial coordinates from Stellina alt/az
    if (!calculateImageCoordinates(imageData)) {
        logMessage(QString("Failed to calculate coordinates for %1").arg(filename), "orange");
        return false;
    }
    
    // Find matching dark frame
    DarkFrame bestDark;
    if (!m_darkFrames.isEmpty() && matchDarkFrame(imageData, bestDark)) {
        // Apply dark calibration
        QString outputName = QString("dark_calibrated_%1").arg(filename);
        QString outputPath = QDir(m_calibratedDirectory).absoluteFilePath(outputName);
        
        if (applyMasterDark(currentFile, bestDark.filepath, outputPath)) {
            // Update image data with new path
            imageData.currentFitsPath = outputPath;
            imageData.processingStage = "DARK_CALIBRATED";
            
            // Write enhanced metadata to calibrated file
            if (writeStellinaMetadataWithCoordinates(outputPath, imageData)) {
                m_darkCalibratedFiles.append(outputPath);
                m_darkCalibratedCount++;
                
                // Store in our tracking list
                m_stellinaImageData.append(imageData);
                
                logMessage(QString("Dark calibration successful: %1").arg(outputName), "green");
                return true;
            } else {
                logMessage(QString("Failed to write metadata to %1").arg(outputName), "orange");
            }
        } else {
            logMessage(QString("Dark calibration failed for %1").arg(filename), "red");
        }
    } else {
        // No dark frame available - just copy with coordinate metadata
        QString outputName = QString("dark_calibrated_%1").arg(filename);
        QString outputPath = QDir(m_calibratedDirectory).absoluteFilePath(outputName);
        
        if (QFile::copy(currentFile, outputPath)) {
            imageData.currentFitsPath = outputPath;
            
            if (writeStellinaMetadataWithCoordinates(outputPath, imageData)) {
                m_darkCalibratedFiles.append(outputPath);
                m_stellinaImageData.append(imageData);
                
                logMessage(QString("Coordinates added (no dark applied): %1").arg(outputName), "blue");
                return true;
            }
        }
        
        logMessage(QString("Failed to process %1").arg(filename), "red");
    }
    
    return false;
}

bool StellinaProcessor::finalizeDarkCalibrationStage() {
    logMessage("Finalizing dark calibration stage...", "blue");
    
    if (m_darkCalibratedCount == 0) {
        logMessage("No images were successfully dark calibrated", "orange");
        return false;
    }
    
    logMessage(QString("Dark calibration complete: %1 images processed")
              .arg(m_darkCalibratedCount), "green");
    
    // Update processing statistics
    logMessage(QString("Dark calibrated files saved to: %1").arg(m_calibratedDirectory), "blue");
    
    return true;
}

// ============================================================================
// Dark Frame Management
// ============================================================================

bool StellinaProcessor::loadDarkFrames() {
    m_darkFrames.clear();
    
    if (m_darkDirectory.isEmpty()) {
        logMessage("No dark frame directory specified", "orange");
        return true; // Not an error, just no darks available
    }
    
    QDir darkDir(m_darkDirectory);
    if (!darkDir.exists()) {
        logMessage(QString("Dark frame directory does not exist: %1").arg(m_darkDirectory), "orange");
        return true;
    }
    
    // Look for dark frame files
    QStringList darkFiles = darkDir.entryList(
        QStringList() << "*.fits" << "master_dark_*.fits" << "dark_*.fits",
        QDir::Files);
    
    if (darkFiles.isEmpty()) {
        logMessage("No dark frames found in dark directory", "orange");
        return true;
    }
    
    // Process each dark frame
    for (const QString &darkFile : darkFiles) {
        QString darkPath = darkDir.absoluteFilePath(darkFile);
        DarkFrame dark;
        
        if (extractDarkFrameMetadata(darkPath, dark)) {
            m_darkFrames.append(dark);
            logMessage(QString("Loaded dark frame: %1 (exp=%2s, temp=%3°C)")
                      .arg(darkFile)
                      .arg(dark.exposure)
                      .arg(dark.temperature), "gray");
        }
    }
    
    logMessage(QString("Loaded %1 dark frames").arg(m_darkFrames.size()), "blue");
    return true;
}

bool StellinaProcessor::matchDarkFrame(const StellinaImageData &imageData, DarkFrame &bestMatch) {
    if (m_darkFrames.isEmpty()) {
        return false;
    }
    
    // Extract light frame parameters
    int lightExposure = imageData.exposureSeconds;
    int lightTemperature = imageData.temperatureKelvin;
    
    // Find best matching dark frame
    int bestScore = -1;
    int bestIndex = -1;
    
    for (int i = 0; i < m_darkFrames.size(); ++i) {
        const DarkFrame &dark = m_darkFrames[i];
        
        // Calculate matching score
        int exposureDiff = abs(dark.exposure - lightExposure);
        int tempDiff = abs(dark.temperature - lightTemperature);
        
        // Apply tolerances
        if (exposureDiff > m_exposureTolerance || tempDiff > m_temperatureTolerance) {
            continue; // Outside tolerance
        }
        
        // Score: lower is better
        int score = exposureDiff + tempDiff;
        
        if (bestIndex == -1 || score < bestScore) {
            bestScore = score;
            bestIndex = i;
        }
    }
    
    if (bestIndex >= 0) {
        bestMatch = m_darkFrames[bestIndex];
        
        if (m_debugMode) {
            logMessage(QString("Matched dark frame: exp=%1s, temp=%2°C (score=%3)")
                      .arg(bestMatch.exposure)
                      .arg(bestMatch.temperature)
                      .arg(bestScore), "gray");
        }
        
        return true;
    }
    
    return false;
}

bool StellinaProcessor::createMasterDark(const QList<DarkFrame> &darkFrames, const QString &outputPath) {
    if (darkFrames.isEmpty()) {
        logMessage("No dark frames provided for master creation", "red");
        return false;
    }
    
    logMessage(QString("Creating master dark from %1 frames").arg(darkFrames.size()), "blue");
    
    // Load first dark frame to get dimensions
    cv::Mat masterDark;
    bool firstFrame = true;
    int successfulFrames = 0;
    
    for (const DarkFrame &dark : darkFrames) {
        cv::Mat darkImage;
        if (!loadFITSImage(dark.filepath, darkImage)) {
            logMessage(QString("Failed to load dark frame: %1").arg(dark.filepath), "orange");
            continue;
        }
        
        // Convert to float for averaging
        cv::Mat darkFloat;
        darkImage.convertTo(darkFloat, CV_32F);
        
        if (firstFrame) {
            masterDark = darkFloat.clone();
            firstFrame = false;
        } else {
            // Check dimensions match
            if (darkFloat.size() != masterDark.size()) {
                logMessage(QString("Dark frame size mismatch: %1").arg(dark.filepath), "orange");
                continue;
            }
            
            // Accumulate
            masterDark += darkFloat;
        }
        
        successfulFrames++;
    }
    
    if (successfulFrames == 0) {
        logMessage("No dark frames could be loaded", "red");
        return false;
    }
    
    // Average the accumulated frames
    masterDark /= static_cast<float>(successfulFrames);
    
    // Convert back to 16-bit
    cv::Mat masterDark16;
    masterDark.convertTo(masterDark16, CV_16U);
    
    // Save the master dark
    if (!saveFITSImage(outputPath, masterDark16)) {
        logMessage("Failed to save master dark frame", "red");
        return false;
    }
    
    logMessage(QString("Master dark created successfully from %1 frames").arg(successfulFrames), "green");
    return true;
}

bool StellinaProcessor::applyMasterDark(const QString &lightFrame, const QString &masterDark, const QString &outputFrame) {
    logMessage(QString("Applying master dark to %1").arg(QFileInfo(lightFrame).fileName()), "blue");
    
    // Load images
    cv::Mat lightImage, darkImage;
    
    if (!loadFITSImage(lightFrame, lightImage)) {
        logMessage("Failed to load light frame", "red");
        return false;
    }
    
    if (!loadFITSImage(masterDark, darkImage)) {
        logMessage("Failed to load master dark", "red");
        return false;
    }
    
    // Check dimensions match
    if (lightImage.size() != darkImage.size()) {
        logMessage("Light and dark frame dimensions don't match", "red");
        return false;
    }
    
    // Subtract dark from light
    cv::Mat calibratedImage;
    if (lightImage.type() == CV_16U && darkImage.type() == CV_16U) {
        // 16-bit processing
        cv::subtract(lightImage, darkImage, calibratedImage, cv::noArray(), CV_16U);
    } else {
        // Convert to float for safer subtraction
        cv::Mat lightFloat, darkFloat;
        lightImage.convertTo(lightFloat, CV_32F);
        darkImage.convertTo(darkFloat, CV_32F);
        
        cv::Mat calibratedFloat;
        cv::subtract(lightFloat, darkFloat, calibratedFloat);
        
        // Clamp negative values to zero
        cv::threshold(calibratedFloat, calibratedFloat, 0, 0, cv::THRESH_TOZERO);
        
        // Convert back to original type
        calibratedFloat.convertTo(calibratedImage, lightImage.type());
    }
    
    // Save the calibrated image
    if (!saveFITSImage(outputFrame, calibratedImage)) {
        logMessage("Failed to save calibrated image", "red");
        return false;
    }
    
    // Copy FITS headers from original
    copyFITSHeaders(lightFrame, outputFrame);
    
    return true;
}

// ============================================================================
// Dark Frame Utilities
// ============================================================================

bool StellinaProcessor::extractDarkFrameMetadata(const QString &darkPath, DarkFrame &dark) {
    dark.filepath = darkPath;
    
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = darkPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Read exposure time
    if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &dark.exposure, nullptr, &status)) {
        // Try alternative keywords
        status = 0;
        if (fits_read_key(fptr, TDOUBLE, "EXPOSURE", &dark.exposure, nullptr, &status)) {
            dark.exposure = 10; // Default exposure
            status = 0;
        }
    }
    
    // Read temperature
    if (fits_read_key(fptr, TINT, "CCD-TEMP", &dark.temperature, nullptr, &status)) {
        status = 0;
        if (fits_read_key(fptr, TINT, "TEMP", &dark.temperature, nullptr, &status)) {
            dark.temperature = 15; // Default temperature
            status = 0;
        }
    }
    
    // Read binning
    char binning[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "BINNING", binning, nullptr, &status) == 0) {
        dark.binning = QString(binning);
    } else {
        dark.binning = "1x1";
        status = 0;
    }
    
    // Read Bayer pattern
    char pattern[FLEN_VALUE];
    if (fits_read_key(fptr, TSTRING, "BAYERPAT", pattern, nullptr, &status) == 0) {
        dark.bayerPattern = QString(pattern);
    } else {
        dark.bayerPattern = "RGGB";
        status = 0;
    }
    
    fits_close_file(fptr, &status);
    return true;
}

bool StellinaProcessor::validateDarkFrameCompatibility(const DarkFrame &dark, const StellinaImageData &light) {
    // Check exposure time tolerance
    int exposureDiff = abs(dark.exposure - light.exposureSeconds);
    if (exposureDiff > m_exposureTolerance) {
        return false;
    }
    
    // Check Bayer pattern compatibility
    if (!dark.bayerPattern.isEmpty() && !light.bayerPattern.isEmpty()) {
        if (dark.bayerPattern != light.bayerPattern) {
            return false;
        }
    }
    
    return true;
}

QString StellinaProcessor::generateMasterDarkPath(const DarkFrame &darkTemplate) {
    return QDir(m_darkDirectory).absoluteFilePath(
        QString("master_dark_%1s_%2C_%3.fits")
        .arg(darkTemplate.exposure)
        .arg(darkTemplate.temperature)
        .arg(darkTemplate.bayerPattern));
}