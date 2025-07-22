// StellinaProcessor_PlateSolving.cpp - Fixed version
// Key fixes:
// 1. Changed hasStellinaMetadata to readStellinaMetadataFromFits
// 2. Fixed hasValidWCS member name
// 3. Added pixelScale member to ProcessedImageData
// 4. Fixed convertStellinaToEquatorial function signature
// 5. Fixed undefined variables in runSolveField call
// 6. Fixed Qt signal/slot connection

#include <fitsio.h>
#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
#include <QDir>
#include <QProcess>
#include <QThread>
#include <QFileInfo>
#include <QApplication>

// ============================================================================
// Plate Solving Stage Functions (Symmetric with other stages)
// ============================================================================

bool StellinaProcessor::setupPlatesolvingStage() {
    m_currentStage = STAGE_PLATE_SOLVING;
    
    // Plate solving uses dark-calibrated images if available, otherwise fail
    QDir calibratedDir(m_calibratedDirectory);
    if (m_calibratedDirectory.isEmpty() || !calibratedDir.exists()) {
        logMessage("Plate solving requires calibrated images directory to be set", "red");
        return false;
    }
    
    QStringList calibratedFiles = calibratedDir.entryList(
        QStringList() << "dark_calibrated_*.fits" << "*calibrated*.fits", 
        QDir::Files);
    
    if (calibratedFiles.isEmpty()) {
        logMessage("No calibrated images found for plate solving", "red");
        logMessage("Run dark calibration first or select a directory with calibrated images", "red");
        return false;
    }
    
    // Build full paths and validate each file has coordinates
    m_imagesToProcess.clear();
    int validFiles = 0;
    
    for (const QString &calibratedFile : calibratedFiles) {
        QString fullPath = calibratedDir.absoluteFilePath(calibratedFile);
        
        // Check if file has valid Stellina metadata with coordinates
        // FIXED: Use readStellinaMetadataFromFits instead of hasStellinaMetadata
        StellinaImageData tempData;
        if (readStellinaMetadataFromFits(fullPath, tempData) && tempData.hasValidCoordinates) {
            m_imagesToProcess.append(fullPath);
            validFiles++;
        } else {
            logMessage(QString("Skipping %1: no coordinate metadata").arg(calibratedFile), "orange");
        }
    }
    
    if (validFiles == 0) {
        logMessage("No calibrated images with coordinate metadata found", "red");
        return false;
    }
    
    logMessage(QString("Plate solving stage: found %1 calibrated images with coordinates")
              .arg(validFiles), "blue");
    
    // Ensure plate solved directory exists
    QDir plateSolvedDir(m_plateSolvedDirectory);
    if (!plateSolvedDir.exists()) {
        if (!plateSolvedDir.mkpath(".")) {
            logMessage("Failed to create plate solved directory", "red");
            return false;
        }
    }
    
    return true;
}

bool StellinaProcessor::processImagePlatesolving(const QString &currentFile) {
    QString filename = QFileInfo(currentFile).fileName();
    logMessage(QString("Plate solving: %1").arg(filename), "blue");
    
    // Validate this is a calibrated file
    if (!currentFile.contains(m_calibratedDirectory)) {
        logMessage(QString("ERROR: Plate solving received non-calibrated file: %1").arg(filename), "red");
        return false;
    }
    
    // Read Stellina metadata from calibrated FITS file
    StellinaImageData imageData;
    if (!readStellinaMetadataFromFits(currentFile, imageData)) {
        logMessage(QString("No Stellina metadata in calibrated file: %1").arg(filename), "red");
        return false;
    }
    
    if (!imageData.hasValidCoordinates) {
        logMessage("No valid coordinates in calibrated file metadata", "red");
        return false;
    }
    
    // Use pre-calculated coordinates from calibrated file
    double ra, dec;
    if (imageData.hasPreCalculatedCoords()) {
        ra = imageData.calculatedRA;
        dec = imageData.calculatedDec;
        logMessage(QString("Using coordinates from calibrated file: RA=%1°, Dec=%2°")
                  .arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
    } else {
        logMessage("ERROR: Calibrated file missing pre-calculated coordinates", "red");
        return false;
    }
    
    // Create output file name
    QString baseName = QFileInfo(currentFile).baseName();
    if (baseName.startsWith("dark_calibrated_")) {
        baseName = baseName.mid(16); // Remove "dark_calibrated_" prefix
    }
    QString outputName = QString("plate_solved_%1.fits").arg(baseName);
    QString outputPath = QDir(m_plateSolvedDirectory).absoluteFilePath(outputName);
    
    // Perform plate solving
    logMessage(QString("Plate solving with solve-field: RA=%1°, Dec=%2°").arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
    
    // FIXED: Pass the required parameters to runSolveField
    if (!runSolveField(currentFile, outputPath, ra, dec)) {
        logMessage("Plate solving failed", "red");
        return false;
    }
    
    logMessage("Plate solving succeeded!", "green");
    
    // Write metadata to plate-solved file
    StellinaImageData plateSolvedImageData = imageData;
    plateSolvedImageData.currentFitsPath = outputPath;
    
    if (!writeStellinaMetadataWithCoordinates(outputPath, plateSolvedImageData)) {
        logMessage("Warning: Failed to write metadata to plate-solved file", "orange");
    } else {
        updateProcessingStage(outputPath, "PLATE_SOLVED");
    }
    
    m_plateSolvedFiles.append(outputPath);
    
    // Update our tracking data
    for (auto &data : m_stellinaImageData) {
        if (data.originalFitsPath == imageData.originalFitsPath) {
            data.currentFitsPath = outputPath;
            // FIXED: Use correct member name
            data.hasValidWCS = true;
            break;
        }
    }
    
    return true;
}

bool StellinaProcessor::finalizePlatesolvingStage() {
    logMessage("Finalizing plate solving stage...", "blue");
    
    if (m_plateSolvedFiles.isEmpty()) {
        logMessage("No images were successfully plate solved", "orange");
        return false;
    }
    
    logMessage(QString("Plate solving complete: %1 images solved")
              .arg(m_plateSolvedFiles.size()), "green");
    
    // Validate solve results
    int validSolves = 0;
    for (const QString &solvedPath : m_plateSolvedFiles) {
        if (validateSolveFieldResult(solvedPath)) {
            validSolves++;
        }
    }
    
    logMessage(QString("Valid plate solutions: %1/%2").arg(validSolves).arg(m_plateSolvedFiles.size()), 
               validSolves > 0 ? "green" : "orange");
    
    return validSolves > 0;
}

// Helper function to read solve-field results from WCS headers
bool StellinaProcessor::readSolveFieldResults(const QString &fitsPath, ProcessedImageData &data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Read WCS parameters
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &data.solvedRA, nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    status = 0;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &data.solvedDec, nullptr, &status) != 0) {
        fits_close_file(fptr, &status);
        return false;
    }
    
    // Read CD matrix for pixel scale
    double cd11, cd12, cd21, cd22;
    status = 0;
    if (fits_read_key(fptr, TDOUBLE, "CD1_1", &cd11, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CD1_2", &cd12, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CD2_1", &cd21, nullptr, &status) == 0 &&
        fits_read_key(fptr, TDOUBLE, "CD2_2", &cd22, nullptr, &status) == 0) {
        
        // Calculate pixel scale from CD matrix
        // FIXED: Add pixelScale member to ProcessedImageData
        data.pixelScale = sqrt(cd11*cd11 + cd12*cd12) * 3600.0; // arcsec/pixel
        data.hasValidWCS = true;
    }
    
    fits_close_file(fptr, &status);
    return data.hasValidWCS;
}

// FIXED: Match the function signature from the header file
bool StellinaProcessor::convertStellinaToEquatorial(const StellinaImageData &stellinaData, double &ra, double &dec) {
    if (!stellinaData.hasCalculatedCoords) {
        return false;
    }
    
    ra = stellinaData.calculatedRA;
    dec = stellinaData.calculatedDec;
    return true;
}

// Stellar Solver integration callbacks
void StellinaProcessor::onStellarSolverImageSolved(const QString& filename, double ra, double dec, double pixelScale) {
    if (m_stellarSolverManager) {
        m_processedCount++;
        logMessage(QString("Plate solved: %1 (RA=%2°, Dec=%3°, Scale=%4 arcsec/px)")
                   .arg(QFileInfo(filename).baseName())
                   .arg(ra, 0, 'f', 4)
                   .arg(dec, 0, 'f', 4)
                   .arg(pixelScale, 0, 'f', 2), "green");
    } else {
        m_errorCount++;
        logMessage(QString("Plate solve failed: %1").arg(QFileInfo(filename).baseName()), "red");
    }
    
    // Move to next image automatically
    m_currentImageIndex++;
    updateProcessingStatus();
    
    // Continue processing or finish
    if (m_currentImageIndex < m_imagesToProcess.length()) {
        QTimer::singleShot(10, this, &StellinaProcessor::processNextImage);
    } else {
        finishProcessing();
    }
}

// FIXED: Updated signal/slot connection to match parameter count
void StellinaProcessor::onStellarSolverImageSkipped(const QString& filename, const QString& reason) {
    m_skippedCount++;
    logMessage(QString("Skipped %1: %2").arg(QFileInfo(filename).baseName()).arg(reason), "orange");
}

void StellinaProcessor::onStellarSolverBatchComplete() {
    logMessage("StellarSolver batch processing complete", "green");
    finishProcessing();
}

void StellinaProcessor::onStellarSolverError(const QString& error) {
    logMessage(QString("StellarSolver error: %1").arg(error), "red");
    finishProcessing();
}

bool StellinaProcessor::validateSolveFieldResult(const QString &solvedPath) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = solvedPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
        return false;
    }
    
    // Check for required WCS keywords
    char value[FLEN_VALUE];
    bool hasWCS = (fits_read_key(fptr, TSTRING, "CRVAL1", value, nullptr, &status) == 0);
    status = 0;
    hasWCS = hasWCS && (fits_read_key(fptr, TSTRING, "CRVAL2", value, nullptr, &status) == 0);
    status = 0;
    hasWCS = hasWCS && (fits_read_key(fptr, TSTRING, "CD1_1", value, nullptr, &status) == 0);
    
    fits_close_file(fptr, &status);
    
    if (hasWCS && m_debugMode) {
        logMessage(QString("WCS validation passed for %1").arg(QFileInfo(solvedPath).fileName()), "gray");
    }
    
    return hasWCS;
}

void StellinaProcessor::onStellarSolverProgressUpdated(int current, int total, const QString& status) {
    m_progressBar->setValue(current);
    m_progressBar->setMaximum(total);
    logMessage(status, "blue");
}

void StellinaProcessor::onStellarSolverImageSolved(const QString& filename, bool success, double ra, double dec, double pixelScale) {
    if (success) {
        m_processedCount++;
        logMessage(QString("✓ Solved: %1 - RA: %2°, Dec: %3°, Scale: %4 arcsec/px")
                   .arg(QFileInfo(filename).baseName())
                   .arg(ra, 0, 'f', 4)
                   .arg(dec, 0, 'f', 4)
                   .arg(pixelScale, 0, 'f', 2), "green");
    } else {
        m_errorCount++;
        logMessage(QString("✗ Failed: %1").arg(QFileInfo(filename).baseName()), "red");
    }
    
    updateProcessingStatus();
}
