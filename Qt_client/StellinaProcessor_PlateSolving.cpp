// StellinaProcessor_PlateSolving.cpp - Plate solving operations module
// Extracted from StellinaProcessor_Core.cpp

#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
#include <QDir>
#include <QProcess>
#include <QThread>
#include <QFileInfo>
#include <QApplication>
#include <fitsio.h>

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
        if (hasStellinaMetadata(fullPath)) {
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

// ============================================================================
// Plate Solving Core Operations
// ============================================================================

bool StellinaProcessor::runSolveField(const QString &inputPath, const QString &outputPath, double ra, double dec) {
    // Convert coordinates to solve-field format
    QString raStr = CoordinateUtils::formatRaAsHMS(ra);
    QString decStr = CoordinateUtils::formatDecAsDMS(dec);
    
    // Build solve-field command
    QStringList arguments;
    arguments << inputPath;
    arguments << "--ra" << QString::number(ra, 'f', 6);
    arguments << "--dec" << QString::number(dec, 'f', 6);
    arguments << "--radius" << "2.0";  // 2 degree search radius
    arguments << "--scale-units" << "arcsecperpix";
    arguments << "--scale-low" << "1.0";
    arguments << "--scale-high" << "2.0";
    arguments << "--pixel-error" << "2";
    arguments << "--no-plots";
    arguments << "--overwrite";
    arguments << "--new-fits" << outputPath;
    
    // Add timeout and other solve-field options
    arguments << "--cpulimit" << "60";  // 60 second timeout
    arguments << "--no-verify";
    arguments << "--crpix-center";
    
    if (m_debugMode) {
        logMessage(QString("solve-field command: %1").arg(arguments.join(" ")), "gray");
    }
    
    // Run solve-field
    QProcess solveProcess;
    solveProcess.setProgram("solve-field");
    solveProcess.setArguments(arguments);
    
    // Set timeout for the process
    QTimer::singleShot(90000, &solveProcess, &QProcess::kill); // 90 second total timeout
    
    solveProcess.start();
    
    if (!solveProcess.waitForStarted(5000)) {
        logMessage("Failed to start solve-field process", "red");
        return false;
    }
    
    // Allow UI updates while solve-field runs
    while (solveProcess.state() == QProcess::Running) {
        QApplication::processEvents();
        QThread::msleep(100);
        
        if (!m_processing) {
            solveProcess.kill();
            return false; // User cancelled
        }
    }
    
    if (!solveProcess.waitForFinished(5000)) {
        logMessage("solve-field process timeout", "red");
        solveProcess.kill();
        return false;
    }
    
    int exitCode = solveProcess.exitCode();
    QString output = solveProcess.readAllStandardOutput();
    QString errorOutput = solveProcess.readAllStandardError();
    
    if (m_debugMode && !output.isEmpty()) {
        logMessage(QString("solve-field output: %1").arg(output.trimmed()), "gray");
    }
    
    if (exitCode != 0) {
        logMessage(QString("solve-field failed with exit code %1").arg(exitCode), "red");
        if (!errorOutput.isEmpty()) {
            logMessage(QString("Error: %1").arg(errorOutput.trimmed()), "red");
        }
        return false;
    }
    
    // Check if output file was created
    if (!QFile::exists(outputPath)) {
        logMessage("solve-field did not create output file", "red");
        return false;
    }
    
    return true;
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

bool StellinaProcessor::readSolveFieldResults(const QString &solvedPath, ProcessedImageData &data) {
    fitsfile *fptr = nullptr;
    int status = 0;
    
    QByteArray pathBytes = solvedPath.toLocal8Bit();
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
        data.pixelScale = sqrt(cd11*cd11 + cd12*cd12) * 3600.0; // arcsec/pixel
        data.hasValidWCS = true;
    }
    
    fits_close_file(fptr, &status);
    return data.hasValidWCS;
}

// ============================================================================
// Coordinate Calculation Functions
// ============================================================================

bool StellinaProcessor::calculateImageCoordinates(StellinaImageData &imageData) {
    if (!imageData.hasValidCoordinates) {
        return false;
    }
    
    // Parse observation time
    QDateTime obsDateTime = QDateTime::fromString(imageData.dateObs, Qt::ISODate);
    if (!obsDateTime.isValid()) {
        // Try alternative format
        obsDateTime = QDateTime::fromString(imageData.dateObs, "yyyy-MM-ddThh:mm:ss.zzz");
        if (!obsDateTime.isValid()) {
            logMessage(QString("Invalid observation time format: %1").arg(imageData.dateObs), "orange");
            return false;
        }
    }
    
    // Convert to UTC if needed
    if (obsDateTime.timeSpec() != Qt::UTC) {
        obsDateTime = obsDateTime.toUTC();
    }
    
    // Extract date/time components
    QDate obsDate = obsDateTime.date();
    QTime obsTime = obsDateTime.time();
    
    int year = obsDate.year();
    int month = obsDate.month();
    int day = obsDate.day();
    int hour = obsTime.hour();
    int minute = obsTime.minute();
    int second = obsTime.second();
    
    // Observer location (default to a reasonable location if not set)
    double latitude = 45.0;  // Default latitude
    double longitude = 2.0;  // Default longitude
    
    if (!m_observerLocation.isEmpty()) {
        // Parse observer location if provided (format: "lat,lon")
        QStringList coords = m_observerLocation.split(',');
        if (coords.size() == 2) {
            bool latOk, lonOk;
            double parsedLat = coords[0].trimmed().toDouble(&latOk);
            double parsedLon = coords[1].trimmed().toDouble(&lonOk);
            if (latOk && lonOk) {
                latitude = parsedLat;
                longitude = parsedLon;
            }
        }
    }
    
    // Convert Stellina alt/az to equatorial coordinates
    auto result = CoordinateUtils::calculateRaDec(
        year, month, day, hour, minute, second,
        imageData.altitude, imageData.azimuth,
        latitude, longitude
    );
    
    // Extract results from tuple
    double ra = std::get<0>(result);
    double dec = std::get<1>(result);
    
    // Store calculated coordinates
    imageData.calculatedRA = ra;
    imageData.calculatedDec = dec;
    imageData.hasCalculatedCoords = true;
    
    if (m_debugMode) {
        logMessage(QString("Calculated coordinates: Alt=%1°, Az=%2° -> RA=%3°, Dec=%4°")
                  .arg(imageData.altitude, 0, 'f', 2)
                  .arg(imageData.azimuth, 0, 'f', 2)
                  .arg(ra, 0, 'f', 4)
                  .arg(dec, 0, 'f', 4), "gray");
    }
    
    return true;
}

bool StellinaProcessor::convertStellinaToEquatorial(const StellinaImageData &stellinaData, double &ra, double &dec) {
    if (!stellinaData.hasCalculatedCoords) {
        return false;
    }
    
    ra = stellinaData.calculatedRA;
    dec = stellinaData.calculatedDec;
    return true;
}

// ============================================================================
// Plate Solving Utilities
// ============================================================================

QString StellinaProcessor::generatePlateSolvedPath(const QString &inputPath) {
    QString baseName = QFileInfo(inputPath).baseName();
    if (baseName.startsWith("dark_calibrated_")) {
        baseName = baseName.mid(16); // Remove prefix
    }
    
    QString outputName = QString("plate_solved_%1.fits").arg(baseName);
    return QDir(m_plateSolvedDirectory).absoluteFilePath(outputName);
}

bool StellinaProcessor::estimateFieldOfView(const QString &fitsPath, double &fovWidth, double &fovHeight) {
    // Default Stellina field of view
    fovWidth = 1.0;   // degrees
    fovHeight = 0.75; // degrees
    
    // Try to read image dimensions and calculate from pixel scale
    fitsfile *fptr = nullptr;
    int status = 0;
    long naxes[2];
    
    QByteArray pathBytes = fitsPath.toLocal8Bit();
    if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status) == 0) {
        int naxis;
        if (fits_get_img_dim(fptr, &naxis, &status) == 0 && naxis == 2) {
            if (fits_get_img_size(fptr, 2, naxes, &status) == 0) {
                // Calculate FOV based on Stellina's known pixel scale (~1.25 arcsec/pixel)
                double pixelScale = 1.25 / 3600.0; // degrees per pixel
                fovWidth = naxes[0] * pixelScale;
                fovHeight = naxes[1] * pixelScale;
            }
        }
        fits_close_file(fptr, &status);
    }
    
    return true;
}

QStringList StellinaProcessor::getSupportedSolveFieldOptions() {
    // Return list of solve-field options suitable for Stellina images
    QStringList options;
    
    options << "--scale-units" << "arcsecperpix";
    options << "--scale-low" << "1.0";     // Stellina pixel scale range
    options << "--scale-high" << "2.0";
    options << "--pixel-error" << "2";
    options << "--no-plots";
    options << "--overwrite";
    options << "--cpulimit" << "60";
    options << "--no-verify";
    options << "--crpix-center";
    
    return options;
}
