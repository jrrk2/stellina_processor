#include "StellinaProcessor.h"
#include "WcsAstrometricStacker.h"
#include <QApplication>
#include <QDir>
#include <QProcess>
#include <QFileInfo>
#include <QJsonObject>
#include <QJsonDocument>
#include <QDateTime>
#include <QMessageBox>
// Redesigned clean stage-based processing
// Replace the startStellinaProcessing function in StellinaProcessor_Processing.cpp

void StellinaProcessor::startStellinaProcessing() {
    if (!validateProcessingInputs()) {
        return;
    }
    
    logMessage(QString("Starting %1 processing...")
                  .arg(m_processingModeCombo->currentText()), "blue");
    
    // Initialize processing state
    m_processing = true;
    m_currentImageIndex = 0;
    m_processedCount = 0;
    m_errorCount = 0;
    m_skippedCount = 0;
    m_darkCalibratedCount = 0;
    m_registeredCount = 0;
    m_processingStartTime = QDateTime::currentMSecsSinceEpoch();
    
    // Clear previous results
    m_darkCalibratedFiles.clear();
    m_plateSolvedFiles.clear();
    m_registeredFiles.clear();
    m_finalStackedImage.clear();
    
    // Determine processing stage and input files based on mode
    bool success = false;
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        success = setupPlatesolvingStage();
        break;
    case MODE_DARK_CALIBRATION:
        success = setupDarkCalibrationStage();
        break;
    case MODE_ASTROMETRIC_STACKING:
        success = setupStackingStage();
        break;
    case MODE_FULL_PIPELINE:
        success = setupDarkCalibrationStage(); // Start with first stage
        break;
    }
    
    if (!success) {
        m_processing = false;
        return;
    }
    
    // Set up progress tracking
    m_progressBar->setMaximum(m_imagesToProcess.length());
    m_progressBar->setValue(0);
    
    updateProcessingStatus();
    updateUI();
    
    logMessage(QString("Ready to process %1 images in %2 stage")
                  .arg(m_imagesToProcess.length())
                  .arg(getStageDescription()), "green");
    
    // Start processing timer
    m_processingTimer->start();
}

// New setup functions for each stage
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
            m_stellarSolverManager->addJob(fullPath);
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
    // Start batch processing
    m_stellarSolverManager->setOutputDirectory(m_plateSolvedDirectory);
    m_stellarSolverManager->startBatchSolving();
    return true;
}

bool StellinaProcessor::setupIntegrationStage() {
    m_currentStage = STAGE_INTEGRATION;
    m_currentIntegrationRow = 0;
    m_currentImageIndex = 0;
    m_progressBar->setMaximum( m_wcsStacker->getOutputHeight() );
    m_progressBar->setValue(0);
    updateProcessingStatus();
    return true;
}

bool StellinaProcessor::setupStackingStage() {
        m_currentStage = STAGE_STACKING;
        
    QDir solvedDir(m_plateSolvedDirectory);
    if (!solvedDir.exists()) {
        logMessage(QString("Solved directory does not exist: '%1'").arg(m_sourceDirectory), "red");
        return false;
    }
    
    logMessage(QString("Scanning directory: %1").arg(solvedDir.absolutePath()), "blue");
    
    QStringList fitsFiles = solvedDir.entryList(QStringList() << "*.fits" << "*.fit" << "*.FITS" << "*.FIT", QDir::Files);
    logMessage(QString("Found %1 FITS files").arg(fitsFiles.count()), "blue");
    
    int loaded = 0;
    int qualityRejected = 0;
    int reversedImages = 0;
    // Build full paths
    m_imagesToProcess.clear();
    // Set up WCS stacker with current plate-solved files
    m_wcsStacker->setProgressWidgets(m_subTaskProgressBar, m_currentTaskLabel);
    
    for (const QString &fitsFile : fitsFiles) {
        StellinaImageData imageData;
        imageData.originalFitsPath = solvedDir.absoluteFilePath(fitsFile);
        imageData.currentFitsPath = imageData.originalFitsPath;
        
        // Check if this is a reversed stellina image
        imageData.isReversedImage = isReversedStellinaImage(fitsFile);
        if (imageData.isReversedImage) {
            reversedImages++;
            logMessage(QString("Found reversed Stellina image: %1").arg(fitsFile), "blue");
        }
        
        // Get base name (without 'r' suffix for reversed images)
        imageData.baseName = getBaseName(fitsFile);
        
        // Detect bayer pattern
        imageData.bayerPattern = detectBayerPattern(imageData.originalFitsPath);
        
        QString baseName = QFileInfo(fitsFile).baseName();
        
        imageData.hasValidCoordinates = true;
        
        // Extract FITS metadata
        imageData.exposureSeconds = extractExposureTime(imageData.originalFitsPath);
        imageData.temperatureKelvin = extractTemperature(imageData.originalFitsPath);
        imageData.binning = extractBinning(imageData.originalFitsPath);
        imageData.dateObs = extractDateObs(imageData.originalFitsPath);
        
        // Quality filtering
        if (m_qualityFilter && !checkStellinaQuality(imageData.metadata)) {
            qualityRejected++;
            if (m_debugMode) {
                logMessage(QString("Rejected due to invalid coordinates: %1").arg(fitsFile), "orange");
            }
            qualityRejected++;
            continue;
        }
        
        m_stellinaImageData.append(imageData);
        m_imagesToProcess.append(imageData.originalFitsPath);
        if (m_wcsStacker->addImageWithMetadata(imageData.originalFitsPath, imageData)) {
            loaded++;
        }
    }
    
    qDebug() << "image to process" << m_imagesToProcess.size();
    logMessage(QString("Astrometric stacking stage: found %1 plate-solved images")
               .arg(m_imagesToProcess.size()), "blue");
    
    if (m_imagesToProcess.size() < 3) {
        logMessage("Warning: Need at least 3 images for effective stacking", "orange");
    }
    
    if (loaded < 3) {
        logMessage(QString("Only loaded %1 images for WCS stacking").arg(loaded), "red");
        return false;
    }
    
    // Start WCS stacking
    if (!m_wcsStacker->beginStackWCSImages()) {
        logMessage("Begin WCS stacking failed", "red");
        return false;
    }

    
    return true;
}

bool StellinaProcessor::accumulateStacking() {
    
    // accumulate images
    for (size_t i = 0; i < m_wcsStacker->getImageCount(); ++i) m_wcsStacker->imageAccumWCS(i);

    // accumulate mapped pixels
    for (int y = 0; y < m_wcsStacker->getOutputHeight(); ++y) m_wcsStacker->pixelAccumWCS(y);
    // End WCS stacking
    if (!m_wcsStacker->endStackWCSImages()) {
        logMessage("End WCS stacking failed", "red");
        return false;
    }

    // Save result automatically
    QString outputName = QString("wcs_stacked_%1.fits")
                           .arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
    QString outputPath = QDir(getOutputDirectoryForCurrentStage()).absoluteFilePath(outputName);
    
    if (m_wcsStacker->saveResult(outputPath)) {
        m_finalStackedImage = outputPath;
        logMessage(QString("WCS stacking completed: %1").arg(outputName), "green");
        return true;
    } else {
        logMessage("Failed to save WCS stacked result", "red");
        return false;
    }
    
    return true;
}

void StellinaProcessor::processNextImage() {
    bool complete = m_currentStage == STAGE_INTEGRATION ? m_currentIntegrationRow >= m_wcsStacker->getOutputHeight() :
    m_currentImageIndex >= m_imagesToProcess.length();
    
    if (complete) {
        // Current stage complete
        if (m_processingMode == MODE_FULL_PIPELINE) {
            // Handle pipeline stage transitions
            handlePipelineStageTransition();
        } else if (m_processingMode == MODE_ASTROMETRIC_STACKING && m_currentStage == STAGE_STACKING) {
            setupIntegrationStage();
        } else {
            // Single stage complete
            finishProcessing();
        }
        return;
    }
    
    if (m_currentStage != STAGE_INTEGRATION)
        logMessage(QString("Processing %1 of %2: %3")
                   .arg(m_currentImageIndex)
                   .arg(m_imagesToProcess.length())
                   .arg(QFileInfo(m_imagesToProcess[m_currentImageIndex]).fileName()), "blue");
    
    bool success = false;
    qint64 elapsed;
    double avgTimePerImage, avgTimePerChunk;
    qint64 remaining;
    
    switch (m_currentStage) {
        case STAGE_DARK_CALIBRATION:
            success = processImageDarkCalibration();
            m_progressBar->setValue(++m_currentImageIndex);
            // Update time estimate
            elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
            avgTimePerImage = static_cast<double>(elapsed) / m_currentImageIndex;
            remaining = static_cast<qint64>((m_imagesToProcess.length() - m_currentImageIndex) * avgTimePerImage);
            m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                            .arg(formatProcessingTime(remaining)));
            break;
        case STAGE_PLATE_SOLVING:
            success = processImagePlatesolving();
            m_progressBar->setValue(m_currentImageIndex);
            // Update time estimate
            elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
            avgTimePerImage = static_cast<double>(elapsed) / m_currentImageIndex;
            remaining = static_cast<qint64>((m_imagesToProcess.length() - m_currentImageIndex) * avgTimePerImage);
            m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                            .arg(formatProcessingTime(remaining)));
            break;
        case STAGE_STACKING:
            success = processImageStacking();
            m_progressBar->setValue(++m_currentImageIndex);
            // Update time estimate
            elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
            avgTimePerImage = static_cast<double>(elapsed) / m_currentImageIndex;
            remaining = static_cast<qint64>((m_imagesToProcess.length() - m_currentImageIndex) * avgTimePerImage);
            
            m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                            .arg(formatProcessingTime(remaining)));
            break;
        case STAGE_INTEGRATION:
            success = processImageIntegration();
            m_progressBar->setValue(m_currentIntegrationRow);
            // Update time estimate
            elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
            avgTimePerChunk = static_cast<double>(elapsed) / m_currentIntegrationRow;
            remaining = static_cast<qint64>((m_wcsStacker->getOutputHeight() - m_currentIntegrationRow) * avgTimePerChunk);
            m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                            .arg(formatProcessingTime(remaining)));
            break;
        default:
            logMessage(QString("Unknown processing stage: %1").arg(m_currentStage), "red");
            success = false;
            break;
    }
    
    // Update counters
    if (success) {
        m_processedCount++;
    } else {
        m_errorCount++;
    }
    
    // Update progress
    updateProcessingStatus();
    
}

// Handle full pipeline stage transitions
void StellinaProcessor::handlePipelineStageTransition() {
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        logMessage("Dark calibration stage complete, setting up plate solving...", "blue");
        if (setupPlatesolvingStage()) {
            m_currentImageIndex = 0;
            m_progressBar->setMaximum(m_imagesToProcess.length());
            m_progressBar->setValue(0);
            updateProcessingStatus();
        } else {
            logMessage("Failed to set up plate solving stage", "red");
            finishProcessing();
        }
        break;
        
        case STAGE_PLATE_SOLVING:
            logMessage("Plate solving stage complete, setting up astrometric stacking...", "blue");
            if (setupStackingStage()) {
                m_currentImageIndex = 0;
                m_progressBar->setMaximum(m_imagesToProcess.length());
                m_progressBar->setValue(0);
                updateProcessingStatus();
            } else {
                logMessage("Failed to set up stacking stage", "red");
                finishProcessing();
            }
            break;
            
        case STAGE_STACKING:
            logMessage("Stacking complete, setting up image integration...", "blue");
            if (setupIntegrationStage()) {
            } else {
                logMessage("Failed to set up stacking stage", "red");
                finishProcessing();
            }
            break;

    default:
        finishProcessing();
        break;
    }
}

// Replace your solve-field processing with:
bool StellinaProcessor::processImagePlatesolving() {
    return true; // Return immediately, completion handled by signals
}

bool StellinaProcessor::validateProcessingInputs() {
    if (m_sourceDirectory.isEmpty()) {
        QMessageBox::warning(this, "Directory Error", "Please select the raw light frames directory.");
        return false;
    }
    
    // Check required directories based on processing mode
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        if (m_plateSolvedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", "Please select the plate-solved output directory.");
            return false;
        }
        break;
        
    case MODE_DARK_CALIBRATION:
        if (m_calibratedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", "Please select the calibrated lights output directory.");
            return false;
        }
        if (m_autoMatchDarks && m_darkFrames.isEmpty()) {
            int ret = QMessageBox::question(this, "No Dark Frames", 
                                           "Dark calibration is enabled but no dark frames were found. Continue without dark calibration?",
                                           QMessageBox::Yes | QMessageBox::No);
            if (ret == QMessageBox::No) {
                return false;
            }
        }
        break;
        
    case MODE_ASTROMETRIC_STACKING:
        if (m_stackedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", "Please select the stacked images output directory.");
            return false;
        }
        break;
        
    case MODE_FULL_PIPELINE:
        if (m_calibratedDirectory.isEmpty() || m_plateSolvedDirectory.isEmpty() || m_stackedDirectory.isEmpty()) {
            QMessageBox::warning(this, "Directory Error", 
                                 "Full pipeline requires all output directories:\n"
                                 "- Calibrated Lights\n"
                                 "- Plate-Solved Lights\n" 
                                 "- Stacked Images");
            return false;
        }
        if (m_autoMatchDarks && m_darkFrames.isEmpty()) {
            int ret = QMessageBox::question(this, "No Dark Frames", 
                                           "Dark calibration is enabled but no dark frames were found. Continue without dark calibration?",
                                           QMessageBox::Yes | QMessageBox::No);
            if (ret == QMessageBox::No) {
                return false;
            }
        }
        break;
    }
    
    return true;
}

void StellinaProcessor::updateProcessingStatus() {
    QString stageText = getStageDescription();
    
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1").arg(
            m_processing ? "In Progress" : "Ready"));
        if (m_processing) {
                m_progressLabel->setText(QString("Calibrating %1 of %2 images (Stage: %3) - Success: %4, Errors: %5")
                                            .arg(m_currentImageIndex)
                                            .arg(m_imagesToProcess.length())
                                            .arg(stageText)
                                            .arg(m_processedCount)
                                            .arg(m_errorCount));
            }
        break;
    case STAGE_PLATE_SOLVING:
        // Keep existing registration/stacking status as ready for now
        if (m_processing) {
                m_progressLabel->setText(QString("Solving %1 of %2 images (Stage: %3) - Success: %4, Errors: %5")
                                            .arg(m_currentImageIndex)
                                            .arg(m_imagesToProcess.length())
                                            .arg(stageText)
                                            .arg(m_processedCount)
                                            .arg(m_errorCount));
            }
        break;
    case STAGE_STACKING:
        m_stackingStatusLabel->setText("Stacking: In Progress");
        if (m_processing) {
                m_progressLabel->setText(QString("Stacking %1 of %2 images (Stage: %3) - Success: %4, Errors: %5")
                                            .arg(m_currentImageIndex)
                                            .arg(m_imagesToProcess.length())
                                            .arg(stageText)
                                            .arg(m_processedCount)
                                            .arg(m_errorCount));
            }
        break;
    case STAGE_INTEGRATION:
        m_registrationStatusLabel->setText("Integration: In Progress");
        if (m_processing) {
                m_progressLabel->setText(QString("Integrating %1 of %2 rows (Stage: %3) - Success: %4, Errors: %5")
                                            .arg(m_currentIntegrationRow)
                                            .arg(m_wcsStacker->getOutputHeight())
                                            .arg(stageText)
                                            .arg(m_processedCount)
                                            .arg(m_errorCount));
            }
        else m_progressLabel->setText(QString("Integrated"));
        break;
    case STAGE_COMPLETE:
        m_darkCalibrationStatusLabel->setText("Dark Calibration: Complete");
        m_registrationStatusLabel->setText("Registration: Complete");
        m_stackingStatusLabel->setText("Stacking: Complete");
        m_progressLabel->setText(QString("Complete"));
        break;
    }
    
    if (m_processing && m_darkCalibratedCount > 0) {
            m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1 completed").arg(m_darkCalibratedCount));
        }
}

QString StellinaProcessor::getStageDescription() const {
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION: return "Dark Calibration";
    case STAGE_PLATE_SOLVING: return "Plate Solving";
    case STAGE_INTEGRATION: return "Integration";
    case STAGE_STACKING: return "Stacking";
    case STAGE_COMPLETE: return "Complete";
    default: return "Unknown";
    }
}

void StellinaProcessor::finishProcessing() {
    m_processing = false;
    m_processingTimer->stop();
    
    qint64 totalTime = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    
    QString completionMessage = QString("Processing complete! Total time: %1\n\n")
                                   .arg(formatProcessingTime(totalTime));
    
    switch (m_processingMode) {
    case MODE_BASIC_PLATESOLVE:
        completionMessage += QString("Plate solved: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        break;
        
    case MODE_DARK_CALIBRATION:
        completionMessage += QString("Dark calibrated: %1, Skipped: %2, Errors: %3, Total: %4")
                                .arg(m_darkCalibratedCount)
                                .arg(m_skippedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        break;
        
    case MODE_ASTROMETRIC_STACKING:
        endStacking();
        completionMessage += QString("Images processed: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
        
    case MODE_FULL_PIPELINE:
        endStacking();
        completionMessage += QString("Dark calibrated: %1, Plate solved: %2, Errors: %3, Total: %4")
                                .arg(m_darkCalibratedCount)
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
    }
    
    logMessage(completionMessage, "green");
    m_progressBar->setMaximum(1);
    m_progressBar->setValue(1);

    // Save processing report
    saveProcessingReport();
    
    updateUI();
    
    // Show completion dialog
    if (false) QMessageBox::information(this, "Processing Complete", completionMessage);
}
