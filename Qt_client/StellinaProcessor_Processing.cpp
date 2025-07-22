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

// Simplified processNextImage that doesn't handle stage transitions
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
            m_progressBar->setValue(++m_currentImageIndex);
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
    const QString &calibratedFitsPath = m_imagesToProcess[m_currentImageIndex];
//    if (!m_stellarSolverManager) m_stellarSolverManager = new StellarSolverManager;
    // Use StellarSolver instead of solve-field
    if (m_currentImageIndex == 0) {
            // Initialize batch on first image
            if (!m_stellarSolverManager->initializeBatch(m_imagesToProcess)) {
                return false;
            }
        }
        
    // Start processing current image
    m_stellarSolverManager->processNextImage();
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

bool StellinaProcessor::checkSolveFieldInstalled() {
    QProcess process;
    process.start("/opt/homebrew/bin/solve-field", QStringList() << "--help");
    bool finished = process.waitForFinished(3000);
    return finished && process.exitCode() == 0;
}

QStringList StellinaProcessor::getAstrometryPaths() {
    QStringList astrometryPaths;
    
    // Common directories where astrometry.net tools are installed
    QStringList searchDirs = {
        "/opt/homebrew/bin",              // Homebrew Apple Silicon
        "/opt/homebrew/libexec/astrometry.net/bin",  // Homebrew support tools
        "/usr/local/bin",                 // Homebrew Intel / manual install
        "/usr/local/astrometry/bin",      // Custom install location
        "/usr/bin",                       // System package (Linux)
        "/opt/astrometry.net/bin",        // Alternative install
        "/usr/local/libexec/astrometry.net/bin"  // Support tools location
    };
    
    // Check which directories actually exist and contain astrometry tools
    for (const QString &dir : searchDirs) {
        QDir directory(dir);
        if (directory.exists()) {
            // Check for common astrometry.net tools
            QStringList requiredTools = {"solve-field", "fits2fits", "pnmfile", "image2xy"};
            bool hasTools = false;
            
            for (const QString &tool : requiredTools) {
                if (QFile::exists(directory.absoluteFilePath(tool))) {
                    hasTools = true;
                    break;
                }
            }
            
            if (hasTools) {
                astrometryPaths.append(dir);
                if (m_debugMode) {
                    logMessage(QString("Found astrometry tools in: %1").arg(dir), "gray");
                }
            }
        }
    }
    
    return astrometryPaths;
}

QProcessEnvironment StellinaProcessor::createSolveFieldEnvironment() {
    // Start with current environment
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    
    // Get astrometry.net tool directories
    QStringList astrometryPaths = getAstrometryPaths();
    
    // Get current PATH
    QString currentPath = env.value("PATH");
    QStringList pathList = currentPath.split(QDir::listSeparator(), Qt::SkipEmptyParts);
    
    // Add astrometry paths to the beginning (higher priority)
    for (const QString &astroPath : astrometryPaths) {
        if (!pathList.contains(astroPath)) {
            pathList.prepend(astroPath);
        }
    }
    
    // Set the enhanced PATH
    QString enhancedPath = pathList.join(QDir::listSeparator());
    env.insert("PATH", enhancedPath);
    
    if (m_debugMode) {
        logMessage(QString("Enhanced PATH for solve-field: %1").arg(enhancedPath), "gray");
    }
    
    // Add other useful environment variables for astrometry.net
    
    // Set data directory hints (if not already set)
    if (!env.contains("ASTROMETRY_NET_DATA_DIR")) {
        QStringList dataDirs = {
            "/opt/homebrew/share/astrometry",
            "/usr/local/share/astrometry",
            "/usr/share/astrometry",
            QDir::homePath() + "/.local/share/astrometry"
        };
        
        for (const QString &dataDir : dataDirs) {
            if (QDir(dataDir).exists()) {
                env.insert("ASTROMETRY_NET_DATA_DIR", dataDir);
                if (m_debugMode) {
                    logMessage(QString("Set ASTROMETRY_NET_DATA_DIR: %1").arg(dataDir), "gray");
                }
                break;
            }
        }
    }
    
    // Ensure solve-field can find its configuration
    if (!env.contains("SOLVE_FIELD_CONFIG")) {
        QStringList configPaths = {
            "/opt/homebrew/etc/astrometry.cfg",
            "/usr/local/etc/astrometry.cfg",
            "/etc/astrometry.cfg"
        };
        
        for (const QString &configPath : configPaths) {
            if (QFile::exists(configPath)) {
                env.insert("SOLVE_FIELD_CONFIG", configPath);
                if (m_debugMode) {
                    logMessage(QString("Found astrometry config: %1").arg(configPath), "gray");
                }
                break;
            }
        }
    }
    
    return env;
}

bool StellinaProcessor::runSolveField(const QString &fitsPath, const QString &outputPath, double ra, double dec) {
    logMessage(QString("Running solve-field with coordinate hints: RA=%1°, Dec=%2°").arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
    
    // Check if solve-field is available
    if (!checkSolveFieldInstalled()) {
        logMessage("solve-field not found - install astrometry.net", "red");
        return false;
    }
    
    QProcess solveProcess;
    QStringList arguments;
    
    // solve-field arguments with coordinate hints and binning
    arguments << fitsPath;
    arguments << "--overwrite";           // Overwrite existing files
    arguments << "--no-plots";            // Don't create plot files
    arguments << "--new-fits" << outputPath;  // Output solved FITS
    arguments << "--scale-units" << "arcsecperpix";
    arguments << "--scale-low" << "1.2";   // Original pixel scale range
    arguments << "--scale-high" << "1.3";  // solve-field will auto-adjust for downsampling
    arguments << "--ra" << QString::number(ra, 'f', 6);     // Use calculated RA
    arguments << "--dec" << QString::number(dec, 'f', 6);   // Use calculated Dec
    arguments << "--radius" << "10.0";      // Small search radius since we have good hints
    arguments << "--cpulimit" << "60";     // 1 minute timeout (faster with hints)
    arguments << "--no-verify";           // Skip verification step
    arguments << "--crpix-center";        // Set reference pixel at center
    arguments << "-z" << "2";             // Debayer the CFA image (CRITICAL for Stellina)
    
    if (m_debugMode) {
        logMessage(QString("solve-field command: solve-field %1").arg(arguments.join(" ")), "gray");
    }
    
    // Start solve-field process
    QProcessEnvironment enhancedEnv = createSolveFieldEnvironment();
    solveProcess.setProcessEnvironment(enhancedEnv);
    
    // Start solve-field process with absolute path and enhanced environment
    solveProcess.start("/opt/homebrew/bin/solve-field", arguments);
    
    if (!solveProcess.waitForStarted(5000)) {
        logMessage("Failed to start solve-field process", "red");
        return false;
    }
    
    // Monitor progress and keep UI responsive
    QElapsedTimer timer;
    timer.start();
    bool finished = false;
    
    while (!finished && timer.elapsed() < 600000) { // 1 minute timeout with hints
        finished = solveProcess.waitForFinished(1000);
        QApplication::processEvents(); // Keep UI responsive
        
        // Check if user wants to stop
        if (!m_processing) {
            solveProcess.kill();
            return false;
        }
        
        // Show some progress
        if (timer.elapsed() % 10000 == 0) { // Every 10 seconds
            logMessage(QString("solve-field running... (%1s elapsed)").arg(timer.elapsed() / 1000), "gray");
        }
    }
    
    if (!finished) {
        logMessage("solve-field timed out after 1 minute", "orange");
        solveProcess.kill();
        return false;
    }
    
    int exitCode = solveProcess.exitCode();
    QString output = solveProcess.readAllStandardOutput();
    QString errors = solveProcess.readAllStandardError();
    
    if (m_debugMode && !output.isEmpty()) {
        logMessage(QString("solve-field output: %1").arg(output.left(200)), "gray");
    }
    
    if (exitCode == 0) {
        // Check if output file was created
        if (QFile::exists(outputPath)) {
            logMessage("solve-field succeeded - image plate solved!", "green");
            
            // Extract coordinates from solve-field output
            QRegularExpression raRegex(R"(Field center: \(RA,Dec\) = \(([\d\.]+), ([\d\.]+)\))");
            QRegularExpressionMatch match = raRegex.match(output);
            if (match.hasMatch()) {
                double ra = match.captured(1).toDouble();
                double dec = match.captured(2).toDouble();
                logMessage(QString("solve-field found coordinates: RA=%1°, Dec=%2°").arg(ra, 0, 'f', 4).arg(dec, 0, 'f', 4), "blue");
            }
            
            return true;
        } else {
            logMessage(QString("solve-field completed but output file not found: %1").arg(outputPath), "red");
            return false;
        }
    } else {
        logMessage(QString("solve-field failed with exit code %1").arg(exitCode), "red");
        if (!errors.isEmpty()) {
            logMessage(QString("solve-field errors: %1").arg(errors.left(200)), "red");
        }
        return false;
    }
}
