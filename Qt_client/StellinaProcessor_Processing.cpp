#include "StellinaProcessor.h"
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
// Simplified processNextImage that doesn't handle stage transitions
void StellinaProcessor::processNextImage() {
    if (m_currentImageIndex >= m_imagesToProcess.length()) {
        // Current stage complete
        if (m_processingMode == MODE_FULL_PIPELINE) {
            // Handle pipeline stage transitions
            handlePipelineStageTransition();
        } else {
            // Single stage complete
            finishProcessing();
        }
        return;
    }
    
    QString currentFile = m_imagesToProcess[m_currentImageIndex];
    
    logMessage(QString("Processing %1 of %2: %3")
                  .arg(m_currentImageIndex + 1)
                  .arg(m_imagesToProcess.length())
                  .arg(QFileInfo(currentFile).fileName()), "blue");
    
    bool success = false;
    
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION:
        success = processImageDarkCalibration(currentFile);
        break;
    case STAGE_PLATE_SOLVING:
        success = processImagePlatesolving(currentFile);
        break;
    case STAGE_REGISTRATION:
    case STAGE_STACKING:
        // These stages will be handled by performAstrometricStacking()
        success = true; // Just advance through the list
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
    m_currentImageIndex++;
    m_progressBar->setValue(m_currentImageIndex);
    updateProcessingStatus();
    
    // Update time estimate
    qint64 elapsed = QDateTime::currentMSecsSinceEpoch() - m_processingStartTime;
    double avgTimePerImage = static_cast<double>(elapsed) / m_currentImageIndex;
    qint64 remaining = static_cast<qint64>((m_imagesToProcess.length() - m_currentImageIndex) * avgTimePerImage);
    
    m_timeEstimateLabel->setText(QString("Estimated time remaining: %1")
                                    .arg(formatProcessingTime(remaining)));
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
            m_currentStage = STAGE_REGISTRATION;
            if (performAstrometricStackingEnhanced()) {
                m_currentStage = STAGE_COMPLETE;
                finishProcessing();
            } else {
                logMessage("Astrometric stacking failed", "red");
                finishProcessing();
            }
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
        break;
    case STAGE_PLATE_SOLVING:
        // Keep existing registration/stacking status as ready for now
        break;
    case STAGE_REGISTRATION:
        m_registrationStatusLabel->setText("Registration: In Progress");
        break;
    case STAGE_STACKING:
        m_stackingStatusLabel->setText("Stacking: In Progress");
        break;
    case STAGE_COMPLETE:
        m_darkCalibrationStatusLabel->setText("Dark Calibration: Complete");
        m_registrationStatusLabel->setText("Registration: Complete");
        m_stackingStatusLabel->setText("Stacking: Complete");
        break;
    }
    
    if (m_processing) {
        m_progressLabel->setText(QString("Processing %1 of %2 images (Stage: %3) - Success: %4, Errors: %5")
                                    .arg(m_currentImageIndex + 1)
                                    .arg(m_imagesToProcess.length())
                                    .arg(stageText)
                                    .arg(m_processedCount)
                                    .arg(m_errorCount));
        
        if (m_darkCalibratedCount > 0) {
            m_darkCalibrationStatusLabel->setText(QString("Dark Calibration: %1 completed").arg(m_darkCalibratedCount));
        }
    }
}

QString StellinaProcessor::getStageDescription() const {
    switch (m_currentStage) {
    case STAGE_DARK_CALIBRATION: return "Dark Calibration";
    case STAGE_PLATE_SOLVING: return "Plate Solving";
    case STAGE_REGISTRATION: return "Registration";
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
        completionMessage += QString("Images processed: %1, Errors: %2, Total: %3")
                                .arg(m_processedCount)
                                .arg(m_errorCount)
                                .arg(m_imagesToProcess.length());
        if (!m_finalStackedImage.isEmpty()) {
            completionMessage += QString("\nFinal stack: %1").arg(QFileInfo(m_finalStackedImage).fileName());
        }
        break;
        
    case MODE_FULL_PIPELINE:
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
    
    // Save processing report
    saveProcessingReport();
    
    updateUI();
    
    // Show completion dialog
    QMessageBox::information(this, "Processing Complete", completionMessage);
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
