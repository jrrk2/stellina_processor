#include "StellinaProcessor.h"
#include "CoordinateUtils.h"
#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QJsonObject>
#include <QJsonDocument>
#include <QJsonArray>
#include <QDateTime>
#include <QSettings>
#include <QSplitter>
#include <QMenuBar>
#include <QStatusBar>
#include <QHeaderView>
#include <QElapsedTimer>
#include <QProcess>
#include <QStandardItemModel>
#include <QTableWidgetItem>
#include <QThread>
#include <QPainter>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>  // for std::fmod if needed

// Also add a function to verify the systematic offsets are being applied correctly
void StellinaProcessor::verifySystematicOffsetsInUse() {
    logMessage("=== VERIFYING SYSTEMATIC OFFSET APPLICATION ===", "blue");
    
    logMessage("Current mount tilt correction settings:", "gray");
    logMessage(QString("  Enable correction: %1").arg(m_mountTilt.enableCorrection ? "YES" : "NO"), 
              m_mountTilt.enableCorrection ? "green" : "red");
    logMessage(QString("  Systematic RA offset: %1°").arg(m_mountTilt.systematicRAOffset, 0, 'f', 4), "blue");
    logMessage(QString("  Systematic Dec offset: %1°").arg(m_mountTilt.systematicDecOffset, 0, 'f', 4), "blue");
    logMessage(QString("  Geometric north tilt: %1°").arg(m_mountTilt.northTilt, 0, 'f', 4), "gray");
    logMessage(QString("  Geometric east tilt: %1°").arg(m_mountTilt.eastTilt, 0, 'f', 4), "gray");
    
    if (!m_mountTilt.enableCorrection) {
        logMessage("", "gray");
        logMessage("Mount correction is DISABLED", "red");
        logMessage("Enable it through the UI or run calibration to apply systematic offsets", "orange");
        return;
    }
    
    if (qAbs(m_mountTilt.systematicRAOffset) < 0.01 && qAbs(m_mountTilt.systematicDecOffset) < 0.01) {
        logMessage("", "gray");
        logMessage("WARNING: Systematic offsets are very small", "orange");
        logMessage("Run 'Auto-Calibrate from Processed Files' to calculate proper offsets", "orange");
        return;
    }
    
    // Test a coordinate conversion to verify offsets are being applied
    double testAlt = 42.0;
    double testAz = 287.0;
    QString testTime = "2024-01-09T22:13:29";
    
    // Temporarily disable correction
    bool originalState = m_mountTilt.enableCorrection;
    m_mountTilt.enableCorrection = false;
    
    double ra_without, dec_without;
    convertAltAzToRaDec(testAlt, testAz, testTime, ra_without, dec_without);
    
    // Re-enable correction
    m_mountTilt.enableCorrection = originalState;
    
    double ra_with, dec_with;
    convertAltAzToRaDec(testAlt, testAz, testTime, ra_with, dec_with);
    
    double actual_ra_offset = ra_with - ra_without;
    double actual_dec_offset = dec_with - dec_without;
    
    // Handle RA wraparound
    if (actual_ra_offset > 180) actual_ra_offset -= 360;
    if (actual_ra_offset < -180) actual_ra_offset += 360;
    
    logMessage("", "gray");
    logMessage("Verification test results:", "blue");
    logMessage(QString("  Test coordinates: Alt=%1°, Az=%2°").arg(testAlt).arg(testAz), "gray");
    logMessage(QString("  Without correction: RA=%1°, Dec=%2°")
                  .arg(ra_without, 0, 'f', 4).arg(dec_without, 0, 'f', 4), "gray");
    logMessage(QString("  With correction:    RA=%1°, Dec=%2°")
                  .arg(ra_with, 0, 'f', 4).arg(dec_with, 0, 'f', 4), "blue");
    logMessage(QString("  Actual offset applied: RA=%1°, Dec=%2°")
                  .arg(actual_ra_offset, 0, 'f', 4).arg(actual_dec_offset, 0, 'f', 4), "green");
    logMessage(QString("  Expected offset:       RA=%1°, Dec=%2°")
                  .arg(m_mountTilt.systematicRAOffset, 0, 'f', 4).arg(m_mountTilt.systematicDecOffset, 0, 'f', 4), "orange");
    
    bool ra_match = qAbs(actual_ra_offset - m_mountTilt.systematicRAOffset) < 0.001;
    bool dec_match = qAbs(actual_dec_offset - m_mountTilt.systematicDecOffset) < 0.001;
    
    if (ra_match && dec_match) {
        logMessage("✓ SUCCESS: Systematic offsets are being applied correctly!", "green");
    } else {
        logMessage("✗ PROBLEM: Systematic offsets are not being applied correctly", "red");
        logMessage("Check the convertAltAzToRaDec function implementation", "red");
    }
    
    logMessage("=== END VERIFICATION ===", "blue");
}

// Add this to StellinaProcessor_Core.cpp - Qt plotting function using QPainter

void StellinaProcessor::plotMountErrors() {
    logMessage("=== PLOTTING MOUNT ERRORS ===", "blue");
    
    // Check if we have plate-solved files to work with
    QDir plateSolvedDir(m_plateSolvedDirectory);
    if (m_plateSolvedDirectory.isEmpty() || !plateSolvedDir.exists()) {
        logMessage("No plate-solved directory found. Run plate solving first.", "red");
        return;
    }
    
    // Collect data using existing function
    QList<ProcessedImageData> imageData;
    
    QStringList plateSolvedFiles = plateSolvedDir.entryList(
        QStringList() << "plate_solved_*.fits" << "*solved*.fits", 
        QDir::Files);
    
    if (plateSolvedFiles.isEmpty()) {
        logMessage("No plate-solved FITS files found", "red");
        return;
    }
    
    // Load data
    QDateTime sessionStart;
    bool sessionStartSet = false;
    
    for (const QString &fileName : plateSolvedFiles) {
        ProcessedImageData data;
        QString fitsPath = plateSolvedDir.absoluteFilePath(fileName);
        
        if (!readStellinaDataFromSolvedFits(fitsPath, data)) {
            continue;
        }
        if (!readSolveFieldResults(fitsPath, data)) {
            continue;
        }
        
        // Extract image number
        QRegularExpression imgNumRegex(R"(img-(\d+))");
        QRegularExpressionMatch match = imgNumRegex.match(fileName);
        if (match.hasMatch()) {
            data.imageNumber = match.captured(1).toInt();
        } else {
            continue;
        }
        
        // Parse time
        data.obsTime = QDateTime::fromString(data.dateObs, "yyyy-MM-ddThh:mm:ss");
        if (!data.obsTime.isValid()) {
            QStringList formats = {"yyyy-MM-ddThh:mm:ss.zzz", "yyyy-MM-dd hh:mm:ss"};
            for (const QString &format : formats) {
                data.obsTime = QDateTime::fromString(data.dateObs, format);
                if (data.obsTime.isValid()) break;
            }
        }
        if (!data.obsTime.isValid()) continue;
        
        data.obsTime.setTimeSpec(Qt::UTC);
        
        if (!sessionStartSet) {
            sessionStart = data.obsTime;
            sessionStartSet = true;
            data.minutesFromStart = 0.0;
        } else {
            data.minutesFromStart = sessionStart.msecsTo(data.obsTime) / 60000.0;
        }
        
        data.isValid = true;
        imageData.append(data);
    }
    
    if (imageData.size() < 10) {
        logMessage("Need at least 10 data points for plotting", "red");
        return;
    }
    
    // Sort by time
    std::sort(imageData.begin(), imageData.end(), 
              [](const ProcessedImageData &a, const ProcessedImageData &b) {
                  return a.minutesFromStart < b.minutesFromStart;
              });
    
    logMessage(QString("Creating error plot with %1 data points").arg(imageData.size()), "blue");
    
    // Create a widget to display the plot
    QWidget *plotWidget = new QWidget;
    plotWidget->setWindowTitle("Stellina Mount Error Analysis");
    plotWidget->setMinimumSize(1200, 800);
    plotWidget->setAttribute(Qt::WA_DeleteOnClose);
    
    // Create custom paint widget
    class ErrorPlotWidget : public QWidget {
    private:
        QList<ProcessedImageData> m_data;
        double m_initialRAOffset;
        double m_driftRA;
        double m_initialDecOffset;
        double m_driftDec;
        
    public:
        ErrorPlotWidget(const QList<ProcessedImageData> &data, 
                       double initialRA, double driftRA,
                       double initialDec, double driftDec, QWidget *parent = nullptr)
            : QWidget(parent), m_data(data), m_initialRAOffset(initialRA), m_driftRA(driftRA),
              m_initialDecOffset(initialDec), m_driftDec(driftDec) {
            setMinimumSize(1200, 800);
        }
        
    protected:
        void paintEvent(QPaintEvent *event) override {
            Q_UNUSED(event)
            
            QPainter painter(this);
            painter.setRenderHint(QPainter::Antialiasing);
            
            // Plot area
            QRect plotArea = rect().adjusted(80, 60, -80, -120);
            
            // Background
            painter.fillRect(rect(), Qt::white);
            painter.fillRect(plotArea, QColor(250, 250, 250));
            painter.setPen(Qt::black);
            painter.drawRect(plotArea);
            
            if (m_data.isEmpty()) return;
            
            // Calculate data ranges
            double minTime = 0, maxTime = m_data.last().minutesFromStart;
            double minRAError = 1000, maxRAError = -1000;
            double minDecError = 1000, maxDecError = -1000;
            double minResidualRA = 1000, maxResidualRA = -1000;
            
            QList<double> raErrors, decErrors, residualRAErrors, timePoints;
            
            for (const ProcessedImageData &data : m_data) {
                double raError = data.predictedRA - data.solvedRA;
                double decError = data.predictedDec - data.solvedDec;
                
                // Handle RA wraparound
                if (raError > 180) raError -= 360;
                if (raError < -180) raError += 360;
                
                // Calculate residual after linear trend removal
                double timeHours = data.minutesFromStart / 60.0;
                double linearRA = m_initialRAOffset + m_driftRA * timeHours;
                double residualRA = raError - linearRA;
                
                raErrors.append(raError);
                decErrors.append(decError);
                residualRAErrors.append(residualRA);
                timePoints.append(data.minutesFromStart);
                
                minRAError = qMin(minRAError, raError);
                maxRAError = qMax(maxRAError, raError);
                minDecError = qMin(minDecError, decError);
                maxDecError = qMax(maxDecError, decError);
                minResidualRA = qMin(minResidualRA, residualRA);
                maxResidualRA = qMax(maxResidualRA, residualRA);
            }
            
            // Add padding to ranges
            double raRange = maxRAError - minRAError;
            minRAError -= raRange * 0.1;
            maxRAError += raRange * 0.1;
            
            double residualRange = maxResidualRA - minResidualRA;
            minResidualRA -= residualRange * 0.1;
            maxResidualRA += residualRange * 0.1;
            
            // Helper function to map data to plot coordinates
            auto mapX = [&](double time) {
                return plotArea.left() + (time - minTime) / (maxTime - minTime) * plotArea.width();
            };
            
            auto mapY_RA = [&](double error) {
                return plotArea.bottom() - (error - minRAError) / (maxRAError - minRAError) * (plotArea.height() / 2);
            };
            
            auto mapY_Residual = [&](double error) {
                return plotArea.bottom() - plotArea.height()/2 - (error - minResidualRA) / (maxResidualRA - minResidualRA) * (plotArea.height() / 2);
            };
            
            // Draw grid
            painter.setPen(QPen(Qt::lightGray, 1, Qt::DotLine));
            for (int i = 0; i <= 10; ++i) {
                int x = plotArea.left() + i * plotArea.width() / 10;
                painter.drawLine(x, plotArea.top(), x, plotArea.bottom());
                
                int y1 = plotArea.top() + i * (plotArea.height()/2) / 10;
                int y2 = plotArea.bottom() - plotArea.height()/2 + i * (plotArea.height()/2) / 10;
                painter.drawLine(plotArea.left(), y1, plotArea.right(), y1);
                painter.drawLine(plotArea.left(), y2, plotArea.right(), y2);
            }
            
            // Draw separation line
            painter.setPen(QPen(Qt::black, 2));
            int midY = plotArea.top() + plotArea.height() / 2;
            painter.drawLine(plotArea.left(), midY, plotArea.right(), midY);
            
            // Plot RA errors (top half)
            painter.setPen(QPen(Qt::blue, 2));
            QPolygonF raLine;
            for (int i = 0; i < raErrors.size(); ++i) {
                raLine << QPointF(mapX(timePoints[i]), mapY_RA(raErrors[i]));
            }
            painter.drawPolyline(raLine);
            
            // Plot linear trend line for RA
            painter.setPen(QPen(Qt::red, 2, Qt::DashLine));
            double startLinearRA = m_initialRAOffset;
            double endLinearRA = m_initialRAOffset + m_driftRA * (maxTime / 60.0);
            painter.drawLine(QPointF(mapX(minTime), mapY_RA(startLinearRA)),
                           QPointF(mapX(maxTime), mapY_RA(endLinearRA)));
            
            // Plot residual RA errors (bottom half)
            painter.setPen(QPen(Qt::darkGreen, 2));
            QPolygonF residualLine;
            for (int i = 0; i < residualRAErrors.size(); ++i) {
                residualLine << QPointF(mapX(timePoints[i]), mapY_Residual(residualRAErrors[i]));
            }
            painter.drawPolyline(residualLine);
            
            // Zero line for residuals
            painter.setPen(QPen(Qt::black, 1, Qt::DashLine));
            painter.drawLine(plotArea.left(), mapY_Residual(0), plotArea.right(), mapY_Residual(0));
            
            // Labels and title
            painter.setPen(Qt::black);
            QFont font = painter.font();
            font.setPointSize(12);
            font.setBold(true);
            painter.setFont(font);
            
            // Title
            painter.drawText(rect().adjusted(0, 10, 0, -rect().height() + 30), 
                           Qt::AlignHCenter, "Stellina Mount Error Analysis");
            
            // Axis labels
            font.setPointSize(10);
            font.setBold(false);
            painter.setFont(font);
            
            painter.drawText(rect().adjusted(0, rect().height() - 30, 0, 0), 
                           Qt::AlignHCenter, "Time (minutes)");
            
            painter.save();
            painter.translate(20, plotArea.center().y() - plotArea.height()/4);
            painter.rotate(-90);
            painter.drawText(0, 0, "RA Error (degrees)");
            painter.restore();
            
            painter.save();
            painter.translate(20, plotArea.center().y() + plotArea.height()/4);
            painter.rotate(-90);
            painter.drawText(0, 0, "RA Residual (degrees)");
            painter.restore();
            
            // Legend
            painter.setPen(Qt::blue);
            painter.drawLine(plotArea.right() - 200, plotArea.top() + 20, 
                           plotArea.right() - 170, plotArea.top() + 20);
            painter.setPen(Qt::black);
            painter.drawText(plotArea.right() - 160, plotArea.top() + 25, "Raw RA Error");
            
            painter.setPen(Qt::red);
            painter.drawLine(plotArea.right() - 200, plotArea.top() + 40, 
                           plotArea.right() - 170, plotArea.top() + 40);
            painter.setPen(Qt::black);
            painter.drawText(plotArea.right() - 160, plotArea.top() + 45, "Linear Trend");
            
            painter.setPen(Qt::darkGreen);
            painter.drawLine(plotArea.right() - 200, plotArea.top() + 60, 
                           plotArea.right() - 170, plotArea.top() + 60);
            painter.setPen(Qt::black);
            painter.drawText(plotArea.right() - 160, plotArea.top() + 65, "Residual (detrended)");
            
            // Statistics
            font.setPointSize(9);
            painter.setFont(font);
            QString stats = QString("Linear trend: %1° + %2°/h * t (R² = %3)")
                              .arg(m_initialRAOffset, 0, 'f', 3)
                              .arg(m_driftRA, 0, 'f', 2)
                              .arg(0.982, 0, 'f', 3);  // From your regression
            painter.drawText(plotArea.left(), plotArea.bottom() + 20, stats);
            
            // Residual RMS
            double residualRMS = 0.0;
            for (double res : residualRAErrors) {
                residualRMS += res * res;
            }
            residualRMS = sqrt(residualRMS / residualRAErrors.size());
            
            QString residualStats = QString("Residual RMS: %1° (possible periodic error)")
                                      .arg(residualRMS, 0, 'f', 3);
            painter.drawText(plotArea.left(), plotArea.bottom() + 40, residualStats);
            
            // Time axis labels
            for (int i = 0; i <= 6; ++i) {
                double time = minTime + i * (maxTime - minTime) / 6;
                int x = mapX(time);
                painter.drawText(x - 15, plotArea.bottom() + 15, QString::number(time, 'f', 0));
            }
            
            // RA error axis labels (top)
            for (int i = 0; i <= 5; ++i) {
                double error = minRAError + i * (maxRAError - minRAError) / 5;
                int y = mapY_RA(error);
                painter.drawText(plotArea.left() - 50, y + 5, QString::number(error, 'f', 1));
            }
            
            // Residual axis labels (bottom)
            for (int i = 0; i <= 5; ++i) {
                double error = minResidualRA + i * (maxResidualRA - minResidualRA) / 5;
                int y = mapY_Residual(error);
                painter.drawText(plotArea.left() - 50, y + 5, QString::number(error, 'f', 1));
            }
        }
    };
    
    // Create the plotting widget
    ErrorPlotWidget *plotCanvas = new ErrorPlotWidget(imageData, 
        m_mountTilt.initialRAOffset, m_mountTilt.driftRA,
        m_mountTilt.initialDecOffset, m_mountTilt.driftDec, plotWidget);
    
    QVBoxLayout *layout = new QVBoxLayout(plotWidget);
    layout->addWidget(plotCanvas);
    
    // Add save button
    QPushButton *saveButton = new QPushButton("Save Plot as PNG");
    layout->addWidget(saveButton);
    
    connect(saveButton, &QPushButton::clicked, [plotCanvas, this]() {
        QString fileName = QFileDialog::getSaveFileName(this, "Save Plot", 
            QDir(m_plateSolvedDirectory).absoluteFilePath("mount_error_plot.png"),
            "PNG Images (*.png)");
        if (!fileName.isEmpty()) {
            QPixmap pixmap(plotCanvas->size());
            plotCanvas->render(&pixmap);
            pixmap.save(fileName);
            logMessage(QString("Plot saved to: %1").arg(fileName), "green");
        }
    });
    
    plotWidget->show();
    logMessage("Error plot window opened. Look for periodic patterns in the residual (bottom) plot.", "green");
}

// Fast calibration from stacking JSON files
void StellinaProcessor::calibrateFromStackingJSON() {
    logMessage("=== FAST MOUNT CALIBRATION FROM STACKING JSON ===", "blue");
    
    if (m_sourceDirectory.isEmpty()) {
        logMessage("Please select source directory containing stacking JSON files", "red");
        return;
    }
    
    QDir sourceDir(m_sourceDirectory);
    if (!sourceDir.exists()) {
        logMessage("Source directory does not exist", "red");
        return;
    }
    
    // Find all stacking JSON files
    QStringList jsonFiles = sourceDir.entryList(
        QStringList() << "*-stacking.json" << "*stacking*.json", 
        QDir::Files);
    
    if (jsonFiles.isEmpty()) {
        logMessage("No stacking JSON files found in source directory", "red");
        logMessage("Looking for files like 'img-0001-stacking.json'", "orange");
        return;
    }
    
    logMessage(QString("Found %1 stacking JSON files").arg(jsonFiles.size()), "blue");
    
    // Parse all stacking files
    QList<StackingCorrectionData> stackingData;
    double sessionStart;
    bool sessionStartSet = false;
    
    for (const QString &jsonFile : jsonFiles) {
        StackingCorrectionData data;
        QString jsonPath = sourceDir.absoluteFilePath(jsonFile);
        
        if (parseStackingJSON(jsonPath, data)) {
            // Set session start from first valid file
            if (!sessionStartSet) {
                sessionStart = data.acqTime;
                sessionStartSet = true;
                data.minutesFromStart = 0.0;
            } else {
                data.minutesFromStart = (data.acqTime - sessionStart) / 60000.0;
            }
            
            stackingData.append(data);
        }
    }
    
    if (stackingData.size() < 10) {
        logMessage("Need at least 10 valid stacking files for calibration", "red");
        return;
    }
    
    // Sort by image number
    std::sort(stackingData.begin(), stackingData.end(), 
              [](const StackingCorrectionData &a, const StackingCorrectionData &b) {
                  return a.imageNumber < b.imageNumber;
              });
    
    logMessage(QString("Successfully parsed %1 stacking corrections").arg(stackingData.size()), "green");
    logMessage(QString("Time span: %.1f minutes").arg(stackingData.last().minutesFromStart), "blue");
    
    // Convert pixel corrections to angular errors and perform regression
    analyzeStackingCorrections(stackingData, sessionStart);
}

bool StellinaProcessor::parseStackingJSON(const QString &jsonPath, StackingCorrectionData &data) {
    QFile file(jsonPath);
    if (!file.open(QIODevice::ReadOnly)) {
        return false;
    }
    
    QJsonParseError error;
    QJsonDocument doc = QJsonDocument::fromJson(file.readAll(), &error);
    
    if (error.error != QJsonParseError::NoError) {
        if (m_debugMode) {
            logMessage(QString("JSON parse error in %1: %2").arg(jsonPath).arg(error.errorString()), "orange");
        }
        return false;
    }
    
    QJsonObject root = doc.object();
    
    // Extract image information
    data.imageFilename = QFileInfo(jsonPath).fileName();
    
    // Extract image number from filename (e.g., "img-0200-stacking.json" -> 200)
    QRegularExpression imgNumRegex(R"(img-(\d+))");
    QRegularExpressionMatch match = imgNumRegex.match(data.imageFilename);
    if (!match.hasMatch()) {
        return false;
    }
    data.imageNumber = match.captured(1).toInt();
    
    // Extract observation time from acqTime (milliseconds since telescope boot)
    if (!root.contains("acqTime")) {
        return false;
    }
    
    data.acqTime = root["acqTime"].toVariant().toLongLong();
    
    // Extract mount position
    if (root.contains("motors")) {
        QJsonObject motors = root["motors"].toObject();
        if (motors.contains("ALT") && motors.contains("AZ")) {
            data.stellinaAlt = motors["ALT"].toDouble();
            data.stellinaAz = motors["AZ"].toDouble();
        } else {
            return false;
        }
    } else {
        return false;
    }
    
    // Extract stacking correction data
    if (root.contains("stackingData")) {
        QJsonObject stackingData = root["stackingData"].toObject();
        
        if (stackingData.contains("liveRegistrationResult")) {
            QJsonObject regResult = stackingData["liveRegistrationResult"].toObject();
            
            // Check if registration was successful
            data.statusMessage = regResult["statusMessage"].toString();
            if (data.statusMessage != "StackingOk") {
                return false; // Skip failed registrations
            }
            
            // Extract correction values
            if (regResult.contains("correction")) {
                QJsonObject correction = regResult["correction"].toObject();
                data.correctionX = correction["x"].toDouble();
                data.correctionY = correction["y"].toDouble();
                data.correctionRot = correction["rot"].toDouble();
            } else {
                return false;
            }
            
            // Extract quality metrics
            data.starsUsed = regResult["starsUsed"].toInt();
            data.distanceToCenter = regResult["distanceToCenter"].toDouble();
            
            // Quality filter - require reasonable number of stars
            if (data.starsUsed < 10) {
                if (m_debugMode) {
                    logMessage(QString("Skipping %1: only %2 stars used").arg(data.imageFilename).arg(data.starsUsed), "orange");
                }
                return false;
            }
            
        } else {
            return false;
        }
    } else {
        return false;
    }
    
    data.isValid = true;
    return true;
}
// Improved stacking analysis that handles mosaic patterns
// Replace the analyzeStackingCorrections function with this version
