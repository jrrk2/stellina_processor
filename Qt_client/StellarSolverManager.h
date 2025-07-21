// StellarSolver Integration - Replace solve-field with embedded StellarSolver
// This implementation processes images one at a time with progress bar updates

#include "StellinaProcessor.h"
#include <stellarsolver.h>
#include <parameters.h>
#include <structuredefinitions.h>
#include <QElapsedTimer>

class StellarSolverManager : public QObject {
    Q_OBJECT

public:
    StellarSolverManager(QObject* parent = nullptr) 
        : QObject(parent), m_solver(nullptr), m_currentImageIndex(0), m_totalImages(0) {
        setupParameters();
        findIndexFiles();
    }

    ~StellarSolverManager() {
        if (m_solver) {
            if (m_solver->isRunning()) {
                m_solver->abortAndWait();
            }
            m_solver->deleteLater();
        }
    }

    bool initializeBatch(const QStringList& imageFiles) {
        qDebug() << "=== initializeBatch DEBUG ===";
        
        // Check input parameter
        qDebug() << "Input parameter valid?" << (&imageFiles != nullptr);
        qDebug() << "Input size:" << imageFiles.size();
        
        if (imageFiles.isEmpty()) {
            qDebug() << "WARNING: Empty imageFiles list provided!";
        } else {
            qDebug() << "First file:" << imageFiles.first();
            qDebug() << "Last file:" << imageFiles.last();
        }
        
        // Check member variable state before assignment
        qDebug() << "m_imageFiles before assignment - size:" << m_imageFiles.size();
        
        try {
            // Try the assignment with error checking
            qDebug() << "Attempting assignment...";
            m_imageFiles = imageFiles;
            qDebug() << "Assignment successful!";
            
            // Verify assignment worked
            qDebug() << "m_imageFiles after assignment - size:" << m_imageFiles.size();
            
            if (m_imageFiles.size() != imageFiles.size()) {
                qDebug() << "ERROR: Size mismatch after assignment!";
                return false;
            }
            
        } catch (const std::exception& e) {
            qDebug() << "Exception during assignment:" << e.what();
            return false;
        } catch (...) {
            qDebug() << "Unknown exception during assignment!";
            return false;
        }
        
        // Rest of your function
        m_totalImages = imageFiles.size();
        m_currentImageIndex = 0;
        m_results.clear();
        
        if (m_indexPaths.isEmpty()) {
            qDebug() << "ERROR: No astrometry index files found!";
            emit errorOccurred("No astrometry index files found!");
            return false;
        }
        
        qDebug() << "initializeBatch completed successfully";
        return true;
    }

    // Start processing the next image in sequence
    void processNextImage() {
        if (m_currentImageIndex >= m_totalImages) {
            emit batchComplete();
            return;
        }

        const QString& currentFile = m_imageFiles[m_currentImageIndex];
        emit progressUpdated(m_currentImageIndex, m_totalImages, 
                           QString("Processing: %1").arg(QFileInfo(currentFile).baseName()));

        // Clean up previous solver
        if (m_solver) {
            m_solver->deleteLater();
        }

        // Create new solver for this image
        m_solver = new StellarSolver(this);
        
        // Connect signals for this solve
        connect(m_solver, &StellarSolver::finished, this, &StellarSolverManager::onSolveFinished);
        connect(m_solver, &StellarSolver::logOutput, this, &StellarSolverManager::onLogOutput);
        
        // Configure solver
        m_solver->setProperty("ProcessType", SSolver::SOLVE);
        m_solver->setProperty("ExtractorType", SSolver::EXTRACTOR_INTERNAL);
        m_solver->setProperty("SolverType", SSolver::SOLVER_STELLARSOLVER);
        m_solver->setIndexFolderPaths(m_indexPaths);
        m_solver->setParameters(m_params);

        // Load image and start solving
        if (loadFITSImage(currentFile)) {
            m_solveStartTime.start();
            m_solver->start();
        } else {
            // Skip this image and continue
            emit imageSkipped(currentFile, "Failed to load FITS image");
            m_currentImageIndex++;
            QTimer::singleShot(10, this, &StellarSolverManager::processNextImage);
        }
    }

    // Abort current processing
    void abortProcessing() {
        if (m_solver && m_solver->isRunning()) {
            m_solver->abortAndWait();
        }
    }

signals:
    void imageSolved(const QString& filename, double ra, double dec, double pixelScale);
    void progressUpdated(int current, int total, const QString& status);
    void imageProcessed(const QString& filename, bool success, double ra, double dec, double pixelScale);
    void imageSkipped(const QString& filename, const QString& reason);
    void batchComplete();
    void errorOccurred(const QString& error);
    void logOutput(const QString& message);

private slots:
    void onSolveFinished() {
        if (!m_solver) return;

        const QString& currentFile = m_imageFiles[m_currentImageIndex];
        bool success = false;
        double ra = -1, dec = -91, pixelScale = 0;

        if (m_solver->solvingDone() && m_solver->hasWCSData()) {
            FITSImage::Solution solution = m_solver->getSolution();
            success = true;
            ra = solution.ra;
            dec = solution.dec;
            pixelScale = solution.pixscale;
            
            emit logOutput(QString("✓ Solved: %1 - RA: %2°, Dec: %3°, Scale: %4 arcsec/pix")
                          .arg(QFileInfo(currentFile).baseName())
                          .arg(ra, 0, 'f', 4)
                          .arg(dec, 0, 'f', 4)
                          .arg(pixelScale, 0, 'f', 2));
        } else {
            QString reason = "Unknown failure";
            if (m_solver->failed()) {
                reason = "Solver failed";
            } else if (!m_solver->solvingDone()) {
                reason = "Solving incomplete";
            } else {
                reason = "No WCS data";
            }
            
            emit logOutput(QString("✗ Failed: %1 - %2 (%.1fs)")
                          .arg(QFileInfo(currentFile).baseName())
                          .arg(reason)
                          .arg(m_solveStartTime.elapsed() / 1000.0));
        }

        emit imageProcessed(currentFile, success, ra, dec, pixelScale);

        // Move to next image
        m_currentImageIndex++;
        
        // Process next image after a brief delay to keep UI responsive
        QTimer::singleShot(10, this, &StellarSolverManager::processNextImage);
    }

    void onLogOutput(const QString& message) {
        emit logOutput(QString("StellarSolver: %1").arg(message));
    }

private:
    bool convertStellinaToEquatorial(const StellinaImageData &stellinaData, double &ra, double &dec);
    void onStellarSolverImageSolved(const QString& filename, double ra, double dec, double pixelScale);
    void setupParameters() {
        // Get built-in profiles
        QList<Parameters> profiles = StellarSolver::getBuiltInProfiles();
        if (!profiles.isEmpty()) {
            m_params = profiles.at(0); // Use first profile as base
        }

        // Configure parameters for your use case
        m_params.multiAlgorithm = SSolver::MULTI_AUTO;
        m_params.search_radius = 10.0;  // Match your solve-field --radius
        m_params.minwidth = 0.1;
        m_params.maxwidth = 10.0;
        m_params.resort = true;
        m_params.autoDownsample = false;
        m_params.downsample = 1;
        m_params.inParallel = false;  // Single-threaded for progress control
        m_params.solverTimeLimit = 60;  // 60 second timeout

        // Star extraction parameters (tune these for your images)
        m_params.initialKeep = 2000;
        m_params.keepNum = 500;
        m_params.r_min = 1.0;
        m_params.removeBrightest = 0;
        m_params.removeDimmest = 50;
        m_params.saturationLimit = 65000;
        m_params.minarea = 5;
        m_params.threshold_offset = 0;
        m_params.threshold_bg_multiple = 2.0;
    }

    void findIndexFiles() {
        QStringList searchDirs = {
            "/opt/homebrew/share/astrometry",
            "/usr/local/share/astrometry", 
            "/usr/share/astrometry",
            QDir::homePath() + "/.local/share/astrometry"
        };

        for (const QString& dir : searchDirs) {
            QDir indexDir(dir);
            if (indexDir.exists()) {
                QStringList filters;
                filters << "index-*.fits";
                auto files = indexDir.entryList(filters, QDir::Files);
                if (!files.isEmpty()) {
                    m_indexPaths << indexDir.absolutePath();
                    break;
                }
            }
        }
    }

    bool loadFITSImage(const QString& filename) {
        if (!m_solver) return false;

        // Load FITS file using CFITSIO
        fitsfile *fptr;
        int status = 0;
        
        if (fits_open_file(&fptr, filename.toLocal8Bit().data(), READONLY, &status)) {
            return false;
        }

        // Read coordinate hints if available (STELLRA/STELLDEC)
        double fitsRA = -1, fitsDec = -91;
        char comment[FLEN_COMMENT];
        
        status = 0;
        bool hasSTELLRA = (fits_read_key(fptr, TDOUBLE, "STELLRA", &fitsRA, comment, &status) == 0);
        status = 0;
        bool hasSTELLDEC = (fits_read_key(fptr, TDOUBLE, "STELLDEC", &fitsDec, comment, &status) == 0);
        
        // Set coordinate hints if available
        if (hasSTELLRA && hasSTELLDEC && 
            fitsRA >= 0 && fitsRA <= 360 && 
            fitsDec >= -90 && fitsDec <= 90) {
            m_solver->setSearchPositionInDegrees(fitsRA, fitsDec);
        }

        // Get image dimensions
        int naxis, bitpix;
        long naxes[2];
        status = 0;
        if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) || naxis != 2) {
            fits_close_file(fptr, &status);
            return false;
        }

        int width = static_cast<int>(naxes[0]);
        int height = static_cast<int>(naxes[1]);
        long npixels = width * height;

        // Read image data
        m_imageBuffer.resize(npixels);
        long firstPix[2] = {1, 1};
        if (fits_read_pix(fptr, TFLOAT, firstPix, npixels, nullptr, 
                         m_imageBuffer.data(), nullptr, &status)) {
            fits_close_file(fptr, &status);
            return false;
        }

        fits_close_file(fptr, &status);

        // Calculate statistics
         float minVal = *std::min_element(m_imageBuffer.begin(), m_imageBuffer.end());
         float maxVal = *std::max_element(m_imageBuffer.begin(), m_imageBuffer.end());
         
         // Downsample with Bayer-aware processing
         const int downsampleFactor = 2;
         int outputWidth = width / downsampleFactor;
         int outputHeight = height / downsampleFactor;
         
         std::vector<uint8_t> buffer(outputWidth * outputHeight);
         float range = maxVal - minVal;
         
         if (range > 0) {
             for (int y = 0; y < outputHeight; y++) {
                 for (int x = 0; x < outputWidth; x++) {
                     int srcY = y * downsampleFactor;
                     int srcX = x * downsampleFactor;
                     
                     if (srcY + 1 >= height || srcX + 1 >= width) continue;
                     
                     // 2x2 Bayer pattern with proper weighting
                     float r = m_imageBuffer[srcY * width + srcX];
                     float g1 = m_imageBuffer[srcY * width + srcX + 1];
                     float g2 = m_imageBuffer[(srcY + 1) * width + srcX];
                     float b = m_imageBuffer[(srcY + 1) * width + srcX + 1];
                     
                     float avg = (r + (g1 + g2) * 0.5f + b) / 3.0f;
                     float normalized = (avg - minVal) / range;
                     buffer[y * outputWidth + x] = static_cast<uint8_t>(std::clamp(normalized * 255.0f, 0.0f, 255.0f));
                 }
             }
         } else {
             std::fill(buffer.begin(), buffer.end(), 128);
         }

        // Create statistics
         FITSImage::Statistic stats{};
         stats.width = outputWidth;
         stats.height = outputHeight;
         stats.channels = 1;
         stats.dataType = TBYTE;
         stats.bytesPerPixel = 1;
         
         for (int i = 0; i < 3; i++) {
             stats.min[i] = (i == 0) ? minVal : 0.0;
             stats.max[i] = (i == 0) ? maxVal : 0.0;
             stats.mean[i] = (i == 0) ? (minVal + maxVal) / 2.0 : 0.0;
             stats.stddev[i] = 0.0;
             stats.median[i] = (i == 0) ? (minVal + maxVal) / 2.0 : 0.0;
         }
         stats.SNR = 1.0;
         
         // Store buffer in solver-specific storage to keep it alive
         m_imageBuffers[m_solver] = std::move(buffer);
         
         // Load into solver
         return m_solver->loadNewImageBuffer(stats, m_imageBuffers[m_solver].data());
        
        return true;
    }

private:
    StellarSolver* m_solver;
    Parameters m_params;
    QStringList m_indexPaths;
    
    QStringList m_imageFiles;
    int m_currentImageIndex;
    int m_totalImages;
    
    std::vector<float> m_imageBuffer;
    QHash<StellarSolver*, std::vector<uint8_t>> m_imageBuffers; // Keep buffers alive per solver
    std::vector<uint8_t> m_imageBuffer8;
    QElapsedTimer m_solveStartTime;
    
    QList<QPair<QString, bool>> m_results;
};
