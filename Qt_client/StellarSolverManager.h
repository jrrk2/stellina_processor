#include <fitsio.h>
#include <stellarsolver.h>
#include <parameters.h>
#include <structuredefinitions.h>
#include <QElapsedTimer>
#include <QTimer>
// ===============================================================================
// CRITICAL FIX: StellarSolverManager Constructor Missing Member Initialization
// ===============================================================================

// The problem is in StellarSolverManager.h - the constructor doesn't initialize m_imageFiles!

// REPLACE the constructor in StellarSolverManager.h with this:

class StellarSolverManager : public QObject {
    Q_OBJECT

public:
    // FIXED CONSTRUCTOR - Initialize ALL member variables
    StellarSolverManager(QObject* parent = nullptr) 
        : QObject(parent)
        , m_solver(nullptr)
        , m_currentImageIndex(0)
        , m_totalImages(0)
        , m_imageFilesVector()          // ADD THIS LINE - Initialize QStringList
        , m_indexPaths()          // ADD THIS LINE - Initialize QStringList
        , m_results()             // ADD THIS LINE - Initialize QList
        , m_imageBuffer()         // ADD THIS LINE - Initialize vector
        , m_imageBuffer8()        // ADD THIS LINE - Initialize vector
    {
        setupParameters();
        findIndexFiles();
        // Explicitly initialize all member variables
        m_indexPaths = QStringList();
        m_results = QList<QPair<QString, bool>>();
        m_imageBuffer = std::vector<float>();
        m_imageBuffer8 = std::vector<uint8_t>();
        m_currentImageIndex = 0;
        m_totalImages = 0;
        
        qDebug() << "StellarSolverManager: Members explicitly initialized";
        
        // Verify initialization worked
        qDebug() << "StellarSolverManager initialized:";
        qDebug() << "  m_imageFiles size:" << m_imageFilesVector.size();
        qDebug() << "  m_indexPaths size:" << m_indexPaths.size();
        qDebug() << "  m_results size:" << m_results.size();
    }

    ~StellarSolverManager() {
        if (m_solver) {
            if (m_solver->isRunning()) {
                m_solver->abortAndWait();
            }
            m_solver->deleteLater();
        }
    }
    
    void processNextImage() {
        if (m_currentImageIndex >= m_totalImages) {
            emit batchComplete();
            return;
        }

        const QString& currentFile = m_imageFilesVector[m_currentImageIndex];
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
    
    // SAFE VERSION of initializeBatch
    bool initializeBatch(const QStringList& imageFiles) {
        qDebug() << "initializeBatch called with" << imageFiles.size() << "files";
        
        if (imageFiles.isEmpty()) {
            emit errorOccurred("No image files provided!");
            return false;
        }
        
        try {
	  m_imageFilesVector.clear();
	  for (const QString& file : imageFiles) {
            m_imageFilesVector.push_back(file);
	  }
            m_totalImages = m_imageFilesVector.size();
            m_currentImageIndex = 0;
            m_results.clear();
            
            qDebug() << "Assignment successful, m_imageFiles now has" << m_imageFilesVector.size() << "files";
            
        } catch (const std::exception& e) {
            qDebug() << "Exception during assignment:" << e.what();
            return false;
        } catch (...) {
            qDebug() << "Unknown exception during assignment";
            return false;
        }
        
        if (m_indexPaths.isEmpty()) {
            emit errorOccurred("No astrometry index files found!");
            return false;
        }
        
        return true;
    }

signals:
    void imageSolved(const QString& filename, double ra, double dec, double pixelScale);
    void progressUpdated(int current, int total, const QString& status);
    void imageProcessed(const QString& filename, bool success, double ra, double dec, double
 pixelScale);
    void imageSkipped(const QString& filename, const QString& reason);
    void batchComplete();
    void errorOccurred(const QString& error);
    void logOutput(const QString& message);

private slots:
    void onSolveFinished() {
        if (!m_solver) return;

        const QString& currentFile = m_imageFilesVector[m_currentImageIndex];
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
    std::vector<QString> m_imageFilesVector;
    StellarSolver* m_solver;
    Parameters m_params;
    
    QStringList m_indexPaths;                              // This was missing initialization!
    QList<QPair<QString, bool>> m_results;                 // This was missing initialization!
    
    int m_currentImageIndex;
    int m_totalImages;
    
    std::vector<float> m_imageBuffer;                      // This was missing initialization!
    std::vector<uint8_t> m_imageBuffer8;                   // This was missing initialization!
    QHash<StellarSolver*, std::vector<uint8_t>> m_imageBuffers;
    QElapsedTimer m_solveStartTime;
    
    // Add these methods if they don't exist:
    void setupParameters() {
        // Initialize solver parameters
        m_params.listName = "Default";
        // ... set other parameters
    }
    
    void findIndexFiles() {
        // Search for astrometry index files
        QStringList searchPaths = {
            "/usr/local/astrometry/data",
            "/opt/homebrew/share/astrometry",
            "/usr/share/astrometry",
            QDir::homePath() + "/.local/share/astrometry"
        };
        
        m_indexPaths.clear();
        
        for (const QString& path : searchPaths) {
            QDir dir(path);
            if (dir.exists()) {
                QStringList indexFiles = dir.entryList(QStringList() << "index-*.fits", QDir::Files);
                for (const QString& file : indexFiles) {
                    m_indexPaths << dir.absoluteFilePath(file);
                }
            }
        }
        
        qDebug() << "Found" << m_indexPaths.size() << "astrometry index files";
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
        return true;
    }
};

/*
// ===============================================================================
// ALTERNATIVE: If you can't modify the header, create a new constructor
// ===============================================================================

// If the header file is auto-generated or you can't modify it,
// add this initialization in the .cpp file:

void StellarSolverManager::initializeMembers() {
    // Call this immediately after construction
    m_imageFiles = QStringList();
    m_indexPaths = QStringList(); 
    m_results = QList<QPair<QString, bool>>();
    m_imageBuffer = std::vector<float>();
    m_imageBuffer8 = std::vector<uint8_t>();
    
    qDebug() << "Manually initialized member variables";
}

// Then call it like this:
// StellarSolverManager* manager = new StellarSolverManager(this);
// manager->initializeMembers();  // Add this line

// ===============================================================================
// IMMEDIATE WORKAROUND - Add this to your StellinaProcessor where you create the manager
// ===============================================================================

void StellinaProcessor::initializeStellarSolver() {
    if (m_stellarSolverManager) {
        delete m_stellarSolverManager;
    }
    
    m_stellarSolverManager = new StellarSolverManager(this);
    
    // WORKAROUND: Force proper initialization by calling a dummy operation
    QStringList emptyList;
    try {
        // This will initialize the internal QStringList properly
        bool result = m_stellarSolverManager->initializeBatch(emptyList);
        // We expect this to fail due to empty list, but it initializes the member
    } catch (...) {
        qDebug() << "Expected failure during dummy initialization";
    }
    
    qDebug() << "StellarSolver initialized and member variables reset";
}

// ===============================================================================
// ROOT CAUSE SUMMARY
// ===============================================================================
*/

/*
 * THE PROBLEM:
 * In StellarSolverManager.h, the constructor was:
 * 
 * StellarSolverManager(QObject* parent = nullptr) 
 *     : QObject(parent), m_solver(nullptr), m_currentImageIndex(0), m_totalImages(0) {
 * 
 * Notice: m_imageFiles, m_indexPaths, m_results are NOT in the initializer list!
 * This means they are default-constructed but may be in an undefined state.
 * 
 * THE FIX:
 * Add all member variables to the constructor initializer list:
 * 
 * StellarSolverManager(QObject* parent = nullptr) 
 *     : QObject(parent)
 *     , m_solver(nullptr)
 *     , m_currentImageIndex(0)
 *     , m_totalImages(0)
 *     , m_imageFiles()        // FIX: Add this
 *     , m_indexPaths()        // FIX: Add this
 *     , m_results()           // FIX: Add this
 *     , m_imageBuffer()       // FIX: Add this
 *     , m_imageBuffer8()      // FIX: Add this
 * 
 * This ensures all member variables are properly initialized before use.
 */


