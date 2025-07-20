#pragma once

#include <QMainWindow>
#include <QMenuBar>
#include <QLabel>
#include <QScrollArea>
#include <QSlider>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <QProgressBar>
#include <QStatusBar>
#include <QApplication>
#include <QKeySequence>
#include <QAction>
#include <QMenu>
#include <fitsio.h>
#include <vector>
#include <algorithm>
#include <cmath>

/**
 * @brief Stretched Image Viewer for astronomical images
 * 
 * This class provides a comprehensive viewer for FITS astronomical images
 * with various stretching algorithms to enhance visibility of faint details.
 * 
 * Features:
 * - Multiple stretch algorithms (Linear, Gamma, Log, Sinh, Auto)
 * - Interactive parameter adjustment
 * - Zoom controls
 * - Statistics calculation
 * - Save functionality
 * 
 * Usage:
 * ```cpp
 * StretchedImageViewer *viewer = new StretchedImageViewer(parent);
 * viewer->loadStackedImage("/path/to/stacked.fits");
 * viewer->showViewer();
 * ```
 */
class StretchedImageViewer : public QMainWindow
{
    Q_OBJECT

public:
    /**
     * @brief Construct a new Stretched Image Viewer
     * @param parent Parent widget
     */
    explicit StretchedImageViewer(QWidget *parent = nullptr);
    
    /**
     * @brief Destructor
     */
    ~StretchedImageViewer();

    /**
     * @brief Load a FITS image for viewing
     * @param fitsPath Path to the FITS file
     */
    void loadStackedImage(const QString &fitsPath);
    
    /**
     * @brief Show the viewer window
     */
    void showViewer();

    /**
     * @brief Check if an image is currently loaded
     * @return true if image is loaded, false otherwise
     */
    bool hasImageLoaded() const { return !m_imageData.empty(); }
    
    /**
     * @brief Get the current image path
     * @return Path to currently loaded image
     */
    QString getCurrentImagePath() const { return m_currentImagePath; }

public slots:
    /**
     * @brief Apply auto-stretch to the current image
     */
  //    void autoStretch();
    
    /**
     * @brief Reset zoom to 100%
     */
  //    void resetZoom();

private slots:
    void onLoadImage();
    void onStretchTypeChanged();
    void onStretchParametersChanged();
    void onZoomChanged(int value);
    void onResetZoom();
    void onSaveStretchedImage();
    void onAutoStretch();
    //    void onShowStretchedViewer();
    //    void onLoadImageInViewer();

private:
    // UI Components
    QWidget *m_centralWidget;
    QScrollArea *m_scrollArea;
    QLabel *m_imageLabel;
    
    // Controls
    QGroupBox *m_stretchControlsGroup;
    QComboBox *m_stretchTypeCombo;
    QDoubleSpinBox *m_gammaSpinBox;
    QDoubleSpinBox *m_blackPointSpinBox;
    QDoubleSpinBox *m_whitePointSpinBox;
    QSlider *m_zoomSlider;
    QLabel *m_zoomLabel;
    QPushButton *m_loadButton;
    QPushButton *m_saveButton;
    QPushButton *m_autoStretchButton;
    QPushButton *m_resetZoomButton;
    QProgressBar *m_progressBar;
    
    // Image data
    QString m_currentImagePath;
    QImage m_originalImage;
    QImage m_stretchedImage;
    std::vector<float> m_imageData;
    long m_width, m_height;
    float m_minValue, m_maxValue;
    double m_currentZoom;

    // Setup methods
    void setupUI();
    void setupMenuBar();
    void connectSignals();
    
    // Image processing methods
    bool loadFITSImage(const QString &filePath);
    void calculateImageStatistics();
    void applyStretch();
    void updateImageDisplay();
    void updateZoom();
//    void autoStretch();
//    void resetZoom();
    void logMessage(const QString &message);
    
    QImage convertToQImage(const std::vector<float> &data, long width, long height);
    
    // Stretch algorithms
    std::vector<float> applyLinearStretch(const std::vector<float> &data);
    std::vector<float> applyGammaStretch(const std::vector<float> &data);
    std::vector<float> applyLogStretch(const std::vector<float> &data);
    std::vector<float> applySinhStretch(const std::vector<float> &data);
    std::vector<float> applyAutoStretch(const std::vector<float> &data);
    
    // Utility methods
    float percentile(const std::vector<float> &data, float percent);
//    void logMessage(const QString &message);
};
