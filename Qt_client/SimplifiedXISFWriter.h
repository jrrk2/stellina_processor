#ifndef SIMPLIFIED_XISF_WRITER_H
#define SIMPLIFIED_XISF_WRITER_H

#include <QString>
#include <QVariant>
#include <QJsonObject>
#include <memory>

// Forward declarations - no PCL includes in header!
class SimplifiedXISFWriterPrivate;

// Your project types (adjust paths as needed)
struct StellinaImageData;
namespace FITSImage { struct Solution; }

enum class CompressionType {
    None,
    ZLib,
    LZ4,
    ZSTD
};

enum class ChecksumType {
    None,
    SHA1,
    SHA256,
    SHA512
};

class SimplifiedXISFWriter
{
public:
    explicit SimplifiedXISFWriter(const QString& filePath, 
                                  CompressionType compression = CompressionType::ZLib);
    ~SimplifiedXISFWriter();

    // Image management
    bool addImage(const QString& id, 
                  const float* pixels, 
                  int width, 
                  int height, 
                  int channels = 1);

    // Property management  
    void addProperty(const QString& id, 
                     const QString& type, 
                     const QVariant& value, 
                     const QString& comment = QString());

    void addImageProperty(const QString& id, 
                          const QString& type, 
                          const QVariant& value, 
                          const QString& comment = QString());

    // Convenience methods for specific metadata
    void addWCSData(const FITSImage::Solution& solution);
    void addStellinaMetadata(const StellinaImageData& imageData);
    void addProcessingHistory(const QString& step, const QJsonObject& parameters);

    // Configuration
    void setCompression(CompressionType compression, int level = 6);
    void setChecksum(ChecksumType checksum);
    void setCreatorApplication(const QString& appName);
    void setVerbosity(int level); // 0=silent, 1=normal, 2=verbose

    // File operations
    bool write();
    void clear();

    // Status
    QString lastError() const;
    bool hasImages() const;
    int imageCount() const;

private:
    // PIMPL idiom - all PCL details hidden in implementation
    std::unique_ptr<SimplifiedXISFWriterPrivate> d;
    
    // Non-copyable
    SimplifiedXISFWriter(const SimplifiedXISFWriter&) = delete;
    SimplifiedXISFWriter& operator=(const SimplifiedXISFWriter&) = delete;
};

#endif // SIMPLIFIED_XISF_WRITER_H
