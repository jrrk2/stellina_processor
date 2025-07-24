// Standalone XISF Format Checker
// Compile with: g++ -std=c++17 xisf_checker.cpp -o xisf_checker `pkg-config --cflags --libs Qt5Core Qt5Xml`

#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QXmlStreamReader>
#include <QDebug>
#include <QTextStream>
#include <QCryptographicHash>
#include <QtEndian>
#include <iostream>

class XISFFormatChecker {
public:
    struct ValidationResult {
        bool isValid = false;
        QStringList errors;
        QStringList warnings;
        QStringList info;
        
        void addError(const QString& msg) { errors << msg; }
        void addWarning(const QString& msg) { warnings << msg; }
        void addInfo(const QString& msg) { info << msg; }
    };
    
    struct ImageInfo {
        QString id;
        int width = 0;
        int height = 0;
        int channels = 0;
        QString sampleFormat;
        QString colorSpace;
        qint64 location = 0;
        qint64 size = 0;
        QString checksum;
    };
    
    struct PropertyInfo {
        QString id;
        QString type;
        QString value;
        QString comment;
    };

private:
    QString m_filePath;
    QFile m_file;
    ValidationResult m_result;
    QList<ImageInfo> m_images;
    QList<PropertyInfo> m_properties;
    qint64 m_xmlHeaderSize = 0;
    qint64 m_blockAlignment = 4096;
    
public:
    XISFFormatChecker(const QString& filePath) : m_filePath(filePath), m_file(filePath) {}
    
    ValidationResult validate() {
        m_result = ValidationResult();
        
        if (!m_file.exists()) {
            m_result.addError("File does not exist: " + m_filePath);
            return m_result;
        }
        
        if (!m_file.open(QIODevice::ReadOnly)) {
            m_result.addError("Cannot open file for reading: " + m_file.errorString());
            return m_result;
        }
        
        m_result.addInfo(QString("File size: %1 bytes").arg(m_file.size()));
        
        // Step 1: Validate XML header
        if (!validateXMLHeader()) {
            return m_result;
        }
        
        // Step 2: Validate block alignment
        if (!validateBlockAlignment()) {
            return m_result;
        }
        
        // Step 3: Validate image data blocks
        if (!validateImageBlocks()) {
            return m_result;
        }
        
        // Step 4: Final validation
        if (m_result.errors.isEmpty()) {
            m_result.isValid = true;
            m_result.addInfo("✓ File passes all XISF monolithic format checks");
        }
        
        return m_result;
    }
    
    void printReport() {
        QTextStream out(stdout);
        
        out << "=== XISF Format Checker Report ===" << Qt::endl;
        out << "File: " << QFileInfo(m_filePath).fileName() << Qt::endl;
        out << "Path: " << m_filePath << Qt::endl;
        out << Qt::endl;
        
        // Overall status
        if (m_result.isValid) {
            out << "✓ VALID: File is a properly formatted monolithic XISF file" << Qt::endl;
        } else {
            out << "✗ INVALID: File has format issues" << Qt::endl;
        }
        out << Qt::endl;
        
        // Errors
        if (!m_result.errors.isEmpty()) {
            out << "ERRORS:" << Qt::endl;
            for (const QString& error : m_result.errors) {
                out << "  ✗ " << error << Qt::endl;
            }
            out << Qt::endl;
        }
        
        // Warnings
        if (!m_result.warnings.isEmpty()) {
            out << "WARNINGS:" << Qt::endl;
            for (const QString& warning : m_result.warnings) {
                out << "  ⚠ " << warning << Qt::endl;
            }
            out << Qt::endl;
        }
        
        // Information
        if (!m_result.info.isEmpty()) {
            out << "INFORMATION:" << Qt::endl;
            for (const QString& info : m_result.info) {
                out << "  ℹ " << info << Qt::endl;
            }
            out << Qt::endl;
        }
        
        // Properties summary
        if (!m_properties.isEmpty()) {
            out << "PROPERTIES (" << m_properties.size() << "):" << Qt::endl;
            for (const auto& prop : m_properties) {
                out << QString("  %1 (%2) = %3").arg(prop.id, prop.type, prop.value);
                if (!prop.comment.isEmpty()) {
                    out << " // " << prop.comment;
                }
                out << Qt::endl;
            }
            out << Qt::endl;
        }
        
        // Images summary
        if (!m_images.isEmpty()) {
            out << "IMAGES (" << m_images.size() << "):" << Qt::endl;
            for (const auto& img : m_images) {
                out << QString("  %1: %2x%3x%4 %5 %6")
                       .arg(img.id)
                       .arg(img.width).arg(img.height).arg(img.channels)
                       .arg(img.sampleFormat, img.colorSpace);
                out << QString(" @ location %1 (%2 bytes)").arg(img.location).arg(img.size);
                if (!img.checksum.isEmpty()) {
                    out << " checksum: " << img.checksum;
                }
                out << Qt::endl;
            }
        }
    }

private:
    bool validateXMLHeader() {
        m_file.seek(0);
        QByteArray xmlData;
        
        // Read until we find the end of XML or reasonable limit
        const qint64 maxXMLSize = 1024 * 1024; // 1MB max for XML header
        char buffer[8192];
        bool foundXMLEnd = false;
        
        while (!foundXMLEnd && xmlData.size() < maxXMLSize) {
            qint64 bytesRead = m_file.read(buffer, sizeof(buffer));
            if (bytesRead <= 0) break;
            
            xmlData.append(buffer, bytesRead);
            
            // Look for </xisf> to find end of XML
            if (xmlData.contains("</xisf>")) {
                // Find the exact position
                int endPos = xmlData.indexOf("</xisf>") + 7; // 7 = length of "</xisf>"
                xmlData = xmlData.left(endPos);
                foundXMLEnd = true;
            }
        }
        
        if (!foundXMLEnd) {
            m_result.addError("Could not find end of XML header (</xisf> tag)");
            return false;
        }
        
        m_xmlHeaderSize = xmlData.size();
        m_result.addInfo(QString("XML header size: %1 bytes").arg(m_xmlHeaderSize));
        
        // Parse XML
        QXmlStreamReader xml(xmlData);
        
        if (!xml.readNextStartElement()) {
            m_result.addError("No XML start element found");
            return false;
        }
        
        if (xml.name() != "xisf") {
            m_result.addError("Root element is not 'xisf', found: " + xml.name().toString());
            return false;
        }
        
        // Check version
        QString version = xml.attributes().value("version").toString();
        if (version != "1.0") {
            m_result.addWarning("XISF version is not '1.0', found: " + version);
        }
        
        // Check namespace
        QString xmlns = xml.attributes().value("xmlns").toString();
        if (xmlns != "http://www.pixinsight.com/xisf") {
            m_result.addError("Incorrect xmlns: " + xmlns);
            return false;
        }
        
        m_result.addInfo("✓ Valid XISF root element with correct namespace");
        
        // Parse elements
        while (xml.readNextStartElement()) {
            if (xml.name() == "Metadata") {
                if (!parseMetadata(xml)) {
                    return false;
                }
            } else if (xml.name() == "Image") {
                if (!parseImage(xml)) {
                    return false;
                }
            } else {
                xml.skipCurrentElement();
            }
        }
        
        if (xml.hasError()) {
            m_result.addError("XML parsing error: " + xml.errorString());
            return false;
        }
        
        return true;
    }
    
    bool parseMetadata(QXmlStreamReader& xml) {
        while (xml.readNextStartElement()) {
            if (xml.name() == "Property") {
                PropertyInfo prop;
                prop.id = xml.attributes().value("id").toString();
                prop.type = xml.attributes().value("type").toString();
                prop.value = xml.attributes().value("value").toString();
                prop.comment = xml.attributes().value("comment").toString();
                
                if (prop.id.isEmpty() || prop.type.isEmpty()) {
                    m_result.addWarning("Property missing id or type");
                } else {
                    m_properties.append(prop);
                    
                    // Check for required properties
                    if (prop.id == "XISF:BlockAlignmentSize") {
                        bool ok;
                        m_blockAlignment = prop.value.toInt(&ok);
                        if (!ok || m_blockAlignment <= 0) {
                            m_result.addError("Invalid BlockAlignmentSize: " + prop.value);
                            return false;
                        }
                        m_result.addInfo(QString("Block alignment: %1 bytes").arg(m_blockAlignment));
                    }
                }
                
                xml.skipCurrentElement();
            } else {
                xml.skipCurrentElement();
            }
        }
        return true;
    }
    
    bool parseImage(QXmlStreamReader& xml) {
        ImageInfo img;
        img.id = xml.attributes().value("id").toString();
        
        QString geometry = xml.attributes().value("geometry").toString();
        QStringList geomParts = geometry.split(':');
        if (geomParts.size() >= 3) {
            img.width = geomParts[0].toInt();
            img.height = geomParts[1].toInt();
            img.channels = geomParts[2].toInt();
        }
        
        img.sampleFormat = xml.attributes().value("sampleFormat").toString();
        img.colorSpace = xml.attributes().value("colorSpace").toString();
        
        QString location = xml.attributes().value("location").toString();
        if (location.startsWith("attachment:")) {
            QStringList locParts = location.mid(11).split(':'); // Remove "attachment:"
            if (locParts.size() >= 2) {
                img.location = locParts[0].toLongLong();
                img.size = locParts[1].toLongLong();
            }
        }
        
        img.checksum = xml.attributes().value("checksum").toString();
        
        m_images.append(img);
        xml.skipCurrentElement();
        return true;
    }
    
    bool validateBlockAlignment() {
        // Check if XML header is properly aligned
        qint64 alignedHeaderSize = ((m_xmlHeaderSize + m_blockAlignment - 1) / m_blockAlignment) * m_blockAlignment;
        
        if (alignedHeaderSize != m_xmlHeaderSize) {
            // Check if there's zero padding
            m_file.seek(m_xmlHeaderSize);
            QByteArray padding = m_file.read(alignedHeaderSize - m_xmlHeaderSize);
            
            bool isZeroPadded = true;
            for (int i = 0; i < padding.size(); ++i) {
                if (padding[i] != 0) {
                    isZeroPadded = false;
                    break;
                }
            }
            
            if (isZeroPadded) {
                m_result.addInfo(QString("✓ Header properly zero-padded from %1 to %2 bytes")
                               .arg(m_xmlHeaderSize).arg(alignedHeaderSize));
            } else {
                m_result.addError(QString("Header not properly zero-padded (size: %1, should be: %2)")
                                .arg(m_xmlHeaderSize).arg(alignedHeaderSize));
                return false;
            }
        } else {
            m_result.addInfo("✓ Header size is already block-aligned");
        }
        
        return true;
    }
    
    bool validateImageBlocks() {
        for (const auto& img : m_images) {
            // Check if location is properly aligned
            if (img.location % m_blockAlignment != 0) {
                m_result.addWarning(QString("Image %1 location %2 is not block-aligned")
                                  .arg(img.id).arg(img.location));
            }
            
            // Check if we can read the expected amount of data
            if (!m_file.seek(img.location)) {
                m_result.addError(QString("Cannot seek to image %1 location %2").arg(img.id).arg(img.location));
                continue;
            }
            
            QByteArray imageData = m_file.read(img.size);
            if (imageData.size() != img.size) {
                m_result.addError(QString("Image %1: expected %2 bytes, read %3 bytes")
                                .arg(img.id).arg(img.size).arg(imageData.size()));
                continue;
            }
            
            // Validate checksum if present
            if (!img.checksum.isEmpty()) {
                QString expectedChecksum = img.checksum;
                if (expectedChecksum.startsWith("sha1:")) {
                    expectedChecksum = expectedChecksum.mid(5);
                }
                
                QCryptographicHash hash(QCryptographicHash::Sha1);
                hash.addData(imageData);
                QByteArray actualChecksum = hash.result().toHex();
                
                if (actualChecksum.toLower() != expectedChecksum.toLower()) {
                    m_result.addError(QString("Image %1 checksum mismatch: expected %2, got %3")
                                    .arg(img.id, expectedChecksum, QString(actualChecksum)));
                } else {
                    m_result.addInfo(QString("✓ Image %1 checksum verified").arg(img.id));
                }
            }
            
            // Basic sanity checks
            qint64 expectedSize = static_cast<qint64>(img.width) * img.height * img.channels;
            if (img.sampleFormat == "Float32") {
                expectedSize *= 4;
            } else if (img.sampleFormat == "UInt16") {
                expectedSize *= 2;
            }
            
            if (img.size != expectedSize) {
                m_result.addWarning(QString("Image %1 size mismatch: declared %2 bytes, expected %3 bytes")
                                  .arg(img.id).arg(img.size).arg(expectedSize));
            }
            
            m_result.addInfo(QString("✓ Image %1 data block validated").arg(img.id));
        }
        
        return true;
    }
};

int main(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <xisf_file>" << std::endl;
        std::cout << "Validates XISF file format and reports issues" << std::endl;
        return 1;
    }
    
    QString filePath = argv[1];
    
    XISFFormatChecker checker(filePath);
    XISFFormatChecker::ValidationResult result = checker.validate();
    checker.printReport();
    
    return result.isValid ? 0 : 1;
}
