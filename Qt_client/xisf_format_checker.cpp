// Fixed XISF Format Checker that works around QXmlStreamReader namespace issues
// Compile with: g++ -std=c++17 fixed_xisf_checker.cpp -o fixed_xisf_checker `pkg-config --cflags --libs Qt5Core Qt5Xml`

#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QXmlStreamReader>
#include <QDebug>
#include <QTextStream>
#include <QRegularExpression>
#include <iostream>

class FixedXISFChecker {
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

private:
    const QString m_filePath;
    QFile m_file;
    ValidationResult m_result;
    
public:
    FixedXISFChecker(const QString& filePath) : m_filePath(filePath), m_file(filePath) {}
    
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
        
        // Read XML header completely
        const QByteArray xmlData = readCompleteXMLHeader();
        if (xmlData.isEmpty()) {
            m_result.addError("Could not read complete XML header (</xisf> end tag not found within limit)");
            return m_result;
        }
        
        m_result.addInfo(QString("XML header size: %1 bytes").arg(xmlData.size()));
        
        // Validate XML structure using text/regex instead of QXmlStreamReader
        if (!validateXMLStructure(xmlData)) {
            return m_result;
        }
        
        // Validate block alignment with padding
        if (!validateBlockAlignment(xmlData.size())) {
            return m_result;
        }
        
        // Validate image blocks referenced in XML header
        if (!validateImageBlocks(xmlData)) {
            return m_result;
        }
        
        if (m_result.errors.isEmpty()) {
            m_result.isValid = true;
            m_result.addInfo("✓ File passes all XISF monolithic format checks");
        }
        
        return m_result;
    }
    
    void printReport() const {
        QTextStream out(stdout);
        
        out << "=== Fixed XISF Format Checker Report ===" << Qt::endl;
        out << "File: " << QFileInfo(m_filePath).fileName() << Qt::endl;
        out << "Path: " << m_filePath << Qt::endl;
        out << Qt::endl;
        
        if (m_result.isValid) {
            out << "✓ VALID: File is a properly formatted monolithic XISF file" << Qt::endl;
        } else {
            out << "✗ INVALID: File has format issues" << Qt::endl;
        }
        out << Qt::endl;
        
        if (!m_result.errors.isEmpty()) {
            out << "ERRORS:" << Qt::endl;
            for (const QString& error : m_result.errors) {
                out << "  ✗ " << error << Qt::endl;
            }
            out << Qt::endl;
        }
        
        if (!m_result.warnings.isEmpty()) {
            out << "WARNINGS:" << Qt::endl;
            for (const QString& warning : m_result.warnings) {
                out << "  ⚠ " << warning << Qt::endl;
            }
            out << Qt::endl;
        }
        
        if (!m_result.info.isEmpty()) {
            out << "INFORMATION:" << Qt::endl;
            for (const QString& info : m_result.info) {
                out << "  ℹ " << info << Qt::endl;
            }
        }
    }

private:
    // Read the entire XML header by scanning for "</xisf>" tag.
    // Limit read size to avoid infinite loops on corrupted files.
    QByteArray readCompleteXMLHeader() {
        m_file.seek(0);
        QByteArray data;
        
        constexpr int chunkSize = 4096;
        constexpr qint64 maxReadSize = 100 * 1024 * 1024; // 100 MB max to avoid runaway
        char buffer[chunkSize];
        
        while (data.size() < maxReadSize) {
            qint64 bytesRead = m_file.read(buffer, chunkSize);
            if (bytesRead <= 0) {
                break;
            }
            
            data.append(buffer, bytesRead);
            
            int endPos = data.indexOf("</xisf>");
            if (endPos != -1) {
                endPos += 7; // Include the length of "</xisf>"
                return data.left(endPos);
            }
        }
        
        return QByteArray(); // Not found or too big
    }
    
    bool validateBlockAlignment(qint64 xmlSize) {
        constexpr qint64 blockAlignment = 4096;
        
        // Calculate expected padded size aligned to blockAlignment
        const qint64 expectedPaddedSize = ((xmlSize + blockAlignment - 1) / blockAlignment) * blockAlignment;
        
        m_result.addInfo(QString("Raw XML size: %1 bytes").arg(xmlSize));
        m_result.addInfo(QString("Expected padded size (block aligned): %1 bytes").arg(expectedPaddedSize));
        
        if (expectedPaddedSize > m_file.size()) {
            m_result.addError(QString("File too small to contain padded XML header: file size %1 < expected %2")
                              .arg(m_file.size()).arg(expectedPaddedSize));
            return false;
        }
        
        if (expectedPaddedSize == xmlSize) {
            m_result.addInfo("✓ XML size is already block-aligned");
            return true;
        }
        
        // Check zero padding bytes after raw XML data up to padded size
        if (!m_file.seek(xmlSize)) {
            m_result.addError(QString("Failed to seek to XML end at byte %1").arg(xmlSize));
            return false;
        }
        
        const qint64 paddingSize = expectedPaddedSize - xmlSize;
        QByteArray padding = m_file.read(paddingSize);
        
        if (padding.size() != paddingSize) {
            m_result.addError(QString("Could not read full padding bytes, expected %1 got %2")
                              .arg(paddingSize).arg(padding.size()));
            return false;
        }
        
        for (int i = 0; i < padding.size(); ++i) {
            if (padding[i] != 0) {
                m_result.addError(QString("Non-zero byte found in XML header padding at offset %1").arg(xmlSize + i));
                return false;
            }
        }
        
        m_result.addInfo(QString("✓ Header properly zero-padded from %1 to %2 bytes")
                        .arg(xmlSize).arg(expectedPaddedSize));
        
        return true;
    }
    
    bool validateXMLStructure(const QByteArray& xmlData) {
        const QString xmlStr = QString::fromUtf8(xmlData);
        
        // Check XML declaration
        if (!xmlStr.contains(QLatin1String("<?xml version=\"1.0\""))) {
            m_result.addError("Missing or invalid XML declaration");
            return false;
        }
        m_result.addInfo("✓ Valid XML declaration found");
        
        // Check XISF root element with version attribute
        const QRegularExpression xisfRegex(R"(<xisf\s+version="1\.0"[^>]*>)");
        if (!xisfRegex.match(xmlStr).hasMatch()) {
            m_result.addError("Invalid or missing XISF root element");
            return false;
        }
        m_result.addInfo("✓ Valid XISF root element found");
        
        // Check namespace
        if (!xmlStr.contains(QLatin1String("xmlns=\"http://www.pixinsight.com/xisf\""))) {
            m_result.addError("Missing or incorrect XISF namespace");
            return false;
        }
        m_result.addInfo("✓ Correct XISF namespace found");
        
        // Check for required Metadata section
        if (!xmlStr.contains(QLatin1String("<Metadata>"))) {
            m_result.addError("Missing Metadata section");
            return false;
        }
        m_result.addInfo("✓ Metadata section found");
        
        // Check for BlockAlignmentSize property
        if (!xmlStr.contains(QLatin1String("XISF:BlockAlignmentSize"))) {
            m_result.addError("Missing XISF:BlockAlignmentSize property");
            return false;
        }
        m_result.addInfo("✓ BlockAlignmentSize property found");
        
        // Check Image element presence
        const QRegularExpression imageRegex(R"(<Image\s+[^>]*id="[^"]*"[^>]*>)");
        if (!imageRegex.match(xmlStr).hasMatch()) {
            m_result.addError("Missing or invalid Image element");
            return false;
        }
        m_result.addInfo("✓ Image element found");
        
        return true;
    }
    
    bool validateImageBlocks(const QByteArray& xmlData) {
        const QString xmlStr = QString::fromUtf8(xmlData);
        
        // Extract the first image location and size using regex
        // Example: location="attachment:4096:123456"
        QRegularExpression locationRegex(R"(location=\"attachment:(\d+):(\d+)\")");
        QRegularExpressionMatch locationMatch = locationRegex.match(xmlStr);
        
        if (!locationMatch.hasMatch()) {
            m_result.addError("Could not parse image location");
            return false;
        }
        
        const qint64 imageLocation = locationMatch.captured(1).toLongLong();
        const qint64 imageSize = locationMatch.captured(2).toLongLong();
        
        m_result.addInfo(QString("Image location: %1, size: %2 bytes").arg(imageLocation).arg(imageSize));
        
        if (imageLocation % 4096 != 0) {
            m_result.addWarning("Image location is not 4096-byte aligned");
        } else {
            m_result.addInfo("✓ Image location is properly block-aligned");
        }
        
        if (!m_file.seek(imageLocation)) {
            m_result.addError(QString("Cannot seek to image location %1").arg(imageLocation));
            return false;
        }
        
        QByteArray imageData = m_file.read(imageSize);
        if (imageData.size() != imageSize) {
            m_result.addError(QString("Expected %1 bytes of image data, but read %2 bytes")
                              .arg(imageSize).arg(imageData.size()));
            return false;
        }
        
        m_result.addInfo("✓ Image data block validated");
        return true;
    }
};

int main(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <xisf_file>" << std::endl;
        std::cout << "Fixed XISF validator that works around QXmlStreamReader namespace issues" << std::endl;
        return 1;
    }
    
    const QString filePath = argv[1];
    
    FixedXISFChecker checker(filePath);
    const FixedXISFChecker::ValidationResult result = checker.validate();
    checker.printReport();
    
    return result.isValid ? 0 : 1;
}
