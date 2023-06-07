#include "GETDecoder.h"

GETDecoder::GETDecoder() {
}
void GETDecoder::Open(const std::string &fileName) {
    this->inFile_.open(fileName.c_str(), std::ios::binary);
    if (!this->inFile_.is_open()) {
        std::cerr << "! Error in GETDecoder::Open" << std::endl;
        std::cerr << "File " << fileName.c_str() << " does not exist!" << std::endl;
        exit(-1);
    }
    std::cout << "File open: " << fileName.c_str() << std::endl;
}
void GETDecoder::OpenFromList(const std::string &fileListName) {
    this->inFileList_.open(fileListName.c_str());
    if (!this->inFileList_.is_open()) {
        std::cerr << "! Error in GETDecoder::OpenFromList" << std::endl;
        std::cerr << "File " << fileListName.c_str() << " does not exist!" << std::endl;
        exit(-1);
    }
    while (!this->inFileList_.eof()) {
        std::string fileName;
        std::getline(this->inFileList_, fileName);
        if (this->inFileList_.eof()) {
            this->inFileList_.close();
            break;
        }
        this->fileQueue_.push(fileName);
    }
}
void GETDecoder::Close() {
    this->inFile_.close();
}
void GETDecoder::ReadHeader() {
    this->inFile_.read((char *)&(this->header_b), sizeof(Header_b));
    this->MetaType = (int)pow(2, (int)this->header_b.metaType);
    this->FrameSize = (int)((this->header_b.frameSize[0] << 16) | (this->header_b.frameSize[1] << 8) | (this->header_b.frameSize[2]));
    this->DataSource = (int)this->header_b.dataSource;
    this->FrameType = (int)((this->header_b.frameType[0] << 8) | this->header_b.frameType[1]);
    this->Revision = (int)this->header_b.revision;
    this->HeaderSize = (int)((this->header_b.headerSize[0] << 8) | this->header_b.headerSize[1]) * this->MetaType;
    this->ItemSize = (int)((this->header_b.itemSize[0] << 8) | this->header_b.itemSize[1]);
    this->nItems = (int)((this->header_b.nItems[0] << 24) | (this->header_b.nItems[1] << 16) |
                         (this->header_b.nItems[2] << 8) | (this->header_b.nItems[3]));
    this->EventTime = (((uint64_t)this->header_b.eventTime[0] << 40) | ((uint64_t)this->header_b.eventTime[1] << 32) |
                       ((uint64_t)this->header_b.eventTime[2] << 24) | ((uint64_t)this->header_b.eventTime[3] << 16) |
                       ((uint64_t)this->header_b.eventTime[4] << 8) | ((uint64_t)this->header_b.eventTime[5]));
    this->EventID = ((this->header_b.eventId[0] << 24) | (this->header_b.eventId[1] << 16) |
                     (this->header_b.eventId[2] << 8) | (this->header_b.eventId[3]));
    this->CoboID = (int)this->header_b.coboId;
    this->AsadID = (int)this->header_b.asadId;
    this->ReadOffset = (int)((this->header_b.readOffset[0] << 8) | this->header_b.readOffset[1]);
    this->Status = (int)this->header_b.status;

    int chanID = 0;
    for (int i = 0; i < 9; i++) {
        for (int shift = 7; shift >= 0; shift--) {
            if (chanID == 0) shift = 3;
            this->IsHit[0][chanID] = this->header_b.hitPat_0[i] & (1 << shift) ? true : false;
            this->IsHit[1][chanID] = this->header_b.hitPat_1[i] & (1 << shift) ? true : false;
            this->IsHit[2][chanID] = this->header_b.hitPat_2[i] & (1 << shift) ? true : false;
            this->IsHit[3][chanID] = this->header_b.hitPat_3[i] & (1 << shift) ? true : false;
            chanID++;
        }
    }
    //
    uint8_t waste;
    for (int i = 0; i < (this->HeaderSize - 87); i++) this->inFile_.read((char *)&waste, sizeof(uint8_t));
}
void GETDecoder::ReadItems() {
    std::fill_n(&this->ADC[0][0][0], nADCs, 0);  // Clear all ADCs
    if (this->FrameType == 1) {
        // Partial readout mode
        // Frame Item Format:     aacc cccc | cbbb bbbb | bb00 ssss | ssss ssss
        // Read in reverse order: ssss ssss | bb00 ssss | cbbb bbbb | aacc cccc
        // a: agetID, c: chanID, b: buckID, s: sample, 0: empty
        int agetID, chanID, buckID, sample;
        int chanIDLowerBits, chanIDUpperBits;
        int buckIDLowerBits, buckIDUpperBits;
        int sampleLowerBits, sampleUpperBits;
        for (int itemId = 0; itemId < this->nItems; itemId++) {
            this->inFile_.read((char *)&this->item_, sizeof(uint32_t));
            // Read agetID
            agetID = (this->item_ >> 6) & 0x03;
            // Read chanID
            chanIDLowerBits = (this->item_ >> 15) & 0x01;
            chanIDUpperBits = this->item_ & 0x3f;
            chanID = (chanIDUpperBits << 1) | chanIDLowerBits;
            // Read buckID
            buckIDLowerBits = (this->item_ >> 22) & 0x03;
            buckIDUpperBits = (this->item_ >> 8) & 0x7f;
            buckID = (buckIDUpperBits << 2) | buckIDLowerBits;
            // Read sample
            sampleLowerBits = (this->item_ >> 24) & 0xff;
            sampleUpperBits = (this->item_ >> 16) & 0x0f;
            sample = (sampleUpperBits << 8) | sampleLowerBits;
            this->ADC[agetID][chanID][buckID] = sample;
        }
    } else if (this->FrameType == 2) {
        // Full readout mode
        // Frame Item Format:     aa00 ssss | ssss ssss | aa00 ssss | ssss ssss
        // Read in reverse order: ssss ssss | aa00 ssss | ssss ssss | aa00 ssss
        // a: agetID, s: sample, 0: empty
        int sample, sampleLowerBits, sampleUpperBits;
        for (int buckID = 0; buckID < nBUCK; buckID++) {
            for (int chanID = 0; chanID < nCHAN; chanID += 2) {
                for (int agetID = 0; agetID < nAGET; agetID++) {
                    this->inFile_.read((char *)&this->item_, sizeof(uint32_t));
                    // Read sample 1
                    sampleLowerBits = (this->item_ >> 8) & 0xff;
                    sampleUpperBits = this->item_ & 0x0f;
                    sample = (sampleUpperBits << 8) | sampleLowerBits;
                    this->ADC[agetID][chanID][buckID] = sample;
                    // Read sample 2
                    sampleLowerBits = (this->item_ >> 24) & 0xff;
                    sampleUpperBits = (this->item_ >> 16) & 0x0f;
                    sample = (sampleUpperBits << 8) | sampleLowerBits;
                    this->ADC[agetID][chanID + 1][buckID] = sample;
                }
            }
        }
    }
}
bool GETDecoder::Run() {
    if (!this->inFile_.is_open() && this->fileQueue_.empty()) {
        std::cerr << "! Error in GETDecoder::Run" << std::endl;
        std::cerr << "File is not opened!" << std::endl;
        exit(-1);
    } else if (!this->inFile_.is_open() && !this->fileQueue_.empty()) {
        this->Open(this->fileQueue_.front());
        this->fileQueue_.pop();
    }
    this->ReadHeader();
    if (this->inFile_.eof() && this->fileQueue_.empty()) {
        this->inFile_.close();
        std::cout << "[End Of Decoding]" << std::endl;
        return false;
    } else if (this->inFile_.eof() && !this->fileQueue_.empty()) {
        this->inFile_.close();
        std::cout << "[End Of File]" << std::endl;
        this->Run();
        return true;
    }
    this->ReadItems();
    if (this->EventID == 0) {
        this->prevTime_ = EventTime;
        this->DiffTime = 0;
    } else {
        this->DiffTime = EventTime - this->prevTime_;
        this->prevTime_ = EventTime;
    }
    return true;
}
void GETDecoder::PrintHeader() {
    std::cout << std::endl;
    std::cout << "eventId: " << this->EventID << std::endl;
    std::cout << "eventTime: " << (double)EventTime / 1.E8 << " s" << std::endl;
    std::cout << "***************Header Status********************" << std::endl;
    std::cout << "metaType: " << this->MetaType << " bytes" << std::endl;
    std::cout << "frameSize: " << this->FrameSize << std::endl;
    std::cout << "dataSource: " << this->DataSource << std::endl;
    std::cout << "frameType: " << this->FrameType << std::endl;
    std::cout << "revision: " << this->Revision << std::endl;
    std::cout << "headerSize: " << this->HeaderSize << " bytes" << std::endl;
    std::cout << "itemSize: " << this->ItemSize << " bytes" << std::endl;
    std::cout << "nItems: " << this->nItems << std::endl;
    std::cout << "coboIdx: " << this->CoboID << std::endl;
    std::cout << "asadIdx: " << this->AsadID << std::endl;
    std::cout << "readOffset: " << this->ReadOffset << std::endl;
    std::cout << "status: " << this->Status << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << std::endl;
}
bool GETDecoder::IsFPN(int chanID) {
    if (chanID < 0 || chanID > 67) {
        std::cerr << "! Error in GETDecoder::IsFPN" << std::endl;
        std::cerr << "Out of range: chanID must be from 0 to 67" << std::endl;
        exit(-1);
    }
    return (chanID == 11) || (chanID == 22) || (chanID == 45) || (chanID == 56);
}
bool GETDecoder::IsDead(int chanID) {
    if (chanID < 0 || chanID > 67) {
        std::cerr << "! Error in GETDecoder::IsDead" << std::endl;
        std::cerr << "Out of range: chanID must be from 0 to 67" << std::endl;
        exit(-1);
    }
    return (chanID == 33) || (chanID == 35);
}
int GETDecoder::FPNChanId(int chanID) {
    if (0 <= chanID && chanID < 17) {
        return 11;
    } else if (chanID < 34) {
        return 22;
    } else if (chanID < 51) {
        return 45;
    } else if (chanID < 68) {
        return 56;
    } else {
        std::cerr << "! Error in GETDecoder::FPNChanID" << std::endl;
        std::cerr << "Out of range: chanID must be from 0 to 67" << std::endl;
        exit(-1);
    }
}
