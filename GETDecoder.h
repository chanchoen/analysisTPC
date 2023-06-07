#ifndef __GETDecoder_h
#define __GETDecoder_h

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>

#define nAGET 4       // The number of AGETs
#define nCHAN 68      // The number of AGET's channels
#define nBUCK 512     // The number of channel's timebuckets
#define nADCs 139264  // The number of ADC(12 bits); nAGET * nCHAN * nBUCK = 139264

class GETDecoder {
   private:
    struct Header_b {
        uint8_t metaType;       // 1
        uint8_t frameSize[3];   // 4
        uint8_t dataSource;     // 5
        uint8_t frameType[2];   // 7
        uint8_t revision;       // 8
        uint8_t headerSize[2];  // 10
        uint8_t itemSize[2];    // 12
        uint8_t nItems[4];      // 16
        uint8_t eventTime[6];   // 22
        uint8_t eventId[4];     // 26
        uint8_t coboId;         // 27
        uint8_t asadId;         // 28
        uint8_t readOffset[2];  // 30
        uint8_t status;         // 31
        uint8_t hitPat_0[9];    // 40
        uint8_t hitPat_1[9];    // 49
        uint8_t hitPat_2[9];    // 58
        uint8_t hitPat_3[9];    // 67
        uint8_t multip_0[2];    // 69
        uint8_t multip_1[2];    // 71
        uint8_t multip_2[2];    // 73
        uint8_t multip_3[2];    // 75
        uint8_t windowOut[4];   // 79
        uint8_t lastCell_0[2];  // 81
        uint8_t lastCell_1[2];  // 83
        uint8_t lastCell_2[2];  // 85
        uint8_t lastCell_3[2];  // 87 bytes
    };
    Header_b header_b;
    uint64_t prevTime_;
    uint32_t item_;

    std::ifstream inFile_;
    std::ifstream inFileList_;
    std::queue<std::string> fileQueue_;

   public:
    // Header
    int MetaType;           // [bytes]
    int FrameSize;          //
    int DataSource;         //
    int FrameType;          // 1: partial readout or zero suppressed mode, 2: full readout mode
    int Revision;           //
    int HeaderSize;         // [bytes]
    int ItemSize;           // [bytes]
    int nItems;             //
    int EventID;            //
    int CoboID;             //
    int AsadID;             //
    int ReadOffset;         //
    int Status;             //
    int multip[nAGET];      //
    int windowOut;          //
    int lastCell[nAGET];    //
    bool IsHit[nAGET][72];  //
    uint64_t EventTime;     // [10 ns]
    uint64_t DiffTime;      // [10 ns]
    // Items
    uint16_t ADC[nAGET][nCHAN][nBUCK];
    //
    GETDecoder();
    void Open(const std::string& fileName);
    void OpenFromList(const std::string& fileListName);
    void Close();

    void ReadHeader();
    void ReadItems();
    bool Run();
    void PrintHeader();

    bool IsFPN(int chanID);
    bool IsDead(int chanID);
    int FPNChanId(int chanID);
};
#endif
