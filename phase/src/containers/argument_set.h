#ifndef _ARGUMENT_SET_H
#define _ARGUMENT_SET_H

#include <utils/otools.h>

struct argument_set
{
    static std::vector<std::string> mLegalOptions;

    uint32_t mSeed;
    uint32_t mNumThreads;

    std::string mBamFileFilename;
    std::string mBamListFilename;
    std::string mInputGLFilename;
    std::string mRefHapFilename;
    std::string mMapFilename;
    std::string mInputRegion;
    std::string mOutputRegion;
    float mSparseMaf;
    std::string mSamplesFilename;
    std::string mIndName;

    uint32_t mBurnIn;
    uint32_t mMain;
    uint32_t mNe;
	float mMinGL;
    float mErrImp;
    float mErrPhase;

    uint32_t mMaxPbwtDepth;
	uint32_t mMinPbwtDepth;
    float mPbwtCM;
    uint32_t Kpbwt;
    uint32_t Kinit;
    std::string mStateListFilename;

    std::string mCallModel;
    std::string mFastaFilename;
    uint32_t mMapQ;
    uint32_t mBaseQ;
    uint32_t mMaxDepth;

    std::string mOutputFilename;
    std::string mContigFaiFilename;
    bool mPrintGP;
    bool mPrintDS;
    bool mPrintAP;
    unsigned int mBgenNbits;
    std::string mBgenCompression;
    std::string mLogFilename;

    bool mKeepMonomorphicRefSites;
    bool mImputeReferenceOnlyVariants;
    bool mInputFieldGL;
    bool mUseGLIndels;
    bool mCallIndels;
    bool mKeepOrphanReads;
    bool mIgnoreOrientation;
    bool mCheckProperPairing;
    bool mKeepFailedQC;
    bool mKeepDuplicates;
    bool mIllumina13;
    bool mOutputIndex;


/*
    bool mPrintGP;
    bool mPrintDS;
    bool mPrintAP;
    bool mPrintSurrogateFamDosage;
    std::string mSurfbatTest;
    float mSurfbatMaf;
    float mSurfbatInfo;
    //bool mStandardSelect;
    int mPloidyMainSamples;
    //uint8_t mBgenNbits;
*/


    std::string mTargetFilename;
	InputFormat mInputFormat;
    OutputFormat mOutputFormat;
    OutputCompression mOutputCompression;

};

#endif
