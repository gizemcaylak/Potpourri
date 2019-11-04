/**
 * @author Tyler Cowman
 * 
 * These structs serve as containers for the information about a given run of LinDen
 * and are primarily used to provide an organized location for this data when printing
 * to standard out or a file later.
 */


#ifndef COMMON_STRUCTS_H
#define COMMON_STRUCTS_H

#include <iostream>
#include <fstream>

using namespace std;

struct ParameterInfo
{
    int maxThreadUsage_;
    int permuteSamples_;
    int snpBeginIndex_;
    int snpEndIndex_;
    float maxUnknownFraction_;

    float minimumMinorAlleleFrequency_;
    float maxMarginalSignificance_;
    
    string outputFileName_;
    
    int noTrees_;

    
    void printSummaryRelevant(ofstream & ofs)
    {
        ofs <<maxMarginalSignificance_;
    }
    
};

struct DatasetSizeInfo
{
    int snps_;
    int cases_;
    int controls_;
    
    //Varies based on trial
    int mafRemoved_;
    int marginalSignificanceRemoved_;
    int filterFileRemoved_;
    
    int passingSnps_;
    int mergedTreesFormed_;
    
    void printDimensions(ofstream & ofs)
    {
        ofs<<snps_<<"\t"<<cases_<<"\t"<<controls_;
    }
    
    void printSubsetDimensions(ofstream & ofs)
    {
        ofs<<marginalSignificanceRemoved_;
    }
    
};


#endif //COMMON_STRUCTS_H
