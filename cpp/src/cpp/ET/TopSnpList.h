/**
 * @author Tyler Cowman
 * Stores the current most significant interaction for each SNP. 
 */

#ifndef TOP_SNP_LIST_H
#define TOP_SNP_LIST_H

#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>

#define MAX_CUTOFF 100
#define PREFIX_SUM_ROLLOVER 1000

using namespace std;

//Used in the analysis and output of top pairs, not in discovery and cutoff updating
struct TopPairing
{
    TopPairing(int snpIndex1, int snpIndex2, float score)
    {
        score_ = score;
        if(snpIndex1<snpIndex2)
        {
            indexes_.first = snpIndex1;
            indexes_.second = snpIndex2;
        }
        else
        {
            indexes_.first = snpIndex2;
            indexes_.second = snpIndex1;
        }
    }
      
    bool operator <(const TopPairing &other )const
    {
        return indexes_<other.indexes_;
    }
    
    static bool orderByScore(const TopPairing & a, const TopPairing &  b)
    {
        return a.score_ < b.score_;
    }
    
    friend ostream& operator<< (ostream &out, const TopPairing & topPairing)
    {
        out<<topPairing.score_<<" "<<topPairing.indexes_.first<<" "<<topPairing.indexes_.second;
        
        return out;
    }
    
    float score_;
    pair<int, int> indexes_;
};

class TopSnpList
{
    public:
        TopSnpList();
        TopSnpList(const TopSnpList & cpy);
        TopSnpList(int topK, int numberSnps, float cutoff);
        
        bool attemptInsert(int snpIndex1, int snpIndex2, float score);
        
        float getCutoff()const;
        
        void incrementInternalTestsCounter(long long int testsDone);
        void incrementLeafTestsCounter(long long int testsDone);
        
        long long int getInternalTests()const;
        long long int getLeafTests()const;
        vector<TopPairing> getReciprocalPairs()const;
        
        void calculateFormattedResults();
        float calculateRecall(string groundTruthFileName);
        float calculatePrecision(string groundTruthFileName);
        float calculateRecallCutoffPairs(string groundTruthFileName);
        
        void printReciprocalPairs(ofstream & ofs)const;
        void printCutoffPairs(ofstream & ofs)const;
        
        void printReciprocalPairScores(ofstream & ofs)const;
        void printCutoffPairScores(ofstream & ofs)const;
        
        int getNumberOfReciprocalPairs()const;
        
        friend ostream& operator<< (ostream &out, const TopSnpList & topSnpList);
        
    private:
        void calculateReciprocalPairs();
        
        vector<float> topScores_;
        vector<int> topPartners_;
        long long int internalTestsCounter_;
        long long int leafTestsCounter_;
        
        float cutoff_;
        int topK_;

        array<int, MAX_CUTOFF+1> pairwiseSignificanceCounts_;
        array<int, MAX_CUTOFF+1> pairwiseSignificanceCountsPrefixSum_;
        int prefixSumTimer_;
    
        //Formatted results
        vector<TopPairing> reciprocalPairs_;
        vector<TopPairing> cutoffPairs_;
    
};
#endif //TOP_SNP_LIST_H
