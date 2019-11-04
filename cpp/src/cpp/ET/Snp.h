/**
 * @author Tyler Cowman
 * Stores the case and control information for a single SNP, also stores the orginal index in the full dataset.
 * Constructed from plaintext reperesentations of each genotype, i,e. a one byte 0,1,2 per genotype. During
 * construction the genotypes are compressed using an implementation of the BOOST approach from Wan et al.
 */
#ifndef SNP_H
#define SNP_H

#include <iostream>

#include <vector>
#include <array>
#include <stdint.h>
#include <cmath>

#define GENOTYPE_LEVELS 3
#define GENOTYPE_PAIRINGS 9
#define CONTINGENCY_COLUMNS 2

#define PACK_SIZE 64
#define PACK_TYPE uint64_t
#define POPCOUNT_FUNCTION __builtin_popcountll
using namespace std;

static int popCount(const vector<PACK_TYPE> & v, int begin, int distance)
{
    int retVal = 0;
    for(int i=begin; i<begin+distance; ++i)
        retVal += POPCOUNT_FUNCTION(v[i]);

    return retVal;
}

static vector<PACK_TYPE> bitMerge(const vector<PACK_TYPE> & v1, const vector<PACK_TYPE> & v2)
{
    vector<PACK_TYPE>  retVal;
    retVal.reserve(v1.size());
    
    for(int i=0; i<v1.size(); ++i)
        retVal.push_back(v1[i] & v2[i]);
    
    return retVal;
}

static int popCountAnd(const vector<PACK_TYPE> & v1, const vector<PACK_TYPE> & v2, int begin, int distance)
{
    int retVal = 0;

    for(int i=begin; i<begin+distance; ++i)
    {
        if(v1[begin] | v2[begin] !=0)
        retVal += POPCOUNT_FUNCTION(v1[i] & v2[i]);
    }
    return retVal;
}

static int popCountAnd(const vector<PACK_TYPE> & v1, const vector<PACK_TYPE> & v2, int begin1, int begin2, int distance)
{
    int retVal = 0;
    //#pragma acc parallel reduction(+:retVal)
    for(int i=0; i<distance; ++i)
    {
        //There are many instnaces of long strings of zeroes, primarily due to homozygous minor this saves popcount operations
        if(v1[begin1+i] | v2[begin2+i] !=0)
            retVal += POPCOUNT_FUNCTION(v1[begin1 +i] & v2[begin2 +i]);
    }
    return retVal;
}

struct SmallContingencyTable
{
    array<int,GENOTYPE_LEVELS*CONTINGENCY_COLUMNS> M_;
    array<int,GENOTYPE_LEVELS> RT_;
    array<int,CONTINGENCY_COLUMNS> CT_;
    
    int & a(int x, int y)
    {
        return M_[x+(GENOTYPE_LEVELS*y)];
    }
    
    void addOne()
    {
        for(int i=0; i<M_.size(); ++i)
            M_[i]++;
    }
    
    void calculateTotals()
    {
        CT_.fill(0);
        RT_.fill(0);
        
        for(int i=0; i<GENOTYPE_LEVELS; ++i)
        {
            RT_[i] = M_[i] + M_[i+GENOTYPE_LEVELS];
            CT_[0] += M_[i];
            CT_[1] += M_[i+GENOTYPE_LEVELS];
        }
    }
    
};

struct ContingencyTable
{
    array<int,(GENOTYPE_LEVELS*GENOTYPE_LEVELS+1)*(CONTINGENCY_COLUMNS+1)-1> M_;

    int & a(int x, int y)
    {
        return M_[x+(GENOTYPE_LEVELS*GENOTYPE_LEVELS*y)];
    }
    
    void printM()
    {
            
        cout<<M_[0]<<" "<<M_[1]<<endl;
        cout<<M_[2]<<" "<<M_[3]<<endl;
        cout<<M_[4]<<" "<<M_[5]<<endl;
        cout<<M_[6]<<" "<<M_[7]<<endl;
        cout<<M_[8]<<" "<<M_[9]<<endl;
        cout<<M_[10]<<" "<<M_[11]<<endl;
        cout<<M_[12]<<" "<<M_[13]<<endl;
        cout<<M_[14]<<" "<<M_[15]<<endl;
        cout<<M_[16]<<" "<<M_[17]<<endl;
        cout<<endl;
    }
    
    
    void calculateTotals()
    {

        for(int i=0; i<GENOTYPE_LEVELS*GENOTYPE_LEVELS; ++i)
        {

            M_[GENOTYPE_LEVELS*GENOTYPE_LEVELS*2+i] = M_[i] + M_[i+GENOTYPE_LEVELS*GENOTYPE_LEVELS];
            M_[GENOTYPE_LEVELS*GENOTYPE_LEVELS*3] += M_[i];
            M_[GENOTYPE_LEVELS*GENOTYPE_LEVELS*3+1] += M_[i+GENOTYPE_LEVELS*GENOTYPE_LEVELS];
        }
        
    }
};

class Snp
{
    public:
        Snp(int index, const vector<char> & controls, const vector<char> & cases, int weight);
        Snp(const Snp & cpy);
        
        //Getters
        int getIndex()const;
        int getWeight()const;
        int getCaseNo()const;
        
        //Calculations
        float computeMinorAlleleFrequency()const;
        int computePopCoverAnd(const Snp & other)const;
        //Comparisons
        int computeDifferences(const Snp & other)const;
        float computeUnknownRatio()const;
        
        float marginalTest()const;
        float epistasisTest(const Snp & other)const;
        
        friend ostream& operator<< (ostream &out, const Snp & snp);
        //Will consist of controls 0 , 1 , 2 then cases 0 , 1 , 2
        vector<PACK_TYPE> allSamples_;
    private:
        
        array<vector<PACK_TYPE>, GENOTYPE_LEVELS> packGenotypesHelper(const vector<char> & genotypes); 
        void packGenotypes(vector<PACK_TYPE> & allPackedGenotypes, const vector<char> & controlGenotypes, const vector<char> & caseGenotypes);

        int index_;
        int weight_;

        static int numCases_;
        static int numControls_;
        static int CONR_;
        static int CASR_;
        static int CASS_;        
};
#endif //SNP_H