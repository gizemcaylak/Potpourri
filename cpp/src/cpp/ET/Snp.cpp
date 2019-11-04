#include "Snp.h"

//Define static members, will be set correctly in the Snp constructors
int Snp::numControls_;
int Snp::numCases_;
int Snp::CONR_;
int Snp::CASR_;
int Snp::CASS_;

Snp::Snp(int index, const vector<char> & controls, const vector<char> & cases, int weight)
{
    packGenotypes(allSamples_, controls, cases);
    
    index_ = index;
    weight_ = weight;
    
    numControls_ = controls.size();
    numCases_ = cases.size();
    CONR_ = ceil(numControls_/(double)PACK_SIZE);
    CASR_ = ceil(numCases_/(double)PACK_SIZE);
    CASS_ = (GENOTYPE_LEVELS*CONR_) ;//+ 1;
}

Snp::Snp(const Snp & cpy)
{
    index_ = cpy.index_;
    weight_ = cpy.weight_;
    allSamples_ = cpy.allSamples_;
}

int Snp::getIndex()const
{
    return index_;
}

int Snp::getCaseNo()const
{
    return numCases_;
}

int Snp::getWeight()const
{
    return weight_;
}

int Snp::computePopCoverAnd(const Snp & other)const
{
    int hetero = popCountAnd(allSamples_, other.allSamples_, CASS_+CASR_, CASR_)- popCountAnd(allSamples_, other.allSamples_, CONR_, CONR_) ; 
    int homoMinor =  popCountAnd(allSamples_, other.allSamples_, CASS_+(2*CASR_), CASR_) -popCountAnd(allSamples_, other.allSamples_, 2*CONR_, CONR_);
   
    return hetero+homoMinor;
}

float Snp::computeMinorAlleleFrequency()const
{
    int homoMajor = popCount(allSamples_, 0, CONR_) + popCount(allSamples_, CASS_, CASR_) ;
    int hetero = popCount(allSamples_, CONR_, CONR_) + popCount(allSamples_, CASS_+CASR_, CASR_); 
    int homoMinor = popCount(allSamples_, 2*CONR_, CONR_) + popCount(allSamples_, CASS_+(2*CASR_), CASR_) ;
    
    return (2*homoMinor+hetero)/(float)(2*homoMajor + hetero + 2*homoMinor);
}

int Snp::computeDifferences(const Snp & other)const
{
    return ((numControls_ + numCases_)- popCountAnd(allSamples_, other.allSamples_, 0, allSamples_.size()));
}

float Snp::computeUnknownRatio()const
{
    int numberKnown = popCount(allSamples_, 0, allSamples_.size());

    return ( numberKnown/(float)(numControls_ + numCases_) );
}

float Snp::marginalTest()const
{
    float retVal=0.0;
    SmallContingencyTable t;
    t.M_.fill(0);
    for(int i=0; i<GENOTYPE_LEVELS; ++i)
    {
        t.M_[i] = popCount(allSamples_, i*CONR_, CONR_);
        t.M_[i+GENOTYPE_LEVELS] = popCount(allSamples_, CASS_ + (i*CASR_), CASR_);
    }
       
    t.addOne(); //For correction    
    t.calculateTotals();
        
		
    int factors = GENOTYPE_LEVELS;
    int levels = 2;
    
    for(int y=0; y< factors; ++y)
    {

        for(int x=0; x< levels; ++x)
        {            
            //Use margin values
            float c1 = (float)t.RT_[y];
            float c2 = (float)t.CT_[x];
            float c3 = (float)(t.CT_[0] + t.CT_[1]);          
            
            float c4 = (float)t.a(y,x);

            float expected  = ( c1 * c2 ) / c3 ;

            retVal = retVal +  (pow(c4-expected,2)/expected);

        }
    }
    return retVal;
}

float Snp::epistasisTest(const Snp & other)const
{
    float retVal=0.0;

    ContingencyTable t;
    t.M_.fill(0);
    
    //PLUS ONE FOR CORRECTION
    for(int i=0; i<GENOTYPE_LEVELS; ++i)
    {
        //Count controls
        for(int j=0; j<GENOTYPE_LEVELS; ++j)
        {
            t.M_[i+(GENOTYPE_LEVELS*j)] = popCountAnd(allSamples_, other.allSamples_, i*CONR_, j*CONR_, CONR_)+1;
            
        }
        //Count Cases
        for(int j=0; j<GENOTYPE_LEVELS; ++j)
        {
            t.M_[i+(GENOTYPE_LEVELS*j)+(GENOTYPE_LEVELS*GENOTYPE_LEVELS)] = popCountAnd(allSamples_, other.allSamples_, CASS_+(i*CASR_), CASS_+(j*CASR_), CASR_ )+1;
        }
    }   
    
    t.calculateTotals();

    //18 and 27 offset for row and column totals respectively
    for(int y=0; y< GENOTYPE_PAIRINGS; ++y)
    {

        for(int x=0; x< CONTINGENCY_COLUMNS; ++x)
        {            
            //Use margin values
            int c1 = t.M_[y+GENOTYPE_LEVELS*GENOTYPE_LEVELS*2];
            int c2 = t.M_[x+GENOTYPE_LEVELS*GENOTYPE_LEVELS*3];
            float c3 = (t.M_[GENOTYPE_LEVELS*GENOTYPE_LEVELS*3] + t.M_[GENOTYPE_LEVELS*GENOTYPE_LEVELS*3+1]);


            float c4 = t.M_[y+(GENOTYPE_LEVELS*GENOTYPE_LEVELS*x)];

            float expected  = ( c1 * c2 ) / c3 ;

            retVal = retVal +  (pow(c4-expected,2)/expected);
           
        }
    }

    return retVal;
}

array<vector<PACK_TYPE>, GENOTYPE_LEVELS> Snp::packGenotypesHelper(const vector<char> & genotypes)
{
    array<vector<PACK_TYPE>, GENOTYPE_LEVELS> packedGenotypes;

    array<PACK_TYPE,GENOTYPE_LEVELS> sectionCode;
    sectionCode.fill(0);
    for(int i=0; i<genotypes.size(); ++i)
    {
  
        PACK_TYPE mask = (PACK_TYPE)1<<(PACK_SIZE-(i+1)%PACK_SIZE);
        
        //Also convert from ascii value to numerical
        switch(genotypes[i])
        {
            case '0':
                sectionCode[0] = sectionCode[0] | mask;
                break;
            case '1':
                sectionCode[1] = sectionCode[1] | mask;
                break;
            case '2':
                sectionCode[2] = sectionCode[2] | mask;
                break;
        }
        
        //When a section is full, push them back and reset the sections to empty;
        if( ((i+1) % PACK_SIZE == 0) || ( (i+1) == genotypes.size()) )
        {
            for(int l = 0; l < GENOTYPE_LEVELS; l++)
                packedGenotypes[l].push_back(sectionCode[l]);
            sectionCode.fill(0);
        }
    }

    return packedGenotypes;
}

void Snp::packGenotypes(vector<PACK_TYPE> & allPackedGenotypes, const vector<char> & controlGenotypes, const vector<char> & caseGenotypes)
{
    array<vector<PACK_TYPE>, GENOTYPE_LEVELS> packedControls = packGenotypesHelper(controlGenotypes);
    array<vector<PACK_TYPE>, GENOTYPE_LEVELS> packedCases = packGenotypesHelper(caseGenotypes);

    for(int coding=0; coding<GENOTYPE_LEVELS; ++coding)
        for(int i=0; i< packedControls[0].size(); ++i)
        {
            allPackedGenotypes.push_back(packedControls[coding][i]);
        }

    for(int coding=0; coding<GENOTYPE_LEVELS; ++coding)
        for(int i=0; i< packedCases[0].size(); ++i)
        {
            allPackedGenotypes.push_back(packedCases[coding][i]);
        }
}

