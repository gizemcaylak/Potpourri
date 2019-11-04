#include "TopSnpList.h"

TopSnpList::TopSnpList()
{
    
}

TopSnpList::TopSnpList(const TopSnpList & cpy)
{
    topK_ = cpy.topK_;
    cutoff_ = cpy.cutoff_;
    
    topScores_ = cpy.topScores_;
    topPartners_ = cpy.topPartners_;
    internalTestsCounter_ = cpy.internalTestsCounter_;
    leafTestsCounter_ = cpy.leafTestsCounter_;
    
    pairwiseSignificanceCounts_ = cpy.pairwiseSignificanceCounts_;
    pairwiseSignificanceCountsPrefixSum_ = cpy.pairwiseSignificanceCountsPrefixSum_;
    prefixSumTimer_ = cpy.prefixSumTimer_;
}

TopSnpList::TopSnpList(int topK, int numberSnps, float cutoff)
{
    topK_ = topK;
    cutoff_ = cutoff;
    
    topScores_ = vector<float>(numberSnps, 0.0);
    topPartners_ = vector<int>(numberSnps, -1);
    internalTestsCounter_ = 0;
    leafTestsCounter_ = 0;
    
    pairwiseSignificanceCounts_.fill(0);
    pairwiseSignificanceCountsPrefixSum_.fill(0);
    prefixSumTimer_ = 0;
}

bool TopSnpList::attemptInsert(int snpIndex1, int snpIndex2, float score)
{ 
    bool retVal = false;

    #pragma omp critical
    {

        if( prefixSumTimer_ >= PREFIX_SUM_ROLLOVER )
        {
            
            pairwiseSignificanceCountsPrefixSum_ = pairwiseSignificanceCounts_;
            for(int i=MAX_CUTOFF; i>getCutoff(); --i)
            {
                pairwiseSignificanceCountsPrefixSum_[i-1] = pairwiseSignificanceCountsPrefixSum_[i-1]+pairwiseSignificanceCounts_[i];
            }
            
            prefixSumTimer_ = 0;
        }
        
        if(score > getCutoff())
        {
            int index = min((int)score, MAX_CUTOFF);
            pairwiseSignificanceCounts_[index]++;
            
            
            
            if(pairwiseSignificanceCounts_[index] >= topK_ || pairwiseSignificanceCountsPrefixSum_[index] >= topK_)
                cutoff_ = index;
        } 
        
        if(score > topScores_[snpIndex1])
        {
            topScores_[snpIndex1] = score;
            topPartners_[snpIndex1] = snpIndex2;
            retVal = true;
        }
        
        if(score > topScores_[snpIndex2])
        {
            topScores_[snpIndex2] = score;
            topPartners_[snpIndex2] = snpIndex1;
            retVal = true;
        }
        prefixSumTimer_++;
    }

    return retVal;
}

float TopSnpList::getCutoff()const
{
    return cutoff_;
}

void TopSnpList::incrementInternalTestsCounter(long long int testsDone)
{
    #pragma omp critical
    internalTestsCounter_ += testsDone;
}

void TopSnpList::incrementLeafTestsCounter(long long int testsDone)
{
    #pragma omp critical
    leafTestsCounter_ += testsDone;
}

long long int TopSnpList::getInternalTests()const
{
    return internalTestsCounter_;
}

long long int TopSnpList::getLeafTests()const
{
    return leafTestsCounter_;
}

void TopSnpList::calculateFormattedResults()
{
    
    for(int i=0; i< topPartners_.size(); ++i)
        if(topScores_[i] >= getCutoff()-1 && topScores_[i] >0)
            cutoffPairs_.push_back(  TopPairing(i, topPartners_[i], topScores_[i] )       );

    
    calculateReciprocalPairs();
    
    sort(reciprocalPairs_.begin(), reciprocalPairs_.end(), TopPairing::orderByScore);
    sort(cutoffPairs_.begin(), cutoffPairs_.end(), TopPairing::orderByScore);
}

float TopSnpList::calculateRecall(string groundTruthFileName)
{
    ifstream ifs(groundTruthFileName);
    if(!ifs.is_open())
        return -2.0;
    
    //Holds the pairs that should exist
    map<TopPairing, int> groundTruth;
    int tempIndex1; int tempIndex2; float tempScore;
    while(true)
    {
      
        ifs>>tempScore;
        if(ifs.eof())
            break;
        ifs>>tempIndex1;
        ifs>>tempIndex2;
        groundTruth.insert(make_pair(TopPairing(tempIndex1, tempIndex2, tempScore), 0));
    }
    ifs.close();

    int foundCounter = 0;
    for(int i=0; i<reciprocalPairs_.size(); ++i)
    {
        if(groundTruth.count(reciprocalPairs_[i]) == 1)
            foundCounter++;
    }
    
    return foundCounter/(float)groundTruth.size();
}

float TopSnpList::calculatePrecision(string groundTruthFileName)
{
    ifstream ifs(groundTruthFileName);
    if(!ifs.is_open())
        return -2.0;
    
    //Holds the pairs that should exist
    map<TopPairing, int> groundTruth;
    int tempIndex1; int tempIndex2; float tempScore;
    while(true)
    {
        ifs>>tempScore;
        if(ifs.eof())
            break;
        ifs>>tempIndex1;
        ifs>>tempIndex2;
        groundTruth.insert(make_pair(TopPairing(tempIndex1, tempIndex2, tempScore), 0));
    }
    ifs.close();

    int foundCounter = 0;
    for(int i=0; i<reciprocalPairs_.size(); ++i)
    {
        if(groundTruth.count(reciprocalPairs_[i]) == 1)
            foundCounter++;
    }
    
    return foundCounter/(float)reciprocalPairs_.size();
}

float TopSnpList::calculateRecallCutoffPairs(string groundTruthFileName)
{
    ifstream ifs(groundTruthFileName);
    if(!ifs.is_open())
        return -2.0;
    
    //Holds the pairs that should exist
    map<TopPairing, int> groundTruth;
    int tempIndex1; int tempIndex2; float tempScore;
    while(true)
    {
      
        ifs>>tempScore;
        if(ifs.eof())
            break;
        ifs>>tempIndex1;
        ifs>>tempIndex2;
        groundTruth.insert(make_pair(TopPairing(tempIndex1, tempIndex2, tempScore), 0));
    }
    ifs.close();
    
    
    //for(auto it=groundTruth.begin(); it!=groundTruth.end(); ++it)
      //  cout<<it->first<<endl;
    
    int foundCounter = 0;
    for(int i=0; i<cutoffPairs_.size(); ++i)
    {
        //Only count each pairing once
        if(groundTruth.count(cutoffPairs_[i]) == 1)
            if(groundTruth[cutoffPairs_[i]] == 0)
            {
                groundTruth[cutoffPairs_[i]] = 1;
                foundCounter++;
            }    
    }
    
    return foundCounter/(float)groundTruth.size();
}

void TopSnpList::printReciprocalPairs(ofstream & ofs)const
{
    for(int i=0; i<reciprocalPairs_.size(); ++i)
        ofs<<reciprocalPairs_[i]<<endl;
}

vector<TopPairing> TopSnpList::getReciprocalPairs()const {
    return reciprocalPairs_;
}
void TopSnpList::printCutoffPairs(ofstream & ofs)const
{
    for(int i=0; i<cutoffPairs_.size(); ++i)
        ofs<<cutoffPairs_[i]<<endl;
}

void TopSnpList::printReciprocalPairScores(ofstream & ofs)const
{
    for(int i=0; i<reciprocalPairs_.size(); ++i)
        ofs<<reciprocalPairs_[i].score_<<"\t";
    ofs<<endl;
    for(int i=0; i<reciprocalPairs_.size(); ++i)
        ofs<<reciprocalPairs_[i].indexes_.first<<"\t";
    ofs<<endl;
    for(int i=0; i<reciprocalPairs_.size(); ++i)
        ofs<<reciprocalPairs_[i].indexes_.second<<"\t";
    ofs<<endl;
}

void TopSnpList::printCutoffPairScores(ofstream & ofs)const
{
    for(int i=0; i<cutoffPairs_.size(); ++i)
        ofs<<cutoffPairs_[i].score_<<"\t";
    ofs<<endl;
    for(int i=0; i<cutoffPairs_.size(); ++i)
        ofs<<cutoffPairs_[i].indexes_.first<<"\t";
    ofs<<endl;
    for(int i=0; i<cutoffPairs_.size(); ++i)
        ofs<<cutoffPairs_[i].indexes_.second<<"\t";
    ofs<<endl;
}

int TopSnpList::getNumberOfReciprocalPairs()const
{
    return reciprocalPairs_.size();
}

ostream& operator<< (ostream &out, const TopSnpList & topSnpList)
{
    vector< pair<float, pair<int,int> > > passingPairs;
    
    
    for(int i=0; i< topSnpList.topPartners_.size(); ++i)
    {
        
        if(topSnpList.topScores_[i] >= topSnpList.getCutoff()-1)
        {
            //cout<< topSnpList.topScores_[i]<<endl;
            passingPairs.push_back(   make_pair( topSnpList.topScores_[i], make_pair(i, topSnpList.topPartners_[i])       )       );
        }


        
    }
    
    sort(passingPairs.rbegin(), passingPairs.rend());

    map<TopPairing, int> reciprocalCheck;
    
    
    
    
    //for(int i=0; i< topSnpList.topK_; ++i)
    for(int i=0; i< passingPairs.size(); ++i)
    {   
        TopPairing tempPair(passingPairs[i].second.first, passingPairs[i].second.second, passingPairs[i].first);
        
        if(reciprocalCheck.count(tempPair)>0)
            reciprocalCheck[tempPair]++;
        else
        {
            reciprocalCheck.insert(make_pair(tempPair, 1));
            //reciprocalCheck[tempPair]=1;
        } 
    }
    
    vector<float> lazySort;
    int recips =0;
    for(auto it=reciprocalCheck.begin(); it!=reciprocalCheck.end(); ++it)
    {
        if(it->second == 2)
        {
            
            recips++;
            //cout<<it->first.score_<<endl;
            lazySort.push_back(it->first.score_);
        }
    }
    
    sort(lazySort.rbegin(), lazySort.rend());
    for(int i=0; i< lazySort.size(); ++i)
        cout<<lazySort[i]<<endl;
    
    
    out<<"INTERNAL TESTS:   "<<topSnpList.internalTestsCounter_<<endl;
    out<<"LEAF TESTS:       "<<topSnpList.leafTestsCounter_<<endl;
    out<<"PASSING:          "<<passingPairs.size()<<endl;
    out<<"RECIPROCALS:      "<<recips<<endl;
    
    return out;
}

void TopSnpList::calculateReciprocalPairs()
{
    map<TopPairing, int> reciprocalCheck;
    for(int i=0; i< topPartners_.size(); ++i)
    {
        if(topScores_[i] >= getCutoff()-1)
        {
            TopPairing tempPair(i, topPartners_[i], topScores_[i]);
            
            if(reciprocalCheck.count(tempPair)>0)
                reciprocalCheck[tempPair]++;
            else
                reciprocalCheck.insert(make_pair(tempPair, 1));
        }
    }
    
    for(auto it=reciprocalCheck.begin(); it!=reciprocalCheck.end(); ++it)
    {
        if(it->second == 2)
            reciprocalPairs_.push_back(it->first);
    }
}