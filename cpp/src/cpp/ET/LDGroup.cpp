#include "LDGroup.h" 

LDGroup::LDGroup(const Snp & snp , char chromosome, int basePair)
{
   nodes_.push_back(snp);
   genomeLocations_.push_back(GenomeLocation(chromosome, basePair));
}

LDGroup::LDGroup(const LDGroup & cpy)
{
    nodes_ = cpy.nodes_;
    genomeLocations_ = cpy.genomeLocations_;
}
//create LD GROUP
LDGroup::LDGroup(const LDGroup & t1, const LDGroup & t2)
{
    nodes_.reserve(t1.size() + t2.size());
    genomeLocations_.reserve(t1.size() + t2.size());
    
    for(int i=0; i<t1.size(); i++) {
        nodes_.push_back(t1.nodes_[i]);
        genomeLocations_.push_back(t1.genomeLocations_[i]);
    }
    for(int i=0; i<t2.size(); i++) {
        nodes_.push_back(t2.nodes_[i]);
        genomeLocations_.push_back(t2.genomeLocations_[i]);
    }
}

bool LDGroup::empty()const
{
    return nodes_.empty();
}

int LDGroup::size()const
{
    return nodes_.size();
}

const Snp & LDGroup::getRoot()const
{
    return nodes_[0];
}

int LDGroup::computeDifferences(const LDGroup & other)const
{
    return nodes_[0].computeDifferences(other.nodes_[0]);
}

void LDGroup::clear()
{
    nodes_.clear();
}

void LDGroup::epistasisTest(const LDGroup & other, TopSnpList & topSnpList)const
{
   
    int top_k = std::min(10, (int)(other.nodes_.size()*nodes_.size()));
    long long int localLeaftTestsDone = 0;

    PopcoverQueue popcoverQueue;
    for (unsigned int i = 0; i < nodes_.size(); i++) {
        for(unsigned int j = 0; j < other.nodes_.size(); j++) {
            if((nodes_[i].getWeight()) && (other.nodes_[j].getWeight()))
                popcoverQueue.push(HeapData((nodes_[i].computePopCoverAnd(other.nodes_[j])+nodes_[0].getCaseNo()), i, j));
            else
                popcoverQueue.push(HeapData(nodes_[i].computePopCoverAnd(other.nodes_[j]), i, j));
        }
    }

    // fill the stack based on pop cover 
    for(int i = 0; i < top_k; i++) {
            auto a = popcoverQueue.top();
            int Pmax = a.key;
            unsigned int maxIndex1 = a.index1;
            unsigned int maxIndex2 = a.index2;
            popcoverQueue.pop();

            float score = nodes_[maxIndex1].epistasisTest(other.nodes_[maxIndex2]);
            //check to make sure not estimated as being in LD
            if(!genomeLocations_[maxIndex1].inLinkageDisequilibrium(other.genomeLocations_[maxIndex2]) )
            {
                topSnpList.attemptInsert(nodes_[maxIndex1].getIndex(), other.nodes_[maxIndex2].getIndex(), score);
                localLeaftTestsDone += 2;
            }
    }
    topSnpList.incrementLeafTestsCounter(localLeaftTestsDone);
    
}

