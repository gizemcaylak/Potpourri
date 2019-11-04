/**
 * @author Tyler Cowman
 * 
 * Class representing a single LD-Group. The SNP nodes are seperated from the genome locations 
 * to reduce the memory footprint when calculating contingecy tables.
 */

#ifndef LDGROUP_H
#define LDGROUP_H

#include "Snp.h"
#include "TopSnpList.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <queue>
#include <utility>
#include <functional>
#include <boost/heap/binomial_heap.hpp>
#include <queue>
#include <ctime>
#include <cstdlib>

#define ESTIMATED_LD_RANGE 1000000

struct GenomeLocation
{
    GenomeLocation(char chromosome, int basePair)
    {
        chromosome_ = chromosome;
        basePair_ = basePair;
    }
    
    GenomeLocation(const GenomeLocation & l1, const GenomeLocation & l2)
    {
        chromosome_ = '!';
        basePair_ = -1;
    }
    
    bool inLinkageDisequilibrium(const GenomeLocation & other)const
    {
        if(chromosome_ != other.chromosome_)
            return false;
        else if(abs(basePair_ - other.basePair_) > ESTIMATED_LD_RANGE)
            return false;
        else
            return true;
            
    }
    
    char chromosome_;
    int basePair_;
};

struct HeapData {
    int key; // pop cover
    unsigned int index1;
    unsigned int index2;
    HeapData(int _key, unsigned int _index1, unsigned int _index2) : key(_key), index1(_index1), index2(_index2) { }
    bool operator<(HeapData const & rhs) const {
        return key == rhs.key ? index1 < rhs.index1 : key < rhs.key;
    }
};

using namespace std;
class LDGroup
{
    public:
        
        LDGroup(const Snp & snp, char chromosome, int basePair);
        LDGroup(const LDGroup & cpy);
        LDGroup(const LDGroup & t1, const LDGroup & t2);
    
        bool empty()const;
        int size()const;
        
        const Snp & getRoot()const;
        
        int computeDifferences(const LDGroup & other)const;
        
        void clear();
        
        void epistasisTest(const LDGroup & other, TopSnpList & topSnpList)const;

        friend ostream& operator<< (ostream &out, const LDGroup & ldgroup);
        typedef boost::heap::binomial_heap<HeapData> PopcoverQueue;
        typedef PopcoverQueue::handle_type PopcoverHandle;

        vector<Snp> nodes_;
        vector<GenomeLocation> genomeLocations_;

};
#endif //LDGROUP_H
