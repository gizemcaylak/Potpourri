/**
 * @author Tyler Cowman
 *  
 * Acts as a wrapper for a collection of LD-Groups and associated functionality
 * such as merging and testing the set of groups.
 */


#ifndef LDFOREST_H
#define LDFOREST_H

#include "CommonStructs.h"
#include "Snp.h"
#include "LDGroup.h"
#include "TopSnpList.h"
#include "FilePath.h"

#include <limits.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <bits/stdc++.h> 

using namespace std;
class LDForest
{
    public:
        LDForest(int topKSnps, int numberControlSamples, int numberCaseSamples, int numberSnps);
        
        void insert(const Snp & snp, char chromosome, int basePair);
    
        int size()const;
        void createGroups(vector<vector<int>> clusterIndices);
        void testGroups(int maxThreadUsage,  ParameterInfo parameterInfo);
        TopSnpList writeResults(string fileName, ParameterInfo parameterInfo, DatasetSizeInfo datasetSizeInfo);
    
    private:        
        //Results output
        void writeGroundTruthList();
        
        vector<LDGroup> ldgroups_;
        
        TopSnpList topSnpList_;
        
        FilePath outputDirectory_;
};
#endif //LDFOREST_H