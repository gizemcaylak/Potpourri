#include "LDForest.h"

LDForest::LDForest( int topKSnps, int numberControlSamples, int numberCaseSamples, int numberSnps)
{
    topSnpList_ = TopSnpList(topKSnps, numberSnps, 0);
    
    outputDirectory_ = FilePath();
    outputDirectory_.moveDown("output");
}

void LDForest::insert(const Snp & snp, char chromosome, int basePair)
{
    ldgroups_.push_back(LDGroup(snp, chromosome, basePair));
}

int LDForest::size()const
{
    return ldgroups_.size();
}

void LDForest::createGroups(vector<vector<int>> clusterIndices)
{
    cout<<"---Creating LD Groups"<<endl;
    vector<LDGroup> createdGroups;
    for(int c = 0; c < clusterIndices.size(); c++) {
        if(clusterIndices[c].size() > 1) {
            LDGroup prevGroup(ldgroups_[clusterIndices[c][0]], ldgroups_[clusterIndices[c][1]]);
            for(int i = 2; i < clusterIndices[c].size(); i++) {
                LDGroup potentialGroup(prevGroup, ldgroups_[clusterIndices[c][i]]);
                prevGroup.nodes_ = potentialGroup.nodes_;
                prevGroup.genomeLocations_ = potentialGroup.genomeLocations_;
            }
            LDGroup finalgroup(prevGroup);
            {
                createdGroups.push_back(finalgroup);
            }
        }
        else if(clusterIndices[c].size() == 1) {
            createdGroups.push_back(ldgroups_[clusterIndices[c][0]]);
        }
    }
    
    ldgroups_ = createdGroups;
    cout<<"The number of LD Groups: "<<size()<<endl;
}

void LDForest::testGroups(int maxThreadUsage,  ParameterInfo parameterInfo)
{
    cout<<"---Testing Pairs"<<endl;
    cout<<"\tFinished 0/"<<size()<<"               \r"<<flush;
    
    int testsFinished = 0;
    #pragma omp parallel for num_threads(maxThreadUsage) schedule(dynamic, 20)
    for(int i=0; i<size(); ++i)
    {
		
        for(int j=i+1; j<size(); ++j)
            ldgroups_[i].epistasisTest(ldgroups_[j], topSnpList_);
        
        if(i % 100 == 0)
        {   
            #pragma omp critical
            {
            testsFinished +=100;
            cout<<"\tFinished "<<min(testsFinished, size())<<"/"<<size()<<"               \r"<<flush;
            }
                
        }
            
    }
    cout<<endl;
}

TopSnpList LDForest::writeResults(string fileName, ParameterInfo parameterInfo, DatasetSizeInfo datasetSizeInfo)
{
    topSnpList_.calculateFormattedResults();
    
    outputDirectory_.makeActive();
    ofstream ofs;
    
    if(parameterInfo.maxUnknownFraction_ == 0.0)
    {
        ofs.open("reciprocalPairs.gt");
        topSnpList_.printReciprocalPairs(ofs);
        ofs.close();
        
        ofs.open("cutoffPairs.gt");
        topSnpList_.printCutoffPairs(ofs);
        ofs.close();
    }
    
    //Check if the summary file is currently empty
    bool needsHeader = false;
    string tempLine;
    ifstream ifs(fileName + ".summary");
    getline(ifs, tempLine);
    if(tempLine.size() == 0)
        needsHeader = true;
    ifs.close();
    
    ofs.open(fileName + ".summary", ofstream::out | ofstream::app);
    
    if(needsHeader)
        ofs <<"output_file\tsnps_input\tcase_samples\tcontrol_samples\tmax_marginal_significance\tsnps_ignored_marginal_significance\tno_of_tests\tpairs_reported"<<endl;

    ofs<<parameterInfo.outputFileName_<<"\t";
    datasetSizeInfo.printDimensions(ofs); ofs<<"\t";
    parameterInfo.printSummaryRelevant(ofs); ofs<<"\t";
    datasetSizeInfo.printSubsetDimensions(ofs); ofs<<"\t";

    long long totalTestsPerformed = (topSnpList_.getLeafTests())/2; //Divide by two because each test was counted in each 
    ofs<<totalTestsPerformed<<"\t";
    
    ofs<<topSnpList_.getNumberOfReciprocalPairs();
    ofs<<endl;
    ofs.close();
    
    ofs.open(fileName + ".reciprocalPairs", ofstream::out | ofstream::app);
    topSnpList_.printReciprocalPairScores(ofs);
    ofs.close();
    
    ofs.open(fileName + ".cutoffPairs", ofstream::out | ofstream::app);
    topSnpList_.printCutoffPairScores(ofs);
    ofs.close();
    
    outputDirectory_.returnActive();
    
    cout<<"---Done"<<endl;
    cout<<"\tLeaf Tests: "<<topSnpList_.getLeafTests()<<endl;
    cout<<"\tInternal Tests: "<<topSnpList_.getInternalTests()<<endl;
    return topSnpList_;
}