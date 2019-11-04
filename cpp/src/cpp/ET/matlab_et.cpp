/*
* @author Gizem Caylak
* This file is part of Potpourri
*
* Potpourri is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Potpourri is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mex.h"
#include "matrix.h"
#include "CommonStructs.h"
#include "Matrix.h"
#include "MatrixMath.h"
#include "Snp.h"
#include "TopSnpList.h"
#include "FilePath.h"
#include "LDForest.h"

#include <string>
#include <vector>
#include <iostream>
#include <omp.h>
#include <fstream>
#include <map>
#include <ctime>


const static array<float, 7> chi2DegreesFreedomTable{0, 2.71, 6.63, 10.82, 15.14, 19.57, 24.87};
/*
	Input								Type						Size						
	------------------------------		-------						-------------------
	Case Features (X_1)					double						c x n matrix
	Control Features (X_2)				double						m-c x n matrix
	SNP info (I)						string						n x 3 sparse matrix
	Group ids (G) 						integer						vector
	Weights (W)							logical						vector
    OutputFileName (filename)           string                      1
    MaxMarginalSignificance             float                       1
*/


void displayUsageMessage() {
	mexPrintf("Usage: epistasis_test(\n \
		CaseFeatures X_1<Vector, cxn> ,\n \
		ControlFeatures X_2<Vector, (m-c)xn> ,\n \
		SNPInfo I <Matrix, nx3>,\n \
		Group ids G<vector>,\n \
		Weights W<vector>) \n \
        OutputFileName fileName\n \
        MaxMarginalSignificance maxMarginalSignificance \
		\n");
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double* w_ = mxGetPr(prhs[4]);
	double* g_ = mxGetPr(prhs[3]);
	mwSize nRowWeights = mxGetM(prhs[4]);
	std::vector<double> weights;
	weights.reserve(nRowWeights);
	std::vector<double> regionIndices;
	regionIndices.reserve(nRowWeights);
	for (mwSize i = 0; i < nRowWeights; i++) {
		weights.push_back(w_[i]);
		regionIndices.push_back(g_[i]);
	}

	// Read in the data
    Matrix<string> infoMatrix;
	mwSize col_outputBuff     = mxGetN(prhs[2]);
	mwSize row_outputBuff     = mxGetM(prhs[2]);

    const mxArray *cell_element_ptr;
    char* c_array;
    mwSize buflen;
    int status;
    // Read in the data
	for (mwSize row=0; row < row_outputBuff; row++) {
		vector<string> parsedRow;
		for (mwSize col=0; col < col_outputBuff; col++) {
			cell_element_ptr = mxGetCell(prhs[2],row+col*row_outputBuff);
			buflen = mxGetN(cell_element_ptr)*sizeof(mxChar)+1;
		    c_array = (char*)mxMalloc(buflen);
		    status = mxGetString(cell_element_ptr,c_array,buflen);
 			parsedRow.push_back(std::string(c_array));
    		mxFree(c_array);
	  	}
	  	infoMatrix.addRow(parsedRow);
	}

	Matrix<char> casesMatrix; 
	Matrix<char> controlsMatrix;

	double* caseBuff = mxGetPr(prhs[0]);
	col_outputBuff     = mxGetN(prhs[0]);
	row_outputBuff     = mxGetM(prhs[0]);
	// Read in the data
	for (mwSize row=0; row < row_outputBuff; row++) {
		vector<char> parsedRow;
		for (mwSize col=0; col < col_outputBuff; col++) {
			parsedRow.push_back(((int)caseBuff[row+col*row_outputBuff])+'0');
      }
	  casesMatrix.addRow(parsedRow);
	}

	double* controlBuff = mxGetPr(prhs[1]);
	col_outputBuff     = mxGetN(prhs[1]);
	row_outputBuff     = mxGetM(prhs[1]);
	// Read in the data
	for (mwSize row=0; row < row_outputBuff; row++) {
		vector<char> parsedRow;
		for (mwSize col=0; col < col_outputBuff; col++) {
			parsedRow.push_back(((int)controlBuff[row+col*row_outputBuff])+'0');
	  	}
	  controlsMatrix.addRow(parsedRow);
	}


    //If any of the files were read incorrectly exit
    if(infoMatrix.size() == 0 || casesMatrix.size() == 0 || controlsMatrix.size() == 0)
        exit(-1);
    //Recored the initial size of the dataset read in
    DatasetSizeInfo datasetSizeInfo;
    datasetSizeInfo.snps_ = infoMatrix.dim(0);
    datasetSizeInfo.cases_ = casesMatrix.dim(1);
    datasetSizeInfo.controls_ = controlsMatrix.dim(1);
    datasetSizeInfo.marginalSignificanceRemoved_ = 0;
    datasetSizeInfo.filterFileRemoved_ = 0;
    ParameterInfo parameterInfo;
    parameterInfo.snpBeginIndex_ = 0;
    parameterInfo.snpEndIndex_ = infoMatrix.dim(0);
    parameterInfo.maxUnknownFraction_ = 0;
    parameterInfo.maxThreadUsage_ = 20;
    int topK = 2000;

    int maxMarginalSignificance = (int)mxGetScalar(prhs[6]); 
    parameterInfo.maxMarginalSignificance_ = maxMarginalSignificance;

    char *filename;
    string outputFileName;
     /* Find out how long the input string array is. */
    col_outputBuff     = mxGetN(prhs[5]);
    row_outputBuff     = mxGetM(prhs[5]);
    // Read in the data
    for (mwSize row=0; row < row_outputBuff; row++) {
        for (mwSize col=0; col < col_outputBuff; col++) {
            cell_element_ptr = mxGetCell(prhs[5],row+col*row_outputBuff);
            buflen = mxGetN(cell_element_ptr)*sizeof(mxChar)+1;
            filename = (char*)mxMalloc(buflen);
            status = mxGetString(cell_element_ptr,filename,buflen);
            outputFileName=std::string(filename);
            mxFree(filename);
        }
    }
    parameterInfo.outputFileName_ = outputFileName;
    //Call snp constructors to create bitwise snp representations
    vector<Snp> snps;
    
    for(int i= parameterInfo.snpBeginIndex_; ceil(i)<min( parameterInfo.snpEndIndex_, infoMatrix.dim(0)); ++i)
    {
        snps.push_back(Snp(i, controlsMatrix.getRow(i), casesMatrix.getRow(i), weights[i]));
    }
    

    LDForest ldforest( infoMatrix.dim(0) , controlsMatrix.dim(1), casesMatrix.dim(1), infoMatrix.dim(0));

    //Only create from snps with low enough marginal significance
    for(int i=0; i<snps.size(); ++i)
    {  
        if(snps[i].marginalTest() > chi2DegreesFreedomTable[parameterInfo.maxMarginalSignificance_])
        {
            datasetSizeInfo.marginalSignificanceRemoved_++;
            regionIndices[i] = -1;
        }
        else
        {
            if(isdigit(infoMatrix.a(snps[i].getIndex(),1)[0]));
            {
                ldforest.insert(snps[i], (char)stoi(infoMatrix.a(snps[i].getIndex(),1)), stoi(infoMatrix.a(snps[i].getIndex(),2))); 
            }   
        }    
    }

    datasetSizeInfo.passingSnps_ = ldforest.size();
    cout<<"---SNPs in LDForest = "<<ldforest.size()<<endl;
    cout<<"\tSNPs Removed Marginal Significance: "<<datasetSizeInfo.marginalSignificanceRemoved_<<endl;
    // assign SNPs to the regions they belong
    vector<int> tmpRegionInd;
    tmpRegionInd.reserve(ldforest.size());
    for(int i = 0 ; i < snps.size(); i++) {
    	if(regionIndices[i]!=-1)
    		tmpRegionInd.push_back(regionIndices[i]);
    }

    vector<vector<int>> snp_ind_per_region;
    int no_of_regions = tmpRegionInd[tmpRegionInd.size()-1];
    snp_ind_per_region.resize(no_of_regions);
    
    for(int i = 0 ; i < ldforest.size(); i++) {
    	snp_ind_per_region[tmpRegionInd[i]-1].push_back(i);
    }
    TopSnpList topSnpList_;
    if(ldforest.size() > 1)
    {   
        ldforest.createGroups(snp_ind_per_region);
        ldforest.testGroups(parameterInfo.maxThreadUsage_, parameterInfo);
        topSnpList_ = ldforest.writeResults(parameterInfo.outputFileName_, parameterInfo, datasetSizeInfo);
    }
    
    topSnpList_.calculateFormattedResults();
    vector<TopPairing> topPairs = topSnpList_.getReciprocalPairs();
    
    ofstream ofs("output/" + parameterInfo.outputFileName_ + ".reciprocalPairs.formatted");
   
    for(int j=min((size_t)topK, topPairs.size() )-1; j>=0; --j )
    {
        ofs<<topPairs[j].score_<<"\t"<<infoMatrix.a(topPairs[j].indexes_.first, 0)<<"\t"<<infoMatrix.a(topPairs[j].indexes_.second, 0)<<"\t"
                        <<infoMatrix.a(topPairs[j].indexes_.first, 1)<<"\t"<<infoMatrix.a(topPairs[j].indexes_.second, 1)<<"\t"
                        <<infoMatrix.a(topPairs[j].indexes_.first, 2)<<"\t"<<infoMatrix.a(topPairs[j].indexes_.second, 2)<<"\t"<<endl;
    }

    ofs.close();
}