function demo_potpourri
load('data/data.mat'); 	% load feature matrix (X), phenotypes (Y), SNP information(SNP_info), regulatory/coding information (R)
load('data/networks.mat');		% load network information
phenotype_index = 5;			% Index of phenotype selected
k = 1000;						% Number of features to be selected
omega = 1;						% Omega parameter of Potpourri to reward regulatory region
b = 9; 							% number of neighbors (left:9 and right:9) total 18
GS_network = networks{1};		% use GS network
outputFileName = "toy";			% Prefix for the output files
maxMarginalSignificance = 6;	% integer value between 1-6 maximum marginal significance of loci
potpourri( X, Y(:,phenotype_index), GS_network, k, R, SNP_info, b, omega, maxMarginalSignificance, outputFileName);
end

