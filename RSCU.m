function [RSCU_table] = RSCU(ConcatenatedCDSs)
%RSCU is a function that has one input: 
%|CDSsconcat| - a string of concatenated CDSs 

%RSCU is a function that has one output: 
%|RSCU_table| - a 59x2 table where the first column holds each codon (excluding
%TGG, TGA, TAG, TAA, ATG) and the second column holds the RSCU value of the
%codon 

%RSCU:
%1) Creates a map object which associates values with corresponding unique
%   keys - the keys in the created map object are the one letter abbreviation
%   of amino acids (excluding Methionine (M) and Tryptophan (W)), the values
%   in the created map are the codons which encode each amino acid
%2) Next, the function counts the number of occurrences of each codon in
%   the concatenated CDS string (AllCDSconcat) but does not retain the count
%   of codons (TGG, TGA, TAG, TAA and ATG)
%3) The function declares variables:
%   TotalObserved - For each codon, TotalObserved holds the total number 
%                   of occurrences of synonymous codons which encode the 
%                   same amino acid 
%   Expected_under_noCB - TotalObserved divided by the number of synonymous
%                         codons 
%   RSCU - to hold the RSCU values of the 59 codons 
%   codon - to hold the name of each codon 
%4) Next, the function loops through each of the 59 codons and computes the
%   RSCU value of the codon through a series of steps: 
%   i) translates each codon to an amino acid using nt2aa ignoring alternative
%      start codons
%   ii) obtaining the synonymous codons that encode the amino acid
%   iii) works out the number of synonymous codons for the amino acid
%   iv) works out the total number of occurrences of all synonymous codons
%   v) divides this value by the number of synonymous codons to obtain the 
%      expected occurrence of each codon 
%   vi) divides the observed occurrence of each codon (obtained from
%       codoncount) by the expected occurrence of each codon to get the RSCU
%       value of the codon 


%% 1) Creates a map object 
%keySet is a set of unique 'key' elements (i.e. the amino acids) - each key
%is mapped to a correpsonding value from the valueSet
%valueSet is the set of codons that encode each amino acid
Map = rmfield(revgeneticcode(1),{'Name','Starts','Stops','M','W'});
%revgeneticcode is a built-in MATLAB structure that contains the reverse
%mapping of amino acids to codons (default genetic code = 1 (Standard))
%rmfield removes the field 'Name','Starts','Stops','M' and 'W' from the structure revgeneticcode
keySet = fieldnames(Map); %amino acids are the fieldnames of Map defined above
valueSet = struct2cell(Map); %a cell array of codons that encode each amino acid
M = containers.Map(keySet, valueSet, 'UniformValues', false); %creates a Map object

%% 2) Counts for a nucleotide sequence 
%codoncount counts the number of each codon in a nucleotide sequence
CC = codoncount(ConcatenatedCDSs);
CC = rmfield(CC, {'TGG','TGA','TAG','TAA','ATG'});

%% 3) Declaring variables
Total_Observed = zeros(59,1); %to hold the total observed occurences of synonymous codons for each amino acid - Total observed values will be the same for synonymous codons
Expected_under_no_CB = zeros(59,1); %to hold expected occurence of each codon if no codon bias were to exist
RSCU = zeros(59,1); %to hold RSCU values for each codon 
codon = fieldnames(CC); %to hold each codon

%% 4) Main FOR loop
for i=1:59 %loops over each codon 
    
    specific_codon = string(codon(i,1)); 
    key = cellstr(nt2aa(specific_codon,'GeneticCode',1,'AlternativeStartCodons', false)); %works out what aa the codon encodes

    syn_codons = values(M,key); %returns all the codons which encode the aa given by the variable key
    NumberOfCodons = width(syn_codons{1,1}); %computes the number of synonymous codons for the aa 

    %Computes the total observed occurrences of synonymous codons
    for j = 1:NumberOfCodons
        Total_Observed(i,1) = Total_Observed(i,1) + CC.(syn_codons{1,1}{1,j});
    end

    %Computes expected frequency of the specific codon under no codon bias 
    Expected_under_no_CB(i,1) = Total_Observed(i,1)/NumberOfCodons;

    %Computes RSCU value for the specific codon 
    RSCU(i,1) = CC.(specific_codon)/Expected_under_no_CB(i,1);
    
    if isnan(RSCU(i,1))
        RSCU(i,1) = 0;
    end 
end

RSCU = num2cell(RSCU);
RSCU_table = [codon RSCU];
RSCU_table = cell2table(RSCU_table);
RSCU_table.Properties.VariableNames = {'codon', 'RSCUvalues'};



end

