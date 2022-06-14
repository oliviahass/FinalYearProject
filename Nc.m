function [Nc] = Nc(Seq)
%Nc is a function that has one input: 
%|Seq| - a CDS in the form of a string 

%Nc is a function that has one output: 
%|Nc| - the 'effective number of codons' used in a gene 


%Nc: 
%What the funtion does: Nc calculates the effective number of codons - 
%the effective number of codons, Nc, measures the extent to which codon 
%usage of a gene departsfrom equal usage of synonymous codons (Wright 1990). 
%Nc values range from 20 to 61 - the closer the Nc value is to 20 the more 
%extreme the codon bias is. 
%1) Built-in MATLAB function, aacount, counts how many of each amino acid
%is present in the CDS - ignoring Met (M) and Trp (W) because they are not
%degenerate. Built-in MATLAB function, codoncount, counts how many of each
%codon appears in the CDS. 
%2) Creates a Map object, M, for which codons encode which amino acids - the
%keys are the amino acids (single-letter code) and the values are the
%codon(s) which encode each amino acid 
%3) For each amino acid (excluding M and W), codon homozygozity is
%calculated using the formula from Wright's paper 
%4) Average homozygosity is then computed for each 'synonymous family' (SF)
%(synonymous families are groups of amino acids which have the same number of
%synonymous codons). Average homozygosity is computed by computing the mean 
%of codon homozgyozitiies for amino acids belonging to the same SF. In this
%function, average homozygosities are denoted by F2 (for amino acids which
%have two synonymous codons), F3 (for amino acids which have three
%synonymous codons), F4 (for amino acids which have four synonymous codons)
%and F6 (for amino acids which have six synonymous codons) 
%5) Finally, the Nc value of a gene is computed by the following equation: 
% Nc = 2 + 9/F2 + 1/F3 +5/F4 + 3/F6
%Notes: In the original paper by Wright, two adjustments to the method of
%computing Nc were mentioned. i)Rare or missing amino acids: If the codon 
%homozygosity of an amino acid is equal to 0 (the case for rare missing amino acids)
%the amino acid should be treated as absent and the average homozygosity
%should be computed only for present amino acids. However, the SF F3 only
%has one amino acid in it (Isoleucine is the only amino acid to be encoded
%by three synonymous codons) and if Isoleucine is absent, F3 should be
%comptued as the average of F2 and F4. ii) It is possible that the value of
%Nc can be greater than 61 (this occurs if codon usage is more uniform than
%expected by chance). In this scenario, Wright suggests capping the value
%of Nc at 61. Both of these adjustments were incorporated into this fucntion

%% Counts amino acids in sequence 
AAStruct = aacount(nt2aa(Seq));

%% Count codons in sequence 
CC = codoncount(Seq);

%% Remove fields for Methionine and Tryptophan - because these amino acids are not degenerate
AAStruct = rmfield(AAStruct,{'M','W'});

%% Creates a Map object for which codons encode which amino acids 
Map = rmfield(revgeneticcode(1),{'Name','Starts','Stops','M','W'});
keySet = fieldnames(Map); %amino acids are the fieldnames of Map defined above
valueSet = struct2cell(Map); %a cell array of codons that encode each amino acid
M = containers.Map(keySet, valueSet, 'UniformValues', false); %creates a Map object
%%
AminoAcidNames = fieldnames(Map);

%% Compute Codon homozygosities 
F_values = zeros(18,1);
for i=1:18
    A = 0;
    key = AminoAcidNames{i,1};
    syn_codons = values(M,cellstr(key));
    NumberOfCodons = width(syn_codons{1,1});
    for j=1:NumberOfCodons
        codonfreq = (CC.(syn_codons{1,1}{1,j}))/(AAStruct.(AminoAcidNames{i,1}));
        A = A + (codonfreq)^2;
    end 
    n = AAStruct.(AminoAcidNames{i,1});
    F_values(i,1) = ((n*A) - 1)/(n-1);
end 

F_values(isnan(F_values))=0;
F_values = [AminoAcidNames num2cell(F_values) valueSet];


%% Compute average homozygosities for each SF 
F2 = 0;
F3 = 0;
F4 = 0;
F6 = 0;

denom2 = 0;
denom3 = 0;
denom4 = 0;
denom6 = 0;

for i=1:18 
    if width(F_values{i,3}) == 2
        F2 = F2 + F_values{i,2};
        if F_values{i,2} ~= 0
           denom2 = denom2+1;
        end 
    end  
    
    if width(F_values{i,3}) == 3
        F3 = F3 + F_values{i,2};
        if F_values{i,2} ~= 0
        denom3 = denom3+1;
        end 
    end    
    
    if width(F_values{i,3}) == 4
        F4 = F4 + F_values{i,2};
        if F_values{i,2} ~= 0
        denom4 = denom4+1;
        end 
    end 
    
    if width(F_values{i,3}) == 6
        F6 = F6 + F_values{i,2};
        if F_values{i,2} ~= 0
        denom6 = denom6+1;
        end 
    end        
    
end 

F2 = F2/denom2;
F3 = F3/denom3;
F4 = F4/denom4;
F6 = F6/denom6;

%If F3 is equal to 0 (i.e. if the codon homozyogisty of Isoleucine is equal
%to zero), F3 should be computed as the average of F2 and F4.
if isnan(F3)
    F3 = (F2 + F4)/2;
end 

%% Effective number of codons in a gene 
Nc = 2 + 9/F2 + 1/F3 +5/F4 + 3/F6;

if Nc > 61
    Nc = 61;
end 
end

