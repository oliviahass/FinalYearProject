function [GC3_content] = GC3(CDS)
%GC3 is a function that has one input: 
%|CDS| - a CDS or several concateanted CDSs in the form of a string

%GC3 is a function that has one output:
%|GC3_content| - fraction of GC-ending codons excluding ATG, TGG and the
%three stop codons


%GC3 calculates the GC3 content at the third (silent) position of codons
%excluding Methionine (ATG), Tyrptophan (TGG), and the three stop codons
%(TAA, TGA, TAG). It first takes the CDS input (either a single CDS or
%several concatenated CDSs)and uses the built-in function in MATLAB called 
%codoncount to count the number of occurrences of different codons in a 
%sequence. Next, the total number of codons (variable called TotalCodon) 
%is calculated, exclduing the number of Methionine, Tryptophan and stop 
%codons. Finally, the number of GC ending codons is calculated (again 
%excluding any of the GC ending codons mentioned above) and divided by 
%TotalCodon to obtain the GC3 content. 

%% Compute TotalCodon and the codon count of each codon
Codon_occurrences = codoncount(CDS);
TotalCodon = strlength(CDS)/3 - Codon_occurrences.ATG - Codon_occurrences.TGG - Codon_occurrences.TGA - Codon_occurrences.TAG - Codon_occurrences.TAA ;

%% AA_GC3
AA_GC3 = sum(Codon_occurrences.AAC + Codon_occurrences.AAG);

%% AC_GC3
AC_GC3 = sum(Codon_occurrences.ACC + Codon_occurrences.ACG);

%% AG_GC3
AG_GC3 = sum(Codon_occurrences.AGC + Codon_occurrences.AGG);

%% AT_GC3
AT_GC3 = sum(Codon_occurrences.ATC);

%% CA_GC3
CA_GC3 = sum(Codon_occurrences.CAC + Codon_occurrences.CAG);

%% CC_GC3
CC_GC3 = sum(Codon_occurrences.CCC + Codon_occurrences.CCG);

%% CG_GC3
CG_GC3 = sum(Codon_occurrences.CGC + Codon_occurrences.CGG);

%% CT_GC3
CT_GC3 = sum(Codon_occurrences.CTC + Codon_occurrences.CTG);

%% GA_GC3
GA_GC3 = sum(Codon_occurrences.GAC + Codon_occurrences.GAG);

%% GC_CG3
GC_GC3 = sum(Codon_occurrences.GCC + Codon_occurrences.GCG);

%% GG_GC3
GG_GC3 = sum(Codon_occurrences.GGC + Codon_occurrences.GGG);

%% GT_GC3
GT_GC3 = sum(Codon_occurrences.GTC + Codon_occurrences.GTG);

%% TA_GC3
TA_GC3 = sum(Codon_occurrences.TAC);

%% TC_GC3
TC_GC3 = sum(Codon_occurrences.TCC + Codon_occurrences.TCG);

%% TG_GC3
TG_GC3 = sum(Codon_occurrences.TGC);

%% TT_GC3
TT_GC3 = sum(Codon_occurrences.TTC + Codon_occurrences.TTG);


%% GC3
GC3_COUNT = AA_GC3 + AC_GC3 + AG_GC3 + AT_GC3 + CA_GC3 + CC_GC3 + CG_GC3 + CT_GC3 + GA_GC3 + GC_GC3 +GG_GC3 + GT_GC3 + TA_GC3 + TC_GC3 + TG_GC3 + TT_GC3;
GC3_content = GC3_COUNT/TotalCodon;

end

