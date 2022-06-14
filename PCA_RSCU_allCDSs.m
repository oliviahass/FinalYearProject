function [score] = PCA_RSCU_allCDSs(AllChromosomes,CAI_Scores)

%PCA_RSCU_allCDSs 
CDSs = AllChromosomes.Sequence;
RSCUvalues_forallgenes = zeros(59,height(AllChromosomes));

%%
for i=1:height(AllChromosomes)
    RSCUvalues_t = RSCU(string(CDSs(i,1)));
    RSCUvalues_forallgenes(:,i) = RSCUvalues_t.RSCUvalues; 
end

%%
[coeff,score,latent,~,explained,~] = pca(RSCUvalues_forallgenes');


for i=1:height(AllChromosomes)
    GC3_vals(1,i) = GC3(string(CDSs(i,1)));
end 

meanGC3 = mean(GC3_vals);
High = find(GC3_vals > meanGC3); %high GC3 content 
Low = find(GC3_vals < meanGC3); %high GC3 content 
figure;
scatter(score(High,1),score(High,2),50,'.','blue'); %high GC3 content
hold on
scatter(score(Low,1),score(Low,2),70,'.','red');
xlabel(['PC1: ',num2str(explained(1,1)),'%'])
ylabel(['PC2: ',num2str(explained(2,1)),'%'])
title('Prinicpal Component Analysis (PCA) plot of RSCU values')
legend('Higher than average GC3_{s} content','Lower than average GC3_{s} content')

%%
figure;
scatter(coeff(:,1),coeff(:,2),70,'.','red');
xlabel(['PC loading 1 (',num2str(explained(1,1)),'%)'])
ylabel(['PC loading 2 (',num2str(explained(2,1)),'%)'])
title('Prinicpal Component Analysis (PCA) loading plot of RSCU values')
labels1 = {'AAA','AAT','ACA','ACT','AGA','AGT','ATA','ATT','CAA','CAT','CCA','CCT','CGA','CGT','CTA','CTT','GAA','GAT','GCA','GCT','GGA','GGT','GTA','GTT','TAT','TCA','TCT','TGT','TTA','TTT'};
labelpoints(coeff([1,4,5,8,9,12,13,15,16,19,20,23,24,27,28,31,32,35,36,39,40,43,44,47,49,50,53,55,56,59],1),coeff([1,4,5,8,9,12,13,15,16,19,20,23,24,27,28,31,32,35,36,39,40,43,44,47,49,50,53,55,56,59],2),labels1,'FontSize', 4, 'Color', 'r');
hold on 
labels2 = {'AAC','AAG','ACC','ACG','AGC','AGG','ATC','CAC','CAG','CCC','CCG','CGC','CGG','CTC','CTG','GAC','GAG','GCC','GCG','GGC','GGG','GTC','GTG','TAC','TCC','TCG','TGC','TTC','TTG',};
scatter(coeff([2,3,6,7,10,11,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,45,46,48,51,52,54,57,58],1),coeff([2,3,6,7,10,11,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,45,46,48,51,52,54,57,58],2),70,'.','blue');
labelpoints(coeff([2,3,6,7,10,11,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,45,46,48,51,52,54,57,58],1),coeff([2,3,6,7,10,11,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,45,46,48,51,52,54,57,58],2),labels2,'FontSize', 4, 'Color','b');

%% Correlation between corrdinates of first PC axis and GC3 content 
figure;
scatter(GC3_vals,score(:,1),70,'.','green');
xlabel('GC3_{s}');
ylabel('PC1');
[R,P] = corrcoef(GC3_vals, score(:,1)');
str = ['Correlation coefficient = ' string(R(1,2))];
dim = [0.6 0.6 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Correlation between position on PC1 and GC3_{s}')

%% Correlation between corrdinates of second PC axis and GC3 content 
figure;
scatter(GC3_vals,score(:,2),70,'.','green');
xlabel('GC3_{s}');
ylabel('PC2');
[R,p] = corrcoef(GC3_vals, score(:,2));
str = ['r = ' string(R(1,2))];
dim = [0.6 0.6 0.4 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Correlation between position on PC2 and GC3_{s}')

%% Correlation between coordinates on first PC axis and Nc 
Nc_vals = zeros(1,height(AllChromosomes));

for i=1:height(AllChromosomes)
    Nc_vals(1,i) = Nc(string(AllChromosomes.Sequence(i)));
end 

figure;
scatter(score(:,1),Nc_vals,70,'.','black');
xlabel(['PC1: ',num2str(explained(1,1)),'%'])
ylabel('Nc value');
title('Effective number of codons against the first PC axis')


%% CAI vs PC1 - CAI scores obtained from codonW 
CAI = readtable(CAI_Scores);
figure;
scatter(score(:,1),table2array(CAI(:,2)),70,'.','red');
xlabel(['PC1: ',num2str(explained(1,1)),'%']);
ylabel('CAI scores');
[R,p] = corrcoef(score(:,1), table2array(CAI(:,2)));
str = ['r = ' string(R(1,2))];
dim = [0.6 0.6 0.4 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Correlation between position on PC1 and CAI scores')

%% CAI vs PC2 - CAI scores obtained from codonW 
figure;
scatter(score(:,2),table2array(CAI(:,2)),70,'.','red');
xlabel(['PC2: ',num2str(explained(2,1)),'%']);
ylabel('CAI scores');
[R,p] = corrcoef(score(:,2), table2array(CAI(:,2)));
str = ['r = ' string(R(1,2))];
dim = [0.6 0.6 0.4 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Correlation between position on PC2 and CAI scores')

end

