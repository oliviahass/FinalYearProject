function [] = Nc_plot(AllChromosomes)
%%
Nc_vals = zeros(1,height(AllChromosomes));
GC3_vals= zeros(1,height(AllChromosomes));

for i=1:height(AllChromosomes)
    Nc_vals(1,i) = Nc(string(AllChromosomes.Sequence(i)));
    GC3_vals(1,i) = GC3(string(AllChromosomes.Sequence(i)));
end 

figure;
scatter(GC3_vals,Nc_vals,60,'.','black')
xlim([0 1.0])
hold on
fplot(@(x) 2+x+29*(x^2 +(1-x)^2)^(-1), [0 1],'b')
xlabel('GC3_{s}')
ylabel('Nc')



end

