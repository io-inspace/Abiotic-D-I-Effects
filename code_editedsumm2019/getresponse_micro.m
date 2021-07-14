%getting microsite-level results
results_micro = zeros(1,2);
for tls = 1:20

load(['dim_65_',num2str(tls),'_neutRneutGdiffC_uni_1_5.mat'])

micrositerichnessoverIC = zeros(dim,dim,length(endallinds));
micrositerichness = zeros(dim,dim);
for ic = 1:length(endallinds)
for x = 1:dim
    for y = 1:dim
abundances = endallinds{1,ic}(x,y,:);
micrositerichness(x,y)=4-length(find(abundances==0));
    end
end
micrositerichnessoverIC(:,:,ic)=micrositerichness;
end

micrositeS=mean(micrositerichnessoverIC,3);

S = reshape(micrositeS,[dim^2,1]);
landscape = repelem(tls,length(S))';
datmat = horzcat(S,landscape);

results_microNEW = vertcat(results_micro,datmat);
results_micro = results_microNEW;
end

resultslist_micro = results_micro(2:84501,:);

results_micro = array2table(resultslist_micro(:,1:2),'VariableNames',{'S','landscape'});
writetable(results_micro,'results_micro.csv','Delimiter',',','QuoteStrings',true)
