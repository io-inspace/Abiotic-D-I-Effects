load('setsofuncorrelatedlandscapes65.mat','all_landscapes')

dim = 65;

saveit = 1;

divisibles = zeros(1,dim);
for K=1:dim;
divisibles(K) = rem(dim,K);
end
scale = find(divisibles == 0);
scale = scale(2:length(scale)-1);

predictors = cell(1,length(scale));

for s = 1:length(scale)  
for ls = 1:length(all_landscapes)
   working = all_landscapes{ls,1};
   A = working{1,1};
   B = working{1,2};
   %mean and variance:
    [mean,variance,cv] = computemeanvarcv(dim,A,scale);
    meanA = mean;
    varA = variance;
    cvA = cv;
    [mean,variance,cv] = computemeanvarcv(dim,B,scale);
    meanB = mean;
    varB = variance;
    cvB = cv;

   %fractal dimension
    [allfractalsA,allfractalsB] = computeD(dim,A,B,scale);
    fractalA = allfractalsA;
    fractalB = allfractalsB; 
  
    meanR = reshape(meanA{1,s},size(meanA{1,s},1)^2,1);
    varR = reshape(varA{1,s},size(varA{1,s},1)^2,1);
    cvR = reshape(cvA{1,s},size(cvA{1,s},1)^2,1);
    fractalR = reshape(fractalA{1,s},size(fractalA{1,s},1)^2,1);
    
    meanC = reshape(meanB{1,s},size(meanB{1,s},1)^2,1);
    varC = reshape(varB{1,s},size(varB{1,s},1)^2,1);
    cvC = reshape(cvB{1,s},size(cvB{1,s},1)^2,1);
    fractalC = reshape(fractalB{1,s},size(fractalB{1,s},1)^2,1);
    
    labels = repelem(ls,size(meanA{1,s},1)^2)';
    
    if ls == 1
    predictors{1,s} = horzcat(labels,meanR,varR,cvR,fractalR,meanC,varC,cvC,fractalC);
    else
        temp = horzcat(labels,meanR,varR,cvR,fractalR,meanC,varC,cvC,fractalC);
        predictors{1,s} = vertcat(predictors{1,s},temp);
    end
 end
 
end

colNames = {'landscape','meanR','varR','cvR','fractalR','meanC','varC','cvC','fractalC'};
smallscale_pred = array2table(predictors{1,1},'VariableNames',colNames);
largescale_pred = array2table(predictors{1,2},'VariableNames',colNames);

if saveit == 1
filename = ['predictors' num2str(dim) '.mat'];
save(filename,'smallscale_pred','largescale_pred');
end

writetable(smallscale_pred,'smallscale_pred.csv','Delimiter',',','QuoteStrings',true)
writetable(largescale_pred,'largescale_pred.csv','Delimiter',',','QuoteStrings',true)

