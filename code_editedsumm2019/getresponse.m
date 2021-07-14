%EDITED 12-14-2018
%you were computing mean competition felt over species, when really you
%want the average competition felt in a subunit by a species, and then the
%variation over all species. Less variability in competition experienced
%(usually) means increased fitness variance.

files = cd();
cd('/Users/s_catella/Documents/submissions/Ecology2019/uncorrelatedforspatialstats_editedsumm2019');

load('setsofuncorrelatedlandscapes65.mat','all_landscapes','dim')
all_landscapes = cat(1,all_landscapes{:,1});
%determine what size sampling units to use
%touse = [1 2 3 4 5 6 8 9 10 13 15 17 18 19 20];
divisibles = zeros(1,dim);
for K=1:dim;
divisibles(K) = rem(dim,K);
end
scale = find(divisibles == 0);
scale = scale(2:length(scale)-1);

data = cell(length(all_landscapes),length(scale));
for i = 1:length(all_landscapes)

cd(files);
load(['dim_' num2str(dim) '_' num2str(i) '_neutRdiffG#diffC_uni_1_5.mat'])
cd('/Users/s_catella/Documents/MATLAB/simulation2019Jan');

%richness landscape:

[outaverageSoverIC] = averagerichness(endrichness,maxsp,scale,initcond);
averageSoverIC = outaverageSoverIC;

%(spatio) fitness landscape:

%[outaveragefitoverIC] = averagefitness(endfitness,maxsp,scale,initcond,dim);
%averagefitoverICbysps = outaveragefitoverIC;

%get the average fitness across sps:
%averagefitoverIC = cell(length(scale),1);
%for sampunit = 1:length(scale)
%thing = averagefitoverICbysps(:,sampunit);
%entries = (dim/scale(sampunit))^2;
%averagefit = zeros(entries,1);
%for elem = 1:(dim/scale(sampunit))^2
%howmanysp = length(find(cellfun(@(thing) ~isnan(thing(elem)),thing)));
%averagefit(elem,1) = sum(cellfun(@(thing) nansum(thing(elem)),thing))/howmanysp; %this adds across 0s subbed for NaNs, but it doesn't matter bc the previous line gets the correct no. of species to divide by.
%end
%averagefitoverIC{sampunit,1} = averagefit;
%end

for m = 1:length(scale)
    %richness:
        avgfitoverIC=mean(cat(3,endfitness{:,m}),3);
        S = reshape(averageSoverIC{1,m},1,length(averageSoverIC{1,m})^2)';
        fit = reshape(avgfitoverIC,1,length(avgfitoverIC)^2)';
        %average fitness heterogeneity across initial conditions and species        
        %fit_het = averagefitoverIC{m,1};
       %corr_bysps_ratio = zeros((dim/scale(m))^2,4);
       %corr_bysps_inter = zeros((dim/scale(m))^2,4);
       %corr_bysps_intra = zeros((dim/scale(m))^2,4);
       rbysps_avg = zeros((dim/scale(m))^2,4);
       %rbysps_var = zeros((dim/scale(m))^2,4);
       %rbysps_cv = zeros((dim/scale(m))^2,4);
       %comp_totalr = zeros((dim/scale(m))^2,4);
       for n = 1:maxsp
        %sps_ratio=cell(initcond,1);
        %sps_inter=cell(initcond,1);     
        %sps_intra = cell(initcond,1);
        %ratio = cell(initcond,1);
        totalr = cell(initcond,1);
        for j = 1:initcond
            %sps_ratio{j,1} = endcorr_ratio{j,1}{1,m}{1,n};
            %sps_inter{j,1} = endcorr_inter{j,1}{1,m}{1,n};
            %sps_intra{j,1} = endcorr_intra{j,1}{1,m}{1,n};
            %ratio{j,1} = endcomp_inter{j,1}(:,:,n)./endcomp_intra{j,1}(:,:,n);
            %totalr{j,1} = endcomp_inter{j,1}(:,:,n)+endcomp_intra{j,1}(:,:,n);
            totalr{j,1} = endprobdeath{j,1}(:,:,n);
        end
        %spsavg_ratio = nanmean(cat(3,sps_ratio{:}),3);
        %spsavg_inter = nanmean(cat(3,sps_inter{:}),3);
        %spsavg_intra = nanmean(cat(3,sps_intra{:}),3);
        %avg_ratio = nanmean(cat(3,ratio{:}),3);
        avg_totalr = nanmean(cat(3,totalr{:}),3); %average r for species n over ic.
        %then get the scale different average ratios.
        %ratioavg = calculateRatio(avg_ratio,scale,dim,m); 
        [record_rAVG] = calculateTotalr(avg_totalr,scale,dim,m); 

        %corr_bysps_ratio(:,n) = reshape(spsavg_ratio,1,length(spsavg_ratio)^2)';
        %corr_bysps_inter(:,n) = reshape(spsavg_inter,1,length(spsavg_inter)^2)';
        %corr_bysps_intra(:,n) = reshape(spsavg_intra,1,length(spsavg_intra)^2)';
        %comp_bysps_ratio(:,n) = reshape(ratioavg,1,length(ratioavg)^2)';
        rbysps_avg(:,n) = reshape(record_rAVG,1,length(record_rAVG)^2)';
        %rbysps_var(:,n) = reshape(record_rVAR,1,length(record_rVAR)^2)';
        %rbysps_cv(:,n) = reshape(record_rCV,1,length(record_rCV)^2)';

      end
        %meanrvar = mean(rbysps_var,2,'omitnan');
        comphet = zeros((dim/scale(m))^2,1);
        for h = 1:(dim/scale(m))^2
            if sum(~isnan(rbysps_avg(h,:))) <= 1
            comphet(h,1) = NaN;
            else
                comphet(h,1) = nanvar(rbysps_avg(h,:));
            end
        end
        
        %meanrcv = mean(rbysps_cv,2,'omitnan'); %THIS GIVES A VALUE EVEN WHEN THERE'S JUST 1 SPS PRESENT....
        %allics=cat(1,endcorr{:});
        %avgcorroverIC = nanmean(cat(3,allics{:,m}),3);
        %comp_Coptcorr = reshape(avgcorroverIC,1,length(avgcorroverIC)^2)';
        landscapepair = repelem(i,length(S))';  
        %data{i,m} = [S,fit_het,landscapepair];
        %data{i,m} = [S,fit,meantotalr,comp_totalr,meancomp,comp_bysps_ratio,corr_bysps_intra,corr_bysps_inter,corr_bysps_ratio,landscapepair];
        data{i,m} = [S,fit,comphet,landscapepair];
        %I THINK ALL YOU WANT IS COMP_TOTALR!!!!!
end
end

cd(files)
%load('predictors65.mat')
resultslist_small = vertcat(data{:,1});
results_small = array2table(resultslist_small(:,1:4),'VariableNames',{'S','fit','comphet','landscape'});
writetable(results_small,'results_small.csv','Delimiter',',','QuoteStrings',true)

resultslist_large = vertcat(data{:,2});
results_large = array2table(resultslist_large(:,1:4),'VariableNames',{'S','fit','comphet','landscape'});
writetable(results_large,'results_large.csv','Delimiter',',','QuoteStrings',true)

%results_small = horzcat(smallscale_pred,array2table(resultslist_small(:,1:10),'VariableNames',{'S','fit','R1','R2','R3','R4','R1var','R2var','R3var','R4var'}));
%resultslist_large = vertcat(data{:,2});
%results_large = horzcat(largescale_pred,array2table(resultslist_large(:,1:10),'VariableNames',{'S','fit','R1','R2','R3','R4','R1var','R2var','R3var','R4var'}));

%writetable(results_large,'results_large_neutRdiffGneutC.csv','Delimiter',',','QuoteStrings',true)
