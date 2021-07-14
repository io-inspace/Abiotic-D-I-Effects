for b = 1:20
load(['dim_65_',num2str(b),'_simRdiffGdiffC_adj_0_5.mat'])
abundgen = sum(endtotalinds,3)/initcond;
figure(b)
clf
plot(abundgen)
legend('a1','p1','a2','p2')
end

for t = 1:20
load(['simulationv2018_dim65uncorrelatedsets',num2str(t),'neutRdiffGdiffC_adj_5.mat'])
gauscurves(npC,min(min(B)),max(max(B)),sigsq,Ma,B,t)
end

for b = 1:20
load(['simulationv2018_dim65uncorrelatedsets',num2str(b),'neutRdiffGdiffC_adjd_5.mat'])
figure(b)
clf
imagesc(endallinds{1,1}(:,:,3))
colorbar
end

for b = 1:20
load(['simulationv2018_dim65uncorrelatedsets',num2str(b),'neutRdiffGdiffC_adjd_5.mat'])
figure(b)
clf
imagesc(sum(endallinds{1,1},3))
colorbar
end