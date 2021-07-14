%returns the mean and variance of the abiotic landscape at each scale.

function [mean,variance,cv] = computemeanvarcv(dim,ls,scale)
%across different scales (just do two to start, small and large)
%if you do the grid method, do this:

%you need to modify the landscape for cv:
modified_ls = ls+abs(min(min(ls)));

mean = cell(1,length(scale)); 
variance = cell(1,length(scale));
cv = cell(1,length(scale));
for s = 1:length(scale)
    numgrid = dim/scale(s);
    recordmean = zeros(numgrid,numgrid);
    recordvar = zeros(numgrid,numgrid);
    recordcv = zeros(numgrid,numgrid);
    indexi = 0;
    for i = 1:numgrid
        indexj = 0;
        for j = 1:numgrid
              block = ls((indexi+1):((indexi)+scale(s)),(indexj+1):(j*scale(s)));  
              recordmean(i,j) = mean2(block);
              recordvar(i,j) = var(reshape(block,1,length(block)^2));
              cvblock = modified_ls((indexi+1):((indexi)+scale(s)),(indexj+1):(j*scale(s)));
              recordcv(i,j) = std2(cvblock)/mean2(cvblock);
              indexj = (j*scale(s));
        end
        indexi = (i*scale(s));
    end
    mean{s} = recordmean;
    variance{s} = recordvar;
    cv{s} = recordcv;
end