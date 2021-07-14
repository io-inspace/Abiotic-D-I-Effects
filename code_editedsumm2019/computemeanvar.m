%edited 10-18

%returns the mean and variance of the abiotic landscape at each scale.

function [mean,variance] = computemeanvar(dim,ls)
%across different scales (just do two to start, small and large)
%if you do the grid method, do this:
divisibles = zeros(1,dim);
for K=1:dim;
divisibles(K) = rem(dim,K);
end
scale = find(divisibles == 0);
scale = scale(2:length(scale)-1);


mean = cell(1,length(scale)); 
variance = cell(1,length(scale));
for s = 1:length(scale)
    numgrid = dim/scale(s);
    recordmean = zeros(numgrid,numgrid);
    recordvar = zeros(numgrid,numgrid);
    indexi = 0;
    for i = 1:numgrid
        indexj = 0;
        for j = 1:numgrid
              block = ls((indexi+1):((indexi)+scale(s)),(indexj+1):(j*scale(s)));  
              recordmean(i,j) = mean2(block);
              recordvar(i,j) = var(reshape(block,1,length(block)^2));
             
              indexj = (j*scale(s));
        end
        indexi = (i*scale(s));
    end
    mean{s} = recordmean;
    variance{s} = recordvar;
end