%fitness output
function [fitness] = calculatefit(scale,fitness,seedfrac,dim,maxsp,ic)

seeds = seedfrac;
seeds(isnan(seeds))=0;
for z = 1:length(scale)
    numgrid = dim/scale(z);
    recordfit = zeros(numgrid,numgrid,maxsp);
    for eachsp = 1:maxsp
            indexv = 0;
            for v = 1:numgrid
                indexw = 0;
                for w = 1:numgrid
                    block = seeds((indexv+1):((indexv)+scale(z)),(indexw+1):(w*scale(z)),eachsp);
                    %recordfit(v,w) = std2(block)/mean2(block);
                    %index=block==Inf;
                    %if find(index==1)>0
                    %    recordfit(v,w,eachsp) = var(block(index==0));
                    %else
                        recordfit(v,w,eachsp) = var(reshape(block,1,length(block)^2));
                    %end
                    %if isnan(recordfit(v,w))==1
                    %    recordfit(v,w) = 0;
                    %end
                    indexw = (w*scale(z)); 
                end
                indexv = (v*scale(z));    
            end
    end
    avgfit = sum(recordfit,3)./4;
    fitness{ic,z} = avgfit;
end
end
