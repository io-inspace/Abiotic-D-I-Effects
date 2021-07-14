function [outputrichness] = calculateS(richness,scale,spsmatrix,dim,maxsp,ic)


%this goes through each sampling scale and calculates species richness at
%the end of the simulation.

for eachsp = 1:maxsp
        for z = 1:length(scale)
            numgrid = dim/scale(z);
            recordS = zeros(numgrid,numgrid);
            indexv = 0;
            for v = 1:numgrid
                indexw = 0;
                for w = 1:numgrid
                    block = spsmatrix((indexv+1):((indexv)+scale(z)),(indexw+1):(w*scale(z)),eachsp);
                    recordS(v,w) = sum(sum(block));
                    indexw = (w*scale(z)); 
                end
                indexv = (v*scale(z));    
            end
            richness{eachsp,z,ic} = recordS;
        end
end
outputrichness = richness;
end
