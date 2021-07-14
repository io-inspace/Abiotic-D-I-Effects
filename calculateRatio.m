function [outputratio] = calculateRatio(avg_ratio,scale,dim,m)


%this goes through each sampling scale and calculates species richness at
%the end of the simulation.

            numgrid = dim/scale(m);
            record_rat = zeros(numgrid,numgrid);
            indexv = 0;
            for v = 1:numgrid
                indexw = 0;
                for w = 1:numgrid
                    block = avg_ratio((indexv+1):((indexv)+scale(m)),(indexw+1):(w*scale(m)));
                    record_rat(v,w) = mean(mean(block,'omitnan'));
                    indexw = (w*scale(m)); 
                end
                indexv = (v*scale(m));    
            end    

outputratio = record_rat;
end
