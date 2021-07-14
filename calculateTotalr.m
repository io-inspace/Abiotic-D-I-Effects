function [record_rAVG] = calculateTotalr(avg_totalr,scale,dim,m)

            numgrid = dim/scale(m);
            record_rAVG = zeros(numgrid,numgrid);
            %record_rVAR = zeros(numgrid,numgrid);
            %record_rCV = zeros(numgrid,numgrid);
            indexv = 0;
            for v = 1:numgrid
                indexw = 0;
                for w = 1:numgrid
                    block = avg_totalr((indexv+1):((indexv)+scale(m)),(indexw+1):(w*scale(m)));
                    record_rAVG(v,w) = mean(block(~isnan(block)));
                    %record_rVAR(v,w) = var(block(~isnan(block)));
                    %record_rCV(v,w) = mean(block(~isnan(block)))/std(block(~isnan(block)));
                    indexw = (w*scale(m)); 
                end
                indexv = (v*scale(m));    
            end    
end