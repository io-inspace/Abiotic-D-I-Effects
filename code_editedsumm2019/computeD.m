

function [allfractalsA,allfractalsB] = computeD(dim,A,B,scale)

allfractalsA = cell(1,length(scale));
allfractalsB = cell(1,length(scale));

for z = 1:length(scale)
            numgrid = dim/scale(z);
            recordfractalsA = zeros(numgrid,numgrid);
            recordfractalsB = zeros(numgrid,numgrid);
            indexv = 0;
            for v = 1:numgrid
                indexw = 0;
                for w = 1:numgrid
                    [outfractalsA] = diagnostics(scale(z),A((indexv+1):((indexv)+scale(z)), (indexw+1):(w*scale(z))));
                    fractalsA = 3-(outfractalsA(1)/2);
                    recordfractalsA(v,w) = fractalsA;
                    [outfractalsB] = diagnostics(scale(z),B((indexv+1):((indexv)+scale(z)), (indexw+1):(w*scale(z))));
                    fractalsB = 3-(outfractalsB(1)/2);
                    recordfractalsB(v,w) = fractalsB;
                    indexw = (w*scale(z));
                end
                indexv = (v*scale(z));
            end
            allfractalsA{z} = recordfractalsA;
            allfractalsB{z} = recordfractalsB;
end

