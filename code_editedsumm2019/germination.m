function [outputgermind] = germination(dim,maxsp,seedind,G,npG,beta,sigsq,gen,diffG,savg_germs)

totalseed = sum(seedind,3);
occupiedseeds = find(totalseed>0);%gives you the indices to these microsites that have seeds

%save a version of newmatrix, you'll need it to do the beginning of the o
%loop.  Call it...oldind so that you can do all the other old ind stuff
%with it???  No call it something different for now you can always rename
%it.

germind = zeros(dim,dim,maxsp);
rng(gen*0.0034);
for o = 1:length(occupiedseeds)    
    getdim = size(seedind);
    [x,y]=ind2sub(getdim(1:2),occupiedseeds(o));
    
        whichsps = find(seedind(x,y,:));
        %next, find the probability of germinating.
        Rbysps = zeros(1,length(whichsps));
        for bysps = 1:length(whichsps);
            if diffG == 1
            Rbysps(1,bysps) = beta*exp((-.5)*(G(x,y)-(npG(whichsps(bysps))))^2/(sigsq));
            else
                Rbysps(1,bysps) = savg_germs(1,whichsps(bysps));
            end
        end
            %ok new deal -- get the germination probability.  If they
            %germinate then they compete, and that's just the standard code
            %getting us down to k.  Yay.
        for bysps = 1:length(whichsps)                
            germfate = rand(1,seedind(x,y,whichsps(bysps)));
            for germ = 1:length(germfate)
                if germfate(germ) < Rbysps(bysps) %if this is true, add an individual to germind
                    germind(x,y,whichsps(bysps))=germind(x,y,whichsps(bysps))+1;
                end
            end
        end
end
outputgermind = germind;
end
