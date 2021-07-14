function [outnewentry] = bigloop(remove,germind_entry,spsmatrix_entry,whichsps,alpha,C_entry,npC,sigsq,openspace,k,diffC)
           
allnumindsps = germind_entry+spsmatrix_entry; %sum(germind(x,y,inds))+sum(spsmatrix(x,y,inds));
                   
allwhichsps = find(allnumindsps);
                
Rbysps = zeros(1,length(allwhichsps)); %need the competitive effect of all species,
                                                           %including adults
    for bysps = 1:length(allwhichsps);
        if diffC == 1
        Rbysps(1,bysps) = alpha*exp((-.5)*(C_entry-(npC(allwhichsps(bysps))))^2/(sigsq));
        else
            Rbysps(1,bysps) = C_entry(1,allwhichsps(bysps));
        end
    end

probabilities = zeros(1,length(whichsps));
     for j = 1:length(whichsps); %calculate probability of dying only for germinating species
        compeffect = zeros(1,length(allwhichsps)); %but there's a compeffect of all species
           for m = 1:length(allwhichsps);
              if whichsps(j) == allwhichsps(m)
                 compeffect(m) = 0;
              else     
                 compeffect(m) = Rbysps(m)*alpha*allnumindsps(allwhichsps(m));
              end
           end
        probabilities(j) = allnumindsps(whichsps(j))*(allnumindsps(whichsps(j))+sum(compeffect));
     end
     normedprob = probabilities/sum(probabilities);
        
        
    [newentry] = whodies(whichsps,normedprob,germind_entry,remove,spsmatrix_entry,openspace,k);
                    
                    outnewentry = newentry;
end