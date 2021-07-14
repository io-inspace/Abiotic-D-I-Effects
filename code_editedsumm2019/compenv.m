%edited 1-26
%the comphet values were all over the place because you weren't using het
%in the normalized probability of dying (which is the competitive pressure
%felt),you were using the raw inter and intra comp effects, which could be
%huge and tiny

function [dyingprob] = compenv(germind,spsmatrix,dim,maxsp,C,npC,alpha,sigsq,diffC,savg_comps,k)
intra_compeffect = NaN(dim,dim,maxsp);
inter_compeffect = NaN(dim,dim,maxsp);
dyingprob = NaN(dim,dim,maxsp);
%competitive_ratio = NaN(dim,dim,maxsp);
totalgerm = sum(germind,3);
occugerm = find(totalgerm>0);   
for o = 1:length(occugerm) %do this microsite by microsite
    getdim = [dim dim];
    [x,y]=ind2sub(getdim,occugerm(o)); %get the x,y coords that have germinated seedlings
    openspace = k-sum(spsmatrix(x,y,:)); %calculate how many new plants can fit (0-k)
        %calculate the number of germinated seedlings and which species
        %they are that could potentially go in that microsite
        
                whichsps = find(germind(x,y,:));
                    if sum(germind(x,y,:))<openspace %if there are fewer germinating inds
                        %do nothin'                  %than open spaces, there's no comp.
                    else %if there are more germinating seedlings than open spaces, there is comp.
                         if length(whichsps) == 1 %if there's only one type of sps germinating, 
                                                  %the probability of them
                                                  %dying isn't
                                                  %interesting, so skip?
                         %do nothin'
                         else %multiple sps types present
                            allnumindsps = zeros(1,maxsp);
                            %allnumindsps = zeros(1,10); %WTF IS THIS 10?!?! >:(
                                            %since adults have a competitive 
                                            %effect, get the TOTAL #
                                            %individuals and what species
                                            %they are in the microsite (so
                                            %old perennials and germinated
                                            %seedlings).
                             for inds = 1:maxsp;
                               allnumindsps(1,inds) = sum(germind(x,y,inds))+sum(spsmatrix(x,y,inds));
                             end

                             allwhichsps = find(allnumindsps);

                             Rbysps = zeros(1,length(allwhichsps)); %need the competitive effect of all species,
                                                                               %including adults
                             for bysps = 1:length(allwhichsps);
                                 if diffC == 1
                                Rbysps(1,bysps) = alpha*exp((-.5)*(C(x,y)-(npC(allwhichsps(bysps))))^2/(sigsq));
                                 else
                                     Rbysps(1,bysps) = savg_comps(1,allwhichsps(bysps));
                                 end
                             end

                             probabilities = zeros(1,length(whichsps));
                             for j = 1:length(whichsps); %calculate probability of dying only for germinating species
                                intra_compeffect(x,y,whichsps(j)) = allnumindsps(whichsps(j))*(allnumindsps(whichsps(j)));
                                compeffect = zeros(1,length(allwhichsps)); %but there's a compeffect of all species
                                                for m = 1:length(allwhichsps);
                                                    if whichsps(j) == allwhichsps(m)
                                                        compeffect(m) = 0;
                                                    else  
                                                        compeffect(m) = Rbysps(m)*alpha*allnumindsps(allwhichsps(m));
                                                    end
                                                end 
                                probabilities(j) = allnumindsps(whichsps(j))*(allnumindsps(whichsps(j))+sum(compeffect));       
                                inter_compeffect(x,y,whichsps(j)) = allnumindsps(whichsps(j))*sum(compeffect);
                                %competitive_ratio(x,y,whichsps(j)) = inter_compeffect(x,y,whichsps(j))/intra_compeffect(x,y,whichsps(j));
                             end
                             normedprob = probabilities/sum(probabilities);
                             for j = 1:length(whichsps); %
                                 dyingprob(x,y,whichsps(j)) = normedprob(j);
                             end
                         end
                    end

end
end
