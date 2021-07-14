%EDITED 11-13, as in competition_sym.m, thinning now occurs according to
%initial densities, not progressively updated densities.

%competition within the loop is asymmetric:
function [outputcompmatrix] = competition_asym(germind,spsmatrix,dim,k,alpha,sigsq,C,npC,gen,diffC,savg_comps)
totalgerm = sum(germind,3);
occugerm = find(totalgerm>0); 
rng(gen*0.9058);
for o = 1:length(occugerm)%do this microsite by microsite
    getdim = [dim dim];
    [x,y]=ind2sub(getdim,occugerm(o)); %get the x,y coords that have germinated seedlings
    openspace = k-sum(spsmatrix(x,y,:)); %calculate how many new plants can fit (0-k)
    if openspace > 0 %if there's any space
        %calculate the number of germinated seedlings and which species
        %they are that could potentially go in that microsite
        whichsps = find(germind(x,y,:));
        if sum(germind(x,y,:)) <= openspace %and if there are more open spaces than 
                                                             %germinated seedlings
            %then all seedlings that germinate get to fill up the rest of the microsite:                                                 
            for place = 1:length(whichsps) %...of each species
                spsmatrix(x,y,whichsps(place))=spsmatrix(x,y,whichsps(place))+germind(x,y,whichsps(place));
            end
        else%if there are more germinated seedlings than open spaces
           if length(whichsps) == 1  %and if they're all the same species, 
                                     %no need for competition, just fill
                                     %the empty spaces.
               spsmatrix(x,y,whichsps)=spsmatrix(x,y,whichsps)+openspace;
            
           else %if they're NOT all the same species, a competition loop will weed out the extras.
                %this is tricky because this will be asymmetrical competition:
                %adults have an effect on seedlings but seedlings don't
                %have an effect on adults.
                remove = sum(germind(x,y,:))-openspace;
                germind_entry = germind(x,y,:);
                spsmatrix_entry = spsmatrix(x,y,:);
                if diffC == 1
                C_entry = C(x,y);
                else
                    C_entry = savg_comps;
                end
                [outnewentry] = bigloop(remove,germind_entry,spsmatrix_entry,whichsps,alpha,C_entry,npC,sigsq,openspace,k,diffC);
                spsmatrix(x,y,:) = outnewentry;
           end
        end
    end
end
outputcompmatrix=spsmatrix;
end
