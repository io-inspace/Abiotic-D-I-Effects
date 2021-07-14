%EDITED 1-26:
%this now calculates the correct number of seeds created for each sps in a
%given microsite before they're dispersed.

function [seedsout_final] = seedoutput(maxsp,dim,individuals_in,Ma,Mp,R_a,R_p,npR,sigsqR,gen,simR,savg_seeds)
seedsout_final = zeros(dim,dim,maxsp); 
total = sum(individuals_in,3);
%go through each occupied microsite:
occupied = find(total>0); %because you have your check calculating the new total from above.
rng(gen*0.8980);
for o = 1:length(occupied); %do this microsite by microsite
    %get the x,y coordinate for the occupied(o) you're on.
    getdim = size(individuals_in);
    [x,y]=ind2sub(getdim(1:2),occupied(o));
    
    whichsps = find(individuals_in(x,y,:)); %which sps are reproducing
            
    for p = 1:length(whichsps); %do for each species in microsite
        %make an if loop here for annuals v. perennials where the only
        %difference is Ma or Mp
        if simR == 1
            if mod(whichsps(p),2) == 1 %if the species is odd, it's an annual  
                fecundity = Ma*exp(-.5*(R_a(x,y)-(npR(whichsps(p))))^2/sigsqR);
                %CENSUS ANNUALS
                %allinds(x,y,whichsps(p)) = allinds(x,y,whichsps(p))+spsmatrix(x,y,(whichsps(p)));
            else %otherwise it's a perennial
                fecundity = Mp*exp(-.5*(R_p(x,y)-(npR(whichsps(p))))^2/sigsqR);
                %CENSUS PERENNIALS
                %allinds(x,y,whichsps(p)) = allinds(x,y,whichsps(p))+spsmatrix(x,y,(whichsps(p))); %# individuals present from previous generations, + # in this generation.  For perennials, we don't care whether they're the same or different. 
            end
        else
          fecundity = savg_seeds(1,whichsps(p));
        end
        integ=floor(fecundity);
        fract=fecundity-integ;

        fate = rand(1,individuals_in(x,y,(whichsps(p))));
        numberseeds = zeros(1,length(fate));
        for seed = 1:length(fate)
            if fate(seed)<fract
            numberseeds(1,seed) = integ+1;
            else
                numberseeds(1,seed) = integ;
            end
        end
        seedsout_final(x,y,whichsps(p)) = sum(numberseeds);
    end
end
end

        