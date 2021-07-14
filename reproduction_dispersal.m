%EDITED 1-26 to allow annuals to disperse further than perennials

%EDITED 1-9 to make fecundity stochastic. Integer number of seeds are taken
%for granted, but the fate of the "less viable" fractional seeds are
%determined by drawing a random number. Thus viability is partly determined
%by the environment, and partly by chance (or, in the case of the neutR
%sims, entirely by chance).

function [outputseedind] = reproduction_dispersal(dist_a,dist_p,dim,maxsp,total,spsmatrix,Ma,Mp,R_a,R_p,npR,prop,sigsqR,gen,simR,savg_seeds)

seedind = zeros(dim,dim,maxsp);  
%go through each occupied microsite:
occupied = find(total>0); %because you have your check calculating the new total from above.
rng(gen*0.1270);
for o = 1:length(occupied); %do this microsite by microsite
    %get the x,y coordinate for the occupied(o) you're on.
    getdim = size(spsmatrix);
    [x,y]=ind2sub(getdim(1:2),occupied(o));
    
    whichsps = find(spsmatrix(x,y,:)); %which sps are reproducing
            
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

        fate = rand(1,spsmatrix(x,y,(whichsps(p))));
        numberseeds = zeros(1,length(fate));
        for seed = 1:length(fate)
            if fate(seed)<fract
            numberseeds(1,seed) = integ+1;
            else
                numberseeds(1,seed) = integ;
            end
        end
        
        stay = ceil(prop*sum(numberseeds));       
        leave = sum(numberseeds)-stay; %MAKE SURE THESE SUM TO FECUNDITY!!!
        seedind(x,y,whichsps(p))=stay; 
       
    
%DISPERSAL
    if leave>0 %doesn't specify species
            if mod(whichsps(p),2) == 1 %if the species is odd, it's an annual
                dist = dist_a;
            else %the species is a perennial
                dist = dist_p;
            end
            
             xsub = x-dist;
             xadd = x+dist;
             ysub = y-dist;
             yadd = y+dist;
             
             %for a non-reflecting boundary
             %get a block of potential dispersal sites within
             %dist of the focal microsite
             dxsub = xsub;
             dxadd = xadd+1;%(you have to do +1 with floor)
             dispersex = floor((dxadd-dxsub).*rand(leave,1)+dxsub); %draw a random location
             dysub = ysub;
             dyadd = yadd+1;
             dispersey = floor((dyadd-dysub).*rand(leave,1)+dysub); %draw a random location
            
             %dispersing seeds aren't allowed to stay:
             for d = 1:leave
                 while dispersex(d) == x && dispersey(d) == y
                   dispersex = floor((dxadd-dxsub).*rand(leave,1)+dxsub); %draw a random location
                   dispersey = floor((dyadd-dysub).*rand(leave,1)+dysub); %draw a random location
                 end
             end
             
                for d = 1:leave
                    if dispersex(d,1) <= 65 && dispersex(d,1) > 0 && dispersey(d,1) <= 65 && dispersey(d,1) > 0
                    seedind(dispersex(d,1),dispersey(d,1),whichsps(p))=seedind(dispersex(d,1),dispersey(d,1),whichsps(p))+1;
                    end
                end
        
            
    end
    end
    
end
    
outputseedind = seedind;
%outputallinds = allinds;
 
    
end

