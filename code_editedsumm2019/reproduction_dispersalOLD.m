
%EDITED 1-9 to make fecundity stochastic. Integer number of seeds are taken
%for granted, but the fate of the "less viable" fractional seeds are
%determined by drawing a random number. Thus viability is partly determined
%by the environment, and partly by chance (or, in the case of the neutR
%sims, entirely by chance).

function [outputseedind] = reproduction_dispersal(dim,maxsp,total,spsmatrix,Ma,Mp,R_a,R_p,npR,prop,sigsqR,gen,adj,uni,dist,simR,savg_seeds)

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
        if adj == 1
    todisperseto = ceil((9-1)*rand(1,leave)); %draws the number of which adjacent cells randomly
    %the potential locations for dispersal need to be modified
    %for corners and edges (so that adding individuals to seedind is
    %skipped if they fall off the edge.
    %corners
    if x == 1 && y == 1; %you're in the upper left, and only cells 3 - 5 are available.
        disperse = todisperseto(3<todisperseto & todisperseto<7);
        elseif x == 1 && y == dim; %you're in the upper right
            disperse = todisperseto(5<todisperseto);
        elseif x == dim && y == 1 %you're in the lower left
            disperse = todisperseto(1<todisperseto & todisperseto<5);
        elseif x == dim && y == dim %you're in the lower right
            disperse = todisperseto(todisperseto == 8 | todisperseto<3);
    
   %edges
        elseif x == 1;
            disperse = todisperseto(3<todisperseto);
        elseif x == dim;
            disperse = todisperseto(todisperseto == 8 | todisperseto<5);
        elseif y == 1;
            disperse = todisperseto(1<todisperseto & todisperseto<7);
        elseif y == dim;
            disperse = todisperseto(5<todisperseto & todisperseto<3);
    else
            disperse = todisperseto;
    end
                
        
    %now that we have modified dispersal corners we can update
    %seedind: 
    for d = 1:length(disperse);
        if disperse(d) == 1;
            seedind(x-1,y-1,whichsps(p))=seedind(x-1,y-1,whichsps(p))+1;
            elseif disperse(d) == 2;
                seedind(x-1,y,whichsps(p))=seedind(x-1,y,whichsps(p))+1;
            elseif disperse(d) == 3;
                seedind(x-1,y+1,whichsps(p))=seedind(x-1,y+1,whichsps(p))+1;
            elseif disperse(d) == 4;
                seedind(x,y+1,whichsps(p))=seedind(x,y+1,whichsps(p))+1;                    
            elseif disperse(d) == 5;
                seedind(x+1,y+1,whichsps(p))=seedind(x+1,y+1,whichsps(p))+1;                    
            elseif disperse(d) == 6;
                seedind(x+1,y,whichsps(p))=seedind(x+1,y,whichsps(p))+1;
            elseif disperse(d) == 7;
                seedind(x+1,y-1,whichsps(p))=seedind(x+1,y-1,whichsps(p))+1;
        else 
                seedind(x,y-1,whichsps(p))=seedind(x,y-1,whichsps(p))+1;
        end
    end   
        elseif uni == 1 %with universal dispersal, you just want to pick
                              %random x,y coordinates for each seed leaving.
            a = 1;
            b = dim+1;%(you have to do +1 with floor)
            disperse = floor((b-a).*rand(leave,2)+a); %draw a random location
                for d = 1:leave
                    seedind(disperse(d,1),disperse(d,2),whichsps(p))=seedind(disperse(d,1),disperse(d,2),whichsps(p))+1;
                end
        else %this leaves int == 1, wherein you specify a restricted subset
             %of x,y coordinates around the microsite and pick random
             %coordinates to disperse to.
             
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
            
                for d = 1:leave
                    if dispersex(d,1) <= 65 && dispersex(d,1) > 0 && dispersey(d,1) <= 65 && dispersey(d,1) > 0
                    seedind(dispersex(d,1),dispersey(d,1),whichsps(p))=seedind(dispersex(d,1),dispersey(d,1),whichsps(p))+1;
                    end
                end
             
             
             %for a reflecting boundary:
             %if xsub > 65 
             %    xsub = 65;
             %elseif xsub <= 0
             %    xsub = 1;
             %end
             
             %if xadd > 65 
             %    xadd = 65;
             %elseif xadd <= 0
             %    xadd = 1;
             %end
             
             %if ysub > 65 
             %    ysub = 65;
             %elseif ysub <= 0
             %    ysub = 1;
             %end
             
             %if yadd > 65 
             %    yadd = 65;
             %elseif yadd <= 0
             %    yadd = 1;
             %end
             
            
        end
    end
    
    end
    
    
end
outputseedind = seedind;
%outputallinds = allinds;
end
