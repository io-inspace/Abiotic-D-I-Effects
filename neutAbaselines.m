function[neutspatialavg_seeds] = neutAbaselines(hab,initcond,maxsp,Ma,Mp,sigsq,dim,numls,breakpoints)

load(['setsofuncorrelatedlandscapes' num2str(dim)],'all_landscapes')

if hab == 1;
propBsuitable = zeros(numls,maxsp);
end

neutspatialavg_seeds = zeros(numls,maxsp); %AVERAGED OVER INITIAL CONDITIONS FOR A LANDSCAPE


for t = 1:numls %1:length(all_landscapes);

landsc = all_landscapes{t,1};

A = landsc{1,1};
B = landsc{1,2};

if hab == 1
     propBsuitable(t,1) = length(find(B<=breakpoints(t,1)))/(dim^2);
     propBsuitable(t,2) = length(find(B>breakpoints(t,1) & B<=breakpoints(t,2)))/(dim^2);
     propBsuitable(t,3) = length(find(B>breakpoints(t,2) & B<=breakpoints(t,3)))/(dim^2);
     propBsuitable(t,4) = length(find(B>breakpoints(t,3)))/(dim^2);
end

%define how many individuals total to start with 
Ntot = dim^2; 
%split up the total evenly among all species.
NI = repmat(floor(Ntot/maxsp),1,maxsp); 

rng(1);
seeds = cell(initcond,maxsp);


for ic = 1:initcond
    
%STEP 1: RANDOM PLACEMENT
spsmatrix = zeros(dim,dim,maxsp);

for n = 1:length(NI); %do for each species
    for i = 1:NI(n); %do for each individual
        a = 1;
        b = dim+1;%(you have to do +1 with floor)
        r = floor((b-a).*rand(1,2)+a); %draw a random location
        spsmatrix(r(1),r(2),n) = spsmatrix(r(1),r(2),n)+1;
    end
end

for n = 1:length(NI)
     
    seedsprod = zeros(dim,dim);
    occupied = find(spsmatrix(:,:,n)>0);%gives you the indices to these microsites that have seeds
        for o = 1:length(occupied)    
                getdim = size(seedsprod);
                [x,y]=ind2sub(getdim(1:2),occupied(o));
     
            if mod(n,2) == 1 %if the species is odd, it's an annual
            fecundity = Ma*exp(-.5*(A(x,y)-(0))^2/sigsq);
            else
            fecundity = Mp*exp(-.5*(A(x,y)-(0))^2/sigsq);
            end
        
        seedsprod(x,y) = fecundity*spsmatrix(x,y,n);  
        end
        
        if hab == 1
        %propoccupied = length(find(seedsprod>0))/(dim^2);
        if mod(n,2) == 1 %if it's an annual
            if mean2(seedsprod)/propBsuitable(t,n)<=Ma %and has seeds <=Ma, all good.
                seeds{ic,n} = mean2(seedsprod)/propBsuitable(t,n);
            else
                seeds{ic,n} = Ma;
            end
        else %if it's a perennial
            if mean2(seedsprod)/propBsuitable(t,n)<=Mp %and has seeds <=Mp, all good.
                seeds{ic,n} = mean2(seedsprod)/propBsuitable(t,n);
            else
                seeds{ic,n} = Mp;
            end
        end
        else
        seeds{ic,n} = mean2(seedsprod);
        end
end
end

for n = 1:maxsp
    neutspatialavg_seeds(t,n) = mean(cat(1,seeds{:,n}));
end

end
 
end
