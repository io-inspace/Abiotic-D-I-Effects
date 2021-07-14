%EDITED 1-27
function [neutspatialavg_germs,neutspatialavg_compeff] = neutBbaselines(dim,initcond,numls,maxsp,dist_a,dist_p,prop,hab,beta,alpha,dmech,ddstrength,models,rootdir)

if hab == 1;
propBsuitable = zeros(numls,maxsp);
%propRsuitable = zeros(numls,maxsp); R is held constant! So there is no
%need to account for it.
end

rando = 0.1563;
rng(rando)

neutspatialavg_germs = zeros(numls,maxsp);
neutspatialavg_compeff = zeros(numls,maxsp);
for t = 1:numls
%le sigh. you didn't save seedind, so just take this final one and extract
%the info you want for year 301.
%for each entry, calculate a seeds-produced matrix
seeds = cell(initcond,maxsp);
germs = cell(initcond,maxsp);
seedlings = cell(initcond,maxsp);
compeff = cell(initcond,maxsp);


%get the CORRECT full simulation
filepathString = char(['/Users/s_catella/Documents/MATLAB/simulation2019Jan/neutral/' dmech '/' ddstrength '/' char(models(1))]);
cd(filepathString)
load(['dim_' num2str(dim) '_' num2str(t) '_' char(models(1)) '_' dmech, ddstrength, '.mat'])
cd(rootdir);

if hab == 1
     propBsuitable(t,1) = length(find(B<=breakpoints(t,1)))/(dim^2);
     propBsuitable(t,2) = length(find(B>breakpoints(t,1) & B<=breakpoints(t,2)))/(dim^2);
     propBsuitable(t,3) = length(find(B>breakpoints(t,2) & B<=breakpoints(t,3)))/(dim^2);
     propBsuitable(t,4) = length(find(B>breakpoints(t,3)))/(dim^2);
     %propRsuitable(t,1) = length(find(A>=seedbreak_a(1,1) & A<=seedbreak_a(1,2)))/(dim^2);
     %propRsuitable(t,2) = length(find(A>=seedbreak_p(1,1) & A<=seedbreak_p(1,2)))/(dim^2);
     %propRsuitable(t,3) = length(find(A>=seedbreak_a(1,1) & A<=seedbreak_a(1,2)))/(dim^2);
     %propRsuitable(t,4) = length(find(A>=seedbreak_p(1,1) & A<=seedbreak_p(1,2)))/(dim^2);
end



for j = 1:initcond
for n = 1:maxsp

     species = endallinds{1,j}(:,:,n);
     seedind = zeros(dim,dim);
     germind = zeros(dim,dim);
     seedsprod = zeros(dim,dim);
     %seedindbeforedispersal = zeros(dim,dim);
     occupied = find(species>0);%gives you the indices to these microsites that have seeds
        for o = 1:length(occupied)    
                getdim = size(seedind);
                [x,y]=ind2sub(getdim(1:2),occupied(o));
                
            fecundity = neutspatialavg_seeds(t,n);
            %if mod(n,2) == 1 %if the species is odd, it's an annual
            %fecundity = Ma*exp(-.5*(A(x,y)-(0))^2/sigsq);
            %else
            %fecundity = Mp*exp(-.5*(A(x,y)-(0))^2/sigsq);
            %end
        
        seedsprod(x,y) = fecundity;    
            
        integ=floor(fecundity);
        fract=fecundity-integ;

        fate = rand(1,species(x,y));
        numberseeds = zeros(1,length(fate));
        for seed = 1:length(fate)
            if fate(seed)<fract
            numberseeds(1,seed) = integ+1;
            else
                numberseeds(1,seed) = integ;
            end
        end
        
        %seedindbeforedispersal(x,y) = stay;%save a copy for later
        
        stay = ceil(prop*sum(numberseeds));       
        leave = sum(numberseeds)-stay; %MAKE SURE THESE SUM TO FECUNDITY!!!
        seedind(x,y)=stay; %this gets modified in disperse.

        %DISPERSE
       
         if leave>0 %doesn't specify species
            if mod(n,2) == 1 %if the species is odd, it's an annual
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
                    if dispersex(d) <= 65 && dispersex(d) > 0 && dispersey(d) <= 65 && dispersey(d) > 0
                    seedind(dispersex(d),dispersey(d))=seedind(dispersex(d),dispersey(d))+1;
                    end
                end
            
         end  
          
        end
       
     
     germinating = zeros(dim,dim);
            occupiedseeds = find(seedind>0);%gives you the indices to these microsites that have seeds

            for o = 1:length(occupiedseeds)    
                getdim = size(seedind);
                [x,y]=ind2sub(getdim(1:2),occupiedseeds(o));
                germinating(x,y) = beta*exp((-.5)*(B(x,y)-(npG(n)))^2/(sigsq));
                germfate = rand(1,seedind(x,y));
                germind(x,y) = sum(germfate < germinating(x,y)); %how many germinate?
            end
    
    if hab == 1
       %propoccupied = length(find(seedsprod>0))/(dim^2);
       seeds{j,n} = mean2(seedsprod)/propBsuitable(t,n);
    else
    seeds{j,n} = mean2(seedsprod);
    end
    
    germs{j,n} = mean2(germinating);
    
    
    seedlings{j,n} = germind;

end

lings = zeros(dim,dim,4);
sps = zeros(dim,dim,4);
comps = zeros(dim,dim,4);
    for x = 1:dim
        for y = 1:dim
            for n = 1:4
                lings(x,y,n) = seedlings{j,n}(x,y);
                sps(x,y,n) = endallinds{1,j}(x,y,n);
            end
            if sum(lings(x,y,:)) > 0 
    
            %whichsps = find(lings(x,y,:));
                 
       
            allnumindsps = lings(x,y,:)+sps(x,y,:); %sum(germind(x,y,inds))+sum(spsmatrix(x,y,inds));
                   
            allwhichsps = find(allnumindsps);
                
            Rbysps = zeros(1,length(allwhichsps)); %need the competitive effect of all species,
                                                           %including adults
                for bysps = 1:length(allwhichsps);
                   Rbysps(1,bysps) = alpha*exp((-.5)*(B(x,y)-(npC(allwhichsps(bysps))))^2/(sigsq));
                   comps(x,y,allwhichsps(bysps)) = Rbysps(bysps);
                end
            end
        end
    end
                for n = 1:4
                compeff{j,n} = mean2(comps(:,:,n));
                end
               

end
for n = 1:4          
%fullspatialavg_seeds(t,n) = mean(cat(1,seeds{:,n}));
neutspatialavg_germs(t,n) = mean(cat(1,germs{:,n}));
neutspatialavg_compeff(t,n) = mean(cat(1,compeff{:,n}));
end
end

