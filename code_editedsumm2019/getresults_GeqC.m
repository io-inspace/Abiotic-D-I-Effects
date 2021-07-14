%edited 7-2019

%saveit = 0; %I already have versions used for analysis saved, don't overwrite them!
saveit = 1; %FILENAME IS DIFFERENT.

%all necessary functions/files are stored in this directory.
rootdir = cd();

%specify how many replicate landscapes to use.
numls = 20; %20

%specify how many initial conditions (initial random placement of
%individuals) to use.
initcond = 1; 

%specify how many generations to simulate.
generations = 200; 

%choose whether you're doing niche neutral or equivalent reproduction
reproduction = {'neutral'}; %'equivalent'};  

%specify dispersal treatment
dispersal = {'int'}; %{'adj','int','uni'};

%specify per-capita density dependent treatment
interdd = {'_1_0'}; %{'_0_5','_1_0','_1_5'};
a = sqrt([.5 1 1.5]); %in later functions alpha gets multiplied by itself, 
                      %so to get a coefficient = interdd, the square root 
                      %needs to be entered into those functions.
%TREATMENT LOOPS
for repro = 1:length(reproduction)
    for d = 1:length(dispersal)
        for dd = 1:length(interdd)

%specify simulations (N, EF, AMC, or EF + AMC)
if strcmp(reproduction(repro),'neutral') == 1
    %models = {'neutRdiffGdiffC' 'neutRdiffGneutC' 'neutRneutGdiffC' 'neutRneutGneutC'};
    models =  {'neutRdiffGdiffC' 'neutRneutGdiffC' 'neutRneutGneutC'}; %do 
    %these to show that reducing abiotic het doesn't affect neutral but does 
    %effect indirect simulations.
else
    models = {'simRdiffGdiffC' 'simRdiffGneutC' 'simRneutGdiffC' 'simRneutGneutC'};
end

%SIMULATION LOOP
for mods = 1:length(models)

%LANDSCAPE PARAMETERS

%specify the standard deviation of the range of conditions on the landscape.
STD = 2.5; 

%specify landscape dimension.
dim = 65;

%SPECIES PARAMETERS

%max no. seeds produced by annuals (Fa in manuscript)
Ma = 20;

%max no. seeds produced by perennials (Fp in manuscript)
Mp = 5;

%perennial mortality, 1/m*no.inds individuals die every year
m = 6;

%specify the competition coefficient              
alpha = a(dd);

%specify the maximum germination probability
beta = 1;

%specify habitat breadth.
sigsq = 2; 
sigsqR = sigsq; %(specified here because I'm too lazy to go change them in the functions).

%specify how many species should be initialized on the landscape.
maxsp = 4; 

%specify microsite carrying capacity (max established individuals/pixel), c
%in manuscript
k = 3;

%specify what proportion of seeds stay in the microsite
prop = .5; 

%sps identities
perennials = [2,4]; %REMEMBER TO CHANGE WITH MAXSP.

%specify models to simulate.
if strcmp(reproduction(repro),'neutral') == 1 %niche neutral reproduction
    if mods == 1 %EF+AMC
        diffG = 1;
        diffC = 1;                  
        simR = 0;
    elseif mods == 2 %EF
        diffG = 1;
        diffC = 0;                  
        simR = 0;
    elseif mods == 3 %AMC
            diffG = 0;
            diffC = 1;                  
            simR = 0;
    else %NEUTRAL
                diffG = 0;
                diffC = 0;                  
                simR = 0;
    end

else %niche equivalent reproduction
    if mods == 1 %FULL
        diffG = 1;
        diffC = 1;                  
        simR = 1;
    elseif mods == 2 %EF
        diffG = 1;
        diffC = 0;                  
        simR = 1;
    elseif mods == 3 %AMC
            diffG = 0;
            diffC = 1;                  
            simR = 1;
    else %NEUTRAL
                diffG = 0;
                diffC = 0;                  
                simR = 1;
    end
end

%specify dispersal abilities

%standardized to size of VIBI intensives
wsize = 20; %meters

%empirical basis for perennial seed dispersal distance
pdist = 3; %meters

%how much further can annuals disperse compared to perennials?
diffdist = 2; %2x further

if strcmp(dispersal(d),'adj') == 1
    dist_p = 1;
    dist_a = dist_p*diffdist;
elseif strcmp(dispersal(d),'uni') == 1
    dist_p = dim;
    dist_a = dist_p;
elseif strcmp(dispersal(d),'int') == 1
    dist_p = ceil(pdist/(wsize/dim)); 
    dist_a = ceil(dist_p*diffdist);    
end

%load replicate landscapes
load(['setsofuncorrelatedlandscapes' num2str(dim)],'all_landscapes')

%determine what size sampling units to use
divisibles = zeros(1,dim);
for K=1:dim;
divisibles(K) = rem(dim,K);
end
scale = find(divisibles == 0);
scale = scale(2:length(scale)-1);

%loads 20 initial conditions (manuscript only uses 5)
load('originallandscapes65.mat','originallandscapes')

%calculate spatial averages for neutral reproduction, competition, and
%germination
dmech = char(dispersal(d));
ddstrength = char(interdd(dd));

%set the cutoff for bad habitat. this approach keeps the amt of "good habitat" = 
              %for annuals and perennials  
badhab = .25; %e.g. as an annual,if you make >=5 seeds you're in bad habitat. 
              %For perennials it's >=1.25. 
Rsuit_a = Ma*badhab;
Rsuit_p = Mp*badhab; 

%calculate breakpoints (for differentiation = breakpoints.csv. for equivalence = +-2.63)
[breakpoints,reprohab] = getbreakpoints(dim,maxsp,sigsq,Mp,Ma,Rsuit_a,Rsuit_p,numls);

%specify whether to account for the amount of available habitat when
%calculating spatial averages.
hab = 1;

%calculate spatial averages for # seeds produced, germination probability,
%and competitive ability from EF+AMC simulations to be used in N, EF, and AMC
%simulations.
if simR == 0 
    if diffG == 1 && diffC == 1
    neutspatialavg_seeds = neutAbaselines(hab,initcond,maxsp,Ma,Mp,sigsq,dim,numls,breakpoints);
    else
         neutspatialavg_seeds = neutAbaselines(hab,initcond,maxsp,Ma,Mp,sigsq,dim,numls,breakpoints);
        [neutspatialavg_germs,neutspatialavg_compeff] = neutBbaselines(dim,initcond,numls,maxsp,dist_a,dist_p,prop,hab,beta,alpha,dmech,ddstrength,models,rootdir);
    end
else
    if diffG == 1 && diffC == 1
        %there are no baseline values needed in this full-full model
        %EF+AMC+ER (equivalent reproduction simulations were not used in
        %the manuscript).
    else
       [fullspatialavg_germs,fullspatialavg_compeff] = fullbaselines(dim,initcond,numls,maxsp,dist_a,dist_p,prop,hab,beta,alpha,dmech,ddstrength,models,rootdir,breakpoints,reprohab);    
    end
end

%REPLICATE LANDSCAPE LOOP
for t = 1:numls; %1:length(all_landscapes);

%specify landscapes
landsc = all_landscapes{t,1}; %FOR REGULAR SIMS

A = landsc{1,1};
B = landsc{1,2}; 

%specificy baselines to use for that landscape
if simR == 0 
    if diffG == 1 && diffC == 1
    savg_comps = zeros(t,4); %PLACEHOLDER (baseline not used in this simulation)
    savg_germs = zeros(t,4); %PLACEHOLDER (baseline not used in this simulation)
    savg_seeds = neutspatialavg_seeds(t,:);
    else
        savg_comps = neutspatialavg_compeff(t,:);
        savg_germs = neutspatialavg_germs(t,:);
        savg_seeds = neutspatialavg_seeds(t,:);
    end
else
    if diffG == 1 && diffC == 1
        savg_comps = zeros(t,4); %PLACEHOLDER (baseline not used in this simulation)
        savg_germs = zeros(t,4); %PLACEHOLDER (baseline not used in this simulation)
        savg_seeds = zeros(t,4); %PLACEHOLDER (baseline not used in this simulation)
    else
        savg_comps = fullspatialavg_compeff(t,:);
        savg_germs = fullspatialavg_germs(t,:);
        savg_seeds = zeros(t,4); %PLACEHOLDER (baseline not used in this simulation)
    end
end

%I split these up in a previous simulation, and instead of finding all the
%R_as and R_ps I just call them both A here.
R_a = A;
R_p = A;

%specify niche optima for the different types of niche paritioning (e.g. differentiation v.
%equivalence):
npsim = repelem(0,maxsp); 

%Rrange = ((max(max(R)))-(min(min(R))));
%insetR = Rrange/maxsp;
%rangefordiff = (max(max(R))-insetR)-(min(min(R))+insetR); 
%goby = rangefordiff/(maxsp-1);
%npdiffR = (min(min(R))+insetR):goby:(max(max(R))-insetR);

%argh OK this is scaling niches to the lower range of abiotic conditions.
%le sigh.


C = B;  
Crange = ((max(max(C)))-(min(min(C))));
insetC = Crange/maxsp;
rangefordiff = (max(max(C))-insetC)-(min(min(C))+insetC); 
goby = rangefordiff/(maxsp-1);
npdiffC = (min(min(C))+insetC):goby:(max(max(C))-insetC);    

G = B;
Grange = ((max(max(G)))-(min(min(G))));
insetG = Grange/maxsp;
rangefordiff = (max(max(G))-insetG)-(min(min(G))+insetG); 
goby = rangefordiff/(maxsp-1);
npdiffG = (min(min(G))+insetG):goby:(max(max(G))-insetG);    

%so at this point, spatial averages and niches are set based on the high
%het landscape. Now, we're switching to a low het landscape:

%you modified A, so make sure B is the first one (so that it's the one with low STD).
%these are for the N v. I species loss maps:
A = lowerSTD_landscapes{t,2};
B = lowerSTD_landscapes{t,1};

C = B;
G = B;

%how do species partition resources? 


npR = npsim;

npG = npdiffG;

npC = npdiffC;  



%run the simulation
[endallinds,endtotalinds,endrichness,endfitness,endprobdeath] = simulationv2018(dist_a,dist_p,R_a,R_p,C,G,dim,npR,npC,npG,generations,initcond,scale,maxsp,sigsq,sigsqR,Ma,Mp,alpha,beta,m,k,prop,perennials,originallandscapes,diffC,savg_comps,savg_seeds,savg_germs,diffG,simR);

filepathString = char(['/Users/s_catella/Documents/MATLAB/simulation2019Jan/' char(reproduction(repro)) '/' dmech '/' ddstrength '/' char(models(mods))]);
cd(filepathString)
if saveit == 1
        filename = char(['dim_' num2str(dim) '_' num2str(t) '_' char(models(mods)) '_' dmech ddstrength 'lowSTD.mat']);
        save(filename);
end
cd(rootdir);

end

end
        end
    end
end
