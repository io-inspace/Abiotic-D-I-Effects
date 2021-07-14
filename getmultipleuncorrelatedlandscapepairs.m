%edited 7-19

%specify the number of landscape pairs
t = 20; 

saveit=0; %don't overwrite the version used in the analysis, which is already saved as
          %setsofuncorrelatedlandscapes65.mat

%specify landscape dimension
dim = 65;
dimension = dim;

%specify how many landscapes to make at a time (i.e. could be multiple if
%we don't care about correlation structure or Dout).
iterations = 1;

%specify whether to calculate Dout or not 
dodiagnostics = 1;

%specify input fractal dimension, Din
fractaldims = 2.5; %could also be a vector if you want to look across landscapes 
                   %with different D

%specify how close Dout must be to Din
Dcutoff = .1; %how close do you want output landscape D to be to fractaldims?

%specify whether to normalize landscapes
norm = 1;

%specify what standard deviation normalized landscapes should have
STD = 2.5;

%specify whether to do successive additions, sa (sa = 1 adds wiggle to all values
%after every rotation during midpoint displacement)
sa = 1; 

%store output
all_landscapes = cell(t,1); 
landscapes = cell(iterations,length(fractaldims));
fractalD = zeros(iterations,length(fractaldims));

for rep = 1:t

stopper1 = 0;
stopper2 = 0;

[outA, outB, outcorrelationmat,outfractalD_B,outfractalD_A] = uncorrelatedlandscapes(stopper1, stopper2, dodiagnostics,Dcutoff,dimension,iterations,fractaldims,fractalD,landscapes,STD,sa,norm,dim);
A = outA;
B = outB;
correlationmat = outcorrelationmat;

fractalD_A = outfractalD_A;
fractalD_B = outfractalD_B;
all_landscapes{rep,1} = {A, B, correlationmat,fractalD_B,fractalD_A};

end

if saveit == 1
filename = ['setsofuncorrelatedlandscapes' num2str(dimension) '.mat'];
save(filename);
end


