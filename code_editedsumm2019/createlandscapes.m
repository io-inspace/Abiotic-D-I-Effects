%edited 7-19

%a function for getting landscapes and running diagnostics on them:
%via diagnostics.m

function [outlandscapes,outfractalD] = createlandscapes(dimension,iterations,fractaldims,fractalD,landscapes,STD,sa,norm,dodiagnostics,Dcutoff)

again = 1; %make landscapes until Dout is within Dcutoff of Din.

while again == 1
dim = dimension;
k = dim - 1; %k is the number of diagonals in the matrix

%how many steps are needed to fill it all in?
maxlevel = log(dimension-1)/log(2); 

%landscape border:
border = (2^maxlevel); %borders will allow the midpoint averages to calculate
inset = (2^maxlevel)/2; %insetting builds the matrix within the borders, these will be removed at the end
rows = dim + border; %matrix + border = inset on all four sides
cols = dim + border; %matrix + border = inset on all four sides

for z = 1:iterations
for d = 1:length(fractaldims)

%PARAMETERS
D = fractaldims(d); %set the fractal dimension D = 2:3.  D = 2 should be a perfect gradient
H = 3 - D; %H, the scaling factor, H = 0:1.
%sigmasq = 1; %this is always 1, and I can just adjust the SD of the landscape
           %later via Palmer.  That way spatial autocorrelation is not dependent
           %on sd in some obscure way.
           
%the landscape:
landscape = zeros(rows,cols);

%CORNERS:
%corners (1,1) and (last,last) are the intermediate values, while corners
%(1,last) and (last,1) are the extreme values: [1 3 2 2]

%in this version, sigmasq is 1 so just put corners down then define the initial delta.
%corners = sort(sigmasq*randn(1,3));%this is to get close to a pretty gradient look
corners = sort([-1.225,1.225,.01]);
landscape(rows-inset,1+inset) = corners(1); %corner (last,1)
landscape(1+inset,cols-inset) = corners(3); %corner (1,last)
landscape(1+inset,1+inset) = corners(2); %corner (1,1)
landscape(rows-inset,cols-inset) = corners(2); %corner (last,last)

delta = var([-1.225,1.225,.01,.01]);

%LANDSCAPE LOOP:


for n = 1:maxlevel; %repeat for each level.
    f = 2^(maxlevel-n); %define the distance to (not points between) neighbors, which is a function of level.
    
    %DIAGONALS
    %first, do the rotation delta:
    if D == 2;
        delta=0;
    else
    delta = delta*(.5)^(H/2); 
    end
    
    [outlandscape] = diagonals(landscape,delta,f,k,sa);
    landscape = outlandscape;

%OFF DIAGONALS
   
   %define a new delta to account for this rotation
    delta = delta*(.5)^(H/2);
    
    [outlandscape] = offdiagonals(landscape,f,dim,inset,delta,sa);
    landscape = outlandscape;
   
%index = n; %do this for all the levels until the landscape is finished

end

%trim off the borders:
landscape = landscape(1+inset:rows-inset,1+inset:cols-inset);

if norm == 1;
%modify sensu Palmer:

A = reshape(landscape,numel(landscape),1);
A = zscore(A); %standardize for mean = 0 and unit SD
A = STD.*A; %adjust for SD, NOT SIGMA

%now you have to put this back into matrix form
normedsd = reshape(A,[dim,dim]);
else
    normedsd = landscape;
end

landscapes{z,d} = normedsd;

%THIS IS THE CHUNK OF CODE THAT MAKES THE SEMIVARIOGRAM AND GETS THE SLOPE
%(i.e. Dout)

if dodiagnostics == 1;

    [outputslope] = diagnostics(dim,normedsd);
    fractalD(z,d) = 3-(outputslope(1)/2);
    if abs(fractaldims - fractalD)>Dcutoff
        again = 1;
    else
        again = 0;
    end
    
end

end
end

outlandscapes = landscapes;
outfractalD = fractalD;

end

