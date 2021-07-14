%edited 10-18

%off diagonals
function [outlandscape] = offdiagonals(landscape,f,dim,inset,delta,sa)

 for x = 1:length(landscape); %go by row instead of by diagonal
    points = find(landscape(x,:));
    for y = 1:length(points)-1; %get the number of midpoints.
    mp = (points(y)+points(y+1))/2; %define where the midpoints are.
    if x == (inset+1) || x == (inset+dim) %these x are edges, so divide by 3
    eavg = (landscape(x-f,mp)+landscape(x,mp+f)+landscape(x+f,mp)+landscape(x,mp-f))/3;
    else %these x are not, so divide by 4
    eavg = (landscape(x-f,mp)+landscape(x,mp+f)+landscape(x+f,mp)+landscape(x,mp-f))/4;
    end
    
    landscape(x,mp) = eavg + delta*randn;
    end
end

%the previous step gets every entry except for the first and last column
landscape = transpose(landscape); %turn and do the first column (now row)

points = find(landscape(1+inset,:)); 
for y = 1:length(points)-1;
    mp = (points(y)+points(y+1))/2;
    eavg = (landscape(1+inset-f,mp)+landscape(1+inset,mp+f)+landscape(1+inset+f,mp)+landscape(1+inset,mp-f))/3;
    landscape(1+inset,mp) = eavg + delta*randn;
end

landscape = flipud(landscape); %flip up down to make the last column the first row 

points = find(landscape(1+inset,:));
for y = 1:length(points)-1;
    mp = (points(y)+points(y+1))/2;
    eavg = (landscape(1+inset-f,mp)+landscape(1+inset,mp+f)+landscape(1+inset+f,mp)+landscape(1+inset,mp-f))/3;
    landscape(1+inset,mp) = eavg + delta*randn;
end

landscape = flipud(landscape); 
landscape = transpose(landscape); %return to original position

%do random additions
if sa == 1;
for i = 1:length(find(landscape));
    totpnts = find(landscape);
    landscape(totpnts(i)) = landscape(totpnts(i))+delta*randn;
end
end

outlandscape = landscape;

end
