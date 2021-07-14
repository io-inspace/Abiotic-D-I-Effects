%edited 10-18

%trying to get diagnostics to run faster/at all.  modified from
%olddiagnostics.m

%calculate D from the semivariogram
function [outputslope] = diagnostics(dim,normedsd)

entries = (dim^2)-1; %entries minus the first one
output = sum(1:1:entries); %this is the number of values you'll get depending on the size of the landscape

h = zeros(1,output); %initiate the distance vector, remember to take this off.
sqrs = zeros(1,output); %initiate the sum of squares vector, the vector for (z2-z1)^2 -> the difference between their values.
maxh = sqrt(((1-dim)^2)+((1-dim)^2));

index = 0;
for e = 1:entries
    %get the x,y coordinates of the current point:
    [x,y]=ind2sub(size(normedsd),e);
    if(e == 1);
        for next = (e+1):(dim^2)
            %get the coordinates of the other point
            [a,b]=ind2sub(size(normedsd),next);
            newh = sqrt(((x-a)^2)+((y-b)^2)); %Euclidean distance
            newsqrs = (normedsd(next)-normedsd(e))^2; %the difference between the two points, squared
            if newh<(maxh/2)
            h(1,index+1) = newh;
            sqrs(1,index+1) = newsqrs;
            index = index+1;
            end
        end
    else
        for next = (e+1):(dim^2) %repeat for all other e
            [a,b]=ind2sub(size(normedsd),next);
            newh = sqrt(((x-a)^2)+((y-b)^2)); %Euclidean distance
            newsqrs = (normedsd(next)-normedsd(e))^2; %the difference between the two points, squared
            if newh<(maxh/2)
            h(1,index+1) = newh;
            sqrs(1,index+1) = newsqrs;
            index = index+1;
            end
        end
        
    end
end

%get the semivariogram:
uniq = unique(h); %all the unique distances 
%uniq = uniq(uniq<(max(h)/2)); %use only those h < maxdistance/2 according to someone somewhere at one point.

semivar = zeros(1,length(uniq));
for u = 1:length(uniq)
    f = find(h == uniq(u));%find the position of all the entries equal to that particular distance
    if length(f) > 30; %this is one of the criteria for calculating a useful semivariogram via someone somewhere at one point.
    semivar(1,u) = (1/2)*(1/length(f))*sum(sqrs(f)); 
    else
        semivar(1,u) = NaN; %put NaN in for non-useful entries as a placeholder.
    end
end

%get the slope:
outputslope = regress(log(semivar(semivar~=0)'),[log(uniq(semivar~=0)') ones(size(uniq(semivar~=0)'))]);
%dis = fractalD(1,z);

%if saveit == 1
%filename = ['landscapevariables_fractaldim' num2str(fractaldims(d)) '_iteration#' num2str(z) '.mat'];
%save(filename);
end


