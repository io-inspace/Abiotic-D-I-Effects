%get the proportion of suitable habitat for each sps in each landscape:
function [breakpoints,reprohab] = getbreakpoints(dim,maxsp,sigsq,Mp,Ma,Rsuit_a,Rsuit_p,numls)

load(['setsofuncorrelatedlandscapes' num2str(dim)],'all_landscapes')

breakpoints = zeros(numls,(maxsp-1));

for t = 1:numls;

landsc = all_landscapes{t,1};

B = landsc{1,2};
  
Brange = ((max(max(B)))-(min(min(B))));
insetB = Brange/maxsp;
rangefordiff = (max(max(B))-insetB)-(min(min(B))+insetB); 
goby = rangefordiff/(maxsp-1);
npdiffB = (min(min(B))+insetB):goby:(max(max(B))-insetB);    

%how do species partition resources? 
for bp = 1:(maxsp-1)
syms g
eqn = exp((-.5)*(npdiffB(bp)-(g))^2/(sigsq)) == exp((-.5)*(npdiffB(bp+1)-(g))^2/(sigsq));
breakpoints(t,bp) = double(solve(eqn));
end

end

%this is just a sanity check that the good habitat is always the same
%number
seedbreak_a = zeros(1,2);
seedbreak_p = zeros(1,2);
syms r
eqn_a = Ma*exp((-.5)*(0-(r))^2/(sigsq))==Rsuit_a;
eqn_p = Mp*exp((-.5)*(0-(r))^2/(sigsq))==Rsuit_p;
abreak = double(solve(eqn_a,r));
seedbreak_a(1,1) = abreak(1);
seedbreak_a(1,2) = abreak(2);
pbreak = double(solve(eqn_p,r));
seedbreak_p(1,1) = pbreak(1);
seedbreak_p(1,2) = pbreak(2);
reprohab = seedbreak_p;

end