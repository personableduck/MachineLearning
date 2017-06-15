function [ E ] = entropy( prof )
%entropy: calculates the Renyi entropy of the input profile
%         When r=2, Renyi entropy is a measure of the profile sharpness.
%         It increases when noise increases and when the sharpness gets
%         lower (e.g. as a consequence of inflated lobes due to resolution loss)
%Input:   profile/plot vector;
%Output:  Renyi entropy (default parameter value: r=2)
r=2;
E=1/(1-r)*log(sum((prof.^(2*r))))-(r/(1-r)*log(sum((prof.^2))));

end
