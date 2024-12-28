function [xi, eta, w] = GaussTri(NI)
% Gaussian quadrature for triangle
% Return the quadrature points & weights on (0,1)
% Can only do 1,3,4 points quadrature
if NI == 1
    xi  = 1/3;
    eta = 1/3;
    w   = 1;
elseif NI == 3
    xi  = [1/6; 2/3; 1/6];
    eta = [1/6; 1/6; 2/3];
    w   = [1/3; 1/3; 1/3];
elseif NI == 4
    xi  = [1/3; 1/5; 1/5; 3/5];
    eta = [1/3; 3/5; 1/5; 1/5];
    w   = [-27/48; 25/48; 25/48; 25/48];
else
    error('Error: value of intergral points should be 1,3 or 4.');
end

% EOF