function [D,DW]= get_my_D(W);

if nargin < 2
    nbEigenValues = 8;
end
if nargin < 3
    dataNcut.offset = 5e-1;
    dataNcut.verbose = 0;
    dataNcut.maxiterations = 300;
    dataNcut.eigsErrorTolerance = 1e-8;
    dataNcut.valeurMin=1e-6;
end
% if nargin < 3
%     dataNcut.offset = 5e-1;
%     dataNcut.verbose = 0;
%     dataNcut.maxiterations = 100;
%     dataNcut.eigsErrorTolerance = 1e-6;
%     dataNcut.valeurMin=1e-6;
% end

% make W matrix sparse
%fprintf('Sparsifying W\n');
%W = sparsifyc(W,dataNcut.valeurMin);

% check for matrix symmetry
if max(max(abs(W-W'))) > 1e-10 %voir (-12) 
    disp(max(max(abs(W-W'))));
    error('W not symmetric');
end

n = size(W,1);
offset = dataNcut.offset;


% degrees and regularization
fprintf('Degrees & Regularization\n');
d = sum(abs(W),2);
dr = 0.5 * (d - sum(W,2));
d = d + offset * 2;
dr = dr + offset;
W = W + spdiags(dr,0,n,n);
fprintf('Done\n');

D=spdiags(d,0,n,n);
DW=D-0.99*W;

end
