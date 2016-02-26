function W = computeW(is_old,imageX,dataW,emag,ephase)
% W = computeW(imageX,dataW,emag,ephase)
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
[p,q] = size(imageX);

fprintf('Getting indices...\n');


if is_old==1
    [w_i,w_j] = cimgnbmap([p,q],dataW.sampleRadius,dataW.sample_rate);
else
    [w_i,w_j] = my_cimgnbmap([p,q],dataW.sampleRadius,dataW.sample_rate,dataW.innerRadius);
end


fprintf('Computing affinities...\n');
W = affinityic(emag,ephase,w_i,w_j,max(emag(:)) * dataW.edgeVariance);
%W = affinityic(emag,ephase,w_i,w_j);
%W = my_affinityic(emag,ephase,w_i,w_j);
%W = my_affinityic(emag,ephase,w_i,w_j,max(emag(:)) * dataW.edgeVariance);
fprintf('Done building affinity matrix\n');

%W = W/max(W(:));