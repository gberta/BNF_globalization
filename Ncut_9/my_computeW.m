function W = my_computeW(imageX,dataW,emag,ephase)
% W = computeW(imageX,dataW,emag,ephase)
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
%[p,q] = size(imageX);

fprintf('Getting local pairs...\n');
[ii,jj]=my_getLocalPairs(size(emag),[],[],[]);
fprintf('Size local pairs %d\n',size(ii,1));

emag=permute(emag,[2 1 3]);
ephase=permute(ephase,[2 1 3]);

preds = my_affinityic(emag,ephase,ii,jj,max(emag(:)) * dataW.edgeVariance);

%TEST CASE
%temp=permute(imageX,[2 1 3]);
%preds=temp(ii)==temp(jj);
%%%%%%

fprintf('Creating sparse matrix \n');
n_pixels=size(emag,1)*size(emag,2);
W=sparse(double(ii),double(jj),double(preds),n_pixels,n_pixels);

%W = W/max(W(:));
fprintf('Making matrix symmetric...\n');
W=(W+W');
W=W/max(W(:));

fprintf('Done\n')
