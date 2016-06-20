function [W,imageEdges,ephase] = ICgraph(I,edge_im,dataW,dataEdgemap);
% [W,imageEdges] = ICgraph(I,dataW,dataEdgemap);
% Input:
% I = gray-level image
% optional parameters: 
% dataW.sampleRadius=10;
% dataW.sample_rate=0.3;
% dataW.edgeVariance = 0.1;
% 
% dataEdgemap.parametres=[4,3, 21,3];%[number of filter orientations, number of scales, filter size, elongation]
% dataEdgemap.threshold=0.02;
% 
% Output: 
% W: npixels x npixels similarity matrix based on Intervening Contours
% imageEdges: image showing edges extracted in the image
%
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.


[p,q] = size(I);

%if (nargin< 4) | isempty(dataW),
    %dataW.sampleRadius=12;
    %dataW.sample_rate=0.1;
    
    %dataW.sampleRadius=15;
    %dataW.sample_rate=0.07;
    
    %dataW.innerRadius=7;
    %dataW.sampleRadius=10;
    %dataW.sample_rate=1;   

  
    
    dataW.edgeVariance = 0.1;
%end

%if (nargin<5) | isempty(dataEdgemap),
    dataEdgemap.parametres=[4,3, 21,3];%[number of filter orientations, number of scales, filter size, elongation]
    dataEdgemap.threshold=0.02;
%end



%edgemap = computeEdges(I,dataEdgemap.parametres,dataEdgemap.threshold);
%imageEdges = edgemap.imageEdges;

imageEdges=edge_im;




edgemap.emag=edge_im;
edgemap.ephase=randi(2,size(edge_im));


%disp(max(edgemap.emag(:)));
%disp(xy);

%imshow(uint8(edgemap.emag*255));
%disp(max(imageEdges(:)));
%disp(xy);
%edgemap.emag=double(edge_im);

% 
% dataW.innerRadius=3;
% dataW.sampleRadius=30;
% dataW.sample_rate=0.01;
% W = computeW(0,I,dataW,edgemap.emag,edgemap.ephase);
%[i1,j1,w1]=find(W1);
%W=W1;


dataW.sampleRadius=20;
dataW.sample_rate=0.1;

%dataW.sampleRadius=5;
%dataW.sample_rate=1;

ephase=edgemap.ephase;
W = computeW(1,I,dataW,edgemap.emag,edgemap.ephase);
%W=W2;
