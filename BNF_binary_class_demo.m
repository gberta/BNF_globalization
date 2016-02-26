function BNF_binary_class_demo()

    close all
    addpath('Ncut_9');
    
    
    lib_path='./libs/';

    setup_path=[lib_path '/vlfeat-0.9.18/toolbox/vl_setup'];

    run(setup_path);

    run('vl_setupnn.m');

    %% INPUT PARAMS %%

    root_dir='binary_class_data/';
                       
    edge_dir=[root_dir 'edges/'];
    img_dir=[root_dir 'images/'];
    unary_dir=[root_dir 'unary_data/']; %unary potential files must be stored as 
                             %HxW matrix with foreground probabilities
                             %H-height, W-width

        
    sp_size=10;
    TH=1; %descides how aggressively to cut
            %higher threshold meeans less cutting
  


    % Globalization parameters 
    tol=10^-6;
    max_it=200;
    mu=0.025; %lower mu puts more weight on pairwise cost, 
              %higher mu puts more cost on unary cost
    alpha=1/(1+mu);
    beta=mu/(1+mu);   

    %%%%%%%%%%%%%%%%%%%%%%%%%%


    %% MAIN CODE


    files=dir([img_dir '*.jpg']);

    for file_no=4
    %for file_no=1:numel(files)

        cur_name=files(file_no).name;

        [r,nm,ext]=fileparts(cur_name);

        fprintf('Processing file %s %d / %d\n',nm,file_no,numel(files));

        im_file=strcat(img_dir,nm,'.jpg');
        edge_file=strcat(edge_dir,nm,'.jpg');
        unary_file=strcat(unary_dir,nm,'.mat');



        if exist(unary_file) && exist(edge_file)


            fprintf('File no %d File nm %s\n',file_no,nm);


            %im=imread(im_file);  
            im=read_img_rgb(im_file);
            
            sp_info=gen_supperpixel_info(im, sp_size);
            sp_map=sp_info.sp_ind_map;
        

            h=size(im,1);
            w=size(im,2);

            n=h*w;

            edge_im=im2double(imread(edge_file));

            load(unary_file);
    

            %% ASSUMING THAT THE LOADED FILE WILL SHOP UP WITH THE FIELD 'data'
            I=data>0.5;

            fprintf('Constructing Affinity Matrix...\n');

            [W,~]=get_my_W(im,edge_im);

            [ii,jj,v_ic]=find(W);
            
            v=ones(size(ii,1),1);
            v(sp_map(ii)~=sp_map(jj))=v_ic(sp_map(ii)~=sp_map(jj));


            %% Building single affinity matrix

            
            %I_ii=I(ii);
            %I_jj=I(jj);
            
            v_unary=get_fc8_w(ii,jj,data);
      


            %% GLOBALIZATION STEP

            %Combining Unary + Edge affinities; 
            %shouldnt be done if unaries are spatially disjoint !!!!

            %v=exp(v_unary).*v;
            
            W=sparse(ii,jj,v,n,n);           
            d = sum(abs(W),2);
            D=spdiags(d,0,h*w,h*w);


            class_fg=double(data);
            class_bg=imcomplement(class_fg);


            A=D-alpha*W;

            b=beta*class_fg(:);                      

            x=pcg(A,b,tol,max_it); 

    
            b=beta*class_bg(:);  
            xc=pcg(A,b,tol,max_it);  
            
            X=[xc TH*x];
            [~,S]=max(X,[],2); 
            


            %% Extracting the Predictions
            
            S=reshape(S,[h w]);
            
            
            fprintf('Done\n');
            
            figure()
            imagesc(I)

            figure()
            imagesc(S-1)
            
            
            %disp(xy);


        end
    end

end


function v=get_fc8_w(ii,jj,fc8)
    v=max(abs(fc8(ii)-fc8(jj)),0.001);
    sigma=0.12;
    v=exp(-v/sigma);
end
