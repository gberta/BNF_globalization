function BNF_multi_class_demo()

    close all
    addpath('Ncut_9');

        
    lib_path='./libs/';

    setup_path=[lib_path '/vlfeat-0.9.18/toolbox/vl_setup'];

    run(setup_path);
 
    
    %% INPUT PARAMS %%

    type='aggressive'; % aggressive or non-agressive; agressive- slower, 
                       % assigns more pixels to the background, which often
                       % leads to better accuracy
    root_dir='multi_class_data/';
    

    img_dir=[root_dir 'images/'];
    unary_dir=[root_dir 'unary_data/']; %unary potential files must be stored as HxWxC matrix 
                             %H-height, W-width, C-number of classes where
                             %the matrix encodes probabilities for each
                             %class

    if strcmp(type,'aggressive')
        TH=1.5; %descides how aggressively to cut
                %higher threshold meeans less cutting
    end
    
    % Globalization parameters 
    tol=10^-6;
    max_it=200;
    mu=0.025; %lower mu puts more weight on pairwise cost, 
              %higher mu puts more cost on unary cost
    alpha=1/(1+mu);
    beta=mu/(1+mu);   

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sigma=30.0;


    %% MAIN CODE

    files=dir([img_dir '*.jpg']);


    for file_no=1:numel(files)

        cur_name=files(file_no).name;

        [r,nm,ext]=fileparts(cur_name);

        fprintf('Processing file %s %d / %d\n',nm,file_no,numel(files));

        im_file=strcat(img_dir,nm,'.jpg');
        unary_file=strcat(unary_dir,nm,'.mat');



        if exist(unary_file)


            fprintf('File no %d File nm %s\n',file_no,nm);

            im=read_img_rgb(im_file);

            h=size(im,1);
            w=size(im,2);

            n=h*w;
            
            edge_im=randi(2,[h w]);

            load(unary_file);

            %% ASSUMING THAT THE LOADED FILE WILL SHOW UP WITH THE FIELD 'data'
            ch=size(data,3);

            [V,I]=max(data,[],3);


            fprintf('Constructing Affinity Matrix...\n');

            [W,~]=get_my_W(im,edge_im);

            [ii,jj,~]=find(W);
            
            %% RGB AFFINITIES
            
            im_r=double(im(:,:,1))/sigma;
            im_g=double(im(:,:,2))/sigma;
            im_b=double(im(:,:,3))/sigma;

            v_r=(im_r(ii)-im_r(jj)).^2;
            v_g=(im_g(ii)-im_g(jj)).^2;
            v_b=(im_b(ii)-im_b(jj)).^2;
            
            v=exp(-(v_r+v_g+v_b));
            
            %% Building single affinity matrix

            v_unary=zeros(size(ii));

            I_ii=I(ii);
            I_jj=I(jj);
            pos_ind=find(I_ii==I_jj);

            pos_labels=I_ii(pos_ind);

            v_unary(I_ii~=I_jj)=0.001;

            y0_list=cell(ch,1);
            count=0;


            for c=1:ch


                [rr,cc]=find(I==c);

                if size(rr,1)>0

                    temp_ind=pos_ind(pos_labels==c);

                    class_fg=double(data(:,:,c));


                    y0_list{c}=class_fg(:);
                    count=count+1;

                    v_unary(temp_ind)=get_fc8_w(ii(temp_ind),jj(temp_ind),class_fg);

                end


            end


            %% GLOBALIZATION STEP

            %Combining Unary + Edge affinities; 
            %shouldnt be done if unaries are spatially disjoint !!!!

            v=exp(v_unary).*v;
            %v=v_ic;
            
            W=sparse(ii,jj,v,n,n);           
            d = sum(abs(W),2);
            D=spdiags(d,0,h*w,h*w);

            
            
            X=zeros(n,ch);
            XC=X;

            for c=1:ch


                [rr,cc]=find(I==c);

                %if 0
                if size(rr,1)>0    

                    fprintf('Channel %d\n',c);

                    class_fg=double(data(:,:,c));
                    class_bg=imcomplement(class_fg);


                    A=D-alpha*W;

                    b=beta*class_fg(:);                      

                    x=pcg(A,b,tol,max_it);             
                    X(:,c)=x; 

                    if strcmp(type,'aggressive')
                        b=beta*class_bg(:);  
                        xc=pcg(A,b,tol,max_it);           
                        XC(:,c)=xc;
                    end
                end
            end       






            %% Extracting the Predictions

            [~,S]=max(X,[],2); 
            if ~strcmp(type,'aggressive')
  
                S=reshape(S,[h w]);
                S=S-1;
            else
                bg_ind=find(S==1);
                F=zeros(n,ch);

                %% Aggressive Mode
                for c=2:ch
                    temp=[XC(:,c) TH*X(:,c)];
                    [~,pred]=max(temp,[],2);
                    pred=pred-1;
                    pos_ind=find(pred==1);
                    F(pos_ind,c)=TH*X(pos_ind,c);
                end

                [~,S]=max(F,[],2);
                S=reshape(S,[h w]);
                S(bg_ind)=0;
                S=max(S-1,0);
            end
            
            fprintf('Done\n');

            temp=cat(2,I-1,S);
            imshow(temp,colormap);
            %pause(3)
            %output_path=[output_dir nm '.png'];
            %imwrite(temp,colormap,output_path);


        end
    end

end


function v=get_fc8_w(ii,jj,fc8)
    v=max(abs(fc8(ii)-fc8(jj)),0.001);
    sigma=0.12;
    v=exp(-v/sigma);
end
