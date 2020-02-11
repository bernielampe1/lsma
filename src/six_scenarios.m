% Load signature for simulation scenarios
% Enter 1 for Reflectance signature; 2 for Radiance signature
% ImageCub1 for Scenario 1
% ImageCub2 for Scenario 2
% ImageCub3 for Scenario 3
% ImageCub4 for Scenario 4
% ImageCub5 for Scenario 5
% ImageCub6 for Scenario 6
function [TI1, TI2, TI3, TE1, TE2, TE3, M, GTmap] = six_scenarios()

    ind1=input('Load simulated signature: 1:Reflectance, 2:Radiance?');
    while ind1~=1 && ind1~=2
        ind1=input('Load simulated signature: 1:Reflectance, 2:Radiance?');
    end
    if ind1==1
        load Cuprite_simulate_reflectance
    else
        load Cuprite_simulate_radiance
    end
    GTmap=zeros(200);   % Ground Truth Map

    % simulated white Gaussian noise 
    % snr is 50% of average signal divide by standard deviation of noise
    snr=20;
    sigma=0.5*BGsig/snr;
    noise=reshape(repmat(sigma,1,200*200)',200,200,size(MM,1)).*randn(200,200,size(MM,1));

    count=1;

    % scenario 1: clean target , clean background
    ImageCub=reshape(repmat(BGsig,1,200*200)',200,200,size(MM,1));
    for k=1:5        % Each row, one endmember will be the main endmember,
        d=MM(:,k);   % 1:Alunite, 2:Buddingtonite, 3:Calcite, 4:Kaolinite, 5:Muscovite
        % pure pixels in column 1 (4x4)
        j=1;
        for m=1:4
            for n=1:4
                ImageCub(30+30*(k-1)+m,30+30*(j-1)+n,:)=d;
                GTmap(30+30*(k-1)+m,30+30*(j-1)+n,:)=1;
                P(count,:)=[30+30*(j-1)+n,30+30*(k-1)+m];
                str=sprintf('LABEL{count}=''p%d_%d%d%d'';',k,j,m,n);
                eval(str)
                count=count+1;
            end
        end

        % pure pixels in column 2 (2x2)
        j=2;
        for m=1:2
            for n=1:2
                ImageCub(30+30*(k-1)+m,30+30*(j-1)+n,:)=d;
                GTmap(30+30*(k-1)+m,30+30*(j-1)+n,:)=1;
                P(count,:)=[30+30*(j-1)+n,30+30*(k-1)+m];
                str=sprintf('LABEL{count}=''p%d_%d%d%d'';',k,j,m,n);
                eval(str)
                count=count+1;
            end
        end

        % mixed pixel in column 3 (50% of 2 endmembers, 2x2)
        U=MM;
        U(:,k)=[];
        MM1=0.5*(repmat(d,1,4)+U);
        j=3;
        for m=1:2
            for n=1:2
                ImageCub(30+30*(k-1)+m,30+30*(j-1)+n,:)=MM1(:,2*(m-1)+n);
                GTmap(30+30*(k-1)+m,30+30*(j-1)+n,:)=1;
                P(count,:)=[30+30*(j-1)+n,30+30*(k-1)+m];
                str=sprintf('LABEL{count}=''p%d_%d%d%d'';',k,j,m,n);
                eval(str)
                count=count+1;
            end
        end

        % mixed pixel in column 4 (50% of endmember + 50% of Background, 1x1)
        j=4;
        ImageCub(30+30*(k-1)+1,30+30*(j-1)+1,:)=0.5*d+0.5*BGsig;
        GTmap(30+30*(k-1)+1,30+30*(j-1)+1,:)=1;
        P(count,:)=[30+30*(j-1)+1,30+30*(k-1)+1];
        str=sprintf('LABEL{count}=''p%d_%d1'';',k,j);
        eval(str)
        count=count+1;

        % mixed pixel in column 5 (25% of endmember + 75% of Background, 1x1)
        j=5;
        ImageCub(30+30*(k-1)+1,30+30*(j-1)+1,:)=0.25*d+0.75*BGsig;
        GTmap(30+30*(k-1)+1,30+30*(j-1)+1,:)=1;
        P(count,:)=[30+30*(j-1)+1,30+30*(k-1)+1];
        str=sprintf('LABEL{count}=''p%d_%d1'';',k,j);
        eval(str)
        count=count+1;
    end
    ImageCub1=ImageCub;

    % scenario 2: clean target, noisy background
    % remove noise from target area (noise1)
    for i=1:size(MM,1)
        noise1(:,:,i)=noise(:,:,i).*(1-GTmap);
    end
    ImageCub2=ImageCub1+noise1;  % add noise

    % scenario 3: noisy target, noisy background
    ImageCub3=ImageCub1+noise;   % add noise to whole scene

    % scenario 4: embeded clean target, clean background
    % similar to scenario 1, add target to background instead of insert to it
    ImageCub=reshape(repmat(BGsig,1,200*200)',200,200,size(MM,1));
    for k=1:5
        d=MM(:,k);
        j=1;
        for m=1:4
            for n=1:4
                ImageCub(30+30*(k-1)+m,30+30*(j-1)+n,:)=BGsig+d;
                GTmap(30+30*(k-1)+m,30+30*(j-1)+n,:)=1;
            end
        end

        j=2;
        for m=1:2
            for n=1:2
                ImageCub(30+30*(k-1)+m,30+30*(j-1)+n,:)=BGsig+d;
                GTmap(30+30*(k-1)+m,30+30*(j-1)+n,:)=1;
            end
        end

        U=MM;
        U(:,k)=[];
        MM1=0.5*(repmat(d,1,4)+U);
        j=3;
        for m=1:2
            for n=1:2
                ImageCub(30+30*(k-1)+m,30+30*(j-1)+n,:)=BGsig+MM1(:,2*(m-1)+n);
                GTmap(30+30*(k-1)+m,30+30*(j-1)+n,:)=1;
            end
        end

        j=4;
        ImageCub(30+30*(k-1)+1,30+30*(j-1)+1,:)=BGsig+0.5*d+0.5*BGsig;
        GTmap(30+30*(k-1)+1,30+30*(j-1)+1,:)=1;

        j=5;
        ImageCub(30+30*(k-1)+1,30+30*(j-1)+1,:)=BGsig+0.25*d+0.75*BGsig;
        GTmap(30+30*(k-1)+1,30+30*(j-1)+1,:)=1;

    end
    ImageCub4=ImageCub;

    % scenario 5: embeded clean target, noisy background
    for i=1:size(MM,1)
        noise1(:,:,i)=noise(:,:,i).*(1-GTmap);
    end
    ImageCub5=ImageCub4+noise1;

    % scenario 6: embeded noisy target, noisy background
    ImageCub6=ImageCub4+noise;

    % % show ground truth map and image at band 80 for all scenarios
    % figure,imagesc(GTmap);colormap(gray);title('Ground Truth Map')
    % figure,imagesc(ImageCub1(:,:,80));colormap(gray);title('Scenario 1')
    % figure,imagesc(ImageCub2(:,:,80));colormap(gray);title('Scenario 2')
    % figure,imagesc(ImageCub3(:,:,80));colormap(gray);title('Scenario 3')
    % figure,imagesc(ImageCub4(:,:,80));colormap(gray);title('Scenario 4')
    % figure,imagesc(ImageCub5(:,:,80));colormap(gray);title('Scenario 5')
    % figure,imagesc(ImageCub6(:,:,80));colormap(gray);title('Scenario 6')
    TI1 = ImageCub1;
    TI2 = ImageCub2;
    TI3 = ImageCub3;
    TE1 = ImageCub4;
    TE2 = ImageCub5;
    TE3 = ImageCub6;
    M = MM;
end