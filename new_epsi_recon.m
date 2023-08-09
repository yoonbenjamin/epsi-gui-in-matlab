function img = new_epsi_recon(study_details,file_path,ixofim1,lb)

% Parameters
% obtain number of echosimage size and image size

Ne = readprocpar(file_path,'ne'); Ne = Ne(2);
Nx = readprocpar(file_path,'np'); Nx = Nx(2)/2;
Ny = readprocpar(file_path,'nv'); Ny = Ny(2);
dte = readprocpar(file_path,'te2'); sw = 1/dte(2);
Neproc = Ne;

%% setup parameters
n1 = str2double(file_path(end-1:end))+(ixofim1-1);
% PE order
if study_details.centric
    pe = [0 -1 1 -2 2 -3 3 -4 4 -5 5 -6]+7;
else 
    pe = 1:Ny;
end 
% sets variables to fill later
kspace1 = zeros(Ny,Nx,Ne,study_details.nimg_to_process);
kspace1_2=kspace1;
kspace2 = kspace1;
kspace = zeros(Ny,Nx,2*Neproc,study_details.nimg_to_process);
img = kspace;
img2 = img;

%% Reconstruct EPSI Images
% file name loop
for ii = 1:study_details.nimg_to_process
    if n1+ii-1 < 10
       path = strcat(file_path(1:end-1),num2str(n1+ii-1));
    else
       path = strcat(file_path(1:end-2),num2str(n1+ii-1));       
    end
    %path = strcat(file_path(1:end-1),num2str(0+ii));
    % load data

    [RE,IM,~,~,~] = varianloadfid(path);

    % arrange echoes
    for jj = 1:Ne
        idx = 1+(jj-1):Ne:Ny*Ne;
        kspace1(:,:,jj,ii) = [RE(:,idx) + 1i*IM(:,idx)]';

        % this is because real x imaginary doesn't have right orientation
    end
    
    % reorient kspace
    for jj=1:Ny
        kspace1_2(pe(jj),:,:,:)=kspace1(jj,:,:,:);
    end 

    % apply line broadening
    for jj = 1:Ny
        for kk = 1:Nx
            line = squeeze(kspace1_2(jj,kk,:,ii));
            kspace2(jj,kk,:,ii) = line.*exp(-[1:Ne]'*lb/sw);
        end
    end
    
    kspace(:,:,1:Neproc,ii) = kspace2(:,:,1:Neproc,ii);
    img3 = fftshift(fftn(squeeze(kspace(:,:,:,ii))));
    img1=img3;

    % flip spectral x axis 
    for jj = 1:Ny
        for kk = 1:Nx
            img2(jj,kk,:,ii) = img1(jj,Nx-kk+1,end:-1:1);
        end     
    end
end
img=abs(img2);