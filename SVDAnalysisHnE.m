function [Fraction,Err]=SVDAnalysisHnE(FImage,lambda,ProgRep,Ssmooth)

% ProgRep='t';         % type of progress report:  t=text, g=graph, other= none
% Ssmooth=1;           % Smoothing on WL basis befor analysis, the values is the sigma of a gaussian convoluted with image. can take any value. 0= no smoothing. 

%===================================
% SVD analysis
%===================================
if ProgRep=='n'
    disp('calculating SVD analysis for HnE');
end

load('REF_Absorption_Spectra.mat')
clear REF_Spec LA 

Absorption=-log10(FImage);

% PARAMETERS:
% focus data on interesting part ofg spectrum, e.g. (500-700nm):
LambdaMin=475;      % min/max wavelengths of interest for the analysis
LambdaMax=700;

Ns=4;                % number of input spectra, including Hematoxylin, Eosin and two for background
dn=1;                % sub-sampling in X,Y  to reduce amount of computation by dn^2, if needed
Nr=50000;            % progress report rate
%...........................................................................

Ny=size(Absorption,1);
Nx=size(Absorption,2);
NL=size(Absorption,3);


LIdx=(lambda>LambdaMin).*(lambda<LambdaMax);
Lmin=find(LIdx, 1 ); 
Lmax=find(LIdx, 1, 'last' );             
L=lambda(Lmin:Lmax)';
NL=Lmax-Lmin+1;          


REF_Spec(:,1)=L;  %wavelength
% Actuall absorption curves of Eosin:
REF_Spec(:,2)=interp1(Eosin_Abs_Rel(:,1),Eosin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Eosin:
% WLm2=520;  Sig2=10; 
% REF_Spec(:,2)=exp(-(REF_Spec(:,1)-WLm2).^2/(2*Sig2^2));

% Actuall absorption curves of Hematoxylin:
REF_Spec(:,3)=interp1(Hematoxylin_Abs_Rel(:,1),Hematoxylin_Abs_Rel(:,2),REF_Spec(:,1));
% Gaussian approximaltion to Hematoxylin:
% WLm3=590;  Sig3=40; 
% REF_Spec(:,3)=exp(-(REF_Spec(:,1)-WLm3).^2/(2*Sig3^2));

% Additional curves to account for background and unknown materials
WLm4=700;  Sig4=50; 
REF_Spec(:,4)=exp(-(REF_Spec(:,1)-WLm4).^2/(2*Sig4^2));
WLm5=460;  Sig5=30; 
REF_Spec(:,5)=exp(-(REF_Spec(:,1)-WLm5).^2/(2*Sig5^2));

% figure(5); plot(REF_Spec(:,1),REF_Spec(:,2),REF_Spec(:,1),REF_Spec(:,3))
% .............................................

A=REF_Spec(:,2:Ns+1);               % This is the transfer matrix from fractions to combined spectrum
[U,S,V] = svd(A);                        % SVD decomposition
S1=[inv(S(1:Ns,1:Ns)) zeros(Ns,NL-Ns)];                % this is the "inverse" of S

UnMix=inv(V')*S1*inv(U);            % this is generalized inv(A), allowing to go from spectra to fractions

%========================================================================================
% PREFILTERING:

LA=zeros(Ny,Nx,NL);

% Optional smoothing of the data (each WL seperately) in (X,Y) before analysis :
for j=1:NL
    if Ssmooth>0
        LA(:,:,j)=SpatialFilter(Absorption(:,:,Lmin+j-1),Ssmooth);    
    else
        LA(:,:,j)=Absorption(:,:,Lmin+j-1);
    end
end

% Optional subsampling:
if dn~=1   
    LA=LA(1:dn:end,1:dn:end,:);
end

%========================================================================================
% Running SVD analysis per spectrum:

Fraction=zeros(Ny,Nx,Ns);
Err=zeros(Ny,Nx);

n=0;
Tstart=now;

for j=1:Ny
    for k=1:Nx
        n=n+1;
        S=squeeze(LA(j,k,:));
        fr=UnMix*S;
        Fraction(j,k,:)=fr;                    % Fraction of each component is the output of analysis
        Err(j,k)=std(S-A*fr)^2;
        
        if fix(n/Nr)==n/Nr                     % progress report
            pFinish=100*n/(Nx*Ny);
            tPass= (now-Tstart)*24*60;
            tLeft = tPass*(100-pFinish)/pFinish;
            if ProgRep=='t'
                % progress report by text:
                disp([num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  ']);
            elseif ProgRep=='g'
                % progress report by graph:
                figure(22); clf;
                C1=A(:,1)*fr(1);
                C2=A(:,2)*fr(2);
                C3=A(:,3)*fr(3);
                C4=A(:,4)*fr(4);
                plot(L,S,'or', L,A*fr,'-b',L,C1,'--c',L,C2,'--c',L,C3,'--c',L,C4,'--c')
                drawnow
                title(['Error (SSE): ' num2str( Err(Ny,Nx)) ' ||  Finished: ' num2str(pFinish) ' % ||  Time from start: ' num2str(tPass ) ' min. ||  Remaining Time: ' num2str(tLeft)  ' min.  ']);
                pause(0.1)
%            elseif ProgRep=='n'
%                disp('calculating SVD analysis for HnE');
            end
        end
    end
end
end
