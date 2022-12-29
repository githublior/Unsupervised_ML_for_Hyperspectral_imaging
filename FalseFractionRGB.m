function FRGB=FalseFractionRGB(Fraction,W)

% W = weights for RGB false coloring
% e.g. W=[0,0,2,0; 0,0.5,0,0; 1,0,0,1]


FRGB(:,:,1)=W(1,1)*Fraction(:,:,1)+W(1,2)*Fraction(:,:,2)+W(1,3)*Fraction(:,:,3)+W(1,4)*Fraction(:,:,4);      % R
FRGB(:,:,2)=W(2,1)*Fraction(:,:,1)+W(2,2)*Fraction(:,:,2)+W(2,3)*Fraction(:,:,3)+W(2,4)*Fraction(:,:,4);      % G
FRGB(:,:,3)=W(3,1)*Fraction(:,:,1)+W(3,2)*Fraction(:,:,2)+W(3,3)*Fraction(:,:,3)+W(3,4)*Fraction(:,:,4);      % B
end
function XRGB=Spec2RGB(spec,lambda)

R=find(lambda>600);
G=find((lambda>500).* (lambda<600));
B=find(lambda<500);
Dim=size(spec);

XRGB(:,:,1)=sum(spec(:,:,R),3);
Map=squeeze(XRGB(:,:,1));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,1)=XRGB(:,:,1)./prctile(Sig,98);

XRGB(:,:,2)=sum(spec(:,:,G),3);
Map=squeeze(XRGB(:,:,2));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,2)=XRGB(:,:,2)./prctile(Sig,98);

XRGB(:,:,3)=sum(spec(:,:,B),3);
Map=squeeze(XRGB(:,:,3));
Sig=reshape(Map,Dim(1)*Dim(2),1);
XRGB(:,:,3)=XRGB(:,:,3)./prctile(Sig,98);
end