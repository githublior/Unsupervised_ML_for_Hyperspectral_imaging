function XRGB=Spec2RGB(spec,lambda);

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


