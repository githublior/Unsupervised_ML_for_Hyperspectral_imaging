function Af=SpatialFilter(A,Nf)

% Nf - filter half-size:
x=-Nf:Nf;
y=x';
X=ones(2*Nf+1,1)*x;
Y=y*ones(1,2*Nf+1);
F=exp(-2*(X.^2+Y.^2)/Nf^2);        % Gaussian filter definition
F=F/sum(sum(F));         

Af=conv2(A,F,'same');        % applying filter to Mask
end