function Z= spectoRGB(mymat,lambda, save_SpecToRgb,name);
%{
- mymat : mat of spectral img.
 lambda- wavelength list associated to mymat values.
-save_SpecToRgb: boolean. 1 if you want to save the output, 0 else
-name : string. name when saving output
%}
[a,b,c] = size(mymat);
Y = reshape(mymat, [a*b, c]);
sRGB = spectrumRGB(lambda); 
size(sRGB);
X = reshape(sRGB, [1*c,3]);

z =cast(Y,"double" )*X;
Z = reshape(z,[a,b, 3]);

rgbImg =  uint8(255 * mat2gray(Z));

if save_SpecToRgb==1
    nm =append(name,'_SpecToRgb', '.png' );
    imwrite(Z, nm);
end
%imshow(rgbImg);


fprintf('spectoRGB function terminated. \n');
