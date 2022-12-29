function grayscale_out = read_unique_wl(wl_mat , lambda_mat, chosen_wl);

wl1 = wl_mat(:,:,chosen_wl);
%[a,b] = size(wl1);
lb = lambda_mat(chosen_wl);
Y = spectrumRGB(lb);
B = reshape(Y, 1, 1,[]);
%size(B)
y = rgb2gray(B);
%size(y)
z = wl1 * y ;
grayscale_out = 255 * mat2gray(z);

