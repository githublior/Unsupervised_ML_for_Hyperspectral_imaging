function [IndAll,Spectra]=ShowSpec2(RGBHE,spec,lambda,n,Fraction)
% RGBHE - RGB Image for the original H&E transmission
% spec - data cube
% labmbda ...
% n - number of nearest neibors for averaging when displaying spectrum.

FRGB=FalseFractionRGB(Fraction,[0,0,3,0; 0,0.7,0,0; 1,0,0,1]);      % false RGB representing SVD analysis results (Fraction)

f=figure(6);
if isempty(f.CurrentAxes)               % if figure exist no need to redraw so magnif remains as before
    subplot(2,2,1); 
    imshow(RGBHE)  
    a=f.CurrentAxes;
    subplot(2,2,2);
    imshow(FRGB)
    b=f.CurrentAxes;
end

a.Position=[0.0    0.48    0.6    0.5];
b.Position=[0.5    0.48    0.6    0.5];

IndAll=[];
Spectra=[];

stop=0;

while stop~=1
    [x,y,button]=ginput(1);
    X=round(x);
    Y=round(y);
    Ind=[X Y];
    if button=='s' || button=='S'       % 's' for stopping
        stop=1; 
    elseif button=='p' || button=='P'   % 'p' may be used for pausing and changing magnification, hit Enter when done 
        pause
        % adjust magnif similar for both views:
        c=f.CurrentAxes;    % read imag limits from last touched subplot
        XL=c.XLim;          
        YL=c.YLim;
        a.XLim=XL;          % apply same limits to both plots
        a.YLim=YL;
        b.XLim=XL;
        b.YLim=YL;
    elseif button=='n' || button=='N'   % 'n' for changing the averaging area
        n=input('Averaging area size? : ');
    else                                % mouse-click - choosing a point a drawing spectrum
        subplot(2,2,3); 
        a3=f.CurrentAxes;
        a3.Position=[0.1300    0.0800    0.3347    0.3];
        if n==0
            Spectrum=squeeze(spec(Y,X,:));
            F1=Fraction(Y,X,1);
            FE=Fraction(Y,X,2);
            FH=Fraction(Y,X,3);
            F4=Fraction(Y,X,4);
        else
            Spectrum=squeeze(mean(spec(Y-n:Y+n,X-n:X+n,:),[1 2]));
            F1=mean(Fraction(Y-n:Y+n,X-n:X+n,1),[1 2]);
            FE=mean(Fraction(Y-n:Y+n,X-n:X+n,2),[1 2]);
            FH=mean(Fraction(Y-n:Y+n,X-n:X+n,3),[1 2]);           
            F4=mean(Fraction(Y-n:Y+n,X-n:X+n,4),[1 2]);           
        end
        plot(lambda,Spectrum)
        axis([400 800 0 1.2])
        title(['X:  ' num2str(X) ' , Y:  ' num2str(Y)])
        IndAll=[IndAll; Ind];
        Spectra=[Spectra; Spectrum];
        
        subplot(2,2,4)
        a4=f.CurrentAxes;
        a4.Position=[0.5703    0.1100    0.3347    0.3412];

        bar(1:4,[F1 FE FH F4])
        axis([0 5 -.2 1.2])
        title(['F1:  ' num2str(F1,'%.3f')   '      , Eosin:  ' num2str(FE,'%.3f') '       ,  Hematoxylin:  ' num2str(FH,'%.3f') '      , F4:  ' num2str(F4,'%.3f') ])
        text(1,1,['n= ' num2str(n)])
  
    end
   
end

end

