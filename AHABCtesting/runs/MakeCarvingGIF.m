clear functions
f1 = @Function2DRealAsBoxesPlot;

f = figure(60)
clf(f);
%axis tight
set(gca,'NextPlot','replacechildren');


filename = 'carving.gif';

for k = 1:62 
    
    boxesFileNameBase = 'CarverSEBQueueState_';
    boxesFileName = strcat(boxesFileNameBase, int2str(k),'.txt');
    
    clf(f)
    h1 = gca;
    set(h1,'NextPlot','replace');
    
    p = f1(boxesFileName, h1);
    %set(get(h1,'XLabel'),'String',texlabel('theta'));
    %set(get(h1,'YLabel'),'String',texlabel('h'));
    %set(h1,'View',[61 18]); % good for untransformed data
    set(h1,'View',[21 32]); 
    da = daspect;
    if da(1) > da(2)
        da(2) = da(1);
    elseif da(2) > da(1)
        da(1) = da(2); 
    end

    daspect(h1, da);
    %pause(1);
    
   frame = getframe(gcf);
   
   im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
    
