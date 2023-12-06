% Redefing maxInt filters for the B1AR data
function [B1AR,maxInt] = DensityPlotting(B1AR,dz)
    %Definition of MaxInt and ZPOS! :)
    maxInt = imgaussfilt(B1AR.outputs.maxInt,1);
    zPos = B1AR.results.zValueAvg;
    
    %making masks filters based on fitting parameters
    % G and I decided that this has to do xD
    sigMask = maxInt < median(maxInt) & B1AR.outputs.sig.*sqrt(8*log(2)).*dz > 1200 | maxInt > median(maxInt) & B1AR.outputs.sig.*sqrt(8*log(2)).*dz > 900 | B1AR.outputs.sig.*sqrt(8*log(2)).*dz < 200 ;
    maxInt(sigMask) = NaN;
    zPos(sigMask) = NaN;
% %     h1 = figure 
% %     %imagesc(B1AR.B1AR.outputs.sig.*sqrt(8*log(2)).*dz)
% %     imagesc(sigMask)
% %     xlim([157 205]), ylim([246 294])
    
    
    %Raquare mask
    Rsquare = B1AR.outputs.adjrsquare < 0.9 | B1AR.outputs.adjrsquare > 1;
    maxInt(Rsquare) = NaN;
    zPos(Rsquare) = NaN;
% %     h1 = figure 
% %     %imagesc(B1AR.B1AR.outputs.sig.*sqrt(8*log(2)).*dz)
% %     imagesc(Rsquare)
% %     xlim([157 205]), ylim([246 294])
% %     
    %zPosStd
    zPosStd = B1AR.outputs.zPosStd.*dz < 0 | B1AR.outputs.zPosStd.*dz > 120;
    maxInt(zPosStd) = NaN;
    zPos(zPosStd) = NaN;
% %     h1 = figure 
% %     imagesc(B1AR.outputs.zPosStd.*dz)
% %     %imagesc(zPosStd)
% %     xlim([157 205]), ylim([246 294])
    
    
    %maxIntError
    MaxIntError = B1AR.outputs.maxIntStd./B1AR.outputs.maxInt;
    maxIntError = MaxIntError < 0 | MaxIntError > 0.3;
    maxInt(maxIntError) = NaN;
    zPos(maxIntError) = NaN;
% %     h1 = figure 
% %     %imagesc(MaxIntError)
% %     imagesc(maxIntError)
% %     xlim([157 205]), ylim([246 294])
    
    B1AR.normDen = maxInt./B1AR.confocalAreas;
    B1AR.zValue = zPos;
    
    
% %     X = (B1AR.results.x.*30)./1000;
% %     Y = (B1AR.results.y.*30)./1000;
% %     Z = B1AR.results.zValueAvg;%-1400;
% %     %colormaps
% %     m=100;
% %     cm_inferno=plasma(m);%B1ARdensity
% %     CT=cbrewer('div', 'RdYlBu', 150);%Curvatures
% %     Nb80map = load('Nb80cm');
    
% %     
% %     h1 = figure 
% %     surf(X,Y,Z-1150,B1AR.B1AR.normMaxInt.*100,'EdgeColor','interp')
% %     xlim([157*30/1000 205*30/1000]), ylim([246*30/1000 294*30/1000]), zlim([0 1350-1150]),
% %     view(47,51)
% %     colormap(cm_inferno),
% %     c = colorbar; c.LineWidth = 0.5; c.TickDirection = 'out';
% %     ylabel(c,'Density B1AR [A.U.]'),
% %     set(gca,'FontSize',9,'linewidth',0.5,'TickDir','out'),%,'Ydir','normal','Xdir','reverse')
% %     set(gcf,'color','w', 'Units', 'centimeters', 'Position', [0, 0, 10, 8], 'PaperUnits', 'centimeters', 'PaperSize', [10, 8])
% %     xlabel('x [\mum]'), ylabel('y [\mum]'), zlabel('z [nm]')
% %     ax = gca;
% %     ax.Position = ax.Position + [0 0 0.1 0];
% %     c.Position = c.Position + [-0.05 0 -0.05 0];
% %     caxis([3 6.5]),
% %     
% %     h1 = figure 
% %     surf(X,Y,Z-1150,maxInt,'EdgeColor','interp')
% %     xlim([157*30/1000 205*30/1000]), ylim([246*30/1000 294*30/1000]), zlim([0 1350-1150]),
% %     view(47,51)
% %     colormap(cm_inferno),
% %     c = colorbar; c.LineWidth = 0.5; c.TickDirection = 'out';
% %     ylabel(c,'Density B1AR [A.U.]'),
% %     set(gca,'FontSize',9,'linewidth',0.5,'TickDir','out'),%,'Ydir','normal','Xdir','reverse')
% %     set(gcf,'color','w', 'Units', 'centimeters', 'Position', [0, 0, 10, 8], 'PaperUnits', 'centimeters', 'PaperSize', [10, 8])
% %     xlabel('x [\mum]'), ylabel('y [\mum]'), zlabel('z [nm]')
% %     ax = gca;
% %     ax.Position = ax.Position + [0 0 0.1 0];
% %     c.Position = c.Position + [-0.05 0 -0.05 0];
% %     caxis([30 65]),
    
    
end
