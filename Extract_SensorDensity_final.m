 %% Extracting density conformational biosensor

%% Initialize pixel sizes in nm
dx=30;
dy=30; 
dz=30;

%% Loading of data
APO.B1AR=B1AR;
APO.GFP=miniG;
clear B1AR miniG
noCells=length(APO.B1AR);
%% Find the Cytosolic contribution

for p=1:noCells
    figure, 
    BG=APO.GFP(p).avgImageData(:,:,70:end);
    histogram(BG(:),'BinWidth',1,'Normalization','probability'), 
    hold on, 
    line([mode(round(BG(BG>mean(BG)-0.1*mean(BG)))) mode(round(BG(BG>mean(BG)-0.1*mean(BG))))], [0 , 0.2]),
    [xi,yi]=getpts;        
    Icyt(p)=round(2*xi)/2; %Rounding to 0.5
end

%% Determine the FitType
for p=1:noCells
    data = APO.GFP(p).avgImageData;
    dataList = nan(size(data,3), size(data,1)*size(data,2));
    BG=APO.GFP(p).avgImageData(:,:,70:end);
    IcytStd=Icyt(p)*0.1; % 10% on both sides
    
    for frame = 1:size(data,3)
        frameData = data(:,:,frame);
        dataList(frame, :) = frameData(:);
    end
    zData(:,1) = 1:size(data,3);
    
    zValueR=round(APO.B1AR(p).results.zValueAvg(:)./30);
    zValue=APO.B1AR(p).results.zValueAvg(:)./30;
    fitType=zeros(size(zValue,1),5);
    index=nan(size(zValue));
    slopeAll=nan(size(zValue));
    slopeLast=nan(size(zValue));
    sigma=APO.B1AR(p).outputs.sig(:);
    
    rsq=nan(size(zValue));
    aCoef=nan(size(zValue));
    bCoef=nan(size(zValue));
    cCoef=nan(size(zValue));
    dCoef=nan(size(zValue));
    eCoef=nan(size(zValue));
    fCoef=nan(size(zValue));
    gCoef=nan(size(zValue));
    kCoef=nan(size(zValue));
    
    aCoefErr=nan(size(zValue));
    bCoefErr=nan(size(zValue));
    cCoefErr=nan(size(zValue));
    dCoefErr=nan(size(zValue));
    eCoefErr=nan(size(zValue));
    fCoefErr=nan(size(zValue));
    gCoefErr=nan(size(zValue));
    kCoefErr=nan(size(zValue));
    
    for i=1:size(dataList,2)
        if isnan(zValueR(i))==1 || isnan(sigma(i))==1
            continue
        end
        
        sub=[];
        norm_trace=[];
        norm_trace=(dataList(:,i)-(mean(mink(dataList(:,i),10))))./(mean(maxk(dataList(:,i),15))-std(maxk(dataList(:,i),15))-(mean(mink(dataList(:,i),10))));
        std_norm_trace=movstd(norm_trace,10);
        mean_norm_trace=movmean(norm_trace,10);
        [~,I]=mink(abs(0.5-mean_norm_trace(1:zValueR(i)+5)),7);
        sub=[(zValueR(i)+10 ) > I];
        index(i)=I(find(sub,1,'first'));
        
        Bmax(i)=mean(dataList(85:end,i));
        
        if zValue(i)+sigma(i)>length(zData)
            continue
        end 

        z_thresh=6; 
        [pFitLast,SLast] = polyfit(85:100,mean_norm_trace(85:end)',1); 
        [pFitAll,SAll] = polyfit(round(zValue(i)+1.5.*sigma(i)):85,mean_norm_trace(round(zValue(i)+1.5*sigma(i)):85)',1); 
        [y_fitLast,~] = polyval(pFitLast,85:100,SLast);
        [y_fitAll,~] = polyval(pFitAll,round(zValue(i)+1.5*sigma(i)):85,SAll);
        slopeLast(i)=pFitLast(1);       
        slopeAll(i)=pFitAll(1);         


if  abs(zValue(i)-index(i))<z_thresh && abs(slopeLast(i))<(0.02/15) && abs(slopeAll(i))<(0.02/15) && (Icyt(p)-IcytStd)<Bmax(i) && Bmax(i)<(Icyt(p)+sqrt(Icyt(p)))
    fitType(i,1)=1;
else
    fitType(i,1)=0;
end

if (zValue(i)-index(i))>z_thresh && abs(slopeLast(i))<(0.1/15) && (Icyt(p)-IcytStd)<Bmax(i) && Bmax(i)<(Icyt(p)+sqrt(Icyt(p)))
    fitType(i,2)=2;
    else
    fitType(i,2)=0;
end

if (zValue(i)-index(i))<z_thresh && Bmax(i)<=(Icyt(p)+sqrt(Icyt(p)))  && (abs(slopeAll(i))>(0.02/15) || abs(slopeLast(i))>(0.02/15))
    fitType(i,3)=3;
     else
    fitType(i,3)=0;
end

if  (abs(slopeAll(i))>(0.02/15) || abs(slopeLast(i))>(0.02/15)) && Bmax(i)<=(Icyt(p)+sqrt(Icyt(p)))
    fitType(i,4)=4;
     else
    fitType(i,4)=0;
end

if fitType(i,3)==3 && fitType(i,4)==4
   fitType(i,4)=0;
end 

if fitType(i,2)==2 && fitType(i,4)==4
   fitType(i,4)=0;
   %fitType(i,2)=0;
end

if sum(fitType(i,1:4),2)==0 
   fitType(i,5)=5;
end

fitTypeAll=sum(fitType,2);

    end

APO.Fit(p).fitTypeAll=sum(fitType,2);

end
     
%% Removing single pixels and cleaning up the data.

se=strel('square',1);
for p=1:noCells  
APO.Fit(p).MapFitType=reshape(APO.Fit(p).fitTypeAll, [size(APO.B1AR(p).results.zValueAvg,1), size(APO.B1AR(p).results.zValueAvg,2)]);
mask1=[APO.Fit(p).MapFitType==2 | APO.Fit(p).MapFitType==4]; % choose all pixels with a Gaussian contribution
mask2=bwareaopen(mask1,16);
FitTypeAdj=APO.Fit(p).MapFitType.*mask2;
FitTypeAdj(FitTypeAdj==0)=NaN;
mask3=movmean(FitTypeAdj, [2 2],'omitnan');
mask4=zeros(size(APO.B1AR(p).results.zValueAvg));
mask4((mask3-3)>=0)=4;
mask4((mask3-3)<0)=2;

APO.Fit(p).AdjfitTypeAll=mask4(:);

clear mask1 mask2 mask3 
end


%% Fitting after FiType has been determined

% Condition 2
Eqn2='a./(1+exp(-(-0.0002278*b.^3+0.008699*b.^2-0.1226*b+0.7494).*(x-c))) + d*exp(-(1.656/(2.355*b))^2*0.25*(x-c)^2)/(1+(1.656/(2.355*b))^2*((x-c)^2)) + k';
% Condition 4
Eqn4='a./(1+exp(-(-0.0002278*b.^3+0.008699*b.^2-0.1226*b+0.7494).*(x-c))) + d*exp(-(1.656/(2.355*b))^2*0.25*(x-c)^2)/(1+(1.656/(2.355*b))^2*((x-c)^2)) + e.*exp(-(((x-f).^2)./(2.*(b).^2).^g)) +k';

alfa=0.95;

for p=1:noCells 
    vector2=[];
    vector4=[];
    data = APO.GFP(p).avgImageData;
    dataList = nan(size(data,3), size(data,1)*size(data,2));
    BG=APO.GFP(p).avgImageData(:,:,70:end);
    
    for frame = 1:size(data,3)
        frameData = data(:,:,frame);
        dataList(frame, :) = frameData(:);
    end
    zData(:,1) = 1:size(data,3);
    
    zValueR=round(APO.B1AR(p).results.zValueAvg(:)./30);
    zValue=APO.B1AR(p).results.zValueAvg(:)./30;
    sigma=APO.B1AR(p).outputs.sig(:);
    Bmax=nan(size(zValue));
    
    rsq=nan(size(zValue));
    aCoef=nan(size(zValue));
    bCoef=nan(size(zValue));
    cCoef=nan(size(zValue));
    dCoef=nan(size(zValue));
    eCoef=nan(size(zValue));
    fCoef=nan(size(zValue));
    gCoef=nan(size(zValue));
    kCoef=nan(size(zValue));
    
    aCoefErr=nan(size(zValue));
    bCoefErr=nan(size(zValue));
    cCoefErr=nan(size(zValue));
    dCoefErr=nan(size(zValue));
    eCoefErr=nan(size(zValue));
    fCoefErr=nan(size(zValue));
    gCoefErr=nan(size(zValue));
    kCoefErr=nan(size(zValue));
    
    vector2=find((APO.Fit(p).AdjfitTypeAll==2));
    vector4=find((APO.Fit(p).AdjfitTypeAll==4));
    
    parfor i=1:size(dataList,2)      
        if sum(i==vector2)==0 || isnan(zValue(i))
            continue
        else
            
            Bmax(i)=mean(dataList(85:end,i));
            
            sP1=Icyt(p);
            sP2=sigma(i);
            sP3=zValue(i);
            sP4=Icyt(p)-Bmax(i);
            sP5=mean(dataList(1:10,i));
            
            startPoints = [sP1 sP2 sP3 sP4 sP5];


             LowerBound = [Icyt(p)-0.1*Icyt(p), sP2, zValue(i), 0 ,sP5];
             UpperBound = [Icyt(p)+0.1*Icyt(p), sP2, zValue(i), max(dataList(:,i)),sP5];
             
            e=dataList(:,i);
            [f, gof] = fit(zData,e,Eqn2,'Start',startPoints,'Lower',LowerBound,'Upper',UpperBound);
            
            rsq(i)=gof.adjrsquare;
            aCoef(i)=f.a;
            bCoef(i)=f.b;
            cCoef(i)=f.c;
            dCoef(i)=f.d;
            kCoef(i)=f.k;
            
            ci=confint(f, alfa);
            t = tinv((1+alfa)/2, gof.dfe);
            se = (ci(2,:)-ci(1,:)) ./ (2*t);
            
            aCoefErr(i)=se(1);
            bCoefErr(i)=se(2);
            cCoefErr(i)=se(3);
            dCoefErr(i)=se(4);
            kCoefErr(i)=se(5);
     
        end
    end
    
    parfor i=1:size(dataList,2) 
        if sum(i==vector4)==0 || isnan(zValue(i)) || zValueR(i)>89 || isnan(zValueR(i))
            continue
        else
            
            Bmax(i)=mean(dataList(85:end,i));
            sP1=Icyt(p); 
            sP2=sigma(i);
            sP3=zValue(i);
            sP4=max(dataList(:,i))/2;
            sP5=-(max(dataList(:,i))-mean(dataList(80:end,i)));
            sP6=zData(end);
            sP7=1;
            sP8=mean(dataList(1:10,i));
            
           
            startPoints = [sP1 sP2 sP3 sP4 sP5 sP6 sP7 sP8];
            
            LowerBound = [Icyt(p)-0.1*Icyt(p), sP2, zValue(i), 0, -2*max(dataList(:,i)), zValueR(i)+20, 0.8, sP8];
            UpperBound = [Icyt(p)+0.1*Icyt(p), sP2, zValue(i), sP4*2, 0, 150, 1.2, sP8];
            %weights=[0.5*ones(1,zValueR(i)-10),  ones(1, 20),  0.5*ones(1,100-zValueR(i)-10)];
            
            e=dataList(:,i);
            [f, gof] = fit(zData,e,Eqn4,'Start',startPoints,'Lower',LowerBound,'Upper',UpperBound)%,'Weights',weights);
            
            rsq(i)=gof.adjrsquare;
            aCoef(i)=f.a;
            bCoef(i)=f.b;
            cCoef(i)=f.c;
            dCoef(i)=f.d;
            eCoef(i)=f.e;
            fCoef(i)=f.f;
            gCoef(i)=f.g;
            kCoef(i)=f.k;
            
            ci=confint(f, alfa);
            t = tinv((1+alfa)/2, gof.dfe);
            se = (ci(2,:)-ci(1,:)) ./ (2*t);
            
            aCoefErr(i)=se(1);
            bCoefErr(i)=se(2);
            cCoefErr(i)=se(3);
            dCoefErr(i)=se(4);
            eCoefErr(i)=se(5);
            fCoefErr(i)=se(6);
            gCoefErr(i)=se(7);
            kCoefErr(i)=se(8);
            
        end
    end
    
    APO.GFPFit(p).aCoef=reshape(aCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).bCoef=reshape(bCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).cCoef=reshape(cCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).dCoef=reshape(dCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).eCoef=reshape(eCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).fCoef=reshape(fCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).gCoef=reshape(gCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).kCoef=reshape(kCoef, [size(data,1), size(data,2)]);
    APO.GFPFit(p).rsq=reshape(rsq, [size(data,1), size(data,2)]);
    
    APO.GFPFit(p).aCoefErr=reshape(aCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).bCoefErr=reshape(bCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).cCoefErr=reshape(cCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).dCoefErr=reshape(dCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).eCoefErr=reshape(eCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).hCoefErr=reshape(fCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).gCoefErr=reshape(gCoefErr, [size(data,1), size(data,2)]);
    APO.GFPFit(p).kCoefErr=reshape(kCoefErr, [size(data,1), size(data,2)]);
    
    disp(p); 

end


%% Filtering on intensity data
for p=1:noCells
    
APO.results(p).filtdCoef=APO.GFPFit(p).dCoef;
APO.results(p).filtdCoef(APO.GFPFit(p).rsq<0.9)=NaN; % Filtering of R-squared 
APO.results(p).filtdCoef((APO.GFPFit(p).dCoefErr./APO.results(p).filtdCoef)>0.3)=NaN; % Filtering on relative error

    r=3; % 3x3 block of moving pixels
    f = @(A)median(A(~isnan(A)));
    f_error = @(A)std(A(~isnan(A)));
    f_sum = @(A)sum(A(~isnan(A)));
    
    APO.results(p).filteredRes = nlfilter(APO.results(p).filtdCoef, [r r], f);
    APO.results(p).filteredResError = nlfilter(APO.results(p).filtdCoef, [r r], f_error);
    
    min_pixel_mask=nlfilter(~isnan(APO.results(p).filtdCoef), [r r], f_sum); % mask that shows how many pixels in vicinity are resolved
    min_pixel=3; % minimal number of pixels needed is (3+1)=4
    min_pixel_mask=min_pixel_mask>min_pixel; % larger than 3, not equal to and larger than 3
    
    APO.results(p).filteredRes=APO.results(p).filteredRes.*min_pixel_mask;
    APO.results(p).filteredRes(APO.results(p).filteredRes==0)=NaN;
    
    APO.results(p).filteredResError=APO.results(p).filteredResError.* ~isnan(APO.results(p).filteredRes);
    APO.results(p).filteredResError(APO.results(p).filteredResError==0)=NaN;
    
    APO.results(p).maxInt=APO.results(p).filteredRes;
    APO.results(p).maxIntErr=APO.results(p).filteredResError;


end

%% Save data

    varFit=APO.GFPFit;
    varFitType=APO.Fit;
    B1ARdata=APO.B1AR;
    results=APO.results;
    save('Fit.mat','varFit');
    save('FitType.mat','varFitType');
    save('Icyt.mat','Icyt');
    save('nom.mat','nom');
    save('B1ARdata.mat','B1ARdata');
    save('results.mat','results');
%% Plot data

Cell=APO.B1AR(p);
[B1AR,maxInt] = DensityPlotting(Cell,dz);
Cell=[];
a=APO.results(p).maxInt;
ratio=a./maxInt;
ExpressionLevelFactor(p)=Icyt(p)./nanmedian(APO.B1AR(p).outputs.maxInt(:));
ratio=ratio./ExpressionLevelFactor(p);
ratio(isnan(ratio) & ~isnan(maxInt))=0;
ratio=100.*ratio./1.07;

figure
surface(ratio,'EdgeColor','k')
axis image
colormap([parula])
caxis([15 40])
cbh = colorbar ; %Create Colorbar
cbh.Ticks = [15 30 45]; %Create 8 ticks from zero to 1
cbh.TickLabels = {'15','30','45'} ;
cbh.Label.String = 'Activation Prob. wrt ISO [%]';
