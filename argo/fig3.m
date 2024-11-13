%% Figure 3 -  POC anomalies during MHWs
% modified 11/12/2024  by Mari Bif, marianabif@gmail.com

% To run this code, make sure bbp_mhw_processed.mat is in the same directory of this .m file.

% This figure uses a dataset that was extracted from the Argo Global Data Assembly Center (USGODAE; usgodae.org/ftp/outgoing/argo/) and processed for small and large particles. 
% Only quality-controlled files (Sprof) and data flagged as “good” were considered. 
% We focused the analysis in the upper ocean between the surface and 400 m, with data every 5 to 10 days and acquired between June 2010 and October 2022, 
% with a gap between February 2016 and August 2018. The study region was constrained between 47˚-52.5˚N and 137.5-146˚W. For more details on how to estimate 
% small and large particles, chlorophyll concentrations and mixed layer depths please see: Bif et al. and references wherein.

% NOTE: The user may need to adjust the fig position for proper screen
%       fig.Position = fig.Position + [0 0 0 0]; ~ line 244

% To reproduce the color scheme of these figures, please download
% ColorBrewer colorschemes:
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps

clc
clear all
close all

load("bbp_mhw_processed.mat"); % Processed bbp signals, separating small and large particles
%% Converting bbp into POC (Johnson et al., 2017), matrixes
poc_small = ((31200*bbp_small)+3.04);
poc_big = ((31200*bbp_big)+3.04);

datematrix=datevec(datet); month = datematrix(:,2); month=month'; year= datematrix(:,1);year=year';day= datematrix(:,3);day=day';

yr1 = [2011,2012,2013];
yr2 = [2018,2021,2022];
%% mean POC of non-warm years, PB1 (pre-blob1) and PB2 (pre-blob2)

poc_small_med_pb1 = zeros(size(poc_small));
poc_small_med_pb2 = zeros(size(poc_small));
mld_mean1 = zeros(size(mld));mld_mean2 = zeros(size(mld));
month_mean1 = zeros(size(month));month_mean2 = zeros(size(month));
press_pb1 = zeros(size(press));press_pb2 = zeros(size(press));

for j=1:size(poc_small,2) 
    for uniqueyear = yr1
            for uniquemonth = 1:12
            [~,n]= find(ismember(year,yr1)); 
            [~,o]= find(ismember(month,uniquemonth));
            [val,~]=intersect(n,o);        
                mld_mean1(1,val(1)) =mean(mld(1,val));
                month_mean1(1,val(1))= uniquemonth;
                poc_small_med_pb1(:,val(1)) =mean(poc_small(:,val),2,"omitmissing");    
                press_pb1(:,val(1)) =mean(press(:,val),2,"omitmissing");
             
            end
    end
end

% Removing zero columns
month_mean1(:,all(~month_mean1,1))=[]; mld_mean1(:,all(~mld_mean1,1))=[];
poc_small_med_pb1(:,~all(poc_small_med_pb1,1))=[];press_pb1(:,~all(press_pb1,1))=[];



% Pre-blob2
for j=1:size(poc_small,2) 
    for uniqueyear = yr2
            for uniquemonth = 1:12
            [~,n]= find(ismember(year,yr2)); 
            [~,o]= find(ismember(month,uniquemonth));
            [val,~]=intersect(n,o);        
                mld_mean2(1,val(1)) =mean(mld(1,val));
                month_mean2(1,val(1))= uniquemonth;
                poc_small_med_pb2(:,val(1)) =mean(poc_small(:,val),2,"omitmissing");    
                press_pb2(:,val(1)) =mean(press(:,val),2,"omitmissing");
             
            end
    end
end

% Removing zero columns
month_mean2(:,all(~month_mean2,1))=[]; mld_mean2(:,all(~mld_mean2,1))=[];
poc_small_med_pb2(:,~all(poc_small_med_pb2,1))=[];press_pb2(:,~all(press_pb2,1))=[];

%% Smoothing data every 20m during non-blob years, pb1, bp2

% non-blob1
for j=1:size(press_pb1,2)
    i = 0:20:max(press_pb1(:,j));  
    for m=1:size(i,2)
         int2 = find(press_pb1(:,j)>=i(1,m),1,'first');
         int3 = find(press_pb1(:,j)>=i(1,m)+20,1,'first');

        poc_smooth_pb1(m,j) =mean(poc_small_med_pb1(int2:int3,j));    
        press_smooth_pb1(m,j) =i(1,m);

    end
end     
poc_smooth_pb1(poc_smooth_pb1==0) = NaN;press_smooth_pb1(press_smooth_pb1==0) = NaN;press_smooth_pb1(1,:)=0;
poc_smooth_pb1=poc_smooth_pb1(1:50,:);press_smooth_pb1=press_smooth_pb1(1:50,:); % nothing matter below this point. some profiles >900m dont have data

% non-blob2 
for j=1:size(press_pb2,2)
    i = 0:20:max(press_pb2(:,j));  
    for m=1:size(i,2)
         int2 = find(press_pb2(:,j)>=i(1,m),1,'first');
         int3 = find(press_pb2(:,j)>=i(1,m)+20,1,'first');

        poc_smooth_pb2(m,j) =mean(poc_small_med_pb2(int2:int3,j));    
        press_smooth_pb2(m,j) =i(1,m);

    end
end
poc_smooth_pb2(poc_smooth_pb2==0) = NaN;press_smooth_pb2(press_smooth_pb2==0) = NaN;press_smooth_pb2(1,:)=0;
poc_smooth_pb2=poc_smooth_pb2(1:50,:);press_smooth_pb2=press_smooth_pb2(1:50,:);


%% Selecting and smoothing data every 20m for 2014, 2015, 2019, 2020

%2014, 2015 blob1

poc_small_med_2014 = zeros(size(poc_small));
poc_small_med_2015 = zeros(size(poc_small));
month_2014 = zeros(size(month));month_2015 = zeros(size(month));
datet_2014 = zeros(size(datet));datet_2015 = zeros(size(datet));
press_2014 = zeros(size(press));press_2015 = zeros(size(press));

for j=1:size(poc_small,2) 
    for uniqueyear = 2014
            [~,n]= find(ismember(year,2014));       
                month_2014 = month(1,n);
                datet_2014 = datet(1,n);
                poc_small_med_2014 = poc_small(:,n);    
                press_2014 = press(:,n);
    end
end
for j=1:size(poc_small_med_2014,2) 
    i = 0:20:max(press_2014(:,j));  
    for      m=1:size(i,2)
         int2 = find(press_2014(:,j)>=i(1,m),1,'first');
         int3 = find(press_2014(:,j)>=i(1,m)+20,1,'first');

        poc_smooth_2014(m,j) =mean(poc_small_med_2014(int2:int3,j));    
        press_smooth_2014(m,j) =i(1,m);
    end
end
poc_smooth_2014 = poc_smooth_2014(1:50,:); press_smooth_2014 = press_smooth_2014(1:50,:);


for j=1:size(poc_small,2) 
    for uniqueyear = 2015
            [~,n]= find(ismember(year,2015));       
                month_2015 = month(1,n);
                datet_2015 = datet(1,n);
                poc_small_med_2015 = poc_small(:,n);    
                press_2015 = press(:,n);
    end
end
for j=1:size(poc_small_med_2015,2) 
    i = 0:20:max(press_2015(:,j));  
    for      m=1:size(i,2)
         int2 = find(press_2015(:,j)>=i(1,m),1,'first');
         int3 = find(press_2015(:,j)>=i(1,m)+20,1,'first');

        poc_smooth_2015(m,j) =mean(poc_small_med_2015(int2:int3,j));    
        press_smooth_2015(m,j) =i(1,m);
    end
end
poc_smooth_2015 = poc_smooth_2015(1:50,:); press_smooth_2015 = press_smooth_2015(1:50,:);





%2019, 2020 blob2
poc_small_med_2019 = zeros(size(poc_small));
poc_small_med_2020 = zeros(size(poc_small));
month_2019 = zeros(size(month));month_2020 = zeros(size(month));
datet_2019 = zeros(size(datet));datet_2020 = zeros(size(datet));
press_2019 = zeros(size(press));press_2020 = zeros(size(press));

for j=1:size(poc_small,2) 
    for uniqueyear = 2019
            [~,n]= find(ismember(year,2019));       
                month_2019 = month(1,n);
                datet_2019 = datet(1,n);
                poc_small_med_2019 = poc_small(:,n);    
                press_2019 = press(:,n);
    end
end
for j=1:size(poc_small_med_2019,2) 
    i = 0:20:max(press_2019(:,j));  
    for      m=1:size(i,2)
         int2 = find(press_2019(:,j)>=i(1,m),1,'first');
         int3 = find(press_2019(:,j)>=i(1,m)+20,1,'first');

        poc_smooth_2019(m,j) =mean(poc_small_med_2019(int2:int3,j));    
        press_smooth_2019(m,j) =i(1,m);
    end
end
poc_smooth_2019 = poc_smooth_2019(1:50,:); press_smooth_2019 = press_smooth_2019(1:50,:);


for j=1:size(poc_small,2) 
    for uniqueyear = 2020
            [~,n]= find(ismember(year,2020));       
                month_2020 = month(1,n);
                datet_2020 = datet(1,n);
                poc_small_med_2020 = poc_small(:,n);    
                press_2020 = press(:,n);
    end
end
for j=1:size(poc_small_med_2020,2) 
    i = 0:20:max(press_2020(:,j));  
    for      m=1:size(i,2)
         int2 = find(press_2020(:,j)>=i(1,m),1,'first');
         int3 = find(press_2020(:,j)>=i(1,m)+20,1,'first');

        poc_smooth_2020(m,j) =mean(poc_small_med_2020(int2:int3,j));    
        press_smooth_2020(m,j) =i(1,m);
    end
end
poc_smooth_2020 = poc_smooth_2020(1:50,:); press_smooth_2020 = press_smooth_2020(1:50,:);

%% Anomaly sections, pb1 - non-blob and pb2 - non-blob

for j=1:size(month_mean1,2)

    [~,m]= find(ismember(month_2014,month_mean1(1,j)));
    [~,n]= find(ismember(month_2015,month_mean1(1,j)));
    anom_2014(:,m) = bsxfun(@minus,poc_smooth_2014(:,m), poc_smooth_pb1(:,j));
    anom_2015(:,n) = bsxfun(@minus,poc_smooth_2015(:,n), poc_smooth_pb1(:,j));

end


for j=1:size(month_mean2,2)

    [~,m]= find(ismember(month_2019,month_mean2(1,j)));
    [~,n]= find(ismember(month_2020,month_mean2(1,j)));
    anom_2019(:,m) = bsxfun(@minus,poc_smooth_2019(:,m), poc_smooth_pb2(:,j));
    anom_2020(:,n) = bsxfun(@minus,poc_smooth_2020(:,n), poc_smooth_pb2(:,j));

end

% Figures
fig = figure;
%fig.Position = fig.Position + [0 -400 0 0];
fig.Position(4) = fig.Position(4) * 1.5;  

[m,~] = size(anom_2014);
                date_tseries = repmat(datet_2014,[m,1]);
datetint = date_tseries; bbpint =  anom_2014; pressint = press_smooth_2014;
pressint(isnan(datetint))=[];bbpint(isnan(datetint))=[];datetint(isnan(datetint))=[];
datetint(isnan(bbpint))=[];pressint(isnan(bbpint))=[];bbpint(isnan(bbpint))=[];
bbpint(isnan(pressint))=[];datetint(isnan(pressint))=[];pressint(isnan(pressint))=[];
[O,P] = ndgrid(min(datetint(:)):1:max(datetint(:)), 0:5:1000);
V = griddata(datetint(:),pressint(:),bbpint(:),O,P);
datetint=datetint';pressint=pressint';bbpint=bbpint';
F = scatteredInterpolant(datetint,pressint,bbpint); %ndgrid needs to be employed
Vq = F(O,P);

t = tiledlayout(4,1);
nexttile(1);
hold on
set(gca,'YDir','reverse')
contourf(O,P,Vq,30,'LineStyle','none')
box off; grid off; datetick('x','mmm','keeplimits');
title('POC anomalies - 2014');
set(gca,'FontSize',14,'LineWidth',2); cbh=colorbar; ylabel(cbh,'POC (mg m^-^3)'); colormap (brewermap(16,'PRGn')); 
set(gca,'tickdir','out');  set(cbh,'tickdir','out'); set(cbh,'Ytick',[-15, -10,-5,0,5, 10, 15]); caxis([-15 15]);
ylim([0 400]); ylabel('Depth (m)');xtickangle(45);
xlim([min(datet_2014) max(datet_2014)]);
xL=xlim;
text(xL(2),-50,'a','FontSize', 14)
hold off

[m,~] = size(anom_2015);
                date_tseries = repmat(datet_2015,[m,1]);
datetint = date_tseries; bbpint =  anom_2015; pressint = press_smooth_2015;
pressint(isnan(datetint))=[];bbpint(isnan(datetint))=[];datetint(isnan(datetint))=[];
datetint(isnan(bbpint))=[];pressint(isnan(bbpint))=[];bbpint(isnan(bbpint))=[];
bbpint(isnan(pressint))=[];datetint(isnan(pressint))=[];pressint(isnan(pressint))=[];
[O,P] = ndgrid(min(datetint(:)):1:max(datetint(:)), 0:5:1000);
V = griddata(datetint(:),pressint(:),bbpint(:),O,P);
datetint=datetint';pressint=pressint';bbpint=bbpint';
F = scatteredInterpolant(datetint,pressint,bbpint); %ndgrid needs to be employed
Vq = F(O,P);

nexttile(2);
hold on
set(gca,'YDir','reverse')
contourf(O,P,Vq,30,'LineStyle','none')
box off; grid off; datetick('x','mmm','keeplimits');
title('2015');
set(gca,'FontSize',14,'LineWidth',2); cbh=colorbar; ylabel(cbh,'POC (mg m^-^3)'); colormap (brewermap(16,'PRGn')); 
set(gca,'tickdir','out');  set(cbh,'tickdir','out'); set(cbh,'Ytick',[-15, -10,-5,0,5, 10, 15]); caxis([-15 15]);
ylim([0 400]); ylabel('Depth (m)');xtickangle(45);
xlim([min(datet_2015) max(datet_2015)]);
xL=xlim;
text(xL(2),-50,'b','FontSize', 14)
hold off




[m,~] = size(anom_2019);
                date_tseries = repmat(datet_2019,[m,1]);
datetint = date_tseries; bbpint =  anom_2019; pressint = press_smooth_2019;
pressint(isnan(datetint))=[];bbpint(isnan(datetint))=[];datetint(isnan(datetint))=[];
datetint(isnan(bbpint))=[];pressint(isnan(bbpint))=[];bbpint(isnan(bbpint))=[];
bbpint(isnan(pressint))=[];datetint(isnan(pressint))=[];pressint(isnan(pressint))=[];
[O,P] = ndgrid(min(datetint(:)):1:max(datetint(:)), 0:5:1000);
V = griddata(datetint(:),pressint(:),bbpint(:),O,P);
datetint=datetint';pressint=pressint';bbpint=bbpint';
F = scatteredInterpolant(datetint,pressint,bbpint); %ndgrid needs to be employed
Vq = F(O,P);

nexttile(3);
hold on
set(gca,'YDir','reverse')
contourf(O,P,Vq,30,'LineStyle','none')
box off; grid off; datetick('x','mmm','keeplimits');
title('2019');
set(gca,'FontSize',14,'LineWidth',2); cbh=colorbar; ylabel(cbh,'POC (mg m^-^3)'); colormap (brewermap(16,'PRGn')); 
set(gca,'tickdir','out');  set(cbh,'tickdir','out'); set(cbh,'Ytick',[-15, -10,-5,0,5, 10, 15]); caxis([-15 15]);
ylim([0 400]); ylabel('Depth (m)');xtickangle(45);
xlim([min(datet_2019) max(datet_2019)]);
xL=xlim;
text(xL(2),-50,'c','FontSize', 14)
hold off


[m,~] = size(anom_2020);
                date_tseries = repmat(datet_2020,[m,1]);
datetint = date_tseries; bbpint =  anom_2020; pressint = press_smooth_2020;
pressint(isnan(datetint))=[];bbpint(isnan(datetint))=[];datetint(isnan(datetint))=[];
datetint(isnan(bbpint))=[];pressint(isnan(bbpint))=[];bbpint(isnan(bbpint))=[];
bbpint(isnan(pressint))=[];datetint(isnan(pressint))=[];pressint(isnan(pressint))=[];
[O,P] = ndgrid(min(datetint(:)):1:max(datetint(:)), 0:5:1000);
V = griddata(datetint(:),pressint(:),bbpint(:),O,P);
datetint=datetint';pressint=pressint';bbpint=bbpint';
F = scatteredInterpolant(datetint,pressint,bbpint); %ndgrid needs to be employed
Vq = F(O,P);

nexttile(4);
hold on
set(gca,'YDir','reverse')
contourf(O,P,Vq,30,'LineStyle','none')
box off; grid off; datetick('x','mmm','keeplimits');
title('2020');
set(gca,'FontSize',14,'LineWidth',2); cbh=colorbar; ylabel(cbh,'POC (mg m^-^3)'); colormap (brewermap(16,'PRGn')); 
set(gca,'tickdir','out');  set(cbh,'tickdir','out'); set(cbh,'Ytick',[-15, -10,-5,0,5, 10, 15]); clim([-15 15]);
ylim([0 400]); ylabel('Depth (m)');xtickangle(45);
xlim([min(datet_2020) max(datet_2020)]);
xL=xlim;
text(xL(2),-50,'d','FontSize', 14)
xL=xlim;
hold off

%saveas(gcf,'Fig3','eps'); saveas(gcf,'Fig3','tiff'); saveas(gcf,'Fig3','svg'); saveas(gcf,'Fig3.fig'); 