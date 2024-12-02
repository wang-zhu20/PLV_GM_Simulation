% Program for generating artificial ground motions by stochastic ground
% motion model using wavelet packets
% 10/31/2010
% coded by Yoshi Yamamoto (Stanford University)
% yama4423@stanford.edu

% load('E:\PulseResearch_NGA2\PredictPara.mat')
% PredictPara = PredictPara(10,:);
% M = PredictPara(:,4);R = PredictPara(:,6);Vs30 = PredictPara(:,7);
% SorD = PredictPara(:,8);


clear;
close all
clc;

addpath('./subroutines');
% Rrup:     Rupture distance (km)
% Rhyp:     Hypocentral distance (km)
% inVs30:   Average shear wave velocity with 30m in surface
% Mw:       Minimum Moment Magnitude
% nsmpl:    # of samples to generate
% iflg:     =1  1 realization                %%% with variability G.W.
%           =2  wave from median parameters  %%% no variability   G.W.
% offset:   # of offset for output filename

% data file that contains the result of regression analysis
filereg='.//regdata//regression_equation.txt';

iflg = 1;
offset=0;

% % conditions for each case
% ind=1;
% Mw(ind)=6.53;	Rrup(ind)= 0.07;	Rhyp(ind)= 0.07;	Vs30(ind)= 264.57; 
% SorD(ind) = 19.5;    nsmpl(ind)= 100;  	ind=ind+1;
% Mw(ind)=8;	Rrup(ind)= 10.0000;	Rhyp(ind)= 10.0000;	Vs30(ind)= 400.0;	
% nsmpl(ind)= 3;	ind=ind+1;
% 10^(-2.57+0.62*M)
% total of case 

ind=1;
Mw(ind)=6.5;	Rrup(ind)=20;	Rhyp(ind)= 20;	Vs30(ind)= 270; 
SorD(ind) = 26;    nsmpl(ind)= 100;  	ind=ind+1;


% Mw(ind)=6.5;	Rrup(ind)=10;	Rhyp(ind)= 10;	Vs30(ind)= 400; 
% SorD(ind) = 26;    nsmpl(ind)= 100;  	ind=ind+1;


Mw(ind)=7.0;	Rrup(ind)=20;	Rhyp(ind)= 20;	Vs30(ind)= 270; 
SorD(ind) = 53;    nsmpl(ind)= 100;  	ind=ind+1;


% Mw(ind)=7.0;	Rrup(ind)=10;	Rhyp(ind)= 10;	Vs30(ind)= 400; 
% SorD(ind) = 53;    nsmpl(ind)= 100;  	ind=ind+1;


Mw(ind)=7.5;	Rrup(ind)=20;	Rhyp(ind)= 20;	Vs30(ind)= 270; 
SorD(ind) = 108;    nsmpl(ind)= 100;  	ind=ind+1;

% Mw(ind)=7.5;	Rrup(ind)=10;	Rhyp(ind)= 10;	Vs30(ind)= 400; 
% SorD(ind) = 108;    nsmpl(ind)= 100;  	ind=ind+1;

ncase=ind-1;

% initializing random number generater
rs.TotalSign = RandStream('mt19937ar');
rs.MinorRand = RandStream('mt19937ar');
rs.MajorLoc  = RandStream('mt19937ar');
rs.MajorAmp  = RandStream('mt19937ar');
inrand=ceil(rand(1)*100);
for i=1:1:inrand
    rand(rs.TotalSign);
    rand(rs.MinorRand);
    rand(rs.MajorLoc);
    rand(rs.MajorAmp);
end

% loop for each case
for j=1:1:ncase
    % making directory for output files
    filename=sprintf('M%03.1f_Rr%08.4f_Rh%08.4f_Vs%06.1f',Mw(j),Rrup(j),Rhyp(j),Vs30(j));
    system(sprintf('mkdir %s',filename));

    % loop for each sample
    for i=1:1:nsmpl(j)
        disp(['computing... smpl' num2str(i,'% 6d') '/' num2str(nsmpl(j),'% 6d') '| case' num2str(j,'% 6d') '/' num2str(ncase,'% 6d')]);
        % generating each sample
%         [th dt] = fn_get1Sim(Mw(j),Rrup(j),Vs30(j),rs,iflg,filereg,Rhyp(j));

        [th dt] = fn_get1Sim_ver(Mw(j),Rrup(j),Vs30(j),rs,iflg,filereg,Rhyp(j));
        % storing each sample in NGA format
        
        [Pulse_record, dt,Vp,Tp,Et,Eacc] = sim_single(Mw(j),Rrup(j),Vs30(j),SorD(j));
        
        acc_raw = extract_major(th);

        acc_raw = acc_raw'; % 转为列
        High = 0.1; Low = 30;nrl = 2;
        [acc2,dt,HPused,LPused,nrolused,baselineOption] = Process(acc_raw,dt,High,Low,nrl);
         acc2 = acc2'; % 转为行
         
         
        acc_temp = extract_major(acc2);
        acc = [zeros(1,756),acc_temp]; % 补7.56 秒的0
        vel = cumsum(acc) .* dt .* 981;
        
        
           a = 0.05;
             Pulse_temp = taper(Pulse_record',a);
        Pulse_record = Pulse_temp';
        % 相同长度相加
%         [acc_total,vel_total] = add_pulse_res(vel,Pulse_record);
        sim1 = Pulse_record;dt = 0.01;
        
        %  直接调整位置
        time=[1:length(acc)]*dt;
        Ecum=cumsum(acc.^2);
        total=Ecum(length(Ecum));
        a=0.1+rand*0.4;
        t15=min(time(Ecum>=a*total));
        if Et+7.56<t15
            sim1=[zeros(1,round(100*(t15-Et-7.56))),sim1];
        else
            sim1 = sim1(round(100*(Et+7.56-t15))+1:end);
        end
        
        
        
        if length(sim1)<length(vel)
            sim1=[sim1,zeros(1,length(vel)-length(sim1))];
        else
            sim1 = sim1(1:length(vel));
        end
        vel_total = vel + sim1;
        acc_total = [0,diff(vel_total)./981./dt];
        
        
        
        acc_temp = acc_total';
        [acc_temp,~] = Process(acc_temp,dt,High,Low,nrl);
        acc_end = acc_temp';
        vel_end = cumsum(acc_end) .* dt .* 981;

        
% figure(999)
% subplot(4,1,1)
% plot(0.01*(1:length(vel)),vel,'color','b');set(gca,'FontSize',13);
% subplot(4,1,2)
% plot(0.01*(1:length(vel)),sim1,'color','b');set(gca,'FontSize',13);
% subplot(4,1,3)
% plot(0.01*(1:length(vel)),vel_total,'color','b');set(gca,'FontSize',13);
% subplot(4,1,4)
% plot(0.01*(1:length(vel)),vel_end,'color','b');set(gca,'FontSize',13);

% set(gca,'FontSize',13);
% set(gcf,'Units','centimeters','Position',[12 5 12 10])

% figure(1000)
% subplot(3,1,1)
% plot(0.01*(1:length(vel)),acc,'color','b');
% subplot(3,1,2)
% plot(0.01*(1:length(vel)),acc_total,'color','b');
% subplot(3,1,3)
% plot(0.01*(1:length(vel)),acc_end,'color','b');
figure(999)
t = tiledlayout(3,1) ;t.Padding = 'compact';t.TileSpacing = 'compact';

nexttile(1);
plot((1:length(vel))*dt,sim1,'k','Linewidth',0.5);hold on

set(gca,'xtick',[]);

nexttile(2);
plot((1:length(vel))*dt,vel,'k','Linewidth',0.5);hold on
set(gca,'xtick',[]);

nexttile(3);
plot((1:length(vel))*dt,vel_total,'k','Linewidth',0.5);hold on
ly = ylim;
xlabel( 'Time (s)');
ylabel(t,'Vel (cm/s)');
for iii = 1:3
    nexttile(iii);set(gca,'fontsize',12);xlim([0 length(vel_total)*dt]);
    %     xlim([0 50])
    if iii ~= 2
        ylim(ly);
    end
end
set(gcf,'Units','centimeters','Position',[12 8 12 9]); % 图片大小

% exportgraphics(gcf,[filename,'/Figure/figure',num2str(i),'.jpg'],'Resolution',300)

Simdata.acc = acc;
Simdata.vel = vel;
Simdata.pulse = sim1;
Simdata.vel_total = vel_total;
Simdata.acc_total = acc_total;
Simdata.vel_end = vel_end;
Simdata.acc_end = acc_end;

Simdata.dt = 0.01;
Simdata.pulsepara = [Vp,Tp,Et,Eacc];
DataName = [filename,'/Simdata',num2str(i),'.mat'];
save(DataName,'Simdata');

%         fn_writeACCinNGA(sprintf('.\\%s\\%s_#%07d.txt',filename,filename,(i+offset)),length(th),dt,th, Mw(j), Rrup(j), Rhyp(j), Vs30(j))
%     
%         AccFileName=sprintf('.\\%s\\%s_#%07d.txt',filename,filename,(i+offset));
%         % [T, PSA(i,:)]=PlotResponseSpectrum(AccFileName);   % added by Gang Wang
%         PlotResponseSpectrum(AccFileName);   % added by Gang Wang
    
    end
    

end

rmpath('./subroutines');
function [th] = extract_major(acc)

% time=[1:length(acc)]*dt;
Ecum=cumsum(acc.^2);
total=Ecum(length(Ecum));

% tRange=time((Ecum>=0.0001*total)&(Ecum<=0.9999*total));
% t1=min(tRange);
% t99=max(tRange);

th = acc((Ecum>=0.0001*total)&(Ecum<=0.9999*total));

end

% function [acc_total,vel_total] = add_pulse_res(vel,Pulse_record)
% 
%         sim1 = Pulse_record;dt = 0.01;
%         if length(sim1)<length(vel)
%             sim1=[sim1,zeros(1,length(vel)-length(sim1))];
%         else
%             sim1 = sim1(1:length(vel));
%         end
%         vel_total = vel + sim1;
%         acc_total = [0,diff(vel_total)./981./dt];
% end


%         acc = acc';
%         a = 0.05;LP = 30;HP = 0.1; nroll = 2;
%         acc2 = taper(acc,a);
%         [acc2,tpad] = zeropad(acc2,dt,nroll,HP);
%         acc2 = acausal(LP,HP,nroll,acc2,dt); %filter
     
%         acc2 = Removepad(acc2,tpad,dt);