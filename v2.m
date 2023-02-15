clear all; close all; clc; % clear workspace
c = readNPY('cenogrid_data.npy'); % load data

%% 1. show GEV works well

b = 20; % block size - factors of 24321 (33 & 121 yield same results)
T = -c(:,3); % temperature fluctuations are proportional to delta-18-O fluctuations
t = c(:,1);
%Trl = smoothdata(T(t>-2.6),'rloess',10); % filter quaternary cycles
%T(t>-2.6) = T(t>-2.6)-Trl+median(Trl);
T(c(:,1)>-2) = [];
for i = 1:(length(T)./b);
    Tb = T((1+b.*(i-1)):b.*i); % block
    Tm(i) = max(Tb-nanmean(Tb))./nanstd(Tb); % standardized block maximum
end

gevf = mle(Tm,'distribution','gev'); % fit GEV
q = cdf('gev',unique(sort(Tm)),gevf(1),gevf(2),gevf(3)); % define cdf
figure; [y,x] = ecdf(Tm); % define ecdf
%KS = max(abs(y(2:end)-q')) % kolmogorov-smirnov statistic
KS = max(y(2:end)-q')+max(q'-y(2:end)) % kolmogorov-smirnov statistic

plot(x,y,'k','linewidth',2); % plotting
hold on;
plot(unique(sort(Tm)),q,'color',(1/200)*[110 23 78],'linewidth',2);
box on;
set(gca,'ticklabelinterpreter','latex','fontsize',16,'xtick',1:4,'ytick',0:.2:1);
axis([1 4 0 1])
ylabel('CDF','interpreter','latex');
xlabel('$\delta^{18}$O minima [$z$-score]','interpreter','latex')

lgnd = legend('Data','GEV');
set(lgnd,'interpreter','latex','fontsize',16)

axes('Position',[.475 .2 .4 .45])
[y,x] = ksdensity(Tm,linspace(min(Tm),max(Tm),sqrt(length(Tm))));
scatter(x,y,'k','filled');
hold on;
qp = pdf('gev',unique(sort(Tm)),gevf(1),gevf(2),gevf(3));
plot(unique(sort(Tm)),qp,'color',(1/200)*[110 23 78],'linewidth',2);
set(gca,'ticklabelinterpreter','latex','fontsize',16,'xtick',1:4,'ytick',0:.2:1);
axis([1 4 0 1.05])
ylabel('PDF','interpreter','latex');
box on;

%% show trend with delta-18-O

nboot = 100; % bootstrap iterations for uncertainty

t = c(:,1);
%Trl = smoothdata(T(t>-2.6),'rloess',10); % filter quaternary cycles
%T(t>-2.6) = T(t>-2.6)-Trl+median(Trl);
T(t>-2) = NaN; % remove quaternary period (trend the same when including)
%d13c = c(:,2);

b = 20; % minimum block size -- taken from fischer et al
for i = 1:(length(T)./b);
    Tb = T((1+b.*(i-1)):b.*i); % block
    tb = t((1+b.*(i-1)):b.*i); % block
    %d13cb = d13c((1+b.*(i-1)):b.*i);
    Tbin(i) = nanmean(Tb);
    tbin(i) = nanmean(tb);
    Tm(i) = max((Tb-nanmean(Tb))./nanstd(Tb)); % standardized block maximum
    %Tm(i) = max((d13cb-nanmean(d13cb))./nanstd(d13cb));
end

Ts_i = -4:0;
for i = 1:length(Tbin); % sort blocks by delta-18-O
    [~,indx(i)] = min(abs(Tbin(i)-Ts_i));
end

for i = 1:length(Ts_i); % fit GEV to each metablock
    qq = Tm(indx==i);
    qq = qq(~isnan(qq));
    gevft(:,i) = mle(qq,'distribution','gev');
    for j = 1:nboot; % bootstrap for uncertainty
        boot = randi(length(qq),1,length(qq));
        qqboot = qq(boot);
        gevftboot(:,i,j) = mle(qqboot,'distribution','gev');
    end
    i
end

for i = 1:length(Ts_i);
    qqx(i) = nanmean(Tbin(indx==i));
end

gevft = median(gevftboot,3);
%MU = squeeze(gevftboot(2,:,:));
%SIG = squeeze(gevftboot(3,:,:));
gevftboot = squeeze(gevftboot(1,:,:));
xi_u = mad(gevftboot');

xi = gevft(1,:);
d18O = -Ts_i;
wghts = 1./xi_u.^2;

figure;
subplot(131)
for i = 1:length(d18O);
    plot(linspace(d18O(i),d18O(i)),linspace(xi(i)-xi_u(i),xi(i)+xi_u(i)),'k','linewidth',5)
    hold on;
    scatter(d18O(i),xi(i),400,(1/200)*[110 23 78],'filled');
end

box on;
set(gca,'ticklabelinterpreter','latex','fontsize',16);%
set(gca,'xtick',0:4,'ytick',-.5:.1:0.1);
axis([-.25 4.25 -.5 .1])
ylabel('shape parameter $\xi$','interpreter','latex');
xlabel('$\delta^{18}$O [$^\circ/_{\circ\circ}$]','interpreter','latex')
axis square
text(0.25,-0.45,'A','interpreter','latex','fontsize',16)
%}

%%
%
subplot(132)
q1 = pdf('gev',0.5:.01:4,gevft(1,5),gevft(2,5),gevft(3,5));
q2 = pdf('gev',0.5:.01:4,gevft(1,1),gevft(2,1),gevft(3,1));
p1 = plot(0.5:.01:4,q1,'color','k','linewidth',2);
hold on;
yyaxis right;
qc1 = cdf('gev',0.5:.01:4,gevft(1,5),gevft(2,5),gevft(3,5));
qc2 = cdf('gev',0.5:.01:4,gevft(1,1),gevft(2,1),gevft(3,1));
rl = (1-qc1)./(1-qc2);
p3 = plot(0.5:.01:4,rl,'linewidth',2);
ylim([0 5])
set(gca,'ticklabelinterpreter','latex','fontsize',16,'ytick',0:5)
ylabel('$>z$ event relative probability','interpreter','latex')
yyaxis left;
ylabel('PDF','interpreter','latex')
p2 = plot(0.5:.01:4,q2,'--k','linewidth',2);
ylim([0 1.1])
set(gca,'ytick',0:.5:1,'xtick',1:5)
xlabel('$\delta^{18}$O minima [$z$-score]','interpreter','latex')
lgnd = legend([p1 p2 p3],'$\delta^{18}$O $\sim$ 0','$\delta^{18}$O $\sim$ 4','CCDF ratio');
set(lgnd,'interpreter','latex','fontsize',16)
xlim([0.5 4])
axis square
text(1,1,'B','interpreter','latex','fontsize',16)


%%

clear all; clc; % clear workspace
c = readNPY('cenogrid_data.npy'); % load data
T = -c(:,3); % temperature fluctuations are proportional to delta-18-O fluctuations
nboot = 100; % bootstrap iterations for uncertainty
t = c(:,1);
%Trl = smoothdata(T(t>-2.6),'rloess',10); % filter quaternary cycles
%T(t>-2.6) = T(t>-2.6)-Trl+median(Trl);
T(t>-2) = NaN;

b = 20; % minimum block size -- taken from fischer et al
for i = 1:(length(T)./b);
    Tb = T((1+b.*(i-1)):b.*i); % block
    tb = t((1+b.*(i-1)):b.*i); % block
    %d13cb = d13c((1+b.*(i-1)):b.*i);
    Tbin(i) = nanmean(Tb);
    tbin(i) = nanmean(tb);
    Tm(i) = max((Tb-nanmean(Tb))./nanstd(Tb)); % standardized block maximum
    %Tm(i) = max((d13cb-nanmean(d13cb))./nanstd(d13cb));
end

Ts_i = -4:.5:-1.5;
for i = 1:length(Tbin); % sort blocks by delta-18-O
    [~,indx(i)] = min(abs(Tbin(i)-Ts_i));
end

for i = 1:length(Ts_i); % fit GEV to each metablock
    qq = Tm(indx==i);
    qq = qq(~isnan(qq));
    gevft(:,i) = mle(qq,'distribution','gev');
    for j = 1:nboot; % bootstrap for uncertainty
        boot = randi(length(qq),1,length(qq));
        qqboot = qq(boot);
        gevftboot(:,i,j) = mle(qqboot,'distribution','gev');
    end
    i
end

for i = 1:length(Ts_i);
    qqx(i) = nanmean(Tbin(indx==i));
end

%gevft = median(gevftboot,3);
xi_u = squeeze(gevftboot(1,:,:));
xi_u = mad(xi_u');
mu_u = squeeze(gevftboot(2,:,:));
mu_u = mad(mu_u');
sig_u = squeeze(gevftboot(3,:,:));
sig_u = mad(sig_u');

xi = gevft(1,2:end-1);
mu = gevft(2,2:end-1);
sig = gevft(3,2:end-1);

xi_u = xi_u(2:end-1);
mu_u = mu_u(2:end-1);
sig_u = sig_u(2:end-1);

xi_w = 1./xi_u.^2;
mu_w = 1./mu_u.^2;
sig_w = 1./sig_u.^2;

dd18O = [0 -.5 -1 -1.5];
dT = -dd18O./.22;

%cftool
% xi = .01912(.00678)dT - .1745(.0304) from cftool
% mu -- no trend
% sig -- no trend

dTp = linspace(0,8,1000);
stddevs = 1:.001:4;

xim = .01912;
xib = -.1745;
for i = 1:length(dTp);
   qc1 = cdf('gev',stddevs,xim.*dTp(i)+xib,mu(1),sig(1));
   std3(i) = 1-qc1(stddevs==3);
end
std3 = std3./std3(1);

subplot(133)
plot(dTp,std3,'color',(1/200)*[110 23 78],'linewidth',3);
axis square
box on;
ylabel('$z>3$ event relative probability','interpreter','latex','fontsize',16)
set(gca,'ticklabelinterpreter','latex','fontsize',16,'TickLength',[0 0])
axis([0 8 1 7])
xlabel('Global temperature increase [$^\circ$C]','interpreter','latex')
hold on;
plot(linspace(2.2727,2.2727),linspace(6.8,7),'k')
plot(linspace(4.5455,4.5455),linspace(6.8,7),'k')
plot(linspace(6.8182,6.8182),linspace(6.8,7),'k')
plot(linspace(1.35,1.35),linspace(0,1.2),'k')
plot(linspace(2.7,2.7),linspace(0,1.2),'k')
plot(linspace(4.05,4.05),linspace(0,1.2),'k')
plot(linspace(5.4,5.4),linspace(0,1.2),'k')
plot(linspace(6.75,6.75),linspace(0,1.2),'k')
text(1.4,1.2,'1EgC','interpreter','latex','fontsize',14)
text(2.75,1.2,'2EgC','interpreter','latex','fontsize',14)
text(4.1,1.2,'3EgC','interpreter','latex','fontsize',14)
text(5.45,1.2,'4EgC','interpreter','latex','fontsize',14)
text(6.8,1.2,'5EgC','interpreter','latex','fontsize',14)
title('$\delta^{18}$O [$^\circ/_{\circ\circ}$]','interpreter','latex')
text(0.1,6.8,'3.5','interpreter','latex','fontsize',16)
text(2.3727,6.8,'3','interpreter','latex','fontsize',16)
text(4.6455,6.8,'2.5','interpreter','latex','fontsize',16)
text(6.9182,6.8,'2','interpreter','latex','fontsize',16)
text(1,5,'C','interpreter','latex','fontsize',16)
