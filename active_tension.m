csvf1 = readtable('407037echo/Tmax.csv');
csvf2 = readtable('407037mri/Tmax.csv');
swine = 407037;
if swine== 396669
    total_t = 0.66925; 
end
if swine== 396670
    total_t = 0.4845;
end
if swine== 399018
    total_t = 0.5383;
end
if swine== 399019
    total_t = 0.5383;
end
if swine== 407034
    total_t = 0.5383;
end
if swine== 407036
    total_t = 0.48445;
end
if swine== 407037
    total_t = 0.5921;
end
%Constant values
lr = 1.85;
lo = 1.28;
CaoMax = 4.35;
Cao = 4.35;
B = 4.75;
to = 275; %450];
ttrans = 300;
tau = 25;
BCL = 750; %Length of beats 60/80
Ttotal(:,1) = csvf1{:,2};
Ttotal(:,2) = csvf2{:,2};
ls_t=[];
ls_P=[];
for TT = 1:size(Ttotal,2)
    if TT==1
        T = Ttotal(:,1);
    else
        T = Ttotal(:,2);
    end
for j = 1:length(T)
    Tmax = T(j) *0.007519;
    timepoint =linspace(0, BCL, BCL/2);
    sigma_active = zeros(1, length(timepoint));
for i = 1:length(timepoint)
    ta = timepoint(i);
    lambda = 1.5;
    lso = ActiveLength(lambda, lr, lo);    
    deno = sqrt(exp((B*lso)-1));
    ECa50 = CaoMax/deno;    
    if ta<ttrans
         Ct = 0.5*(1-cos(pi*ta/to));
    else
         Ct = 0.5*(1-cos(pi*ttrans/to))*exp(-((ta-ttrans)/tau));
    CaTerm = Cao^2 /(Cao^2 + ECa50^2);
    Pact = Tmax*CaTerm*Ct;
    sigma_active(i) = Pact;
    end
end
ls_t(end+1)=j;
ls_P(end+1)=max(sigma_active);

end
t1=ls_t(1:30);
t2=ls_t(31:end);
p1=ls_P(1:30);
p2=ls_P(31:end);
end
figure(1)
hold on
plot(t1*total_t/30, p1,'LineWidth',1.5,'DisplayName','Active Tension ECHO')
hold on
plot(t2*total_t/30, p2,'LineWidth',1.5,'DisplayName','Active Tension MRI')
actmax(j) = max(sigma_active)/1000.0;


grid on
% legend({'k = 0','k = 0.14','k = 0.22','k = 0.265'}, 'FontSize',12)%,'17','18','19')%,'change','pp') %,'0.35,0.05','0.35,0.09','0.35,0.5')
legend({'Active Tension ECHO','Active Tension MRI'}, 'FontSize',24, 'Fontname','Times New Roman', 'Location','south')
%title('Isometric tension plot','FontSize',24)
xlabel('Time (s)','FontSize',24, 'fontWeight','bold', 'Fontname','Times New Roman')
ylabel('Active Tension (mmHg)','FontSize',24, 'fontWeight','bold', 'Fontname','Times New Roman')
hold off
saveas(gcf, 'activetension_', 'png')

function lso = ActiveLength(lamda, lr, lo)
ls = lamda*lr;
if ls<=lo
    lso = 0.002;
else 
    lso = ls-lo;

end
end