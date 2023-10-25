% initial conditions
numSpecies = 11;
y0 = zeros(numSpecies,1);
y0(1,1) = 130000; %C8 initial 
y0(3,1) = 21000; %C3
y0(5,1) = 40000; %IAP
y0(7,1) = 40000; %BAR
y0(9,1) = 10000; %receptor
y0(11,1) = 6000; %ligand %100 %1

% y0(1,1) = findParams(1,i); %C8 initial
% y0(3,1) = findParams(2,i); %C3
% y0(5,1) = findParams(3,i); %IAP
% y0(7,1) = findParams(4,i); %BAR
% y0(9,1) = findParams(5,i); %receptor
% y0(11,1) = 100; %ligand

speciesNames = {'C8'; 'C8a'; 'C3'; 'C3a'; 'IAP'; 'C3a-IAP';...
    'BAR'; 'C8a-BAR'; 'receptor'; 'complex'; 'ligand'};
%%
numCells = 10000; %10000;
%%
% pd= makedist('Normal', 130000, 214500);
% 
% t= truncate(pd,13000, 1300000);
% 
% findParams = random(t);
%   
% t= makedist('Lognormal', 130000, 214500);
%%
for i=1:numCells %10000 
    
    %findParams(:,i) = [findParamValue(130000,10);findParamValue(21000,10);findParamValue(40000,10);findParamValue(40000,10);findParamValue(10000,10)]; %uniform distribution
    
    %findParams(:,i) = [random('Normal',130000,3); random('Normal',21000,3); random('Normal',40000,3); random('Normal',40000,3); random('Normal',10000,3)];
    %findParams(:,i) = [random('Normal',130000,214500); random('Normal',21000,34650); random('Normal',40000,66000); random('Normal',40000,66000); random('Normal',10000,16500)];
    findParams(:,i) = [random(truncate(makedist('Normal', 130000, 214500),13000, 1300000)); random(truncate(makedist('Normal', 21000, 34650),2100, 210000)); random(truncate(makedist('Normal', 40000, 66000),4000, 400000)); random(truncate(makedist('Normal', 40000, 66000),4000, 400000)); random(truncate(makedist('Normal', 10000, 16500),1000, 100000))];
    %findParams(:,i) = makedist('Lognormal', 130000, 214500); makedist('Lognormal', 21000, 34650); makedist('Lognormal', 40000, 66000); makedist('Lognormal', 40000, 66000); makedist('Lognormal', 10000, 16500);

    y0(1,1) = findParams(1,i); %C8 initial
    y0(3,1) = findParams(2,i); %C3
    y0(5,1) = findParams(3,i); %IAP
    y0(7,1) = findParams(4,i); %BAR
    y0(9,1) = findParams(5,i); %receptor
    
 
end

figure(1)
subplot(2,5,1); histogram(findParams(1,:),'FaceColor','black'); title('C8');
subplot(2,5,2); histogram(findParams(2,:),'FaceColor','black'); title('C3');
subplot(2,5,3); histogram(findParams(3,:),'FaceColor','black'); title('IAP');
subplot(2,5,4); histogram(findParams(4,:),'FaceColor','black'); title('BAR');
subplot(2,5,5); histogram(findParams(5,:),'FaceColor','black'); title('receptor');

sgtitle('10,000 Normal Distribution Samples for Non-Zero Initial Protein Concentrations');
%%

%     0, %C8a initial 
%     21000, %C3
%          0, %C3a
%       40000, %IAP
%       0,%C3aIAPinitial 
%       40000%BAR inital concentration
%       0,%C8aBAR inital concentration
%       % initial receptor and ligand from Hua 2005
%           10, % receptor
%            0,  % complex
%            1]; % ligand
  
% parameters
p = zeros(32,1);
p(1) = 5.8e-5; %cell.min-1
p(2) = 0;
p(3) = 1e-5;
p(4) = 0;
p(5) = 5e-4;
p(6) = 2.1e-1;
p(7) = 3e-3;
p(8) = 0;
p(9) = 5.8e-3;
p(10) = 0;
p(11) = 5.8e-3;
p(12) = 0;
p(13)= 1.73e-2;
p(14)= 0;
p(15) = 1.16e-2;
p(16) = 464;
p(17) = 3.9e-3;
p(18) = 507;
p(19) = 3.9e-3;
p(20) = 81.9;
p(21) = 5e-4;
p(22) = 2.1e-1;
p(23) = 1e-3;
p(24) = 40;
p(25) = 1.16e-2;
p(26) = 0;

% p27 and 28 based on FasL binding:
%           Kd = 0.4 nM; assume koff = 1.2e-2/min
p(27) = 1e-6; % ligand-receptor association **** does this generate a delay in C3a?
p(28) = 1.2e-2; % ligand-receptor dissociation %assumption -change

% p29 based on FasL signaling from Wu and Finley 2017
p(29) = 8.04e-5; % ligand-receptor complex signaling strength  **** does this generate a delay in C3a?

p(30) = 0; % receptor production

% p31 and p32 based on FasL signaling from Wu and Finley 2017
p(31) = 1.3e-3; % receptor degredation %assumption  -change    **** does this generate a delay in C3a?
p(32) = 4.67e-6;%4.67e-6 original; % receptor-ligand complex internalization -change    **** does this generate a delay in C3a?

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

tstep = 1;
endTime = 2000; %2000
tspan = [0:tstep:endTime];

%% - plot baseline results 
numSpecies = 11;
y0 = zeros(numSpecies,1);
% dydt = zeros(length(tspan), numSpecies, 10000);
for i=1:numCells
    y0(1,1) = findParams(1,i); %C8 initial
    y0(3,1) = findParams(2,i); %C3
    y0(5,1) = findParams(3,i); %IAP
    y0(7,1) = findParams(4,i); %BAR
    y0(9,1) = findParams(5,i); %receptor
    y0(11,1) = 6000; %ligand was 100
    [t, dydt_i] = ode15s(@newEISSINGodeModel_,tspan, y0,options,p);
    dydt{i,1} = dydt_i;
 
end
%%
for i = 1:numCells
    C3a_cells(:,i)= dydt{i,1}(:,4);
%     hold on
%     [m,time_i] = max(dydt{i,1}(:,4));
%     max_C3a_low(i,1) = m;
%     time_max_C3a_low(i,1) = time_i;
end
%%
for i = 1:numCells
    plot(t,dydt{i,1}(:,4),'LineWidth',1, 'Color', '#EDB120'); %#EDB120 #D95319
    hold on
    title('Uniform Distribution: C3a over time (50 cells)');  
    xlabel('time (min)');
    xlim([0 200]);
    ylim([0 200000])
    ylabel('C3a (molecules)');
end
%%

        plot(t,dydt{:,1}(:,4),'LineWidth',1, 'Color', 'black')
        %plot(t,dydt{i,1}(:,j),'-g','LineWidth',1)
        %hold on
        title('C3 concentration');  
        xlabel('time (min)');
        xlim([0 200]);
        ylabel('concentration of species');
%%
%% find number of apoptotic cells per threshold 
timesOfInterest = [5 10 15 20 30 60 90 120 150 180 210 240 270 300 330 360 390 420 450 480 510 540 570 600 700 800 900 1000 1500 2000]';
numApop = zeros(size(timesOfInterest,1),1);
for j = 1:size(timesOfInterest,1)
    apop = 0;
    for i=1:numCells
        timetoCheck = timesOfInterest(j);
        C3alevel = dydt{i,1}(timetoCheck+1,4);
        if C3alevel > 40000 %40000 %100000 %160000
             apop = apop+1;
        end
    end
    numApop(j,1) = apop;
end
%% find max OR index where threshold is first reached FIGURE 5
%% C8
idx = zeros(numCells,1);
% C3nonapop = zeros(findParams(1,:),1);
% C3apop = zeros(findParams(1,:),1);

for i =1:numCells
    C3a_cells = dydt{i,1}(:,4); %C3a (rows, columns)
    try
    idx(i)= find(C3a_cells>40000,1,"first");
    end
    if idx(i) == 0
        C8nonapop(i) = findParams(1,i);
    else
        C8apop(i) = findParams(1,i);
    end
%   C3apop(i,1)= apop;
%   C3nonapop(i,1)= nonapop;
%     [m,time_i] = max(dydt{i,1}(:,4));
%     max_C3a(i,1) = m;
%     time_max_C3a(i,1) = time_i;
end

%% C3
idx = zeros(numCells,1);
% C3nonapop = zeros(findParams(1,:),1);
% C3apop = zeros(findParams(1,:),1);

for i =1:numCells
    C3a_cells = dydt{i,1}(:,4); %C3a (rows, columns)
    try
    idx(i)= find(C3a_cells>40000,1,"first");
    end
    if idx(i) == 0
        C3nonapop(i) = findParams(2,i) ;
    else 
        C3apop(i) = findParams(2,i);
    end

end
%%
subplot(2,2,1); histogram(C8nonapop,'FaceColor','black'); title('nonapop');
subplot(2,2,2); histogram(C8apop,'FaceColor','black'); title('apop');
%% IAP
idx = zeros(numCells,1);
% C3nonapop = zeros(findParams(1,:),1);
% C3apop = zeros(findParams(1,:),1);

for i =1:numCells
    C3a_cells = dydt{i,1}(:,4); %C3a (rows, columns)
    try
    idx(i)= find(C3a_cells>40000,1,"first");
    end
    if idx(i) == 0
        IAPnonapop(i) = findParams(3,i) ;
    else 
        IAPapop(i) = findParams(3,i);
    end

end

%% BAR
idx = zeros(numCells,1);

for i =1:numCells
    C3a_cells = dydt{i,1}(:,4); %C3a (rows, columns)
    try
    idx(i)= find(C3a_cells>40000,1,"first");
    end
    if idx(i) == 0
        BARnonapop(i) = findParams(4,i);
    else 
        BARapop(i) = findParams(4,i);
    end

end

%% Receptor
idx = zeros(numCells,1);

for i =1:numCells
    C3a_cells = dydt{i,1}(:,4); %C3a (rows, columns)
    try
    idx(i)= find(C3a_cells>40000,1,"first");
    end
    if idx(i) == 0
        ReceptorNonapop(i) = findParams(5,i);
    else 
        ReceptorApop(i) = findParams(5,i);
    end

end

%% ATTEMPT FOR FIGURE 3
% [t, dydt_i] = ode15s(@newEISSINGodeModel_,[0:tstep:10], y0,options,p);
% dydt{i,1} = dydt_i;
% for i = 1:numCells
%     y0(11,1) = 600;
%     C3a_cells(:,i) = dydt{i,1}(:,4); %C3a
%     [m,time_i] = max(dydt{i,1}(:,4));
%     max_C3a(i,1) = m;
%     time_max_C3a(i,1) = time_i;
% end

%% -Figure 2: New Experiment
d= mean(max_C3a);
x=1;
errlow = mean(max_C3a)-min(max_C3a);
errhigh = max(max_C3a)-mean(max_C3a);

figure()

bar(d,'yellow')
hold on

er = errorbar(x,d,errlow,errhigh);    
er.Color = [0 0 0];   
er.LineStyle = 'none';  
ylabel('Level of C3a');
title('Mean Maximum level of C3a')
hold off

%% Figure 2b
n= mean(time_max_C3a);
errlowt= mean(time_max_C3a)-min(time_max_C3a);
errhight = max(time_max_C3a)-mean(time_max_C3a);

figure()
bar(x,n,'FaceColor',[0.929411764705882 0.694117647058824 0.125490196078431]);
hold on 

er = errorbar(x,n,errlowt,errhight);    
er.Color = [0 0 0];   
e.CapSize = 15;
er.LineStyle = 'none';  
title('Mean Time to Reach Maximum level of C3a')
hold off

%% - plot baseline results with lowest ligand
numSpecies = 11;
y0 = zeros(numSpecies,1);
% dydt = zeros(length(tspan), numSpecies, 10000);
for i=1:numCells 
    y0(1,1) = findParams(1,i); %C8 initial
    y0(3,1) = findParams(2,i); %C3
    y0(5,1) = findParams(3,i); %IAP
    y0(7,1) = findParams(4,i); %BAR
    y0(9,1) = findParams(5,i); %receptor
    y0(11,1) = 600; %ligand
    [t, dydt_i] = ode15s(@newEISSINGodeModel_,tspan, y0,options,p);
    dydt{i,1} = dydt_i;
end
%%
%Low Ligand 
for i = 1:numCells
     [t, dydt_i] = ode15s(@newEISSINGodeModel_,tspan, y0,options,p);
    dydt{i,1} = dydt_i;
    y0(11,1) = 600;
    C3a_cells(:,i) = dydt{i,1}(:,4); %C3a
    [m,time_i] = max(dydt{i,1}(:,4));
    max_C3a_low(i,1) = m;
    time_max_C3a_low(i,1) = time_i;
end

%%
for i = 1:numCells
    y0(11,1) = 100;
    C8_cells(:,i) = dydt{i,1}(:,1); %C8
    [m,time_i] = max(dydt{i,1}(:,1));
    max_C8(i,1) = m;
    time_max_C8(i,1) = time_i;
end
%%
for i = 1:numCells
    C8a_cells(:,i) = dydt{i,1}(:,2); %C8a
    [m,time_i] = max(dydt{i,1}(:,2));
    max_C8a(i,1) = m;
    time_max_C8a(i,1) = time_i;
end

%% 
for i = 1:numCells
    C3_cells(:,i) = dydt{i,1}(:,3); %C3
    [m,time_i] = max(dydt{i,1}(:,3));
    max_C3(i,1) = m;
    time_max_C3(i,1) = time_i;
end
%%
for i = 1:numCells
    IAP_cells(:,i) = dydt{i,1}(:,5); %IAP
    [m,time_i] = max(dydt{i,1}(:,5));
    max_IAP(i,1) = m;
    time_max_IAP(i,1) = time_i;
end
%%
for i = 1:numCells
    C3aIAP_cells(:,i) = dydt{i,1}(:,6); %C3aIAP
    [m,time_i] = max(dydt{i,1}(:,6));
    max_C3aIAP(i,1) = m;
    time_max_C3aIAP(i,1) = time_i;
end


for i = 1:numCells
    BAR_cells(:,i) = dydt{i,1}(:,7); %BAR
    [m,time_i] = max(dydt{i,1}(:,7));
    max_BAR(i,1) = m;
    time_max_BAR(i,1) = time_i;
end


for i = 1:numCells
    C8aBAR_cells(:,i) = dydt{i,1}(:,8); %C8aBAR
    [m,time_i] = max(dydt{i,1}(:,8));
    max_C8aBAR(i,1) = m;
    time_max_C8aBAR(i,1) = time_i;
end

for i = 1:numCells
    receptor_cells(:,i) = dydt{i,1}(:,9); %receptor
    [m,time_i] = max(dydt{i,1}(:,9));
    max_receptor(i,1) = m;
    time_max_receptor(i,1) = time_i;
end

for i = 1:numCells
    complex_cells(:,i) = dydt{i,1}(:,10); %complex
    [m,time_i] = max(dydt{i,1}(:,10));
    max_complex(i,1) = m;
    time_max_complex(i,1) = time_i;
end

for i = 1:numCells
    ligand_cells(:,i) = dydt{i,1}(:,11); %ligand
    [m,time_i] = max(dydt{i,1}(:,11));
    max_ligand(i,1) = m;
    time_max_ligand(i,1) = time_i;
end

%% create histograms
figure(2); %C3a
subplot(2,1,1);
histogram(max_C3a,100,'FaceColor','[0.8500 0.3250 0.0980]');
xlabel('Level of Activated caspase 3');
ylabel('Number of Simulations');
title('Maximum C3a level')

subplot(2,1,2);
histogram(time_max_C3a,100,'FaceColor','[0.9290 0.6940 0.1250]');
xlabel('Level of Activated caspase 3');
ylabel('Number of Simulations');
title('Time to Reach Maximum C3a level')
%%
figure(3);%C3
subplot(2,1,1);
histogram(max_C3,100);

title('Maximum C3 level')

subplot(2,1,2);
histogram(time_max_C3,100);
title('Time to Reach Maximum C3 level')
%%
figure(4); %C8
subplot(2,1,1);
histogram(max_C8,100);
title('Maximum C8 level')

subplot(2,1,2);
histogram(time_max_C8,100);
title('Time to Reach Maximum C8 level')
%% 

figure(5);%C8a
subplot(2,1,1);
histogram(max_C8a,100);
title('Maximum C8a level')

subplot(2,1,2);
histogram(time_max_C8a,100);
title('Time to Reach Maximum C8a level')
%% 
figure(6);%IAP
subplot(2,1,1);
histogram(max_IAP,100);
title('Maximum IAP level')

subplot(2,1,2);
histogram(time_max_IAP,100);
title('Time to Reach Maximum IAP level')
%% 
figure(7);%C3aIAP
subplot(2,1,1);
histogram(max_C3aIAP,100);
title('Maximum C3aIAP level')

subplot(2,1,2);
histogram(time_max_C3aIAP,100);
title('Time to Reach Maximum C3aIAP level')
%% 

figure(8);%BAR
subplot(2,1,1);
histogram(max_BAR,100);
title('Maximum BAR level')

subplot(2,1,2);
histogram(time_max_BAR,100);
title('Time to Reach Maximum BAR level')
%% 
figure(9);%C8aBAR
subplot(2,1,1);
histogram(max_C8aBAR,100);
title('Maximum C8aBAR level')

subplot(2,1,2);
histogram(time_max_C8aBAR,100);
title('Time to Reach Maximum C8aBAR level')
%%

figure(10);%receptor
subplot(2,1,1);
histogram(max_receptor,100);
title('Maximum receptor level')

subplot(2,1,2);
histogram(time_max_receptor,10);
title('Time to Reach Maximum receptor level')
%%

figure(11);%complex
subplot(2,1,1);
histogram(max_complex,100);
title('Maximum complex level')

subplot(2,1,2);
histogram(time_max_complex,100);
title('Time to Reach Maximum complex level')
%% 
figure(12);%ligand
subplot(2,1,1);
histogram(max_ligand,100);
title('Maximum ligand level')

subplot(2,1,2);
histogram(time_max_ligand,100);
title('Time to Reach Maximum ligand level')
%%
%%
C3a_cells(:,i) = dydt{i,1}(:,4); %C3a
%% plotting the species concentrations (INCLUDES FIGURE 1)
 y0(11,1) = 600;
for j = 1:11 % number of species
%     subplot(4,3,j)
    figure(13+j);
    for i=1:numCells %numCells
        plot(t,dydt{i,1}(:,j),'LineWidth',1, 'Color', 'black')
        %plot(t,dydt{i,1}(:,j),'-g','LineWidth',1)
        hold on
    end
    title(speciesNames{j,1});  
    xlabel('time (min)');
    xlim([0 800]);
    ylabel('concentration of species');
%     sgtitle('Protein Distribution over Time')
end

% plot(t,dydt(i,j))
% hold on
% plot(t,dydt(i,j))
% hold on 
% plot(t,dydt(i,j))
% hold on 
% plot(t,dydt(i,j))
% hold on 
% plot(t,dydt(i,j))

%% - varying receptor number with [Ligand] = 1
tstep = 1;
endTime = 1000; %2000;
tspan = [0:tstep:endTime];
%receptor = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];
%receptor = linspace(100,1000000,20);
receptor = logspace(2,5,20); %keep it in logspace

figure(25)
% subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('varying receptor number with [Ligand] = 1')
% subplot(2,2,2)
% xlabel('time(min)')
% hold on
% title('L:R Complex')
% subplot(2,2,3)
% xlabel('time (min)')
% hold on
% title('Receptor')
% subplot(2,2,4)
% xlabel('time(min)')
% ylim([0 100])
% hold on
% title('Ligand')

y0(11,1) = 1;
for i = 1:size(receptor,2)
%     disp(i)
    y0(9,1) = receptor(i);
    [t, dydt] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(receptor,2);
%     subplot(2,2,1)
    plot(t,dydt(:,4), 'color', [c, 0, 1-c], 'LineWidth', 2)
%     subplot(2,2,2)
%     plot(t,dydt(:,10), 'color', [c, 0, 1-c],'LineWidth', 2 )
%     subplot(2,2,3)
%     plot(t,dydt(:,9), 'color', [c, 0, 1-c],'LineWidth', 2)
%     subplot(2,2,4)
%     plot(t,dydt(:,11), 'color', [c, 0, 1-c],'LineWidth', 2)
end
%sgtitle('varying receptor number with [Ligand] = 1');
xlim([0 800]); 
digits(7);
legend(num2str(receptor',digits));
colormap(jet(256));

%% - varying receptor number with [Ligand] = 10

%receptor = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];
receptor = linspace(100,1000000,20);

figure(26)
subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('C3*')
subplot(2,2,2)
xlabel('time(min)')
hold on
title('L:R Complex')
subplot(2,2,3)
xlabel('time (min)')
hold on
title('Receptor')
subplot(2,2,4)
xlabel('time(min)')
ylim([0 100])
hold on
title('Ligand')

y0(11,1) = 10;
for i = 1:size(receptor,2)
%     disp(i)
    y0(9,1) = receptor(i);
    [t, dydt] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(receptor,2);
    subplot(2,2,1)
    plot(t,dydt(:,4), 'color', [1-c, 0, c], 'LineWidth', 2)
    subplot(2,2,2)
    plot(t,dydt(:,10), 'color', [1-c, 0, c],'LineWidth', 2 )
    subplot(2,2,3)
    plot(t,dydt(:,9), 'color', [1-c, 0, c],'LineWidth', 2)
    subplot(2,2,4)
    plot(t,dydt(:,11), 'color', [1-c, 0, c] ,'LineWidth', 2)
end
sgtitle('varying receptor number with [Ligand] = 10');
xlim([0 800]); 
digits(7);
legend(num2str(receptor',digits));
%% - varying receptor number with [Ligand] = 100

%receptor = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];
%receptor = linspace(100,1000000,20);
receptor = logspace(2,5,20);

figure(27)
%subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('varying receptor number with [Ligand] = 100')
% subplot(2,2,2)
% xlabel('time(min)')
% hold on
% title('L:R Complex')
% subplot(2,2,3)
% xlabel('time (min)')
% hold on
% title('Receptor')
% subplot(2,2,4)
% xlabel('time(min)')
% ylim([0 100])
% hold on
% title('Ligand')

y0(11,1) = 100;
for i = 1:size(receptor,2)
%     disp(i)
    y0(9,1) = receptor(i);
    [t, dydt] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(receptor,2);
    %subplot(2,2,1)
    plot(t,dydt(:,4), 'color', [c, 0, 1-c], 'LineWidth', 2)
%     subplot(2,2,2)
%     plot(t,dydt(:,10), 'color', [c, 0, 1-c],'LineWidth', 2 )
%     subplot(2,2,3)
%     plot(t,dydt(:,9), 'color', [c, 0, 1-c],'LineWidth', 2)
%     subplot(2,2,4)
%     plot(t,dydt(:,11), 'color', [c, 0, 1-c] ,'LineWidth', 2)
end
%sgtitle('varying receptor number with [Ligand] = 100');
xlim([0 800]); 
digits(7);
legend(num2str(receptor',digits));
%% - varying ligand with [receptor] = 100

%ligand = linspace(60,60000,20);
ligand = logspace(2,5,20)
%ligand = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];

figure(28)
%subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('varying ligand number with [Receptor] = 100')
% subplot(2,2,2)
% xlabel('time(min)')
% hold on
% title('L:R Complex')
% subplot(2,2,3)
% xlabel('time (min)')
% hold on
% title('Receptor')
% subplot(2,2,4)
% xlabel('time(min)')
% hold on
% title('Ligand')

y0(9,1) = 100;
for i = 1:size(ligand,2)
%     disp(i)
    y0(11,1) = ligand(i);
    [t, dydt_i] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(ligand,2);
    %subplot(2,2,1)
%     plot(t,dydt(:,4), 'color', [c, 0, 1-c], 'LineWidth', 2)
%     subplot(2,2,2)
%     plot(t,dydt(:,10), 'color', [1-c, 0, c],'LineWidth', 2 )
%     subplot(2,2,3)
%     plot(t,dydt(:,9), 'color', [1-c, 0, c],'LineWidth', 2)
%     subplot(2,2,4)
%     plot(t,dydt(:,11), 'color', [1-c, 0, c] ,'LineWidth', 2)
% end
%sgtitle('varying ligand number with [Receptor] = 100');
xlim([0 800]); 
% ylim([0 40000]); 
% digits(7);
% legend(num2str(ligand',digits));

% for i = 1:numCells
%      [t, dydt_i] = ode15s(@newEISSINGodeModel_,tspan, y0,options,p);
    dydt{i,1} = dydt_i;
%     y0(11,1) = ligand(i);
    C3a_cells(:,i) = dydt{i,1}(:,4); %C3a
    [m,time_i] = max(dydt{i,1}(:,4));
    max_C3a_low(i,1) = m;
    time_max_C3a_low(i,1) = time_i;
end
figure()
bar(time_i);
%%

%% - varying ligand with [receptor] = 1,000

%ligand = linspace(60,60000,20);
ligand = logspace(2,5,20)
%ligand = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];

figure(29)
subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('C3*')
subplot(2,2,2)
xlabel('time(min)')
hold on
title('L:R Complex')
subplot(2,2,3)
xlabel('time (min)')
hold on
title('Receptor')
subplot(2,2,4)
xlabel('time(min)')
hold on
title('Ligand')

y0(9,1) = 1000;
for i = 1:size(ligand,2)
%     disp(i)
    y0(11,1) = ligand(i);
    [t, dydt] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(ligand,2);
    subplot(2,2,1)
    plot(t,dydt(:,4), 'color', [1-c, 0, c], 'LineWidth', 2)
    subplot(2,2,2)
    plot(t,dydt(:,10), 'color', [1-c, 0, c],'LineWidth', 2 )
    subplot(2,2,3)
    plot(t,dydt(:,9), 'color', [1-c, 0, c],'LineWidth', 2)
    subplot(2,2,4)
    plot(t,dydt(:,11), 'color', [1-c, 0, c] ,'LineWidth', 2)
end
sgtitle('varying ligand number with [Receptor] = 1,000');
xlim([0 800]); 
ylim([0 40000]); 
digits(7);
legend(num2str(ligand',digits));
%% - varying ligand with [receptor] = 10,000

%ligand = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];
%ligand = linspace(60,60000,20);
ligand = logspace(2,5,20);

figure(30)
subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('C3*')
subplot(2,2,2)
xlabel('time(min)')
hold on
title('L:R Complex')
subplot(2,2,3)
xlabel('time (min)')
hold on
title('Receptor')
subplot(2,2,4)
xlabel('time(min)')
hold on
title('Ligand')

y0(9,1) = 10000;
for i = 1:size(ligand,2)
%     disp(i)
    y0(11,1) = ligand(i);
    [t, dydt] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(ligand,2);
    subplot(2,2,1)
    plot(t,dydt(:,4), 'color', [1-c, 0, c], 'LineWidth', 2)
    subplot(2,2,2)
    plot(t,dydt(:,10), 'color', [1-c, 0, c],'LineWidth', 2 )
    subplot(2,2,3)
    plot(t,dydt(:,9), 'color', [1-c, 0, c],'LineWidth', 2)
    subplot(2,2,4)
    plot(t,dydt(:,11), 'color', [1-c, 0, c] ,'LineWidth', 2)
end
sgtitle('varying ligand number with [Receptor] = 10,000');
xlim([0 800]); 
ylim([0 40000]); 
digits(7);
legend(num2str(ligand',digits));
%% - varying ligand with [receptor] = 100,000

%ligand = linspace(60,60000,20);
ligand= logspace(2,5,20);
%ligand = [0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,2,3,4,5,6,7,8,9,10];

figure(31)
%subplot(2,2,1)
xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('varying ligand number with [Receptor] = 100,000')
% subplot(2,2,2)
% xlabel('time(min)')
% hold on
% title('L:R Complex')
% subplot(2,2,3)
% xlabel('time (min)')
% hold on
% title('Receptor')
% subplot(2,2,4)
% xlabel('time(min)')
% hold on
% title('Ligand')

y0(9,1) = 100000;
for i = 1:size(ligand,2)
%     disp(i)
    y0(11,1) = ligand(i);
    [t, dydt] = ode15s(@newEISSINGodeModel_,tspan,y0,options,p);
    c = i/size(ligand,2);
    %subplot(2,2,1)
    plot(t,dydt(:,4), 'color', [c, 0, 1-c], 'LineWidth', 2)
%     subplot(2,2,2)
%     plot(t,dydt(:,10), 'color', [1-c, 0, c],'LineWidth', 2 )
%     subplot(2,2,3)
%     plot(t,dydt(:,9), 'color', [1-c, 0, c],'LineWidth', 2)
%     subplot(2,2,4)
%     plot(t,dydt(:,11), 'color', [1-c, 0, c] ,'LineWidth', 2)
end
%sgtitle('varying ligand number with [Receptor] = 100,000');
xlim([0 800]);   
digits(7);
legend(num2str(ligand',digits));
%%
%kk=size(receptor);
%blueGradientfixed = [linspace(lightBLUE(1), darkBLUE(1),kk)',linspace(lightBLUE(2), darkBLUE(2),kk)',linspace(lightBLUE(3), darkBLUE(3),kk),'linspace(lightBLUE(4), darkBLUE(4),kk)']
%plot(t,dydt(:,10),'color',blueGradientfixed(i,:))

%%
% figure(32)
% bar(max_C3a,100,'FaceColor','b');
% hold on
% bar(max_C3a_low,100,'FaceColor','g');
% res = [max_C3a_low,max_C3a]; 
% bar(res);
% sgtitle('Activated Caspase 3 Maximums at Peak over 100 Cells' );

% figure(33)
% res2 = [time_max_C3a_low,time_max_C3a];
% bar(res2);
% sgtitle('Time to Reach Activated Caspase 3 Maximums Over 100 Cells');


% figure(33)
% histogram([time_max_C3a;time_max_C3a_low],100,'FaceColor','g');

% figure(35)
% histogram([h2,h4],100, 'FaceColor', 'g');

% xlabel('Number of Simulations');
% ylabel('Time to Reach Activated Caspase 3')
%%

% % subplot(2,1,1);
% % histogram(res,100,'FaceColor','[0.8500 0.3250 0.0980]');

% 
% subplot(2,1,2);
% % histogram(res2,100,'FaceColor','[0.9290 0.6940 0.1250]');
% bar(time_max_C3a_low(i,1));
% hold on
% bar(time_max_C3a(i,1));






