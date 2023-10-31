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
%% Specify the Number of Cells You Want to Simulate
numCells = 10; %10000;

%% Finding Non-Zero Initial Protein Amounts with Uniform Distribution OR Normal Distribution (Supplemental Figures)
for i=1:numCells %10000 
    
   %For a uniform Distribution uncomment the next line
   findParams(:,i) = [findParamValue(130000,10);findParamValue(21000,10);findParamValue(40000,10);findParamValue(40000,10);findParamValue(10000,10)]; %uniform distribution
   %For a normal Distribution uncomment the next line and comment the prior
   %findParams(:,i) = [random(truncate(makedist('Normal', 130000,214500),13000, 1300000)); random(truncate(makedist('Normal', 21000,34650),2100, 210000)); random(truncate(makedist('Normal', 40000,66000),4000, 400000)); random(truncate(makedist('Normal', 40000, 66000),4000, 400000)); random(truncate(makedist('Normal', 10000,16500),1000, 100000))];

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
%% Assign Parameter Values and Time Steps

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

%% Run ODE with randomly selected initial protein amounts 
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
    [t, dydt_i] = ode15s(@ApoptosisODEModelPub,tspan, y0,options,p);
    dydt{i,1} = dydt_i;
 
end
%% Finding Time Where C3a Reaches Maximum and The Time it Occurs (Figure 4)
for i = 1:numCells
    C3a_cells(:,i)= dydt{i,1}(:,4);
    [m,time_i] = max(dydt{i,1}(:,4));
    max_C3a_low(i,1) = m;
    time_max_C3a_low(i,1) = time_i;
end
%% Plotting C3a concentration over Time for numCells (Figure 3)
for i = 1:numCells
    plot(t,dydt{i,1}(:,4),'LineWidth',1, 'Color', '#EDB120'); %#EDB120 #D95319
    hold on
    title('C3a over time');  
    xlabel('time (min)');
    xlim([0 200]);
    ylim([0 200000])
    ylabel('C3a (molecules)');
end
%% Find number of apoptotic cells per threshold (Figure 6)
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
%% Find max OR index where caspase threshold is first reached (Figure 5,7, and 8)
%% C8
idx = zeros(numCells,1);

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

end

%% C3
idx = zeros(numCells,1);

for i =1:numCells
    C3a_cells = dydt{i,1}(:,4); %C3a (rows, columns)
    try
    idx(i)= find(C3a_cells>40000,1,"first"); %Change depending on threshold to 40000, 100000, or 160000
    end
    if idx(i) == 0
        C3nonapop(i) = findParams(2,i) ;
    else 
        C3apop(i) = findParams(2,i);
    end

end

%% IAP
idx = zeros(numCells,1);

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


%% - Varying receptor with constant ligand 

receptor = logspace(2,5,20);

xlabel('time (min)')
ylabel('activated caspase 3 concentration')
hold on
title('varying receptor number with [Ligand] = 100')

y0(11,1) = 100;
%y0(11,1) = 100000;
for i = 1:size(receptor,2)
%     disp(i)
    y0(9,1) = receptor(i);
    [t, dydt] = ode15s(@ApoptosisODEModelPub,tspan,y0,options,p);
    c = i/size(receptor,2);

    plot(t,dydt(:,4), 'color', [c, 0, 1-c], 'LineWidth', 2)

end

xlim([0 800]); 
digits(7);
legend(num2str(receptor',digits));

%% - Varying ligand with [receptor] = 100 (Figure 2)

ligand = logspace(2,5,20);

y0(9,1) = 100;
for i = 1:size(ligand,2)
    y0(11,1) = ligand(i);
    [t, dydt_i] = ode15s(@ApoptosisODEModelPub,tspan,y0,options,p);

    dydt{i,1} = dydt_i;
    C3a_cells(:,i) = dydt{i,1}(:,4); %C3a
    [m,time_i] = max(dydt{i,1}(:,4));
    max_C3a_low(i,1) = m;
    time_max_C3a_low(i,1) = time_i;
end

%% - varying ligand with [receptor] = 100000 (Figure 2)

ligand = logspace(2,5,20);

y0(9,1) = 100000;
for i = 1:size(ligand,2)
    y0(11,1) = ligand(i);
    [t, dydt_i] = ode15s(@ApoptosisODEModelPub,tspan,y0,options,p);

    dydt{i,1} = dydt_i;
    C3a_cells(:,i) = dydt{i,1}(:,4); %C3a
    [m,time_i] = max(dydt{i,1}(:,4));
    max_C3a(i,1) = m;
    time_max_C3a(i,1) = time_i;
end
