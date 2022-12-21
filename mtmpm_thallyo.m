%% MODELO DE TRAJETORIA MASSA-PONTO MODIFICADO - THALLYO
% codigo desenvolvido e testado na versao r2021b do matlab
% para simulacao de trajetoria da municao 114.3mm inerte
% e com propelente base bleed na versao 2021
clc
clear
tic

%% DADOS 

%--- constantes
g0 = 9.80665;               %aceleracao padrao da gravidade
RT = 6.356766e6;            %raio da Terra (STANAG)
omega = 7.292115e-5;        %velocidade angular da Terra

%--- condicoes atmosfericas
P0 = 1.01325e5;             %pressao ao nivel do mar
T0 = 288.15;                %temperatura do ar ao nivel do mar
rho0 = 1.225;               %densidade do ar ao nivel do mar
R = 287.05;                 %constante de gas ideal (ar)
K = 1.4;                    %constante de expansao adiabatica (ar)
B = -6.5e-3;                %gradiente de temperatura
Pb = 22.632e3;              %pressao a 11000m de altitude 

%--- condicoes balisticas
d = 0.154700;               %diametro de referencia
db = 0.133731;              %diametro de boattail
m0 = 43.5450;               %massa inicial do projetil
dm1_0015 = -0.015213252;    %taxa de queima do propelente 1 - 1pol
dm1_0030 = -0.030426503;    %taxa de queima do propelente 2 - 1pol
dm2_0030 = dm1_0030;        %taxa de queima do propelente 1 - 2pol
dm2_0060 = 0.060853006;     %taxa de queima do propelente 2 - 2pol
Ix0 = 0.142453;             %momento de inercia inicial em X
Iy0 = 1.2699;               %momento de inercia inicial em Y
CG0 = 0.45835;              %localizacao inicial do centro de gravidade
m1 = 42.9850;               %massa final do projetil
Ix1 = 0.13730;              %momento de inercia final em X
Iy1 = 1.2699;               %momento de inercia final em Y
CG1 = 0.44645;              %localizacao final do centro de gravidade

%--- condicoes iniciais de tiro
u0 = 878;                   %velocidade de boca
Az = 0;                     %azimute
lat = -23*pi/180;           %latitude do Rio de Janeiro
tc = 25;                    %taxa de torcao no tubo
QD = 1.2;                   %fator de ajuste da forca de arrasto
QM = 1.2;                   %fator de ajuste da forca de Magnus

%% CONDICOES AERODINAMICAS - PRODAS + CFD
tab_coef_in = readtable('coef_aero_in.xlsx');
tab_coef_bb_1pol_0015 = readtable('coef_aero_bb_2306K_1pol_0015.xlsx');
tab_coef_bb_1pol_0030 = readtable('coef_aero_bb_2306K_1pol_0030.xlsx');
tab_coef_bb_2pol_0030 = readtable('coef_aero_bb_2306K_2pol_0030.xlsx');
tab_coef_bb_2pol_0060 = readtable('coef_aero_bb_2306K_2pol_0060.xlsx');

%--- condicoes iniciais 
pos0 = zeros(1,3);                                             
spn0 = (2*pi*u0)/(tc*d);
pmd1pol_0015 = [m0, Ix0, CG0];
pmd1pol_0030 = pmd1pol_0015;
pmd2pol_0030 = pmd1pol_0015;
pmd2pol_0060 = pmd2pol_0030;
acc = pos0;

%% SOLUCAO NUMERICA
% Para QE = 40 graus (711 mil)
% Primeira conta admite que seria sem base bleed (inerte)
% Para 1 polegada, considerar 2 vazoes: 0.015 kg/s e 0.030 kg/s
% Para 2 polegada, considerar 2 vazoes: 0.030 kg/s e 0.060 kg/s

% convertendo de mil para radianos
QE = 0.0009817477*711.1;
% velocidade inicial de disparo
vel0 = [u0*cos(QE)*cos(Az), u0*sin(QE), u0*cos(QE)*sin(Az)];

%% IMPLEMENTACAO DAS FUNCOES RK4 PARA A TRAJETORIA - QE = 40 GRAUS:
% inerte
[VELin_711, POSin_711, SPNin_711, Tin_711] = RK4_IN(vel0, pos0, spn0);
% com base bleed - 1 polegada na saida
[VELbb1pol0015_711, POSbb1pol0015_711, SPNbb1pol0015_711, PMDbb1pol0015_711, Tbb1pol0015_711] = RK4_BB_1pol_0015(vel0, pos0, spn0, pmd1pol_0015);
[VELbb1pol0030_711, POSbb1pol0030_711, SPNbb1pol0030_711, PMDbb1pol0030_711, Tbb1pol0030_711] = RK4_BB_1pol_0030(vel0, pos0, spn0, pmd1pol_0030);
% com base bleed - 2 polegadas na saída
[VELbb2pol0030_711, POSbb2pol0030_711, SPNbb2pol0030_711, PMDbb2pol0030_711, Tbb2pol0030_711] = RK4_BB_2pol_0030(vel0, pos0, spn0, pmd2pol_0030);
[VELbb2pol0060_711, POSbb2pol0060_711, SPNbb2pol0060_711, PMDbb2pol0060_711, Tbb2pol0060_711] = RK4_BB_2pol_0060(vel0, pos0, spn0, pmd2pol_0060);

%% CALCULANDO O ALCANCE HORIZONTAL PARA GERAR TABELA DE RESULTADOS FINAIS - QE = 40 GRAUS 
% inerte
alcin_711 = interp1(POSin_711(:,2),POSin_711(:,1),0);
% com base bleed - 1 polegada na saída
alcbb1pol0015_711 = interp1(POSbb1pol0015_711(:,2),POSbb1pol0015_711(:,1),0);
alcbb1pol0030_711 = interp1(POSbb1pol0030_711(:,2),POSbb1pol0030_711(:,1),0);
aumnt1pol0015_711 = (alcbb1pol0015_711-alcin_711)/alcin_711*100;
aumnt1pol0030_711 = (alcbb1pol0030_711-alcin_711)/alcin_711*100;
% com base bleed - 2 polegadas na saída
alcbb2pol0030_711 = interp1(POSbb2pol0030_711(:,2),POSbb2pol0030_711(:,1),0);
alcbb2pol0060_711 = interp1(POSbb2pol0060_711(:,2),POSbb2pol0060_711(:,1),0);
aumnt2pol0030_711 = (alcbb2pol0030_711-alcin_711)/alcin_711*100;
aumnt2pol0060_711 = (alcbb2pol0060_711-alcin_711)/alcin_711*100;

%% CALCULANDO O APOGEU PARA GERAR TABELA DE RESULTADOS FINAIS - QE = 40 GRAUS
% inerte
apgin_711 = max(POSin_711(:,2));
% com base bleed - 1 polegada na saída
apgbb1pol0015_711 = max(POSbb1pol0015_711(:,2));
apgbb1pol0030_711 = max(POSbb1pol0030_711(:,2));
amtapg1pol0015_711 = (apgbb1pol0015_711-apgin_711)/apgin_711*100;
amtapg1pol0030_711 = (apgbb1pol0030_711-apgin_711)/apgin_711*100;
t_in_711 = max(Tin_711);
t_bb_1pol0015_711 = max(Tbb1pol0015_711);
t_bb_1pol0030_711 = max(Tbb1pol0030_711);
% com base bleed - 2 polegada na saída
apgbb2pol0030_711 = max(POSbb2pol0030_711(:,2));
apgbb2pol0060_711 = max(POSbb2pol0060_711(:,2));
amtapg2pol0030_711 = (apgbb2pol0030_711-apgin_711)/apgin_711*100;
amtapg2pol0060_711 = (apgbb2pol0060_711-apgin_711)/apgin_711*100;
t_in_711 = max(Tin_711);
t_bb_2pol0030_711 = max(Tbb2pol0030_711);
t_bb_2pol0060_711 = max(Tbb2pol0060_711);

% nova conversao de mils para radianos
QE = 0.0009817477*800;
% mais uma vez programando a condicao de disparo
vel0 = [u0*cos(QE)*cos(Az), u0*sin(QE), u0*cos(QE)*sin(Az)];

%% IMPLEMENTACAO DAS FUNCOES RK4 PARA A TRAJETORIA - QE = 45 GRAUS:
% inerte
[VELin_800, POSin_800, SPNin_800, Tin_800] = RK4_IN(vel0, pos0, spn0);
% com base bleed - 1 polegada na saida
[VELbb1pol0015_800, POSbb1pol0015_800, SPNbb1pol0015_800, PMDbb1pol0015_800, Tbb1pol0015_800] = RK4_BB_1pol_0015(vel0, pos0, spn0, pmd1pol_0015);
[VELbb1pol0030_800, POSbb1pol0030_800, SPNbb1pol0030_800, PMDbb1pol0030_800, Tbb1pol0030_800] = RK4_BB_1pol_0030(vel0, pos0, spn0, pmd1pol_0030);
% com base bleed - 2 polegadas na saída
[VELbb2pol0030_800, POSbb2pol0030_800, SPNbb2pol0030_800, PMDbb2pol0030_800, Tbb2pol0030_800] = RK4_BB_2pol_0030(vel0, pos0, spn0, pmd2pol_0030);
[VELbb2pol0060_800, POSbb2pol0060_800, SPNbb2pol0060_800, PMDbb2pol0060_800, Tbb2pol0060_800] = RK4_BB_2pol_0060(vel0, pos0, spn0, pmd2pol_0060);

%% CALCULANDO O ALCANCE HORIZONTAL PARA GERAR TABELA DE RESULTADOS FINAIS - QE = 45 GRAUS 
% inerte
alcin_800 = interp1(POSin_800(:,2),POSin_800(:,1),0);
% com base bleed - 1 polegada na saída
alcbb1pol0015_800 = interp1(POSbb1pol0015_800(:,2),POSbb1pol0015_800(:,1),0);
alcbb1pol0030_800 = interp1(POSbb1pol0030_800(:,2),POSbb1pol0030_800(:,1),0);
aumnt1pol0015_800 = (alcbb1pol0015_800-alcin_800)/alcin_800*100;
aumnt1pol0030_800 = (alcbb1pol0030_800-alcin_800)/alcin_800*100;
% com base bleed - 2 polegadas na saída
alcbb2pol0030_800 = interp1(POSbb2pol0030_800(:,2),POSbb2pol0030_800(:,1),0);
alcbb2pol0060_800 = interp1(POSbb2pol0060_800(:,2),POSbb2pol0060_800(:,1),0);
aumnt2pol0030_800 = (alcbb2pol0030_800-alcin_800)/alcin_800*100;
aumnt2pol0060_800 = (alcbb2pol0060_800-alcin_800)/alcin_800*100;

%% CALCULANDO O APOGEU PARA GERAR TABELA DE RESULTADOS FINAIS - QE = 40 GRAUS
% inerte
apgin_800 = max(POSin_800(:,2));
% com base bleed - 1 polegada na saída
apgbb1pol0015_800 = max(POSbb1pol0015_800(:,2));
apgbb1pol0030_800 = max(POSbb1pol0030_800(:,2));
amtapg1pol0015_800 = (apgbb1pol0015_800-apgin_800)/apgin_800*100;
amtapg1pol0030_800 = (apgbb1pol0030_800-apgin_800)/apgin_800*100;
t_in_800 = max(Tin_800);
t_bb_1pol0015_800 = max(Tbb1pol0015_800);
t_bb_1pol0030_800 = max(Tbb1pol0030_800);
% com base bleed - 2 polegada na saída
apgbb2pol0030_800 = max(POSbb2pol0030_800(:,2));
apgbb2pol0060_800 = max(POSbb2pol0060_800(:,2));
amtapg2pol0030_800 = (apgbb2pol0030_800-apgin_800)/apgin_800*100;
amtapg2pol0060_800 = (apgbb2pol0060_800-apgin_800)/apgin_800*100;
t_in_800 = max(Tin_800);
t_bb_2pol0030_800 = max(Tbb2pol0030_800);
t_bb_2pol0060_800 = max(Tbb2pol0060_800);

%% PLOTAGEM
% CANTO SUPERIOR ESQUERDO - QE = 40 GRAUS P/ 1 POLEGADA
subplot(2,2,1);
plot(POSin_711(:,1),POSin_711(:,2),'linewidth',2);
xlabel('\fontname{Times} Downrange (m)');
ylabel('\fontname{Times} Height (m)');
title('Quadrant Elevation = 711.1 mil ($\phi_{bb}$ = 25.4mm)','interpreter','latex');
axis([0, 30000, 0, 12500]);
hold on
plot(POSbb1pol0015_711(:,1),POSbb1pol0015_711(:,2),'--','linewidth',2);
plot(POSbb1pol0030_711(:,1),POSbb1pol0030_711(:,2),':','linewidth',2);
legend('Sem base bleed','$\dot{m}_{bb}$ = 0.015 kg/s','$\dot{m}_{bb}$ = 0.030 kg/s','interpreter','latex');
% CANTO SUPERIOR DIREITO - QE = 45 GRAUS P/ 1 POLEGADA
subplot(2,2,2);
plot(POSin_800(:,1),POSin_800(:,2),'linewidth',2);
xlabel('\fontname{Times} Downrange (m)');
ylabel('\fontname{Times} Height (m)');
title('Quadrant Elevation = 800 mil ($\phi_{bb}$ = 25.4mm)','interpreter','latex');
axis([0, 30000, 0, 12500]);
hold on
plot(POSbb1pol0015_800(:,1),POSbb1pol0015_800(:,2),'--','linewidth',2);
plot(POSbb1pol0030_800(:,1),POSbb1pol0030_800(:,2),':','linewidth',2);
legend('Sem base bleed','$\dot{m}_{bb}$ = 0.015 kg/s','$\dot{m}_{bb}$ = 0.030 kg/s','Interpreter','latex');
% CANTO SUPERIOR ESQUERDO - QE = 40 GRAUS P/ 2 POLEGADAS
subplot(2,2,3);
plot(POSin_711(:,1),POSin_711(:,2),'linewidth',2);
xlabel('\fontname{Times} Downrange (m)');
ylabel('\fontname{Times} Height (m)');
title('Quadrant Elevation = 711.1 mil ($\phi_{bb}$ = 50.8mm)','interpreter','latex');
axis([0, 30000, 0, 12500]);
hold on
plot(POSbb2pol0030_711(:,1),POSbb2pol0030_711(:,2),'--','linewidth',2);
plot(POSbb2pol0060_711(:,1),POSbb2pol0060_711(:,2),':','linewidth',2);
legend('Sem base bleed','$\dot{m}_{bb}$ = 0.030 kg/s','$\dot{m}_{bb}$ = 0.060 kg/s','interpreter','latex');
% CANTO SUPERIOR DIREITO - QE = 45 GRAUS P/ 1 POLEGADA
subplot(2,2,4);
plot(POSin_800(:,1),POSin_800(:,2),'linewidth',2);
xlabel('\fontname{Times} Downrange (m)');
ylabel('\fontname{Times} Height (m)');
title('Quadrant Elevation = 800 mil ($\phi_{bb}$ = 50.8mm)','interpreter','latex');
axis([0, 30000, 0, 12500]);
hold on
plot(POSbb2pol0030_800(:,1),POSbb2pol0030_800(:,2),'--','linewidth',2);
plot(POSbb2pol0060_800(:,1),POSbb2pol0060_800(:,2),':','linewidth',2);
legend('Sem base bleed','$\dot{m}_{bb}$ = 0.030 kg/s','$\dot{m}_{bb}$ = 0.060 kg/s','Interpreter','latex');

%% FORMACAO DA TABELA
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20s %1s %20s %1s %20s %1s %20s %1s %20s %1s %20s %1s\n','|','Elevacao','|','Alcance(Inerte)','|','Alcance(BB 1pol)','|','Aumento(%)','|','Apogeu(Inerte)','|','Apogeu(BB 1pol)','|','Aumento(%)','|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','711 mil (15 g/s)','|',alcin_711,'|',alcbb1pol0015_711,'|',aumnt1pol0015_711,'|',apgin_711,'|',apgbb1pol0015_711,'|',amtapg1pol0015_711,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','711 mil (30 g/s)','|',alcin_711,'|',alcbb1pol0030_711,'|',aumnt1pol0030_711,'|',apgin_711,'|',apgbb1pol0030_711,'|',amtapg1pol0030_711,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','800 mil (15 g/s)','|',alcin_800,'|',alcbb1pol0015_800,'|',aumnt1pol0015_800,'|',apgin_800,'|',apgbb1pol0015_800,'|',amtapg1pol0015_800,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','800 mil (30 g/s)','|',alcin_800,'|',alcbb1pol0030_800,'|',aumnt1pol0015_800,'|',apgin_800,'|',apgbb1pol0030_800,'|',amtapg1pol0030_800,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20s %1s %20s %1s %20s %1s %20s %1s %20s %1s %20s %1s\n','|','Elevacao','|','Alcance(Inerte)','|','Alcance(BB 2pol)','|','Aumento(%)','|','Apogeu(Inerte)','|','Apogeu(BB 2po)','|','Aumento(%)','|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','711 mil (30 g/s)','|',alcin_711,'|',alcbb2pol0030_711,'|',aumnt2pol0030_711,'|',apgin_711,'|',apgbb2pol0030_711,'|',amtapg2pol0030_711,'|');2
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','711 mil (60 g/s)','|',alcin_711,'|',alcbb2pol0060_711,'|',aumnt2pol0060_711,'|',apgin_711,'|',apgbb2pol0060_711,'|',amtapg2pol0060_711,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','800 mil (30 g/s)','|',alcin_800,'|',alcbb2pol0030_800,'|',aumnt2pol0030_800,'|',apgin_800,'|',apgbb2pol0030_800,'|',amtapg2pol0030_800,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%1s %10s %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s %20.1f %1s\n','|','800 mil (60 g/s)','|',alcin_800,'|',alcbb2pol0060_800,'|',aumnt2pol0060_800,'|',apgin_800,'|',apgbb2pol0060_800,'|',amtapg2pol0060_800,'|');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------\n');    
toc

%% FUNCOES
%----- MODELO DE TRAJETORIA MASSA-PONTO MODIFICADO (INERTE)
function TRAJ = MTMPM_IN(PROJ)
    
    g0 RT omega P0 T0 R K B Pb d m m1 Ix Ix1 Az lat QD QM tab_coef_in rho acc

    Vel = [PROJ(1), PROJ(2), PROJ(3)];
    Pos = [PROJ(4), PROJ(5), PROJ(6)];
    Spn =  PROJ(7);    

    %--- modulos de velocidade e guinada de repouso
    v = norm(Vel);
    Ix = Ix1;
    m = m1;

    %--- modelo de atmosfera
    if (Pos(2) <= 1.1e4)
        Temp = T0 + B*Pos(2);
        Pres = P0*(1 + (B*Pos(2))/T0)^(-g0/(B*R));
    else
        Temp = 216.65;
        Pres = Pb*exp((-g0*(Pos(2)-11000))/(R*Temp));
    end
    rho = Pres/(R*Temp);
    
    %--- coeficientes aerodinamicos
    Ma = v/sqrt(K*R*Temp);

    CD0 = interp1(tab_coef_in.mach,tab_coef_in.cd0,Ma);
    CD2 = interp1(tab_coef_in.mach,tab_coef_in.cd2,Ma);
    CLa = interp1(tab_coef_in.mach,tab_coef_in.cla,Ma);
    CMa = interp1(tab_coef_in.mach,tab_coef_in.cma,Ma);
    Cmf = interp1(tab_coef_in.mach,tab_coef_in.cmag_f,Ma);
    Csp = interp1(tab_coef_in.mach,tab_coef_in.cxb,Ma);

    %--- guinada de repouso
    Yaw = yor(Vel, Spn, CMa);
    y = norm(Yaw);

    %--- forcas
    w = [omega*cos(lat)*cos(Az), omega*sin(lat), -omega*cos(lat)*sin(Az)];
    CF = -2*cross(w,Vel);

    DF = -((pi*rho*(d^2))/(8*m))*(CD0+CD2*((QD*y)^2))*v*Vel;
    
    g = g0*(1 - 0.0026*cos(2*lat));
    GF = -g*[Pos(1)/RT, 1-2*Pos(2)/RT, Pos(3)/RT];

    LF = (pi*rho*(d^2)/(8*m))*CLa*(v^2)*Yaw;

    MF = -((pi*rho*(d^3)*Spn*Cmf*QM)/(8*m))*cross(Yaw,Vel);
    
    %--- aceleracao, guinada de repouso e aceleracao de spin

    acc = CF+DF+GF+LF+MF;
    vel = Vel;
    spr = (pi*rho*(d^4)*Spn*v*Csp)/(8*Ix);

    TRAJ = [acc, vel, spr];

end

%----- MODELO DE TRAJETORIA MASSA-PONTO MODIFICADO (BASE BLEED)
function TRAJ = MTMPM_BB_1pol_0015(PROJ)
    
    g0 RT omega P0 T0 R K B Pb d db m0 dm1_0015 Ix0 CG0 m1 Ix1 CG1 Az lat QD QM rho acc tab_coef_bb_1pol_0015
  
    Vel = [PROJ(1), PROJ(2), PROJ(3)];
    Pos = [PROJ(4), PROJ(5), PROJ(6)];
    Spn =  PROJ(7); 
    Pmd = [PROJ(8), PROJ(9), PROJ(10)];

    %--- modulos de velocidade e guinada de repouso
    v = norm(Vel);

    %--- caracteristicas do projetil
    if Pmd(1) > m1
        m = Pmd(1);
        Ix = Pmd(2);
        CG = Pmd(3);
    else
        m = m1;
        Ix = Ix1;
        CG = CG1;
    end

    %--- modelo de atmosfera
    if (Pos(2) <= 1.1e4)
        Temp = T0 + B*Pos(2);
        Pres = P0*(1 + (B*Pos(2))/T0)^(-g0/(B*R));
    else
        Temp = 216.65;
        Pres = Pb*exp((-g0*(Pos(2)-11000))/(R*Temp));
    end
    rho = Pres/(R*Temp);
    
    %--- coeficientes aerodinamicos
    Ma = v/sqrt(K*R*Temp);

    CD0 = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cd0,Ma);
    CD2 = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cd2,Ma);
    CLa = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cla,Ma);
    CMa = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cma,Ma);
    Cmf = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cmag_f,Ma);
    Csp = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cspin,Ma);
    CDb = interp1(tab_coef_bb_1pol_0015.mach,tab_coef_bb_1pol_0015.cxb,Ma);

    %--- guinada de repouso
    Yaw = yor(Vel, Spn, CMa);
    y = norm(Yaw);

    %--- forcas
    w = [omega*cos(lat)*cos(Az), omega*sin(lat), -omega*cos(lat)*sin(Az)];
    CF = -2*cross(w,Vel);

    DF = -((pi*rho*(d^2))/(8*m))*(CD0+CD2*((QD*y)^2))*v*Vel;
    
    g = g0*(1 - 0.0026*cos(2*lat));
    GF = -g*[Pos(1)/RT, 1-2*Pos(2)/RT, Pos(3)/RT];

    LF = (pi*rho*(d^2)/(8*m))*CLa*(v^2)*Yaw;

    MF = -((pi*rho*(d^3)*Spn*Cmf*QM)/(8*m))*cross(Yaw,Vel);

    dm1pol = dm1_0015;

    inj = par_inj(rho, v, db, dm1pol);
    if Pmd(1) > m1
        BBF = (pi*rho*((d*v)^2)*CDb*inj)*(Vel*cos(y)/v + Yaw)/(8*m);
    else
        BBF = 0;
    end

    %--- aceleracao, guinada de repouso e aceleracao de spin
    acc = CF+DF+GF+LF+MF+BBF;
    vel = Vel;
    spr = (pi*rho*(d^4)*Spn*v*Csp)/(8*Ix);
    mrt = [dm1_0015, (Ix0-Ix1)*(m-m0)/(m0-m1), (CG0-CG1)*(m-m0)/(m0-m1)];

    TRAJ = [acc, vel, spr, mrt];

end

function TRAJ = MTMPM_BB_1pol_0030(PROJ)
    
    g0 RT omega P0 T0 R K B Pb d db m0 dm1_0030 Ix0 CG0 m1 Ix1 CG1 Az lat QD QM rho acc tab_coef_bb_1pol_0030
  
    Vel = [PROJ(1), PROJ(2), PROJ(3)];
    Pos = [PROJ(4), PROJ(5), PROJ(6)];
    Spn =  PROJ(7); 
    Pmd = [PROJ(8), PROJ(9), PROJ(10)];

    %--- modulos de velocidade e guinada de repouso
    v = norm(Vel);

    %--- caracteristicas do projetil
    if Pmd(1) > m1
        m = Pmd(1);
        Ix = Pmd(2);
        CG = Pmd(3);
    else
        m = m1;
        Ix = Ix1;
        CG = CG1;
    end

    %--- modelo de atmosfera
    if (Pos(2) <= 1.1e4)
        Temp = T0 + B*Pos(2);
        Pres = P0*(1 + (B*Pos(2))/T0)^(-g0/(B*R));
    else
        Temp = 216.65;
        Pres = Pb*exp((-g0*(Pos(2)-11000))/(R*Temp));
    end
    rho = Pres/(R*Temp);
    
    %--- coeficientes aerodinamicos
    Ma = v/sqrt(K*R*Temp);

    CD0 = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cd0,Ma);
    CD2 = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cd2,Ma);
    CLa = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cla,Ma);
    CMa = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cma,Ma);
    Cmf = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cmag_f,Ma);
    Csp = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cspin,Ma);
    CDb = interp1(tab_coef_bb_1pol_0030.mach,tab_coef_bb_1pol_0030.cxb,Ma);

    %--- guinada de repouso
    Yaw = yor(Vel, Spn, CMa);
    y = norm(Yaw);

    %--- forcas
    w = [omega*cos(lat)*cos(Az), omega*sin(lat), -omega*cos(lat)*sin(Az)];
    CF = -2*cross(w,Vel);

    DF = -((pi*rho*(d^2))/(8*m))*(CD0+CD2*((QD*y)^2))*v*Vel;
    
    g = g0*(1 - 0.0026*cos(2*lat));
    GF = -g*[Pos(1)/RT, 1-2*Pos(2)/RT, Pos(3)/RT];

    LF = (pi*rho*(d^2)/(8*m))*CLa*(v^2)*Yaw;

    MF = -((pi*rho*(d^3)*Spn*Cmf*QM)/(8*m))*cross(Yaw,Vel);

    dm1pol = dm1_0030;

    inj = par_inj(rho, v, db, dm1pol);
    if Pmd(1) > m1
        BBF = (pi*rho*((d*v)^2)*CDb*inj)*(Vel*cos(y)/v + Yaw)/(8*m);
    else
        BBF = 0;
    end

    %--- aceleracao, guinada de repouso e aceleracao de spin
    acc = CF+DF+GF+LF+MF+BBF;
    vel = Vel;
    spr = (pi*rho*(d^4)*Spn*v*Csp)/(8*Ix);
    mrt = [dm1_0030, (Ix0-Ix1)*(m-m0)/(m0-m1), (CG0-CG1)*(m-m0)/(m0-m1)];

    TRAJ = [acc, vel, spr, mrt];

end

%----- MODELO DE TRAJETORIA MASSA-PONTO MODIFICADO (BASE BLEED)
function TRAJ = MTMPM_BB_2pol_0030(PROJ)

    g0 RT omega P0 T0 R K B Pb d db m0 dm2_0030 Ix0 CG0 m1 Ix1 CG1 Az lat QD QM rho acc tab_coef_bb_2pol_0030
  
    Vel = [PROJ(1), PROJ(2), PROJ(3)];
    Pos = [PROJ(4), PROJ(5), PROJ(6)];
    Spn =  PROJ(7); 
    Pmd = [PROJ(8), PROJ(9), PROJ(10)];

    %--- modulos de velocidade e guinada de repouso
    v = norm(Vel);

    %--- caracteristicas do projetil
    if Pmd(1) > m1
        m = Pmd(1);
        Ix = Pmd(2);
        CG = Pmd(3);
    else
        m = m1;
        Ix = Ix1;
        CG = CG1;
    end

    %--- modelo de atmosfera
    if (Pos(2) <= 1.1e4)
        Temp = T0 + B*Pos(2);
        Pres = P0*(1 + (B*Pos(2))/T0)^(-g0/(B*R));
    else
        Temp = 216.65;
        Pres = Pb*exp((-g0*(Pos(2)-11000))/(R*Temp));
    end
    rho = Pres/(R*Temp);
    
    %--- coeficientes aerodinamicos
    Ma = v/sqrt(K*R*Temp);

    CD0 = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cd0,Ma);
    CD2 = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cd2,Ma);
    CLa = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cla,Ma);
    CMa = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cma,Ma);
    Cmf = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cmag_f,Ma);
    Csp = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cspin,Ma);
    CDb = interp1(tab_coef_bb_2pol_0030.mach,tab_coef_bb_2pol_0030.cxb,Ma);

    %--- guinada de repouso
    Yaw = yor(Vel, Spn, CMa);
    y = norm(Yaw);

    %--- forcas
    w = [omega*cos(lat)*cos(Az), omega*sin(lat), -omega*cos(lat)*sin(Az)];
    CF = -2*cross(w,Vel);

    DF = -((pi*rho*(d^2))/(8*m))*(CD0+CD2*((QD*y)^2))*v*Vel;
    
    g = g0*(1 - 0.0026*cos(2*lat));
    GF = -g*[Pos(1)/RT, 1-2*Pos(2)/RT, Pos(3)/RT];

    LF = (pi*rho*(d^2)/(8*m))*CLa*(v^2)*Yaw;

    MF = -((pi*rho*(d^3)*Spn*Cmf*QM)/(8*m))*cross(Yaw,Vel);

    dm2pol = dm2_0030;

    inj = par_inj(rho, v, db, dm2pol);
    if Pmd(1) > m1
        BBF = (pi*rho*((d*v)^2)*CDb*inj)*(Vel*cos(y)/v + Yaw)/(8*m);
    else
        BBF = 0;
    end

    %--- aceleracao, guinada de repouso e aceleracao de spin
    acc = CF+DF+GF+LF+MF+BBF;
    vel = Vel;
    spr = (pi*rho*(d^4)*Spn*v*Csp)/(8*Ix);
    mrt = [dm2_0030, (Ix0-Ix1)*(m-m0)/(m0-m1), (CG0-CG1)*(m-m0)/(m0-m1)];

    TRAJ = [acc, vel, spr, mrt];

end

%----- MODELO DE TRAJETORIA MASSA-PONTO MODIFICADO (BASE BLEED)
function TRAJ = MTMPM_BB_2pol_0060(PROJ)

    g0 RT omega P0 T0 R K B Pb d db m0 dm2_0060 Ix0 CG0 m1 Ix1 CG1 Az lat QD QM rho acc tab_coef_bb_2pol_0060
  
    Vel = [PROJ(1), PROJ(2), PROJ(3)];
    Pos = [PROJ(4), PROJ(5), PROJ(6)];
    Spn =  PROJ(7); 
    Pmd = [PROJ(8), PROJ(9), PROJ(10)];

    %--- modulos de velocidade e guinada de repouso
    v = norm(Vel);

    %--- caracteristicas do projetil
    if Pmd(1) > m1
        m = Pmd(1);
        Ix = Pmd(2);
        CG = Pmd(3);
    else
        m = m1;
        Ix = Ix1;
        CG = CG1;
    end

    %--- modelo de atmosfera
    if (Pos(2) <= 1.1e4)
        Temp = T0 + B*Pos(2);
        Pres = P0*(1 + (B*Pos(2))/T0)^(-g0/(B*R));
    else
        Temp = 216.65;
        Pres = Pb*exp((-g0*(Pos(2)-11000))/(R*Temp));
    end
    rho = Pres/(R*Temp);
    
    %--- coeficientes aerodinamicos
    Ma = v/sqrt(K*R*Temp);

    CD0 = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cd0,Ma);
    CD2 = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cd2,Ma);
    CLa = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cla,Ma);
    CMa = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cma,Ma);
    Cmf = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cmag_f,Ma);
    Csp = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cspin,Ma);
    CDb = interp1(tab_coef_bb_2pol_0060.mach,tab_coef_bb_2pol_0060.cxb,Ma);

    %--- guinada de repouso
    Yaw = yor(Vel, Spn, CMa);
    y = norm(Yaw);

    %--- forcas
    w = [omega*cos(lat)*cos(Az), omega*sin(lat), -omega*cos(lat)*sin(Az)];
    CF = -2*cross(w,Vel);

    DF = -((pi*rho*(d^2))/(8*m))*(CD0+CD2*((QD*y)^2))*v*Vel;
    
    g = g0*(1 - 0.0026*cos(2*lat));
    GF = -g*[Pos(1)/RT, 1-2*Pos(2)/RT, Pos(3)/RT];

    LF = (pi*rho*(d^2)/(8*m))*CLa*(v^2)*Yaw;

    MF = -((pi*rho*(d^3)*Spn*Cmf*QM)/(8*m))*cross(Yaw,Vel);

    dm2pol = dm2_0060;

    inj = par_inj(rho, v, db, dm2pol);
    if Pmd(1) > m1
        BBF = (pi*rho*((d*v)^2)*CDb*inj)*(Vel*cos(y)/v + Yaw)/(8*m);
    else
        BBF = 0;
    end

    %--- aceleracao, guinada de repouso e aceleracao de spin
    acc = CF+DF+GF+LF+MF+BBF;
    vel = Vel;
    spr = (pi*rho*(d^4)*Spn*v*Csp)/(8*Ix);
    mrt = [dm2_0060, (Ix0-Ix1)*(m-m0)/(m0-m1), (CG0-CG1)*(m-m0)/(m0-m1)];

    TRAJ = [acc, vel, spr, mrt];

end

%----- GUINADA DE REPOUSO
function Yaw = yor(Vel, Spn, CMa)

    Ix d rho acc
    
    v = norm(Vel);

    Yaw = -((8*(Ix*Spn))/(pi*rho*(d^3)*CMa*(v^4)))*cross(Vel,acc);

end

%--- PARÂMETRO DE INJEÇÃO
function inj = par_inj(rho, vel, db, dm)
    
    I0 = 0.005;

    I = (4*dm)/(pi*rho*vel*(db^2));
    
    if I < I0
        inj = I/I0;
    else
        inj = 1;
    end

end


%----- RUNGE-KUTTA (INERTE)
function [VEL, POS, SPN, T] = RK4_IN(Vel, Pos, Spn)
    
    %--- vetor caracteristicas do projetil
    proj = [Vel, Pos, Spn];
    
    %--- iteracoes e tempo
    i = 1;
    
    t = 0;
    dt = 0.006;
    
    while(proj(5) >= 0)
        
        %calculo das variacoes
        K1 = MTMPM_IN(proj);
        K2 = MTMPM_IN(proj + K1*dt/2);
        K3 = MTMPM_IN(proj + K2*dt/2);
        K4 = MTMPM_IN(proj + K3*dt);
    
        %solucao numerica pela funcao incremental
        proj = proj + (K1 + 2*K2 + 2*K3 + K4)*dt/6;
            
        %caracteristicas do projetil
        VEL(i,:) = [proj(1), proj(2), proj(3)];
        POS(i,:) = [proj(4), proj(5), proj(6)];
        SPN(i,:) =  proj(7);
        
        T(i,:) = t;
    
        %--- proximo passo
        t = t+dt;
        i = i+1;
    
    end
    
end

%----- RUNGE-KUTTA (BASE BLEED)
function [VEL, POS, SPN, PMD, T] = RK4_BB_1pol_0015(Vel, Pos, Spn, Pmd)

    %--- vetor caracteristicas do projetil
    proj = [Vel, Pos, Spn, Pmd];

    %--- iteracoes e tempo
    i = 1;
    
    t = 0;
    dt = 0.006;
    
    while(proj(5) >= 0)
        
        %calculo das variacoes
        K1 = MTMPM_BB_1pol_0015(proj);
        K2 = MTMPM_BB_1pol_0015(proj + K1*dt/2);
        K3 = MTMPM_BB_1pol_0015(proj + K2*dt/2);
        K4 = MTMPM_BB_1pol_0015(proj + K3*dt);
    
        %solucao numerica pela funcao incremental
        proj = proj + (K1 + 2*K2 + 2*K3 + K4)*dt/6;
    
        %caracteristicas do projetil
        V = [proj(1), proj(2), proj(3)];
        P = [proj(4), proj(5), proj(6)];
        S =  proj(7);
        M = [proj(8), proj(9), proj(10)];
    
        %--- armazenamento de dados
        VEL(i,:) = V;   
        POS(i,:) = P;
        SPN(i,:) = S;
        PMD(i,:) = M;

        T(i,:) = t;
    
        %--- proximo passo
        t = t+dt;
        i = i+1;
    
    end
    
end

%----- RUNGE-KUTTA (BASE BLEED)
function [VEL, POS, SPN, PMD, T] = RK4_BB_1pol_0030(Vel, Pos, Spn, Pmd)

    %--- vetor caracteristicas do projetil
    proj = [Vel, Pos, Spn, Pmd];

    %--- iteracoes e tempo
    i = 1;
    
    t = 0;
    dt = 0.006;
    
    while(proj(5) >= 0)
        
        %calculo das variacoes
        K1 = MTMPM_BB_1pol_0030(proj);
        K2 = MTMPM_BB_1pol_0030(proj + K1*dt/2);
        K3 = MTMPM_BB_1pol_0030(proj + K2*dt/2);
        K4 = MTMPM_BB_1pol_0030(proj + K3*dt);
    
        %solucao numerica pela funcao incremental
        proj = proj + (K1 + 2*K2 + 2*K3 + K4)*dt/6;
    
        %caracteristicas do projetil
        V = [proj(1), proj(2), proj(3)];
        P = [proj(4), proj(5), proj(6)];
        S = proj(7);
        M = [proj(8), proj(9), proj(10)];
    
        %--- armazenamento de dados
        VEL(i,:) = V;   
        POS(i,:) = P;
        SPN(i,:) = S;
        PMD(i,:) = M;

        T(i,:) = t;
    
        %--- proximo passo
        t = t+dt;
        i = i+1;
    
    end
    
end

%----- RUNGE-KUTTA (BASE BLEED)
function [VEL, POS, SPN, PMD, T] = RK4_BB_2pol_0030(Vel, Pos, Spn, Pmd)

    %--- vetor caracteristicas do projetil
    proj = [Vel, Pos, Spn, Pmd];

    %--- iteracoes e tempo
    i = 1;
    
    t = 0;
    dt = 0.006;
    
    while(proj(5) >= 0)
        
        %calculo das variacoes
        K1 = MTMPM_BB_2pol_0030(proj);
        K2 = MTMPM_BB_2pol_0030(proj + K1*dt/2);
        K3 = MTMPM_BB_2pol_0030(proj + K2*dt/2);
        K4 = MTMPM_BB_2pol_0030(proj + K3*dt);
    
        %solucao numerica pela funcao incremental
        proj = proj + (K1 + 2*K2 + 2*K3 + K4)*dt/6;
    
        %caracteristicas do projetil
        V = [proj(1), proj(2), proj(3)];
        P = [proj(4), proj(5), proj(6)];
        S =  proj(7);
        M = [proj(8), proj(9), proj(10)];
    
        %--- armazenamento de dados
        VEL(i,:) = V;   
        POS(i,:) = P;
        SPN(i,:) = S;
        PMD(i,:) = M;

        T(i,:) = t;
    
        %--- proximo passo
        t = t+dt;
        i = i+1;
    
    end
    
end

function [VEL, POS, SPN, PMD, T] = RK4_BB_2pol_0060(Vel, Pos, Spn, Pmd)

    %--- vetor caracteristicas do projetil
    proj = [Vel, Pos, Spn, Pmd];

    %--- iteracoes e tempo
    i = 1;
    
    t = 0;
    dt = 0.006;
    
    while(proj(5) >= 0)
        
        %calculo das variacoes
        K1 = MTMPM_BB_2pol_0060(proj);
        K2 = MTMPM_BB_2pol_0060(proj + K1*dt/2);
        K3 = MTMPM_BB_2pol_0060(proj + K2*dt/2);
        K4 = MTMPM_BB_2pol_0060(proj + K3*dt);
    
        %solucao numerica pela funcao incremental
        proj = proj + (K1 + 2*K2 + 2*K3 + K4)*dt/6;
    
        %caracteristicas do projetil
        V = [proj(1), proj(2), proj(3)];
        P = [proj(4), proj(5), proj(6)];
        S =  proj(7);
        M = [proj(8), proj(9), proj(10)];
    
        %--- armazenamento de dados
        VEL(i,:) = V;   
        POS(i,:) = P;
        SPN(i,:) = S;
        PMD(i,:) = M;

        T(i,:) = t;
    
        %--- proximo passo
        t = t+dt;
        i = i+1;
    
    end
    
end