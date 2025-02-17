%% ----------passaggi---------------
% 1 - Caricare input_atmosfera.mat e file hitran
% 2 - Calcolare l'altezza degli strati
% 3 - Calcolare la densita colonnare dei gas
% 4 - Calcolare il trasferimento radiativo

%% ------1------
% Costanti fisiche di riferimento
c = 2.998e8; % Velocità della luce (m/s)
h = 6.626e-34; % Costante di Planck (J·s)
k = 1.381e-23; % Costante di Boltzmann (J/K)
g0 = 9.80665; % Accelerazione di gravita'
Na = 6.022142e23; % Numero di Avogadro (mol^-1)
R = k * Na; % Costante dei gas (J/(mol·K))
Re = 6371; % Raggio terra
mair = 28.964; % Massa molare aria

% Caricamento dell'atmosfera e dei dati spettrali
load('input_atmosfera.mat'); % carica p, t e mr

% Definizione dell'intervallo di frequenze (numeri d'onda)
w_min = 900;
w_max = 1100;
w = (w_min:0.01:w_max)';

%% ------2------
% Calcolo dell'altezza degli strati
Z = calcola_altezze_strati(p,t);

%% ------3------
% Calcolo delle densita' colonnari dei gas
% (H2O, CO2, O3, N2O, CO, CH4, SO2, HNO3, NH3, HDO, OCS, CF4).

vars = whos; % Ottiene tutte le variabili nel workspace
for i = 1:length(vars)
    var_name = vars(i).name;
    % Controlla se la variabile termina con '_mr'
    if endsWith(var_name, '_mr')
        gas_mr = eval(var_name); % Ottiene il valore della variabile
        new_var_name = strrep(var_name, '_mr', '_col'); % Sostituisce '_mr' con '_col'
        % Calcola la densità colonnare e assegna alla nuova variabile
        assignin('base', new_var_name, calcola_gas_col(t, p, gas_mr, Z));
    end
end
% Calcolo della densità colonnare per l'acqua (con specifica massa molecolare)
H2O_col = calcola_gas_col(t, p, H2O_mr, Z, 18.01);

%% -------4-------
% Caricamento dei dati spettrali
h2o  = load('H2O_hitran.txt');
co2  = load('CO2_hitran.txt');
o3   = load('O3_hitran.txt');
nh3  = load('NH3_hitran.txt');
hno3 = load('HNO3_hitran.txt');

% Calcolo delle pressioni e temperature medie di ogni strato
[pm, tm] = calcola_medie(p, t);
n_layers = length(pm);  % Numero di strati
n_w = length(w);

% Parametri angolari e superficie
angolo_zenitale = 45 * pi/180;      % 45 gradi in radianti
angolo_osservazione = 0 * pi/180;   % nadir
mu_s = cos(angolo_zenitale);        % coseno angolo solare
mu = cos(angolo_osservazione);      % coseno angolo di osservazione
emissivita_superficie = 0.7;        % emissività superficiale
temperatura_superficie = 300;       % temperatura superficiale in K

% Inizializzazione della matrice delle trasmittanze
tau = ones(n_w, n_layers + 1);
tau_mu = ones(n_w, n_layers + 1);  % trasmittanza per il percorso solare

% Calcolo dello spessore ottico dalla parte alta dell'atmosfera verso il basso
for l = n_layers:-1:1
    TT = tm(l); % Temperatura dello strato l-esimo
    PP = pm(l); % Pressione dello strato l-esimo
    
    % Calcolo dei coefficienti di assorbimento (k) per i gas
    H2O  = lorentzian(h2o(:,3), h2o(:,4), h2o(:,5), w, 0, TT, PP);
    CO2  = lorentzian(co2(:,3), co2(:,4), co2(:,5), w, 0, TT, PP);
    O3   = lorentzian(o3(:,3), o3(:,4), o3(:,5), w, 0, TT, PP);
    NH3  = lorentzian(nh3(:,3), nh3(:,4), nh3(:,5), w, 0, TT, PP);
    HNO3 = lorentzian(hno3(:,3), hno3(:,4), hno3(:,5), w, 0, TT, PP);

    % Calcolo dello spessore ottico complessivo (k*rho*dz) dovuto a tutte le righe spettrali
    opd = (sum(H2O, 1) .* H2O_col(l) + ...
           sum(CO2, 1) .* CO2_col(l) + ...
           sum(O3, 1) .* O3_col(l) + ...
           sum(NH3, 1) .* NH3_col(l) + ...
           sum(HNO3, 1) .* HNO3_col(l))';
    
    % Calcolo trasmittanze per vista nadir e percorso solare
    tau(:,l) = tau(:,l+1) .* exp(-opd);
    tau_mu(:,l) = tau_mu(:,l+1) .* exp(-opd/mu_s);
end

% Calcolo della radianza up-welling (combinando superficie ed atmosfera)
radianza_upwelling = emissivita_superficie * blackbodyn(temperatura_superficie, w) .* tau(:,1);
for l = n_layers:-1:1
    dtau = tau(:,l+1) - tau(:,l);  % variazione discreta della trasmittanza
    radianza_upwelling = radianza_upwelling + blackbodyn(tm(l), w) .* dtau;
end


% Calcolo della radianza down-welling atmosferica
radianza_downwelling = zeros(n_w, 1);
for l = 1:n_layers
    dtau_inv = 1./tau(:,l+1) - 1./tau(:,l);  % variazione discreta
    radianza_downwelling = radianza_downwelling + blackbodyn(tm(l), w) .* dtau_inv;
end
radianza_downwelling = (emissivita_superficie-1) * tau(:,1).^2 .*radianza_downwelling;

% Calcolo del contributo solare diretto
T_sole = 5800;
dim_ang_sole = 6.8e-5;
E_sigma = blackbodyn(T_sole, w) * dim_ang_sole;

radianza_solare = (1-emissivita_superficie)/pi * tau(:,1) .* ...
    tau_mu(:,1) .* mu_s .* E_sigma;

% Calcolo delle componenti riflesse
riflettanza_superficie = 1 - emissivita_superficie;

% Calcolo della radianza totale includendo tutte le componenti
radianza_totale = radianza_upwelling + ...  % Componente diretta up-welling
                  radianza_downwelling + ...  % Componente diretta down-welling
                  radianza_solare ;      % Componente solare diretta
                  

%% Visualizzazione delle radianze
figure;
subplot(3,1,1); 
plot(w, radianza_upwelling, 'b', ...
     w, radianza_downwelling, 'r', ...
     w, radianza_solare, 'g');
legend('Up-welling', 'Down-welling', 'Contributo Solare Diretto');
xlabel('Numero d''onda (cm^{-1})');
ylabel('Radianza');
title('Componenti Radiative');
grid on;

subplot(3,1,2);
plot(w, radianza_solare, 'c', ...
     w, radianza_downwelling, 'm');
legend('Solare Riflessa', 'Atmosferica Riflessa');
xlabel('Numero d''onda (cm^{-1})');
ylabel('Radianza');
title('Componenti Riflesse');
grid on;

subplot(3,1,3);
plot(w, radianza_totale, 'k', 'LineWidth', 1.5);
xlabel('Numero d''onda (cm^{-1})');
ylabel('Radianza');
title('Radianza Totale (Tutte le Componenti)');
grid on;

%% Visualizzazione dello spettro di assorbimento dei gas sovrapposti
figure;
hold on;
plot(w, H2O, 'b', 'LineWidth', 1.5, 'DisplayName', 'H2O');
plot(w, CO2, 'r', 'LineWidth', 1.5, 'DisplayName', 'CO2');
plot(w, O3, 'g', 'LineWidth', 1.5, 'DisplayName', 'O3');
plot(w, NH3, 'm', 'LineWidth', 1.5, 'DisplayName', 'NH3');
plot(w, HNO3, 'c', 'LineWidth', 1.5, 'DisplayName', 'HNO3');
title('Spettri di Assorbimento dei Gas');
xlabel('Numero d''onda (cm^{-1})');
ylabel('Coefficiente di Assorbimento');
legend('Location', 'best');
grid on;
hold off;

%% Grafico Pressione vs Temperatura per ciascun strato atmosferico
figure;
plot(tm, pm, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Temperatura media [K]');
ylabel('Pressione media [hPa]');
title('Pressione vs Temperatura per ciascun strato atmosferico');
grid on;
set(gca, 'YDir','reverse');  % In atmosfera la pressione diminuisce con l'altezza

%% Grafico del Mixing Ratio dei gas nell'atmosfera
figure;
semilogy(H2O_mr, pm, 'Color','b', 'LineWidth',1.5); hold on;
semilogy(CO2_mr, pm, 'Color','r', 'LineWidth',1.5);
semilogy(O3_mr, pm, 'Color','g', 'LineWidth',1.5);
semilogy(NH3_mr, pm, 'Color','m', 'LineWidth',1.5);
semilogy(HNO3_mr, pm, 'Color','c', 'LineWidth',1.5);
semilogy(CF4_mr, pm, 'Color','y', 'LineWidth',1.5);
semilogy(CH4_mr, pm, 'Color',[0.2 0.7 0.8], 'LineWidth',1.5);
semilogy(HDO_mr, pm, 'Color',[0.5 0 0.5], 'LineWidth',1.5);
semilogy(N2O_mr, pm, 'Color',[0.5 0.5 0], 'LineWidth',1.5);
semilogy(OCS_mr, pm, 'Color',[0 0.5 0.5], 'LineWidth',1.5);
semilogy(SO2_mr, pm, 'Color',[0.25 0.25 0.25], 'LineWidth',1.5);
semilogy(CO_mr, pm, 'Color',[0.8 0.8 0.8], 'LineWidth',1.5);

xlabel('Mixing Ratio');
ylabel('Pressione [hPa]');
title('Pressione vs Mixing Ratio per gas');

set(gca, 'YDir','reverse', 'XScale','log', 'YScale','log');
grid on;
legend('H_2O','CO_2','O_3','NH_3','HNO_3','CF_4','CH_4','HDO','N_2O','OCS','SO_2','CO',...
       'Location','best');

%% FUNCTIONS UTILIZZATE     
function [BB] = blackbodyn(T,numerionda)
    c = 2.998e8; % Velocità della luce
    h = 6.626e-34; % Costante di Planck
    k = 1.381e-23; % % Costante di Boltzmann
    c1 = 2 * h * c^2 * 1e8;
    c2 = h * c / k * 1e2;

    BB = c1 .* numerionda.^3;
    BB = BB .* exp(-c2 .* numerionda ./ T) ./ (1.0 - exp(-c2 .* numerionda ./ T));
end

%VERSIONE OTTIMIZZATA
function [f] = lorentzian(w0, S, a0, w, dwm, T, p)
    % Implementazione della funzione di Lorentz per il calcolo dello spettro di assorbimento
    % dei gas in atmosfera
    %
    % Parametri di input:
    % w0  - numeri d'onda centrali delle linee spettrali [cm^-1]
    % S   - intensità delle linee spettrali
    % a0  - parametri di allargamento delle linee a temperatura e pressione di riferimento
    % w   - array dei numeri d'onda su cui calcolare lo spettro [cm^-1]
    % dwm - (opzionale) distanza massima dal centro della riga da considerare
    % T   - (opzionale) temperatura [K]
    % p   - (opzionale) pressione [hPa]
    % 
    T0 = 296;  % Temperatura di riferimento [K]
    p0 = 1013.25; % Pressione di riferimento [hPa]
    
    % Valori di default
    if ~exist('dwm','var'), dwm = 0; end
    if ~exist('T','var'), T = T0; end
    if ~exist('p','var'), p = p0; end
    
    % Riorganizzazione degli array di input per garantire dimensioni corrette
    a0 = reshape(a0, length(a0), 1);
    S = reshape(S, length(S), 1) ./ pi;
    w0 = reshape(w0, length(w0), 1);
    w = reshape(w, 1, length(w));
    
    % Calcolo dei parametri di allargamento corretti per temperatura e pressione attuali
    a = (T0 / T).^0.5 * p / p0 .* a0;
    
    % Inizializzazione array di output
    f = zeros(1, length(w));
    
    % Elaborazione a blocchi per ottimizzare l'uso della memoria
    chunk_size = 1000;  % Elabora 1000 linee alla volta
    num_chunks = ceil(length(w0) / chunk_size);
    
    for i = 1:num_chunks
        % Calcolo degli indici per il blocco corrente
        idx_start = (i-1)*chunk_size + 1;
        idx_end = min(i*chunk_size, length(w0));
        current_indices = idx_start:idx_end;
        
        % Estrazione dei dati per il blocco corrente
        w0_chunk = w0(current_indices);
        S_chunk = S(current_indices);
        a_chunk = a(current_indices);
        
        % Calcolo della distanza dal centro per ogni riga
        DW = w - w0_chunk * ones(1, length(w));
        if dwm > 0
            mask = abs(DW) <= dwm;
        else
            mask = true(size(DW));
        end
        
        % Calcolo del profilo di Lorentz per il blocco corrente
        f_chunk = (S_chunk * ones(1, length(w))) .* (a_chunk * ones(1, length(w))) ...
                 ./ (DW.^2 + (a_chunk * ones(1, length(w))).^2);
        
        % Applicazione del cut-off se specificato
        f_chunk = f_chunk .* mask;
        
        % Somma del contributo del blocco corrente al totale
        f = f + sum(f_chunk, 1);
    end
end

function [Z] = calcola_altezze_strati(p,t)
    % COSTANTI
    g0 = 9.80665; % Accelerazione di gravita'
    Na = 6.022142e+23; % Numero di Avogadro
    Re = 6371; % Raggio Terra
    k = 1.38065e-23; % Costante di Boltzmann
    mair = 28.964; % Massa molare aria
    R = k * Na; % Costante dei gas
    
    % Pressioni e temperature superiori
    Pu = p(2:end) * 100; % conversione da hPa a Pa
    Tu = t(2:end);
    % Pressioni e temperature inferiori
    Pl = p(1:end-1) * 100; % conversione da hPa a Pa
    Tl = t(1:end-1);
    
    Z  = zeros(numel(p), 1); % altezza
    G = zeros(numel(p), 1); % gravità 
    Pm = zeros(numel(p)-1, 1); % pressione media dello strato
    Tm = zeros(numel(p)-1,1); % temperatura media
    Z(1) = 0.0; % altezza al livello del mare
    
    % Loop
    for i = 2:numel(p)
        tup = Tu(i-1); tlo = Tl(i-1);
        pup = Pu(i-1); plo = Pl(i-1); % indicizzati da 1 a n-1
        Tm(i-1) = tup + (tlo - tup)/log(plo/pup)*((plo*(log(plo)-1)-pup*(log(pup)-1))/(plo-pup)-log(pup));
        G(i-1) = g0 * Re^2 / (Re + Z(i-1))^2;
        % G(i-1) = g0; % se g costante, anche nella parte bassa dell'atm
        Z(i) = Z(i-1) - log(pup/plo) * R * Tm(i-1) / (mair * G(i-1)); % in Km
        Pm(i-1) = (pup - plo) / (log(pup/plo)); % in Pascal
    end
end 

function [pm, tm] = calcola_medie(p,t)
    n = length(p) - 1;
    pm = zeros(n, 1);
    tm = zeros(n, 1);

    for i = 1:n
        pm(i) = (p(i+1) + p(i)) / 2; % media delle pressioni
        tm(i) = (t(i+1) + t(i)) / 2; % media delle temperature
    end
end


function gas_col = calcola_gas_col(t, p, gas_mr, Z, Mgas)
    % Calcola la densità colonnare di un gas in atmosfera.
    %
    % INPUT:
    %   t      - vettore delle temperature (K) per ciascun livello
    %   p      - vettore delle pressioni (hPa) per ciascun livello
    %   gas_mr - vettore del mixing ratio; se Mgas è specificato, si assume
    %            che sia espresso in g/kg, altrimenti in frazione molare.
    %   Z      - vettore delle altezze (in metri) dei bordi degli strati
    %   Mgas   - (opzionale) massa molecolare del gas in g/mol.
    %
    % OUTPUT:
    %   gas_col - vettore della densità colonnare (molecole/m^2) per ciascuno strato.
    
    % Calcola le medie (per pressione e temperatura) per ciascuno strato.
    [pm, tm] = calcola_medie(p,t);
    
    % Conversione da hPa a Pa per la pressione media
    pm = pm * 100;   
    
    % Costanti fisiche
    R = 8.314;              % J/(mol·K)
    Na = 6.022e23;          % molecole/mol
    
    % Numero di livelli (strati = livelli - 1)
    n_levels = length(p);
    gas_col = zeros(n_levels-1, 1);
    
    % Se Mgas è specificato, si assume che gas_mr sia in g/kg e va convertito
    if nargin > 4
        % Conversione da g/kg a frazione molare
        gas_kg = gas_mr * 0.001;
        M_air = 28.97;  % g/mol
        epsilon = Mgas / M_air;
        chi = gas_kg ./ (gas_kg + epsilon);
    else
        % Si assume che gas_mr sia già in frazione molare
        chi = gas_mr;
    end

    % Loop su ogni strato per calcolare la densità colonnare
    for i = 1:n_levels-1
        % Spessore dello strato (in metri)
        dz = Z(i+1) - Z(i);
        % Numero di moli d'aria per unità di area (mol/m^2)
        moli_aria_per_m2 = (pm(i) / (R * tm(i))) * dz * 1000 * 1e-4;
        % Calcola la densità colonnare del gas
        gas_col(i) = chi(i) * moli_aria_per_m2 * Na;
    end
end