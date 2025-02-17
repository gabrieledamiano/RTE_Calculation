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