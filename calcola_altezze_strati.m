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