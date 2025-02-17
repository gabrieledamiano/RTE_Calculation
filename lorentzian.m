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