function [pm, tm] = calcola_medie(p,t)
    pm = zeros(60, 1);
    tm = zeros(60, 1);

    for i = 1:60
        pm(i) = (p(i+1) + p(i)) / 2; % Media pressioni
        tm(i) = (t(i+1) + t(i)) / 2; % Media temperature
    end
end