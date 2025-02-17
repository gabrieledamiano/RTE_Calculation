function [BB] = blackbodyn(T,numerionda)
    c = 2.998e8; % Velocità della luce
    h = 6.626e-34; % Costante di Planck
    k = 1.381e-23; % % Costante di Boltzmann
    c1 = 2 * h * c^2 * 1e8;
    c2 = h * c / k * 1e2;

    BB = c1 .* numerionda.^3;
    BB = BB .* exp(-c2 .* numerionda ./ T) ./ (1.0 - exp(-c2 .* numerionda ./ T));
end