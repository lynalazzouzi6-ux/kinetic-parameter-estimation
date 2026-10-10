function r = global_residual(p, Tvec, tcell, Ccell, CA0, R, Tref)
% p(1) = ln(k_ref), p(2) = Ea en kJ/mol
% Résidu en log car le bruit est relatif
r = [];
for i = 1:numel(Tvec)
    k = exp(p(1)) * exp(-p(2)*1e3/R * (1/Tvec(i) - 1/Tref));
    C = model_first_order(k, tcell{i}, CA0);
    r = [r; log(C) - log(Ccell{i})];
end
end