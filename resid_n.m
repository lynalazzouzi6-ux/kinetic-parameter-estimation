function r = resid_n(p, t, CA0s, data)
% Résidus (en log) pour l'ajustement de p = [k, n] sur plusieurs séries
%   p    : [k, n]
%   t    : vecteur colonne des temps
%   CA0s : vecteur des concentrations initiales (une par série)
%   data : cell array, data{i} = mesures de la série i
% Le log est utilisé car le bruit est relatif.
r = [];
for i = 1:numel(CA0s)
    C = model_order_n(p, t, CA0s(i));
    r = [r; log(max(C, 1e-8)) - log(data{i})];
end
end