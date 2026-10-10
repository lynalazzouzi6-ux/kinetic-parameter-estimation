function CA = model_first_order(k, t, CA0)
% Résout dCA/dt = -k*CA et renvoie CA aux instants t (vecteur colonne)
odefun = @(~, C) -k*C;
opts   = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~, CA] = ode45(odefun, t, CA0, opts);
end