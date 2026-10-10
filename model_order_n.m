function CA = model_order_n(p, t, CA0)   % p = [k, n]
odefun = @(~,C) -p(1)*max(C,0)^p(2);
[~, CA] = ode45(odefun, t, CA0, odeset('RelTol',1e-8,'AbsTol',1e-10));
end