% generate zero-inlfated SV
gamma = 0.1;
mu = -1;
phi = 0.98;
sigma2 = 0.02;
sigma = sqrt(sigma2);

T = 10000;

h_true = zeros(T,1); 

h_true(1) = mu + sqrt(sigma2/(1-phi^2))*randn;
u = rand(T,1);

for ii = 2:T
    h_true(2) = mu + phi*(h_true(2) - mu) + sigma*randn;
end

y = exp(h_true/2).*randn(T,1);
y0 = y;
y0(u<=gamma) = 0;

