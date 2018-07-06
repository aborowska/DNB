clear all;

arg1 = 0; % 1; % DA Kim
arg2 = 0; % 2; % DA RW eff
arg3 = 0; % semi DA
arg4 = 0; % semi DA shifting
arg5 = 0; % semi DA eff
arg6 = 0; % semi DA adaptive
arg7 = 1; % semi DA adaptive eff

for ii = [0,1,2,4]
% ii=0; % 0: simulation; 1-6 different data sets
    for bins = [15,25] %[10,20,30]
        try_SV_param_linux(ii,arg1,arg2,0,0,0,0,bins)
    end
end