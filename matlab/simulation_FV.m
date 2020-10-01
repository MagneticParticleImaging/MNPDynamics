function [t, exp, y] = simulation_FV(B, t_vec, tr, varargin)

% Returns the solution of the Neel rotation model for given parameters
% using the Finite Volume method. For Brownian rotation, use the spherical
% harmonics method.
%       INPUTS:
%       B: function of one scalar variable (time) that returns a 3D column
%           vector, describes the magnetic field over time
%       t_vec: vector of time points where the solution is to be evaluated
%       tr: struct containing the information about the triangulation as
%           returned by the FV_mesh script.
%       (optional): struct containing parameters; if any are not provided,
%           standard values are used.
%          Possible parameters:
%           M_S: Core magnetization of a particle
%           D_core: Core diameter of a particle in meters
%           Temp: Temperature in K
%           kAnis: Neel anisotropy constant
%           tau_N: scalar, Neel relaxation time constant, usually defined as
%               M_S*V_C/(2*alphha*gamma_tilde*k_B*T)
%           alpha: damping coefficient, usually 0.1
%           n: Neel easy axis, to be specified as a 3D column vector or a 
%               function handle that returns a 3D column vector given a 1D 
%               point in time
%           beta: upwind interpolation parameters, 0 corresponds to central difference
%               averaging of the advection term, 1 corresponds to taking
%               the value of the upstream node, 0 < beta < 1 corresponds to
%               a weighted sum of the two.
%           p1 to p4: coefficients of the advection term. It takes the form
%               p1* (H x m) + p2 (m x H) x m + p4 (n*m)n x m + p4 (n*m)(m x n)
%               x m. These can be set explicitly; if not, they will be
%               calculated with the other given or assumed parameters. Setting
%               p1 and p3 to zero explicitly may speed up the computation
%               significantly with little to no error made.
%
%       OUTPUTS:
%       t: vector of time points where the solution was evaluated. Is equal
%          to input t_vec if the integration of the ODE succeeded.
%       exp: matrix of calculated mean magnetic moment.
%           Dimension: length(t) x 3.
%       y: complete probability distribution, given in terms of the mean
%           value of the distribution on each triangle.
%           Dimension: length(t) x (number of triangles) 




k_B = 1.38064852e-23;
gam_gyro=1.75*10^11;
if nargin > 3
    params = varargin{1};
else
    params.i = {};
end

if isfield(params,'M_S')
    M_S = params.M_S;
else
    M_S = 474000;
end
if isfield(params, 'D_core')
    D_core = params.D_core;
else
    D_core = 20e-9;
end
V_core = pi/6 * D_core^3;
if isfield(params,'Temp')
    Temp = params.Temp;
else
    Temp = 293;
end
if isfield(params, 'alpha')
    alpha = params.alpha;
else
    alpha = 0.1;
end
if isfield(params, 'kAnis')
    kAnis = params.kAnis;
else
    kAnis = 625;
end
if isfield(params,'tau_N')
    tau_N = params.tau_N;
else
    tau_N = M_S*V_core/(k_B*Temp*gam_gyro)*(1+alpha^2)/(2*alpha);
end
if isfield(params,'n')
    n = params.n;
else
    n = [0;0;1];
end
if isfield(params, 'RelTol')
    RelTol = params.RelTol;
else
    RelTol = 1e-3;
end
if isfield(params, 'beta')
    beta = params.beta;
else
    beta = 0;
end
if isfield(params, 'p1')
    p1 = params.p1;
else
    p1 = gam_gyro/(1+alpha^2);
end
if isfield(params,'p2')
    p2 = params.p2;
else
    p2 = alpha*gam_gyro/(1+alpha^2);
end
if isfield(params, 'p3')
    p3 = params.p3;
else
    p3 = 2*gam_gyro/(1+alpha^2)*kAnis/M_S;
end
if isfield(params, 'p4')
    p4 = params.p4;
else
    p4 = 2*alpha*gam_gyro/(1+alpha^2)*kAnis/M_S;
end

N = size(tr.fMat,1);
tr.C = 1/(2*tau_N)*tr.areasi.*tr.C;
y0 = 1/(4*pi)*ones(1,N);
rhs = @(t,y) FV_matrix(t, B, N, p1, p2, p3, p4, beta, tr.mids, tr.ds, n, tr.iis, tr.valcs, tr.C, tr.C, tr.a_ijs, tr.e_is, tr.areasidil, tr.tr2edge, tr.flow_signs)*y;
jac = @(t,y) FV_matrix(t, B, N, p1, p2, p3, p4, beta, tr.mids, tr.ds, n, tr.iis, tr.valcs, tr.C, tr.C, tr.a_ijs, tr.e_is, tr.areasidil, tr.tr2edge, tr.flow_signs);
opts = odeset('Jacobian', jac, 'RelTol', RelTol);
[t,y] = ode15s(rhs, t_vec, y0, opts);

exp = zeros(length(t), 3);
for i=1:N
    exp = exp + y(:,i).*tr.areas(i).*tr.centers(i,:);
end


end
