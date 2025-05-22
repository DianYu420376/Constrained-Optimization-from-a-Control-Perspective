clc
clear all
% Solve the OPF problem
define_constants;
mpc = case39;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
     MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
 [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc);

Va   = bus(:, VA) * (pi/180);
global dum_idx
dum_idx = find(Va == 0);
Vm   = bus(:, VM);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;


mpopt = mpoption();%mpoption('opf.use_vg',1);
%mpopt = mpoption('opf.flow_lim','P');

om = opf_setup(mpc, mpopt); %create OPF model

x = [Va; Vm; Pg; Qg]; % Notice that in Va, there's a dum variable 0, so also in our setting, we are not optimizing the Pg of


Nvar = om.getN("var");
Nnli = om.getN("nli");
Nnle = om.getN("nle");

x0 = om.params_var;
x0(dum_idx) = [];
global Al Au idxl idxu vl vu
vl = om.var.params.vl;
vu = om.var.params.vu;
idxl = (vl > -inf);
idxu = (vu < inf);
Al = -sparse(eye(Nvar));
Au = sparse(eye(Nvar));
Al = Al(idxl,:);
Au = Au(idxu,:);
Al(:,dum_idx) = [];
Au(:,dum_idx) = [];
Nl = size(Al,1);
Nu = size(Au,1)


Df1 = @(x) Df(om,x);
Dhi1 = @(x) Dhi(om,x);
Dhe1 = @(x) Dhe(om,x);
eta = 5;
k = 1/eta * eye(Nnli + Nnle + Nl + Nu);
max_iters = 2000;

f_lst = [];
df_lst = [];
x_lst = [];
x=x0;
% for iter = 1:max_iters
%     [f_val,grad_f] = om.eval_nln_constraint(x,0);
%     idx = 10;
%     f_val = f_val(idx);
%     grad_f = grad_f(idx,:)';
%     x = x- eta*grad_f;
%     f_lst = [f_lst, f_val];
%     df_lst = [df_lst, norm(grad_f)];
%     x_lst = [x_lst, x];
% end

[x, f_vals, hi_vals, he_vals, KKT_gaps] = FL(Df1, Dhi1, Dhe1, x0, eta, k, max_iters);
plot(log10(KKT_gaps));
save('case39_history')
function [f_x, grad_f] = Df(om, x)
global dum_idx
x1 = insertElement(x, 0, dum_idx);
[f_x, grad_f] = om.eval_quad_cost(x1);
grad_f(dum_idx) = [];
f_x = f_x/2000;
grad_f = grad_f/2000;
% [f_x,grad_f] = om.eval_quad_cost(x);
end

function [hi_x, grad_hi] = Dhi(om, x)
global dum_idx
x1 = insertElement(x, 0, dum_idx);
[hi1_x, grad_hi1] = om.eval_nln_constraint(x1,0);
grad_hi1(:,dum_idx) = [];
global Al Au idxl idxu vl vu
hi2_x = vl(idxl) - x1(idxl) ;
hi3_x = x1(idxu) - vu(idxu);
hi_x = [hi1_x;hi2_x;hi3_x];
grad_hi = [grad_hi1;Al;Au];
end

function [he_x, grad_he] = Dhe(om, x)
global dum_idx
x1 = insertElement(x, 0, dum_idx);
[he_x, grad_he] = om.eval_nln_constraint(x1,1);
grad_he(:,dum_idx) = [];
end