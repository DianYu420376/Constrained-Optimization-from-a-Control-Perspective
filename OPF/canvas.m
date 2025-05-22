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
Vm   = bus(:, VM);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;


mpopt = mpoption();%mpoption('opf.use_vg',1);
mpopt = mpoption('opf.flow_lim','P');

om = opf_setup(mpc, mpopt); %create OPF model

x = [Va; Vm; Pg; Qg]; % Notice that in Va, there's a dum variable 0, so also in our setting, we are not optimizing the Pg of

% Okay, some thoughts on what needs to be done
% 1. Given x, find corresponding index for dummy varaible, optimization
% variable (OV) and dependent variable (DV)
% 2. Call run_pf to solve for the dependent variable
% 3. Calculate jacobian DV/OV (using implicit function and eval_nln_constraint(x,1))
% 4. Calculate the jacobian of the inequality constraints (using eval_nln_constraint(x,0) and chain rule)
% 5. Calculate the gradient of the quadratic cost
% 6. SQP problem


Nvar = om.getN("var");
Nnli = om.getN("nli");
Nnle = om.getN("nle");

[G, DG]=om.eval_nln_constraint(x,0);


%% Runyu Note
% "help opt_model" might be helpful

% https://matpower.app/manual/matpower/matpowerFunctions.html

%% Copied from PF_function2.m
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;
gencost = mpc.gencost;

Pload = bus(:, PD); 
Qload = bus(:, QD);
BusType = bus(:, BUS_TYPE);
BusPV_Ref = BusType(BusType > 1);
BusPV = find(BusPV_Ref == PV);
BusRef = find(BusPV_Ref == REF);

Pgen = gen(:, PG);
Qgen = gen(:, QG);
Pgen_max = gen(:, PG) * 1.5;
Pgen_min = gen(:, PG) * 0.5;
Qgen_max = gen(:, QMAX);
Qgen_min = gen(:, QMIN);

Prating = branch(:, RATE_A) / 1.5;

% 设置发电机输出
n1 = length(mpc.gen(BusPV, PG));
n2 = length(mpc.gen(:,VG));
assert(length(Action) == n1+n2)
Action1 = Action(1:n1);
Action2 = Action(n1+1:end);

temp = Pgen_min(BusPV) + (Pgen_max(BusPV) - Pgen_min(BusPV)) .* Action1;
temp2 = 0.95 + 0.1*Action2;
mpc.gen(BusPV, PG) = temp;
mpc.gen(:,VG) = temp2;

% 设置潮流计算选项
mpopt = mpoption('PF_ALG', 1, 'PF_TOL', 1e-4, 'PF_MAX_IT', 100, 'OUT_ALL', 0);

% 运行潮流计算
results = runpf(mpc, mpopt);
res_success = results.success;  % if success == 1, power flow converge

% 结果处理
res_Pgen = results.gen(:, PG);
res_Vbus = results.bus(:, VM);
res_Thetabus = results.bus(:, VA);
res_Pbranch = abs(results.branch(:, PF)); %TODO Also add PT
res_Pbranch_T = abs(results.branch(:,PT));

% 计算奖励和代价
Obj = sum(gencost(:, COST).*res_Pgen.^2 + gencost(:, COST+1).*res_Pgen + gencost(:, COST+1)) / 2000;
Cost1 = res_Vbus - 1.05;
Cost2 = 0.95 - res_Vbus;
Cost3 = (res_Pbranch - Prating)./Prating; % TODO: Add PT constraints and constraints on Q
 
Cost4 = (res_Pgen(BusRef) - Pgen_max(BusRef))/Pgen_max(BusRef);
Costs = [Cost1;Cost2;Cost3;Cost4];
g_val = [Obj;Costs];
