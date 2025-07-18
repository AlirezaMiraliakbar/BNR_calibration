using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using JuMP: @variable, @objective, @constraint, Model, optimize!, termination_status, value
using EOptInterface
using CSV, DataFrames, Distributions
using Random, Dates, Unitful
using EAGO
#=============================================================================================================#
include("./Influent/InfluentGenerator.jl");
include("./Effluent/EffluentGenerator.jl");
# ============================================== INFLUENT DATA ===============================================#
influent = generate_influent_ss();
inlet_temp = influent[1];
inlet_comp = influent[2:end];
# UConn average flowrate in m3 / day -> tau for each reactor = 0.3 days -> 1.5 days for whole BNR system
inlet_flowrate = 3785.42 
# ============================================= EFFLUENT DATA ================================================#
effluent = generate_effluent_ss();
effluent_comp = effluent[2:end];
# ============================================= SIMULATION TIME ==============================================#
# time of simulation in days
t_end = 100;

init_vec = copy(inlet_comp);
init_vec[10] = 5;
init_vec[11] = 5;
init_vec[12] = 5;

init_vec_anoxic = copy(init_vec);
init_vec_anoxic[1] = 0.0;


# ============================================ MODEL DEFINITION ==============================================#
function compile_BNR(p)
    #=param_ranges = [
    [0.1, 1],      # p1: K_O2
    [0.1, 1],      # p2: b_HO2
    [0.1, 1],      # p3: K_AO2
    [0.1, 1],      # p4: Y_HO2
    [0.1, 1],      # p5: Y_A
]=#
    p1 = p[1]  # K_O2
    p2 = p[2]  # b_HO2
    p3 = p[3]  # K_AO2
    p4 = p[4]  # Y_HO2
    p5 = p[5]  # Y_A

    @connector MaterialStream begin
        flow_rate(t) , [input = true];
        x1(t) , [input=true];
        x2(t) , [input=true];
        x3(t) , [input=true];
        x4(t) , [input=true];
        x5(t) , [input=true];
        x6(t) ,[input=true];
        x7(t) , [input=true];
        x8(t) , [input=true];
        x9(t) , [input=true];
        x10(t) , [input=true];
        x11(t) , [input=true];
        x12(t) , [input=true];
        x13(t) , [input=true];
    end


    @mtkmodel InletStream begin
        @components begin
            port = MaterialStream()
        end

        @parameters begin
            comp[1:13]  # Concentrations of the influent stream
            flow_rate  # Flow rate of the influent stream
        end

        @equations begin
            port.x1 ~ comp[1]
            port.x2 ~ comp[2]
            port.x3 ~ comp[3]
            port.x4 ~ comp[4]
            port.x5 ~ comp[5]
            port.x6 ~ comp[6]
            port.x7 ~ comp[7]
            port.x8 ~ comp[8]
            port.x9 ~ comp[9]
            port.x10 ~ comp[10]
            port.x11 ~ comp[11]
            port.x12 ~ comp[12]
            port.x13 ~ comp[13]

            port.flow_rate ~ flow_rate

        end
    end

    @mtkmodel Mixer_2in_1out begin
        @components begin
            In1 = MaterialStream()
            In2 = MaterialStream()
            Out1 = MaterialStream()
        end

        @equations begin
        
            Out1.flow_rate ~ In1.flow_rate + In2.flow_rate

            Out1.x1 ~ (In1.x1*In1.flow_rate + In2.x1*In2.flow_rate) / Out1.flow_rate
            Out1.x2 ~ (In1.x2*In1.flow_rate + In2.x2*In2.flow_rate) / Out1.flow_rate
            Out1.x3 ~ (In1.x3*In1.flow_rate + In2.x3*In2.flow_rate) / Out1.flow_rate
            Out1.x4 ~ (In1.x4*In1.flow_rate + In2.x4*In2.flow_rate) / Out1.flow_rate
            Out1.x5 ~ (In1.x5*In1.flow_rate + In2.x5*In2.flow_rate) / Out1.flow_rate
            Out1.x6 ~ (In1.x6*In1.flow_rate + In2.x6*In2.flow_rate) / Out1.flow_rate
            Out1.x7 ~ (In1.x7*In1.flow_rate + In2.x7*In2.flow_rate) / Out1.flow_rate
            Out1.x8 ~ (In1.x8*In1.flow_rate + In2.x8*In2.flow_rate) / Out1.flow_rate
            Out1.x9 ~ (In1.x9*In1.flow_rate + In2.x9*In2.flow_rate) / Out1.flow_rate
            Out1.x10 ~ (In1.x10*In1.flow_rate + In2.x10*In2.flow_rate) / Out1.flow_rate
            Out1.x11 ~ (In1.x11*In1.flow_rate + In2.x11*In2.flow_rate) / Out1.flow_rate
            Out1.x12 ~ (In1.x12*In1.flow_rate + In2.x12*In2.flow_rate) / Out1.flow_rate
            Out1.x13 ~ (In1.x13*In1.flow_rate + In2.x13*In2.flow_rate) / Out1.flow_rate
            
        end
    end

    @mtkmodel Splitter_1in_3out begin
        @components begin
            In1 = MaterialStream()
            Out1 = MaterialStream()
            Out2 = MaterialStream()
            Out3 = MaterialStream()
        end

        @parameters begin
            Out_factor[1:3] = [0.5,0.4,0.1]
        end

        @equations begin
            # Mass balance equations
            Out1.x1 ~ In1.x1
            Out2.x1 ~ In1.x1
            Out3.x1 ~ In1.x1

            Out1.x2 ~ In1.x2
            Out2.x2 ~ In1.x2
            Out3.x2 ~ In1.x2

            Out1.x3 ~ In1.x3
            Out2.x3 ~ In1.x3
            Out3.x3 ~ In1.x3

            Out1.x4 ~ In1.x4
            Out2.x4 ~ In1.x4
            Out3.x4 ~ In1.x4

            Out1.x5 ~ In1.x5
            Out2.x5 ~ In1.x5
            Out3.x5 ~ In1.x5

            Out1.x6 ~ In1.x6
            Out2.x6 ~ In1.x6
            Out3.x6 ~ In1.x6

            Out1.x7 ~ In1.x7
            Out2.x7 ~ In1.x7
            Out3.x7 ~ In1.x7

            Out1.x8 ~ In1.x8
            Out2.x8 ~ In1.x8
            Out3.x8 ~ In1.x8

            Out1.x9 ~ In1.x9
            Out2.x9 ~ In1.x9
            Out3.x9 ~ In1.x9

            Out1.x10 ~ In1.x10
            Out2.x10 ~ In1.x10
            Out3.x10 ~ In1.x10

            Out1.x11 ~ In1.x11
            Out2.x11 ~ In1.x11
            Out3.x11 ~ In1.x11

            Out1.x12 ~ In1.x12
            Out2.x12 ~ In1.x12
            Out3.x12 ~ In1.x12

            Out1.x13 ~ In1.x13
            Out2.x13 ~ In1.x13
            Out3.x13 ~ In1.x13


            # [Out1.x[i] ~ In1.x[i] for i = 1:13]...  # Concentrations are equal
            # [Out2.x[i] ~ In1.x[i] for i = 1:13]...
            # [Out3.x[i] ~ In1.x[i] for i = 1:13]...
            
            # Flow balance
            Out1.flow_rate ~ In1.flow_rate * Out_factor[1]
            Out2.flow_rate ~ In1.flow_rate * Out_factor[2]
            Out3.flow_rate ~ In1.flow_rate - (Out1.flow_rate + Out2.flow_rate)  # Ensure mass balance
        end
    end

    @mtkmodel ASM3_CSTR_Aerator begin
        
        @components begin
            In = MaterialStream()
            Out = MaterialStream()
        end

        @parameters begin
            vol
            # Parameters for aeration
            KLa = 240
            # Flag to turn on the aeration or shut down, 0 for shut down, 1 for turn on
            switch = 1

            Temp = 15

            # Kinetic parameters
            k_H = 3.0
            K_X = 1.5
            k_STO = 5.0
            ny_NOX = 0.6
            K_O2 = p1
            K_NOX = 0.5
            K_S = 2.0
            K_STO = 1.0
            mu_H = 2.0
            K_NH4 = 0.01
            K_ALK = 0.1
            b_HO2 = p2
            b_HNOX = 0.1
            b_STOO2 = 0.2
            b_STONOX = 0.1
            mu_A = 1.0
            K_ANH4 = 1.0
            K_AO2 = p3
            K_AALK = 0.5
            b_AO2 = 0.15
            b_ANOX = 0.05

            # Stoichiometric parameters
            f_SI = 0.0
            Y_STOO2 = 0.85
            Y_STONOX = 0.80
            Y_HO2 = p4
            Y_HNOX = 0.54
            Y_A = p5
            f_XI = 0.2
            i_NSI = 0.01
            i_NSS = 0.03
            i_NXI = 0.02
            i_NXS = 0.04
            i_NBM = 0.07
            i_SSXI = 0.075
            i_SSXS = 0.075
            i_SSBM = 0.90
            i_SSSTO = 0.60

            # Stoichiometric numbers from Table 1 
            # obtained by \sum_i^12 νji*ikI
            x1 = 1.0-f_SI
            x2 = -1.0+Y_STOO2
            x3 = (-1.0+Y_STONOX)/(64.0/14.0-24.0/14.0)
            x4 = 1.0-1.0/Y_HO2
            x5 = (+1.0-1.0/Y_HNOX)/(64.0/14.0-24.0/14.0)
            x6 = -1.0+f_XI
            x7 = (f_XI-1.0)/(64.0/14.0-24.0/14.0)
            x8 = -1.0
            x9 = -1.0/(64.0/14.0-24.0/14.0)
            x10 = -(64.0/14.0)/Y_A+1.0
            x11 = f_XI-1.0
            x12 = (f_XI-1.0)/(64.0/14.0-24.0/14.0)

            y1 = -f_SI*i_NSI-(1.0-f_SI)*i_NSS+i_NXS
            y2 = i_NSS
            y3 = i_NSS
            y4 = -i_NBM
            y5 = -i_NBM
            y6 = -f_XI*i_NXI+i_NBM
            y7 = -f_XI*i_NXI+i_NBM
            y10 = -1.0/Y_A-i_NBM
            y11 = -f_XI*i_NXI+i_NBM
            y12 = -f_XI*i_NXI+i_NBM

            z1 = y1/14.0
            z2 = y2/14.0
            z3 = y3/14.0-x3/14.0
            z4 = y4/14.0
            z5 = y5/14.0-x5/14.0
            z6 = y6/14.0
            z7 = y7/14.0-x7/14.0
            z9 = -x9/14.0
            z10 = y10/14.0-1.0/(Y_A*14.0)
            z11 = y11/14.0
            z12 = y12/14.0-x12/14.0

            t1 = -i_SSXS
            t2 = Y_STOO2*i_SSSTO
            t3 = Y_STONOX*i_SSSTO
            t4 = i_SSBM-1.0/Y_HO2*i_SSSTO
            t5 = i_SSBM-1.0/Y_HNOX*i_SSSTO
            t6 = f_XI*i_SSXI-i_SSBM
            t7 = f_XI*i_SSXI-i_SSBM
            t8 = -i_SSSTO
            t9 = -i_SSSTO
            t10 = i_SSBM
            t11 = f_XI*i_SSXI-i_SSBM
            t12 = f_XI*i_SSXI-i_SSBM
            
            # Oxygen Mass Transfer Parameters

            # van't Hoff equation for saturation concentration of O2
            SO_sat_temp =  0.9997743214*8.0/10.5*(56.12*6791.5*exp(-66.7354 + 87.4755/((Temp+273.15)/100.0) + 24.4526*log((Temp+273.15)/100.0))) 
            KLa_temp  = KLa*1.024^(Temp-15.0)
        end

        @variables begin
            # State variables 
            # x[1:13] are:
            # S_O2 (O2), S_I (COD), S_s (COD), S_NH4 (N), S_N2 (N), S_NOX (N), S_ALK (Mole), X_I (COD), X_S (COD), X_H (COD), X_STO (COD), X_A (COD), X_SS (SS)
            x(t)[1:13]

            # Kinectic rate expressions Table 2 
            # Hydrolysis
            proc1(t)
            # Heterotrophic organisms, denitrification
            proc2(t); proc3(t); proc4(t); proc5(t); proc6(t); proc7(t); proc8(t); proc9(t)
            # Autotrophic organisms, nitrification
            proc10(t); proc11(t); proc12(t)

            reac1(t); reac2(t); reac3(t); reac4(t); reac5(t)
            reac6(t); reac7(t); reac8(t); reac9(t); reac10(t)
            reac11(t); reac12(t); reac13(t); reac15(t); reac16(t)
            reac17(t); reac18(t); reac19(t)

        end


        # any variables ended with "_c" still represent the same physical variables, but defined as parameter-like variables, because these variables are going to be redefined in the equation part soon.
        @equations begin
            proc1 ~ k_H*(x[9]/x[10])/(K_X+x[9]/x[10])*x[10]
            proc2 ~ k_STO*x[1]/(K_O2+x[1])*x[3]/(K_S+x[3])*x[10]
            proc3 ~ k_STO*ny_NOX*K_O2/(K_O2+x[1])*x[6]/(K_NOX+x[6])*x[3]/(K_S+x[3])*x[10]
            proc4 ~ mu_H*x[1]/(K_O2+x[1])*x[4]/(K_NH4+x[4])*x[7]/(K_ALK+x[7])*(x[11]/x[10])/(K_STO+x[11]/x[10])*x[10]
            proc5 ~ mu_H*ny_NOX*K_O2/(K_O2+x[1])*x[6]/(K_NOX+x[6])*x[4]/(K_NH4+x[4])*x[7]/(K_ALK+x[7])*(x[11]/x[10])/(K_STO+x[11]/x[10])*x[10]
            proc6 ~ b_HO2*x[1]/(K_O2+x[1])*x[10]
            proc7 ~ b_HNOX*K_O2/(K_O2+x[1])*x[6]/(K_NOX+x[6])*x[10]
            proc8 ~ b_STOO2*x[1]/(K_O2+x[1])*x[11]
            proc9 ~ b_STONOX*K_O2/(K_O2+x[1])*x[6]/(K_NOX+x[6])*x[11]
            proc10 ~ mu_A*x[1]/(K_AO2+x[1])*x[4]/(K_ANH4+x[4])*x[7]/(K_AALK+x[7])*x[12]
            proc11 ~ b_AO2*x[1]/(K_AO2+x[1])*x[12]
            proc12 ~ b_ANOX*K_AO2/(K_AO2+x[1])*x[6]/(K_NOX+x[6])*x[12]

            # use kinetic rate expression from table 2 times the stoichiometric matrix in table 1.
            reac1 ~ x2*proc2+x4*proc4+x6*proc6+x8*proc8+x10*proc10+x11*proc11
            reac2 ~ f_SI*proc1
            reac3 ~ x1*proc1-proc2-proc3
            reac4 ~ y1*proc1+y2*proc2+y3*proc3+y4*proc4+y5*proc5+y6*proc6+y7*proc7+y10*proc10+y11*proc11+y12*proc12
            reac5 ~ -x3*proc3-x5*proc5-x7*proc7-x9*proc9-x12*proc12
            reac6 ~ x3*proc3+x5*proc5+x7*proc7+x9*proc9+1.0/Y_A*proc10+x12*proc12
            reac7 ~ z1*proc1+z2*proc2+z3*proc3+z4*proc4+z5*proc5+z6*proc6+z7*proc7+z9*proc9+z10*proc10+z11*proc11+z12*proc12
            reac8 ~ f_XI*proc6+f_XI*proc7+f_XI*proc11+f_XI*proc12
            reac9 ~ -proc1
            reac10 ~ proc4+proc5-proc6-proc7
            reac11 ~ Y_STOO2*proc2+Y_STONOX*proc3-1.0/Y_HO2*proc4-1.0/Y_HNOX*proc5-proc8-proc9
            reac12 ~ proc10-proc11-proc12
            reac13 ~ t1*proc1+t2*proc2+t3*proc3+t4*proc4+t5*proc5+t6*proc6+t7*proc7+t8*proc8+t9*proc9+t10*proc10+t11*proc11+t12*proc12
            # The following equations are the mass balance equations for the CSTR
            D(x[1]) ~ 1.0/vol*(In.flow_rate*(In.x1-x[1])) + reac1 + switch*KLa_temp*(SO_sat_temp-x[1])       # S_O2 (SO2) # S_O2 (SO2)
            D(x[2]) ~ 1.0/vol*(In.flow_rate*(In.x2-x[2])) + reac2     # S_I (Si)
            D(x[3]) ~ 1.0/vol*(In.flow_rate*(In.x3-x[3])) + reac3      # S_S (Ss)
            D(x[4]) ~ 1.0/vol*(In.flow_rate*(In.x4-x[4])) + reac4       # S_NH (Snh4)
            D(x[5]) ~ 1.0/vol*(In.flow_rate*(In.x5-x[5])) + reac5       # S_N2 (Sn2)
            D(x[6]) ~ 1.0/vol*(In.flow_rate*(In.x6-x[6])) + reac6       # S_NO (Snox)
            D(x[7]) ~ 1.0/vol*(In.flow_rate*(In.x7-x[7])) + reac7      # S_HCO (Salk)
            D(x[8]) ~ 1.0/vol*(In.flow_rate*(In.x8-x[8])) + reac8        # X_I (Xi)
            D(x[9]) ~ 1.0/vol*(In.flow_rate*(In.x9-x[9])) + reac9      # X_S (Xs)
            D(x[10]) ~ 1.0/vol*(In.flow_rate*(In.x10-x[10])) + reac10   # X_H (Xbh)
            D(x[11]) ~ 1.0/vol*(In.flow_rate*(In.x11-x[11])) + reac11    # X_STO (Xsto)
            D(x[12]) ~ 1.0/vol*(In.flow_rate*(In.x12-x[12])) + reac12   # X_A (Xba)
            D(x[13]) ~ 1.0/vol*(In.flow_rate*(In.x13-x[13])) + reac13   # X_TS (TSS) 
            
            # Specify the relationship of outflow and the CSTR
            # State variables of outflow are the same as the reactors'
            Out.x1 ~ x[1]
            Out.x2 ~ x[2]
            Out.x3 ~ x[3]
            Out.x4 ~ x[4]
            Out.x5 ~ x[5]
            Out.x6 ~ x[6]
            Out.x7 ~ x[7]
            Out.x8 ~ x[8]
            Out.x9 ~ x[9]
            Out.x10 ~ x[10]
            Out.x11 ~ x[11]
            Out.x12 ~ x[12]
            Out.x13 ~ x[13]

            # inflow = outflow
            Out.flow_rate ~ In.flow_rate
        end
    end


    @mtkmodel BNRPlantModel begin
        @components begin

            Influent = InletStream(comp = inlet_comp, flow_rate = inlet_flowrate)
            reactor1 = ASM3_CSTR_Aerator(vol = 1135.6, switch=0, Temp = 20)
            reactor2 = ASM3_CSTR_Aerator(vol = 1135.6, switch=1, Temp = 20)
            reactor3 = ASM3_CSTR_Aerator(vol = 1135.6, switch=0, Temp = 20)
            reactor4 = ASM3_CSTR_Aerator(vol = 1135.6, switch=1, Temp = 20)
            reactor5 = ASM3_CSTR_Aerator(vol = 1135.6, switch=0, Temp = 20)
            splitter1 = Splitter_1in_3out(Out_factor = [0.5,0.4,0.1])
            mixer1 = Mixer_2in_1out()
            mixer2 = Mixer_2in_1out()

        end
    
        @equations begin
            connect(Influent.port,mixer1.In1)
            connect(splitter1.Out3,mixer1.In2)
            connect(mixer1.Out1,reactor1.In)
            connect(reactor1.Out,mixer2.In1)
            connect(splitter1.Out1,reactor2.In)
            connect(reactor2.Out,mixer2.In2)
            connect(mixer2.Out1,reactor3.In)
            connect(reactor3.Out,reactor4.In)
            connect(reactor4.Out,reactor5.In)
            connect(reactor5.Out,splitter1.In1)
        end
    end

    @mtkcompile BNR_system = BNRPlantModel() 

    return BNR_system
end

function solve_BNR(p)
    BNR_system = compile_BNR(p)
    #TODO: if the number of reactors becomes a decision variables then we need to modify the initial conditions accordingly

    u0 = vcat(
        [BNR_system.reactor1.x[i] => init_vec_anoxic[i] for i = 1:13],
        [BNR_system.reactor2.x[i] => init_vec[i] for i = 1:13],
        [BNR_system.reactor3.x[i] => init_vec[i] for i = 1:13],
        [BNR_system.reactor4.x[i] => init_vec[i] for i = 1:13],
        [BNR_system.reactor5.x[i] => init_vec[i] for i = 1:13]
        ) 

    prob = ODEProblem(BNR_system, u0, (0, t_end); fully_determined = true)

    sol = solve(prob, Rodas5P(), initializealg=BrownFullBasicInit())

    # Extracting plant performance metrics
    SNH_EFF = sol[BNR_system.splitter1.Out2.x4][end]
    
    SNO_EFF = sol[BNR_system.splitter1.Out2.x5][end]

    COD_EFF = sol[BNR_system.splitter1.Out2.x3][end] + 
          sol[BNR_system.splitter1.Out2.x2][end] + 
          sol[BNR_system.splitter1.Out2.x9][end] + 
          sol[BNR_system.splitter1.Out2.x8][end] +
          sol[BNR_system.splitter1.Out2.x10][end] +
          sol[BNR_system.splitter1.Out2.x11][end] +
          sol[BNR_system.splitter1.Out2.x12][end]

    TSS_EFF = sol[BNR_system.splitter1.Out2.x13][end]

    eff_metric = [SNH_EFF, SNO_EFF, COD_EFF, TSS_EFF]

    return eff_metric
end

params = [10 for i in 1:26]  # Initial guess for the 26 parameters
sys = compile_BNR(params)
decision_vars(sys)
# ================================================= Optimization formulation ===============================================#
param_ranges = [
    [0.1, 1],      # p5: K_O2
    [0.1, 1],      # p12: b_HO2
    [0.1, 1],      # p18: K_AO2
    [0.1, 1],      # p24: Y_HO2
    [0.1, 1],      # p26: Y_A
]

pL = [r[1] for r in param_ranges];
pU = [r[2] for r in param_ranges];

function obj_func(p)
    eff_metric = solve_BNR(p)
    # compare with effluent data
    effluent_comp = effluent[2:end];
    # ============ COD ==============#
    COD_model = eff_metric[3]
    COD_data = effluent_comp[4]

    err_COD = (COD_model - COD_data)^2

    # ============ SNH ==============#
    SNH_model = eff_metric[1]
    SNH_data = effluent_comp[2]

    err_SNH = (SNH_model - SNH_data)^2


    # ============ SNO ==============#
    SNO_model = eff_metric[2]
    SNO_data = effluent_comp[3]
    err_SNO = (SNO_model - SNO_data)^2

    # ============ TSS ==============#
    TSS_model = eff_metric[4]
    TSS_data = effluent_comp[6]

    err_TSS = (TSS_model - TSS_data)^2

    # ============ MSE ==============#
    SE = err_COD + err_SNH + err_SNO + err_TSS #TODO: modify after proper program definition

    return SE
end

function obj_func_numeric(p_vals...)
    p = collect(p_vals)
    return obj_func(p)
end
# Define the optimization model
model = Model(EAGO.Optimizer)
@variable(model, pL[i] <= p[i=1:5] <= pU[i])
@objective(model, Min, obj_func_numeric(p...))
# having constraint here depends on how we formulate the problem, if the ODEs are implicily included in the objective function, then we don't need to add constraints here.
optimize!(model)


