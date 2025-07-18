# reference
# [1] Gujer, W., Henze, M., Mino, T., & Van Loosdrecht, M. (1999). Activated sludge model No. 3. Water science and technology, 39(1), 183-193.

# [2] Henze, M., Gujer, W., Mino, T., & Van Loosedrecht, M. (2006). Activated sludge models ASM1, ASM2, ASM2d and ASM3. IWA publishing.

using ModelingToolkit, Plots, DifferentialEquations, Unitful, IfElse, NativeNaNMath
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics

include("./ASM3_Parameters.jl")
params_asm3 = default_asm3_parameters();

p1 = params_asm3[:k_H]
p2 = params_asm3[:K_X]
p3 = params_asm3[:k_STO]
p4 = params_asm3[:ny_NOX]
p5 = params_asm3[:K_O2]
p6 = params_asm3[:K_NOX]
p7 = params_asm3[:K_S]
p8 = params_asm3[:K_STO]
p9 = params_asm3[:mu_H]
p10 = params_asm3[:K_NH4]
p11 = params_asm3[:K_ALK]
p12 = params_asm3[:b_HO2]
p13 = params_asm3[:b_HNOX]
p14 = params_asm3[:b_STOO2]
p15 = params_asm3[:b_STONOX]
p16 = params_asm3[:mu_A]
p17 = params_asm3[:K_ANH4]
p18 = params_asm3[:K_AO2]
p19 = params_asm3[:K_AALK]
p20 = params_asm3[:b_AO2]
p21 = params_asm3[:b_ANOX]
p22 = params_asm3[:f_SI]
p23 = params_asm3[:Y_STOO2]
p24 = params_asm3[:Y_STONOX]
p25 = params_asm3[:Y_HO2]
p26 = params_asm3[:Y_HNOX]
p27 = params_asm3[:Y_A]
p28 = params_asm3[:f_XI]
p29 = params_asm3[:i_NSI]
p30 = params_asm3[:i_NSS]
p31 = params_asm3[:i_NXI]
p32 = params_asm3[:i_NXS]
p33 = params_asm3[:i_NBM]
p34 = params_asm3[:i_SSXI]
p35 = params_asm3[:i_SSXS]
p36 = params_asm3[:i_SSBM]
p37 = params_asm3[:i_SSSTO]

@mtkmodel ASM3 begin

    @parameters begin

        # Kinetic parameters
        k_H = p1
        K_X = p2
        k_STO = p3
        ny_NOX = p4
        K_O2 = p5
        K_NOX = p6
        K_S = p7
        K_STO = p8
        mu_H = p9
        K_NH4 = p10
        K_ALK = p11
        b_HO2 = p12
        b_HNOX = p13
        b_STOO2 = p14
        b_STONOX = p15
        mu_A = p16
        K_ANH4 = p17
        K_AO2 = p18
        K_AALK = p19
        b_AO2 = p20
        b_ANOX = p21

        # Stoichiometric parameters
        f_SI = p22
        Y_STOO2 = p23
        Y_STONOX = p24
        Y_HO2 = p25
        Y_HNOX = p26
        Y_A = p27
        f_XI = p28
        i_NSI = p29
        i_NSS = p30
        i_NXI = p31
        i_NXS = p32
        i_NBM = p33
        i_SSXI = p34
        i_SSXS = p35
        i_SSBM = p36
        i_SSSTO = p37

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
    end

    @variables begin
        # State variables 
        # x[1:13] are:
        # S_O2 (O2), S_I (COD), S_s (COD), S_NH4 (N), S_N2 (N), S_NOX (N), S_ALK (Mole), X_I (COD), X_S (COD), X_H (COD), X_STO (COD), X_A (COD), X_SS (SS)
        x(t)[1:13]                                                                                                                                          # Table 3 Autotrophic organisms, nitrification, X_A

        # Kinectic rate expressions Table 2 
        # Hydrolysis
        proc1(t); 
        # Heterotrophic organisms, denitrification
        proc2(t); proc3(t); proc4(t); proc5(t); proc6(t); proc7(t); proc8(t); proc9(t); 
        # Autotrophic organisms, nitrification
        proc10(t); proc11(t); proc12(t)

        # T used in the interpolation of kinetic parameters 
        reac1(t); reac2(t); reac3(t); reac4(t); reac5(t); 
        reac6(t); reac7(t); reac8(t); reac9(t); reac10(t);
         reac11(t); reac12(t); reac13(t); reac15(t); reac16(t); 
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
        reac1 ~ x2*proc2+x4*proc4+x6*proc6+x8*proc8+x10*proc10+x11*proc11 # g/d/L
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
    end
end

