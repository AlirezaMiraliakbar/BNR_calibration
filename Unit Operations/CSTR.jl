using ModelingToolkit, Plots, DifferentialEquations, Unitful, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics

include("MaterialStream.jl")
include("../Kinetics Models/ASM3.jl")

@mtkmodel CSTR begin
        
        # Choosing the bio-kinetic model by @extend macro (e.g. ASM3)
        @extend ASM3()

        @components begin
            In = MaterialStream()
            Out = MaterialStream()
        end

        @parameters begin
            vol
            # Parameters for aeration
            KLa = 240 # 1 / day
            # Flag to turn on the aeration or shut down, 0 for shut down, 1 for turn on
            switch = 1

            Temp
            SOTE = 0.1
            # air_flow
        end

        @variables begin
            # Varaibles for aeration
            SO_sat_temp(t)
            KLa_temp(t)
        end


        # any variables ended with "_c" still represent the same physical variables, but defined as parameter-like variables, because these variables are going to be redefined in the equation part soon.
        @equations begin
            # The following equations are the mass balance equations for the CSTR
            D(x[1]) ~ 1.0/vol*(In.flow_rate*(In.S_O-x[1])) + reac1 + switch*KLa_temp*(SO_sat_temp-x[1])       # S_O2 (SO2) # S_O2 (SO2)
            D(x[2]) ~ 1.0/vol*(In.flow_rate*(In.S_I-x[2])) + reac2     # S_I (Si)
            D(x[3]) ~ 1.0/vol*(In.flow_rate*(In.S_S-x[3])) + reac3      # S_S (Ss)
            D(x[4]) ~ 1.0/vol*(In.flow_rate*(In.S_NH-x[4])) + reac4       # S_NH (Snh4)
            D(x[5]) ~ 1.0/vol*(In.flow_rate*(In.S_N2-x[5])) + reac5       # S_N2 (Sn2)
            D(x[6]) ~ 1.0/vol*(In.flow_rate*(In.S_NO-x[6])) + reac6       # S_NO (Snox)
            D(x[7]) ~ 1.0/vol*(In.flow_rate*(In.S_ALK-x[7])) + reac7      # S_HCO (Salk)
            D(x[8]) ~ 1.0/vol*(In.flow_rate*(In.X_I-x[8])) + reac8        # X_I (Xi)
            D(x[9]) ~ 1.0/vol*(In.flow_rate*(In.X_S-x[9])) + reac9      # X_S (Xs)
            D(x[10]) ~ 1.0/vol*(In.flow_rate*(In.X_H-x[10])) + reac10   # X_H (Xbh)
            D(x[11]) ~ 1.0/vol*(In.flow_rate*(In.X_STO-x[11])) + reac11    # X_STO (Xsto)
            D(x[12]) ~ 1.0/vol*(In.flow_rate*(In.X_A-x[12])) + reac12   # X_A (Xba)
            D(x[13]) ~ 1.0/vol*(In.flow_rate*(In.X_TS-x[13])) + reac13   # X_TS (TSS) 
            
            # Specify the relationship of outflow and the CSTR
            # State variables of outflow are the same as the reactors'
            Out.S_O ~ x[1]
            Out.S_I ~ x[2]
            Out.S_S ~ x[3]
            Out.S_NH ~ x[4]
            Out.S_N2 ~ x[5]
            Out.S_NO ~ x[6]
            Out.S_ALK ~ x[7]
            Out.X_I ~ x[8]
            Out.X_S ~ x[9]
            Out.X_H ~ x[10]
            Out.X_STO ~ x[11]
            Out.X_A ~ x[12]
            Out.X_TS ~ x[13]

            # inflow = outflow
            Out.flow_rate ~ In.flow_rate
            # Temperature of outflow = temperature of reactor

            # For aeration
            SO_sat_temp ~ 0.9997743214*8.0/10.5*(56.12*6791.5*exp(-66.7354 + 87.4755/((Temp+273.15)/100.0) + 24.4526*log((Temp+273.15)/100.0))) # van't Hoff equation for saturation concentration of O2, contributing to the SO2
            # KLa ~ 0.2967 * SOTE * air_flow / (V * SO_sat_temp) # Chenyu's model
            KLa_temp ~ KLa*1.024^(Temp-15.0)
        end
end