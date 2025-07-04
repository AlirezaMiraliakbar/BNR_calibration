using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

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
        end

        @variables begin
            # Varaibles for aeration
            SO_sat_temp(t)
            KLa_temp(t)
        end


        # any variables ended with "_c" still represent the same physical variables, but defined as parameter-like variables, because these variables are going to be redefined in the equation part soon.
        @equations begin
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
            # Temperature of outflow = temperature of reactor

            # For aeration
            SO_sat_temp ~ 0.9997743214*8.0/10.5*(56.12*6791.5*exp(-66.7354 + 87.4755/((Temp+273.15)/100.0) + 24.4526*log((Temp+273.15)/100.0))) # van't Hoff equation for saturation concentration of O2, contributing to the SO2
            KLa_temp ~ KLa*1.024^(Temp-15.0)
        end
end