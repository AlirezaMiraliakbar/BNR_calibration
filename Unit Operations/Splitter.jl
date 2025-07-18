using ModelingToolkit, Plots, DifferentialEquations, Unitful, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics

include("MaterialStream.jl")

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
        #TODO: change the x2,x3,... to state variable names before anything
        @equations begin
            # Mass balance equations
            Out1.S_O ~ In1.S_O
            Out2.S_O ~ In1.S_O
            Out3.S_O ~ In1.S_O

            Out1.S_I ~ In1.S_I
            Out2.S_I ~ In1.S_I
            Out3.S_I ~ In1.S_I

            Out1.S_S ~ In1.S_S
            Out2.S_S ~ In1.S_S
            Out3.S_S ~ In1.S_S

            Out1.S_NH ~ In1.S_NH
            Out2.S_NH ~ In1.S_NH
            Out3.S_NH ~ In1.S_NH

            Out1.S_N2 ~ In1.S_N2
            Out2.S_N2 ~ In1.S_N2
            Out3.S_N2 ~ In1.S_N2

            Out1.S_NO ~ In1.S_NO
            Out2.S_NO ~ In1.S_NO
            Out3.S_NO ~ In1.S_NO

            Out1.S_ALK ~ In1.S_ALK
            Out2.S_ALK ~ In1.S_ALK
            Out3.S_ALK ~ In1.S_ALK

            Out1.X_I ~ In1.X_I
            Out2.X_I ~ In1.X_I
            Out3.X_I ~ In1.X_I

            Out1.X_S ~ In1.X_S
            Out2.X_S ~ In1.X_S
            Out3.X_S ~ In1.X_S

            Out1.X_H ~ In1.X_H
            Out2.X_H ~ In1.X_H
            Out3.X_H ~ In1.X_H

            Out1.X_STO ~ In1.X_STO
            Out2.X_STO ~ In1.X_STO
            Out3.X_STO ~ In1.X_STO

            Out1.X_A ~ In1.X_A
            Out2.X_A ~ In1.X_A
            Out3.X_A ~ In1.X_A

            Out1.X_TS ~ In1.X_TS
            Out2.X_TS ~ In1.X_TS
            Out3.X_TS ~ In1.X_TS


            # [Out1.x[i] ~ In1.x[i] for i = 1:13]...  # Concentrations are equal
            # [Out2.x[i] ~ In1.x[i] for i = 1:13]...
            # [Out3.x[i] ~ In1.x[i] for i = 1:13]...
            
            # Flow balance
            Out1.flow_rate ~ In1.flow_rate * Out_factor[1]
            Out2.flow_rate ~ In1.flow_rate * Out_factor[2]
            Out3.flow_rate ~ In1.flow_rate - (Out1.flow_rate + Out2.flow_rate)  # Ensure mass balance
        end
    end