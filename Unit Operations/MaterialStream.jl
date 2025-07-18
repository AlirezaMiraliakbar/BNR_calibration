using ModelingToolkit, Plots, DifferentialEquations, Unitful, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics

@connector MaterialStream begin
        flow_rate(t) , [input = true];
        S_O(t) , [input=true];
        S_I(t) , [input=true];
        S_S(t) , [input=true];
        S_NH(t) , [input=true];
        S_N2(t) , [input=true];
        S_NO(t) ,[input=true];
        S_ALK(t) , [input=true];
        X_I(t) , [input=true];
        X_S(t) , [input=true];
        X_H(t) , [input=true];
        X_STO(t) , [input=true];
        X_A(t) , [input=true];
        X_TS(t) , [input=true];
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
            port.S_O ~ comp[1]
            port.S_I ~ comp[2]
            port.S_S ~ comp[3]
            port.S_NH ~ comp[4]
            port.S_N2 ~ comp[5]
            port.S_NO ~ comp[6]
            port.S_ALK ~ comp[7]
            port.X_I ~ comp[8]
            port.X_S ~ comp[9]
            port.X_H ~ comp[10]
            port.X_STO ~ comp[11]
            port.X_A ~ comp[12]
            port.X_TS ~ comp[13]

            port.flow_rate ~ flow_rate

        end
end