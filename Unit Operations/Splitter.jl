using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D


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