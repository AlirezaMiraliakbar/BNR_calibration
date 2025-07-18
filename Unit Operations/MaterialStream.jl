# using ModelingToolkit
# using ModelingToolkit: t_nounits as t, D_nounits as D

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