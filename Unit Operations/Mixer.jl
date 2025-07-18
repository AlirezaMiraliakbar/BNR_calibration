# using ModelingToolkit
# using ModelingToolkit: t_nounits as t, D_nounits as D


include("MaterialStream.jl")

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