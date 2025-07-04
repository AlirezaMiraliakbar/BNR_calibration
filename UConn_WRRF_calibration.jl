using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D


include("./Unit Operations/MaterialStream.jl");
include("./Unit Operations/Mixer.jl");
include("./Unit Operations/Splitter.jl");
include("./Unit Operations/CSTR.jl");
include("./Influent/InfluentGenerator.jl");

influent = generate_influent_ss();
inlet_temp = influent[1];
inlet_comp = influent[2:end];

 
# UConn average flowrate in m3 / day -> tau for each reactor = 0.3 days -> 1.5 days for whole BNR system
inlet_flowrate = 3785.42 

# time of simulation in days
t_end = 100;

@mtkmodel BNRPlantModel begin

    @components begin

        Influent = InletStream(comp = inlet_comp, flow_rate = inlet_flowrate)
        reactor1 = CSTR(vol = 1135.6, switch=0, Temp = inlet_temp)
        reactor2 = CSTR(vol = 1135.6, switch=1, Temp = inlet_temp)
        reactor3 = CSTR(vol = 1135.6, switch=0, Temp = inlet_temp)
        reactor4 = CSTR(vol = 1135.6, switch=1, Temp = inlet_temp)
        reactor5 = CSTR(vol = 1135.6, switch=0, Temp = inlet_temp)
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

@mtkcompile BNR_system = BNRPlantModel();

