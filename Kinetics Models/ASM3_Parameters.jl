function default_asm3_parameters()
    Dict(
        # Kinetic parameters
        :k_H => 3.0,
        :K_X => 1.5,
        :k_STO => 5.0,
        :ny_NOX => 0.6,
        :K_O2 => 0.2, 
        :K_NOX => 0.5, 
        :K_S => 2.0, 
        :K_STO => 1.0,
        :mu_H => 2.0,
        :K_NH4 => 0.01, 
        :K_ALK => 0.1, 
        :b_HO2 => 0.2,
        :b_HNOX => 0.1,
        :b_STOO2 => 0.2,
        :b_STONOX => 0.1,
        :mu_A => 1.0,
        :K_ANH4 => 1.0, 
        :K_AO2 => 0.5, 
        :K_AALK => 0.5, 
        :b_AO2 => 0.15,
        :b_ANOX => 0.05,
        
        # Stoichiometric parameters
        :f_SI => 0.0,
        :Y_STOO2 => 0.85,
        :Y_STONOX => 0.80,
        :Y_HO2 => 0.63,
        :Y_HNOX => 0.54,
        :Y_A => 0.24,
        :f_XI => 0.2,
        :i_NSI => 0.01,
        :i_NSS => 0.03,
        :i_NXI => 0.02,
        :i_NXS => 0.04,
        :i_NBM => 0.07,
        :i_SSXI => 0.75,
        :i_SSXS => 0.75,
        :i_SSBM => 0.90,
        :i_SSSTO => 0.60
    )
end
