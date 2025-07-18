using ModelingToolkit, Plots, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

# [1] Hiatt, W. C., & Grady, C. L. (2008). An updated process model for carbon oxidation, nitrification, and denitrification. Water Environment Research, 80(11), 2145-2156.

@mtkmodel ASMN begin
    @parameters begin
        # ref [1] Table 2 
        μ_B_H = 6.25      # Maximum specific growth rate, heterotrophs (1/day)
        b_L_H = 0.408     # Decay coefficient, heterotrophs (1/day)
        Y_H = 0.6         # Heterotrophic yield (mg biomass COD formed / mg substrate COD removed)
        
        η_Y = 0.9         # Anoxic yield factor (Dimensionless)
        η_g = 0.8         # Anoxic growth factor, R17, R18 (Dimensionless)
        η_g2 = 0.28       # Anoxic growth factor, R2 (Dimensionless)
        η_g3 = 0.16       # Anoxic growth factor, R3 (Dimensionless)
        η_g4 = 0.35       # Anoxic growth factor, R4 (Dimensionless)
        η_g5 = 0.35       # Anoxic growth factor, R5 (Dimensionless)
        
        K_S1 = 20.0       # Half-saturation coefficient for substrate, R1 (mg/L, as COD)
        K_S2 = 20.0       # Half-saturation coefficient for substrate, R2 (mg/L, as COD)
        K_S3 = 20.0       # Half-saturation coefficient for substrate, R3 (mg/L, as COD)
        K_S4 = 20.0       # Half-saturation coefficient for substrate, R4 (mg/L, as COD)
        K_S5 = 40.0       # Half-saturation coefficient for substrate, R5 (mg/L, as COD)
        
        K_O_H1 = 0.1      # Half-saturation coefficient for O₂, heterotrophs, R1 (mg/L, as O₂)
        K_O_H2 = 0.1      # Half-saturation coefficient for O₂, heterotrophs, R2 (mg/L, as O₂)
        K_O_H3 = 0.1      # Half-saturation coefficient for O₂, heterotrophs, R3 (mg/L, as O₂)
        K_O_H4 = 0.1      # Half-saturation coefficient for O₂, heterotrophs, R4 (mg/L, as O₂)
        K_O_H5 = 0.1      # Half-saturation coefficient for O₂, heterotrophs, R5 (mg/L, as O₂)
        
        K_NO3 = 0.2       # Half-saturation coefficient for nitrate-nitrogen, heterotrophs (mg/L, as nitrogen)
        K_NO2 = 0.2       # Half-saturation coefficient for nitrite-nitrogen, heterotrophs (mg/L, as nitrogen)
        K_NO = 0.05       # Half-saturation coefficient for nitric oxide-nitrogen, heterotrophs (mg/L, as nitrogen)
        K_N2O = 0.05      # Half-saturation coefficient for nitrous oxide-nitrogen, heterotrophs (mg/L, as nitrogen)
        
        K_I3NO = 0.5      # Nitric oxide inhibition coefficient, R3 (mg/L, as nitrogen)
        K_I4NO = 0.3      # Nitric oxide inhibition coefficient, R4 (mg/L, as nitrogen)
        K_I5NO = 0.075    # Nitric oxide inhibition coefficient, R5 (mg/L, as nitrogen)
        
        K_INO3 = 0.1      # Nitrate inhibition coefficient for nitrate-nitrogen, ANRA R6 (mg/L, as nitrogen)
        K_INO2 = 0.05     # Nitrite inhibition coefficient, ANRA R6 (mg/L, as nitrogen)
        K_I6NH = 0.05     # Ammonia inhibition coefficient, ANRA R6 (mg/L, as nitrogen)
        K_7NO2   = 0.1       # Half-saturation coefficient for nitrite-nitrogen, ANRA R7 (mg/L, as nitrogen)
        K_I7NH   = 0.05      # Ammonia inhibition coefficient, ANRA R7 (mg/L, as nitrogen)
        
        μ_B_A1   = 0.78      # Maximum specific growth rate, AOB (1/day)
        μ_B_A2   = 0.78      # Maximum specific growth rate, NOB (autotrophic) (1/day)
        β_g      = 6.0       # NOB mixotrophic growth factor (dimensionless)
        
        b_L_A1   = 0.096     # Decay coefficient, AOB (1/day)
        b_L_A2   = 0.096     # Decay coefficient, NOB (1/day)
        
        Y_A1     = 0.18      # Autotrophic yield, AOB (mg biomass COD formed / mg nitrogen removed)
        Y_A2     = 0.06      # Autotrophic yield, NOB (mg biomass COD formed / mg nitrogen removed)
        
        K_FA     = 0.0075    # Half-saturation coefficient for free ammonia (mg/L)
        K_FNA    = 0.0001    # Half-saturation coefficient for free nitrous acid (mg/L)
        
        K_S11    = 20.0      # Half-saturation coefficient for substrate, Ss, R11 (mg/L, as COD)
        
        K_O_A1   = 0.6       # Half-saturation coefficient for O₂, AOB (mg/L, as O₂)
        K_O_A2   = 1.2       # Half-saturation coefficient for O₂, NOB (mg/L, as O₂)
        
        K_I9FA   = 1.0       # Free ammonia inhibition coefficient, R9 (mg/L, as nitrogen)
        K_I9FNA  = 0.1       # Free nitrous acid inhibition coefficient, R9 (mg/L, as nitrogen)
        K_I10FA  = 0.2       # Free ammonia inhibition coefficient, R10 (mg/L, as nitrogen)
        K_I10FNA = 0.04      # Free nitrous acid inhibition coefficient, R10 (mg/L, as nitrogen)
        
        f_prime_D = 0.08     # Fraction active biomass contributing to biomass debris (mg debris COD / mg biomass COD)
        
        i_N_XB   = 0.086     # Mass of nitrogen per mass of COD in active biomass (mg N / mg COD)
        i_N_XD   = 0.06      # Mass of nitrogen per mass of COD in biomass debris (mg N / mg COD)
        
        η_h      = 0.4       # Anoxic hydrolysis factor (dimensionless)
        k_a      = 0.1608    # Ammonification rate coefficient (L / (mg biomass COD · h))
        
        K_X      = 0.15      # Half-saturation coefficient for hydrolysis of slowly biodegradable substrate (mg COD / mg biomass COD)
        
        K_N1     = 0.1       # Half-saturation coefficient for NH₃, heterotrophs (mg/L, as nitrogen)
        k_h      = 2.208     # Hydrolysis coefficient (mg COD / (mg biomass COD · d))
    end

    @variables begin
        # ref [1] Table 1 
        X_I(t) = 0.0 # Inert particulate organic matter (mg/L as COD)
        X_S(t) = 0.0 # Slowly biodegradable substrate (mg/L as COD)
        X_B_H(t) = 0.0 # Active heterotrophic biomass (mg/L as COD)
        X_B_A1(t) = 0.0 # Active AOB biomass (mg/L as COD)
        X_B_A2(t) = 0.0 # Active NOB biomass (mg/L as COD)
        X_D(t) = 0.0 # Debris from biomass death and lysis (mg/L as COD)
        S_I(t) = 0.0 # Inert soluble organic matter (mg/L as COD)
        S_O(t) = 0.0 # Dissolved Oxygen concentration (mg/L as COD)
        S_S(t) = 0.0 # Readily (soluble) Biodegradable substrate concentration (mg/L as COD)
        S_NO3(t) = 0.0 # Nitrate concentration (mg/L as N)
        S_NO2(t) = 0.0 # Nitrite concentration (mg/L as N)
        S_NO(t) = 0.0 # Nitric oxide (mg/L as N)
        S_N2O(t) = 0.0 # Nitrous oxide (mg/L as N)
        S_NH(t) = 0.0 # Ammonia (mg/L as N)
        S_NS(t) = 0.0 # Soluble biodegradable organic nitrogen (mg/L as N)
        X_NS(t) = 0.0 # Particulate biodegradable organic nitrogen (mg/L as N)
        S_ALK(t) = 0.0 # Alkalinity (molar units)
        S_Z(t) = 0.0 # Salt (mg/L)
        S_J(t) = 0.0 # Biodegradable AOB inhibitor (mg/L as COD)
        S_W(t) = 0.0 # Priority pollutant (mg/L as COD)
        # ref [1] Page 4, The variable to represent summation of all the nitrogen oxides variables
        sum_NOX(t) = S_NO3 + S_NO2 +  S_NO + S_N2O # Sum of nitrogen oxides (mg/L as N)
        # ref [1] Table 3
        R1(t)  = 0.0  # Aerobic growth of heterotrophs
        R2(t)  = 0.0  # Anoxic growth of heterotrophs, reducing nitrate to nitrite
        R3(t)  = 0.0  # Anoxic growth of heterotrophs, reducing nitrite to nitric oxide
        R4(t)  = 0.0  # Anoxic growth of heterotrophs, reducing nitric oxide to nitrous oxide
        R5(t)  = 0.0  # Anoxic growth of heterotrophs, reducing nitrous oxide to nitrogen
        R6(t)  = 0.0  # ANRA, nitrate reduction to nitrite
        R7(t)  = 0.0  # ANRA, nitrite reduction to ammonia
        R8(t)  = 0.0  # Death and lysis of heterotrophs
        R9(t)  = 0.0  # Autotrophic growth of AOB
        R10(t) = 0.0  # Autotrophic growth of NOB
        R11(t) = 0.0  # Mixotrophic growth of NOB
        R12(t) = 0.0  # Death and lysis of AOB
        R13(t) = 0.0  # Death and lysis of NOB
        R14(t) = 0.0  # Ammonification of soluble organic nitrogen
        R15(t) = 0.0  # Hydrolysis of particulate organics
        R16(t) = 0.0  # Hydrolysis of particulate organic nitrogen
        R17(t) = 0.0  # Biodegradation of an inhibitory compound
        R18(t) = 0.0  # Biodegradation of a special interest compound
        proc
    end

    @equations begin
        # ref [1] Table 3
        R2  ~ μ_H * η_g2 * X_B_H * (S_S / (K_S2 + S_S)) * (S_NO3 / (K_NO3 + S_NO3)) * (S_O / (K_O_H2 + S_O))

        R3  ~ μ_H * η_g3 * X_B_H * (S_S / (K_S3 + S_S)) * (S_NO2 / (K_NO2 + S_NO2)) * (S_O / (K_O_H3 + S_O)) * (K_I3NO / (K_I3NO + S_NO))
        
        R4  ~ μ_H * η_g4 * X_B_H * (S_S / (K_S4 + S_S)) * (S_NO / (K_NO + S_NO + S_NO^2 / K_I4NO)) * (S_O / (K_O_H4 + S_O))
        
        R5  ~ μ_H * η_g5 * X_B_H * (S_S / (K_S5 + S_S)) * (S_N2O / (K_N2O + S_N2O)) * (S_O / (K_O_H5 + S_O)) * (K_I5NO / (K_I5NO + S_NO))
        
        R6  ~ 1.2 * i_N_XB * (S_NO3 / (K_INO3 + S_NO3)) * (K_I6NH / (K_I6NH + S_NH)) * (K_I6NO2 / (K_I6NO2 + S_NO2)) * (R1 + R2 + R3 + R4 + R5 + R9 + R10 + R11 + R17 + R18 - R14)
        
        R7  ~ 1.2 * i_N_XB * (S_NO2 / (K_7NO2 + S_NO2)) * (K_I7NH / (K_I7NH + S_NH)) * (R1 + R2 + R3 + R4 + R5 + R9 + R10 + R11 + R17 + R18 - R14)
        
        R8  ~ b_L_H * X_B_H
        
        R9  ~ μ_B_A1 * X_B_A1 * (S_FA / (K_FA + S_FA + S_FA^2 / K_I9FA)) * (S_O / (K_O_A1 + S_O)) * (K_I9FNA / (K_I9FNA + S_FNA))
        
        R10 ~ μ_B_A2 * X_B_A2 * (S_FNA / (K_FNA + S_FNA + S_FNA^2 / K_I10FNA)) * (S_O / (K_O_A2 + S_O)) * (K_I10FA / (K_I10FA + S_FA))
        
        R11 ~ β_g * (S_S / (K_S1 + S_S)) * R10

        R12 ~ b_L_A1 * X_B_A1

        R13 ~ b_L_A2 * X_B_A2

        R14 ~ k_a * S_NS * X_B_H

        R15 ~ k_h * X_B_H * ((X_S / X_B_H) / (K_X + X_S / X_B_H)) * (S_O / (K_O_H1 + S_O)) + η_h * (K_O_H1 / (K_O_H1 + S_O)) * (sum_NOX / (K_NO3 + sum_NOX))

        R16 ~ R15 * (X_NS / X_S)

        R17 ~ μ_H * (S_J / S_S) * X_B_H * (S_J / (K_SJ + S_J)) * (S_O / (K_O_H1 + S_O)) + η_g * (K_O_H1 / (K_O_H1 + S_O)) * (sum_NOX / (K_NO3 + sum_NOX))

        R18 ~ μ_H * (S_W / S_S) * X_B_H * (S_W / (K_SW + S_W)) * (S_O / (K_O_H1 + S_O)) + η_g * (K_O_H1 / (K_O_H1 + S_O)) * (sum_NOX / (K_NO3 + sum_NOX))

        sum_NOX ~ S_NO3 + S_NO2 + S_NO + S_N2O

        # ref [1] notation below Table 4
        A ~ (1 - Y_H * η_Y) / (1.143 * Y_H * η_Y)

        B ~ (1 - Y_H * η_Y) / (0.571 * Y_H * η_Y)

        C ~ -(i_N_XB / 14) + (1 - Y_H * η_Y) / (14 * 0.571 * Y_H * η_Y)

        D ~ -(i_N_XB / 14) - (1 / (7 * Y_A1))

        # ref [1] Table 4, ASMN component and process matrix
        reac1  ~ (1 - f_prime_D) * R8 + (1 - f_prime_D) * R12 + (1 - f_prime_D) * R13 - R15        # D(X_S)

        reac2  ~ R1 + R2 + R3 + R4 + R5 - R8 + R17 + R18                                          # D(X_B_H)

        reac3  ~ R9 - R12                                                                         # D(X_B_A1)

        reac4  ~ R10 + R11 - R13                                                                  # D(X_B_A2)

        reac5  ~ f_prime_D * R8 + f_prime_D * R12 + f_prime_D * R13                               # D(X_D)

        reac6  ~ -1 / Y_H * R1 - 1 / (Y_H * η_Y) * R2 - 1 / (Y_H * η_Y) * R3 - 1 / (Y_H * η_Y) * R4 - 1 / (Y_H * η_Y) * R5 - 1.14 * R6 - 3.43 * R7 - 1 / Y_H * R11 + R15              # D(S_S)

        reac7  ~ (1 - Y_H) / Y_H * R1 + (3.43 - Y_A1) / Y_A1 * R9 + (1.14 - Y_A2) / Y_A2 * R10 + (1 - Y_H) / Y_H * R11 + (1 - Y_H) / Y_H * R17 + (1 - Y_H) / Y_H * R18           # D(S_O)
        
        reac8  ~ -A * R2 - R6 + 1 / Y_A2 * R10                                                    # D(S_NO3)

        reac9  ~ A * R2 - B * R3 + R6 - R7 + 1 / Y_A1 * R9 - 1 / Y_A2 * R10                        # D(S_NO2)

        reac10 ~ B * R3 - B * R4                                                                   # D(S_NO)

        reac11 ~ B * R4 - B * R5                                                                   # D(S_N2O)

        reac12 ~ -i_N_XB * R1 - i_N_XB * R2 - i_N_XB * R3 - i_N_XB * R4 - i_N_XB * R5 + R7 + (-i_N_XB - 1 / Y_A1) * R9 - i_N_XB * R10 - i_N_XB * R11 + R14 - i_N_XB * R17 - i_N_XB * R18                                                     # D(S_NH)
        
        reac13 ~ -R14 + R16                                                                        # D(S_NS)
        
        reac14 ~ (i_N_XB - f_prime_D * i_N_XD) * R8 + (i_N_XB - f_prime_D * i_N_XD) * R12 + (i_N_XB - f_prime_D * i_N_XD) * R13 - R16                                       # D(X_NS)
        
        reac15 ~ -i_N_XB / 14 * R1 - i_N_XB / 14 * R2 + C * R3 - i_N_XB / 14 * R4 - i_N_XB / 14 * R5 + 1 / 7 * R8 + D * R9 - i_N_XB / 14 * R10 - i_N_XB / 14 * R11 + 1 / 14 * R14 - i_N_XB / 14 * R17 - i_N_XB / 14 * R18         # D(S_ALK)
        
        reac16 ~ -1 / Y_H * R17                                                                    # D(S_J)
       
        reac17 ~ -1 / Y_H * R18                                                                    # D(S_W)

        # ref [1] Table 4, notation below the table says "Components XI,SI, and SZ are not shown because they are not involved in any of the ASMN reactions."
        reac18 ~ 0                                                                                 # D(X_I)
        
        reac19 ~ 0                                                                                 # D(S_I)
        
        reac20 ~ 0                                                                                 # D(S_Z)
    end
end
