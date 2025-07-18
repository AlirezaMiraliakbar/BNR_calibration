#################################################################
#                                                               #                             
#          This code implements the Fractionation               #
#          problem for ASM-X models. X represents               #
#          the model number                                     #
#                                                               #         
#          Author: Alireza Miraliakbar                          #
#          Date: 06-16-2025                                     #
#          Process Systems and Operations Research Lab (PSOR)   #
#          University of Connecticut                            #
#          Storrs, CT                                           #
#                                                               #
#################################################################
# using CSV, DataFrames, Distributions
# using Random

function generate_influent_df()
    # reading CSV file containing concentrations
    df = CSV.read("./Influent/UConn_Influent.csv", DataFrame)

    #removing any missing values 
    # Drop rows with any missing values
    df = dropmissing(df)



    # Initialize influent_gen DataFrame with specified column names
    influent_gen = DataFrame(Date = String[], Temp = Float64[] ,SO = Float64[], SI = Float64[], SS = Float64[], SNH = Float64[],
    SN2 = Float64[], SNO = Float64[], SHCO = Float64[], XI = Float64[], XS = Float64[], XH = Float64[], XSTO = Float64[], 
    XA = Float64[], XTS = Float64[])

    # Add rows to influent_gen DataFrame with Date and Temp from df
    for row in eachrow(df)
        push!(influent_gen, (
            Date = row.Date,
            Temp = row.Temp,
            SO = 0.0,  # These will be populated later
            SI = 0.0,
            SS = 0.0,
            SNH = 0.0,
            SN2 = 0.0,
            SNO = 0.0,
            SHCO = 0.0,
            XI = 0.0,
            XS = 0.0,
            XH = 0.0,
            XSTO = 0.0,
            XA = 0.0,
            XTS = 0.0
        ))
    end

    # sampling COD/BOD ratio from N(2.3, 0.35)
    # Sample 460 values for COD/BOD ratio from Normal distribution

    Random.seed!(42)  # Set a fixed seed for reproducible results
    COD_BOD_ratio = rand(Normal(2.3, 0.35), length(influent_gen.Date));
    # Save the COD/BOD ratios to a CSV file for later use
    ratios_df = DataFrame(Date = influent_gen.Date, COD_BOD_ratio = COD_BOD_ratio)
    CSV.write("./Influent/COD_BOD_ratios.csv", ratios_df)
    println("COD/BOD ratios saved to ./Influent/COD_BOD_ratios.csv")


    BOD_vals = df[!,"BOD (mg/L)"]
    COD_vals = COD_BOD_ratio .* BOD_vals


    # sampling COD fractions from normal distributions 
    len_data = length(influent_gen.Date)
    for i = 1 : len_data

        F_SI = rand(Normal(0.06, 0.01));
        F_SS = rand(Normal(0.26, 0.0533));
        F_XI = rand(Normal(0.39, 0.0367));
        F_XS = rand(Normal(0.28, 0.0667));

        # Ensure the fractions sum to 1
        total_fraction = F_SI + F_SS + F_XI + F_XS;
        F_SI /= total_fraction;
        F_SS /= total_fraction;
        F_XI /= total_fraction;
        F_XS /= total_fraction;

        influent_gen.SI[i] = F_SI * COD_vals[i]
        influent_gen.SS[i] = F_SS * COD_vals[i]
        influent_gen.XI[i] = F_XI * COD_vals[i]
        influent_gen.XS[i] = F_XS * COD_vals[i]

    end

    # Updating values for SNH, SNO, SHCO, and XTS
    for i in 1:length(influent_gen.Date)
        influent_gen.SNH[i] = df[i, "Ammonia (g N / L)"] * 1000;
        influent_gen.SNO[i] = df[i, "Nitrite (mg N/L)"] + df[i, "Nitrate (mg N/L)"];
        influent_gen.XTS[i] = df[i, "XTSS(g / L)"] * 1000;
        influent_gen.SHCO[i] = df[i, "Alkalanity (mol/L)"];
    end
    return influent_gen
end

function generate_influent_ss()
    # calculate the means of all columns except Date
    # Create a new DataFrame to store means
    influent_means = DataFrame()
    influent_gen = generate_influent_df()
    # Calculate means for all numeric columns (excluding Date)
    for col_name in names(influent_gen)
        if col_name != "Date"
            influent_means[!, col_name] = [mean(influent_gen[!, col_name])]
        end
    end

    # Convert DataFrame to a vector (excluding the Date column)
    influent_means_vector = []
    for col_name in names(influent_means)
        push!(influent_means_vector, influent_means[1, col_name])
    end

    # Convert to a standard Vector{Float64}
    influent_means_vector = convert(Vector{Float64}, influent_means_vector)
    return influent_means_vector
end

#TODO: add a function to read data based on seasonal variations

if abspath(PROGRAM_FILE) == @__FILE__
influent = generate_influent_ss()
end
