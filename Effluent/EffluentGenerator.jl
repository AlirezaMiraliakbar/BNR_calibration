#################################################################
#                                                               #                             
#          This code implements the Fractionation               #
#          procedure for effluent data for ASM-X models.        #
#          X represents the model number                        #
#                                                               #         
#          Author: Alireza Miraliakbar                          #
#          Date: 07-07-2025                                     #
#          Process Systems and Operations Research Lab (PSOR)   #
#          University of Connecticut                            #
#          Storrs, CT                                           #
#                                                               #
#################################################################

using CSV, DataFrames, Distributions
using Random, Dates

function generate_effluent_ss()
    
    # reading CSV file containing concentrations
    df = CSV.read("./Effluent/UConn_Effluent_converted.csv", DataFrame)

    # computing COD from BOD using the previously saved COD/BOD ratios
    ratios_df = CSV.read("./Influent/COD_BOD_ratios.csv", DataFrame)

    

    #removing any missing values 
    # Drop rows with any missing values
    df = dropmissing(df)

    # picking dates that are common in both ratios_df and df
    common_dates = intersect(ratios_df.Date, df.Date)
    println("Number of common dates: ", length(common_dates))
    
    # Filter both dataframes to only include common dates
    df = df[in.(df.Date, Ref(common_dates)), :]
    ratios_df = ratios_df[in.(ratios_df.Date, Ref(common_dates)), :]
    
    # Sort both dataframes by date to ensure alignment
    sort!(df, :Date)
    sort!(ratios_df, :Date)
    
    println("Number of rows in filtered ratios_df: ", nrow(ratios_df))
    println("Number of rows in filtered df: ", nrow(df))
    
    # Verify that both dataframes have the same number of rows after filtering
    if nrow(ratios_df) != nrow(df)
        error("After filtering, the number of rows in ratios_df ($(nrow(ratios_df))) does not match the number of rows in df ($(nrow(df)))")
    end

    

    # Initialize effluent_gen DataFrame with specified column names
    effluent_gen = DataFrame(Date = String[], Temp = Float64[], SHCO = Float64[], SNH = Float64[], SNO= Float64[], COD = Float64[], TKN = Float64[], XTS = Float64[])

    # Add rows to effluent_gen DataFrame with Date and Temp from df
    for row in eachrow(df)
        push!(effluent_gen, (
            Date = row.Date,
            Temp = row.Temp,
            SHCO = 0.0,  # These will be populated later
            SNH = 0.0,
            SNO = 0.0,
            COD = 0.0,
            TKN = 0.0,
            XTS = 0.0
        ))
    end

    
    BOD_vals = df[!,"BOD (mg/L)"]
    COD_vals = ratios_df.COD_BOD_ratio .* BOD_vals

    effluent_gen.COD = COD_vals
    effluent_gen.TKN = df[!, "TKN(mg /L)"]

    # Updating values for SNH, SNO, SHCO, and XTS
    for i in 1:length(effluent_gen.Date)
        effluent_gen.SNH[i] = df[i, "Ammonia (g N / L)"] .* 1000;
        effluent_gen.SNO[i] = df[i, "Nitrite (mg N/L)"] + df[i, "Nitrate (mg N/L)"];
        effluent_gen.XTS[i] = df[i, "XTSS(mg / L)"];
        effluent_gen.SHCO[i] = parse(Float64, strip(df[i, "Alkalanity (mol/L)"]));
    end

    # calculate the means of all columns except Date
    # Create a new DataFrame to store means
    effluent_means = DataFrame()

    # Calculate means for all numeric columns (excluding Date)
    for col_name in names(effluent_gen)
        if col_name != "Date"
            effluent_means[!, col_name] = [mean(effluent_gen[!, col_name])]
        end
    end

    # Convert DataFrame to a vector (excluding the Date column)
    effluent_means_vector = []
    for col_name in names(effluent_means)
        push!(effluent_means_vector, effluent_means[1, col_name])
    end

    # Convert to a standard Vector{Float64}
    effluent_means_vector = convert(Vector{Float64}, effluent_means_vector)
    return effluent_means_vector
end

# Call the function to generate the effluent data
generate_effluent_ss()

