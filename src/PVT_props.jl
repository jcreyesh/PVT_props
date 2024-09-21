
module PVT_props

# importing libraries
using DataFrames

# Gas
export PTsc_Standing, PTsc_Kessler_Lee, PTsc_Matthews_Roland    # 3
export PTsc_Kay, PTsc_Stewart_Burkhardt_Voo, PTsc_Brown_Katz_Oberfell_Alden, PTsc_Sutton, PTsc_Wichert_Aziz   # 5
export z_Sarem, z_Papay, z_Hall_Yarborough, z_Beggs_Brill, z_Dranchuk_Purvis_Robinson, z_Dranchuk_Kassem, z_Gopal # 7
export cg_Mattar_Brar, cg_Sarem, cg_Papay, cg_Hall, cg_Beggs_Brill, cg_Gopal # 6
export ug_Carr_Kobayashi_Burrows, ug_Lee_Gonzalez_Eakin # 2

# Oil
export pb_Standing, pb_Lasater, pb_Vasquez_Beggs_sp, pb_Vasquez_Beggs, pb_Glaso, pb_Total_CFP, pb_AlMarhoun, pb_Dokla_Osman, 
       pb_Petrosky_Farshad, pb_Kartoatmodjo_Schmidt_sp, pb_Kartoatmodjo_Schmidt # 9

export rs_Standing, rs_Lasater, rs_Vasquez_Beggs_sp, rs_Vasquez_Beggs, rs_Glaso, rs_Total_CFP, rs_AlMarhoun, rs_Dokla_Osman,
       rs_Petrosky_Farshad, rs_Kartoatmodjo_Schmidt_sp, rs_Kartoatmodjo_Schmidt # 9

export bo_Standing, bo_Vasquez_Beggs_sp, bo_Vasquez_Beggs, bo_Glaso, bo_Total_CFP, bo_AlMarhoun, bo_Dokla_Osman, 
       bo_Petrosky_Farshad, bo_Kartoatmodjo_Schmidt_sp, bo_Kartoatmodjo_Schmidt # 8

export bt_Glaso, bt_AlMarhoun # 2

# Py > Pb
export co_Vasquez_Beggs_sp, co_Vasquez_Beggs, co_Petrosky_Farshad, co_Kartoatmodjo_Schmidt_sp, 
       co_Kartoatmodjo_Schmidt, co_McCain_Rollins_Villena_1, co_McCain_Rollins_Villena_2, 
       co_McCain_Rollins_Villena_3 # 6

# Py <= Pb
export co_Standing, co_Vasquez_Beggs_sat_sp, co_Vasquez_Beggs_sat, co_Glaso, co_Total, co_AlMarhoun, co_Dokla_Osman, 
       co_Petrosky_Farshad_sat, co_Kartoatmodjo_Schmidt_sat # 8
       
# Py > Pb
export uob_Beal, uob_Beggs_Robinson, uob_Glaso, uob_Egbogad, uob_Kartoatmodjo_Schmidt # 5

# Py = Pb
export uos_Chew_Conally, uos_Beggs_Robinson, uos_Kartoadmodjo_Schmidt # 3

# Py < Pb
export uosb_Beal, uosb_Vasquez_Beggs, uosb_Kartoatmodjo_Schmidt # 3

# Water
export rsw_Culberson_Mcketta, rsw_McCoy  # 2 

export bw_McCain, bw_McCoy # 2

export cwb_Dodson_Standing, cwb_Osif, cws_Ramey # 3

export uw_Van_Wingen, uw_Matthews_Russell, uw_McCain, uw_McCoy # 3

export dw_McCain, tigw_Jennings_Newman # 2

# Gas
include("Gas_properties_c7+.jl")
include("Gas_properties_PTsc.jl")
include("Gas_properties_z.jl")
include("Gas_properties_cg.jl")
include("Gas_properties_ug.jl")

# Oil
include("Oil_properties_pb.jl")
include("Oil_properties_rs.jl")
include("Oil_properties_bo.jl")
include("Oil_properties_bt.jl")
include("Oil_properties_co.jl")
include("Oil_properties_uo.jl")

# Waters
include("Water_properties.jl")

end
