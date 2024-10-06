# ---- Psc y Tsc de la mezcla de hc's ----
# Correlación de Kay para la Presión y Temperatura Pseudocríticas de la mezcla de Hc's
""" Correlación de Standing - 1977 \n
PTsc\\_Standing(Ge_gas) \n
\\- Argumentos - \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64}         \n
\t tipo -> tipo de yacimientos - (gas = "gn", gas y condensado = "gc")     \n
\t Psc -> Presión pseudocrítica - (psia) {Float64}                         \n
\t Tsc -> Temperatura pseudocrítica - (°F) {Float64}                       \n
"""
function PTsc_Standing(Ge_gas::Union{Int64, Float64}, tipo::String)::Tuple{Float64, Float64}

    if tipo == "gn"
        Psc = 677 + 15*Ge_gas - 37.5*Ge_gas^2
        Tsc = 168 + 325*Ge_gas - 12.5*Ge_gas^2
   elseif tipo == "gc"
        Psc = 706 + 517*Ge_gas - 11.1*Ge_gas^2
        Tsc = 187 + 330*Ge_gas - 71.55*Ge_gas^2
    end
    return round(Psc, digits=4), round(Tsc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Kay para la Presión y Temperatura Pseudocríticas de la mezcla de Hc's
""" Correlación de Kay - 1936 \n
PTsc\\_Kay(Pc, Tc, fm) \n
\\- Argumentos - \n
\t Pc -> Presión crítica - (psia) {AbstractVector}     \n
\t Tc -> Temperatura crítica - (°F) {AbstractVector}   \n
\t fm -> - Fracción molar (fracción) {AbstractVector}  \n
\t Psc -> Presión pseudocrítica - (psia) {Float64}     \n
\t Tsc -> Temperatura pseudocrítica - (°F) {Float64}   \n
"""
function PTsc_Kay(Pc::AbstractVector, Tc::AbstractVector, fm::AbstractVector)::Tuple{Float64, Float64}
    Psc, Tsc = 0, 0
    for i in eachindex(fm)
        Psc += Pc[i] * fm[i]
        Tsc += Tc[i] * fm[i]
    end
    return round(Psc, digits=4), round(Tsc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Stewart, W.B., Burkhardt S.F. y Voo para la presión y temperatura pseudocríticas
""" Correlación de Stewart, W.B., Burkhardt S.F. y Voo - 1959 \n
PTsc\\_Stewart\\_Burkhardt\\_Voo(Pc, Tc, fm) \n
\\- Argumentos - \n
\t Pc -> Presión crítica - (psia) {AbstractVector}     \n
\t Tc -> Temperatura crítica - (°F) {AbstractVector}   \n
\t fm -> - Fracción molar (fracción) {AbstractVector}  \n
\t Psc -> Presión pseudocrítica - (psia) {Float64}     \n
\t Tsc -> Temperatura pseudocrítica - (°F) {Float64}   \n
"""
function PTsc_Stewart_Burkhardt_Voo(Pc::AbstractVector, Tc::AbstractVector, fm::AbstractVector)::Tuple{Float64, Float64}
    J1, J2 , K = 0, 0, 0
    for i in eachindex(fm)
        J1 += fm[i]*(Tc[i]/Pc[i])
        J2 += fm[i]*sqrt(Tc[i]/Pc[i])
        K += fm[i] * (Tc[i]/sqrt(Pc[i]))
    end

    J = 1/3*J1 + 2/3*J2^2
    T = K^2/J 
    P = T/J
    Fj = 1/3*(fm[end]*Tc[end]/Pc[end]) + 2/3*(fm[end]*sqrt(Tc[end]/Pc[end]))^2
    Ej = 0.6081*Fj + 1.1325*Fj^2 - 14.004*Fj*fm[end] + 64.434*Fj*fm[end]^2
    Ek = (Tc[end]/Pc[end]^0.5) * (0.3129*fm[end] - 4.8156*fm[end]^2 + 27.3751*fm[end]^3)
    J_corr = J - Ej
    K_corr = K - Ek
    Tsc = K_corr^2/J_corr
    Psc = Tsc/J_corr

    return round(Psc, digits=4), round(Tsc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Brown, Katz, Oberfell y Alden para la Presión y Temperatura pseudocríticas
""" Correlación de Brown, G.G., Katz, D.L., Oberfell, G.G. y Alden, R.C - 1948     \n
PTsc\\_Bronw\\_Katz\\_Oberfell\\_Alden(fm\\_nhc, pm\\_m, kind)                     \n
\\- Argumentos -                                                                   \n
\t fm\\_nhc -> Fracción molar no-hidrocarburos (n2, co2, h2s) - (fracción) {AbstractVector}       \n
\t pm\\_m -> Peso molecular de la mezcla - (lbm/lbm-mol) {Float64}                 \n
\t kind -> Tipo de fluido, ng-natural gas, cg-condensate gas - (adim) {String}     \n
\t Psc\\_M -> Presión pseudocrítica de la mezcla - (psia) {Float64}                \n
\t Tsc\\_M -> Temperatura pseudocrítica de la mezcla - (°F) {Float64}              \n
"""
function PTsc_Brown_Katz_Oberfell_Alden(fm_nhc::AbstractVector, pm_m::Union{Int64,Float64}, kind::String="ng")::Tuple{Float64, Float64}
    Gem = pm_m/28.96
    Gehc = (Gem - 0.967*fm_nhc[1] - 1.52*fm_nhc[2] - 1.18*fm_nhc[3])/(1 - fm_nhc[1] - fm_nhc[2] - fm_nhc[3])
  
    if kind == "gn"
        Psc_hc = 677 + 15*Gehc - 37.5*Gehc^2
        Tsc_hc = 168 + 325*Gehc - 12.5*Gehc^2

    elseif kind == "gc"
        Psc_hc = 706 - 51.7*Gehc - 11.1*Gehc^2
        Tsc_hc = 187 + 330*Gehc - 71.5*Gehc^2
    end

    Psc_M = (1 - fm_nhc[1] - fm_nhc[2] - fm_nhc[3])*Psc_hc + 493*fm_nhc[1] + 1071*fm_nhc[2] + 1306*fm_nhc[3]
    Tsc_M = (1 - fm_nhc[1] - fm_nhc[2] - fm_nhc[3])*Tsc_hc + 227*fm_nhc[1] + 548*fm_nhc[2] + 672*fm_nhc[3]

    return round(Psc_M, digits=4), round(Tsc_M, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Sutton para la Presión y Temperatura pseudocríticas
""" Correlación de Sutton, R.P. - 1985                                         \n
PTsc\\_Sutton(fm\\_nhc, pm\\_m)                                                \n
\t fm_nhc -> Fracción molar no-hidrocarburos (n2, co2, h2s) - (fracción) {AbstractVector}    \n
\t pm_m -> Peso molecular de la mezcla - (lbm/lbm-mol) {Float64}               \n
\t Psc_M -> Presión pseudocrítica de la mezcla - (psia) {Float64}              \n
\t Tsc_M -> Temperatura pseudocrítica de la mezcla - (°F) {Float64}            \n
"""
function PTsc_Sutton(fm_nhc::AbstractVector, pm_m::Float64)::Tuple{Float64, Float64}
    Gem = pm_m/28.96
    Gehc = (Gem - 0.9672*fm_nhc[1] - 1.5196*fm_nhc[2] - 1.1767*fm_nhc[3])/(1 - fm_nhc[1] - fm_nhc[2] - fm_nhc[3])

    Psc_hc = 756.8 - 131*Gehc - 3.6*Gehc^2
    Tsc_hc = 169.2 + 349.5*Gehc - 74*Gehc^2

    Psc_M = (1 - fm_nhc[1] - fm_nhc[2] - fm_nhc[3])*Psc_hc + 493*fm_nhc[1] + 1071*fm_nhc[2] + 1306*fm_nhc[3]
    Tsc_M = (1 - fm_nhc[1] - fm_nhc[2] - fm_nhc[3])*Tsc_hc + 227*fm_nhc[1] + 548*fm_nhc[2] + 672*fm_nhc[3]
    
    return round(Psc_M, digits=4), round(Tsc_M, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Corrección de Wichert & Aziz para la Psc y Tsc de la mezcla por presencia de CO2 y H2S
""" Correlación de Wichert, E. y Aziz, K. - 1971                               \n
PTsr\\_Wichert\\_Aziz(Pc, Tc, fm, fm_nhc, Py, Ty)                              \n
\\- Argumentos -                                                               \n
\t Pc -> Presión crítica - (psia) {AbstractVector}                             \n
\t Tc -> Temperatura crítica - (°F) {AbstractVector}                           \n
\t fm -> Fracción molar hc's - (fracción) {AbstractVector}                     \n
\t fm_nhc -> Fracción molar no-hidrocarburos (n2, co2, h2s) - (fracción) {AbstractVector}     \n
\t Py -> Presión del yacimiento - (psia) {Int64, Float64}                      \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64}                    \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}               \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}             \n
"""
function PTsr_Wichert_Aziz(Pc::AbstractVector,Tc::AbstractVector,fm::AbstractVector,fm_nhc::AbstractVector, Py::Union{Int64,Float64},Ty::Union{Int64,Float64})::Tuple{Float64,Float64}
    # Factor de corrección epsilon
    eps = 120 * ((fm_nhc[1] + fm_nhc[2])^0.9 - (fm_nhc[1] + fm_nhc[2])^1.6) + 15 * (fm_nhc[2]^0.5 - fm_nhc[2]^4)
    
    # Presión y temperatura pseudocríticas
    Psc, Tsc = 0, 0
    for i in eachindex(fm)
        Psc += Pc[i] * fm[i]
        Tsc += Tc[i] * fm[i]
    end

    # Presión y temperatura pseudocríticas corregidas
    Tsc_mc = Tsc - eps
    Psc_mc = (Psc * Tsc_mc)/(Tsc + fm_nhc[2]*eps*(1 - fm_nhc[2]))

    # Presión y temperatura pseudoreducidas
    Psr = Py/Psc_mc
    Tsr = (Ty + 460)/Tsc_mc
    
    return round(Psr, digits=4), round(Tsr, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Corrección de Wichert & Aziz para la Psc y Tsc de la mezcla por presencia de CO2 y H2S
""" Correlación de Wichert, E. y Aziz, K. - 1971                                              \n
PTsc\\_Wichert\\_Aziz(Pc, Tc, fm, fm_nhc)                                                     \n
\\- Argumentos -                                                                              \n
\t Pc -> Presión crítica - (psia) {AbstractVector}                                            \n
\t Tc -> Temperatura crítica - (°F) {AbstractVector}                                          \n
\t fm -> Fracción molar hc's - (fracción) {AbstractVector}                                    \n
\t fm_nhc -> Fracción molar no-hidrocarburos (co2, h2s) - (fracción) {AbstractVector}         \n
\t Psc -> Presión pseudocrítica de la mezcla - (psia) {Float64}                               \n
\t Tsc -> Temperatura pseudocrítica de la mezcla - (°F) {Float64}                             \n
"""
function PTsc_Wichert_Aziz(Pc::AbstractVector,Tc::AbstractVector,fm::AbstractVector,fm_nhc::AbstractVector)::Tuple{Float64,Float64}
    # Factor de corrección epsilon
    eps = 120 * ((fm_nhc[1] + fm_nhc[2])^0.9 - (fm_nhc[1] + fm_nhc[2])^1.6) + 15 * (fm_nhc[2]^0.5 - fm_nhc[2]^4)
    
    # Presión y temperatura pseudocríticas
    Psc, Tsc = 0, 0
    for i in eachindex(fm)
        Psc += Pc[i] * fm[i]
        Tsc += Tc[i] * fm[i]
    end

    # Presión y temperatura pseudocríticas corregidas
    Tsc_mc = Tsc - eps
    Psc_mc = (Psc * Tsc_mc)/(Tsc + fm_nhc[2]*eps*(1 - fm_nhc[2]))
    
    return round(Psc_mc, digits=4), round(Tsc_mc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
function PTsc_Wichert_Aziz(Pc::AbstractVector,Tc::AbstractVector,fm::AbstractVector,fm_co2::Union{Int64, Float64}, fm_h2s::Union{Int64, Float64})::Tuple{Float64,Float64}
    # Factor de corrección epsilon
    eps = 120 * ((fm_co2 + fm_h2s)^0.9 - (fm_co2 + fm_h2s)^1.6) + 15 * (fm_h2s^0.5 - fm_h2s^4)
    
    # Presión y temperatura pseudocríticas
    Psc, Tsc = 0, 0
    for i in eachindex(fm)
        Psc += Pc[i] * fm[i]
        Tsc += Tc[i] * fm[i]
    end

    # Presión y temperatura pseudocríticas corregidas
    Tsc_mc = Tsc - eps
    Psc_mc = (Psc * Tsc_mc)/(Tsc + fm_h2s*eps*(1 - fm_h2s))
    
    return round(Psc_mc, digits=4), round(Tsc_mc, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Corrección de para la Psc y Tsc de la mezcla por presencia de N2 y H2O
""" Corrección para Nitrógeno y Agua                                                          \n
PTsc\\_n2\\_h2o(Psc, Tsc, fm\\_n2, fm\\_h2o)                                                  \n
\\- Argumentos -                                                                              \n
\t Psc -> Presión seudocrítica - (psia) Union{Int64, Float64}                                 \n
\t Tc -> Temperatura seudocrítica - (°F) Union{Int64, Float64}                                \n
\t fm_n2 -> Fracción molar n2 - (fracción) {Float64}                                          \n
\t fm_h2o -> Fracción molar h2o - (fracción) {Float64}                                        \n
\t Psc -> Presión pseudocrítica de la mezcla - (psia) {Float64}                               \n
\t Tsc -> Temperatura pseudocrítica de la mezcla - (°F) {Float64}                             \n
"""
function PTsc_n2_h2o(Psci::Union{Int64, Float64},Tsci::Union{Int64, Float64},fm_n2::Union{Int64, Float64}, fm_h2o::Union{Int64, Float64})::Tuple{Float64,Float64}
    Tsc_corr = -246.1*fm_n2 + 400*fm_h2o
    Psc_corr = -162*fm_n2 + 1270*fm_h2o

    Tsc = (Tsci - 227.2*fm_n2 - 1165*fm_h2o)/(1 - fm_n2 - fm_h2o) + Tsc_corr
    Psc = (Psci - 493.1*fm_n2 - 3200*fm_h2o)/(1 - fm_n2 - fm_h2o) + Psc_corr

    return Psc, Tsc
end

#---------------------------------------------------------------------------------------------------------------------------------
function PTsc_n2_h2o(Psci::Union{Int64, Float64},Tsci::Union{Int64, Float64},fm_nhc::AbstractVector)::Tuple{Float64,Float64}
    Tsc_corr = -246.1*fm_nhc[1] + 400*fm_nhc[2]
    Psc_corr = -162*fm_nhc[1] + 1270*fm_nhc[2]

    Tsc = (Tsci - 227.2*fm_nhc[1] - 1165*fm_nhc[2]/(1 - fm_nhc[1] - fm_nhc[2])) + Tsc_corr
    Psc = (Psci - 493.1*fm_nhc[1] - 3200*fm_nhc[2])/(1 - fm_nhc[1] - fm_nhc[2]) + Psc_corr

    return Psc, Tsc
end