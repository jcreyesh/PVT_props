# ---- Factor de volumen de la mezcla - Bt ----
# Correlación de Glaso, O.
""" Correlación de Glaso, O. - 1980 \n
bt\\_Glaso(Py, Ty, Pb, API, Ge_gas, Co) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Bt -> Factor de volumen de la mezcla - (bl @yac/bl @std) {Float64} \n
"""
function bt_Glaso(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    E = 10^(2.8869 - (14.1811 - 3.3093*log10(Py))^0.5)
    rs = Ge_gas*(E*API^0.989/Ty^0.172)^1.2255
    F = rs*(Ty^0.5/Ge_gas^0.3)*Py^-1.1089*Ge_oil^(2.9*10^(-27e-5*rs))
    bt = 10^(8.0135e-2 + 4.7257e-1*log10(F) + 1.7351e-1*log10(F)^2)

    return round(bt, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Al-Marhoun, M.A.
""" Correlación de Al-Marhoun, M.A. - 1988 \n
bt\\_AlMarhoun(Py, Ty, API, Ge_gas) \n
\t Py -> Presión de yacimiento - (psia) {Int64, Float64} \n
\t Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
\t API -> Gravedad API - (grados) {Int64, Float64} \n
\t Ge_gas -> Gravedad específica del gas - (adim) {Int64, Float64} \n
\t Bt -> Factor de volumen de la mezcla - (bl @yac/bl @std) {Float64} \n
"""
function bt_AlMarhoun(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, API::Union{Int64, Float64}, Ge_gas::Union{Int64, Float64})::Float64

    Ge_oil = 141.5/(API + 131.5)
    rs = (185.84321*Py*(Ge_gas^1.87784)*(Ge_oil^-3.1437)*((Ty + 460)^-1.32657))^1.3984
    F = (rs^0.644516)*(Ge_gas^-1.07934)*(Ge_oil^0.724874)*(Py^-0.76191)*((Ty + 460)^2.00621)
    bt = 0.314693 + 0.106253e-4*F + 0.18883e-10*F^2

    return round(bt, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
