# ---- Viscosidad del gas - ug ----
# Correlación de Carr, N.L., Kobayachi, R.y Burrows, D.B.
""" Correlación de Carr, N.L., Kobayachi, R. y Burrows, D.B. - 1954 \n
ug\\_Carr\\_Kobayashi\\_Burrows(Ty, Psr, Tsr, M, fm, fm_nhc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
Psr -> Presión pseudoreducida - (psia) {Int64, Float64} \n
Tsr -> Temperatura pseudoreducida - (°F) {Int64, Float64} \n
M -> Pesos moleculares de los componentes de la mezcla - (lbm/lbm-mol) {AbstractVector}
fm -> Fracciones molares de los componentes de la mezcla - (fracción) {AbstractVector} 
fm_nhc -> Fracciones molares de los componentes no-hc's de la mezcla - (fracción) {AbstractVector}
ug -> Viscosidad del gas - (cp) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function ug_Carr_Kobayashi_Burrows(Ty::Union{Int64, Float64}, Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, M::AbstractVector, fm::AbstractVector, fm_nhc::AbstractVector)
    a0, a1, a2, a3= -2.46211820, 2.97054714, -2.86264054e-1, 8.05420522e-3
    a4, a5, a6, a7 = 2.80860949, -3.49803305, 3.60373020e-1, -1.04432413e-2
    a8, a9, a10, a11 = -7.93385684e-1, 1.39643306, -1.49144925e-1, 4.41015512e-3
    a12, a13, a14, a15 = 8.39387178e-2, -1.86408848e-1, 2.03367881e-2, -6.09579263e-4
    
    Mm = sum(M .* fm)
    Ge = Mm/28.96
    ug1 = (1.709e-5 - 2.062e-6*Ge)*Ty + 8.188e-3 - 6.15e-3*log10(Ge)

    C_co2 = fm_nhc[1]*(9.08e-3*log10(Ge) + 6.24e-3)
    C_h2s = fm_nhc[2]*(8.49e-3*log10(Ge) + 3.73e-3)
    C_n2 = fm_nhc[3]*(8.48e-3*log10(Ge) + 9.59e-3)
    ug1c = ug1 + C_co2 + C_h2s + C_n2

    ln = a0 + a1*Psr + a2*Psr^2 + a3*Psr^3 + Tsr*(a4 + a5*Psr + a6*Psr^2 + a7*Psr^3) + 
         Tsr^2*(a8 + a9*Psr + a10*Psr^2 + a11*Psr^3) + Tsr^3*(a12 + a13*Psr + a14*Psr^2 + a15*Psr^3)

    ug = ug1c*exp(ln)/Tsr

    return round(ug, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Lee, A.L., González, M.H. y Eakin, B.E.
""" Correlación de Lee, A.L., González, M.H. y Eakin, B.E. - 1966 \n
ug\\_Lee\\_Gonzalez\\_Eakin(Py, Ty, z, M, fm) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Py - Presión de yacimiento - (psia) {Int64, Float64}
Ty -> Temperatura del yacimiento - (°F) {Int64, Float64} \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}
M -> Pesos moleculares de los componentes de la mezcla - (lbm/lbm-mol) {AbstractVector}
fm -> Fracciones molares de los componentes de la mezcla - (fracción) {AbstractVector} 
ug -> Viscosidad del gas - (cp) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function ug_Lee_Gonzalez_Eakin(Py::Union{Int64, Float64}, Ty::Union{Int64, Float64}, z::Union{Int64, Float64}, M::AbstractVector, fm::AbstractVector)
    Ty = Ty + 460
    Mm = sum(M .* fm)
    dg = 1.4935e-3*(Py*Mm)/(z*Ty)
    K = (9.4 + 0.02*Mm)*Ty^1.5/(209 + 19*Mm + Ty)
    X = 3.5 + 986/Ty + 0.01*Mm
    Y = 2.4 - 0.2*X
    ug = K*exp(X*dg^Y)/10^4
    return round(ug, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
