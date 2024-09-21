# ---- Compresibilidad del gas - cg ----
# Correlación de Mattar, L., Brar, G.S., y Aziz, K. - 1975
""" Correlación de Mattar, L., Brar, G.S., y Aziz, K. - 1975 \n
cg\\_Mattar\\_Brar\\_Aziz(Psr, Tsr, z, Psc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
Psc -> Presión pseducrítica - (psia) {Float64} \n
cg -> Compresibilidad del gas - (psia^-1) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cg_Mattar_Brar_Aziz(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, z::Union{Int64, Float64}, Psc::Union{Int64, Float64})::Float64
    A1, A2, A3, A4 = 0.31506237, -1.0467099, -0.57832729, 0.53530771
    A5, A6, A7, A8 = -0.61232032, -0.10488813, 0.68157001, 0.68446549
    dr = (0.27*Psr)/(z*Tsr)
    dz = A1 + A2/Tsr + A3/Tsr^3 + 2*dr*(A4 + A5/Tsr) + (5*A5*A6*dr^4)/Tsr +
         ((2*A7*dr)*(1 + A8*dr^2 - A8^2*dr^4)*exp(-A8*dr^2))/Tsr^3

    cr = 1/Psr - 0.27/(z^2*Tsr) * (dz/(1 + (dr/z)*dz))
    cg = cr/Psc

    return round(cg, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Sarem, A.M. - 1961
""" Correlación de Sarem, A.M. - 1961 \n
cg\\_Sarem(Psr, Tsr, z, Psc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
z -> Factor de comresibilidad del gas -(adim) {Int64, Float64}    \n
Psc -> Presión pseducrítica - (psia) {Float64} \n
cg -> Compresibilidad del gas - (psia^-1) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cg_Sarem(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, z::Union{Int64, Float64}, Psc::Union{Int64, Float64})::Float64
    x = (2*Psr - 15)/14.8
    y = (2*Tsr - 4)/1.9
    
    A = [2.1433504 0.083176184 -0.021467042 -0.00087140318 0.0042846283 -0.0016595343;
    0.33123524 -0.13403614 0.066880961 -0.027174261 0.0088512291 -0.0021520929;
    0.10572871 -0.050393654 0.0050924798 0.010551336 -0.0073181933 0.0026959963;
    -0.0521840 0.0443121 -0.0193294 0.0058973 0.0015367 -0.0028327;
    0.019703980 -0.026383354 0.019262143 -0.01153539 0.0042910089 -0.0081302526;
    -0.0053095900 0.0089178330 -0.010894821 0.009559389 -0.0060114017 0.0031175170]

    Px = [0, 0.16551, 0.641002*x, 0.379221*(5*x^2-1), 0.716652*(7*x^3-3*x), 0.594225*(21*x^4-14*x^2+1)]
    Py = [0.7071068, 1.224745*y, 0.7905695*(3*y^2-1), 0.9354145*(5*y^3-3*y), 0.265165*(35*y^4-30*y^2+3), 
    0.293151*(63*y^5-70*y^3+15*y)]

    dz = 0
    for i in 1:6
        for j in 1:6
        dz += A[i,j] * Px[i] * Py[j]
        end
    end

    cr = 1/Psr - dz/z
    cg = cr/Psc

    return round(cg, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Papay,J.
""" Correlación de Papay, J. - 1968 \n
cg\\_Papay(Psr, Tsr, z, Psc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
z -> Factor de comresibilidad del gas -(adim) {Int64, Float64}    \n
Psc -> Presión pseducrítica - (psia) {Float64} \n
cg -> Compresibilidad del gas - (psia^-1) {Float64} \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cg_Papay(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, z::Union{Int64, Float64}, Psc::Union{Int64, Float64})::Float64
    dz = -3.52/(10^(0.9813*Tsr)) + 0.548*Psr/(10^(0.8157*Tsr))
    cr = 1/Psr - dz/z
    cg = cr/Psc
    return round(cg, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Hall, K.R. y Yarborough, L.
""" Correlación de Hall,K.R. y Yarborough, L. - 1973 \n
cg\\_Hall\\_Yarborough(Psr, Tsr, yi, tol, z, Psc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
yi -> Valor inicial de la variable y - (adim) {Int64, Float64}    \n
tol -> Tolerancia para el error - (adim) {Int64, Float64}         \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
Psc -> Presión pseducrítica - (psia) {Float64}                    \n
cg -> Compresibilidad del gas - (psia^-1) {Float64}               \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cg_Hall_Yarborough(Psr::Union{Int64,Float64}, Tsr::Union{Int64,Float64}, yi::Float64, tol::Union{Int64, Float64}, z::Union{Int64, Float64}, Psc::Union{Int64, Float64})
    t = 1/Tsr
    A = 0.06125*t*exp(-1.2*(1-t)^2)
    B = 14.76*t - 9.76*t^2 + 4.58*t^3
    C = 90.7*t - 242.2*t^2 + 42.4*t^3
    D = 2.18 + 2.82*t
    F = -A*Psr + (yi + yi^2 + yi^3 - yi^4)/(1 - yi)^3 - B*yi^2 + C*yi^D
    dF = (1 + 4*yi + 4*yi^2 - 4*yi^3 + yi^4)/(1 - yi)^4 - 2*B*yi + C*D*yi^(D-1)
    y = yi - F/dF
    
    while abs(F) >= tol
        A = 0.06125*t*exp(-1.2*(1-t)^2)
        B = 14.76*t - 9.76*t^2 + 4.58*t^3
        C = 90.7*t - 242.2*t^2 + 42.4*t^3
        D = 2.18 + 2.82*t
        F = -A*Psr + (y + y^2 + y^3 - y^4)/(1 - y)^3 - B*y^2 + C*y^D
        dF = (1 + 4*y + 4*y^2 - 4*y^3 + y^4)/(1 - y)^4 - 2*B*y + C*D*y^(D-1)
        y = y - F/dF
    end

    dy = (A * (1 - y)^4)/(1 + 4*y + 4*y^2 - 4*y^3 + y^4 - (1-y)^4*(2*B*y - C*D*y^(D-1)))
    dz = A/y - (A*Psr*dy)/y^2
    cr = 1/Psr - dz/z
    cg = cr/Psc

    return round(cg, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Beggs, H.D. y Brill, J.P.
""" Correlación de Beggs, H.D. y Brill, J.P. - 1974 \n
cg\\_Beggs\\_Brill(Psr, Tsr, z, Psc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
Psc -> Presión pseducrítica - (psia) {Float64}                    \n
cg -> Compresibilidad del gas - (psia^-1) {Float64}               \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cg_Beggs_Brill(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, z::Union{Int64, Float64}, Psc::Union{Int64, Float64})::Float64
    A = 1.39*(Tsr - 0.92)^0.5 - 0.36*Tsr - 0.101
    B = (0.62 - 0.23*Tsr)*Psr + (0.066/(Tsr - 0.86) - 0.037)*Psr^2 + (0.32*Psr^6/(10^(9*(Tsr-1))))
    C = 0.132 - 0.32*log10(Tsr)
    D = exp10(0.3106 - 0.49*Tsr + 0.1824*Tsr^2)

    dz = (1 - A)/(((0.62-0.23*Tsr) + (0.132/(Tsr-0.86) - 0.074)*Psr + (1.92*Psr^5)/(10^(9*(Tsr - 1))))*exp(B)) + 
          C*D*Psr^(D - 1)
    
    cr = 1/Psr - dz/z
    cg = cr/Psc
    
    return round(cg, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Gopal, V.N.
""" Correlación de Gopal, V.N. - 1977 \n
cg\\_Gopal(Psr, Tsr, z, Psc) \n
\\--------------------------------------------------------------------------------------------------------------\n
\\- Argumentos - \n
Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
Psc -> Presión pseducrítica - (psia) {Float64}                    \n
cg -> Compresibilidad del gas - (psia^-1) {Float64}               \n
\\--------------------------------------------------------------------------------------------------------------\n
"""
function cg_Gopal(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, z::Union{Int64, Float64}, Psc::Union{Int64, Float64})::Float64
    if 0.2 <= Psr < 1.2
        if 1.05 <= Tsr < 1.2
            dz = 1.6643*Tsr - 2.2114
        elseif 1.2 <= Tsr < 1.4
            dz = 0.0522*Tsr - 0.8511
        elseif 1.4 <= Tsr < 2
            dz = 0.1391*Tsr - 0.2988
        elseif 2 <= Tsr <= 3
            dz = 0.0295*Tsr - 0.0825
        end

    elseif 1.2 <= Psr < 2.8
        if 1.05 <= Tsr < 1.2
            dz = -1.3570*Tsr + 1.4942
        elseif 1.2 <= Tsr < 1.4
            dz = 0.1717*Tsr - 0.3232
        elseif 1.4 <= Tsr < 2
            dz = 0.0984*Tsr - 0.2053
        elseif 2 <= Tsr <= 3
            dz = 0.0211*Tsr - 0.0527
        end

    elseif 2.8 <= Psr < 5.4
        if 1.05 <= Tsr < 1.2
            dz = -0.3278*Tsr + 0.4752
        elseif 1.2 <= Tsr < 1.4
            dz = -0.2521*Tsr + 0.3871
        elseif 1.4 <= Tsr < 2
            dz = -0.0284*Tsr + 0.0625
        elseif 2 <= Tsr <= 3
            dz = 0.0041*Tsr + 0.0039
        end
    
    elseif 5.4 <= Psr < 15
        if 1.05 <= Tsr < 3
            dz = (0.711 + 3.66*Tsr)^-1.4667
        end    
    end

    cr = 1/Psr - dz/z
    cg = cr/Psc

    return round(cg, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
