# ---- Factor de compresibilidad z ----
# Método de Sarem
""" Correlación de Sarem, A.M. - 1961 \n
z\\_Sarem(Psr, Tsr) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Int64, Float64}      \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Int, Float64}      \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}          \n
"""
function z_Sarem(Psr::Union{Int64,Float64}, Tsr::Union{Int64,Float64})::Float64
    x = (2*Psr - 15)/14.8
    y = (2*Tsr - 4)/1.9

    A = [2.1433504 0.083176184 -0.021467042 -0.00087140318 0.0042846283 -0.0016595343;
         0.33123524 -0.13403614 0.066880961 -0.027174261 0.0088512291 -0.0021520929;
         0.10572871 -0.050393654 0.0050924798 0.010551336 -0.0073181933 0.0026959963;
         -0.0521840 0.0443121 -0.0193294 0.0058973 0.0015367 -0.0028327;
         0.019703980 -0.026383354 0.019262143 -0.01153539 0.0042910089 -0.0081302526;
         -0.0053095900 0.0089178330 -0.010894821 0.009559389 -0.0060114017 0.0031175170]

    Px = [0.7071068, 1.224745*x, 0.7905695*(3*x^2-1), 0.9354145*(5*x^3 - 3*x), 
          0.265165*(35*x^4 - 30*x^2 + 3), 0.293151*(63*x^5 - 70*x^3 + 15*x)]
    Py = [0.7071068, 1.224745*y, 0.7905695*(3*y^2-1), 0.9354145*(5*y^3 - 3*y), 
    0.265165*(35*y^4 - 30*y^2 + 3), 0.293151*(63*y^5 - 70*y^3 + 15*y)]

    z = 0
    for i in 1:6
        for j in 1:6
            z += A[i,j] * Px[i] * Py[j]
        end
    end
    
    return round(z, digits=4)

end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Papay
""" Correlación de Papay, J. - 1968 \n
z\\_Papay(Psr, Tsr) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Int64, Float64}     \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Int64, Float64}   \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
"""
function z_Papay(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64})::Float64
    z = 1 - ((3.52*Psr)/(10^(0.9813*Tsr))) + ((0.274*Psr^2)/(10^(0.8157*Tsr)))
    return round(z, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Hall, K.R. y Yarborough L.
""" Correlación de Hall, K.R. y Yarborough, L. - 1973 \n
z\\_Hall\\_Yarborough(Psr, Tsr, yi, tol) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Float64}     \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Float64}   \n
\t yi -> Valor inicial de la variable y - (adim) {Int64, Float64}    \n
\t tol -> Tolerancia para el error - (adim) {Int64, Float64}         \n
\t df -> Resumen de resultados - (DataFrame)                         \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
"""
function z_Hall_Yarborough(Psr::Union{Int64,Float64}, Tsr::Union{Int64,Float64}, yi::Float64, tol::Float64)
    y_val = Vector{Float64}(undef, 0)
    F_val = Vector{Float64}(undef, 0)
    dF_val = Vector{Float64}(undef, 0)
    
    t = 1/Tsr
    A = 0.06125*t*exp(-1.2*(1-t)^2)
    B = 14.76*t - 9.76*t^2 + 4.58*t^3
    C = 90.7*t - 242.2*t^2 + 42.4*t^3
    D = 2.18 + 2.82*t
    F = -A*Psr + (yi + yi^2 + yi^3 - yi^4)/(1 - yi)^3 - B*yi^2 + C*yi^D
    dF = (1 + 4*yi + 4*yi^2 - 4*yi^3 + yi^4)/(1 - yi)^4 - 2*B*yi + C*D*yi^(D-1)
    y = yi - F/dF

    push!(y_val, yi)
    push!(F_val, F)
    push!(dF_val, dF)
    
    while abs(F) >= tol
        A = 0.06125*t*exp(-1.2*(1-t)^2)
        B = 14.76*t - 9.76*t^2 + 4.58*t^3
        C = 90.7*t - 242.2*t^2 + 42.4*t^3
        D = 2.18 + 2.82*t
        F = -A*Psr + (y + y^2 + y^3 - y^4)/(1 - y)^3 - B*y^2 + C*y^D
        dF = (1 + 4*y + 4*y^2 - 4*y^3 + y^4)/(1 - y)^4 - 2*B*y + C*D*y^(D-1)
        y = y - F/dF
        push!(y_val, y)
        push!(F_val, F)
        push!(dF_val, dF)
    end

    z = (0.06125*Psr*t*exp(-1.2*(1 - t)^2))/y
    df = DataFrame(y = y_val, F = F_val, dF = dF_val)
    
    return df, round(z, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Beggs, H.D. y Brill, J.P.
""" Correlación de Beggs, H.D. y Brill, J.P. - 1974 \n
z\\_Beggs\\_Brill(Psr, Tsr) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Int64, Float64}     \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Int64, Float64}   \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
"""
function z_Beggs_Brill(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64})::Float64
    A = 1.39*(Tsr - 0.92)^0.5 - 0.36*Tsr - 0.101
    B = (0.62 - 0.23*Tsr)*Psr + (0.066/(Tsr - 0.86) - 0.037)*Psr^2 + (0.32*Psr^6/(10^(9*(Tsr-1))))
    C = 0.132 - 0.32*log10(Tsr)
    D = exp10(0.3106 - 0.49*Tsr + 0.1824*Tsr^2)
    z = A + (1-A)/(exp(B)) + C*Psr^D
    return round(z, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Correlación de Dranchuk, P.M., Purvis, R.A. y Robinson, D.B.
""" Correlación de Dranchuk, P.M., Purvis, R.A. y Robinson, D.B. - 1974 \n
z\\_Dranchuk\\_Purvis\\_Robinson(Psr, Tsr, zi, tol) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Int64, Float64}     \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Int64, Float64}   \n
\t zi -> Valor inicial del factor z - (adim) {Int64, Float64}    \n
\t tol -> Tolerancia para el error - (adim) {Int64, Float64}         \n
\t df -> Resumen de resultados - (DataFrame)                         \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
"""
function z_Dranchuk_Purvis_Robinson(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, zi::Union{Int64, Float64}, tol::Float64)
    z_val = Vector{Float64}(undef, 0)
    F_val = Vector{Float64}(undef, 0)
    dF_val = Vector{Float64}(undef, 0)
    dr_val = Vector{Float64}(undef, 0)
    
    A1, A2, A3, A4 = 0.31506237, -1.0467099, -0.57832729, 0.53530771
    A5, A6, A7, A8 = -0.61232032, -0.10488813, 0.68157001, 0.68446549
    dr = (0.27*Psr)/(zi*Tsr)
    A = (A1 + A2/Tsr + A3/Tsr^3)*dr
    B = (A4 + A5/Tsr)*dr^2
    C = A5*A6*dr^5/Tsr
    D = A7*(1 + A8*dr^2)*(dr^2/Tsr^3)*exp(-A8*dr^2)
    E = ((2*A7*dr^2)/(zi*Tsr^3)) * (1 + A8*dr^2 - (A8*dr^2)^2) * exp(-A8*dr^2)
    F = zi - (1 + A + B + C + D)
    dF = 1 + A/zi + 2*B/zi + 5*C/zi + E
    
    push!(z_val, zi)
    push!(F_val, F)
    push!(dF_val, dF)
    push!(dr_val, dr)

    z = zi - F/dF

    while abs(F) >= tol
        dr = (0.27*Psr)/(z*Tsr)
        A = (A1 + A2/Tsr + A3/Tsr^3)*dr
        B = (A4 + A5/Tsr)*dr^2
        C = A5*A6*dr^5/Tsr
        D = A7*(1 + A8*dr^2)*(dr^2/Tsr^3)*exp(-A8*dr^2)
        E = ((2*A7*dr^2)/(z*Tsr^3)) * (1 + A8*dr^2 - (A8*dr^2)^2) * exp(-A8*dr^2)
        F = z - (1 + A + B + C + D)
        dF = 1 + A/z + 2*B/z + 5*C/z + E
        z = z - F/dF
        
        push!(z_val, z)
        push!(F_val, F)
        push!(dF_val, dF)
        push!(dr_val, dr)
    end

    df = DataFrame(z = z_val, F = F_val, dF = dF_val, dr = dr_val)
    
    return df, round(z, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Dranchuk, P.M. y Abou-Kassem, J.H.
""" Correlación de Dranchuk, P.M. y Abou-Kassem, J.H. - 1975 \n
z\\_Dranchuk\\_Kassem(Psr, Tsr, zi, tol) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Int64, Float64}     \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Int64, Float64}   \n
\t zi -> Valor inicial del factor z - (adim) {Int64, Float64}    \n
\t tol -> Tolerancia para el error - (adim) {Int64, Float64}         \n
\t df -> Resumen de resultados - (DataFrame)                         \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}  \n
"""
function z_Dranchuk_Kassem(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64}, zi::Union{Int64, Float64}, tol::Float64)
    z_val = Vector{Float64}(undef, 0)
    F_val = Vector{Float64}(undef, 0)
    dF_val = Vector{Float64}(undef, 0)
    dr_val = Vector{Float64}(undef, 0)
    
    A1, A2, A3, A4 = 0.3265, -1.07, -0.5339, 0.01569
    A5, A6, A7, A8 = -0.05165, 0.5475, -0.7361, 0.1844
    A7, A8, A9, A10, A11 = -0.7361, 0.1844, 0.1056, 0.6134, 0.721
    dr = (0.27*Psr)/(zi*Tsr)
    A = (A1 + A2/Tsr + A3/Tsr^3 + A4/Tsr^4 + A5/Tsr^5)*dr
    B = (A6 + A7/Tsr + A8/Tsr^2)*dr^2
    C = A9*(A7/Tsr + A8/Tsr^2)*dr^5
    D = A10*(1 + A11*dr^2)*(dr^2/Tsr^3)*exp(-A11*dr^2)
    E = ((2*A10*dr^2)/(zi*Tsr^3)) * (1 + A11*dr^2 - (A11*dr^2)^2) * exp(-A11*dr^2)
    F = zi - (1 + A + B - C + D)
    dF = 1 + A/zi + 2*B/zi - 5*C/zi + E
    
    push!(z_val, zi)
    push!(F_val, F)
    push!(dF_val, dF)
    push!(dr_val, dr)

    z = zi - F/dF

    while abs(F)>= tol
        dr = (0.27*Psr)/(z*Tsr)
        A = (A1 + A2/Tsr + A3/Tsr^3 + A4/Tsr^4 + A5/Tsr^5)*dr
        B = (A6 + A7/Tsr + A8/Tsr^2)*dr^2
        C = A9*(A7/Tsr + A8/Tsr^2)*dr^5
        D = A10*(1 + A11*dr^2)*(dr^2/Tsr^3)*exp(-A11*dr^2)
        E = ((2*A10*dr^2)/(z*Tsr^3)) * (1 + A11*dr^2 - (A11*dr^2)^2) * exp(-A11*dr^2)
        F = z - (1 + A + B - C + D)
        dF = 1 + A/z + 2*B/z - 5*C/z + E
        z = z - F/dF
        
        push!(z_val, z)
        push!(F_val, F)
        push!(dF_val, dF)
        push!(dr_val, dr)
    end

    df = DataFrame(z = z_val, F = F_val, dF = dF_val, dr = dr_val)
    
    return df, round(z, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
# Método de Gopal, V.N.
""" Correlación de Gopal, V.N. - 1977 \n
z\\_Gopal(Psr, Tsr) \n
\\- Argumentos - \n
\t Psr -> Presión pseudoreducida de la mezcla - (psia) {Int64, Float64}     \n
\t Tsr -> Temperatura pseudoreducida de la mezcla - (°F) {Int64, Float64}   \n
\t z -> Factor de compresibilidad del gas - (adim) {Int64, Float64}         \n
"""
function z_Gopal(Psr::Union{Int64, Float64}, Tsr::Union{Int64, Float64})::Float64
    if 0.2 <= Psr < 1.2
        if 1.05 <= Tsr < 1.2
            z = Psr*(1.6643*Tsr - 2.2114) - 0.3467*Tsr + 1.4385
        elseif 1.2 <= Tsr < 1.4
            z = Psr*(0.0522*Tsr - 0.8511) - 0.0364*Tsr + 1.0490
        elseif 1.4 <= Tsr < 2
            z = Psr*(0.1391*Tsr - 0.2988) + 0.0007*Tsr + 0.9969
        elseif 2 <= Tsr <= 3
            z = Psr*(0.0295*Tsr - 0.0825) + 0.0009*Tsr + 0.9967
        end

    elseif 1.2 <= Psr < 2.8
        if 1.05 <= Tsr < 1.2
            z = Psr*(-1.3570*Tsr + 1.4942) + 4.6315*Tsr - 4.7009
        elseif 1.2 <= Tsr < 1.4
            z = Psr*(0.1717*Tsr - 0.3232) + 0.5869*Tsr + 0.1229
        elseif 1.4 <= Tsr < 2
            z = Psr*(0.0984*Tsr - 0.2053) + 0.0621*Tsr + 0.8580
        elseif 2 <= Tsr <= 3
            z = Psr*(0.0211*Tsr - 0.0527) + 0.0127*Tsr + 0.9549
        end

    elseif 2.8 <= Psr < 5.4
        if 1.05 <= Tsr < 1.2
            z = Psr*(-0.3278*Tsr + 0.4752) + 1.8223*Tsr - 1.9036
        elseif 1.2 <= Tsr < 1.4
            z = Psr*(-0.2521*Tsr + 0.3871) + 1.6087*Tsr - 1.6635
        elseif 1.4 <= Tsr < 2
            z = Psr*(-0.0284*Tsr + 0.0625) + 0.4714*Tsr + 0.0011
        elseif 2 <= Tsr <= 3
            z = Psr*(0.0041*Tsr + 0.0039) + 0.0607*Tsr + 0.7927
        end
    
    elseif 5.4 <= Psr < 15
        if 1.05 <= Tsr < 3
            z = Psr*(0.711 + 3.66*Tsr)^-1.4667 - 1.637/(0.319*Tsr + 0.522) + 2.071
        end    
    end

    return round(z, digits=4)
end

#---------------------------------------------------------------------------------------------------------------------------------
