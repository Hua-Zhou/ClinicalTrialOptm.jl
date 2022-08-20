struct Country{T <: Real}
    m  :: T # mean of Gamma-distributed enrollment rate 
    s² :: T # var of Gamma-distributed enrollment rate 
    l  :: Integer # lower bound of number of centers
    u  :: Integer # upper bound of number of centers
    c₀ :: T # cost of initializing one center
    c  :: T # cost of running one center per unit of time
    q  :: T # cost per one enrolled patient
    d  :: T # probability of an enrolled patient dropping out
    T₀ :: Distribution{Univariate, Continuous} # distribution of center initialization time
    Td :: T # total duration of clinical trial
    end
    
    # constructor
    function Country(
        m  :: T, # mean of Gamma-distributed enrollment rate 
        s² :: T, # var of Gamma-distributed enrollment rate 
        l  :: Integer, # lower bound of number of centers
        u  :: Integer, # upper bound of number of centers
        c₀ :: T, # cost of initializing one center
        c  :: T, # cost of running one center per unit of time
        q  :: T, # cost per one enrolled patient
        d  :: T, # probability of an enrolled patient dropping out
        T₀ :: Distribution{Univariate, Continuous}, # distribution of center initialization time    
        Td :: T
    ) where T <: Real
        Country{T}(m, s², l, u, c₀, c, q, d, T₀, Td)
    end
    
    """ 
        mean(ctry :: Country)
    
    Mean number of patients enrolled by a center in country `ctry`.
    """
    mean(ctry :: Country) = ctry.m * (1 - ctry.d) * (ctry.Td - mean(ctry.T₀))
    
    """ 
        var(ctry :: Country)
    
    Variance of the number of patients enrolled by a center in country `ctry`.
    """
    var(ctry :: Country) = 
        (abs2(ctry.m) + ctry.s²) * abs2(1 - ctry.d) * var(ctry.T₀) + 
        ctry.s² * abs2(1 - ctry.d) * abs2(ctry.Td - mean(ctry.T₀)) + 
        ctry.m * (1 - ctry.d) * (ctry.Td - mean(ctry.T₀))
    
    """
        pgf(ctry, z)
    
    Probability generating function of the number of patients enrolled by one 
    center in country `ctry`.     
    """
    function pgf(
        ctry :: Country{<:Real},
        z :: Union{Real, Complex}
        )
        θ = ctry.s² / ctry.m
        d = ctry.d
        α = ctry.m^2 / ctry.s²
        T₀ = ctry.T₀
        Td = ctry.Td
        I, _ = quadgk(t -> (1 + θ * (1 - d) * (Td - t) * (1 - z))^(-α) * pdf(T₀, t), minimum(ctry.T₀), maximum(ctry.T₀))
        I
    end
    