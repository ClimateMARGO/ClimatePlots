function DICE_EBM(t, Δt, F; physics=nothing, default="DICE-2013r", return_params=false)
    ζ = [0., 0., 0., 0.]
    if isnothing(physics)
        if default == "Geoffroy-2013"
            physics = Dict("Δt" => Δt, "Cu" => 7.3, "Cd" => 106., "κ" => 0.73, "B" => 1.13)
        elseif default == "Geoffroy-NoMixing"
            physics = Dict("Δt" => Δt, "Cu" => 7.3, "Cd" => 106., "κ" => 0.0, "B" => 1.13)
        elseif default == "MARGO"
            physics = Dict("Δt" => Δt, "Cu" => 0., "Cd" => 106., "κ" => 0.73, "B" => 1.13)
        elseif default == "DICE-2013r"
            ζ = [0.098, 3.8/2.9, 0.088, 0.025]
            physics = Dict("Δt" => Δt, "Cu" => 5. /ζ[1], "Cd" =>  ζ[3]*5. /ζ[4], "κ" => ζ[3], "B" => ζ[2])
        elseif default == "DICE-2016R"
            ζ = [0.1005, 3.6813/3.1, 0.088, 0.025]
            physics = Dict("Δt" => Δt, "Cu" => 5. /ζ[1], "Cd" =>  ζ[3]*5. /ζ[4], "κ" => ζ[3], "B" => ζ[2])
        elseif default == "DICE-2013r-NoMixing"
            ζ = [0.098, 3.8/2.9, 0.088, 0.025]
            ζ[3] = 0.
            physics = Dict("Δt" => Δt, "Cu" => 5. /ζ[1], "Cd" =>  ζ[3]*5. /ζ[4], "κ" => ζ[3], "B" => ζ[2])
        end
    end
    
    ζ[1] = Δt/physics["Cu"]
    ζ[2] = physics["B"]
    ζ[3] = physics["κ"]
    ζ[4] = physics["κ"]*physics["Δt"]/physics["Cd"]
    
    T = zeros(size(t))
    T_LO = zeros(size(t))
    
    for i = 1:length(t)-1
        if ζ[1] == Inf
            T[i+1] = (F[i] + ζ[3]*T_LO[i]) / (ζ[2] + ζ[3])
        else
            T[i+1] = T[i] + ζ[1]*(F[i] - ζ[2]*T[i] - ζ[3]*(T[i] - T_LO[i]))
        end
        T_LO[i+1] = T_LO[i] + ζ[4]*(T[i] - T_LO[i])
    end
    
    if return_params
        return physics
    else
        return T
    end
end;