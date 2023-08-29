 # SEIR ode function        
 function seir_ode_log!(du, u, p, t)
    (S, E, I, R, C) = exp.(u)
    (β, γ, ν) = p
    N = S + E + I + R
  
    # infection = β * I * S / N
    infection = β * I * S / N
    
    # progression to I = γ * E
    progression= γ * E
  
    # progression to R = ν * I
    progression_R = ν * I
    
    @inbounds begin
      du[1] = -(infection) / S # S
      du[2] = (infection - progression) / E # E
      du[3] = (progression - progression_R) / I # I
      du[4] = progression_R / R # R
      du[5] = progression/C # C
    end
    nothing
  end
  