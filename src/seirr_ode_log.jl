 # SEIRR ode function        
 function seirr_ode_log!(du, u, p, t)
    # S,E,I R1 is recovered but emitting, R2 is recovered and not emitting (done)
    (S, E, I, R1, R2, C) = exp.(u)
    (β, γ, ν, η) = p
    N = S + E + I + R1 + R2
  
    # infection = β * I * S / N
    infection = β* I* S / N
    
    # progression to I = γ * E
    progression= γ * E
  
    # progression to Re = ν * I
    progression_R1 = ν * I
  
    # progression to Rd = η * Re
    progression_R2 = η * R1
  
    @inbounds begin
      du[1] = -(infection) / S # S
      du[2] = (infection - progression) / E # E
      du[3] = (progression - progression_R1) / I # I
      du[4] = (progression_R1 - progression_R2) / R1 # R1
      du[5] = progression_R2 / R2# R2
      du[6] = infection/C # C
    end
    nothing
  end
  