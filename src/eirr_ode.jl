# EIRRC ODE not on log scale
function eirr_ode!(du, u, p, t)
    # S,E,I R1 is recovered but emitting, R2 is recovered and not emitting (done)
    (E, I, R1, R2, C) = u
    (alpha, γ, ν, η) = p
  
    # infection = alpha * I
    infection = alpha * I
    
    # progression to I = γ * E
    progression= γ * E
  
    # progression to Re = ν * I
    progression_R1 = ν * I
  
    # progression to Rd = η * Re
    progression_R2 = η * R1
  
    @inbounds begin
      du[1] = (infection - progression) # E
      du[2] = (progression - progression_R1)# I
      du[3] = (progression_R1 - progression_R2) # R1
      du[4] = progression_R2 # R2
      du[5] = progression # keeping track of cumulative cases
    end
    nothing
  end
  
  