# ODE models

**Base (TIVF):**

dT/dt = -βTV  
dI/dt = βTV - δI  
dV/dt = pI - cV  
dF/dt = sI - αF  

with p = p̂/(1 + εF).

**Ours (TIVFD):**

dT/dt = -βTV  
dI/dt = βTV - δI  
dV/dt = (1-d)p_eff I - cV  
dD/dt = d·p_eff I - hD  
dF/dt = s(D)I - αF  

with p_eff = p̂/(1 + ε₂F) and s(D) = s₀ + s_d·D/(1 + θ_s·D).
