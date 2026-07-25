```
1   s = s0  (TIVF base, no D)
2   s(D) = s0 + s_d * D / (1 + theta_s * D)
3   s(D) = s0 + s_d * (D/D0) * exp(1 - D/D0)
4   s(D) = s0 + s_d * (D/D0) * exp(1 - D/D0)  [T0 free]
5   s(D) = s0 + s_d * (1 - exp(-D/D0))
6   s(D) = s_d * D^2 / (theta_s^2 + D^2)
7   s(D) = s_d * D^3 / (theta_s^3 + D^3)
8   s(D) = s_d * D^4 / (theta_s^4 + D^4)
9   s(D) = s0 + s_d * exp(-0.5*(x-1)^2)                          x = D/D0
10  s(D) = s0 + s_d * exp(-0.5*(x-mu_g)^2/sigma_g^2)             x = D/D0
11  s(D) = s0 + s_d * exp(-0.5*(x-mu_g)^2)                       x = D/D0
12  s(D) = s0 + s_d * exp(-0.5*(x-1)^2/sigma_g^2)                x = D/D0
13  s(D) = s0 + s_d * x / (1 + x^2)                              x = D/D0
14  s(D) = s0 + s_d * (x/mu_r)^n / (1 + (x/mu_r)^(2n))          x = D/D0
15  s(D) = s0 + s_d * (x/mu_r) / (1 + (x/mu_r)^2)               x = D/D0
16  s(D) = s0 + s_d * x^n / (1 + x^(2n))                         x = D/D0
17  s(D) = s0 + s_d * D / (D_ref/D0 + D)
18  s(D) = s0 + s_d * exp(-|ln(D/D_pk)| / w)
19  s(D) = s0 + s_d / (1 + (ln(D/D_pk)/w)^2)
20  s(D) = s0 + s_d * exp(-(ln(D/D0))^2)
21  s(D) = s0 + s_d * exp(-(ln(D/D*))^2)
22  s(D) = s0 + s_d * exp(-(lnD - lnD_pk)^2 / (2*sigma_log^2))
23  s(D) = s0 + s_d * exp(-|ln(D/D_pk)/w|^p)
24  s(D) = s0 + s_d * log(1 + D0*D/D_ref) / log2
25  s(D) = s0 + s_d * tanh(D/D0)
```
