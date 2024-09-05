library(Rcpp)


cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')
# add works like a regular R function
add
#> function (x, y, z) 
#> .Call(<pointer: 0x107536a00>, x, y, z)
add(1, 2, 3)
#> [1] 6
#> 
#> 
cppFunction('std::map<std::string, double> calculate_rates(double beta_E_colony, double beta_I_colony,
                                              double sigma, double eta, double gamma, double mu,
                                              double zeta_to_colony, double zeta_to_sea, double psi, 
                                              double rho_to_colony, double rho_to_sea,
                                              double prop_dispersal, double prop_prospecting, int dispersal_date,
                                              int hatching_date,
                                              double S_a, double E_a, double I_a, double R_a, double D_a,
                                              double S_sea_a, double E_sea_a, double I_sea_a, double R_sea_a, double D_sea_a,
                                              double S_sea_b, double E_sea_b, double I_sea_b, double R_sea_b, double D_sea_b,
                                              double S_b, double E_b, double I_b, double R_b, double D_b,
                                              double S_a_NB, double E_a_NB, double I_a_NB, double R_a_NB, double D_a_NB,
                                              double S_sea_a_NB, double E_sea_a_NB, double I_sea_a_NB, double R_sea_a_NB, double D_sea_a_NB,
                                              double S_sea_b_NB, double E_sea_b_NB, double I_sea_b_NB, double R_sea_b_NB, double D_sea_b_NB,
                                              double S_b_NB, double E_b_NB, double I_b_NB, double R_b_NB, double D_b_NB,
                                              double S_a_N, double E_a_N, double I_a_N, double R_a_N, double D_a_N,
                                              double S_b_N, double E_b_N, double I_b_N, double R_b_N, double D_b_N) 
{
    std::map<std::string, double> rates;

    // SEIR

    // Nestlings
    // Colony A
    rates["S_a_N_to_E_a_N"] = beta_E_colony * S_a_N * (E_a + E_a_NB + E_a_N) +
                              beta_I_colony * S_a_N * (I_a + I_a_NB + I_a_N);
    rates["E_a_N_to_S_a_N"] = eta * E_a_N;
    rates["E_a_N_to_I_a_N"] = sigma * E_a_N;
    rates["I_a_N_to_R_a_N"] = gamma * I_a_N;
    rates["I_a_N_to_D_a_N"] = mu * I_a_N;

    // Colony B
    rates["S_b_N_to_E_b_N"] = beta_E_colony * S_b_N * (E_b + E_b_NB + E_b_N) +
                              beta_I_colony * S_b_N * (I_b + I_b_NB + I_b_N);
    rates["E_b_N_to_S_b_N"] = eta * E_b_N;
    rates["E_b_N_to_I_b_N"] = sigma * E_b_N;
    rates["I_b_N_to_R_b_N"] = gamma * I_b_N;
    rates["I_b_N_to_D_b_N"] = mu * I_b_N;

    // Non-Breeders
    // Colony A
    rates["S_a_NB_to_E_a_NB"] = beta_E_colony * S_a_NB * (E_a + E_a_NB + E_a_N) +
                                beta_I_colony * S_a_NB * (I_a + I_a_NB + I_a_N);
    rates["E_a_NB_to_S_a_NB"] = eta * E_a_NB;
    rates["E_a_NB_to_I_a_NB"] = sigma * E_a_NB;
    rates["I_a_NB_to_R_a_NB"] = gamma * I_a_NB;
    rates["I_a_NB_to_D_a_NB"] = mu * I_a_NB;

    // Sea A
    rates["S_sea_a_NB_to_E_sea_a_NB"] = 0;
    rates["E_sea_a_NB_to_S_sea_a_NB"] = eta * E_sea_a_NB;
    rates["E_sea_a_NB_to_I_sea_a_NB"] = sigma * E_sea_a_NB;
    rates["I_sea_a_NB_to_R_sea_a_NB"] = gamma * I_sea_a_NB;
    rates["I_sea_a_NB_to_D_sea_a_NB"] = mu * I_sea_a_NB;

    // Sea B
    rates["S_sea_b_NB_to_E_sea_b_NB"] = 0;
    rates["E_sea_b_NB_to_S_sea_b_NB"] = eta * E_sea_b_NB;
    rates["E_sea_b_NB_to_I_sea_b_NB"] = sigma * E_sea_b_NB;
    rates["I_sea_b_NB_to_R_sea_b_NB"] = gamma * I_sea_b_NB;
    rates["I_sea_b_NB_to_D_sea_b_NB"] = mu * I_sea_b_NB;

    // Colony B
    rates["S_b_NB_to_E_b_NB"] = beta_E_colony * S_b_NB * (E_b + E_b_NB + E_b_N) +
                                beta_I_colony * S_b_NB * (I_b + I_b_NB + I_b_N);
    rates["E_b_NB_to_S_b_NB"] = eta * E_b_NB;
    rates["E_b_NB_to_I_b_NB"] = sigma * E_b_NB;
    rates["I_b_NB_to_R_b_NB"] = gamma * I_b_NB;
    rates["I_b_NB_to_D_b_NB"] = mu * I_b_NB;

    // Breeders
    // Colony A
    rates["S_a_to_E_a"] = beta_E_colony * S_a * (E_a + E_a_NB + E_a_N) +
                          beta_I_colony * S_a * (I_a + I_a_NB + I_a_N);
    rates["E_a_to_S_a"] = eta * E_a;
    rates["E_a_to_I_a"] = sigma * E_a;
    rates["I_a_to_R_a"] = gamma * I_a;
    rates["I_a_to_D_a"] = mu * I_a;

    // Sea A
    rates["S_sea_a_to_E_sea_a"] = 0;
    rates["E_sea_a_to_S_sea_a"] = eta * E_sea_a;
    rates["E_sea_a_to_I_sea_a"] = sigma * E_sea_a;
    rates["I_sea_a_to_R_sea_a"] = gamma * I_sea_a;
    rates["I_sea_a_to_D_sea_a"] = mu * I_sea_a;

    // Sea B
    rates["S_sea_b_to_E_sea_b"] = 0;
    rates["E_sea_b_to_S_sea_b"] = eta * E_sea_b;
    rates["E_sea_b_to_I_sea_b"] = sigma * E_sea_b;
    rates["I_sea_b_to_R_sea_b"] = gamma * I_sea_b;
    rates["I_sea_b_to_D_sea_b"] = mu * I_sea_b;

    // Colony B
    rates["S_b_to_E_b"] = beta_E_colony * S_b * (E_b + E_b_NB + E_b_N) +
                          beta_I_colony * S_b * (I_b + I_b_NB + I_b_N);
    rates["E_b_to_S_b"] = eta * E_b;
    rates["E_b_to_I_b"] = sigma * E_b;
    rates["I_b_to_R_b"] = gamma * I_b;
    rates["I_b_to_D_b"] = mu * I_b;

    // Mobility
    // Non-Breeders
    // From colony A to sea A
    rates["S_a_NB_to_S_sea_a_NB"] = rho_to_sea * S_a_NB;
    rates["E_a_NB_to_E_sea_a_NB"] = rho_to_sea * E_a_NB;
    rates["I_a_NB_to_I_sea_a_NB"] = rho_to_sea * I_a_NB;
    rates["R_a_NB_to_R_sea_a_NB"] = rho_to_sea * R_a_NB;

    // From sea A to colony A
    rates["S_sea_a_NB_to_S_a_NB"] = rho_to_colony * S_sea_a_NB;
    rates["E_sea_a_NB_to_E_a_NB"] = rho_to_colony * E_sea_a_NB;
    rates["I_sea_a_NB_to_I_a_NB"] = rho_to_colony * I_sea_a_NB;
    rates["R_sea_a_NB_to_R_a_NB"] = rho_to_colony * R_sea_a_NB;

    // From colony B to sea B
    rates["S_b_NB_to_S_sea_b_NB"] = rho_to_sea * S_b_NB;
    rates["E_b_NB_to_E_sea_b_NB"] = rho_to_sea * E_b_NB;
    rates["I_b_NB_to_I_sea_b_NB"] = rho_to_sea * I_b_NB;
    rates["R_b_NB_to_R_sea_b_NB"] = rho_to_sea * R_b_NB;

    // From sea B to colony B
    rates["S_sea_b_NB_to_S_b_NB"] = rho_to_colony * S_sea_b_NB;
    rates["E_sea_b_NB_to_E_b_NB"] = rho_to_colony * E_sea_b_NB;
    rates["I_sea_b_NB_to_I_b_NB"] = rho_to_colony * I_sea_b_NB;
    rates["R_sea_b_NB_to_R_b_NB"] = rho_to_colony * R_sea_b_NB;

    // Breeders
    // From colony A to sea A
    rates["S_a_to_S_sea_a"] = zeta_to_sea * S_a;
    rates["E_a_to_E_sea_a"] = zeta_to_sea * E_a;
    rates["I_a_to_I_sea_a"] = zeta_to_sea * I_a;
    rates["R_a_to_R_sea_a"] = zeta_to_sea * R_a;

    // From sea A to colony A
    rates["S_sea_a_to_S_a"] = zeta_to_colony * S_sea_a;
    rates["E_sea_a_to_E_a"] = zeta_to_colony * E_sea_a;
    rates["I_sea_a_to_I_a"] = zeta_to_colony * I_sea_a;
    rates["R_sea_a_to_R_a"] = zeta_to_colony * R_sea_a;

    // From colony B to sea B
    rates["S_b_to_S_sea_b"] = zeta_to_sea * S_b;
    rates["E_b_to_E_sea_b"] = zeta_to_sea * E_b;
    rates["I_b_to_I_sea_b"] = zeta_to_sea * I_b;
    rates["R_b_to_R_sea_b"] = zeta_to_sea * R_b;

    // From sea B to colony B
    rates["S_sea_b_to_S_b"] = zeta_to_colony * S_sea_b;
    rates["E_sea_b_to_E_b"] = zeta_to_colony * E_sea_b;
    rates["I_sea_b_to_I_b"] = zeta_to_colony * I_sea_b;
    rates["R_sea_b_to_R_b"] = zeta_to_colony * R_sea_b;

    // Prospecting
    // From A to B
    rates["S_sea_a_NB_to_S_sea_b_NB"] = prop_prospecting * S_sea_a_NB;
    rates["E_sea_a_NB_to_E_sea_b_NB"] = prop_prospecting * E_sea_a_NB;
    rates["I_sea_a_NB_to_I_sea_b_NB"] = prop_prospecting * I_sea_a_NB;
    rates["R_sea_a_NB_to_R_sea_b_NB"] = prop_prospecting * R_sea_a_NB;

    return rates;
}
            
            
            
            ')
