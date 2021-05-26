function water_props = compute_water_properties_v2(T)

T_k = T + 273.15;
P = 1;


%% IAPWS IF97 thourhg XSteam for MATLAB

water_props.thermal_conductivity = XSteam('tc_pT',1,T);
water_props.specific_heat_capacity = XSteam('Cp_pT',1,T)*10^(3);
water_props.density = XSteam('rho_pT',1,T);
water_props.dynamic_viscosity = XSteam('my_pT',1,T);
water_props.kinematic_viscosity = water_props.dynamic_viscosity / water_props.density;
water_props.thermal_diffusivity = water_props.thermal_conductivity /...
    (water_props.specific_heat_capacity * water_props.density);


%% Thermal expansion

% from:  Laura Otero, Antonio D. Molina-García & Pedro D. Sanz (2002)
% Some Interrelated Thermophysical Properties of Liquid Water and Ice. I.
% A User-Friendly Modeling Review for Food High-Pressure Processing,
% Critical Reviews in Food Science and Nutrition, 42:4, 339-352,
% DOI: 10.1080/10408690290825565 

a_1 = 4.78506e1;
a_2 = -8.12847e-2;
a_3 = 8.49849e-5;
a_4 = 5.56047e5;
a_5 = -3.76355e3;
a_6 = 5.56395e0;
a_7 = 5.59682e-3;
a_8 = -2.76522e1;
a_9 = -4.28067e3;
a_10 = -3.39150e1;
a_11 = 3.65873e-1;
a_12 = -5.89617e-4;
a_13 = 3.28892e-4;
a_14 = -2.65933e-8;

PI = P + a_13 * P^2 + a_14 * P^3;
A = a_1 + a_2 * T_k + a_3 * T_k^2;
B = a_4 + a_5 * T_k + a_6 * T_k^2 + a_7 * PI * T_k + a_8 * PI;
C = a_9 + a_10 *T_k + a_11 * T_k^2 + a_12 * T_k^3;

water_props.thermal_expansion = (A + B / (C + PI)) * 10^(-4);

end