
vec equation_motion(float, float, float, float);
mat jacobian_measurement(vec, vec);
mat equation_measurement(vec, vec);
mat inverse_measurement(vec, mat);

double measure(double);
mat inverse(mat q);
