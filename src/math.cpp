CubicSpline::CubicSpline(some_Eigen_vector x, some_Eigen_vector y)
{
    // compute a, b, c, d coffiecients for given input data x and y
    scipy.interpolate.CubicSpline(np.array(time), np.array(bc_values), bc_type = 'periodic')
    last here - write this function to compute the splines. Notice that I can't just solve for a, b, c, d, but I have to solve for a separate piecewise spline for every single interval in my given input data
}

double CubicSpline::evaluate(double x)
{
    return a + b * x + c * x^2 + d * x^3;
}