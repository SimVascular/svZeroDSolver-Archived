// double cubic_spline(double a, double b, double c, double d, double x) 
// {
//     return a + b * x + c * x^2 + d * x^3;
// }
// 
// function_pointer create_cubic_spline()
// {
// 
// }

class CubicSpline
{
    double a;
    double b;
    double c;
    double d;
public:
    // create the spline by solving for the coefficients comprising the spline, for the given input data
    CubicSpline(const some_Eigen_vector& x, const some_Eigen_vector& y);
    
    // evaluate the spline
    double evaulate(double x); // todo: need to make this a function template or whatever, so that this function can work with either int, float or double types // or for now, I can just make this function work only for double inputs and then in the future, if a user ever wants to make this function work for multiple data types, then we can make that change
    
    // then for UnsteadyFlowRef block, Qfunc can just be a pointer to an instance of the CubicSpline class and to evaluate the flow rate at a given time, just call Qfunc->evaluate
    
    last here - clean up this class
}