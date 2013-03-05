float erff(float x)
{
	float gammp(double a, double x);

	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}
/*  eliminates an odyssey conflicting types error - and i don't use this 
    so who cares 

double erff(double x)
{
	double gammp(double a, double x);

	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}
*/
