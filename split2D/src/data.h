#include<math.h>

extern double *xpoints, *ypoints, *zpoints;
extern const double PI;

double hx(int i);
double hy(int j);
double hz(int k);

//extern double arb_a2_coeff(double x, double y, double z);
//extern double arb_a3_coeff(double x, double y, double z);
//extern double arb_a4_coeff(double x, double y, double z);


int Nx;
int Ny;
int Nz;
int external;

double arb_a2_coeff(double x, double y, double z)
{
  return 1.0/(2.0 + cos(3*PI*x)*cos(2*PI*y));
}
double arb_a2_coeff_gradX(double x, double y, double z)
{
  return (3.0*PI)*sin(3*PI*x)*cos(2*PI*y)*arb_a2_coeff(x,y,z)*arb_a2_coeff(x,y,z);
}
double arb_a2_coeff_gradY(double x, double y, double z)
{
  return (2.0*PI)*cos(3*PI*x)*sin(2*PI*y)*arb_a2_coeff(x,y,z)*arb_a2_coeff(x,y,z);
}

double arb_a3_coeff(double x, double y, double z)
{
  if (x <= 0.5)
    return 1.0 + 0.5*sin(5*PI*x) + y*y*y;
  else
    return (1.5/(1 + (x - 0.5)*(x - 0.5))) + y*y*y;
}
double arb_a3_coeff_gradX(double x, double y, double z)
{
  if (x <= 0.5)
    return 0.5*5.0*PI*cos(5*PI*x);
  else
    return -1.5*2.0*(x - 0.5)*(1.0/(1 + (x - 0.5)*(x - 0.5)))*(1.0/(1 + (x - 0.5)*(x - 0.5)));
}
double arb_a3_coeff_gradY(double x, double y, double z)
{
  return 3*y*y;
}
double arb_a4_coeff(double x, double y, double z)
{
  if ((x <= 0.5) && (y <= 0.5))
    return 2.0 + sin(x*y*y) + 32.0*(x - 0.5)*(y - 0.5);
  else if ((x > 0.5) && (y > 0.5))
    return 2.0 + sin(x*y*y) + 8.0*(x - 0.5)*(y - 0.5);
  else
    return 2.0 + sin(x*y*y);
}

double arb_a4_coeff_gradX(double x, double y, double z)
{
  if ((x <= 0.5) && (y <= 0.5))
    return y*y*cos(x*y*y) + 32.0*(y - 0.5);
  else if ((x > 0.5) && (y > 0.5))
    return y*y*cos(x*y*y) + 8.0*(y - 0.5);
  else
    return y*y*cos(x*y*y);
}

double arb_a4_coeff_gradY(double x, double y, double z)
{
  if ((x <= 0.5) && (y <= 0.5))
    return x*2.0*y*cos(x*y*y) + 32.0*(x - 0.5);
  else if ((x > 0.5) && (y > 0.5))
    return x*2.0*y*cos(x*y*y) + 8.0*(x - 0.5);
  else
    return x*2.0*y*cos(x*y*y);
}

double density(int x, int y, int z)
{
	double temp = 0;
	switch(external)
	{
	case 1:
		if (x<Nx/2)
			return 10;
		else
			return 1;
		break;
	case 11:
		if (z<Nz/2)
			return 10;
		else
			return 1;
		break;
	default:
		return 1.0;
		break;
	}
}

double heat_conductivity(int i, int j, int k)
{
	switch(external)
	{
	case 4:
		if (j<Ny/2)
			return 1;
		else
			return 2;
    break;
  case 10:
    //nonconstant "a2" of Arbogast for smooth solutions (numSolution = 108,109)
    return arb_a2_coeff(xpoints[i] + 0.5*hx(i), ypoints[j] + 0.5*hy(j), zpoints[k] + 0.5*hz(k));
  case 11:
    //nonconstant "a3" of Arbogast for smooth solutions (numSolution = 108,109)
    return arb_a3_coeff(xpoints[i] + 0.5*hx(i), ypoints[j] + 0.5*hy(j), zpoints[k] + 0.5*hz(k));
  case 12:
    //nonconstant "a" of Arbogast for nonsmooth solutions (numSolution = 110)
    return arb_a4_coeff(xpoints[i] + 0.5*hx(i), ypoints[j] + 0.5*hy(j), zpoints[k] + 0.5*hz(k));
	default:
		return 1;
	}

}

double heat_conductivity_func(double x, double y, double z)
{
  switch(external)
  {
  case 4:
    if (y<0.5)
      return 1;
    else
      return 2;
    break;
  case 10:
    //non constant "a2" of Arbogast for smooth solutions (numSolution = 108,109)
    return arb_a2_coeff(x,y,z);
    break;
  case 11:
    //non constant "a3" of Arbogast for smooth solutions (numSolution = 108,109)
    return arb_a3_coeff(x,y,z);
    break;
  case 12:
    //non constant "a" of Arbogast for non smooth solutions (numSolution = 110)
    return arb_a4_coeff(x,y,z);
    break;
  default:
    return 1.0;
    break;
  }

}
double heat_capacity(int x, int y, int z)
{
	switch(external)
	{
	case 2:
		//true external 2
		if (y<Ny/4)
			return 1;
		else
			return 10;
		break;
	case 22:
		//external 2 like external 1
		if (y<Ny/2)
			return 10;
		else
			return 1;
		break;
	default:
		return 1;
	}
}

double hx(int i)
{
	switch(external)
	{
	//case 3:
	//	if (Nx == 8)
	//	{
	//		if (i<5)
	//			return 2.0/15;
	//		else
	//			return 1.0/9;
	//	}
	//	if (Nx == 16)
	//	{
	//		if (i<10)
	//			return 1.0/15;
	//		else
	//			return 1.0/18;
	//	}
	//	if (Nx == 32)
	//	{
	//		if (i<20)
	//			return 1.0/30;
	//		else
	//			return 1.0/36;
	//	}
	//	if (Nx == 64)
	//	{
	//		if (i<40)
	//			return 1.0/60;
	//		else
	//			return 1.0/72;
	//	}
	//	if (Nx == 128)
	//	{
	//		if (i<80)
	//			return 1.0/120;
	//		else
	//			return 1.0/144;
	//	}
	//	if (Nx == 256)
	//	{
	//		if (i<160)
	//			return 1.0/240;
	//		else
	//			return 1.0/288;
	//	}
	//	if (Nx == 512)
	//	{
	//		if (i<320)
	//			return 1.0/480;
	//		else
	//			return 1.0/576;
	//	}
	//	break;

	default:
		return 1.0/Nx;
		break;
	}

	return 0.0;
}
double hy(int j)
{
	switch(external)
	{

	case 3:
		if (Ny == 8)
		{
			if (j<5)
				return 2.0/15;
			else
				return 1.0/9;
		}
		if (Ny == 16)
		{
			if (j<10)
				return 1.0/15;
			else
				return 1.0/18;
		}
		if (Ny == 32)
		{
			if (j<20)
				return 1.0/30;
			else
				return 1.0/36;
		}
		if (Ny == 64)
		{
			if (j<40)
				return 1.0/60;
			else
				return 1.0/72;
		}
		if (Ny == 128)
		{
			if (j<80)
				return 1.0/120;
			else
				return 1.0/144;
		}
		if (Ny == 256)
		{
			if (j<160)
				return 1.0/240;
			else
				return 1.0/288;
		}
		if (Ny == 512)
		{
			if (j<320)
				return 1.0/480;
			else
				return 1.0/576;
		}
		break;

	default:
		return 1.0/Ny;
		break;
	}

	return 0.0;
}
double hz(int k)
{
	return 1.0/Nz;
}
