extern double *xpoints, *ypoints, *zpoints;
extern const double PI;

double hx(int i);
double hy(int j);
double hz(int k);
double density(int i, int j, int k);
double heat_conductivity(int i, int j, int k);
double heat_capacity(int i, int j, int k);
double heat_conductivity_func(double x, double y, double z);



int Nx;
int Ny;
int Nz;
int external;


double hx(int i)
{
	switch(external)
	{

	case 3:
		if (Nx == 8)
		{
			if (i<5)
				return 2.0/15;
			else
				return 1.0/9;
		}
		if (Nx == 16)
		{
			if (i<10)
				return 1.0/15;
			else
				return 1.0/18;
		}
		if (Nx == 32)
		{
			if (i<20)
				return 1.0/30;
			else
				return 1.0/36;
		}
		if (Nx == 64)
		{
			if (i<40)
				return 1.0/60;
			else
				return 1.0/72;
		}
		if (Nx == 128)
		{
			if (i<80)
				return 1.0/120;
			else
				return 1.0/144;
		}
		if (Nx == 256)
		{
			if (i<160)
				return 1.0/240;
			else
				return 1.0/288;
		}
		break;

	default:
		return 1.0/Nx;
		break;
	}

	//	return 1.0/Nx;
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
		break;

	default:
		return 1.0/Ny;
		break;
	}

	//	return 1.0/Ny;
}
double hz(int k)
{
	switch(external)
	{

	case 3:
		if (Nz == 8)
		{
			if (k<5)
				return 2.0/15;
			else
				return 1.0/9;
		}
		if (Nz == 16)
		{
			if (k<10)
				return 1.0/15;
			else
				return 1.0/18;
		}
		if (Nz == 32)
		{
			if (k<20)
				return 1.0/30;
			else
				return 1.0/36;
		}
		if (Nz == 64)
		{
			if (k<40)
				return 1.0/60;
			else
				return 1.0/72;
		}
		if (Nz == 128)
		{
			if (k<80)
				return 1.0/120;
			else
				return 1.0/144;
		}
		if (Nz == 256)
		{
			if (k<160)
				return 1.0/240;
			else
				return 1.0/288;
		}
		break;

	default:
		return 1.0/Nz;
		break;
	}
//	return 1.0/Nz;
}
double density(int i, int j, int k)
{
	double temp = 0;
	switch(external)
	{
		/*	case 1:
		if (j<Ny/4)
		return 1;
		else
		return 10;
		*/
	case 1:
		if (i<Nx/2)
			return 1000;
		else
			return 1;
		break;
	default:
		return 1;
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
	default:
		return 1;
	}

}

double heat_capacity(int i, int j, int k)
{

	switch(external)
	{
	case 2:
		if (k<Nz/4)
			return 1;
		else
			return 10;
		break;
	default:
		return 1;
	}

	/*
	if (i<Nx/2)
	return 10;
	else
	return 1;
	*/
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
  default:
    return 1.0;
    break;
  }

}