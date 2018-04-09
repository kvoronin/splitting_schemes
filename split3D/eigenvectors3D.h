
// p^k_DD & p^k_NN
int harmonics_W ( double * Harmon_W, int ind_W, int N, boundaryConditions Condition )
{
	if ( Condition != eNeumann && Condition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_W() \n", Condition );
		return -1;
	}

	if ( ind_W < 0 || ind_W > N )
	{
		printf ( "Wrong number ind_W = %d ( N = %d ) in harmonics_W() \n", ind_W, N );
		return -1;
	}

	for ( int i = 0; i <= N; i++ )
	{
		if ( Condition == eNeumann )
			Harmon_W[i] = sin ( ind_W * i * PI /  N );
		if ( Condition == eDirichlet )
			Harmon_W[i] = cos ( ind_W * i * PI /  N );
	}

	return 0;
}

// p^k_DD & p^k_NN
int harmonics_T ( double * Harmon_T, int ind_T, int N, boundaryConditions Condition )
{
	if ( Condition != eNeumann && Condition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_T() \n", Condition );
		return -1;
	}

	if ( ind_T < 0 || ind_T > N - 1 )
	{
		printf ( "Wrong number ind_T = %d ( N = %d ) in harmonics_T() \n", ind_T, N );
		return -1;
	}

	//double * Harmon_T = ( double * ) malloc ( N * sizeof ( double ) );

	for ( int i = 0; i < N; i++ )
	{
		if ( Condition == eNeumann )
			Harmon_T[i] = cos ( ind_T * PI * ( 2.0 * i + 1.0 ) / ( 2.0 * N ) );
		if ( Condition == eDirichlet )
			Harmon_T[i] = sin ( ind_T * PI * ( 2.0 * i + 1.0 ) / ( 2.0 * N ) );
	}

	return 0;
}

int	harmonics3D_Wx ( double * Vector_x, int ind_Wx, int ind_Ty, int ind_Tz, int Nx, int Ny, int Nz, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition )
{
	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary x-condition = %d in harmonics_W() \n", xCondition );
		return -1;
	}

	if ( ind_Wx < 0 || ind_Wx > Nx )
	{
		printf ( "Wrong number ind_Wx = %d ( Nx = %d ) in harmonics_W() \n", ind_Wx, Nx );
		return -1;
	}

	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary y-condition = %d in harmonics_T() \n", yCondition );
		return -1;
	}

	if ( ind_Ty < 0 || ind_Ty > Ny - 1 )
	{
		printf ( "Wrong number ind_Ty = %d ( Ny = %d ) in harmonics_T() \n", ind_Ty, Ny );
		return -1;
	}

	if ( zCondition != eNeumann && zCondition != eDirichlet )
	{
		printf ( "Wrong boundary z-condition = %d in harmonics_T() \n", zCondition );
		return -1;
	}

	if ( ind_Tz < 0 || ind_Tz > Nz - 1 )
	{
		printf ( "Wrong number ind_Tz = %d ( Nz = %d ) in harmonics_T() \n", ind_Tz, Nz );
		return -1;
	}

	int dim_Wx = ( Nx + 1 ) * Ny * Nz;
	double * Harmon_Wx, * Harmon_Ty, * Harmon_Tz;
	Harmon_Wx = ( double * ) malloc ( (Nx + 1) * sizeof ( double ) );
	Harmon_Ty = ( double * ) malloc ( Ny * sizeof ( double ) );
	Harmon_Tz = ( double * ) malloc ( Nz * sizeof ( double ) );

	harmonics_W ( Harmon_Wx, ind_Wx, Nx, xCondition );
	harmonics_T ( Harmon_Ty, ind_Ty, Ny, yCondition );
	harmonics_T ( Harmon_Tz, ind_Tz, Nz, zCondition );

	for ( int k = 0; k < Nz; k++ )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			for ( int i = 0; i <= Nx; i++ )
			{
				Vector_x [ k * ( Nx + 1 ) * Ny + j * ( Nx + 1 ) + i ] = Harmon_Wx [ i ] * Harmon_Ty [ j ] * Harmon_Tz [ k ];
			}
		}
	}

	free ( Harmon_Wx );
	free ( Harmon_Ty );
	free ( Harmon_Tz );

	return 0;
}


int	harmonics3D_Wy ( double * Vector_y, int ind_Tx, int ind_Wy, int ind_Tz, int Nx, int Ny, int Nz, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition )
{
	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary y-condition = %d in harmonics_W() \n", yCondition );
		return -1;
	}

	if ( ind_Wy < 0 || ind_Wy > Ny )
	{
		printf ( "Wrong number ind_Wy = %d ( Ny = %d ) in harmonics_W() \n", ind_Wy, Ny );
		return -1;
	}

	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary x-condition = %d in harmonics_T() \n", xCondition );
		return -1;
	}

	if ( ind_Tx < 0 || ind_Tx > Nx - 1 )
	{
		printf ( "Wrong number ind_Tx = %d ( Nx = %d ) in harmonics_T() \n", ind_Tx, Nx );
		return -1;
	}

	if ( zCondition != eNeumann && zCondition != eDirichlet )
	{
		printf ( "Wrong boundary z-condition = %d in harmonics_T() \n", zCondition );
		return -1;
	}

	if ( ind_Tz < 0 || ind_Tz > Nz - 1 )
	{
		printf ( "Wrong number ind_Tz = %d ( Nz = %d ) in harmonics_T() \n", ind_Tz, Nz );
		return -1;
	}

	int dim_Wy = Nx * ( Ny + 1 ) * Nz;
	double * Harmon_Wy, * Harmon_Tx, * Harmon_Tz;
	Harmon_Wy = ( double * ) malloc ( ( Ny +1 ) * sizeof ( double ) );
	Harmon_Tx = ( double * ) malloc ( Nx * sizeof ( double ) );
	Harmon_Tz = ( double * ) malloc ( Nz * sizeof ( double ) );

	harmonics_T ( Harmon_Tx, ind_Tx, Nx, xCondition );
	harmonics_W ( Harmon_Wy, ind_Wy, Ny, yCondition );
	harmonics_T ( Harmon_Tz, ind_Tz, Nz, zCondition );

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int k = 0; k < Nz; k++ )
		{
			for ( int j = 0; j <= Ny; j++ )
			{
				Vector_y [ i * (Ny + 1) * Nz + k * ( Ny + 1 ) + j ] = Harmon_Tx [ i ] * Harmon_Wy [ j ] * Harmon_Tz [ k ];
			}
		}
	}

	free (Harmon_Tx);
	free (Harmon_Wy);
	free (Harmon_Tz);

	return 0;
}

int	harmonics3D_Wz ( double * Vector_z, int ind_Tx, int ind_Ty, int ind_Wz, int Nx, int Ny, int Nz, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition )
{
	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary z-condition = %d in harmonics_W() \n", zCondition );
		return -1;
	}

	if ( ind_Wz < 0 || ind_Wz > Nz )
	{
		printf ( "Wrong number ind_Wz = %d ( Nz = %d ) in harmonics_W() \n", ind_Wz, Nz );
		return -1;
	}

	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary x-condition = %d in harmonics_T() \n", xCondition );
		return -1;
	}

	if ( ind_Tx < 0 || ind_Tx > Nx - 1 )
	{
		printf ( "Wrong number ind_Tx = %d ( Nx = %d ) in harmonics_T() \n", ind_Tx, Nx );
		return -1;
	}

	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary y-condition = %d in harmonics_T() \n", yCondition );
		return -1;
	}

	if ( ind_Ty < 0 || ind_Ty > Ny - 1 )
	{
		printf ( "Wrong number ind_Ty = %d ( Ny = %d ) in harmonics_T() \n", ind_Ty, Ny );
		return -1;
	}

	int dim_Wy = Nx * Ny * ( Nz + 1 );
	double * Harmon_Tx, * Harmon_Ty, * Harmon_Wz;

	Harmon_Tx = ( double * ) malloc ( Nx * sizeof ( double ) );
	Harmon_Ty = ( double * ) malloc ( Ny * sizeof ( double ) );
	Harmon_Wz = ( double * ) malloc ( ( Nz +1 ) * sizeof ( double ) );

	harmonics_T ( Harmon_Tx, ind_Tx, Nx, xCondition );
	harmonics_T ( Harmon_Ty, ind_Ty, Ny, yCondition );
	harmonics_W ( Harmon_Wz, ind_Wz, Nz, zCondition );

	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			for ( int k = 0; k <= Nz; k++ )
			{
				Vector_z [ j * (Nz + 1) * Nx + i * ( Nz + 1 ) + k ] = Harmon_Tx [ i ] * Harmon_Ty [ j ] * Harmon_Wz [ k ];
			}
		}
	}

	free (Harmon_Tx);
	free (Harmon_Ty);
	free (Harmon_Wz);

	return 0;
}

int	harmonics3D_T ( double * Vector_T, int ind_Tx, int ind_Ty, int ind_Tz, int Nx, int Ny, int Nz, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition)
{
	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary x-condition = %d in harmonics_T() \n", xCondition );
		return -1;
	}

	if ( ind_Tx < 0 || ind_Tx > Nx - 1 )
	{
		printf ( "Wrong number ind_Wx = %d ( Nx = %d ) in harmonics_T() \n", ind_Tx, Nx );
		return -1;
	}

	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary y-condition = %d in harmonics_T() \n", yCondition );
		return -1;
	}

	if ( ind_Ty < 0 || ind_Ty > Ny - 1 )
	{
		printf ( "Wrong number ind_Ty = %d ( Ny = %d ) in harmonics_T() \n", ind_Ty, Ny );
		return -1;
	}

	if ( zCondition != eNeumann && zCondition != eDirichlet )
	{
		printf ( "Wrong boundary z-condition = %d in harmonics_T() \n", zCondition );
		return -1;
	}

	if ( ind_Tz < 0 || ind_Tz > Nz - 1 )
	{
		printf ( "Wrong number ind_Tz = %d ( Nzy = %d ) in harmonics_T() \n", ind_Tz, Nz );
		return -1;
	}

	int dim_T = Nx * Ny * Nz;
	double * Harmon_Tx, * Harmon_Ty, * Harmon_Tz;
	Harmon_Tx = ( double * ) malloc ( Nx * sizeof ( double ) );
	Harmon_Ty = ( double * ) malloc ( Ny * sizeof ( double ) );
	Harmon_Tz = ( double * ) malloc ( Nz * sizeof ( double ) );

	harmonics_T ( Harmon_Tx, ind_Tx, Nx, xCondition );
	harmonics_T ( Harmon_Ty, ind_Ty, Ny, yCondition );
	harmonics_T ( Harmon_Tz, ind_Tz, Nz, zCondition );

	for ( int k = 0; k < Nz; k++ )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			for ( int i = 0; i < Nx; i++ )
			{
				Vector_T [ k * Nx * Ny + j * Nx + i ] = Harmon_Tx [ i ] * Harmon_Ty [ j ] * Harmon_Tz[k];
			}
		}
	}

	return 0;
}

double gamma ( int index, int N, boundaryConditions Condition )
{
	if ( index < 0 || index > N )
	{
		printf ( "Wrong number index = %d ( N = %d ) in gamma() \n", index, N );
		return -1;
	}

	double res;

	if ( Condition == eNeumann )
		res = 2.0 * sin ( index * PI / ( 2.0 * N ) );
	else
		res = - 2.0 * sin ( index * PI / ( 2.0 * N ) );
	//if ( Condition == eDirichlet )
	//	res = 2.0 * sin ( index * PI / ( 2.0 * N ) );
	//else
	//	res = 2.0 * sin ( index * PI / ( 2.0 * N ) );

	return res;
}

int myprint ( double * Vector, int dim, char * filename )
{
	FILE * f = fopen ( filename, "wt" );
	if ( f == NULL )
	{
		printf ( "Cannot open file %s for writing \n", filename );
		return -1;
	}
	else
	{
		for ( int i = 0; i < dim; i++ )
			fprintf ( f, "%f \n", Vector[i] );
		fclose ( f );
	}

	return 0;
}