double scal_l2_bnd ( char * TorW, char * XorY, double * Vec1, double * Vec2, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition  );

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

	//double * Harmon_W = ( double * ) malloc ( (N + 1) * sizeof ( double ) );

	for ( int i = 0; i <= N; i++ )
	{
		if ( Condition == eNeumann )
			Harmon_W[i] = sin ( ind_W * i * PI /  N );
//			Harmon_Wx[i] = sin ( ind_Wx * PI * ( 2.0 * i + 1.0 ) / ( 2.0 * Nx ) );
		if ( Condition == eDirichlet )
			Harmon_W[i] = cos ( ind_W * i * PI /  N );
//			Harmon_Wx[i] = cos ( ind_Wx * PI * ( 2.0 * i + 1.0 ) / ( 2.0 * Nx ) );
	}

	//* Harmon_W_pt = Harmon_W;
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

	//* Harmon_T_pt = Harmon_T;
	return 0;
}

int	harmonics2D_Wx ( double * Vector_x, int ind_Wx, int ind_Ty, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_W() \n", xCondition );
		return -1;
	}

	if ( ind_Wx < 0 || ind_Wx > Nx )
	{
		printf ( "Wrong number ind_Wx = %d ( Nx = %d ) in harmonics_W() \n", ind_Wx, Nx );
		return -1;
	}

	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_T() \n", yCondition );
		return -1;
	}

	if ( ind_Ty < 0 || ind_Ty > Ny - 1 )
	{
		printf ( "Wrong number ind_Ty = %d ( Ny = %d ) in harmonics_T() \n", ind_Ty, Ny );
		return -1;
	}

	int dim_Wx = ( Nx + 1 ) * Ny;
	double * Harmon_Wx, * Harmon_Ty;
	Harmon_Wx = ( double * ) malloc ( (Nx + 1) * sizeof ( double ) );
	Harmon_Ty = ( double * ) malloc ( Ny * sizeof ( double ) );

	harmonics_W ( Harmon_Wx, ind_Wx, Nx, xCondition );
	harmonics_T ( Harmon_Ty, ind_Ty, Ny, yCondition );

	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i <= Nx; i++ )
		{
			Vector_x [ j * ( Nx + 1 ) + i ] = Harmon_Wx [ i ] * Harmon_Ty [ j ];
		}
	}

	return 0;
}


int	harmonics2D_Wy ( double * Vector_y, int ind_Wy, int ind_Tx, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_W() \n", xCondition );
		return -1;
	}

	if ( ind_Wy < 0 || ind_Wy > Ny )
	{
		printf ( "Wrong number ind_Wy = %d ( Ny = %d ) in harmonics_W() \n", ind_Wy, Ny );
		return -1;
	}

	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_T() \n", yCondition );
		return -1;
	}

	if ( ind_Tx < 0 || ind_Tx > Nx - 1 )
	{
		printf ( "Wrong number ind_Tx = %d ( Nx = %d ) in harmonics_T() \n", ind_Tx, Nx );
		return -1;
	}

	int dim_Wy = Nx * ( Ny + 1 );
	double * Harmon_Wy, * Harmon_Tx;
	Harmon_Wy = ( double * ) malloc ( ( Ny +1 ) * sizeof ( double ) );
	Harmon_Tx = ( double * ) malloc ( Nx * sizeof ( double ) );

	harmonics_W ( Harmon_Wy, ind_Wy, Ny, yCondition );
	harmonics_T ( Harmon_Tx, ind_Tx, Nx, xCondition );

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j <= Ny; j++ )
		{
			Vector_y [ i * ( Ny + 1 ) + j ] = Harmon_Wy [ j ] * Harmon_Tx [ i ];
		}
	}

	return 0;
}
int	harmonics2D_T ( double * Vector_T, int ind_Tx, int ind_Ty, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	if ( xCondition != eNeumann && xCondition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_T() \n", xCondition );
		return -1;
	}

	if ( ind_Tx < 0 || ind_Tx > Nx - 1 )
	{
		printf ( "Wrong number ind_Wx = %d ( Nx = %d ) in harmonics_T() \n", ind_Tx, Nx );
		return -1;
	}

	if ( yCondition != eNeumann && yCondition != eDirichlet )
	{
		printf ( "Wrong boundary condition = %d in harmonics_T() \n", yCondition );
		return -1;
	}

	if ( ind_Ty < 0 || ind_Ty > Ny - 1 )
	{
		printf ( "Wrong number ind_Ty = %d ( Ny = %d ) in harmonics_T() \n", ind_Ty, Ny );
		return -1;
	}

	int dim_T = Nx * Ny;
	double * Harmon_Tx, * Harmon_Ty;
	Harmon_Tx = ( double * ) malloc ( Nx * sizeof ( double ) );
	Harmon_Ty = ( double * ) malloc ( Ny * sizeof ( double ) );

	harmonics_T ( Harmon_Tx, ind_Tx, Nx, xCondition );
	harmonics_T ( Harmon_Ty, ind_Ty, Ny, yCondition );

	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			Vector_T [ j * Nx + i ] = Harmon_Tx [ i ] * Harmon_Ty [ j ];
		}
	}

	return 0;
}






int	coeffs_Wx ( double * Coeff_Vec_x, double * Vector_x, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	int dim_Wx = ( Nx + 1 ) * Ny;
	double norm_sq;

	double * Harmonic_ij = (double *) malloc ( dim_Wx * sizeof ( double ) );
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i <= Nx; i++ )
		{
			Coeff_Vec_x [ j * ( Nx + 1 ) + i ] = 0.0;

			harmonics2D_Wx ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
			
			Coeff_Vec_x[ j * ( Nx + 1 ) + i ] = scal_l2_bnd ( "W", "X", Vector_x, Harmonic_ij, Nx, Ny, xCondition, yCondition  );
			norm_sq = scal_l2_bnd ( "W", "X", Harmonic_ij, Harmonic_ij, Nx, Ny, xCondition, yCondition  );

			//norm_sq = 0.0;
			//for ( int k = 0; k < dim_Wx; k++ )
			//{
			//	Coeff_Vec_x [ j * ( Nx + 1 ) + i ] += Vector_x [ k ] * Harmonic_ij [ k ];
			//	norm_sq += Harmonic_ij [ k ] * Harmonic_ij [ k ];
			//}

			if ( fabs ( norm_sq ) < 1.0e-14 )
				;
			else
				Coeff_Vec_x [ j * ( Nx + 1 ) + i ] /= norm_sq;
		}
	}

	free ( Harmonic_ij );

	return 0;
}

int	coeffs_Wy ( double * Coeff_Vec_y, double * Vector_y, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	int dim_Wy = Nx * ( Ny + 1 );
	double norm_sq;

	double * Harmonic_ij = (double *) malloc ( dim_Wy * sizeof ( double ) );
	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j <= Ny; j++ )
		{
			Coeff_Vec_y [ i * ( Ny + 1 ) + j ] = 0.0;

			harmonics2D_Wy ( Harmonic_ij, j, i, Nx, Ny, xCondition, yCondition );
			
			norm_sq = 0.0;
			//printf ( "i = %d j = %d \n", i, j );

			Coeff_Vec_y[ i * ( Ny + 1 ) + j ] = scal_l2_bnd ( "W", "Y", Vector_y, Harmonic_ij, Nx, Ny, xCondition, yCondition  );
			norm_sq = scal_l2_bnd ( "W", "Y", Harmonic_ij, Harmonic_ij, Nx, Ny, xCondition, yCondition  );

			//for ( int k = 0; k < dim_Wy; k++ )
			//{
			//	Coeff_Vec_y [ i * ( Ny + 1 ) + j ] += Vector_y [ k ] * Harmonic_ij [ k ];
			//	printf ( "k = %d vec[%d] = %f harm[%d] = %f \n", k, k, Vector_y[k], k, Harmonic_ij[k] );
			//	norm_sq += Harmonic_ij [ k ] * Harmonic_ij [ k ];
			//}
			//printf ( "(vec_y, harm_ij) = %f norm_sq = %f \n", Coeff_Vec_y [ i * ( Ny + 1 ) + j ], norm_sq );

			if ( fabs ( norm_sq ) < 1.0e-14 )
				;
			else
				Coeff_Vec_y [ i * ( Ny + 1 ) + j ] /= norm_sq;
		}
	}

	free ( Harmonic_ij );

	return 0;
}

int	coeffs_T ( double * Coeff_Vec_T, double * Vector_T, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	int dim_T = Nx * Ny;
	double norm_sq;
	char cdummy;

	double * Harmonic_ij = (double *) malloc ( dim_T * sizeof ( double ) );
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			Coeff_Vec_T [ j * Nx + i ] = 0.0;

			harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );

			Coeff_Vec_T[ j * Nx + i ] = scal_l2_bnd( "T", &cdummy, Vector_T, Harmonic_ij, Nx, Ny, xCondition, yCondition  );
			norm_sq = scal_l2_bnd ( "T", &cdummy, Harmonic_ij, Harmonic_ij, Nx, Ny, xCondition, yCondition  );

			
			//norm_sq = 0.0;
			//for ( int k = 0; k < dim_T; k++ )
			//{
			//	Coeff_Vec_T [ j * Nx + i ] += Vector_T [ k ] * Harmonic_ij [ k ];
			//	norm_sq += Harmonic_ij [ k ] * Harmonic_ij [ k ];
			//}

			if ( norm_sq == 0 && Coeff_Vec_T [ j * Nx + i ] == 0 )
				;
			else
				Coeff_Vec_T [ j * Nx + i ] /= norm_sq;
		}
	}

	free ( Harmonic_ij );

	return 0;
}


int vector_Wx ( double * Vector_x, double * Coeff_Vec_x, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i <= Nx; i++ )
		{
			Vector_x [j * (Nx + 1) + i] = 0.0;
		}
	}

	int dim_Wx = ( Nx + 1 ) * Ny;
	double * Harmonic_ij = (double *) malloc ( dim_Wx * sizeof ( double ) );
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i <= Nx; i++ )
		{
			harmonics2D_Wx ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
			for ( int k = 0; k < dim_Wx; k++ )
				Vector_x [ k ] += Coeff_Vec_x [ j * ( Nx + 1 ) + i ] * Harmonic_ij [ k ];
		}
	}
		
	free ( Harmonic_ij );
	return 0;
}


int vector_Wy ( double * Vector_y, double * Coeff_Vec_y, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j <= Ny; j++ )
		{
			Vector_y [i * (Ny + 1) + j] = 0.0;
		}
	}

	int dim_Wy = Nx * (Ny + 1);
	double * Harmonic_ij = (double *) malloc ( dim_Wy * sizeof ( double ) );
	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j <= Ny; j++ )
		{
			//printf ( "i = %d j = %d \n", i, j );
			harmonics2D_Wy ( Harmonic_ij, j, i, Nx, Ny, xCondition, yCondition );
			for ( int k = 0; k < dim_Wy; k++ )
			{
				Vector_y [ k ] += Coeff_Vec_y [ i * ( Ny + 1 ) + j ] * Harmonic_ij [ k ];
				//printf ( "Coeff_ij = %f Harmonic_ij[%d] = %f \n", Coeff_Vec_y [ i * (Ny + 1) + j ], k, Harmonic_ij [ k ] );
			}
			//_getch();
		}
	}
		
	free ( Harmonic_ij );
	return 0;
}

int vector_T ( double * Vector_T, double * Coeff_Vec_T, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			Vector_T [j * Nx + i] = 0.0;
		}
	}

	int dim_T = Nx * Ny;
	double * Harmonic_ij = (double *) malloc ( dim_T * sizeof ( double ) );
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			//printf ( "i = %d j = %d \n", i, j );
			harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
			for ( int k = 0; k < dim_T; k++ )
			{
				Vector_T [ k ] += Coeff_Vec_T [ j * Nx + i ] * Harmonic_ij [ k ];
				//printf ( "Coeff_ij = %f Harmonic_ij[%d] = %f \n", Coeff_Vec_T [ j * Nx + i ], k, Harmonic_ij [ k ] );
			}
			//_getch();
		}
	}
		
	free ( Harmonic_ij );
	return 0;
}


double gamma ( int index, int N, boundaryConditions Condition )
//double gamma ( int index, int N )
{
	//if ( Condition != eNeumann && Condition != eDirichlet )
	//{
	//	printf ( "Wrong boundary condition = %d in gamma() \n", Condition );
	//	return -1;
	//}

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

int	coeffs_Aminus1B ( double * Coeffs_Wx0, double * Coeffs_Wy0, double * Coeffs_T0, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	if ( Coeffs_T0 == NULL )
	{
		printf ( "Bad input coeffs vector_T in coeffs_Aminus1B() \n" );
		return -1;
	}

	double gamma_i, gamma_j, koeff_i, koeff_j, h_i, h_j;

	h_i = 1.0 / Nx;
	for ( int j = 0; j < Ny; j++ )
	{
		//gamma_j = gamma ( j, Ny, yCondition );
		for ( int i = 0; i <= Nx; i++ )
		{
			if ( i < Nx )
			{
				gamma_i = gamma ( i, Nx, xCondition );
				koeff_i = 6.0 * gamma_i / ( (6.0 - gamma_i * gamma_i ) * h_i );
			}
			else 
				koeff_i = 0.0;
			Coeffs_Wx0 [ j * ( Nx + 1 ) + i ] = Coeffs_T0 [ j * Nx + i ] * koeff_i ;
		}
	}

	h_j = 1.0 / Ny;
	for ( int i = 0; i < Nx; i++ )
	{
		//gamma_i = gamma ( i, Nx, xCondition );
		for ( int j = 0; j <= Ny; j++ )
		{
			if ( j < Ny )
			{
				gamma_j = gamma ( j, Ny, yCondition );
				koeff_j = 6.0 * gamma_j / ( (6.0 - gamma_j * gamma_j ) * h_j );
			}
			else
				koeff_j = 0.0;
			Coeffs_Wy0 [ i * ( Ny + 1 ) + j ] = Coeffs_T0 [ j * Nx + i ] * koeff_j ;
		}
	}

	return 0;
}


int	coeffs_Btr ( double * Coeffs_Btr, double * Coeff_Vec_x, double * Coeff_Vec_y, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	if ( Coeff_Vec_x == NULL || Coeff_Vec_y == NULL )
	{
		printf ( "Bad input coeffs' vectors_x,y in coeffs_Btr() \n" );
		return -1;
	}

	int dim_T = Nx * Ny;
	double gamma_i, gamma_j;

	for ( int j = 0; j < Ny; j++ )
	{
		gamma_j = gamma ( j, Ny, yCondition );
		//gamma_j = gamma ( j, Ny );
		for ( int i = 0; i < Nx; i++ )
		{
			gamma_i = gamma ( i, Nx, xCondition );
			//gamma_i = gamma ( i, Nx );
			Coeffs_Btr [ j * Nx + i ] = Coeff_Vec_y [ i * ( Ny + 1 ) + j ] * gamma_j / ( 1.0 / Ny ) + Coeff_Vec_x [ j * ( Nx + 1 ) + i ] * gamma_i / ( 1.0 / Nx );
		}
	}

	return 0;
}

int	coeff_applyLambdax ( double * Coeffs_LambdaxBtr, double * Coeffs_Btr, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	int dim_T = Nx * Ny;
	double gamma_i, lambda_i, h;

	h = 1.0 / Nx;
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			gamma_i = gamma ( i, Nx, xCondition );
			lambda_i = 6.0 * gamma_i * gamma_i / ( ( 6.0 - gamma_i * gamma_i ) * h * h ) ;
			Coeffs_LambdaxBtr [ j * Nx + i ] = Coeffs_Btr [ j * Nx + i ] * lambda_i;
		}
	}

	return 0;
}

int	coeff_applyLambday ( double * Coeffs_LambdayBtr, double * Coeffs_Btr, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	int dim_T = Nx * Ny;
	double gamma_j, lambda_j, h;

	h = 1.0 / Ny;
	for ( int j = 0; j < Ny; j++ )
	{
		gamma_j = gamma ( j, Ny, yCondition );
		lambda_j = 6.0 * gamma_j * gamma_j / ( ( 6.0 - gamma_j * gamma_j ) * h * h ) ;
		for ( int i = 0; i < Nx; i++ )
		{
			Coeffs_LambdayBtr [ j * Nx + i ] = Coeffs_Btr [ j * Nx + i ] * lambda_j;
		}
	}

	return 0;
}


int	coeff_applyLambda ( double * Coeffs_LambdaBtr, double * Coeffs_Btr, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition )
{
	int dim_T = Nx * Ny;
	double gamma_j, lambda_j, h_j;
	double gamma_i, lambda_i, h_i;

	h_j = 1.0 / Ny;
	h_i = 1.0 / Nx;
	for ( int j = 0; j < Ny; j++ )
	{
		gamma_j = gamma ( j, Ny, yCondition );
		lambda_j = 6.0 * gamma_j * gamma_j / ( ( 6.0 - gamma_j * gamma_j ) * h_j * h_j ) ;
		for ( int i = 0; i < Nx; i++ )
		{
			gamma_i = gamma ( i, Nx, xCondition );
			lambda_i = 6.0 * gamma_i * gamma_i / ( ( 6.0 - gamma_i * gamma_i ) * h_i * h_i ) ;

			Coeffs_LambdaBtr [ j * Nx + i ] = Coeffs_Btr [ j * Nx + i ] * ( lambda_j + lambda_i );
		}
	}

	return 0;
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