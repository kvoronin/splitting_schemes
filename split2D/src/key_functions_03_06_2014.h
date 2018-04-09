
extern boundaryConditions;
extern double righthand(int number, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

int F_alpha_beta_fill(double * F, int t, double alpha, double beta, int numSolution, int Mx, int My, int Mz)
{
	//usually we use (alpha, beta) = (1, 0)[explicit], (0, 1)[implicit],(0.5,0.5)[Crank-Nicolson]
	for (int j = 0; j < My; j++ )
	{
		for (int i = 0; i < Mx; i++)
		{
			F[i + j * Mx] = alpha * righthand(numSolution, t, i, j, 1, xpoints, ypoints, zpoints) + beta * righthand(numSolution, t + 1, i, j, 1, xpoints, ypoints, zpoints);
		}
	}
	return 0;
}

int ByMminus1_F_Add(double * Output, double * F, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	double m_ij, m_ijminus1;
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My + 1; j++)
		{
			m_ij = heat_capacity(i,j,1)*density(i,j,1)*hy(j)*hx(i);
			m_ijminus1 = heat_capacity(i,j-1,1)*density(i,j-1,1)*hy(j-1)*hx(i);

			if (j == 0)
			{
				if (yCondition == eNeumann)
					Output[i * (Ny + 1) + j] += 0;
				else
				{
					//temp_0 = hx(i1)*(1.0/m_0)*righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints);
					//f_0 = -temp_0 ; //вычислили B* (M(-1) * Fn+1)

					Output[i * (Ny + 1) + j] += - hx(i) * (1.0/m_ij) * F[j*Nx + i];
				}
			}
			else
				if (j == Ny)
				{
					if (yCondition == eNeumann || yCondition == eMixed)
						Output[i * (Ny + 1) + j] += 0;
					else
					{
						//temp_nminus1 = hx(i1)*(1.0/m_nminus1)*righthand(numSolution,t,i1,Ny-1,k1,xpoints,ypoints,zpoints);
						//f_n = temp_nminus1; //вычислили B* (M(-1) * Fn+1)

						Output[i * (Ny + 1) + j] += hx(i) * (1.0/m_ijminus1) * F[(j-1)*Nx + i];
					}
				}
				else //0 < j < Ny
				{
					//temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints);
					//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints);
					//f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 +  Fn)

					Output[i * (Ny + 1) + j] += hx(i) * ((1.0/m_ijminus1) * F[(j - 1)*Nx + i] - (1.0/m_ij) * F[j*Nx + i]);
				}
		}
	return 0;
}
int BxMminus1_F_Add(double * Output, double * F, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	double m_ij, m_iminus1j;
	for (int j = 0; j < My; j++)
		for (int i = 0; i < Mx + 1; i++)
		{
			m_ij = heat_capacity(i,j,1)*density(i,j,1)*hy(j)*hx(i);
			m_iminus1j = heat_capacity(i-1,j,1)*density(i-1,j,1)*hy(j)*hx(i-1);

			if (i == 0)
			{
				if (xCondition == eNeumann)
					Output[j * (Nx + 1) + i] += 0;
				else
				{
					Output[j * (Nx + 1) + i] += - hy(j) * (1.0/m_ij) * F[j*Nx + i];
				}
			}
			else
				if (i == Nx)
				{
					if (xCondition == eNeumann)
						Output[j * (Nx + 1) + i] += 0;
					else
					{
						Output[j * (Nx + 1) + i] += hy(j) * (1.0/m_iminus1j) * F[j*Nx + i - 1];
					}
				}
				else //0 < i < Nx
				{
					Output[j * (Nx + 1) + i] += hy(j) * ((1.0/m_iminus1j) * F[j*Nx + i - 1] - (1.0/m_ij) * F[j*Nx + i]);
				}
		}
		return 0;
}

int Bxtr_Wx_Add(double * Output, double alpha, double * Wx_n, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	for (int j = 0; j < My; j++)
		for (int i = 0; i < Mx; i++)
			Output[j*Nx + i] += alpha * hy(j) * (Wx_n[j * (Nx + 1) + i + 1] - Wx_n[j * (Nx + 1) + i]);
	return 0;
}

int Bytr_Wy_Add(double * Output, double beta, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My; j++)
			Output[j*Nx + i] += beta * hx(i) * (Wy_n[i * (Ny + 1) + j + 1] - Wy_n[i * (Ny + 1) + j]);
	return 0;
}
int C_Wx_Add(double * Output, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	double cx_ij, cx_iplus1j, koeff_ij, a_ij_old;
	for (int j = 0; j < My; j++)
		for (int i = 0; i < Mx; i++)
		{
			a_ij_old = hx(i)*hy(j)/(heat_conductivity(i,j,1)*6.0);
			koeff_ij = heat_capacity(i,j,1)*density(i,j,1)*a_ij_old;

			cx_ij = koeff_ij * (2 * Vx[j * (Mx + 1) + i] + Vx[j * (Mx + 1) + i + 1]);
			cx_iplus1j = koeff_ij * (Vx[j * (Mx + 1) + i] + 2 * Vx[j * (Mx + 1) + i + 1]);

			Output[j*Nx + i] += alpha * hy(j) * (cx_ij*Wx_n[j * (Mx + 1) + i] + cx_iplus1j*Wx_n[j * (Mx + 1) + i + 1]);
		}
		return 0;
}
int C_Wy_Add(double * Output, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	double cy_ij, cy_ijplus1, koeff_ij, a_ij_old;
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My; j++)
		{
			a_ij_old = hx(i)*hy(j)/(heat_conductivity(i,j,1)*6.0);
			koeff_ij = heat_capacity(i,j,1)*density(i,j,1)*a_ij_old;

			cy_ij = koeff_ij * (2 * Vy[i * (My + 1) + j] + Vy[i * (My + 1) + j + 1]);
			cy_ijplus1 = koeff_ij * (Vy[i * (My + 1) + j] + 2 * Vy[i * (My + 1) + j + 1]);

			Output[j*Nx + i] += beta * hx(i) * (cy_ij * Wy_n[i * (My + 1) + j] + cy_ijplus1*Wy_n[i * (My + 1) + j + 1]);
		}
		return 0;
}


int Btr_W_Add(double * Output, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	Bxtr_Wx_Add(Output, alpha, Wx_n, xCondition, Mx, My, Mz);
	//printf("Bxtr_Wx_Add[0] = %f \n", Output[0]);
	Bytr_Wy_Add(Output, beta, Wy_n, yCondition, Mx, My, Mz);
	//printf("Bxtr_Wx + Bytr_Wy[0] = %f \n", Output[0]);
	return 0;
}
int C_W_Add(double * Output, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	C_Wx_Add(Output, alpha, Vx, Wx_n, xCondition, Mx, My, Mz);
	//printf("Bxtr_Wx_Add[0] = %f \n", Output[0]);
	C_Wy_Add(Output, beta, Vy, Wy_n, yCondition, Mx, My, Mz);
	//printf("Bxtr_Wx + Bytr_Wy[0] = %f \n", Output[0]);
	return 0;
}

int ByMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	Btr_W_Add(Temp, alpha, Wx_n, xCondition, beta, Wy_n, yCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	ByMminus1_F_Add(Output, Temp, yCondition, Mx, My, Mz);

	return 0;
}
int BxMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	Btr_W_Add(Temp, alpha, Wx_n, xCondition, beta, Wy_n, yCondition, Mx, My, Mz);
	BxMminus1_F_Add(Output, Temp, xCondition, Mx, My, Mz);

	return 0;
}

int BxMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	C_W_Add(Temp, alpha, Vx, Wx_n, xCondition, beta, Vy, Wy_n, yCondition, Mx, My, Mz);
	BxMminus1_F_Add(Output, Temp, xCondition, Mx, My, Mz);

	return 0;
}
int ByMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	C_W_Add(Temp, alpha, Vx, Wx_n, xCondition, beta, Vy, Wy_n, yCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	ByMminus1_F_Add(Output, Temp, yCondition, Mx, My, Mz);

	return 0;
}
int HeatFlux_Wx0_Init(double * T_0, double * Wx_0, boundaryConditions xCondition, int numSolution, int Mx, int My, int Mz)
{
	double a_0_old, a_i_old, a_iminus1_old, a_nminus1_old;
	double b_0_old, b_i_old, b_n_old = 0;
	double temper_0, temper_k, temper_n = 0;
	double Tx_0_temp, Tx_1_temp, Ty_0_temp, Ty_1_temp, wx_0_temp, wx_1_temp, wy_0_temp, wy_1_temp;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	double * alfax = (double *) malloc ((Mx + 1) * sizeof(double));
	double * betax = (double *) malloc ((Mx + 1) * sizeof(double));

	switch(xCondition)
	{
	case eNeumann:
		//правильное для неоднородных условий Неймана по x
		//инициализация для потока по x - Wx_n
		//Нейман, x
		for (int i = 0; i < My*Mz; i++)
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			wx_0_temp = wx_0_bound(numSolution,0,j1,0*tau);
			wx_1_temp = wx_1_bound(numSolution,Mx-1,j1,0*tau);

			betax[0] = 0 ;
			betax[1] = wx_0_temp ;
			alfax[0] = 0 ;
			alfax[1] = 0 ;
			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = T_0[k - 1 + j1*Mx + k1*Mx*My] - T_0[k + j1*Mx + k1*Mx*My];

				if (k==Mx-1)
				{
					alfax[k + 1] = 0;
					betax[k + 1] = ( ( temper_k - a_i_old*wx_1_temp ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( ( temper_k - a_iminus1_old*wx_0_temp ) ) / (  b_i_old );
				}
				else
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( ( temper_k ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
				}

				//if (i == 0)
				//{
				//	printf("k = %d \n",k);
				//	printf("temper_%d = %f \n",k, temper_k);
				//}
			}

			Wx_0[i*(Mx+1) + Mx ] = wx_1_temp ;
			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				Wx_0[i*(Mx+1) + Mx - j] = Wx_0[i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 - j];
			}

		}
		break;
	case eDirichlet:
		//для неоднородных условий Дирихле
		//Дирихле, x
		{;}

		for (int i = 0; i < My*Mz; i++)
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			Tx_0_temp = Tx_0_bound(numSolution,0,j1,0*tau);
			Tx_1_temp = Tx_1_bound(numSolution,Mx-1,j1,0*tau);

			//к-ты матрицы А
			a_0_old = hx(0)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;

			temper_0 = Tx_0_temp - T_0[0 + j1*Mx + k1*Mx*My];
			betax[0] = 0 ;
			alfax[0] = 0 ;
			alfax[1] = -0.5 ;
			betax[1] =  temper_0 / b_0_old; 
			for ( int k = 1 ; k < Mx ; k++ )
			{
				i1 = k;
				//к-ты матрицы А
				a_i_old = hx(k)/(heat_conductivity(i1,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)/(heat_conductivity(i1-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;
				//правая часть
				temper_k = T_0[k - 1 + j1*Mx + k1*Mx*My] - T_0[k + j1*Mx + k1*Mx*My];

				alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
				betax[k + 1] = ( ( temper_k ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
			}
			a_nminus1_old = hx(Mx-1)/(heat_conductivity(Mx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			temper_n = T_0[Mx-1 + j1*Mx + k1*Mx*My] - Tx_1_temp;

			Wx_0[i*(Mx+1) + Mx] = ( temper_n  -  betax[Mx]*a_nminus1_old ) / (alfax[Mx]*a_nminus1_old + b_n_old);

			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				Wx_0[i*(Mx+1) + Mx - j] = Wx_0 [i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 -j]; 
			}
		}
		break;
	}

	free(alfax);
	free(betax);

	return 0;
}



int HeatFlux_Wy0_Init(double * T_0, double * Wy_0, boundaryConditions yCondition, int numSolution, int Mx, int My, int Mz)
{
	double a_0_old, a_i_old, a_iminus1_old, a_nminus1_old;
	double b_0_old, b_i_old, b_n_old = 0;
	double temper_0, temper_k, temper_n = 0;
	double Tx_0_temp, Tx_1_temp, Ty_0_temp, Ty_1_temp, wx_0_temp, wx_1_temp, wy_0_temp, wy_1_temp;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	double * alfay = (double *) malloc ((Mx + 1) * sizeof(double));
	double * betay = (double *) malloc ((Mx + 1) * sizeof(double));

#ifdef DEBUGE
	int dim_Wy = Mx * (My + 1);
	FILE * temper_y0 = fopen("temper_y0.xls","wt");
	fprintf(temper_y0,"%s \t %s \t %s \t %s \n","rhand","rhand-exact","A_wexact","rhand-A_wexact");
	double * rhand_Wy0 = (double*)malloc(dim_Wy * sizeof(double));
	if (rhand_Wy0 == NULL)
		printf("Cannot allocate rhand_Wy0 \n");
	double * AWy0_exact = (double*)malloc(dim_Wy * sizeof(double));
	double * Wy0_exact = (double*)malloc(dim_Wy * sizeof(double));
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My + 1; j++)
			Wy0_exact[i*(My + 1) + j] = -exact_gradientY(numSolution, xpoints[i] + 0.5*hx(i), ypoints[j],0,0);
	//FILE * AWy0_exact_file = fopen("AWy0_exact.xls","wt");
#endif

	switch(yCondition)
	{
	case eNeumann:
		//инициализация для потока по y - Wy_n
		//Нейман, y
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			wy_0_temp = wy_0_bound(numSolution,i1,0,0*tau);
			wy_1_temp = wy_1_bound(numSolution,i1,My-1,0*tau);

			betay[0] = 0 ;
			betay[1] = wy_0_temp ;
			alfay[0] = 0 ;
			alfay[1] = 0 ;
			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = hx(i1)*(T_0[i1 + (k - 1)*Mx + k1*Mx*My] - T_0[i1 + k*Mx + k1*Mx*My]);

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = ( ( temper_k - a_i_old*wy_1_temp ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k - a_iminus1_old*wy_0_temp ) ) / (  b_i_old );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
			}
			Wy_0[i*(My+1) + My ] = wy_1_temp ;
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_0[i*(My+1) + My - j ] = Wy_0 [i*(My+1) + My + 1 - j ] * alfay[My + 1 - j] + betay[My + 1 - j]; 
			}
		}
		break;
	case eDirichlet:
		//для неоднородных условий Дирихле
		//Дирихле, y

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			Ty_0_temp = Ty_0_bound(numSolution,i1,0,0*tau);
			Ty_1_temp = Ty_1_bound(numSolution,i1,My-1,0*tau);

			//Ty_0_temp = 0.0;
			//Ty_1_temp = 0.0;


			//к-ты матрицы А
			a_0_old = hy(0)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = Ty_0_temp - T_0[i1 + 0*Mx + k1*Mx*My];

			double diff = 0.0;
			diff = temper_0 - (-hy(0)*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[0], 0, 0));
#ifdef DEBUGE
			rhand_Wy0[i*(My + 1) + 0] = temper_0 / hy(0);
			AWy0_exact[i*(My + 1) + 0] = b_0_old * Wy0_exact[i*(My + 1) + 0] + a_0_old * Wy0_exact[i*(My + 1) + 1];
			double diff2 = 0.0;
			diff2 = temper_0 - AWy0_exact[i*(My + 1) + 0];
			//fprintf(temper_y0,"%f \t %f \n",temper_0, diff);
			fprintf(temper_y0,"%f \t %f \t %f \t %f \n",temper_0, diff, AWy0_exact[i*(My + 1) + 0], temper_0 - AWy0_exact[i*(My + 1) + 0]);
			if (i == 0)
			{
				printf("%f \t %f \t %f \t %f \n",temper_0, diff, AWy0_exact[i*(My + 1) + 0], temper_0 - AWy0_exact[i*(My + 1) + 0]);
				printf("boundT = %f \t T[] = %f \n",Ty_0_temp, T_0[i1 + 0*Mx + k1*Mx*My]);
			}
			//printf("rhand[%d] = %f \n",i*(My + 1) + 0,rhand_Wy0[i*(My + 1) + 0]);
#endif


			betay[0] = 0 ;
			alfay[0] = 0 ;
			alfay[1] = -0.5 ;
			//betay[1] = diff / b_0_old; 
			//betay[1] = diff2 / b_0_old; 
			betay[1] =  temper_0 / b_0_old; 

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;
				temper_k = T_0[i1 + (k - 1)*Mx + k1*Mx*My] - T_0[i1 + k*Mx + k1*Mx*My];

#ifdef DEBUGE
				diff = temper_k - (-hy(0)*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[k], 0, 0));
				rhand_Wy0[i*(My + 1) + k] = temper_k / hy(k);
				AWy0_exact[i*(My + 1) + k] = a_iminus1_old * Wy0_exact[i*(My + 1) + k-1] + b_i_old * Wy0_exact[i*(My + 1) + k] + a_i_old * Wy0_exact[i*(My + 1) + k + 1];
				diff2 = temper_k - AWy0_exact[i*(My + 1) + k];
				fprintf(temper_y0,"%f \t %f \t %f \t %f \n",temper_k, diff, AWy0_exact[i*(My + 1) + k],temper_k - AWy0_exact[i*(My + 1) + k]);
				//printf("rhand[%d] = %f \n",i*(My + 1) + k,rhand_Wy0[i*(My + 1) + k]);
#endif

				alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
				//betay[k + 1] = ( ( diff ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				//betay[k + 1] = ( ( diff2 ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );

			}
			a_nminus1_old = hy(My-1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			temper_n = T_0[i1 + (My-1)*Mx + k1*Mx*My] - Ty_1_temp;

#ifdef DEBUGE
			diff = temper_n - (-hy(0)*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[My], 0, 0));
			rhand_Wy0[i*(My + 1) + My] = temper_n / hy(My-1);
			AWy0_exact[i*(My + 1) + My] = a_nminus1_old * Wy0_exact[i*(My + 1) + My-1] + b_n_old * Wy0_exact[i*(My + 1) + My];
			diff2 = temper_n - AWy0_exact[i*(My + 1) + My];
			fprintf(temper_y0,"%f \t %f \t %f \t %f \n",temper_n, diff, AWy0_exact[i*(My + 1) + My], temper_n - AWy0_exact[i*(My + 1) + My]);
			if (i == 0)
			{
				printf("rhand[%d] = %f \n",i*(My + 1) + My,rhand_Wy0[i*(My + 1) + My]);
				printf("boundTy1 = %f T[] = %f \n",Ty_1_temp, T_0[i1 + (My-1)*Mx + k1*Mx*My]);
			}
#endif

			//Wy_n[i*(My+1) + My] = ( diff  -  betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);
			//Wy_n[i*(My+1) + My] = ( diff2  -  betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);
			Wy_0[i*(My+1) + My] = ( temper_n  -  betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);

			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_0[i*(My+1) + My - j ] = Wy_0 [i*(My+1) + My + 1 - j ] * alfay[My + 1 - j] + betay[My + 1 - j]; 
			}
		}
		break;
	case eMixed:
		//Дирихле-Нейман, y
		{;}
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			Ty_0_temp = Ty_0_bound(numSolution,i1,0,0*tau);
			wy_1_temp = wy_1_bound(numSolution,i1,My-1,0*tau);

			//к-ты матрицы А
			a_0_old = hy(0)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = Ty_0_temp - T_0[i1 + 0*Mx + k1*Mx*My];

			betay[0] = 0 ;
			alfay[0] = 0 ;
			alfay[1] = -0.5 ;
			betay[1] =  temper_0 / b_0_old; 

			//printf("temper_0 = %f \n",temper_0);

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = T_0[i1 + (k - 1)*Mx + k1*Mx*My] - T_0[i1 + k*Mx + k1*Mx*My];

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = ( ( temper_k - a_i_old*wy_1_temp ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}

				//printf("k = %d \n",k);
				//printf("temper_k = %f \n",temper_k);
				//_getch();

			}

			Wy_0[i*(My+1) + My] = wy_1_temp;
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_0[i*(My+1) + My - j ] = Wy_0 [i*(My+1) + My + 1 - j ] * alfay[My + 1 - j] + betay[My + 1 - j]; 
			}
		}
		break;
	}

#ifdef DEBUGE
	fclose(temper_y0);
#endif
	free(alfay);
	free(betay);

	return 0;
}

int HeatFlux_Init(double * T_0, double * Wx_0, boundaryConditions xCondition, double * Wy_0, boundaryConditions yCondition, int numSolution, int Mx, int My, int Mz)
{
	//computes w_0 by solving Aw_0 = BT_0 + g_0

	HeatFlux_Wx0_Init(T_0, Wx_0, xCondition, numSolution, Mx, My, Mz);
	HeatFlux_Wy0_Init(T_0, Wy_0, yCondition, numSolution, Mx, My, Mz);

	{;}

#ifdef DEBUGE
	int dim_Wy = Mx * (My + 1);
//	FILE * temper_y0 = fopen("temper_y0.xls","wt");
//	fprintf(temper_y0,"%s \t %s \t %s \t %s \n","rhand","rhand-exact","A_wexact","rhand-A_wexact");
	double * rhand_Wy0 = (double*)malloc(dim_Wy * sizeof(double));
	if (rhand_Wy0 == NULL)
		printf("Cannot allocate rhand_Wy0 \n");
	double * AWy0_exact = (double*)malloc(dim_Wy * sizeof(double));
	double * Wy0_exact = (double*)malloc(dim_Wy * sizeof(double));
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My + 1; j++)
			Wy0_exact[i*(My + 1) + j] = -exact_gradientY(numSolution, xpoints[i] + 0.5*hx(i), ypoints[j],0,0);
	//FILE * AWy0_exact_file = fopen("AWy0_exact.xls","wt");
#endif


	{;}

#ifdef DEBUGE
	{
		char* ext = ".xls";
		char nx_string[5];
		sprintf(nx_string,"_%d",Mx);
		char filename[30];
		strcpy(filename,"rhand_Wy_0");
		strcat(filename,nx_string);
		strcat(filename,ext);
		FILE * rhand_Wy0_file = fopen(filename,"wt");
		for ( int i = 0 ; i < Mx ; i++ )
		{
			for (int j = 0; j < My + 1; j++)
			{
				fprintf(rhand_Wy0_file,"%f \t %f \n",rhand_Wy0[j + i*(My + 1)],-exact_gradientY(numSolution, xpoints[i]+0.5*hx(i), ypoints[j], 0, 0));
				//if (i == 0)
				//	printf("x = %f y = %f \n%f \t %f \n",xpoints[i]+0.5*hx(i), ypoints[j], rhand_Wy0[j + i*(My + 1)],-exact_gradientY(numSolution, xpoints[i]+0.5*hx(i), ypoints[j], 0, 0));
			}
		}
		fclose(rhand_Wy0_file);

	}
#endif

	return 0;
}

int Temperature_standart_step(double * T_n, double * T_nplus1, double * Vx, double * Wx_n, double * Wx_nplus1, double * Vy, double * Wy_n, double * Wy_nplus1, double sigma, double tau, int numSolution, int t, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double m_i, a_i_old, temp, koeff_k, conv_k, conv_kplus1, CONVECT_x, CONVECT_y;
	for (int i = 0; i < My*Mz; i++)
	{
		j1 = i - (i/My)*My;
		k1 = i/My;
		for ( int k = 0 ; k < Mx ; k++ )
		{
			i1 = k;
			m_i = hx(i1)*hy(j1)*density(i1,j1,k1)*heat_capacity(i1,j1,k1);

			temp = 0.5*(righthand(numSolution,(t+1),i1,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,j1,k1,xpoints,ypoints,zpoints));

			a_i_old = 1.0/(heat_conductivity(i1,j1,k1)*6.0);
			koeff_k = a_i_old;
			conv_k = koeff_k*( 2*Vx[i*(Mx+1)+ i1] + Vx[i*(Mx+1)+ i1+1])*Wx_n[i*(Mx+1) + i1];
			conv_kplus1 = koeff_k*( Vx[i*(Mx+1)+ i1] + 2*Vx[i*(Mx+1)+ i1+1])*Wx_n[i*(Mx+1) + i1 + 1];
			CONVECT_x = conv_k + conv_kplus1;

			//правильное
			//T_nplus1[i*Mx + k] = T[i*Mx + k] + tau * (1.0/m_i) * temp - tau * (1.0/m_i)*hy(j1)*( 0.5*(Wx_n[i*(Mx+1) + k + 1] - Wx_n[i*(Mx+1) + k] ) + 0.5*(Wx_nplus1[i*(Mx+1) + k + 1] - Wx_nplus1[i*(Mx+1) + k] )) ;
			T_nplus1[i*Mx + k] = T_n[i*Mx + k] + tau * (1.0/m_i) * temp + tau * CONVECT_x - tau * (1.0/m_i)*hy(j1)*hz(k1)*( sigma*(Wx_n[i*(Mx+1) + k + 1] - Wx_n[i*(Mx+1) + k] ) + (1 - sigma)*(Wx_nplus1[i*(Mx+1) + k + 1] - Wx_nplus1[i*(Mx+1) + k] )) ;
		}

	}
	//		_getch();

	for ( int i = 0 ; i < Mx*Mz ; i++ )
	{
		k1 = i-(i/Mz)*Mz;
		i1 = i/Mz;
		for ( int j = 0 ; j < My ; j++ )
		{
			j1 = j;
			m_i = hx(i1)*hy(j)*density(i1,j1,k1)*heat_capacity(i1,j1,k1);

			a_i_old = 1.0/(heat_conductivity(i1,j1,k1)*6.0);
			koeff_k = a_i_old;
			conv_k = koeff_k*(2*Vy[i*(My+1)+ j1] + Vy[i*(My+1) + j1+1])*Wy_n[i*(My+1) + j1];
			conv_kplus1 = koeff_k*(Vy[i*(My+1)+ j1] + 2*Vy[i*(My+1) + j1+1])*Wy_n[i*(My+1) + j1 + 1];
			CONVECT_y = conv_k + conv_kplus1;

			//T_nplus1[i1 + j1*Mx + k1*Mx*My] +=  - tau * (1.0/m_i)* hx(i1)*( 0.5*(Wy_n[i*(My + 1) + 1 + j] - Wy_n[i*(My + 1) + j] ) + 0.5*(Wy_nplus1[i*(My + 1) + 1 + j] - Wy_nplus1[i*(My + 1) + j] )) ;
			T_nplus1[i1 + j1*Mx + k1*Mx*My] +=  tau * CONVECT_y - tau * (1.0/m_i)* hx(i1)*hz(k1)*( sigma*(Wy_n[i*(My + 1) + 1 + j] - Wy_n[i*(My + 1) + j] ) + (1 - sigma)*(Wy_nplus1[i*(My + 1) + 1 + j] - Wy_nplus1[i*(My + 1) + j] )) ;
		}
	}
	return 0.0;
}

int Righthand_step1_wynplus05_old(double * righthandY, double * Vx, double * Wx_n, double * Vy, double * Wy_n, boundaryConditions yCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old, g_po_t_0, g_po_t_1, m_0, m_i, m_iminus1, m_nminus1, temp_0, temp_k, temp_kminus1, temp_nminus1, f_0, f_k, f_n, BMBy_0, BMBy_k, BMBy_n, ByMBx_0, ByMBx_k, ByMBx_n, koeff_0, koeff_k, koeff_kminus1, koeff_nminus1, BMCy_0, BMCy_k, BMCy_n, ByMCx_0, ByMCx_k, ByMCx_n, cx_minus1_kminus1, cx_0_kminus1, cx_0_k, cx_plus1_k, cy_minus1_kminus1, cy_0_kminus1, cy_0_k, cy_plus1_k, cx_0_0, cx_plus1_0, cy_0_0, cy_plus1_0, BMBGy_k, ByMBGx_k, BMBGy_0, BMBGy_n, ByMBGx_0, ByMBGx_n, cx_minus1_nminus1, cx_0_nminus1, cy_minus1_nminus1, cy_0_nminus1;

#ifdef DEBUGE
	FILE *deleit0 = fopen("right_wy_nplus05.xls","wt");
	FILE *righthandfile = fopen("righthand_wy_nplus05.xls","wt");
	FILE *bmby_n = fopen("bmby_n_wy_nplus05.xls","wt");
	FILE *bymbx_n = fopen("bymbx_n_wy_nplus05.xls","wt");
	FILE *bymcx_n = fopen("bymcx_n_wy_nplus05.xls","wt");
	FILE *bmcy_n = fopen("bmcy_n_wy_nplus05.xls","wt");
	FILE *f_nfile = fopen("f_n_wy_nplus05.xls","wt");
	FILE *g_po_t = fopen("g_po_t_wy_nplus05.xls","wt");
#ifdef NEW_BOUNDARY_CONDITION
	FILE *deleit4 = fopen("rightGnew_wy_nplus05.xls","wt");
#endif
#endif

	switch(yCondition)
	{
	case eNeumann:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ y-компоненты
		//Нейман, y
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//g_po_t_0 = (wy_0_bound(numSolution,i1,0,(t+0.5)*tau) - wy_0_bound(numSolution,i1,0,t*tau)) / (0.5*tau);
			//g_po_t_1 = (wy_1_bound(numSolution,i1,My-1,(t+0.5)*tau) - wy_1_bound(numSolution,i1,My-1,t*tau)) / (0.5*tau);

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau) + exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
			g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,(t+1)*tau) + exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = (wy_0_bound(numSolution,i1,0,(t+1)*tau) - wy_0_bound(numSolution,i1,0,t*tau)) / (tau);
			g_po_t_1 = (wy_1_bound(numSolution,i1,My-1,(t+1)*tau) - wy_1_bound(numSolution,i1,My-1,t*tau)) / (tau);
#endif

			//g_po_t_0 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1))*(exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),t*tau))/tau;
			//g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*(exact_gradientY(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;

			righthandY[i * (My + 1) + 0] = g_po_t_0;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				temp_k = 0;
				temp_kminus1 = 0;
				temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints);
				temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints);
				//temp_k = hx(i1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 +  Fn)

				BMBy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));
				ByMBx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My]));

				BMBGy_k = 0.0;
				ByMBGx_k = 0.0;
#ifdef NEW_BOUNDARY_CONDITION
				BMBGy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Gy_temp[i*(My+1) + k] - Gy_temp[i*(My+1) + k-1])  - (1.0/m_i)*(Gy_temp[i*(My+1) + k+1] - Gy_temp[i*(My+1) + k]));
				ByMBGx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Gx_temp[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Gx_temp[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Gx_temp[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Gx_temp[i1 + k*(Mx+1) + k1*(Mx+1)*My]));
#endif
				koeff_kminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(i1,k,k1)*density(i1,k,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_k = koeff_k*( 2*Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);
				cx_plus1_k = koeff_k*( Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[i*(My+1) + k-1] + Vy[i*(My+1) + k]);
				cy_0_kminus1 = koeff_kminus1*( Vy[i*(My+1) + k-1] + 2*Vy[i*(My+1) + k]);
				cy_0_k = koeff_k*( 2*Vy[i*(My+1) + k] + Vy[i*(My+1) + k + 1]);
				cy_plus1_k = koeff_k*( Vy[i*(My+1) + k] + 2*Vy[i*(My+1) + k + 1]);

				BMCy_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[i*(My+1) + k-1] + cy_0_kminus1*Wy_n[i*(My+1) + k])  - (1.0/m_i)*(cy_0_k*Wy_n[i*(My+1) + k] + cy_plus1_k*Wy_n[i*(My+1) + k+1]));
				ByMCx_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_kminus1*Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - (1.0/m_i)*(cx_0_k*Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_k*Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]));


				if (k==My-1)
				{
					righthandY[i * (My + 1) + k] =  - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k + BMBGy_k + ByMBGx_k) + f_k - a_i_old * g_po_t_1 ;
				}
				else if (k==1)
				{
					righthandY[i * (My + 1) + k] =  - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k + BMBGy_k + ByMBGx_k) + f_k - a_iminus1_old * g_po_t_0 ;
				}
				else
				{
					righthandY[i * (My + 1) + k] =  - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k + BMBGy_k + ByMBGx_k) + f_k;
				}

				//           if (i==0)
				//           {
				//printf("k = %d \n",k);
				//           	printf("right_%d = %15.15f \n",k,- (BMBy_k + ByMBx_k) + f_k);
				//           	printf("BMBy_%d = %15.15f ByMBx_%d = %15.15f f_%d = %15.15f g_1 = %15.15f \n",k,BMBy_k,k, ByMBx_k,k,f_k,g_po_t_1);
				//             printf("temp_k = %15.15f temp_kminus1 = %15.15f \n",temp_k, temp_kminus1);
				//             printf("temp_k_1 = %15.15f temp_k_2 = %15.15f \n",righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints),righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				//             printf("temp_kminus1_1 = %15.15f temp_kminus1_2 = %15.15f \n",righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints),righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));
				//           	printf("alfay[%d+1] = %15.15f  betay[%d+1] = %15.15f \n",k,alfay[k+1],k,betay[k+1]);
				//             _getch();
				//           }
				//fprintf(deleteit1,"%f \n",- (BMBy_k + ByMBx_k) + f_k);

			}
		}
		break;
	case eDirichlet:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ y-компоненты
		//Дирихле, y

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//g_po_t_0 = hx(i1)*hz(k1)*(exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
			//g_po_t_1 = -hx(i1)*hz(k1)*(exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;

			//g_po_t_0 = hx(i1)*hz(k1)*(Ty_0_bound(numSolution,i1,0,(t+0.5)*tau) - Ty_0_bound(numSolution,i1,0,t*tau))/(0.5*tau);
			//g_po_t_1 = -hx(i1)*hz(k1)*(Ty_1_bound(numSolution,i1,0,(t+0.5)*tau) - Ty_1_bound(numSolution,i1,0,t*tau))/(0.5*tau);

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau) + exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
			g_po_t_1 = -hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,(t+1)*tau) + exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*(Ty_0_bound(numSolution,i1,0,(t+1)*tau) - Ty_0_bound(numSolution,i1,0,t*tau))/(tau);
			g_po_t_1 = -hx(i1)*hz(k1)*(Ty_1_bound(numSolution,i1,0,(t+1)*tau) - Ty_1_bound(numSolution,i1,0,t*tau))/(tau);
#endif

			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1);

			temp_0 = 0;
			temp_0 = hx(i1)*(1.0/m_0)*righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints);
			//temp_0 = hx(i1)*(1.0/m_0)*0.5*(righthand(numSolution,t+1,i1,0,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints));


			f_0 = -temp_0 ; //вычислили B* (M(-1) * Fn+1)
			BMBy_0 = hx(i1)*hx(i1)*( - (1.0/m_0)*(Wy_n[i*(My+1) + 1] - Wy_n[i*(My+1) + 0]));
			ByMBx_0 = hx(i1)*hy(0)*( - (1.0/m_0)*(Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			BMBGy_0 = 0.0;
			ByMBGx_0 = 0.0;
#ifdef NEW_BOUNDARY_CONDITION
			BMBGy_0 = hx(i1)*hx(i1)*( - (1.0/m_0)*(Gy_temp[i*(My+1) + 1] - Gy_temp[i*(My+1) + 0]));
			ByMBGx_0 = hx(i1)*hy(0)*( - (1.0/m_0)*(Gx_temp[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My] - Gx_temp[i1 + 0*(Mx+1) + k1*(Mx+1)*My]));
#endif
			koeff_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*a_0_old;

			cx_0_0 = koeff_0*( 2*Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);
			cx_plus1_0 = koeff_0*( Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);

			cy_0_0 = koeff_0*( 2*Vy[i*(My+1) + 0] + Vy[i*(My+1) + 0 + 1]);
			cy_plus1_0 = koeff_0*( Vy[i*(My+1) + 0] + 2*Vy[i*(My+1) + 0 + 1]);

			BMCy_0 = hx(i1)*hz(k1)*( - (1.0/m_0)*(cy_0_0*Wy_n[i*(My+1) + 0] + cy_plus1_0*Wy_n[i*(My+1) + 0+1]));
			ByMCx_0 = hx(i1)*hz(k1)*(- (1.0/m_0)*(cx_0_0*Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_0*Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			righthandY[i*(My + 1) + 0] = (- (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0 + BMBGy_0 + ByMBGx_0) + f_0 + g_po_t_0);

#ifdef DEBUGE
			if (i==0)
			{
				printf("right_0 = %f \n",- (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0) + f_0 + g_po_t_0);
				printf("BMBy_0 = %f  f_0 = %f \n",BMBy_0,f_0);
				printf("ByMBx_0 = %f g_0 = %f \n",ByMBx_0,g_po_t_0);
				printf("special_Bxtr_Wx_0 = %f \n",hy(0)*(Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My]));
				printf("special_Bytr_Wy_0 = %f \n",hx(i1) * (Wy_n[i*(My+1) + 1] - Wy_n[i*(My+1) + 0]));
				//printf("g_po_t_0_1 = %15.15f, g_po_t_0_2 = %15.15f \n",hx(i1)*hz(k1)*exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau),hx(i1)*hz(k1)*exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
				//printf("g_po_t_0_old = %15.15f\n",hx(i1)*hz(k1)*(Ty_0_bound(numSolution,i1,0,(t+1)*tau) - Ty_0_bound(numSolution,i1,0,t*tau))/(tau));
				//printf("f_0_1 = %15.15f, f_0_2 = %15.15f \n",hx(i1)*(1.0/m_0)*righthand(numSolution,t+1,i1,0,k1,xpoints,ypoints,zpoints), hx(i1)*(1.0/m_0)*righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints));
				//printf("alfay[1] = %15.15f  betay[1] = %15.15f \n",alfay[1],betay[1]);
				fprintf(deleit0,"%15.15f \n", - (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0) + f_0 + g_po_t_0);

				fprintf(righthandfile,"%15.15f \n", - (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0 + BMBGy_0 + ByMBGx_0) + f_0 + g_po_t_0);
				fprintf(bmby_n,"%15.15f \n", BMBy_0);
				fprintf(bymbx_n,"%15.15f \n", ByMBx_0);
				fprintf(bmcy_n,"%15.15f \n", BMCy_0);
				fprintf(bymcx_n,"%15.15f \n", ByMCx_0);
				fprintf(f_nfile,"%15.15f \n", f_0);
				fprintf(g_po_t,"%15.15f \n", g_po_t_0);

#ifdef NEW_BOUNDARY_CONDITION
				fprintf(deleit4,"%f \n", BMBGy_0 + ByMBGx_0);
#endif
				//_getch();
			}
#endif

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				temp_k = 0;
				temp_kminus1 = 0;
				temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints);
				temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints);
				//temp_k = hx(i1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 +  Fn)

				BMBy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));
				ByMBx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My]));

				BMBGy_k = 0.0;
				ByMBGx_k = 0.0;
#ifdef NEW_BOUNDARY_CONDITION
				BMBGy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Gy_temp[i*(My+1) + k] - Gy_temp[i*(My+1) + k-1])  - (1.0/m_i)*(Gy_temp[i*(My+1) + k+1] - Gy_temp[i*(My+1) + k]));
				ByMBGx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Gx_temp[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Gx_temp[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Gx_temp[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Gx_temp[i1 + k*(Mx+1) + k1*(Mx+1)*My]));
#endif


				koeff_kminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(i1,k,k1)*density(i1,k,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_k = koeff_k*( 2*Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);
				cx_plus1_k = koeff_k*( Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[i*(My+1) + k-1] + Vy[i*(My+1) + k]);
				cy_0_kminus1 = koeff_kminus1*( Vy[i*(My+1) + k-1] + 2*Vy[i*(My+1) + k]);
				cy_0_k = koeff_k*( 2*Vy[i*(My+1) + k] + Vy[i*(My+1) + k + 1]);
				cy_plus1_k = koeff_k*( Vy[i*(My+1) + k] + 2*Vy[i*(My+1) + k + 1]);

				BMCy_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[i*(My+1) + k-1] + cy_0_kminus1*Wy_n[i*(My+1) + k])  - (1.0/m_i)*(cy_0_k*Wy_n[i*(My+1) + k] + cy_plus1_k*Wy_n[i*(My+1) + k+1]));
				ByMCx_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_kminus1*Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - (1.0/m_i)*(cx_0_k*Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_k*Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]));

				righthandY[i*(My + 1) + k] =  - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k + BMBGy_k + ByMBGx_k) + f_k;

#ifdef DEBUGE
//				if (i==0)
//				{
//					printf("right_%d = %15.15f \n",k,- (BMBy_k + ByMBx_k) + f_k);
//					printf("BMBy_%d = %15.15f ByMBx_%d = %15.15f f_%d = %15.15f g_1 = %15.15f \n",k,BMBy_k,k, ByMBx_k,k,f_k,g_po_t_1);
//					//printf("alfay[%d+1] = %15.15f  betay[%d+1] = %15.15f \n",k,alfay[k+1],k,betay[k+1]);
//					//fprintf(deleit0,"%15.15f \n", - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k);
//
//					fprintf(righthandfile,"%15.15f \n", - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k + BMBGy_k + ByMBGx_k) + f_k);
//					fprintf(bmby_n,"%15.15f \n", BMBy_k);
//					fprintf(bymbx_n,"%15.15f \n", ByMBx_k);
//					fprintf(bmcy_n,"%15.15f \n", BMCy_k);
//					fprintf(bymcx_n,"%15.15f \n", ByMCx_k);
//					fprintf(f_nfile,"%15.15f \n", f_k);
//					fprintf(g_po_t,"%15.15f \n",0.0);
//
//#ifdef NEW_BOUNDARY_CONDITION
//					fprintf(deleit4,"%f \n", BMBGy_k + ByMBGx_k);
//#endif
//					//_getch();
//				}
#endif
				//							fprintf(deleteit1,"%f \n",- (BMBy_k + ByMBx_k) + f_k);
			}

			//к-ты матрицы А
			a_nminus1_old = hy(My-1)*hx(i1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			//к-ты матрицы M
			m_nminus1 = heat_capacity(i1,My-1,k1)*density(i1,My-1,k1)*hy(My-1)*hx(i1);

			temp_nminus1 = 0;
			temp_nminus1 = hx(i1)*(1.0/m_nminus1)*righthand(numSolution,t,i1,My-1,k1,xpoints,ypoints,zpoints);
			//temp_nminus1 = hx(i1)*(1.0/m_nminus1)*0.5*(righthand(numSolution,t+1,i1,My-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,My-1,k1,xpoints,ypoints,zpoints));

			f_n = temp_nminus1; //вычислили B* (M(-1) * Fn+1)
			BMBy_n = hx(i1)*hx(i1)*((1.0/m_nminus1)*(Wy_n[i*(My+1) + My] - Wy_n[i*(My+1) + My-1]));
			ByMBx_n = hx(i1)*hy(My-1)*((1.0/m_nminus1)*(Wx_n[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]) );

			BMBGy_n = 0.0;
			ByMBGx_n = 0.0;
#ifdef NEW_BOUNDARY_CONDITION
			BMBGy_n = hx(i1)*hx(i1)*((1.0/m_nminus1)*(Gy_temp[i*(My+1) + My] - Gy_temp[i*(My+1) + My-1]));
			ByMBGx_n = hx(i1)*hy(My-1)*((1.0/m_nminus1)*(Gx_temp[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] - Gx_temp[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]) );
#endif
			koeff_nminus1 = heat_capacity(i1,My-1,k1)*density(i1,My-1,k1)*a_nminus1_old;

			cx_minus1_nminus1 = koeff_nminus1*( 2*Vx[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]);
			cx_0_nminus1 = koeff_nminus1*( Vx[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]);

			cy_minus1_nminus1 = koeff_nminus1*( 2*Vy[i*(My+1) + My-1] + Vy[i*(My+1) + My]);
			cy_0_nminus1 = koeff_nminus1*( Vy[i*(My+1) + My-1] + 2*Vy[i*(My+1) + My]);

			BMCy_n = hx(i1)*hz(k1)*((1.0/m_nminus1)*(cy_minus1_nminus1*Wy_n[i*(My+1) + My-1] + cy_0_nminus1*Wy_n[i*(My+1) + My]));
			ByMCx_n = hx(i1)*hz(k1)*((1.0/m_nminus1)*(cx_minus1_nminus1*Wx_n[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_nminus1*Wx_n[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]));

			righthandY[i*(My + 1) + My] =  - (BMBy_n + ByMBx_n - BMCy_n - ByMCx_n + BMBGy_n + ByMBGx_n) + f_n + g_po_t_1;

#ifdef DEBUGE
			//							if (i==0)
			//							{
			//								printf("right_n = %15.15f \n",- (BMBy_n + ByMBx_n ) + f_n + g_po_t_1);
			//								//printf("BMBy_n = %15.15f ByMBx_n = %15.15f f_n = %15.15f g_1 = %15.15f \n",BMBy_n,ByMBx_n,f_n, g_po_t_1);
			//								fprintf(deleit0,"%15.15f \n", - (BMBy_n + ByMBx_n - BMCy_n - ByMCx_n) + f_n + g_po_t_1);
			//
			//								fprintf(righthandfile,"%15.15f \n", - (BMBy_n + ByMBx_n - BMCy_n - ByMCx_n + BMBGy_n + ByMBGx_n) + f_n + g_po_t_1);
			//								fprintf(bmby_n,"%15.15f \n", BMBy_n);
			//								fprintf(bymbx_n,"%15.15f \n", ByMBx_n);
			//								fprintf(bmcy_n,"%15.15f \n", BMCy_n);
			//								fprintf(bymcx_n,"%15.15f \n", ByMCx_n);
			//								fprintf(f_nfile,"%15.15f \n", f_n);
			//								fprintf(g_po_t,"%15.15f \n", g_po_t_1);
			//
			//#ifdef NEW_BOUNDARY_CONDITION
			//								fprintf(deleit4,"%f \n", BMBGy_n + ByMBGx_n);
			//#endif
			//								//_getch();
			//							}
#endif
		}
		break;
	case eMixed:
		//Дирихле-Нейман, y
		{;}

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau) + exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
			g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,(t+1)*tau) + exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*(exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
			g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*(exact_gradientY(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
#endif


			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1);

			temp_0 = 0;
			temp_0 = hx(i1)*(1.0/m_0)*0.5*(righthand(numSolution,t+1,i1,0,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints));

			f_0 = -temp_0 ; //вычислили B* (M(-1) * Fn+1)
			BMBy_0 = hx(i1)*hx(i1)*( - (1.0/m_0)*(Wy_n[i*(My+1) + 1] - Wy_n[i*(My+1) + 0]));
			ByMBx_0 = hx(i1)*hy(0)*( - (1.0/m_0)*(Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			koeff_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*a_0_old;

			cx_0_0 = koeff_0*( 2*Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);
			cx_plus1_0 = koeff_0*( Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);

			cy_0_0 = koeff_0*( 2*Vy[i*(My+1) + 0] + Vy[i*(My+1) + 0 + 1]);
			cy_plus1_0 = koeff_0*( Vy[i*(My+1) + 0] + 2*Vy[i*(My+1) + 0 + 1]);

			BMCy_0 = hx(i1)*hz(k1)*( - (1.0/m_0)*(cy_0_0*Wy_n[i*(My+1) + 0] + cy_plus1_0*Wy_n[i*(My+1) + 0+1]));
			ByMCx_0 = hx(i1)*hz(k1)*(- (1.0/m_0)*(cx_0_0*Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_0*Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			righthandY[i*(My + 1) + 0] =  - (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0) + f_0 + g_po_t_0;

			//fprintf(file1,"%f \t %10.10f\n",ypoints[0],f_0);
			//fprintf(file2,"%f \t %10.10f\n",ypoints[0],-BMBy_0);
			//fprintf(file3,"%f \t %10.10f\n",ypoints[0],-ByMBx_0);
			//fprintf(file4,"%f \t %10.10f \n",ypoints[0] + 0.5*hy(0),(1.0/hy(0))*(Wy_n[i*(My+1) + 1] - Wy_n[i*(My+1) + 0]));
			//fprintf(file5,"%f \t %10.10f \n",ypoints[0] + 0.5*hy(0),(n*PI/2)*(n*PI/2)*exact_solution510(xpoints[i1]+0.5*hx(i1), ypoints[0] + 0.5*hy(0), zpoints[k1] + 0.5*hz(k1), t*tau, m, n));

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				temp_k = 0;
				temp_kminus1 = 0;
				temp_k = hx(i1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				temp_kminus1 = hx(i1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 +  Fn)

				BMBy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));

				//							BMBy_k2 = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));

				ByMBx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My]));

				koeff_kminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(i1,k,k1)*density(i1,k,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_k = koeff_k*( 2*Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);
				cx_plus1_k = koeff_k*( Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[i*(My+1) + k-1] + Vy[i*(My+1) + k]);
				cy_0_kminus1 = koeff_kminus1*( Vy[i*(My+1) + k-1] + 2*Vy[i*(My+1) + k]);
				cy_0_k = koeff_k*( 2*Vy[i*(My+1) + k] + Vy[i*(My+1) + k + 1]);
				cy_plus1_k = koeff_k*( Vy[i*(My+1) + k] + 2*Vy[i*(My+1) + k + 1]);

				BMCy_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[i*(My+1) + k-1] + cy_0_kminus1*Wy_n[i*(My+1) + k])  - (1.0/m_i)*(cy_0_k*Wy_n[i*(My+1) + k] + cy_plus1_k*Wy_n[i*(My+1) + k+1]));
				ByMCx_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_kminus1*Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - (1.0/m_i)*(cx_0_k*Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_k*Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]));


				if (k==My-1)
				{
					righthandY[i*(My + 1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k - a_i_old * g_po_t_1;
				}
				else
				{
					righthandY[i*(My + 1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k ;
				}

				//							righthand_vector1[i*(My+1) + k] = - (BMBy_k + ByMBx_k) + f_k;
				/*
				if (i==0)
				{
				printf("k = %d \n",k);
				printf("alfa[k+1] = %f beta[k+1] = %f \n",alfay[k+1],betay[k+1]);
				printf("right_k = %f \n",- (BMBy_k + ByMBx_k) + f_k);
				printf("f_k = %f \n",f_k);
				//printf("tempkminus1 = %f temp_lminus1 = %f \n",temp_kminus1, temp_k);
				printf("BMBy_k = %f \n",BMBy_k);
				printf("Wy_n[...+ k-1] = %f Wy_n[...+ k] = %f Wy_n[...+ k+1] = %f",Wy_n[i*(My+1) + k-1],Wy_n[i*(My+1) + k], Wy_n[i*(My+1)+ k+1]);
				printf("s m_i razniza1 = %f razniza2 = %f \n",(1.0/m_iminus1)*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1]),(1.0/m_i)*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));
				printf("s h   razniza1 = %f razniza2 = %f \n",(1.0/hy(k-1))*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1]),(1.0/hy(k))*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));
				//temp2 = fabs(-heat_conductivity(i1,j1,k1)*exact_gradient510y(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau, m, n) );
				//printf("proizv1 = %f proizv2 = %f \n",-heat_conductivity(i1,j1,k1)*exact_gradient510y(xpoints[i1]+0.5*hx(i1), ypoints[kj1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau, m, n),-heat_conductivity(i1,j1,k1)*exact_gradient510y(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau, m, n));
				printf("proizvodnay_1 = %f proizvodnay2 = %f \n",(n*PI/2)*(n*PI/2)*exact_solution510(xpoints[i1]+0.5*hx(i1), ypoints[k-1] + 0.5*hy(k-1), zpoints[k1] + 0.5*hz(k1), t*tau, m, n),(n*PI/2)*(n*PI/2)*exact_solution510(xpoints[i1]+0.5*hx(i1), ypoints[k] + 0.5*hy(k), zpoints[k1] + 0.5*hz(k1), t*tau, m, n) );
				double slag1 = hx(i1)*hx(i1)*(1.0/m_iminus1)*(Wy_n[i*(My+1) + k] - Wy_n[i*(My+1) + k-1]);
				double slag2 = hx(i1)*hx(i1)*(1.0/m_i)*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]);
				printf("slag1 = %f slag2 = %f \n",slag1, slag2);
				printf("hy_kminus1 = %f hy_k = %f \n",hy(k-1),hy(k));
				printf("ByMBx_k = %f \n",ByMBx_k);
				double slag3 = hx(i1)*hy(k-1)*(1.0/m_iminus1)*(Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				double slag4 = hx(i1)*hy(k)*(1.0/m_i)*(Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My]);
				//printf("slag3 = %f slag4 = %f \n",slag3, slag4);
				printf("\n");
				//printf("i1 = %d k1 = %d \n",i1,k1);
				_getch();
				}
				*/
				//fprintf(file1,"%f \t %10.10f \n",ypoints[k],f_k);
				//fprintf(file2,"%f \t %10.10f \n",ypoints[k],-BMBy_k);
				//fprintf(file3,"%f \t %10.10f \n",ypoints[k],-ByMBx_k);
				//fprintf(file4,"%f \t %10.10f \n",ypoints[k] + 0.5*hy(k),(1.0/hy(k))*(Wy_n[i*(My+1) + k+1] - Wy_n[i*(My+1) + k]));
				//fprintf(file5,"%f \t %10.10f \n",ypoints[k] + 0.5*hy(k),(n*PI/2)*(n*PI/2)*exact_solution510(xpoints[i1]+0.5*hx(i1), ypoints[k] + 0.5*hy(k), zpoints[k1] + 0.5*hz(k1), t*tau, m, n));

			}
			righthandY[i*(My+1) + My] = g_po_t_1 ;
			//						righthand_vector1[i*(My+1) + My] = 0;

			//fprintf(file1,"%f \t %10.10f\n",ypoints[My],0);
			//fprintf(file2,"%f \t %10.10f\n",ypoints[My],0);
			//fprintf(file3,"%f \t %10.10f\n",ypoints[My],0);

		}

		//fclose(file1);
		//fclose(file2);
		//fclose(file3);
		//fclose(file4);
		//fclose(file5);

		break;
	}
	{;}
#ifdef DEBUGE
	fclose(righthandfile);
	fclose(bmby_n);
	fclose(bymbx_n);
	fclose(bymcx_n);
	fclose(bmcy_n);
	fclose(f_nfile);
	fclose(g_po_t);
	fclose(deleit0);
#ifdef NEW_BOUNDARY_CONDITION
	fclose(deleit4);
#endif
#endif
	return 0.0;
}

int Righthand_step1_wxnplus05_old(double * righthandX, double * Vx, double * Wx_n, double * Vy, double * Wy_n, double* Wy_nplus05,  boundaryConditions xCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old, a_0, a_iminus1, a_i, a_nminus1, b_0, b_i, b_n, g_po_t_0, g_po_t_1, m_0, m_i, m_iminus1, m_nminus1, temp_0, temp_k, temp_kminus1, temp_nminus1, f_0, f_k, f_n, BMBx_0, BMBx_k, BMBx_n, BxMBy_0, BxMBy_k, BxMBy_n, koeff_0, koeff_k, koeff_kminus1, koeff_nminus1, BMCx_0, BMCx_k, BMCx_n, BxMCy_0, BxMCy_k, BxMCy_n, cx_minus1_kminus1, cx_0_kminus1, cx_0_k, cx_plus1_k, cy_minus1_kminus1, cy_0_kminus1, cy_0_k, cy_plus1_k, cx_0_0, cx_plus1_0, cy_0_0, cy_plus1_0, BMBGy_k, ByMBGx_k, BMBGy_0, BMBGy_n, ByMBGx_0, ByMBGx_n, cx_minus1_nminus1, cx_0_nminus1, cy_minus1_nminus1, cy_0_nminus1;

#ifdef DEBUGE
	FILE * deleteit = fopen("BxMBy.xls","wt");
	FILE * deleteit2 = fopen("BxMBy2D.xls","wt");
	FILE * deleteit3 = fopen("right_wx_nplus05.xls","wt");

	FILE * deleit0 = fopen("right_wx_nplus05.xls","wt");
	FILE * righthandfile = fopen("righthand_wx_nplus05.xls","wt");
	FILE * bmbx_n = fopen("bmbx_wx_nplus05.xls","wt");
	FILE * bxmby_n = fopen("bxmby_wx_nplus05.xls","wt");
	FILE * bmcx_n = fopen("bmcx_wx_nplus05.xls","wt");
	FILE * bxmcy_n = fopen("bxmcy_wx_nplus05.xls","wt");
	FILE * f_nfile = fopen("f_wx_nplus05.xls","wt");
	FILE * g_po_t = fopen("gt_wx_nplus05.xls","wt");
#endif

	switch(xCondition)
	{
	case eNeumann:
		//understep 2 of Step 1 - обращение верхнего блока, т.е ~x-компоненты
		//правильное для условий Неймана, x
		{;}

		for (int i = 0; i < My*Mz; i++ )
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			//g_po_t_0 = -heat_conductivity_func(0,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1))*(exact_gradientX(numSolution,0,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientX(numSolution,0,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),t*tau))/tau;
			//g_po_t_1 = -heat_conductivity_func(1,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1))*(exact_gradientX(numSolution,1,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientX(numSolution,1,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;

			//g_po_t_0 = (wx_0_bound(numSolution,0,j1,(t+0.5)*tau) - wx_0_bound(numSolution,0,j1,t*tau))/(0.5*tau);
			//g_po_t_1 = (wx_1_bound(numSolution,0,j1,(t+0.5)*tau) - wx_1_bound(numSolution,0,j1,t*tau))/(0.5*tau);

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = -heat_conductivity_func(0,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientX_po_t(numSolution, 0,ypoints[j1] + 0.5*hy(j1),0,(t+1)*tau) + exact_gradientX_po_t(numSolution, 0, ypoints[j1] + 0.5*hy(j1),0,t*tau));
			g_po_t_1 = -heat_conductivity_func(1,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientX_po_t(numSolution, 1,ypoints[j1] + 0.5*hy(j1),0,(t+1)*tau) + exact_gradientX_po_t(numSolution, 0, ypoints[j1] + 0.5*hy(j1),0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = (wx_0_bound(numSolution,0,j1,(t+1)*tau) - wx_0_bound(numSolution,0,j1,t*tau))/(tau);
			g_po_t_1 = (wx_1_bound(numSolution,0,j1,(t+1)*tau) - wx_1_bound(numSolution,0,j1,t*tau))/(tau);
#endif

			righthandX[i * (Mx + 1) + 0] = g_po_t_0 ;

#ifdef DEBUGE
			//if (i==0)
			{
				fprintf(deleit0,"%15.15f \n", g_po_t_0);
				fprintf(righthandfile,"%15.15f \n", g_po_t_0);
				fprintf(bmbx_n,"%15.15f \n", 0.0);
				fprintf(bxmby_n,"%15.15f \n", 0.0);
				fprintf(bmcx_n,"%15.15f \n", 0.0);
				fprintf(bxmcy_n,"%15.15f \n", 0.0);
				fprintf(f_nfile,"%15.15f \n", 0.0);
			}
#endif

			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(k,j1,k1)*density(k,j1,k1)*hx(k)*hy(j1);
				m_iminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*hx(k-1)*hy(j1);

				//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
				a_i = a_i_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hy(j1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				temp_k = 0;
				temp_kminus1 = 0;

				temp_k = hy(j1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,k,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,k,j1,k1,xpoints,ypoints,zpoints));
				temp_kminus1 = hy(j1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,k-1,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,k-1,j1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 + Fn)

				BMBx_k = hy(j1)*hy(j1)*((1.0/m_iminus1)*(Wx_n[i*(Mx+1) + k] - Wx_n[i*(Mx+1) + k-1])  - (1.0/m_i)*(Wx_n[i*(Mx+1) + k+1] - Wx_n[i*(Mx+1) + k]));
				BxMBy_k = hy(j1)*(hx(k-1)*(1.0/m_iminus1)*(Wy_nplus05[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz] - Wy_nplus05[j1 + k1*(My+1) + (k-1)*(My+1)*Mz]) - hx(k)*(1.0/m_i)*(Wy_nplus05[j1 + 1 + k1*(My+1) + k*(My+1)*Mz] - Wy_nplus05[j1 + k1*(My+1) + k*(My+1)*Mz]));

				koeff_kminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(k,j1,k1)*density(k,j1,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i*(Mx+1)+ k-1] + Vx[i*(Mx+1)+ k]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i*(Mx+1)+ k-1] + 2*Vx[i*(Mx+1)+ k]);
				cx_0_k = koeff_k*( 2*Vx[i*(Mx+1)+ k] + Vx[i*(Mx+1)+ k+1]);
				cx_plus1_k = koeff_k*( Vx[i*(Mx+1)+ k] + 2*Vx[i*(Mx+1)+ k+1]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[j1 + k1*(My+1) + (k-1)*(My+1)*Mz] + Vy[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz]);
				cy_0_kminus1 = koeff_kminus1*( Vy[j1 + k1*(My+1) + (k-1)*(My+1)*Mz] + 2*Vy[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz]);
				cy_0_k = koeff_k*( 2*Vy[j1 + k1*(My+1) + k*(My+1)*Mz] + Vy[j1 + 1 + k1*(My+1) + k*(My+1)*Mz]);
				cy_plus1_k = koeff_k*( Vy[j1 + k1*(My+1) + k*(My+1)*Mz] + 2*Vy[j1 + 1 + k1*(My+1) + k*(My+1)*Mz]);

				BMCx_k = hy(j1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i*(Mx+1) + k-1] + cx_0_kminus1*Wx_n[i*(Mx+1) + k])  - (1.0/m_i)*(cx_0_k*Wx_n[i*(Mx+1) + k] + cx_plus1_k*Wx_n[i*(Mx+1) + k+1]));
				BxMCy_k = hy(j1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[j1 + k1*(My+1) + (k-1)*(My+1)*Mz] + cy_0_kminus1*Wy_n[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz]) - (1.0/m_i)*(cy_0_k*Wy_n[j1 + k1*(My+1) + k*(My+1)*Mz] + cy_plus1_k*Wy_n[j1 + 1 + k1*(My+1) + k*(My+1)*Mz]));

				if (k==Mx-1)
				{
					righthandX[i*(Mx + 1) + k] = - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k - a_i * g_po_t_1;
				}
				else if (k==1)
				{
					righthandX[i*(Mx + 1) + k] = - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k - a_iminus1 * g_po_t_0;
				}
				else
				{
					righthandX[i*(Mx + 1) + k] = - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k;
				}

#ifdef DEBUGE
				////if (i==0)
				//{
				//	if (k==Mx-1)
				//	{
				//		fprintf(deleit0,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k - a_i * g_po_t_1);
				//		fprintf(righthandfile,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k - a_i * g_po_t_1);
				//	}
				//	else if (k==1)
				//	{
				//		fprintf(deleit0,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k - a_iminus1 * g_po_t_0);
				//		fprintf(righthandfile,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k - a_iminus1 * g_po_t_0);
				//	}
				//	else
				//	{
				//		fprintf(deleit0,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k );
				//		fprintf(righthandfile,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k );
				//	}
				//	fprintf(bmbx_n,"%15.15f \n", BMBx_k);
				//	fprintf(bxmby_n,"%15.15f \n", BxMBy_k);
				//	fprintf(bmcx_n,"%15.15f \n", BMCx_k);
				//	fprintf(bxmcy_n,"%15.15f \n", BxMCy_k);
				//	fprintf(f_nfile,"%15.15f \n", f_k);
				//	fprintf(g_po_t,"%15.15f \n",0.0);
				//}
#endif
			}
			righthandX[i*(Mx + 1) + Mx] = g_po_t_1 ;

#ifdef DEBUGE
			//if (i==0)
			{
				fprintf(deleit0,"%15.15f \n", g_po_t_1);
				fprintf(righthandfile,"%15.15f \n", g_po_t_1);
				fprintf(bmbx_n,"%15.15f \n", 0.0);
				fprintf(bxmby_n,"%15.15f \n", 0.0);
				fprintf(bmcx_n,"%15.15f \n", 0.0);
				fprintf(bxmcy_n,"%15.15f \n", 0.0);
				fprintf(f_nfile,"%15.15f \n", 0.0);
				fprintf(g_po_t,"%15.15f \n",g_po_t_1);
			}
#endif

		}
		break;
	case eDirichlet:
		//understep 2 of Step 1 - обращение верхнего блока, т.е ~x-компоненты
		//правильное для условий Дирихле, x
		{;}

		for (int i = 0; i < My*Mz; i++ )
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			//g_po_t_0 = hy(j1)*hz(k1)*(exact_solution(numSolution,0,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,0,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
			//g_po_t_1 = -hy(j1)*hz(k1)*(exact_solution(numSolution,1,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,1,ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;

			//g_po_t_0 = hy(j1)*hz(k1)*(Tx_0_bound(numSolution,0,j1,(t+0.5)*tau) - Tx_0_bound(numSolution,0,j1,t*tau))/(0.5*tau) ;
			//g_po_t_1 = -hy(j1)*hz(k1)*(Tx_1_bound(numSolution,0,j1,(t+0.5)*tau) - Tx_1_bound(numSolution,0,j1,t*tau))/(0.5*tau) ;

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, 0, ypoints[j1] + 0.5*hy(j1),0,(t+1)*tau) + exact_solution_po_t(numSolution, 0, ypoints[j1] + 0.5*hy(j1),0,t*tau));
			g_po_t_1 = -hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, 1, ypoints[j1] + 0.5*hy(j1),0,(t+1)*tau) + exact_solution_po_t(numSolution, 1, ypoints[j1] + 0.5*hy(j1),0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hy(j1)*hz(k1)*(Tx_0_bound(numSolution,0,j1,(t+1)*tau) - Tx_0_bound(numSolution,0,j1,t*tau))/(tau) ;
			g_po_t_1 = -hy(j1)*hz(k1)*(Tx_1_bound(numSolution,0,j1,(t+1)*tau) - Tx_1_bound(numSolution,0,j1,t*tau))/(tau) ;
#endif


			//к-ты матрицы А
			a_0_old = hx(0)*hy(j1)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(0,j1,k1)*density(0,j1,k1)*hx(0)*hy(j1);

			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hy(j1)*hy(j1)*(1.0/m_0); 

			temp_0 = 0;
			temp_0 = hy(j1)*(1.0/m_0)*0.5*(righthand(numSolution,t+1,0,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,0,j1,k1,xpoints,ypoints,zpoints));

			f_0 = (- temp_0 ); //вычислили B* M(-1) * 0.5 *(Fn+1 + Fn)

			BMBx_0 = hy(j1)*hy(j1)*(- (1.0/m_0)*(Wx_n[i*(Mx+1) + 0+1] - Wx_n[i*(Mx+1) + 0]));
			BxMBy_0 = hy(j1)*(- hx(0)*(1.0/m_0)*(Wy_nplus05[j1 + 1 + k1*(My+1) + 0*(My+1)*Mz] - Wy_nplus05[j1 + k1*(My+1) + 0*(My+1)*Mz]));

			koeff_0 = heat_capacity(0,j1,k1)*density(0,j1,k1)*a_0_old;

			cx_0_0 = koeff_0*( 2*Vx[i*(Mx+1)+ 0] + Vx[i*(Mx+1)+ 0+1]);
			cx_plus1_0 = koeff_0*( Vx[i*(Mx+1)+ 0] + 2*Vx[i*(Mx+1)+ 0+1]);

			cy_0_0 = koeff_0*( 2*Vy[j1 + k1*(My+1) + 0*(My+1)*Mz] + Vy[j1 + 1 + k1*(My+1) + 0*(My+1)*Mz]);
			cy_plus1_0 = koeff_0*( Vy[j1 + k1*(My+1) + 0*(My+1)*Mz] + 2*Vy[j1 + 1 + k1*(My+1) + 0*(My+1)*Mz]);

			BMCx_0 = hy(j1)*(- (1.0/m_0)*(cx_0_0*Wx_n[i*(Mx+1) + 0] + cx_plus1_0*Wx_n[i*(Mx+1) + 0+1]));
			BxMCy_0 = hy(j1)*(- (1.0/m_0)*(cy_0_0*Wy_n[j1 + k1*(My+1) + 0*(My+1)*Mz] + cy_plus1_0*Wy_n[j1 + 1 + k1*(My+1) + 0*(My+1)*Mz]));

			righthandX[i * (Mx + 1) + 0] = - (BMBx_0 + BxMBy_0 - BMCx_0 - BxMCy_0) + f_0 + g_po_t_0;

#ifdef DEBUGE
			//if (i==0)
			{
				fprintf(deleteit,"%f \n", BxMBy_0);
				fprintf(deleteit2,"%f \t", BxMBy_0);
				fprintf(deleteit3,"%f \n", - (BMBx_0 + BxMBy_0 - BMCx_0 - BxMCy_0) + f_0 + g_po_t_0 );

				fprintf(deleit0,"%15.15f \n", - (BMBx_0 + BxMBy_0 - BMCx_0 - BxMCy_0) + f_0 + g_po_t_0 );
				fprintf(righthandfile,"%15.15f \n", - (BMBx_0 + BxMBy_0 - BMCx_0 - BxMCy_0) + f_0 );
				fprintf(bmbx_n,"%15.15f \n", BMBx_0);
				fprintf(bxmby_n,"%15.15f \n", BxMBy_0);
				fprintf(bmcx_n,"%15.15f \n", BMCx_0);
				fprintf(bxmcy_n,"%15.15f \n", BxMCy_0);
				fprintf(f_nfile,"%15.15f \n", f_0);
				fprintf(g_po_t,"%15.15f \n",g_po_t_0);
			}
#endif

			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(k,j1,k1)*density(k,j1,k1)*hx(k)*hy(j1);
				m_iminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*hx(k-1)*hy(j1);

				//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
				a_i = a_i_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hy(j1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				temp_k = 0;
				temp_kminus1 = 0;

				temp_k = hy(j1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,k,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,k,j1,k1,xpoints,ypoints,zpoints));
				temp_kminus1 = hy(j1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,k-1,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,k-1,j1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 + Fn)

				BMBx_k = hy(j1)*hy(j1)*((1.0/m_iminus1)*(Wx_n[i*(Mx+1) + k] - Wx_n[i*(Mx+1) + k-1])  - (1.0/m_i)*(Wx_n[i*(Mx+1) + k+1] - Wx_n[i*(Mx+1) + k]));
				BxMBy_k = hy(j1)*(hx(k-1)*(1.0/m_iminus1)*(Wy_nplus05[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz] - Wy_nplus05[j1 + k1*(My+1) + (k-1)*(My+1)*Mz]) - hx(k)*(1.0/m_i)*(Wy_nplus05[j1 + 1 + k1*(My+1) + k*(My+1)*Mz] - Wy_nplus05[j1 + k1*(My+1) + k*(My+1)*Mz]));

				koeff_kminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(k,j1,k1)*density(k,j1,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i*(Mx+1)+ k-1] + Vx[i*(Mx+1)+ k]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i*(Mx+1)+ k-1] + 2*Vx[i*(Mx+1)+ k]);
				cx_0_k = koeff_k*( 2*Vx[i*(Mx+1)+ k] + Vx[i*(Mx+1)+ k+1]);
				cx_plus1_k = koeff_k*( Vx[i*(Mx+1)+ k] + 2*Vx[i*(Mx+1)+ k+1]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[j1 + k1*(My+1) + (k-1)*(My+1)*Mz] + Vy[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz]);
				cy_0_kminus1 = koeff_kminus1*( Vy[j1 + k1*(My+1) + (k-1)*(My+1)*Mz] + 2*Vy[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz]);
				cy_0_k = koeff_k*( 2*Vy[j1 + k1*(My+1) + k*(My+1)*Mz] + Vy[j1 + 1 + k1*(My+1) + k*(My+1)*Mz]);
				cy_plus1_k = koeff_k*( Vy[j1 + k1*(My+1) + k*(My+1)*Mz] + 2*Vy[j1 + 1 + k1*(My+1) + k*(My+1)*Mz]);

				BMCx_k = hy(j1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i*(Mx+1) + k-1] + cx_0_kminus1*Wx_n[i*(Mx+1) + k])  - (1.0/m_i)*(cx_0_k*Wx_n[i*(Mx+1) + k] + cx_plus1_k*Wx_n[i*(Mx+1) + k+1]));
				BxMCy_k = hy(j1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[j1 + k1*(My+1) + (k-1)*(My+1)*Mz] + cy_0_kminus1*Wy_n[j1 + 1 + k1*(My+1) + (k-1)*(My+1)*Mz]) - (1.0/m_i)*(cy_0_k*Wy_n[j1 + k1*(My+1) + k*(My+1)*Mz] + cy_plus1_k*Wy_n[j1 + 1 + k1*(My+1) + k*(My+1)*Mz]));

				righthandX[i * (Mx + 1) + k] =  - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k ;
#ifdef DEBUGE
				//if (i==0)
				{
					fprintf(deleteit,"%f \n", BxMBy_k);
					fprintf(deleteit2,"%f \t", BxMBy_k);
					fprintf(deleteit3,"%f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k );

					fprintf(deleit0,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k );
					fprintf(righthandfile,"%15.15f \n", - (BMBx_k + BxMBy_k - BMCx_k - BxMCy_k) + f_k );
					fprintf(bmbx_n,"%15.15f \n", BMBx_k);
					fprintf(bxmby_n,"%15.15f \n", BxMBy_k);
					fprintf(bmcx_n,"%15.15f \n", BMCx_k);
					fprintf(bxmcy_n,"%15.15f \n", BxMCy_k);
					fprintf(f_nfile,"%15.15f \n", f_k);
					fprintf(g_po_t,"%15.15f \n",0.0);
				}
#endif
			}

			//к-ты матрицы А
			a_nminus1_old = hx(Mx-1)*hy(j1)/(heat_conductivity(Mx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			//к-ты матрицы mx
			m_nminus1 = heat_capacity(Mx-1,j1,k1)*density(Mx-1,j1,k1)*hy(j1)*hx(Mx-1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hy(j1)*hy(j1)*(1.0/m_nminus1);

			temp_nminus1 = 0;
			temp_nminus1 = hy(j1)*(1.0/m_nminus1)*0.5*(righthand(numSolution,t+1,Mx-1,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,Mx-1,j1,k1,xpoints,ypoints,zpoints));

			f_n = temp_nminus1; //вычислили B* (mx(-1) * Fn+1)

			BMBx_n = hy(j1)*hy(j1)*((1.0/m_nminus1)*(Wx_n[i*(Mx+1) + Mx] - Wx_n[i*(Mx+1) + Mx-1]));
			BxMBy_n = hy(j1)*(hx(Mx-1)*(1.0/m_nminus1)*(Wy_nplus05[j1 + 1 + k1*(My+1) + (Mx-1)*(My+1)*Mz] - Wy_nplus05[j1 + k1*(My+1) + (Mx-1)*(My+1)*Mz]));

			koeff_nminus1 = heat_capacity(Mx-1,j1,k1)*density(Mx-1,j1,k1)*a_nminus1_old;

			cx_minus1_nminus1 = koeff_nminus1*( 2*Vx[i*(Mx+1)+ Mx-1] + Vx[i*(Mx+1)+ Mx]);
			cx_0_nminus1 = koeff_nminus1*( Vx[i*(Mx+1)+ Mx-1] + 2*Vx[i*(Mx+1)+ Mx]);

			cy_minus1_nminus1 = koeff_nminus1*( 2*Vy[j1 + k1*(My+1) + (Mx-1)*(My+1)*Mz] + Vy[j1 + 1 + k1*(My+1) + (Mx-1)*(My+1)*Mz]);
			cy_0_nminus1 = koeff_nminus1*( Vy[j1 + k1*(My+1) + (Mx-1)*(My+1)*Mz] + 2*Vy[j1 + 1 + k1*(My+1) + (Mx-1)*(My+1)*Mz]);

			BMCx_n = hy(j1)*((1.0/m_nminus1)*(cx_minus1_nminus1*Wx_n[i*(Mx+1) + Mx-1] + cx_0_nminus1*Wx_n[i*(Mx+1) + Mx]));
			BxMCy_n = hy(j1)*((1.0/m_nminus1)*(cy_minus1_nminus1*Wy_n[j1 + k1*(My+1) + (Mx-1)*(My+1)*Mz] + cy_0_nminus1*Wy_n[j1 + 1 + k1*(My+1) + (Mx-1)*(My+1)*Mz]));

			righthandX[i * (Mx + 1) + Mx] =  - (BMBx_n + BxMBy_n - BMCx_n - BxMCy_n) + f_n + g_po_t_1;

#ifdef DEBUGE
			//if (i==0)
			{
				fprintf(deleteit,"%f \n", BxMBy_n);
				fprintf(deleteit2,"%f \n", BxMBy_n);
				fprintf(deleteit3,"%f \n", - (BMBx_n + BxMBy_n - BMCx_n - BxMCy_n) + f_n + g_po_t_1 );

				fprintf(deleit0,"%15.15f \n", - (BMBx_n+ BxMBy_n - BMCx_n - BxMCy_n) + f_n + g_po_t_1 );
				fprintf(righthandfile,"%15.15f \n", - (BMBx_n + BxMBy_n - BMCx_n - BxMCy_n) + g_po_t_1 );
				fprintf(bmbx_n,"%15.15f \n", BMBx_n);
				fprintf(bxmby_n,"%15.15f \n", BxMBy_n);
				fprintf(bmcx_n,"%15.15f \n", BMCx_n);
				fprintf(bxmcy_n,"%15.15f \n", BxMCy_n);
				fprintf(f_nfile,"%15.15f \n", f_n);
				fprintf(g_po_t,"%15.15f \n",g_po_t_1);
			}
#endif
		}
		break;
	}

#ifdef DEBUGE
	fclose(deleteit);
	fclose(deleteit2);
	fclose(deleteit3);

	fclose(deleit0);
	fclose(righthandfile);
	fclose(bmbx_n);
	fclose(bxmby_n);
	fclose(bmcx_n);
	fclose(bxmcy_n);
	fclose(f_nfile);
	fclose(g_po_t);
#endif

	return 0;
}

int Righthand_step2_wynplus1_old(double * righthandY, double * Vx, double * Wx_n, double * Wx_nplus1, double * Vy, double * Wy_n, double * Wy_nplus05, boundaryConditions yCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old, a_0, a_iminus1, a_i, a_nminus1, b_0, b_i, b_n, g_po_t_0, g_po_t_1, m_0, m_i, m_iminus1, m_nminus1, temp_0, temp_k, temp_kminus1, temp_nminus1, f_0, f_k, f_n, BMBy_0, BMBy_k, BMBy_n, ByMBx_0, ByMBx_k, ByMBx_n, koeff_0, koeff_k, koeff_kminus1, koeff_nminus1, BMCy_0, BMCy_k, BMCy_n, ByMCx_0, ByMCx_k, ByMCx_n, cx_minus1_kminus1, cx_0_kminus1, cx_0_k, cx_plus1_k, cy_minus1_kminus1, cy_0_kminus1, cy_0_k, cy_plus1_k, cx_0_0, cx_plus1_0, cy_0_0, cy_plus1_0, BMBGy_k, ByMBGx_k, BMBGy_0, BMBGy_n, ByMBGx_0, ByMBGx_n, cx_minus1_nminus1, cx_0_nminus1, cy_minus1_nminus1, cy_0_nminus1;

	switch(yCondition)
	{
	case eNeumann:
		//understep 2 of Step 2 - обращение среднего блока, т.е ~ y-компоненты
		//Нейман, y
		{;}
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//g_po_t_0 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1))*(exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),t*tau))/tau;
			//g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*(exact_gradientY(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau) + exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
			g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,(t+1)*tau) + exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = (wy_0_bound(numSolution,i1,0,(t+1)*tau) - wy_0_bound(numSolution,i1,0,t*tau)) / (tau);
			g_po_t_1 = (wy_1_bound(numSolution,i1,My-1,(t+1)*tau) - wy_1_bound(numSolution,i1,My-1,t*tau)) / (tau);
#endif

			righthandY[i*(My+1) + 0] = g_po_t_0;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hx(i1)*( 1.0/m_i + 1.0/m_iminus1 );

				temp_k = 0;
				temp_kminus1 = 0;
				temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints);
				temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints);
				//temp_k = hx(i1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B * M(-1) * 0.5* (Fn+1 + Fn)

				BMBy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_nplus05[i*(My+1) + k] - Wy_nplus05[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_nplus05[i*(My+1) + k+1] - Wy_nplus05[i*(My+1) + k]));
				ByMBx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Wx_nplus1[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Wx_nplus1[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + k*(Mx+1) + k1*(Mx+1)*My]));

				koeff_kminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(i1,k,k1)*density(i1,k,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_k = koeff_k*( 2*Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);
				cx_plus1_k = koeff_k*( Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[i*(My+1) + k-1] + Vy[i*(My+1) + k]);
				cy_0_kminus1 = koeff_kminus1*( Vy[i*(My+1) + k-1] + 2*Vy[i*(My+1) + k]);
				cy_0_k = koeff_k*( 2*Vy[i*(My+1) + k] + Vy[i*(My+1) + k + 1]);
				cy_plus1_k = koeff_k*( Vy[i*(My+1) + k] + 2*Vy[i*(My+1) + k + 1]);

				BMCy_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[i*(My+1) + k-1] + cy_0_kminus1*Wy_n[i*(My+1) + k])  - (1.0/m_i)*(cy_0_k*Wy_n[i*(My+1) + k] + cy_plus1_k*Wy_n[i*(My+1) + k+1]));
				ByMCx_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_kminus1*Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - (1.0/m_i)*(cx_0_k*Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_k*Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]));

				if (k==My-1)
				{
					righthandY[i*(My+1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k - a_i * g_po_t_1;
				}
				else if (k==1)
				{
					righthandY[i*(My+1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k - a_iminus1 * g_po_t_0;
				}
				else
				{
					righthandY[i*(My+1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k ;
				}

			}
			righthandY[i*(My+1) +  My] = g_po_t_1 ;
		}
		break;
	case eDirichlet:
		//understep 2 of Step 2 - обращение среднего блока, т.е ~ y-компоненты
		//Дирихле, y
		{;}

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//g_po_t_0 = hx(i1)*hz(k1)*(exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
			//g_po_t_1 = -hx(i1)*hz(k1)*(exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau) + exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
			g_po_t_1 = -hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,(t+1)*tau) + exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*(Ty_0_bound(numSolution,i1,0,(t+1)*tau) - Ty_0_bound(numSolution,i1,0,t*tau))/(tau);
			g_po_t_1 = -hx(i1)*hz(k1)*(Ty_1_bound(numSolution,i1,0,(t+1)*tau) - Ty_1_bound(numSolution,i1,0,t*tau))/(tau);
#endif


			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0); 

			temp_0 = 0;
			temp_0 = hx(i1)*(1.0/m_0)*righthand(numSolution,t+1,i1,0,k1,xpoints,ypoints,zpoints);
			//temp_0 = hx(i1)*(1.0/m_0)*0.5*(righthand(numSolution,t+1,i1,0,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints));

			f_0 = -temp_0 ; //вычислили B* (M(-1) * Fn+1)

			BMBy_0 = hx(i1)*hx(i1)*( - (1.0/m_0)*(Wy_nplus05[i*(My+1) + 1] - Wy_nplus05[i*(My+1) + 0]));
			ByMBx_0 = hx(i1)*hy(0)*( - (1.0/m_0)*(Wx_nplus1[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			koeff_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*a_0_old;

			cx_0_0 = koeff_0*( 2*Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);
			cx_plus1_0 = koeff_0*( Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);

			cy_0_0 = koeff_0*( 2*Vy[i*(My+1) + 0] + Vy[i*(My+1) + 0 + 1]);
			cy_plus1_0 = koeff_0*( Vy[i*(My+1) + 0] + 2*Vy[i*(My+1) + 0 + 1]);

			BMCy_0 = hx(i1)*hz(k1)*( - (1.0/m_0)*(cy_0_0*Wy_n[i*(My+1) + 0] + cy_plus1_0*Wy_n[i*(My+1) + 0+1]));
			ByMCx_0 = hx(i1)*hz(k1)*(- (1.0/m_0)*(cx_0_0*Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_0*Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			righthandY[i*(My+1) +  0] = - (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0) + f_0 + g_po_t_0;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				//к-ты матрицы A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hx(i1)*( 1.0/m_i + 1.0/m_iminus1 );

				temp_k = 0;
				temp_kminus1 = 0;
				temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints);
				temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints);
				//temp_k = hx(i1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B * M(-1) * 0.5* (Fn+1 + Fn)

				BMBy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_nplus05[i*(My+1) + k] - Wy_nplus05[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_nplus05[i*(My+1) + k+1] - Wy_nplus05[i*(My+1) + k]));
				ByMBx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Wx_nplus1[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Wx_nplus1[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + k*(Mx+1) + k1*(Mx+1)*My]));

				koeff_kminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(i1,k,k1)*density(i1,k,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_k = koeff_k*( 2*Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);
				cx_plus1_k = koeff_k*( Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[i*(My+1) + k-1] + Vy[i*(My+1) + k]);
				cy_0_kminus1 = koeff_kminus1*( Vy[i*(My+1) + k-1] + 2*Vy[i*(My+1) + k]);
				cy_0_k = koeff_k*( 2*Vy[i*(My+1) + k] + Vy[i*(My+1) + k + 1]);
				cy_plus1_k = koeff_k*( Vy[i*(My+1) + k] + 2*Vy[i*(My+1) + k + 1]);

				BMCy_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[i*(My+1) + k-1] + cy_0_kminus1*Wy_n[i*(My+1) + k])  - (1.0/m_i)*(cy_0_k*Wy_n[i*(My+1) + k] + cy_plus1_k*Wy_n[i*(My+1) + k+1]));
				ByMCx_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_kminus1*Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - (1.0/m_i)*(cx_0_k*Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_k*Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]));

				righthandY[i*(My+1) +  k] =  - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k ;
			}

			//к-ты матрицы А
			a_nminus1_old = hy(My-1)*hx(i1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			//к-ты матрицы M
			m_nminus1 = heat_capacity(i1,My-1,k1)*density(i1,My-1,k1)*hy(My-1)*hx(i1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hx(i1)*hx(i1)*(1.0/m_nminus1);

			temp_nminus1 = 0;
			temp_nminus1 = hx(i1)*(1.0/m_nminus1)*righthand(numSolution,t+1,i1,My-1,k1,xpoints,ypoints,zpoints);
			//temp_nminus1 = hx(i1)*(1.0/m_nminus1)*0.5*(righthand(numSolution,t+1,i1,My-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,My-1,k1,xpoints,ypoints,zpoints));

			f_n = temp_nminus1; //вычислили B* (M(-1) * Fn+1)
			BMBy_n = hx(i1)*hx(i1)*((1.0/m_nminus1)*(Wy_nplus05[i*(My+1) + My] - Wy_nplus05[i*(My+1) + My-1]));
			ByMBx_n = hx(i1)*hy(My-1)*((1.0/m_nminus1)*(Wx_nplus1[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]) );

			koeff_nminus1 = heat_capacity(i1,My-1,k1)*density(i1,My-1,k1)*a_nminus1_old;

			cx_minus1_nminus1 = koeff_nminus1*( 2*Vx[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]);
			cx_0_nminus1 = koeff_nminus1*( Vx[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]);

			cy_minus1_nminus1 = koeff_nminus1*( 2*Vy[i*(My+1) + My-1] + Vy[i*(My+1) + My]);
			cy_0_nminus1 = koeff_nminus1*( Vy[i*(My+1) + My-1] + 2*Vy[i*(My+1) + My]);

			BMCy_n = hx(i1)*hz(k1)*((1.0/m_nminus1)*(cy_minus1_nminus1*Wy_n[i*(My+1) + My-1] + cy_0_nminus1*Wy_n[i*(My+1) + My]));
			ByMCx_n = hx(i1)*hz(k1)*((1.0/m_nminus1)*(cx_minus1_nminus1*Wx_n[i1 + (My-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_nminus1*Wx_n[i1 + 1 + (My-1)*(Mx+1) + k1*(Mx+1)*My]));

			//right
			righthandY[i*(My+1) +  My] = - (BMBy_n + ByMBx_n - BMCy_n - ByMCx_n) + f_n + g_po_t_1;

		}
		break;
	case eMixed:
		//Дирихле-Нейман, y
		{;}

		//					double* righthand_vector3 = (double*)malloc(dim_Wy*sizeof(double));
		//					for (int i = 0; i < dim_Wy; i++)
		//						righthand_vector3[i] = 0;


		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

#ifdef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*0.5*(exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,(t+1)*tau) + exact_solution_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),0,0,t*tau));
			g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*0.5*(exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,(t+1)*tau) + exact_gradientY_po_t(numSolution, xpoints[i1] + 0.5*hx(i1),1,0,t*tau));
#endif
#ifndef EXACT_BOUNDARY_CONDITION
			g_po_t_0 = hx(i1)*hz(k1)*(exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_solution(numSolution,xpoints[i1] + 0.5*hx(i1),0,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
			g_po_t_1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1))*(exact_gradientY(numSolution,xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),(t+1)*tau) - exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1),1,zpoints[k1] + 0.5*hz(k1),t*tau))/tau ;
#endif


			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0); 

			temp_0 = 0;
			temp_0 = hx(i1)*(1.0/m_0)*0.5*(righthand(numSolution,t+1,i1,0,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints));

			f_0 = -temp_0 ; //вычислили B* (M(-1) * Fn+1)
			BMBy_0 = hx(i1)*hx(i1)*( - (1.0/m_0)*(Wy_nplus05[i*(My+1) + 1] - Wy_nplus05[i*(My+1) + 0]));
			ByMBx_0 = hx(i1)*hy(0)*( - (1.0/m_0)*(Wx_nplus1[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			koeff_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*a_0_old;

			cx_0_0 = koeff_0*( 2*Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);
			cx_plus1_0 = koeff_0*( Vx[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]);

			cy_0_0 = koeff_0*( 2*Vy[i*(My+1) + 0] + Vy[i*(My+1) + 0 + 1]);
			cy_plus1_0 = koeff_0*( Vy[i*(My+1) + 0] + 2*Vy[i*(My+1) + 0 + 1]);

			BMCy_0 = hx(i1)*hz(k1)*( - (1.0/m_0)*(cy_0_0*Wy_n[i*(My+1) + 0] + cy_plus1_0*Wy_n[i*(My+1) + 0+1]));
			ByMCx_0 = hx(i1)*hz(k1)*(- (1.0/m_0)*(cx_0_0*Wx_n[i1 + 0*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_0*Wx_n[i1 + 1 + 0*(Mx+1) + k1*(Mx+1)*My]));

			righthandY[i * (My + 1) + 0] =  - (BMBy_0 + ByMBx_0 - BMCy_0 - ByMCx_0) + f_0 + g_po_t_0;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				//к-ты матрицы A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hx(i1)*( 1.0/m_i + 1.0/m_iminus1 );

				temp_k = 0;
				temp_kminus1 = 0;
				//temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints);
				//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints);
				temp_k = hx(i1)*(1.0/m_i)*0.5*(righthand(numSolution,t+1,i1,k,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints));
				temp_kminus1 = hx(i1)*(1.0/m_iminus1)*0.5*(righthand(numSolution,t+1,i1,k-1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints));

				f_k = ( temp_kminus1 - temp_k ); //вычислили B * M(-1) * 0.5* (Fn+1 + Fn)

				BMBy_k = hx(i1)*hx(i1)*((1.0/m_iminus1)*(Wy_nplus05[i*(My+1) + k] - Wy_nplus05[i*(My+1) + k-1])  - (1.0/m_i)*(Wy_nplus05[i*(My+1) + k+1] - Wy_nplus05[i*(My+1) + k]));
				ByMBx_k = hx(i1)*(hy(k-1)*(1.0/m_iminus1)*(Wx_nplus1[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - hy(k)*(1.0/m_i)*(Wx_nplus1[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My] - Wx_nplus1[i1 + k*(Mx+1) + k1*(Mx+1)*My]));

				koeff_kminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*a_iminus1_old;
				koeff_k = heat_capacity(i1,k,k1)*density(i1,k,k1)*a_i_old;

				cx_minus1_kminus1 = koeff_kminus1*( 2*Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_kminus1 = koeff_kminus1*( Vx[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My]);
				cx_0_k = koeff_k*( 2*Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);
				cx_plus1_k = koeff_k*( Vx[i1 + k*(Mx+1) + k1*(Mx+1)*My] + 2*Vx[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]);

				cy_minus1_kminus1 = koeff_kminus1*( 2*Vy[i*(My+1) + k-1] + Vy[i*(My+1) + k]);
				cy_0_kminus1 = koeff_kminus1*( Vy[i*(My+1) + k-1] + 2*Vy[i*(My+1) + k]);
				cy_0_k = koeff_k*( 2*Vy[i*(My+1) + k] + Vy[i*(My+1) + k + 1]);
				cy_plus1_k = koeff_k*( Vy[i*(My+1) + k] + 2*Vy[i*(My+1) + k + 1]);

				BMCy_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cy_minus1_kminus1*Wy_n[i*(My+1) + k-1] + cy_0_kminus1*Wy_n[i*(My+1) + k])  - (1.0/m_i)*(cy_0_k*Wy_n[i*(My+1) + k] + cy_plus1_k*Wy_n[i*(My+1) + k+1]));
				ByMCx_k = hx(i1)*hz(k1)*((1.0/m_iminus1)*(cx_minus1_kminus1*Wx_n[i1 + (k-1)*(Mx+1) + k1*(Mx+1)*My] + cx_0_kminus1*Wx_n[i1 + 1 + (k-1)*(Mx+1) + k1*(Mx+1)*My])  - (1.0/m_i)*(cx_0_k*Wx_n[i1 + k*(Mx+1) + k1*(Mx+1)*My] + cx_plus1_k*Wx_n[i1 + 1 + k*(Mx+1) + k1*(Mx+1)*My]));

				if (k==My-1)
				{
					righthandY[i * (My + 1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k - a_i*g_po_t_1;
				}
				else
				{
					righthandY[i * (My + 1) + k] = - (BMBy_k + ByMBx_k - BMCy_k - ByMCx_k) + f_k;
				}

			}
			righthandY[i * (My + 1) + My] = g_po_t_1 ;

		}
		break;
	}


	return 0;
}
int Explicit_Wy(double * Wy_nplus05, double * Wy_n, double * righthandY, boundaryConditions yCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;

	double * alfay = (double*) malloc ((My + 1) * sizeof(double));
	double * betay = (double*) malloc ((My + 1) * sizeof(double));

	switch(yCondition)
	{
	case eNeumann:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ y-компоненты
		//Нейман, y
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			betay[0] = 0 ;
			betay[1] = righthandY[i * (My + 1) + 0];
			alfay[0] = 0 ;
			alfay[1] = 0 ;
			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					//right
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1_old ) / (alfay[k]*a_iminus1_old + b_i_old);
				}
				else if (k==1)
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( righthandY[i * (My + 1) + k] ) / ( b_i_old );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1_old ) / (alfay[k]*a_iminus1_old + b_i_old);
				}

			}

			Wy_nplus05[i*(My+1) + My] =  righthandY[i * (My + 1) + My] ;
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus05[i*(My+1) + My - j] = Wy_nplus05[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
			}
		}
		break;
	case eDirichlet:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ y-компоненты
		//Дирихле, y

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;

			betay[0] = 0 ;
			betay[1] = righthandY[i * (My + 1) + 0] / b_0_old ;
			alfay[0] = 0 ;
			alfay[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
				betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1_old ) / (alfay[k]*a_iminus1_old + b_i_old);
			}

			//к-ты матрицы А
			a_nminus1_old = hy(My-1)*hx(i1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			Wy_nplus05[i*(My+1) +  My] = (   righthandY[i * (My + 1) + My]  - betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus05[i*(My+1) + My - j] = Wy_nplus05[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
			}
		}
		break;
	case eMixed:
		//Дирихле-Нейман, y
		{;}

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;

			betay[0] = 0 ;
			betay[1] =  righthandY[i * (My + 1) + 0]/b_0_old ;
			alfay[0] = 0 ;
			alfay[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1_old ) / (alfay[k]*a_iminus1_old + b_i_old);
				}
				else
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1_old ) / (alfay[k]*a_iminus1_old + b_i_old);
				}
			}
			Wy_nplus05[i*(My+1) + My] = righthandY[i*(My + 1) + My] ;

			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus05[i*(My+1) + My - j] = Wy_nplus05[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
			}
		}
		break;
	}
	{;}

	free(alfay);
	free(betay);

	int dim_Wy = Mx * (My + 1);
	for (int i = 0; i < dim_Wy; i++)
	{
		Wy_nplus05[i] = Wy_n[i] + sigma*tau*Wy_nplus05[i];
#ifdef ZEROING
		Wy_nplus05[i] = 0;
#endif
	}

	return 0;
}

int Implicit_Wx(double * Wx_nplus05, double * Wx_n, double * righthandX, boundaryConditions xCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old, a_0, a_iminus1, a_i, a_nminus1, b_0, b_i, b_n, m_0, m_iminus1, m_i, m_nminus1;

	double * alfax = (double*) malloc ((Mx + 1) * sizeof(double));
	double * betax = (double*) malloc ((Mx + 1) * sizeof(double));

	switch(xCondition)
	{
	case eNeumann:
		//understep 2 of Step 1 - обращение верхнего блока, т.е ~x-компоненты
		//правильное для условий Неймана, x
		{;}

		for (int i = 0; i < My*Mz; i++ )
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			betax[0] = 0 ;
			betax[1] = righthandX[i * (Mx + 1) + 0] ;
			alfax[0] = 0 ;
			alfax[1] = 0 ;

			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(k,j1,k1)*density(k,j1,k1)*hx(k)*hy(j1);
				m_iminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*hx(k-1)*hy(j1);

				//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
				a_i = a_i_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hy(j1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==Mx-1)
				{
					alfax[k + 1] = 0;
					betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1 ) / (alfax[k]*a_iminus1 + b_i);
				}
				else if (k==1)
				{
					alfax[k + 1] = (-1)*a_i/(alfax[k]*a_iminus1 + b_i);
					betax[k + 1] = ( righthandX[i * (Mx + 1) + k] ) / ( b_i );
				}
				else
				{
					alfax[k + 1] = (-1)*a_i/(alfax[k]*a_iminus1 + b_i);
					betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1 ) / (alfax[k]*a_iminus1 + b_i);
				}

			}

			Wx_nplus05[i*(Mx+1) + Mx] = righthandX[i * (Mx + 1) + Mx];

			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				Wx_nplus05[i*(Mx+1) + Mx - j] = Wx_nplus05[i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 - j];
			}
		}
		break;
	case eDirichlet:
		//understep 2 of Step 1 - обращение верхнего блока, т.е ~x-компоненты
		//правильное для условий Дирихле, x
		{;}

		for (int i = 0; i < My*Mz; i++ )
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			//к-ты матрицы А
			a_0_old = hx(0)*hy(j1)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(0,j1,k1)*density(0,j1,k1)*hx(0)*hy(j1);

			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hy(j1)*hy(j1)*(1.0/m_0); 

			betax[0] = 0 ;
			betax[1] = righthandX[i * (Mx + 1) + 0] / b_0 ;
			alfax[0] = 0 ;
			alfax[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(k,j1,k1)*density(k,j1,k1)*hx(k)*hy(j1);
				m_iminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*hx(k-1)*hy(j1);

				//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
				a_i = a_i_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hy(j1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				alfax[k + 1] = (-1)*a_i/(alfax[k]*a_iminus1 + b_i);
				betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1 ) / (alfax[k]*a_iminus1 + b_i);

			}

			//к-ты матрицы А
			a_nminus1_old = hx(Mx-1)*hy(j1)/(heat_conductivity(Mx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			//к-ты матрицы mx
			m_nminus1 = heat_capacity(Mx-1,j1,k1)*density(Mx-1,j1,k1)*hy(j1)*hx(Mx-1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hy(j1)*hy(j1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hy(j1)*hy(j1)*(1.0/m_nminus1);

			Wx_nplus05[i*(Mx+1) +  Mx] = (  righthandX[i * (Mx + 1) + Mx]  - betax[Mx]*a_nminus1 ) / (alfax[Mx]*a_nminus1 + b_n);
			
			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				Wx_nplus05[i*(Mx+1) + Mx - j] = Wx_nplus05[i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 - j];
			}

		}
		break;
	}

	{;}

	free(alfax);
	free(betax);

	int dim_Wx = (Mx + 1) * My;
	for (int i = 0; i < dim_Wx; i++)
		Wx_nplus05[i] = Wx_n[i] + sigma*tau*Wx_nplus05[i];


	return 0;
}

int Special_Wx_nplus1(double * Wx_nplus1, double * Wx_nplus05, double * Wx_n, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My;
	for (int i = 0; i < dim_Wx; i++)
		Wx_nplus1[i] = 2.0 * Wx_nplus05[i] - Wx_n[i];
	return 0;
}

int Implicit_Wy(double * Wy_nplus1, double * Wy_nplus05, double * righthandY, boundaryConditions yCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old, a_0, a_iminus1, a_i, a_nminus1, b_0, b_i, b_n, m_0, m_iminus1, m_i, m_nminus1;

	double * alfay = (double*) malloc ((My + 1) * sizeof(double));
	double * betay = (double*) malloc ((My + 1) * sizeof(double));

	switch(yCondition)
	{
	case eNeumann:
		//understep 2 of Step 2 - обращение среднего блока, т.е ~ y-компоненты
		//Нейман, y
		{;}
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			betay[0] = 0 ;
			betay[1] = righthandY[i * (My + 1) + 0];
			alfay[0] = 0 ;
			alfay[1] = 0 ;
			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hx(i1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					//right
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
					//								betay[k + 1] = (  ( righthand_my(k) )  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}
				else if (k==1)
				{
					alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
					//right
					betay[k + 1] = ( righthandY[i * (My + 1) + k] ) / ( b_i );
					//								betay[k + 1] = ( ( righthand_my(k) ) ) / ( b_i );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
					//right	
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
					//								betay[k + 1] = (  ( righthand_my(k) )  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}

			}
			Wy_nplus1[i*(My+1) + My] = righthandY[i * (My + 1) + My] ;
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus1[i*(My+1) + My - j] = Wy_nplus1[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
			}

		}
		break;
	case eDirichlet:
		//understep 2 of Step 2 - обращение среднего блока, т.е ~ y-компоненты
		//Дирихле, y
		{;}

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;


			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0); 

			//
			betay[0] = 0 ;
			//right
			betay[1] = righthandY[i * (My + 1) + 0] / b_0 ;
			//						betay[1] = ( righthand_my(0) )/b_0 ;
			alfay[0] = 0 ;
			alfay[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				//к-ты матрицы A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hx(i1)*( 1.0/m_i + 1.0/m_iminus1 );

				alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
				//right	
				betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				//							betay[k + 1] = (  ( righthand_my(k) )  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
			}

			//к-ты матрицы А
			a_nminus1_old = hy(My-1)*hx(i1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			//к-ты матрицы M
			m_nminus1 = heat_capacity(i1,My-1,k1)*density(i1,My-1,k1)*hy(My-1)*hx(i1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hx(i1)*hx(i1)*(1.0/m_nminus1);

			//right
			Wy_nplus1[i*(My+1) +  My] = (  righthandY[i * (My + 1) + My]  - betay[My]*a_nminus1 ) / (alfay[My]*a_nminus1 + b_n);

			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus1[i*(My+1) + My - j] = Wy_nplus1[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
			}

		}
		break;
	case eMixed:
		//Дирихле-Нейман, y
		{;}

		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hx(i1)*(1.0/m_0); 

			//
			betay[0] = 0 ;
			//right
			betay[1] = righthandY[i * (My + 1) + 0] / b_0 ;
			alfay[0] = 0 ;
			alfay[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1);

				//к-ты матрицы A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hx(i1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hx(i1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}
				else
				{
					alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}

				//							righthand_vector3[i*(My+1) + k] = ( - (BMBy_k + ByMBx_k) + f_k);
			}

			Wy_nplus1[i*(My+1) + My] = righthandY[i * (My + 1) + My] ;
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus1[i*(My+1) + My - j] = Wy_nplus1[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
			}

		}
		break;
	}

	free(alfay);
	free(betay);

	int dim_Wy = Mx * (My + 1);
	for (int i = 0; i < dim_Wy; i++)
	{
		Wy_nplus1[i] = Wy_nplus05[i] + sigma*tau*Wy_nplus1[i];

#ifdef ZEROING
		Wy_nplus1[i] = 0;
#endif

	}

	return 0;
}

int Accuracy_calculate(double * T_n, double * Wx_n, double * Wy_n, int numSolution, int t, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	////////модуль проверки = расчета промежуточной погрешности на очередном шаге
	double eps_max_w = *eps_max_w_pt;
	double eps_max_T = *eps_max_T_pt;
	double eps_l2_w = *eps_l2_w_pt;
	double eps_l2_T = *eps_l2_T_pt;
	double eps_relative_max_w = *eps_relative_max_w_pt;
	double eps_relative_max_T = *eps_relative_max_T_pt;
	double eps_relative_l2_w = *eps_relative_l2_w_pt;
	double eps_relative_l2_T = *eps_relative_l2_T_pt;

	double norm_l2_T = 0;
	double norm_max_T = 0;
	double norm_l2_w = 0;
	double norm_max_w = 0;
	double eps_special_max_w = 0;
	double eps_special_max_T = 0;
	double eps_special_l2_w = 0;
	double eps_special_l2_T = 0;

	double eps_n1_T = 0 ;		//погрешность в равномерной норме на очередном слое
	double eps_n2_T = 0 ;		//погрешность в l2-норме на очередном слое
	double diffz = 0;
	double diffx = 0;
	double diffy = 0;
	double tx = 0;
	double ty = 0;
	double diff_l2x = 0;
	double diff_l2y = 0;
	double diff_l2z = 0;
	double norm_l2_wx = 0;
	double norm_l2_wy = 0;
	double norm_l2_wz = 0;
	double eps_n1_w = 0;
	double eps_wn1 = 0;
	double eps_wn2 = 0;
	double temp = 0, temp2 = 0;


	int i_max_T = 0;
	int j_max_T = 0;
	int k_max_T = 0;
	int i1, j1, k1;
	int dim_T = Mx * My;
	int dim_Wx = (Mx + 1) * My;
	int dim_Wy = Mx * (My + 1);

	for (int ind = 0; ind < My*Mz; ind++)
	{
		j1 = ind - (ind/My)*My;
		k1 = ind/My;
		for ( int i = 0 ; i < Mx ; i++ )
		{
			temp = 0;
			temp2 = 0;

			temp = fabs(T_n[ind*Mx + i] - exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau));
			temp2 = fabs(exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau));

			//if (ind == 0 && t == -1)
			//	printf("temp = %f \n", temp);
			
			if ( fabs(temp) > eps_n1_T )
			{
				eps_n1_T = fabs(temp);
				i_max_T = i;
				j_max_T = j1;
				k_max_T = k1;
			}
			if ( fabs(temp2) > norm_max_T )
				norm_max_T = fabs(temp2);

			eps_n2_T += temp*temp; 
			norm_l2_T += temp2*temp2;
		}
	}
	eps_n2_T = sqrt(eps_n2_T/(dim_T));
	norm_l2_T = sqrt(norm_l2_T/dim_T);

	eps_special_max_T = eps_n1_T / norm_max_T;
	eps_special_l2_T = eps_n2_T / norm_l2_T;

	//модуль проверки погрешности для теплового потока на очередном шаге
	///проверка теплового потока
	{;}

	diff_l2z = 0.0;
	norm_l2_wz = 0.0;

	/////////////////
	int imax_x = 0;
	int jmax_x = 0;
	int kmax_x = 0;
	temp2 = 0;
	for (int i = 0; i < My*Mz; i++)
		for (int k = 0; k < Mx+1; k++)
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			diffx = fabs(Wx_n[i*(Mx+1) + k] - (-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) ));
			temp2 = fabs(-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );

			if (diffx > tx)
			{
				tx = diffx;
				imax_x = k;
				jmax_x = j1;
				kmax_x = k1;
			}
			if (temp2 > norm_max_w)
				norm_max_w = temp2;

			diff_l2x += diffx*diffx; 
			norm_l2_wx += temp2*temp2;
		}
		diff_l2x = diff_l2x/dim_Wx;
		norm_l2_wx = norm_l2_wx/dim_Wx;
		//						printf("diff_l2x = %f \n", diff_l2x);

		if ( tx > eps_n1_w)
			eps_n1_w = tx;

		/////////////////////
		int imax_y = 0;
		int jmax_y = 0;
		int kmax_y = 0;
		temp2 = 0;
		for (int i = 0; i < Mz*Mx; i++)
			for (int k = 0; k < My+1; k++)
			{
				i1 = (i-i%Mz)/Mz;
				j1 = k;
				k1 = i%Mz;

				diffy = fabs(Wy_n[i*(My+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau) ));
				temp2 = fabs(-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );

				if (diffy > ty)
				{
					ty = diffy;
					imax_y = i1;
					jmax_y = j1;
					kmax_y = k1;
				}
				if (temp2 > norm_max_w)
					norm_max_w = temp2;

				diff_l2y += diffy*diffy; 
				norm_l2_wy += temp2*temp2;
			}
			diff_l2y = diff_l2y/dim_Wy;
			norm_l2_wy = norm_l2_wy/dim_Wy;

			if ( ty > eps_n1_w)
				eps_n1_w = ty;

			eps_wn2 = sqrt(diff_l2x + diff_l2y + diff_l2z);
			norm_l2_w = sqrt(norm_l2_wx + norm_l2_wy + norm_l2_wz);

			eps_special_max_w = eps_n1_w / norm_max_w;
			eps_special_l2_w = eps_wn2 / norm_l2_w;

			{;}

			if ( eps_n1_w > eps_max_w)
				eps_max_w = eps_n1_w;
			if ( eps_wn2 > eps_l2_w )
				eps_l2_w = eps_wn2;
			if ( eps_n1_T > eps_max_T)
				eps_max_T = eps_n1_T;
			if ( eps_n2_T > eps_l2_T )
				eps_l2_T = eps_n2_T;

			if ( eps_special_max_w > eps_relative_max_w )
				eps_relative_max_w = eps_special_max_w;
			if ( eps_special_max_T > eps_relative_max_T )
				eps_relative_max_T = eps_special_max_T;

			if ( eps_special_l2_w > eps_relative_l2_w )
				eps_relative_l2_w = eps_special_l2_w;
			if ( eps_special_l2_T > eps_relative_l2_T )
				eps_relative_l2_T = eps_special_l2_T;

			//if ( t == N - 1) 
			//	my_eps_max_T = eps_n1_T;
			//if ( t == N - 1) 
			//	my_eps_l2_T = eps_n2_T;
			//if ( t == N - 1) 
			//	my_eps_max_w = eps_n1_w;
			//if ( t == N - 1) 
			//	my_eps_l2_w = eps_wn2;

			//if ( t == N - 1) 
			//	my_releps_max_T = eps_relative_max_T;
			//if ( t == N - 1) 
			//	my_releps_max_w = eps_relative_max_w;
			//if ( t == N - 1) 
			//	my_releps_l2_T = eps_relative_l2_T;
			//if ( t == N - 1) 
			//	my_releps_l2_w = eps_relative_l2_w;

			//проверка на расходимость
			//////////////////////////////////////////////////
			//if ((eps_max_T > 10000000)||(eps_max_w > 10000000))
			//{
			//	//									printf("slishkom mnoho \n");
			//	//									_getch();
			//	printf("eps_max = %f \n", eps_max_T);
			//	printf("eps_max_w = %f \n", eps_max_w);
			//	fprintf(f1,"too much...\n",t);
			//	goto Exit;
			//}
			////////////////////////конец модуля проверки
			{;}


			//////////////блок выдачи и печати в файл
			if (t == 0)
			{
				printf("t=%d inside function Accuracy \n",t);
				printf("eps_max_T = %e \n",eps_max_T); 
				printf("eps_l2_T = %e \n",eps_l2_T); 
				printf("eps_max_w = %e \n",eps_max_w); 
				printf("eps_l2_w = %e \n",eps_l2_w); 

				printf("eps_wx_l2 = %f \n", sqrt(diff_l2x / dim_Wx));
				printf("eps_wy_l2 = %f \n", sqrt(diff_l2y / dim_Wy));

			}
			if (t+1-((t+1)/print_step)*print_step==0)
			{
				printf("t=%d \n",t);
				printf("eps_max_T = %e \n",eps_max_T); 
				printf("eps_l2_T = %e \n",eps_l2_T); 
				printf("eps_max_w = %e \n",eps_max_w); 
				printf("eps_l2_w = %e \n",eps_l2_w); 

				printf("n-th step \n");
				printf("tx = %f, ty = %f\n",tx, ty); 
				printf("eps_n1_w = %e \n", eps_n1_w);
				printf("eps_n2_w = %e \n", eps_wn2);
				printf("eps_n1_T = %e \n", eps_n1_T);
				printf("eps_n2_T = %e \n", eps_n2_T);

				//if ( eps_n1_w > eps_max_w)
				//	eps_max_w = eps_n1_w;
				//printf("eps_max_w = %e \n",eps_max_w); 
#ifdef GETCH						
				switch(pause)
				{
				case 0:
					break;
				case 1:
					_getch();
					break;
				}
#endif
			}
			*eps_max_w_pt = eps_max_w;
			*eps_max_T_pt = eps_max_T;
			*eps_l2_w_pt = eps_l2_w;
			*eps_l2_T_pt = eps_l2_T;

			*eps_relative_max_w_pt = eps_relative_max_w;
			*eps_relative_max_T_pt = eps_relative_max_T;
			*eps_relative_l2_w_pt = eps_relative_l2_w;
			*eps_relative_l2_T_pt = eps_relative_l2_T;

			return 0;
}



int Temperature_Flux_Update(double * Wx_n, double * Wx_nplus1, double * Wy_n, double * Wy_nplus1, double * T_n, double * T_nplus1, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My;
	int dim_Wy = Mx * (My + 1);
	int dim_T = Mx * My;

	for (int i = 0; i < dim_Wx; i++)
		Wx_n[i] = Wx_nplus1[i];
	for (int i = 0; i < dim_Wy; i++)
		Wy_n[i] = Wy_nplus1[i];
	for (int i = 0; i < dim_T; i++)
		T_n[i] = T_nplus1[i];

	return 0;
}

int Righthand_step1_wynplus05(double * righthandY, double * F, double alpha, double beta, double * Vx, double * Wx_n, double * Vy, double * Wy_n, boundaryConditions yCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wy = Mx * (My + 1);
	int dim_T = Mx * My;
	for (int i = 0; i < dim_Wy; i++)
	{
		righthandY[i] = 0;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	ByMminus1_F_Add(righthandY, F, yCondition, Mx, My, Mz);
	ByMBtr_W_Add(righthandY, F, -1, Wx_n, xCondition, -1, Wy_n, yCondition, Mx, My, Mz);
	ByMC_W_Add(righthandY, F, -1, Vx, Wx_n, xCondition, -1, Vy, Wy_n, yCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wxnplus05(double * righthandX, double * F, double alpha, double beta, double * Vx, double * Wx_n, double * Vy, double * Wy_n, double * Wy_nplus05, boundaryConditions xCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My;
	int dim_T = Mx * My;
	for (int i = 0; i < dim_Wx; i++)
	{
		righthandX[i] = 0;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	BxMminus1_F_Add(righthandX, F, xCondition, Mx, My, Mz);
	//printf("righthandX[0] = %f \n", righthandX[0]);
	BxMBtr_W_Add(righthandX, F, -1, Wx_n, xCondition, -1, Wy_nplus05, yCondition, Mx, My, Mz);
	BxMC_W_Add(righthandX, F, -1, Vx, Wx_n, xCondition, -1, Vy, Wy_n, yCondition, Mx, My, Mz);
	//printf("righthandX[0] = %f \n", righthandX[0]);
	//F is spoiled because it is used as Temp during call of BxMBtr_W_Add
	return 0;

}

int Righthand_step2_wynplus1(double * righthandY, double * F, double alpha, double beta, double * Vx, double * Wx_n, double * Vy, double * Wy_n, boundaryConditions yCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wy = Mx * (My + 1);
	int dim_T = Mx * My;
	for (int i = 0; i < dim_Wy; i++)
	{
		righthandY[i] = 0;
	}
	F_alpha_beta_fill(F, t+1, alpha, beta, numSolution, Mx, My, Mz);
	ByMminus1_F_Add(righthandY, F, yCondition, Mx, My, Mz);
	ByMBtr_W_Add(righthandY, F, -1, Wx_n, xCondition, -1, Wy_n, yCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}
int Splitting2D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus05, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus05, double * Wy_nplus1, double sigma, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	int dim_Wy = Mx * (My + 1);
	int dim_Wx = (Mx + 1) * My;
	int dim_T = Mx * My;
	//computing new FLUX W_nplus1
	//Righthand_step1_wynplus05_old(righthandY, Vx, Wx_n, Vy, Wy_n, yCondition, numSolution, t, Mx, My, Mz);
	Righthand_step1_wynplus05(righthandY, F, 1, 0, Vx, Wx_n, Vy, Wy_n, yCondition, numSolution, t, Mx, My, Mz);
	Explicit_Wy(Wy_nplus05, Wy_n, righthandY, yCondition, tau, 0.5, Mx, My, Mz);

	//Righthand_step1_wxnplus05_old(righthandX, Vx, Wx_n, Vy, Wy_n, Wy_nplus05, xCondition, numSolution, t, Mx, My, Mz);
	Righthand_step1_wxnplus05(righthandX, F, 0.5, 0.5, Vx, Wx_n, Vy, Wy_n, Wy_nplus05, xCondition, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus05, Wx_n, righthandX, xCondition, tau, 0.5, Mx, My, Mz);

	//Righthand_step1_wxnplus05();
	//Explicit_Wx;
	Special_Wx_nplus1(Wx_nplus1, Wx_nplus05, Wx_n, Mx, My, Mz);

	//Righthand_step2_wynplus1_old(righthandY, Vx, Wx_n, Wx_nplus1, Vy, Wy_n, Wy_nplus05, yCondition, numSolution, t, Mx, My, Mz);
	Righthand_step2_wynplus1(righthandY, F, 1, 0, Vx, Wx_nplus1, Vy, Wy_nplus05, yCondition, numSolution, t, Mx, My, Mz);
	Implicit_Wy(Wy_nplus1, Wy_nplus05, righthandY, yCondition, tau, 0.5, Mx, My, Mz);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, Wx_n, Wx_nplus1, Vy, Wy_n, Wy_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, T_n, T_nplus1, Mx, My, Mz);

	Accuracy_calculate(T_n, Wx_n, Wy_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);

	return 0.0;
}

int Uzawa2D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus05, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus05, double * Wy_nplus1, double sigma, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	int dim_Wy = Mx * (My + 1);
	int dim_Wx = (Mx + 1) * My;
	int dim_T = Mx * My;
	//computing new FLUX W_nplus1
	Righthand_step1_wynplus05(righthandY, F, 0.5, 0.5, Vx, Wx_n, Vy, Wy_n, yCondition, numSolution, t, Mx, My, Mz);
	//for (int i = 0; i < Ny + 1; i++)
	//	printf("RhandY[%d] = %f \n",i, righthandY[i]);

	Implicit_Wy(Wy_nplus05, Wy_n, righthandY, yCondition, tau, 1, Mx, My, Mz);
	//for (int i = 0; i < Ny + 1; i++)
	//	printf("Wy_nplus05[%d] = %f \n",i, Wy_nplus05[i]);

	for (int i = 0; i < dim_Wy; i++)
		Temp_Wy[i] = 0.5 * (Wy_nplus05[i] + Wy_n[i]);
	Righthand_step1_wxnplus05(righthandX, F, 0.5, 0.5, Vx, Wx_n, Vy, Wy_n, Temp_Wy, xCondition, numSolution, t, Mx, My, Mz);
	//for (int i = 0; i < Ny + 1; i++)
	//	printf("RhandX[%d] = %f \n",i, righthandX[i]);

	Implicit_Wx(Wx_nplus05, Wx_n, righthandX, xCondition, tau, 1, Mx, My, Mz);
	//for (int i = 0; i < Ny + 1; i++)
	//	printf("Wx_nplus05[%d] = %f \n",i, Wx_nplus05[i]);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, Wx_n, Wx_nplus05, Vy, Wy_n, Wy_nplus05, 0.5, tau, numSolution, t, Mx, My, Mz);
	//for (int i = 0; i < Ny; i++)
	//	printf("T[%d] = %f \n",i, T_nplus1[i]);

	//correction of FLUX Wy
	HeatFlux_Wy0_Init(T_nplus1, Wy_nplus1, yCondition, numSolution, Mx, My, Mz);
	//for (int i = 0; i < Ny + 1; i++)
	//	printf("Wy_nplus1[%d] = %f \n",i, Wy_nplus1[i]);

	Temperature_Flux_Update(Wx_n, Wx_nplus05, Wy_n, Wy_nplus1, T_n, T_nplus1, Mx, My, Mz);

	Accuracy_calculate(T_n, Wx_n, Wy_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);

	return 0.0;
}

int Local1D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus05, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus05, double * Wy_nplus1, double sigma, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	int dim_Wy = Mx * (My + 1);
	int dim_Wx = (Mx + 1) * My;
	int dim_T = Mx * My;
	//computing new FLUX W_nplus1
	Righthand_step1_wynplus05(righthandY, F, 0.5, 0.5, Vx, Wx_n, Vy, Wy_n, yCondition, numSolution, t, Mx, My, Mz);
	Implicit_Wy(Wy_nplus1, Wy_n, righthandY, yCondition, tau, 1.0, Mx, My, Mz);

	for (int i = 0; i < dim_Wy; i++)
		Temp_Wy[i] = 0.5 * (Wy_nplus1[i] + Wy_n[i]);
	Righthand_step1_wxnplus05(righthandX, F, 0.5, 0.5, Vx, Wx_n, Vy, Wy_n, Wy_nplus1, xCondition, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus1, Wx_n, righthandX, xCondition, tau, 1.0, Mx, My, Mz);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, Wx_n, Wx_nplus1, Vy, Wy_n, Wy_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, T_n, T_nplus1, Mx, My, Mz);

	Accuracy_calculate(T_n, Wx_n, Wy_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);

	return 0.0;
}

int TwoCyclic_Local1D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus05, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus05, double * Wy_nplus1, double sigma, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	int dim_Wy = Mx * (My + 1);
	int dim_Wx = (Mx + 1) * My;
	int dim_T = Mx * My;
	//computing new FLUX W_nplus1
	Righthand_step1_wynplus05(righthandY, F, 0, 0, Vx, Wx_n, Vy, Wy_n, yCondition, numSolution, t, Mx, My, Mz);
	Implicit_Wy(Wy_nplus05, Wy_n, righthandY, yCondition, tau, 0.5, Mx, My, Mz);

	Righthand_step1_wxnplus05(righthandX, F, 0, 0, Vx, Wx_n, Vy, Wy_n, Wy_nplus05, xCondition, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus05, Wx_n, righthandX, xCondition, tau, 0.5, Mx, My, Mz);


	//...write it down first
	for (int i = 0; i < dim_Wy; i++)
		Temp_Wy[i] = 0.5 * (Wy_nplus1[i] + Wy_n[i]);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, Wx_n, Wx_nplus1, Vy, Wy_n, Wy_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, T_n, T_nplus1, Mx, My, Mz);

	Accuracy_calculate(T_n, Wx_n, Wy_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);

	return 0.0;
}
