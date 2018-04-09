extern boundaryConditions;


int HeatFlux_Init(double * T_0, double * Wx_0, boundaryConditions xCondition, double * Wy_0, boundaryConditions yCondition, double * Wz_0, boundaryConditions zCondition, int numSolution, int Mx, int My, int Mz);
int HeatFlux_Wx0_Init(double * T_0, double * Wx_0, boundaryConditions xCondition, int numSolution, int Mx, int My, int Mz);
int HeatFlux_Wy0_Init(double * T_0, double * Wy_0, boundaryConditions yCondition, int numSolution, int Mx, int My, int Mz);
int HeatFlux_Wz0_Init(double * T_0, double * Wz_0, boundaryConditions zCondition, int numSolution, int Mx, int My, int Mz);
int Velocity_Init(double *Vx, double *Vy, double *Vz, int Mx, int My, int Mz);

int Accuracy_calculate(double * T_n, double * Wx_n, double * Wy_n, double * Wz_n, int numSolution, int t, int Mx, int My, int Mz, int print_step,
					   double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt);
int	Accuracy_calculateL2Pr(double * T_n, double * Wx_n, double * Wy_n, double * Wz_n,
						   double * T_nAv, double * Wx_Av, double * Wy_Av, double * Wz_Av,
						   boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition,
						   int numSolution, int t, int Mx, int My, int Mz, int print_step, 
						   double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt,
						   double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt);

int Righthand_step1_wynplus13(double * righthandY, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz);
int Righthand_step1_wznplus13(double * righthandZ, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz);
int Righthand_step1_wxnplus1(double * righthandX, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, int numSolution, int t, int Mx, int My, int Mz);
//not needed, included into Righthand_step1_wxnplus1
//int Righthand_step1_wxnplus1_local(double * righthandX, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz);
int Righthand_step1_wynplus1(double * righthandY, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, int numSolution, int t, int Mx, int My, int Mz);
//not needed
//int Righthand_step1_wynplus1_local(double * righthandY, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double * Wx_nplus1, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz);
int Righthand_step1_wznplus1(double * righthandZ, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, int numSolution, int t, int Mx, int My, int Mz);
//not needed
//int Righthand_step1_wznplus1_local(double * righthandZ, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double * Wx_nplus1, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double *Wy_nplus1, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz);
int F_alpha_beta_fill(double * F, int t, double alpha, double beta, int numSolution, int Mx, int My, int Mz);
int BzMminus1_F_Add(double * Output, double * F, boundaryConditions zCondition, int Mx, int My, int Mz);
int ByMminus1_F_Add(double * Output, double * F, boundaryConditions yCondition, int Mx, int My, int Mz);
int BxMminus1_F_Add(double * Output, double * F, boundaryConditions xCondition, int Mx, int My, int Mz);

int ByMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int BzMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int BxMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);

int Btr_W_Add(double * Output, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int Bxtr_Wx_Add(double * Output, double alpha, double * Wx_n, boundaryConditions xCondition, int Mx, int My, int Mz);
int Bytr_Wy_Add(double * Output, double beta, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz);
int Bztr_Wz_Add(double * Output, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int ByMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int BzMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int BxMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int C_W_Add(double * Output, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);
int C_Wx_Add(double * Output, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, int Mx, int My, int Mz);
int C_Wy_Add(double * Output, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz);
int C_Wz_Add(double * Output, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz);

int Explicit_Wy(double * Wy_nplus05, double * Wy_n, double * righthandY, boundaryConditions yCondition, double tau, double sigma, int Mx, int My, int Mz);
int Explicit_Wz(double * Wz_nplus05, double * Wz_n, double * righthandZ, boundaryConditions zCondition, double tau, double sigma, int Mx, int My, int Mz);
int Explicit_Wx(double * Wx_nplus05, double * Wx_n, double * righthandX, boundaryConditions xCondition, double tau, double sigma, int Mx, int My, int Mz);
int Implicit_Wx(double * Wx_nplus05, double * Wx_n, double * righthandX, boundaryConditions xCondition, double tau, double sigma, int Mx, int My, int Mz);
int Implicit_Wy(double * Wy_nplus05, double * Wy_n, double * righthandY, boundaryConditions yCondition, double tau, double sigma, int Mx, int My, int Mz);
int Implicit_Wz(double * Wz_nplus05, double * Wz_n, double * righthandZ, boundaryConditions zCondition, double tau, double sigma, int Mx, int My, int Mz);

int Temperature_standart_step(double * T_n, double * T_nplus1, double * Vx, double coeff_Vx, double * Wx_n, double * Wx_nplus1, double * Vy, double coeff_Vy, double * Wy_n, double * Wy_nplus1, double * Vz, double coeff_Vz, double * Wz_n, double * Wz_nplus1, double sigma, double tau, int numSolution, int t, int Mx, int My, int Mz);
int Temperature_Flux_Update(double * Wx_n, double * Wx_nplus1, double * Wy_n, double * Wy_nplus1, double * Wz_n, double * Wz_nplus1, double * T_n, double * T_nplus1, int Mx, int My, int Mz);

int DouglasGunn_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, 
								double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus1,
								double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus13, double * Wy_nplus1, 
								double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus13, double * Wz_nplus1, 
								int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step,
								double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt,
								double * eps_relative_max_w_pt, double * eps_relative_max_T_pt,	double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt,
								double * diffPr_max_w_pt, double * diffPr_l2_w_pt, double * diffrelPr_max_w_pt, double * diffrelPr_l2_w_pt,
								double * diffPr_max_T_pt, double * diffPr_l2_T_pt, double * diffrelPr_max_T_pt, double * diffrelPr_l2_T_pt );
int PredictorCorrector3D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, 
										 double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus12, double * Wx_nplus1, 
										 double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus12, double * Wy_nplus1, 
										 double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus12, double * Wz_nplus1, 
										 int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step,
										 double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt,
										 double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt,
				 						 double * diffPr_max_w_pt, double * diffPr_l2_w_pt, double * diffrelPr_max_w_pt, double * diffrelPr_l2_w_pt,
										 double * diffPr_max_T_pt, double * diffPr_l2_T_pt, double * diffrelPr_max_T_pt, double * diffrelPr_l2_T_pt );
int Uzawa3D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus12, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus12, double * Wy_nplus1, double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus1, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt);
int Local1D_3D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus1, double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus1, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt);
int	DifferenceProjectors_calc( double * T_av, double * Wx_av, double * Wy_av, double * Wz_av,
							  double time, int Mx, int My, int Mz, int num_sol, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition,
							  double * diffPr_max_w_pt, double * diffPr_l2_w_pt, double * diffrelPr_max_w_pt, double * diffrelPr_l2_w_pt,
							  double * diffPr_max_T_pt, double * diffPr_l2_T_pt, double * diffrelPr_max_T_pt, double * diffrelPr_l2_T_pt);

int fprintdVec ( char * filename, double * vec, int dim );


int HeatFlux_Init(double * T_0, double * Wx_0, boundaryConditions xCondition, double * Wy_0, boundaryConditions yCondition, double * Wz_0, boundaryConditions zCondition, int numSolution, int Mx, int My, int Mz)
{
	//computes w_0 by solving Aw_0 = BT_0 + g_0
	HeatFlux_Wx0_Init(T_0, Wx_0, xCondition, numSolution, Mx, My, Mz);
	HeatFlux_Wy0_Init(T_0, Wy_0, yCondition, numSolution, Mx, My, Mz);
	HeatFlux_Wz0_Init(T_0, Wz_0, zCondition, numSolution, Mx, My, Mz);
	{;}
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

			wx_0_temp = wx_0_bound(numSolution,0,j1,k1,0*tau);
			wx_1_temp = wx_1_bound(numSolution,Mx-1,j1,k1,0*tau);

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

			Tx_0_temp = Tx_0_bound(numSolution,0,j1,k1,0*tau);
			Tx_1_temp = Tx_1_bound(numSolution,Mx-1,j1,k1,0*tau);

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

	double * alfay = (double *) malloc ((My + 1) * sizeof(double));
	double * betay = (double *) malloc ((My + 1) * sizeof(double));

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

			wy_0_temp = wy_0_bound(numSolution,i1,0,k1,0*tau);
			wy_1_temp = wy_1_bound(numSolution,i1,My-1,k1,0*tau);

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

			Ty_0_temp = Ty_0_bound(numSolution,i1,0,k1,0*tau);
			Ty_1_temp = Ty_1_bound(numSolution,i1,My-1,k1,0*tau);

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

			Ty_0_temp = Ty_0_bound(numSolution,i1,0,k1,0*tau);
			wy_1_temp = wy_1_bound(numSolution,i1,My-1,k1,0*tau);

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

int HeatFlux_Wz0_Init(double * T_0, double * Wz_0, boundaryConditions zCondition, int numSolution, int Mx, int My, int Mz)
{
	double a_i_old, a_iminus1_old, b_i_old, a_0_old, a_n_old, b_0_old, b_n_old;
	double temper_0, temper_k, temper_n;
	double Tz_0_temp, Tz_1_temp, wz_0_temp, wz_1_temp;
	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	double * alfaz = (double *) malloc ((Mz + 1) * sizeof(double));
	double * betaz = (double *) malloc ((Mz + 1) * sizeof(double));

	switch(zCondition)
	{
	case eNeumann:
		//инициализация для потока по z - Wz_n
		//Нейман, z (новое, проверено)
		{;}
		for (int i = 0; i < Ny*Nx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/Ny;
			j1 = i - i1*Ny;

			wz_0_temp = wz_0_bound(numSolution,i1,j1,0,0*tau);
			wz_1_temp = wz_1_bound(numSolution,i1,j1,Mz-1,0*tau);


			betaz[0] = 0 ;
			betaz[1] = wz_0 ;
			alfaz[0] = 0 ;
			alfaz[1] = 0 ;
			for ( int k = 1 ; k < Nz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//			temper_k = T[k - 1 + j1*Nx + k1*Nx*Ny] - T[k + j1*Nx + k1*Nx*Ny];
				temper_k = T_0[i1 + j1*Nx + (k-1)*Nx*Ny] - T_0[i1 + j1*Nx + k*Nx*Ny];

				if (k==Nz-1)
				{
					alfaz[k + 1] = 0;
					betaz[k + 1] = ( ( temper_k - a_i_old*wz_1 ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = ( ( temper_k - a_iminus1_old*wz_0 ) ) / (  b_i_old );
				}
				else
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = ( ( temper_k ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				}
			}
			Wz_0[i*(Nz+1) + Nz ] = wz_1 ;
			for ( int j = 1 ; j < Nz + 1 ; j++ )
			{
				Wz_0[i*(Nz+1) + Nz - j] = Wz_0[i*(Nz+1) + Nz + 1 - j] * alfaz[Nz + 1 - j] + betaz[Nz + 1 - j];
			}
		}
		break;
	case eDirichlet:
		//инициализация для потока по z - Wz_n
		//Дирихле, z (новое, проверено)
		{;}

		for (int i = 0; i < Ny*Nx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/Ny;
			j1 = i - i1*Ny;

			Tz_0_temp = Tz_0_bound(numSolution,i1,j1,0,0*tau);
			Tz_1_temp = Tz_1_bound(numSolution,i1,j1,Mz-1,0*tau);

			//к-ты матрицы А
			a_0_old = hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = Tz_0_temp - T_0[i1 + j1*Nx + 0*Nx*Ny];

			betaz[0] = 0 ;
			alfaz[0] = 0 ;
			alfaz[1] = -0.5 ;
			betaz[1] =  temper_0 / b_0_old; 
			//		printf("temper_0 = %f \n",temper_0);
			for ( int k = 1 ; k < Nz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//			temper_k = T[k - 1 + j1*Nx + k1*Nx*Ny] - T[k + j1*Nx + k1*Nx*Ny];
				temper_k = T_0[i1 + j1*Nx + (k-1)*Nx*Ny] - T_0[i1 + j1*Nx + k*Nx*Ny];

				alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
				betaz[k + 1] =  ( ( temper_k ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				//			printf("temper_k = %f \n",temper_k);
			}
			a_n_old = hz(Nz-1)/(heat_conductivity(i1,j1,Nz-1)*6.0);
			b_n_old = 2*a_n_old;
			temper_n = T_0[i1 + j1*Nx + (Nz-1)*Nx*Ny] - Tz_1_temp;
			//		printf("temper_n = %f \n",temper_n);

			Wz_0[i*(Nz+1) + Nz] = ( temper_n  -  betaz[Nz]*a_n_old ) / (alfaz[Nz]*a_n_old + b_n_old);
			for ( int j = 1 ; j < Nz + 1 ; j++ )
			{
				Wz_0[i*(Nz+1) + Nz - j] = Wz_0[i*(Nz+1) + Nz + 1 - j] * alfaz[Nz + 1 - j] + betaz[Nz + 1 - j];
			}
		}
		break;
	default:
		printf("Bad bnd_Z because of wrong Bnd_z\n");
		return -1;
		break;
	}
	free(alfaz);
	free(betaz);

}


int Accuracy_calculate(double * T_n, double * Wx_n, double * Wy_n, double *Wz_n,
					   int numSolution, int t, int Mx, int My, int Mz, int print_step, 
					   double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt,
					   double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
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
	double tz = 0;
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
	int dim_T = Mx * My * Mz;
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);

	for (int ind = 0; ind < My*Mz; ind++)
	{
		j1 = ind - (ind/My)*My;
		k1 = ind/My;
		for ( int i = 0 ; i < Mx ; i++ )
		{
			temp = 0;
			temp2 = 0;

#ifdef SPECIALTEST
			temp = fabs(T_n[ind*Mx + i]);
			temp2 = temp; 
#else
			temp = fabs(T_n[ind*Mx + i] - exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau));
			temp2 = fabs(exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau));
#endif
			//if (ind == 0 && t == -1)
			//printf("temp = %f \n", temp);
			
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

	//diff_l2z = 0.0;
	//norm_l2_wz = 0.0;
	int imax_z = 0;
	int jmax_z = 0;
	int kmax_z = 0;
	temp2 = 0;

	for (int i=0; i<Ny*Nx; i++)
	{
		for (int k=0; k<Nz+1; k++)
		{
			//i1 = i%Nx;
			//j1 = (i-i%Nx)/Nx;
			i1 = i/Ny;
			j1 = i - i1*Ny;

			k1 = k;

#ifdef SPECIALTEST
			diffz = fabs(Wz_n[i*(Nz+1) + k]);
			temp2 = diffz;
#else
			diffz = fabs(Wz_n[i*(Nz+1) + k] - (-heat_conductivity(i1,j1,k)*exact_gradientZ(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1], (t+1)*tau) ));
			temp2 = fabs(-heat_conductivity(i1,j1,k)*exact_gradientZ(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1], (t+1)*tau) );
#endif
			if (diffz > tz)
			{
				tz = diffz;
				imax_z = i1;
				jmax_z = j1;
				kmax_z = k1;
			}

			if (temp2 > norm_max_w)
				norm_max_w = temp2;

			diff_l2z += diffz*diffz; 
			norm_l2_wz += temp2*temp2;
		}
	}
	diff_l2z = diff_l2z/dim_Wz;
	norm_l2_wz = norm_l2_wz/dim_Wz;

	if ( tz > eps_n1_w)
		eps_n1_w = tz;

	//printf("----------------------\n");
	//printf("after wz_loop eps_n1_w = %f \n", eps_n1_w);
	//printf("after wz_loop diff_l2z = %f \n", diff_l2z);
	//printf("after wz_loop norm_l2_wz = %f \n", norm_l2_wz);


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

#ifdef SPECIALTEST
			diffx = fabs(Wx_n[i*(Mx+1) + k]);
			temp2 = diffx;
#else
			diffx = fabs(Wx_n[i*(Mx+1) + k] - (-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) ));
			temp2 = fabs(-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );
#endif
			//double tempval = (-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );
			//diffx = 1.0;

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

		//printf("after wx_loop eps_n1_w = %f \n", eps_n1_w);
		//printf("after wx_loop diff_l2x = %f \n", diff_l2x);
		//printf("after wx_loop norm_l2_wx = %f \n", norm_l2_wx);


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

#ifdef SPECIALTEST
				diffy = fabs(Wy_n[i*(My+1) + k]);
				temp2 = diffy;
#else
				diffy = fabs(Wy_n[i*(My+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau) ));
				temp2 = fabs(-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );
#endif
				if (diffy > ty)
				{
					ty = diffy;
					imax_y = i1;
					jmax_y = j1;
					kmax_y = k1;
				}

//				printf ( "i1 = %d j1 = %d k1 = %d \n", i1,j1,k1);
//				printf ( "val_numer = %3.3e val_ex = %3.3e \n", Wy_n[i*(My+1) + k], (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau) ) ); 

				if (temp2 > norm_max_w)
					norm_max_w = temp2;

				diff_l2y += diffy*diffy; 
				norm_l2_wy += temp2*temp2;
			}
			diff_l2y = diff_l2y/dim_Wy;
			norm_l2_wy = norm_l2_wy/dim_Wy;

			if ( ty > eps_n1_w)
				eps_n1_w = ty;

			//printf("after wy_loop eps_n1_w = %f \n", eps_n1_w);
			//printf("after wy_loop diff_l2z = %f \n", diff_l2y);
			//printf("after wy_loop norm_l2_wz = %f \n", norm_l2_wy);


			eps_wn2 = sqrt(diff_l2x + diff_l2y + diff_l2z);
			norm_l2_w = sqrt(norm_l2_wx + norm_l2_wy + norm_l2_wz);

			//printf("after all eps_wn2 = %f \n", eps_wn2);
			//printf("after all norm_l2_w = %f \n", norm_l2_w);
			//printf("----------------------\n");

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
			if (t == -1)
			{
				printf("t=%d inside function Accuracy \n",t);
				printf("eps_max_T = %e \n",eps_max_T); 
				printf("eps_l2_T = %e \n",eps_l2_T); 
				printf("eps_max_w = %e \n",eps_max_w); 
				printf("eps_l2_w = %e \n",eps_l2_w); 

				printf("eps_wx_l2 = %f \n", sqrt(diff_l2x));
				printf("eps_wy_l2 = %f \n", sqrt(diff_l2y));
				printf("eps_wz_l2 = %f \n", sqrt(diff_l2z));

			}
			if (t+1-((t+1)/print_step)*print_step==0 && t != -1)
			{
				printf("t=%d \n",t);
				printf("eps_max_T = %e \n",eps_max_T); 
				printf("eps_l2_T = %e \n",eps_l2_T); 
				printf("eps_max_w = %e \n",eps_max_w); 
				printf("eps_l2_w = %e \n",eps_l2_w); 

				printf("eps_wx_l2 = %f \n", sqrt(diff_l2x));
				printf("eps_wy_l2 = %f \n", sqrt(diff_l2y));
				printf("eps_wz_l2 = %f \n", sqrt(diff_l2z));

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


int Righthand_step1_wynplus13(double * righthandY, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wy; i++)
	{
		righthandY[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
//	fprintdVec ( "F_n_wy_nplus13_3d.txt", F, dim_T );
	ByMminus1_F_Add(righthandY, F, yCondition, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus13_1part_3d.txt", righthandY, dim_Wy );
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	ByMBtr_W_Add(righthandY, F, coeff_BMBWx, Wx_n, xCondition, coeff_BMBWy, Wy_n, yCondition, coeff_BMBWz, Wz_n, zCondition, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus13_12parts_3d.txt", righthandY, dim_Wy );
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	ByMC_W_Add(righthandY, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus13_123parts_3d.txt", righthandY, dim_Wy );
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wznplus13(double * righthandZ, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wz; i++)
	{
		righthandZ[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	BzMminus1_F_Add(righthandZ, F, zCondition, Mx, My, Mz);

//	fprintdVec ( "righthandZ_1.txt", righthandZ, dim_Wz);

	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	BzMBtr_W_Add(righthandZ, F, coeff_BMBWx, Wx_n, xCondition, coeff_BMBWy, Wy_n, yCondition, coeff_BMBWz, Wz_n, zCondition, Mx, My, Mz);

//	fprintdVec ( "righthandZ_2.txt", righthandZ, dim_Wz);

	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	BzMC_W_Add(righthandZ, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);

//	fprintdVec ( "righthandZ_3.txt", righthandZ, dim_Wz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wxnplus1(double * righthandX, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wx; i++)
	{
		righthandX[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	BxMminus1_F_Add(righthandX, F, xCondition, Mx, My, Mz);
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	BxMBtr_W_Add(righthandX, F, coeff_BMBWx, Temp_Wx, xCondition, coeff_BMBWy, Temp_Wy, yCondition, coeff_BMBWz, Temp_Wz, zCondition, Mx, My, Mz);
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	BxMC_W_Add(righthandX, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wxnplus1_local(double * righthandX, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wx; i++)
	{
		righthandX[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	BxMminus1_F_Add(righthandX, F, xCondition, Mx, My, Mz);
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	BxMBtr_W_Add(righthandX, F, coeff_BMBWx, Wx_n, xCondition, coeff_BMBWy, Wy_n, yCondition, coeff_BMBWz, Wz_n, zCondition, Mx, My, Mz);
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	BxMC_W_Add(righthandX, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}
int Righthand_step1_wynplus1(double * righthandY, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wy; i++)
	{
		righthandY[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
//	fprintdVec ( "F_wy_nplus1_3d.txt", F, dim_T );
	ByMminus1_F_Add(righthandY, F, yCondition, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus1_1part_3d.txt", righthandY, dim_Wy );
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	ByMBtr_W_Add(righthandY, F, coeff_BMBWx, Temp_Wx, xCondition, coeff_BMBWy, Temp_Wy, yCondition, coeff_BMBWz, Temp_Wz, zCondition, Mx, My, Mz);
//	fprintdVec ( "Temp_Wy_wy_nplus1_12parts_3d.txt", Temp_Wy, dim_Wy );
//	fprintdVec ( "rhand_wy_nplus1_12parts_3d.txt", righthandY, dim_Wy );
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	ByMC_W_Add(righthandY, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus1_123parts_3d.txt", righthandY, dim_Wy );
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wynplus1_local(double * righthandY, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double * Wx_nplus1, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wy; i++)
	{
		righthandY[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	ByMminus1_F_Add(righthandY, F, yCondition, Mx, My, Mz);
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	ByMBtr_W_Add(righthandY, F, coeff_BMBWx, Wx_nplus1, xCondition, coeff_BMBWy, Wy_n, yCondition, coeff_BMBWz, Wz_n, zCondition, Mx, My, Mz);
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	ByMC_W_Add(righthandY, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wznplus1(double * righthandZ, double * F, double alpha, double beta, 
							 double * Vx, double coeff_Vx, double * Wx_n, double coeff_BMBWx, boundaryConditions xCondition, 
							 double * Vy, double coeff_Vy, double * Wy_n, double coeff_BMBWy, boundaryConditions yCondition, 
							 double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, 
							 double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wz; i++)
	{
		righthandZ[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	BzMminus1_F_Add(righthandZ, F, zCondition, Mx, My, Mz);
	//fprintdVec ( "Righthand_Wz_nplus12_part1_3d.txt", righthandZ, dim_Wz);
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	BzMBtr_W_Add(righthandZ, F, coeff_BMBWx, Temp_Wx, xCondition, coeff_BMBWy, Temp_Wy, yCondition, coeff_BMBWz, Temp_Wz, zCondition, Mx, My, Mz);
	//fprintdVec ( "Righthand_Wz_nplus12_parts12_3d.txt", righthandZ, dim_Wz);
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	BzMC_W_Add(righthandZ, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}

int Righthand_step1_wznplus1_local(double * righthandZ, double * F, double alpha, double beta, double * Vx, double coeff_Vx, double * Wx_n, double * Wx_nplus1, double coeff_BMBWx, boundaryConditions xCondition, double * Vy, double coeff_Vy, double * Wy_n, double *Wy_nplus1, double coeff_BMBWy, boundaryConditions yCondition, double * Vz, double coeff_Vz, double * Wz_n, double coeff_BMBWz, boundaryConditions zCondition, int numSolution, int t, int Mx, int My, int Mz)
{
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_Wz; i++)
	{
		righthandZ[i] = 0.;
	}
	F_alpha_beta_fill(F, t, alpha, beta, numSolution, Mx, My, Mz);
	BzMminus1_F_Add(righthandZ, F, zCondition, Mx, My, Mz);
	//obi4no coeff_BMBWx = -1, coeff_BMBWy = -1;
	BzMBtr_W_Add(righthandZ, F, coeff_BMBWx, Wx_nplus1, xCondition, coeff_BMBWy, Wy_nplus1, yCondition, coeff_BMBWz, Wz_n, zCondition, Mx, My, Mz);
	//obi4no coeff_Vx = -1, coeff_Vy = -1;
	// UNCHECKED BY NOW
	BzMC_W_Add(righthandZ, F, coeff_Vx, Vx, Wx_n, xCondition, coeff_Vy, Vy, Wy_n, yCondition, coeff_Vz, Vz, Wz_n, zCondition, Mx, My, Mz);
	//F is spoiled because it is used as Temp during call of ByMBtr_W_Add
	return 0;
}
int F_alpha_beta_fill(double * F, int t, double alpha, double beta, int numSolution, int Mx, int My, int Mz)
{
	if (alpha*alpha + beta*beta != 0)
	{
		//usually we use (alpha, beta) = (1, 0)[explicit], (0, 1)[implicit], or (0.5,0.5)[Crank-Nicolson]
		for (int k = 0; k < Mz; k++ )
		{
			for (int j = 0; j < My; j++ )
			{
				for (int i = 0; i < Mx; i++)
				{
					F[i + j * Mx + k * Mx * My] = alpha * righthand(numSolution, t, i, j, k, xpoints, ypoints, zpoints) + beta * righthand(numSolution, t + 1, i, j, k, xpoints, ypoints, zpoints);
				}
			}
		}
	}
	else
		for (int i = 0; i < Mx * My * Mz; i++)
			F[i] = 0.0;
	return 0;
}

int BzMminus1_F_Add(double * Output, double * F, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	double m_ijk, m_ijkminus1;
	int n_block_z;

	//zyx
	for (int i = 0; i < Mx; i++)
	{
		for (int j = 0; j < My; j++)
		{
			n_block_z = j + i * My;

			for (int k = 0; k < Mz + 1; k++)
			{
				m_ijk = heat_capacity(i,j,k)*density(i,j,k)*hy(j)*hx(i)*hz(k);
				m_ijkminus1 = heat_capacity(i,j,k-1)*density(i,j,k-1)*hy(j)*hx(i)*hz(k-1);

				if (k == 0)
				{
					if (zCondition == eNeumann)
						Output[n_block_z * (Mz + 1) + k] += 0;
					else
					{
						Output[n_block_z * (Mz + 1) + k] += - hx(i) * hy(j) * (1.0/m_ijk) * F[k * Mx * My + j*Mx + i];
					}
				}
				else
					if (k == Mz)
					{
						if (zCondition == eNeumann || zCondition == eMixed)
							Output[n_block_z * (Mz + 1) + k] += 0;
						else
						{
							Output[n_block_z * (Mz + 1) + k] += hx(i) * hy(j) * (1.0/m_ijkminus1) * F[(k-1) * Mx * My + j * Mx + i];
						}
					}
					else //0 < k < Mz
					{
						Output[n_block_z * (Mz + 1) + k] += hx(i) * hy(j) * ((1.0/m_ijkminus1) * F[(k - 1) * Mx * My + j * Mx + i] - (1.0/m_ijk) * F[k * Mx * My + j*Mx + i] );
					}

				//if ( n_block_z <= 2 )
				//{
				//	printf ("nblock = %d Outpt[%d] = %f \n", n_block_z, n_block_z * (Mz + 1) + k, Output[n_block_z * (Mz + 1) + k] );
				//	if ( k == 0 )
				//		printf ( "F[%d] = %f \n", k * Mx * My + j*Mx + i, F[k * Mx * My + j*Mx + i] );
				//	else if ( k == Mz )
				//		printf ( "F[%d] = %f \n", (k-1) * Mx * My + j*Mx + i, F[(k-1) * Mx * My + j*Mx + i] );
				//	else
				//		printf ( "diff F = [k-1] - [k] = %f - %f = %f \n", F[(k - 1) * Mx * My + j * Mx + i], F[k * Mx * My + j*Mx + i], F[(k - 1) * Mx * My + j * Mx + i] - F[k * Mx * My + j*Mx + i]);
				//}

			}

		}
	}
	return 0;
}


int ByMminus1_F_Add(double * Output, double * F, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	double m_ijk, m_ijminus1k;
	int n_block_y;

	for (int k = 0; k < Mz; k++)
	{
		for (int i = 0; i < Mx; i++)
		{
			n_block_y = k + i * Mz;

			for (int j = 0; j < My + 1; j++)
			{
				m_ijk = heat_capacity(i,j,k)*density(i,j,k)*hy(j)*hx(i)*hz(k);
				m_ijminus1k = heat_capacity(i,j-1,k)*density(i,j-1,k)*hy(j-1)*hx(i)*hz(k);

				if (j == 0)
				{
					if (yCondition == eNeumann)
						Output[n_block_y * (My + 1) + j] += 0;
					else
					{
						//temp_0 = hx(i1)*(1.0/m_0)*righthand(numSolution,t,i1,0,k1,xpoints,ypoints,zpoints);
						//f_0 = -temp_0 ; //вычислили B* (M(-1) * Fn+1)

						Output[n_block_y * (My + 1) + j] += - hx(i) * hz(k) * (1.0/m_ijk) * F[k * Mx * My + j*Mx + i];
					}
				}
				else
					if (j == My)
					{
						if (yCondition == eNeumann || yCondition == eMixed)
							Output[n_block_y * (My + 1) + j] += 0;
						else
						{
							//temp_nminus1 = hx(i1)*(1.0/m_nminus1)*righthand(numSolution,t,i1,Ny-1,k1,xpoints,ypoints,zpoints);
							//f_n = temp_nminus1; //вычислили B* (M(-1) * Fn+1)

							Output[n_block_y * (My + 1) + j] += hx(i) * hz(k) * (1.0/m_ijminus1k) * F[k * Mx * My + (j-1) * Mx + i];
						}
					}
					else //0 < j < Ny
					{
						//temp_k = hx(i1)*(1.0/m_i)*righthand(numSolution,t,i1,k,k1,xpoints,ypoints,zpoints);
						//temp_kminus1 = hx(i1)*(1.0/m_iminus1)*righthand(numSolution,t,i1,k-1,k1,xpoints,ypoints,zpoints);
						//f_k = ( temp_kminus1 - temp_k ); //вычислили B* M(-1) * 0.5 *(Fn+1 +  Fn)

						Output[n_block_y * (My + 1) + j] += hx(i) * hz(k) * ((1.0/m_ijminus1k) * F[k * Mx * My + (j - 1) * Mx + i] - (1.0/m_ijk) * F[k * Mx * My + j*Mx + i] );
						/*if (n_block_y == 0 && j == 1)
						{
							printf("val_j-1 = %f val_j = %f \n",F[k * Mx * My + (j - 1) * Mx + i], F[k * Mx * My + j * Mx + i]);
							printf("index1 = %d index2 = %d \n",k * Mx * My + (j - 1) * Mx + i, k * Mx * My + j * Mx + i);
						}*/
					}
			}
		}
	}
	return 0;
}

int BxMminus1_F_Add(double * Output, double * F, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	double m_ijk, m_iminus1jk;
	int n_block_x;

	//xyz
	for (int k = 0; k < Mz; k++)
	{
		for (int j = 0; j < My; j++)
		{
			n_block_x = j + k * My;

			for (int i = 0; i < Mx + 1; i++)
			{
				m_ijk = heat_capacity(i,j,k)*density(i,j,k)*hy(j)*hx(i)*hz(k);
				m_iminus1jk = heat_capacity(i-1,j,k)*density(i-1,j,k)*hy(j)*hx(i-1)*hz(k);

				if (i == 0)
				{
					if (xCondition == eNeumann)
						Output[n_block_x * (Mx + 1) + i] += 0;
					else
					{
						Output[n_block_x * (Mx + 1) + i] += - hy(j) * hz(k) * (1.0/m_ijk) * F[k * Mx * My + j*Mx + i];
					}
				}
				else
					if (i == Mx)
					{
						if (xCondition == eNeumann || xCondition == eMixed)
							Output[n_block_x * (Mx + 1) + i] += 0;
						else
						{
							Output[n_block_x * (Mx + 1) + i] += hy(j) * hz(k) * (1.0/m_iminus1jk) * F[k * Mx * My + j * Mx + (i - 1)];
						}
					}
					else //0 < i < Nx
					{
						Output[n_block_x * (Mx + 1) + i] += hy(j) * hz(k) * ((1.0/m_iminus1jk) * F[k * Mx * My + j * Mx + i - 1] - (1.0/m_ijk) * F[k * Mx * My + j*Mx + i] );
					}
			}
		}
	}
	return 0;
}

int ByMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	Btr_W_Add(Temp, alpha, Wx_n, xCondition, beta, Wy_n, yCondition, gamma, Wz_n, zCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	ByMminus1_F_Add(Output, Temp, yCondition, Mx, My, Mz);

	return 0;
}

int BzMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	Btr_W_Add(Temp, alpha, Wx_n, xCondition, beta, Wy_n, yCondition, gamma, Wz_n, zCondition, Mx, My, Mz);

//	fprintdVec ( "Temp.txt", Temp, dim_T );
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	BzMminus1_F_Add(Output, Temp, zCondition, Mx, My, Mz);

	return 0;
}

int BxMBtr_W_Add(double * Output, double * Temp, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	Btr_W_Add(Temp, alpha, Wx_n, xCondition, beta, Wy_n, yCondition, gamma, Wz_n, zCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	BxMminus1_F_Add(Output, Temp, xCondition, Mx, My, Mz);

	return 0;
}

int Btr_W_Add(double * Output, double alpha, double * Wx_n, boundaryConditions xCondition, double beta, double * Wy_n, boundaryConditions yCondition, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	Bxtr_Wx_Add(Output, alpha, Wx_n, xCondition, Mx, My, Mz);
	//printf("Bxtr_Wx_Add[0] = %f \n", Output[0]);
	Bytr_Wy_Add(Output, beta, Wy_n, yCondition, Mx, My, Mz);
	//printf("Bxtr_Wx + Bytr_Wy[0] = %f \n", Output[0]);
	Bztr_Wz_Add(Output, gamma, Wz_n, zCondition, Mx, My, Mz);
	return 0;
}

int Bxtr_Wx_Add(double * Output, double alpha, double * Wx_n, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	for (int k = 0; k < Mz; k++)
		for (int j = 0; j < My; j++)
			for (int i = 0; i < Mx; i++)
				Output[k * Mx * My + j*Mx + i] += alpha * hy(j) * hz(k) * (Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i + 1] - Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i]);
	return 0;
}

int Bytr_Wy_Add(double * Output, double beta, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	for (int i = 0; i < Mx; i++)
		for (int k = 0; k < Mz; k++)
			for (int j = 0; j < My; j++)
				Output[k * Mx * My + j*Mx + i] += beta * hx(i) * hz(k) * (Wy_n[i * (My + 1) * Mz + k * (My + 1) + j + 1] - Wy_n[i * (My + 1) * Mz + k * (My + 1) + j]);
	return 0;
}

int Bztr_Wz_Add(double * Output, double gamma, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	//zxy
	//for (int j = 0; j < My; j++)
	//	for (int i = 0; i < Mx; i++)
	//		for (int k = 0; k < Mz; k++)
	//			Output[k * Mx * My + j*Mx + i] += gamma * hx(i) * hy(j) * (Wz_n[j * Mx * (Mz + 1) + i * (Mz + 1) + k + 1] - Wz_n[j * Mx * (Mz + 1) + i * (Mz + 1) + k]);
	//zyx
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My; j++)
			for (int k = 0; k < Mz; k++)
				Output[k * Mx * My + j*Mx + i] += gamma * hx(i) * hy(j) * (Wz_n[i * My * (Mz + 1) + j * (Mz + 1) + k + 1] - Wz_n[i * My * (Mz + 1) + j * (Mz + 1) + k]);
	return 0;
}

int ByMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	C_W_Add(Temp, alpha, Vx, Wx_n, xCondition, beta, Vy, Wy_n, yCondition, gamma, Vz, Wz_n, zCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	ByMminus1_F_Add(Output, Temp, yCondition, Mx, My, Mz);

	return 0;
}

int BzMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	C_W_Add(Temp, alpha, Vx, Wx_n, xCondition, beta, Vy, Wy_n, yCondition, gamma, Vz, Wz_n, zCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	BzMminus1_F_Add(Output, Temp, yCondition, Mx, My, Mz);

	return 0;
}

int BxMC_W_Add(double * Output, double * Temp, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	int dim_T = Mx * My * Mz;
	for (int i = 0; i < dim_T; i++)
		Temp[i] = 0.0;
	C_W_Add(Temp, alpha, Vx, Wx_n, xCondition, beta, Vy, Wy_n, yCondition, gamma, Vz, Wz_n, zCondition, Mx, My, Mz);
	//printf("Btr_W_Add[0] = %f \n", Temp[0]);
	BxMminus1_F_Add(Output, Temp, xCondition, Mx, My, Mz);

	return 0;
}

int C_W_Add(double * Output, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	C_Wx_Add(Output, alpha, Vx, Wx_n, xCondition, Mx, My, Mz);
	//printf("Bxtr_Wx_Add[0] = %f \n", Output[0]);
	C_Wy_Add(Output, beta, Vy, Wy_n, yCondition, Mx, My, Mz);
	//printf("Bxtr_Wx + Bytr_Wy[0] = %f \n", Output[0]);
	C_Wz_Add(Output, gamma, Vz, Wz_n, zCondition, Mx, My, Mz);
	return 0;
}

int C_Wx_Add(double * Output, double alpha, double * Vx, double * Wx_n, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	double cx_ijk, cx_iplus1jk, koeff_ijk, a_ijk_old;
	for (int k = 0; k < Mz; k++)
		for (int j = 0; j < My; j++)
			for (int i = 0; i < Mx; i++)
			{
				a_ijk_old = hx(i)*hy(j)*hz(k)/(heat_conductivity(i,j,k)*6.0);
				koeff_ijk = heat_capacity(i,j,k)*density(i,j,k)*a_ijk_old;

				cx_ijk = koeff_ijk * (2 * Vx[k * My * (Mx + 1) + j * (Mx + 1) + i] + Vx[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);
				cx_iplus1jk = koeff_ijk * (Vx[k * My * (Mx + 1) + j * (Mx + 1) + i] + 2 * Vx[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);

				//Output[k * Mx * My + j*Mx + i] += alpha * hy(j) * hz(k) * (cx_ijk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i] + cx_iplus1jk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);
				Output[k * Mx * My + j*Mx + i] += alpha * (cx_ijk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i] + cx_iplus1jk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);
				//if (k * Mx * My + j*Mx + i == 0 || k * Mx * My + j*Mx + i == 8)
				//{
				//	printf("Output[%d] = %f \n",k * Mx * My + j*Mx + i, Output[k * Mx * My + j*Mx + i]);
				//	printf("cx_ijk = %f cx_iplus1jk = %f \n",cx_ijk, cx_iplus1jk);
				//	printf("Wx_n[%d] = %f Wx_n[%d] = %f \n",k * My * (Mx + 1) + j * (Mx + 1) + i, Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i],k * My * (Mx + 1) + j * (Mx + 1) + i + 1, Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);
				//	printf("their product1 = %f product2 = %f \n",cx_ijk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i], cx_iplus1jk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);
				//	printf("their sum = %f \n",cx_ijk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i] + cx_iplus1jk*Wx_n[k * My * (Mx + 1) + j * (Mx + 1) + i + 1]);
				//}
			}
	return 0;
}

int C_Wy_Add(double * Output, double beta, double * Vy, double * Wy_n, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	double cy_ijk, cy_ijplus1k, koeff_ijk, a_ijk_old;
	for (int i = 0; i < Mx; i++)
		for (int k = 0; k < Mz; k++)
			for (int j = 0; j < My; j++)
			{
				a_ijk_old = hx(i)*hy(j)*hz(k)/(heat_conductivity(i,j,k)*6.0);
				koeff_ijk = heat_capacity(i,j,k)*density(i,j,k)*a_ijk_old;

				cy_ijk = koeff_ijk * (2 * Vy[i * (My + 1) * Mz + k * (My + 1) + j] + Vy[i * (My + 1) * Mz + k * (My + 1) + j + 1]);
				cy_ijplus1k = koeff_ijk * (Vy[i * (My + 1) * Mz + k * (My + 1) + j] + 2 * Vy[i * (My + 1) * Mz + k * (My + 1) + j + 1]);

				Output[k*Mx*My + j*Mx + i] += beta * (cy_ijk * Wy_n[i * (My + 1) * Mz + k * (My + 1) + j] + cy_ijplus1k * Wy_n[i * (My + 1) * Mz + k * (My + 1) + j + 1]);
			}
		return 0;
}


int C_Wz_Add(double * Output, double gamma, double * Vz, double * Wz_n, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	double cz_ijk, cz_ijplus1k, koeff_ijk, a_ijk_old;
	//zxy
	//for (int j = 0; j < My; j++)
	//	for (int i = 0; i < Mx; i++)
	//		for (int k = 0; k < Mz; k++)
	//		{
	//			a_ijk_old = hx(i)*hy(j)*hz(k)/(heat_conductivity(i,j,k)*6.0);
	//			koeff_ijk = heat_capacity(i,j,k)*density(i,j,k)*a_ijk_old;

	//			cz_ijk = koeff_ijk * (2 * Vz[j * Mx * (Mz + 1) + i * (Mz + 1) + k] + Vz[j * Mx * (Mz + 1) + i * (Mz + 1) + k + 1]);
	//			cz_ijplus1k = koeff_ijk * (Vz[j * Mx * (Mz + 1) + i * (Mz + 1) + k] + 2 * Vz[j * Mx * (Mz + 1) + i * (Mz + 1) + k + 1]);

	//			Output[k*Mx*My + j*Mx + i] += gamma * hx(i) * hy(j) * (cz_ijk * Wz_n[j * Mx * (Mz + 1) + i * (Mz + 1) + k] + cz_ijplus1k * Wz_n[j * Mx * (Mz + 1) + i * (Mz + 1) + k + 1]);
	//		}
	//zyx
	for (int i = 0; i < Mx; i++)
		for (int j = 0; j < My; j++)
			for (int k = 0; k < Mz; k++)
			{
				a_ijk_old = hx(i)*hy(j)*hz(k)/(heat_conductivity(i,j,k)*6.0);
				koeff_ijk = heat_capacity(i,j,k)*density(i,j,k)*a_ijk_old;

				cz_ijk = koeff_ijk * (2 * Vz[i * My * (Mz + 1) + j * (Mz + 1) + k] + Vz[i * My * (Mz + 1) + j * (Mz + 1) + k + 1]);
				cz_ijplus1k = koeff_ijk * (Vz[i * My * (Mz + 1) + j * (Mz + 1) + k] + 2 * Vz[i * My * (Mz + 1) + j * (Mz + 1) + k + 1]);

				Output[k*Mx*My + j*Mx + i] += gamma * (cz_ijk * Wz_n[i * My * (Mz + 1) + j * (Mz + 1) + k] + cz_ijplus1k * Wz_n[i * My * (Mz + 1) + j * (Mz + 1) + k + 1]);
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
				a_i_old = hy(k)*hx(i1)*hz(k1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,k-1,k1)*6.0);
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
			a_0_old = hy(0)*hx(i1)*hz(k1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;

			betay[0] = 0 ;
			betay[1] = righthandY[i * (My + 1) + 0] / b_0_old ;
			alfay[0] = 0 ;
			alfay[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)*hz(k1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
				betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1_old ) / (alfay[k]*a_iminus1_old + b_i_old);
			}

			//к-ты матрицы А
			a_nminus1_old = hy(My-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			Wy_nplus05[i*(My+1) +  My] = (   righthandY[i * (My + 1) + My]  - betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);
//			if (i == 0)
//				printf("delta_Wy_nplus05[%d] = %f \n",i*(My+1) +  My,Wy_nplus05[i*(My+1) +  My]);
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus05[i*(My+1) + My - j] = Wy_nplus05[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
//				if (i == 0)
//					printf("delta_Wy_nplus05[%d] = %f \n",i*(My+1) + My - j,Wy_nplus05[i*(My+1) + My - j]);
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
			a_0_old = hy(0)*hx(i1)*hz(k1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;

			betay[0] = 0 ;
			betay[1] =  righthandY[i * (My + 1) + 0]/b_0_old ;
			alfay[0] = 0 ;
			alfay[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)*hz(k1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,k-1,k1)*6.0);
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

	int dim_Wy = Mx * (My + 1) * Mz;
	for (int i = 0; i < dim_Wy; i++)
	{
		Wy_nplus05[i] = Wy_n[i] + sigma*tau*Wy_nplus05[i];
#ifdef ZEROING
		Wy_nplus05[i] = 0;
#endif
	}

	return 0;
}



int Explicit_Wz(double * Wz_nplus05, double * Wz_n, double * righthandZ, boundaryConditions zCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;

	double * alfaz = (double*) malloc ((Mz + 1) * sizeof(double));
	double * betaz = (double*) malloc ((Mz + 1) * sizeof(double));

	switch(zCondition)
	{
	case eNeumann:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ z-компоненты
		//Нейман, z
		for (int i = 0; i < My*Mx; i++)
		{
			//zxy
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			//zyx
			i1 = i/My;
			j1 = i - i1*My;

			betaz[0] = 0 ;
			betaz[1] = righthandZ[i * (Mz + 1) + 0];
			alfaz[0] = 0 ;
			alfaz[1] = 0 ;
			for ( int k = 1 ; k < Mz ; k++ )
			{
				a_i_old = hx(i1)*hy(j1)*hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hx(i1)*hy(j1)*hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==Mz-1)
				{
					alfaz[k + 1] = 0;
					//right
					betaz[k + 1] = (  righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1_old ) / (alfaz[k]*a_iminus1_old + b_i_old);
				}
				else if (k==1)
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = ( righthandZ[i * (Mz + 1) + k] ) / ( b_i_old );
				}
				else
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1_old ) / (alfaz[k]*a_iminus1_old + b_i_old);
				}

			}

			Wz_nplus05[i*(Mz+1) + Mz] =  righthandZ[i * (Mz + 1) + Mz] ;
			for ( int j = 1 ; j < Mz + 1 ; j++ )
			{
				Wz_nplus05[i*(Mz+1) + Mz - j] = Wz_nplus05[i*(Mz+1) + Mz + 1 - j] * alfaz[Mz + 1 - j] + betaz[Mz + 1 - j];
			}
		}
		break;
	case eDirichlet:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ y-компоненты
		//Дирихле, y
		for (int i = 0; i < My*Mx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/My;
			j1 = i - i1*My;

			//к-ты матрицы А
			a_0_old = hx(i1)*hy(j1)*hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;

			betaz[0] = 0 ;
			betaz[1] = righthandZ[i * (Mz + 1) + 0] / b_0_old ;
			alfaz[0] = 0 ;
			alfaz[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < Mz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(i1)*hy(j1)*hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hx(i1)*hy(j1)*hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
				betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1_old ) / (alfaz[k]*a_iminus1_old + b_i_old);
			}

			//к-ты матрицы А
			a_nminus1_old = hz(Mz-1)*hx(i1)*hy(j1)/(heat_conductivity(i1,j1,Mz-1)*6.0);
			b_n_old = 2*a_nminus1_old;

			Wz_nplus05[i*(Mz+1) +  Mz] = (   righthandZ[i * (Mz + 1) + Mz]  - betaz[Mz]*a_nminus1_old ) / (alfaz[Mz]*a_nminus1_old + b_n_old);

			for ( int j = 1 ; j < Mz + 1 ; j++ )
			{
				Wz_nplus05[i*(Mz+1) + Mz - j] = Wz_nplus05[i*(Mz+1) + Mz + 1 - j] * alfaz[Mz + 1 - j] + betaz[Mz + 1 - j];
//				if (i == 0)
//					printf("delta_Wy_nplus05[%d] = %f \n",i*(My+1) + My - j,Wy_nplus05[i*(My+1) + My - j]);
			}
		}
		break;
	case eMixed:
		//Дирихле-Нейман, y
		{;}

		for (int i = 0; i < My*Mx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/My;
			j1 = i - i1*My;

			//к-ты матрицы А
			a_0_old = hx(i1)*hy(j1)*hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;

			betaz[0] = 0 ;
			betaz[1] =  righthandZ[i * (Mz + 1) + 0]/b_0_old ;
			alfaz[0] = 0 ;
			alfaz[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < Mz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(i1)*hy(j1)*hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hx(i1)*hy(j1)*hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==Mz-1)
				{
					alfaz[k + 1] = 0;
					betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1_old ) / (alfaz[k]*a_iminus1_old + b_i_old);
				}
				else
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1_old ) / (alfaz[k]*a_iminus1_old + b_i_old);
				}
			}
			Wz_nplus05[i*(Mz+1) + Mz] = righthandZ[i*(Mz + 1) + Mz] ;

			for ( int j = 1 ; j < Mz + 1 ; j++ )
			{
				Wz_nplus05[i*(Mz+1) + Mz - j] = Wz_nplus05[i*(Mz+1) + Mz + 1 - j] * alfaz[Mz + 1 - j] + betaz[Mz + 1 - j];
			}
		}
		break;
	}
	{;}

	free(alfaz);
	free(betaz);

	int dim_Wz = Mx * My * (Mz + 1);
	for (int i = 0; i < dim_Wz; i++)
	{
		Wz_nplus05[i] = Wz_n[i] + sigma*tau*Wz_nplus05[i];
#ifdef ZEROING
		Wz_nplus05[i] = 0;
#endif
	}

	return 0;
}


int Explicit_Wx(double * Wx_nplus05, double * Wx_n, double * righthandX, boundaryConditions xCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;

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
				a_i_old = hx(k)*hy(j1)*hz(k1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)*hz(k1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==Mx-1)
				{
					alfax[k + 1] = 0;
					betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1_old ) / (alfax[k]*a_iminus1_old + b_i_old);
				}
				else if (k==1)
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( righthandX[i * (Mx + 1) + k] ) / ( b_i_old );
				}
				else
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1_old ) / (alfax[k]*a_iminus1_old + b_i_old);
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
			a_0_old = hx(0)*hy(j1)*hz(k1)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;

			betax[0] = 0 ;
			betax[1] = righthandX[i * (Mx + 1) + 0] / b_0_old ;
			alfax[0] = 0 ;
			alfax[1] = -a_0_old/b_0_old ;

			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)*hz(k1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)*hz(k1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
				betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1_old ) / (alfax[k]*a_iminus1_old + b_i_old);

			}

			//к-ты матрицы А
			a_nminus1_old = hx(Mx-1)*hy(j1)*hz(k1)/(heat_conductivity(Mx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			Wx_nplus05[i*(Mx+1) +  Mx] = (  righthandX[i * (Mx + 1) + Mx]  - betax[Mx]*a_nminus1_old ) / (alfax[Mx]*a_nminus1_old + b_n_old);
			
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

	int dim_Wx = (Mx + 1) * My * Mz;
	for (int i = 0; i < dim_Wx; i++)
	{
		//if (i < 4)
		//{
		//	printf("Wx_n[i] = %f \nWx_nplus05[i] = %f \n", Wx_n[i], Wx_nplus05[i]);
		//}
		Wx_nplus05[i] = Wx_n[i] + sigma*tau*Wx_nplus05[i];
		//if (i < 4)
		//{
		//	printf("Wx_nplus05_FINAL[i] = %f \n", Wx_nplus05[i]);
		//}
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
				a_i_old = hx(k)*hy(j1)*hz(k1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)*hz(k1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(k,j1,k1)*density(k,j1,k1)*hx(k)*hy(j1)*hz(k1);
				m_iminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*hx(k-1)*hy(j1)*hz(k1);

				//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
				a_i = a_i_old - 0.5*tau*hy(j1)*hz(k1)*hy(j1)*hz(k1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hy(j1)*hz(k1)*hy(j1)*hz(k1)*( 1.0/m_i + 1.0/m_iminus1 );

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
			a_0_old = hx(0)*hy(j1)*hz(k1)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы M
			m_0 = heat_capacity(0,j1,k1)*hz(k1)*density(0,j1,k1)*hx(0)*hy(j1);

			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_0); 

			betax[0] = 0 ;
			betax[1] = righthandX[i * (Mx + 1) + 0] / b_0 ;
			alfax[0] = 0 ;
			alfax[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)*hz(k1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)*hz(k1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы M
				m_i = heat_capacity(k,j1,k1)*density(k,j1,k1)*hx(k)*hy(j1)*hz(k1);
				m_iminus1 = heat_capacity(k-1,j1,k1)*density(k-1,j1,k1)*hx(k-1)*hy(j1)*hz(k1);

				//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
				a_i = a_i_old - 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*( 1.0/m_i + 1.0/m_iminus1 );

				alfax[k + 1] = (-1)*a_i/(alfax[k]*a_iminus1 + b_i);
				betax[k + 1] = (  righthandX[i * (Mx + 1) + k]  - betax[k]*a_iminus1 ) / (alfax[k]*a_iminus1 + b_i);

			}

			//к-ты матрицы А
			a_nminus1_old = hx(Mx-1)*hy(j1)*hz(k1)/(heat_conductivity(Mx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			//к-ты матрицы mx
			m_nminus1 = heat_capacity(Mx-1,j1,k1)*density(Mx-1,j1,k1)*hy(j1)*hx(Mx-1)*hz(k1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hy(j1)*hy(j1)*hz(k1)*hz(k1)*(1.0/m_nminus1);

			Wx_nplus05[i*(Mx+1) +  Mx] = (  righthandX[i * (Mx + 1) + Mx]  - betax[Mx]*a_nminus1 ) / (alfax[Mx]*a_nminus1 + b_n);

			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				Wx_nplus05[i*(Mx+1) + Mx - j] = Wx_nplus05[i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 - j];
				//if (i == 0)
				//	printf("Wx_nplus05[%d] = %f \n",i*(Mx+1) + Mx - j, Wx_nplus05[i*(Mx+1) + Mx - j]);
			}

		}
		break;
	}

	{;}

	free(alfax);
	free(betax);

	int dim_Wx = (Mx + 1) * My * Mz;
	for (int i = 0; i < dim_Wx; i++)
	{
		//if (i < 4)
		//{
		//	printf("Wx_n[i] = %f \nWx_nplus05[i] = %f \n", Wx_n[i], Wx_nplus05[i]);
		//}
		//if (i < 33)
		//{
		//	printf("Wx_n[%d] = %f Wx_nplus05[%d] = %f \n", i, Wx_n[i], i, Wx_nplus05[i]);
		//}
		Wx_nplus05[i] = Wx_n[i] + sigma*tau*Wx_nplus05[i];
		//if (i < 33)
		//	printf("after Wx_nplus05[%d] = %f \n",i, Wx_nplus05[i]);
		//if (i < 4)
		//{
		//	printf("Wx_nplus05_FINAL[i] = %f \n", Wx_nplus05[i]);
		//}
	}


	return 0;
}

int Implicit_Wy(double * Wy_nplus05, double * Wy_n, double * righthandY, boundaryConditions yCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	double a_0, a_iminus1, a_i, a_nminus1, b_0, b_i, b_n;
	double m_0, m_iminus1, m_i, m_nminus1;

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

			//if (i == 0)
			//{
			//	printf("alfay[%d] = %f \n",0, alfay[0]);
			//	printf("betay[%d] = %f \n",0, betay[0]);
			//	printf("alfay[%d] = %f \n",1, alfay[1]);
			//	printf("betay[%d] = %f \n",1, betay[1]);
			//}
			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)*hz(k1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы mx
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1)*hz(k1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1)*hz(k1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					//right
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}
				else if (k==1)
				{
					alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
					betay[k + 1] = ( righthandY[i * (My + 1) + k] ) / ( b_i );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
					betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}

				//if (i == 0)
				//{
				//	printf("alfay[%d] = %f \n",k+1, alfay[k+1]);
				//	printf("betay[%d] = %f \n",k+1, betay[k+1]);
				//}


			}

			Wy_nplus05[i*(My+1) + My] =  righthandY[i * (My + 1) + My] ;
			//if (i == 0)
			//{
			//	printf("Wy_nplus05[%d] = %f \n",i*(My+1) + My,Wy_nplus05[i*(My+1) + My]);
			//}
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
			a_0_old = hy(0)*hx(i1)*hz(k1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы mx
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1)*hz(k1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_0); 

			betay[0] = 0 ;
			betay[1] = righthandY[i * (My + 1) + 0] / b_0 ;
			alfay[0] = 0 ;
			alfay[1] = -a_0/b_0 ;

//			printf ( "alfay[1] = %f betay[1] = %f \n", alfay[1], betay[1] );
			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)*hz(k1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы mx
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1)*hz(k1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1)*hz(k1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*( 1.0/m_i + 1.0/m_iminus1 );

				alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
				betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);

//				printf ( "alfay[%d] = %f betay[%d] = %f rhand[%d] = %f \n", k + 1, alfay[k+1], k+1, betay[k+1], i * (My + 1) + k, righthandY[i * (My + 1) + k] );
			}

			//к-ты матрицы А
			a_nminus1_old = hy(Ny-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,Ny-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			//к-ты матрицы mx
			m_nminus1 = heat_capacity(i1,Ny-1,k1)*density(i1,Ny-1,k1)*hy(Ny-1)*hx(i1)*hz(k1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_nminus1);

			Wy_nplus05[i*(My+1) +  My] = (   righthandY[i * (My + 1) + My]  - betay[My]*a_nminus1 ) / (alfay[My]*a_nminus1 + b_n);

//			printf ( "Wy_nplus05[last] = %f \n", Wy_nplus05[i*(My+1) +  My] );

//			if (i == 0)
//				printf("delta_Wy_nplus05[%d] = %f \n",i*(My+1) +  My,Wy_nplus05[i*(My+1) +  My]);
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				Wy_nplus05[i*(My+1) + My - j] = Wy_nplus05[i*(My+1) + My + 1 - j] * alfay[My + 1 - j] + betay[My + 1 - j];
//				if (i == 0)
//					printf("delta_Wy_nplus05[%d] = %f \n",i*(My+1) + My - j,Wy_nplus05[i*(My+1) + My - j]);
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
			a_0_old = hy(0)*hx(i1)*hz(k1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы mx
			m_0 = heat_capacity(i1,0,k1)*density(i1,0,k1)*hy(0)*hx(i1)*hz(k1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_0); 

			betay[0] = 0 ;
			betay[1] =  righthandY[i * (My + 1) + 0]/b_0 ;
			alfay[0] = 0 ;
			alfay[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)*hz(k1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)*hz(k1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы mx
				m_i = heat_capacity(i1,k,k1)*density(i1,k,k1)*hy(k)*hx(i1)*hz(k1);
				m_iminus1 = heat_capacity(i1,k-1,k1)*density(i1,k-1,k1)*hy(k-1)*hx(i1)*hz(k1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hz(k1)*hx(i1)*hz(k1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = (   righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
				}
				else
				{
					alfay[k + 1] = (-1)*a_i/(alfay[k]*a_iminus1 + b_i);
					betay[k + 1] = (  righthandY[i * (My + 1) + k]  - betay[k]*a_iminus1 ) / (alfay[k]*a_iminus1 + b_i);
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

	int dim_Wy = Mx * (My + 1) * Mz;
	for (int i = 0; i < dim_Wy; i++)
	{
		//if (i < 33)
		//{
		//	printf("Wy_n[%d] = %f Wy_nplus05[%d] = %f \n",i, Wy_n[i],i,Wy_nplus05[i]);
		//}
		Wy_nplus05[i] = Wy_n[i] + sigma*tau*Wy_nplus05[i];
		//if (i < 33)
		//	printf("after Wy_nplus05[%d] = %f \n", i, Wy_nplus05[i]);
#ifdef ZEROING
		Wy_nplus05[i] = 0;
#endif
	}

	return 0;
}


int Implicit_Wz(double * Wz_nplus05, double * Wz_n, double * righthandZ, boundaryConditions zCondition, double tau, double sigma, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	double a_0, a_iminus1, a_i, a_nminus1, b_0, b_i, b_n;
	double m_0, m_iminus1, m_i, m_nminus1;

	double * alfaz = (double*) malloc ((Mz + 1) * sizeof(double));
	double * betaz = (double*) malloc ((Mz + 1) * sizeof(double));

	switch(zCondition)
	{
	case eNeumann:
		//understep 1 of Step 1 - обращение среднего блока, т.е ~ z-компоненты
		//Нейман, z
		for (int i = 0; i < My*Mx; i++)
		{
			//zxy
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			//zyx
			i1 = i/My;
			j1 = i - i1*My;

			betaz[0] = 0 ;
			betaz[1] = righthandZ[i * (Mz + 1) + 0];
			alfaz[0] = 0 ;
			alfaz[1] = 0 ;
			for ( int k = 1 ; k < Mz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(j1)*hx(i1)*hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hy(j1)*hx(i1)*hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы mx
				m_i = heat_capacity(i1,j1,k)*density(i1,j1,k)*hy(j1)*hx(i1)*hz(k);
				m_iminus1 = heat_capacity(i1,j1,k-1)*density(i1,j1,k-1)*hy(j1)*hx(i1)*hz(k-1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==Mz-1)
				{
					alfaz[k + 1] = 0;
					//right
					betaz[k + 1] = (  righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1 ) / (alfaz[k]*a_iminus1 + b_i);
				}
				else if (k==1)
				{
					alfaz[k + 1] = (-1)*a_i/(alfaz[k]*a_iminus1 + b_i);
					betaz[k + 1] = ( righthandZ[i * (Mz + 1) + k] ) / ( b_i );
				}
				else
				{
					alfaz[k + 1] = (-1)*a_i/(alfaz[k]*a_iminus1 + b_i);
					betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1 ) / (alfaz[k]*a_iminus1 + b_i);
				}

			}

			Wz_nplus05[i*(Mz+1) + Mz] =  righthandZ[i * (Mz + 1) + Mz] ;
			for ( int j = 1 ; j < Mz + 1 ; j++ )
			{
				Wz_nplus05[i*(Mz+1) + Mz - j] = Wz_nplus05[i*(Mz+1) + Mz + 1 - j] * alfaz[Mz + 1 - j] + betaz[Mz + 1 - j];
			}
		}
		break;
	case eDirichlet:
		//Дирихле, z
		for (int i = 0; i < My*Mx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/My;
			j1 = i - i1*My;

			//к-ты матрицы А
			a_0_old = hy(j1)*hx(i1)*hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы mx
			m_0 = heat_capacity(i1,j1,0)*density(i1,j1,0)*hy(j1)*hx(i1)*hz(0);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_0); 

			betaz[0] = 0 ;
			betaz[1] = righthandZ[i * (Mz + 1) + 0] / b_0 ;
			alfaz[0] = 0 ;
			alfaz[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < Mz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(j1)*hx(i1)*hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hy(j1)*hx(i1)*hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы mx
				m_i = heat_capacity(i1,j1,k)*density(i1,j1,k)*hy(j1)*hx(i1)*hz(k);
				m_iminus1 = heat_capacity(i1,j1,k-1)*density(i1,j1,k-1)*hy(j1)*hx(i1)*hz(k-1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				alfaz[k + 1] = (-1)*a_i/(alfaz[k]*a_iminus1 + b_i);
				betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1 ) / (alfaz[k]*a_iminus1 + b_i);
			}

			//к-ты матрицы А
			a_nminus1_old = hy(j1)*hx(i1)*hz(Nz-1)/(heat_conductivity(i1,j1,Nz-1)*6.0);
			b_n_old = 2*a_nminus1_old;
			//к-ты матрицы mx
			m_nminus1 = heat_capacity(i1,j1,Nz-1)*density(i1,j1,Nz-1)*hy(j1)*hx(i1)*hz(Nz-1);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_nminus1 = a_nminus1_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_nminus1);
			b_n = b_n_old + 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_nminus1);


			Wz_nplus05[i*(Mz+1) +  Mz] = (   righthandZ[i * (Mz + 1) + Mz]  - betaz[Mz]*a_nminus1 ) / (alfaz[Mz]*a_nminus1 + b_n);

			for ( int j = 1 ; j < Mz + 1 ; j++ )
			{
				Wz_nplus05[i*(Mz+1) + Mz - j] = Wz_nplus05[i*(Mz+1) + Mz + 1 - j] * alfaz[Mz + 1 - j] + betaz[Mz + 1 - j];
//				if (i == 0)
//					printf("delta_Wy_nplus05[%d] = %f \n",i*(My+1) + My - j,Wy_nplus05[i*(My+1) + My - j]);
			}
		}
		break;
	case eMixed:
		//Дирихле-Нейман, z
		{;}

		for (int i = 0; i < My*Mx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/My;
			j1 = i - i1*My;

			//к-ты матрицы А
			a_0_old = hy(j1)*hx(i1)*hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;
			//к-ты матрицы mx
			m_0 = heat_capacity(i1,j1,0)*density(i1,j1,0)*hy(j1)*hx(i1)*hz(0);
			//к-ты матрицы A + tau/2 * BM(-1)B(tr)
			a_0 = a_0_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_0);
			b_0 = b_0_old + 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_0); 

			betaz[0] = 0 ;
			betaz[1] =  righthandZ[i * (Mz + 1) + 0]/b_0 ;
			alfaz[0] = 0 ;
			alfaz[1] = -a_0/b_0 ;

			for ( int k = 1 ; k < Mz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(j1)*hx(i1)*hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hy(j1)*hx(i1)*hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				//к-ты матрицы mx
				m_i = heat_capacity(i1,j1,k)*density(i1,j1,k)*hy(j1)*hx(i1)*hz(k);
				m_iminus1 = heat_capacity(i1,j1,k-1)*density(i1,j1,k-1)*hy(j1)*hx(i1)*hz(k-1);

				//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
				a_i = a_i_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_i);
				a_iminus1 = a_iminus1_old - 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*(1.0/m_iminus1);
				b_i = b_i_old + 0.5*tau*hx(i1)*hy(j1)*hx(i1)*hy(j1)*( 1.0/m_i + 1.0/m_iminus1 );

				if (k==Mz-1)
				{
					alfaz[k + 1] = 0;
					betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1 ) / (alfaz[k]*a_iminus1 + b_i);
				}
				else
				{
					alfaz[k + 1] = (-1)*a_i/(alfaz[k]*a_iminus1 + b_i);
					betaz[k + 1] = (   righthandZ[i * (Mz + 1) + k]  - betaz[k]*a_iminus1 ) / (alfaz[k]*a_iminus1 + b_i);
				}
			}
			Wz_nplus05[i*(Mz+1) + Mz] = righthandZ[i*(Mz + 1) + Mz] ;

			for ( int j = 1 ; j < Mz + 1 ; j++ )
			{
				Wz_nplus05[i*(Mz+1) + Mz - j] = Wz_nplus05[i*(Mz+1) + Mz + 1 - j] * alfaz[Mz + 1 - j] + betaz[Mz + 1 - j];
			}
		}
		break;
	}
	{;}

	free(alfaz);
	free(betaz);

	int dim_Wz = Mx * My * (Mz + 1);
	for (int i = 0; i < dim_Wz; i++)
	{
		Wz_nplus05[i] = Wz_n[i] + sigma*tau*Wz_nplus05[i];
#ifdef ZEROING
		Wz_nplus05[i] = 0;
#endif
	}

	return 0;
}


int Temperature_standart_step(double * T_n, double * T_nplus1, double * Vx, double coeff_Vx, double * Wx_n, double * Wx_nplus1, double * Vy, double coeff_Vy, double * Wy_n, double * Wy_nplus1, double * Vz, double coeff_Vz, double * Wz_n, double * Wz_nplus1, double sigma, double tau, int numSolution, int t, int Mx, int My, int Mz)
{
	int i1, j1, k1;
	double m_i, a_i_old, temp, koeff_k, conv_k, conv_kplus1, CONVECT_x, CONVECT_y, CONVECT_z;
	for (int i = 0; i < My*Mz; i++)
	{
		j1 = i - (i/My)*My;
		k1 = i/My;
		for ( int k = 0 ; k < Mx ; k++ )
		{
			i1 = k;
			m_i = hx(i1)*hy(j1)*hz(k1)*density(i1,j1,k1)*heat_capacity(i1,j1,k1);

			temp = 0.5*(righthand(numSolution,(t+1),i1,j1,k1,xpoints,ypoints,zpoints) + righthand(numSolution,t,i1,j1,k1,xpoints,ypoints,zpoints));

			a_i_old = 1.0/(heat_conductivity(i1,j1,k1)*6.0);
			koeff_k = a_i_old;
			conv_k = koeff_k*( 2*Vx[i*(Mx+1)+ i1] + Vx[i*(Mx+1)+ i1+1])*Wx_n[i*(Mx+1) + i1];
			conv_kplus1 = koeff_k*( Vx[i*(Mx+1)+ i1] + 2*Vx[i*(Mx+1)+ i1+1])*Wx_n[i*(Mx+1) + i1 + 1];
			CONVECT_x = conv_k + conv_kplus1;

			//правильное
			//T_nplus1[i*Mx + k] = T[i*Mx + k] + tau * (1.0/m_i) * temp - tau * (1.0/m_i)*hy(j1)*( 0.5*(Wx_n[i*(Mx+1) + k + 1] - Wx_n[i*(Mx+1) + k] ) + 0.5*(Wx_nplus1[i*(Mx+1) + k + 1] - Wx_nplus1[i*(Mx+1) + k] )) ;
			T_nplus1[i*Mx + k] = T_n[i*Mx + k] + tau * (1.0/m_i) * temp + tau * coeff_Vx * CONVECT_x - tau * (1.0/m_i)*hy(j1)*hz(k1)*( sigma*(Wx_n[i*(Mx+1) + k + 1] - Wx_n[i*(Mx+1) + k] ) + (1 - sigma)*(Wx_nplus1[i*(Mx+1) + k + 1] - Wx_nplus1[i*(Mx+1) + k] )) ;
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
			m_i = hx(i1)*hy(j)*hz(k1)*density(i1,j1,k1)*heat_capacity(i1,j1,k1);

			a_i_old = 1.0/(heat_conductivity(i1,j1,k1)*6.0);
			koeff_k = a_i_old;
			conv_k = koeff_k*(2*Vy[i*(My+1)+ j1] + Vy[i*(My+1) + j1+1])*Wy_n[i*(My+1) + j1];
			conv_kplus1 = koeff_k*(Vy[i*(My+1)+ j1] + 2*Vy[i*(My+1) + j1+1])*Wy_n[i*(My+1) + j1 + 1];
			CONVECT_y = conv_k + conv_kplus1;

			//T_nplus1[i1 + j1*Mx + k1*Mx*My] +=  - tau * (1.0/m_i)* hx(i1)*( 0.5*(Wy_n[i*(My + 1) + 1 + j] - Wy_n[i*(My + 1) + j] ) + 0.5*(Wy_nplus1[i*(My + 1) + 1 + j] - Wy_nplus1[i*(My + 1) + j] )) ;
			T_nplus1[i1 + j1*Mx + k1*Mx*My] +=  tau * coeff_Vy * CONVECT_y - tau * (1.0/m_i)* hx(i1)*hz(k1)*( sigma*(Wy_n[i*(My + 1) + 1 + j] - Wy_n[i*(My + 1) + j] ) + (1 - sigma)*(Wy_nplus1[i*(My + 1) + 1 + j] - Wy_nplus1[i*(My + 1) + j] )) ;
		}
	}

	for (int i = 0; i < Ny*Nx; i++)
	{
		//xyz
		//j1 = i/Nx;
		//i1 = i - (i/Nx)*Nx;
		//zyx
		i1 = i/Ny;
		j1 = i - i1*Ny;

		for ( int k = 0 ; k < Nz ; k++ )
		{
			k1 = k;
			m_i = hx(i1)*hy(j1)*hz(k)*density(i1,j1,k)*heat_capacity(i1,j1,k);

			a_i_old = 1.0/(heat_conductivity(i1,j1,k1)*6.0);
			koeff_k = a_i_old;
			conv_k = koeff_k*(2*Vz[i*(Mz+1)+ k1] + Vz[i*(Mz+1) + k1+1])*Wz_n[i*(Mz+1) + k1];
			conv_kplus1 = koeff_k*(Vz[i*(Mz+1)+ k1] + 2*Vz[i*(Mz+1) + k1+1])*Wz_n[i*(Mz+1) + k1 + 1];
			CONVECT_z = conv_k + conv_kplus1;

			T_nplus1[i1 + j1*Mx + k1*Mx*My] +=  tau * coeff_Vz * CONVECT_z - tau * (1.0/m_i)* hx(i1)*hy(j1)*( sigma*(Wz_n[i*(Mz + 1) + 1 + k] - Wz_n[i*(Mz + 1) + k] ) + (1 - sigma)*(Wz_nplus1[i*(Mz + 1) + 1 + k] - Wz_nplus1[i*(Mz + 1) + k] )) ;
		}


	}

	return 0.0;
}

int Temperature_Flux_Update(double * Wx_n, double * Wx_nplus1, double * Wy_n, double * Wy_nplus1, double * Wz_n, double * Wz_nplus1, double * T_n, double * T_nplus1, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;

	for (int i = 0; i < dim_Wx; i++)
		Wx_n[i] = Wx_nplus1[i];
	for (int i = 0; i < dim_Wy; i++)
		Wy_n[i] = Wy_nplus1[i];
	for (int i = 0; i < dim_Wz; i++)
		Wz_n[i] = Wz_nplus1[i];
	for (int i = 0; i < dim_T; i++)
		T_n[i] = T_nplus1[i];

	return 0;
}

int DouglasGunn_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, 
								double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus1,
								double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus13, double * Wy_nplus1, 
								double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus13, double * Wz_nplus1, 
								int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step,
								double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt,
								double * eps_relative_max_w_pt, double * eps_relative_max_T_pt,	double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt,
								double * diffPr_max_w_pt, double * diffPr_l2_w_pt, double * diffrelPr_max_w_pt, double * diffrelPr_l2_w_pt,
								double * diffPr_max_T_pt, double * diffPr_l2_T_pt, double * diffrelPr_max_T_pt, double * diffrelPr_l2_T_pt )
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;

	//computing new FLUX W_nplus1
	//printf("starting \n");

	//steps
	//with velocity
	Righthand_step1_wynplus13(righthandY, F, 1.0, 0.0, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, numSolution, t, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus13_3d.txt", righthandY, dim_Wy );
//right
	Explicit_Wy(Wy_nplus13, Wy_n, righthandY, yCondition, tau, 1.0, Mx, My, Mz);
//	Explicit_Wy(Wy_nplus13, Wy_n, righthandY, yCondition, tau, 0.5, Mx, My, Mz);

//	fprintdVec ( "Wy_nplus13_3d.txt", Wy_nplus13, dim_Wy);

	Righthand_step1_wznplus13(righthandZ, F, 1.0, 0.0, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, numSolution, t, Mx, My, Mz);
//	fprintdVec ( "righthandZ.txt", righthandZ, dim_Wz);
	Explicit_Wz(Wz_nplus13, Wz_n, righthandZ, zCondition, tau, 1.0, Mx, My, Mz);

//	fprintdVec ( "Wz_nplus13_3d.txt", Wz_nplus13, dim_Wz);

	for (int i = 0; i < dim_Wy; i++)
		Temp_Wy[i] = 0.5 * (Wy_nplus13[i] + Wy_n[i]);
	for (int i = 0; i < dim_Wz; i++)
		Temp_Wz[i] = 0.5 * (Wz_nplus13[i] + Wz_n[i]);

	Righthand_step1_wxnplus1(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_n, Temp_Wy, Temp_Wz, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus1, Wx_n, righthandX, xCondition, tau, 1.0, Mx, My, Mz);

//	fprintdVec ( "Wx_nplus1_3d.txt", Wx_nplus1, dim_Wx);

	for (int i = 0; i < dim_Wx; i++)
		Temp_Wx[i] = 0.5 * (Wx_nplus1[i] - Wx_n[i]);
	for (int i = 0; i < dim_Wy; i++)
//right
		Temp_Wy[i] = 0.5 * (Wy_nplus13[i] - Wy_n[i]);
//		Temp_Wy[i] = Wy_nplus13[i];
	for (int i = 0; i < dim_Wz; i++)
		Temp_Wz[i] = 0.5 * (Wz_nplus13[i] - Wz_n[i]);

//right
	Righthand_step1_wynplus1(righthandY, F, -0.5, 0.5, Vx, 0.0, Wx_n, -1.0, xCondition, Vy, 0.0, Wy_n, -1.0, yCondition, Vz, 0.0, Wz_n, -1.0, zCondition, Temp_Wx, Temp_Wy, Temp_Wz, numSolution, t, Mx, My, Mz);
	Implicit_Wy(Wy_nplus1, Wy_nplus13, righthandY, yCondition, tau, 1.0, Mx, My, Mz);
//	Righthand_step1_wynplus1(righthandY, F, 0.0, 1.0, Vx, 0.0, Wx_n, -1.0, xCondition, Vy, 0.0, Wy_n, -1.0, yCondition, Vz, 0.0, Wz_n, -1.0, zCondition, Temp_Wx, Temp_Wy, Temp_Wz, numSolution, t, Mx, My, Mz);
//	fprintdVec ( "rhand_wy_nplus1_3d.txt", righthandY, dim_Wy );
//	Implicit_Wy(Wy_nplus1, Wy_nplus13, righthandY, yCondition, tau, 0.5, Mx, My, Mz);

//	fprintdVec ( "Wy_nplus1_3d.txt", Wy_nplus1, dim_Wy);

	for (int i = 0; i < dim_Wy; i++)
		Temp_Wy[i] = 0.5 * (Wy_nplus1[i] - Wy_n[i]);

	Righthand_step1_wznplus1(righthandZ, F, -0.5, 0.5, Vx, 0.0, Wx_n, -1.0, xCondition, Vy, 0.0, Wy_n, -1.0, yCondition, Vz, 0.0, Wz_n, -1.0, zCondition, Temp_Wx, Temp_Wy, Temp_Wz, numSolution, t, Mx, My, Mz);
	Implicit_Wz(Wz_nplus1, Wz_nplus13, righthandZ, zCondition, tau, 1.0, Mx, My, Mz);	//computing new TEMPERATURE T_nplus1

//	fprintdVec ( "Wz_nplus1_3d.txt", Wz_nplus1, dim_Wz);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, 1.0, Wx_n, Wx_nplus1, Vy, 1.0, Wy_n, Wy_nplus1, Vz, 1.0, Wz_n, Wz_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

//	fprintdVec ( "T_nplus1_3d.txt", T_nplus1, dim_T);

//	fprintdVec ( "Wy_nplus1_again_3d.txt", Wy_nplus1, dim_Wy);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, Wz_n, Wz_nplus1, T_n, T_nplus1, Mx, My, Mz);

#ifndef ERROR_FOR_L2PROJECTION
	Accuracy_calculate(T_n, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);
#else
	Accuracy_calculateL2Pr(T_n, Wx_n, Wy_n, Wz_n, T_nplus1, Temp_Wx, Temp_Wy, Temp_Wz, xCondition, yCondition, zCondition,
		numSolution, t, Mx, My, Mz, print_step,
		eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, 
		eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);
#endif

#ifdef COMPARE_PROJECTORS
	// T_nplus1 is used as a temporary array for exact averaged T_n
	DifferenceProjectors_calc(T_nplus1, Temp_Wx, Temp_Wy, Temp_Wz, 
		(t+1) * tau, Mx, My, Mz, numSolution, xCondition, yCondition, zCondition,
		diffPr_max_w_pt, diffPr_l2_w_pt, diffrelPr_max_w_pt, diffrelPr_l2_w_pt,
		diffPr_max_T_pt, diffPr_l2_T_pt, diffrelPr_max_T_pt, diffrelPr_l2_T_pt);
#endif
	return 0.0;
}


// something strange:
// 1) seems to be unstable for large Courant numbers although temperature must be stable
// 2) in 2D-case is long way from what POredictorCorrector shows in 2D, although must be the same.
// But 2) is a problem also for DouglasGunn. Solved for DouglasGunn but not for PredictorCorrector.
// 2) Solved also for PredictorCorrector, a bug was fixed.
// 1) now works, not good for temperature but converged with 2nd order for 0.1-32 and 0.05-64.

int PredictorCorrector3D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, 
										 double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus12, double * Wx_nplus1, 
										 double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus12, double * Wy_nplus1, 
										 double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus12, double * Wz_nplus1, 
										 int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step,
										 double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, 
										 double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt,
										 double * diffPr_max_w_pt, double * diffPr_l2_w_pt, double * diffrelPr_max_w_pt, double * diffrelPr_l2_w_pt,
										 double * diffPr_max_T_pt, double * diffPr_l2_T_pt, double * diffrelPr_max_T_pt, double * diffrelPr_l2_T_pt )
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;
	//...write it down first

	//steps
	//--------------------- predictor
	Righthand_step1_wynplus1(righthandY, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, 
		Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz);

	//fprintdVec ( "Righthand_Wy_nplus12_3d.txt", righthandY, dim_Wy);

	Implicit_Wy(Wy_nplus12, Wy_n, righthandY, yCondition, tau, 0.5, Mx, My, Mz);

	//fprintdVec ( "Wy_nplus12_3d.txt", Wy_nplus12, dim_Wy);

	//printf("WTF WTF Wy_nplus12[%d] = %f \n",1, Wy_nplus12[1]);

	//FILE * file1 = fopen("righthand1_pred.xls","wt");
	//for (int i = 0; i < dim_Wy; i++)
	//	fprintf(file1, "%f \n", righthandY[i]);
	//fclose(file1);
	//FILE * file2 = fopen("wy_12222_pred.xls","wt");
	//for (int i = 0; i < dim_Wy; i++)
	//{
	//	//if (i < 33)
	//	//	printf("WTF Wy_nplus12[%d] = %f \n",i, Wy_nplus12[i]);
	//	fprintf(file2, "%f \n", Wy_nplus12[i]);
	//}
	//fclose(file2);

	Righthand_step1_wznplus1(righthandZ, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, 
		Wx_n, Wy_nplus12, Wz_n, numSolution, t, Mx, My, Mz);

	//fprintdVec ( "Righthand_Wz_nplus12_3d.txt", righthandZ, dim_Wz);

	Implicit_Wz(Wz_nplus12, Wz_n, righthandZ, zCondition, tau, 0.5, Mx, My, Mz);

	//fprintdVec ( "Wz_nplus12_3d.txt", Wz_nplus12, dim_Wz);

	//FILE * file3 = fopen("righthand2_pred.xls","wt");
	//for (int i = 0; i < dim_Wz; i++)
	//	fprintf(file3, "%f \n", righthandZ[i]);
	//fclose(file3);
	//FILE * file4 = fopen("wz_12_pred.xls","wt");
	//for (int i = 0; i < dim_Wz; i++)
	//	fprintf(file4, "%f \n", Wz_nplus12[i]);
	//fclose(file4);

	Righthand_step1_wxnplus1(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_n, Wy_nplus12, Wz_nplus12, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus12, Wx_n, righthandX, xCondition, tau, 0.5, Mx, My, Mz);

	//fprintdVec ( "Wx_nplus12_3d.txt", Wx_nplus12, dim_Wx);

	//FILE * file5 = fopen("righthand3_pred.xls","wt");
	//for (int i = 0; i < dim_Wx; i++)
	//	fprintf(file5, "%f \n", righthandX[i]);
	//fclose(file5);
	//FILE * file6 = fopen("wx_1222222_pred.xls","wt");
	//for (int i = 0; i < dim_Wx; i++)
	//	fprintf(file6, "%f \n", Wx_nplus12[i]);
	//fclose(file6);

	//--------------------- corrector

	Righthand_step1_wxnplus1(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_nplus12, Wy_nplus12, Wz_nplus12, numSolution, t, Mx, My, Mz);
	Explicit_Wx(Wx_nplus1, Wx_n, righthandX, xCondition, tau, 1.0, Mx, My, Mz);

	//fprintdVec ( "Wx_nplus1_3d.txt", Wx_nplus1, dim_Wx);

	//FILE * file7 = fopen("righthand4_pred.xls","wt");
	//for (int i = 0; i < dim_Wx; i++)
	//	fprintf(file7, "%f \n", righthandX[i]);
	//fclose(file7);
	//FILE * file8 = fopen("wx_1.xls","wt");
	//for (int i = 0; i < dim_Wx; i++)
	//	fprintf(file8, "%f \n", Wx_nplus1[i]);
	//fclose(file8);

	Righthand_step1_wynplus1(righthandY, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_nplus12, Wy_nplus12, Wz_nplus12, numSolution, t, Mx, My, Mz);
	Explicit_Wy(Wy_nplus1, Wy_n, righthandY, yCondition, tau, 1.0, Mx, My, Mz);

	//fprintdVec ( "Wy_nplus1_3d.txt", Wy_nplus1, dim_Wy);

	//FILE * file9 = fopen("righthand5_pred.xls","wt");
	//for (int i = 0; i < dim_Wy; i++)
	//	fprintf(file9, "%f \n", righthandY[i]);
	//fclose(file9);
	//FILE * file10 = fopen("wy_12_pred.xls","wt");
	//for (int i = 0; i < dim_Wy; i++)
	//	fprintf(file10, "%f \n", Wy_nplus1[i]);
	//fclose(file10);


	Righthand_step1_wznplus1(righthandZ, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_nplus12, Wy_nplus12, Wz_nplus12, numSolution, t, Mx, My, Mz);
	Explicit_Wz(Wz_nplus1, Wz_n, righthandZ, zCondition, tau, 1.0, Mx, My, Mz);

	//fprintdVec ( "Wz_nplus1_3d.txt", Wz_nplus1, dim_Wz);

	//FILE * file11 = fopen("righthand6_pred.xls","wt");
	//for (int i = 0; i < dim_Wz; i++)
	//	fprintf(file11, "%f \n", righthandZ[i]);
	//fclose(file11);
	//FILE * file12 = fopen("wz_1_pred.xls","wt");
	//for (int i = 0; i < dim_Wz; i++)
	//	fprintf(file12, "%f \n", Wz_nplus1[i]);
	//fclose(file12);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, 1.0, Wx_n, Wx_nplus1, Vy, 1.0, Wy_n, Wy_nplus1, Vz, 1.0, Wz_n, Wz_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

	//fprintdVec ( "T_nplus1_3d.txt", T_nplus1, dim_T);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, Wz_n, Wz_nplus1, T_n, T_nplus1, Mx, My, Mz);

//	Accuracy_calculate(T_n, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);
#ifndef ERROR_FOR_L2PROJECTION
	Accuracy_calculate(T_n, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);
#else
	Accuracy_calculateL2Pr(T_n, Wx_n, Wy_n, Wz_n, T_nplus1, Temp_Wx, Temp_Wy, Temp_Wz, xCondition, yCondition, zCondition,
		numSolution, t, Mx, My, Mz, print_step,
		eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, 
		eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);
#endif

#ifdef COMPARE_PROJECTORS
	// T_nplus1 is used as a temporary array for exact averaged T_n
	DifferenceProjectors_calc(T_nplus1, Temp_Wx, Temp_Wy, Temp_Wz, 
		(t+1) * tau, Mx, My, Mz, numSolution, xCondition, yCondition, zCondition,
		diffPr_max_w_pt, diffPr_l2_w_pt, diffrelPr_max_w_pt, diffrelPr_l2_w_pt,
		diffPr_max_T_pt, diffPr_l2_T_pt, diffrelPr_max_T_pt, diffrelPr_l2_T_pt);
#endif

	return 0.0;
}

int Uzawa3D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * Temp_Wx, double * Temp_Wy, double * Temp_Wz, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus12, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus12, double * Wy_nplus1, double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus1, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;
	//...write it down first

	//steps
	//--------------------- predictor
	Righthand_step1_wynplus1(righthandY, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz);
	Implicit_Wy(Wy_nplus12, Wy_n, righthandY, yCondition, tau, 1.0, Mx, My, Mz);

	//printf("WTF WTF Wy_nplus12[%d] = %f \n",1, Wy_nplus12[1]);

	//FILE * file1 = fopen("righthand1_Uz.xls","wt");
	//for (int i = 0; i < dim_Wy; i++)
	//	fprintf(file1, "%f \n", righthandY[i]);
	//fclose(file1);
	//FILE * file2 = fopen("wy_12222_Uz.xls","wt");
	//for (int i = 0; i < dim_Wy; i++)
	//{
	//	fprintf(file2, "%f \n", Wy_nplus12[i]);
	//}
	//fclose(file2);

	for (int i = 0; i < dim_Wy; i++)
		Temp_Wy[i] = 0.5 * (Wy_n[i] + Wy_nplus12[i]);

	Righthand_step1_wxnplus1(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_n, Temp_Wy, Wz_n, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus12, Wx_n, righthandX, xCondition, tau, 1.0, Mx, My, Mz);

	//FILE * file5 = fopen("righthand3_Uz.xls","wt");
	//for (int i = 0; i < dim_Wx; i++)
	//	fprintf(file5, "%f \n", righthandX[i]);
	//fclose(file5);
	//FILE * file6 = fopen("wx_1222222_Uz.xls","wt");
	//for (int i = 0; i < dim_Wx; i++)
	//	fprintf(file6, "%f \n", Wx_nplus12[i]);
	//fclose(file6);

	for (int i = 0; i < dim_Wx; i++)
		Temp_Wx[i] = 0.5 * (Wx_n[i] + Wx_nplus12[i]);

	Righthand_step1_wznplus1(righthandZ, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Temp_Wx, Temp_Wy, Wz_n, numSolution, t, Mx, My, Mz);
	Implicit_Wz(Wz_nplus1, Wz_n, righthandZ, zCondition, tau, 1.0, Mx, My, Mz);

	//FILE * file3 = fopen("righthand2_Uz.xls","wt");
	//for (int i = 0; i < dim_Wz; i++)
	//	fprintf(file3, "%f \n", righthandZ[i]);
	//fclose(file3);
	//FILE * file4 = fopen("wz_1_Uz.xls","wt");
	//for (int i = 0; i < dim_Wz; i++)
	//	fprintf(file4, "%f \n", Wz_nplus1[i]);
	//fclose(file4);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, 1.0, Wx_n, Wx_nplus12, Vy, 1.0, Wy_n, Wy_nplus12, Vz, 1.0, Wz_n, Wz_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

	HeatFlux_Wx0_Init(T_nplus1, Wx_nplus1, xCondition, numSolution, Mx, My, Mz);
	HeatFlux_Wy0_Init(T_nplus1, Wy_nplus1, yCondition, numSolution, Mx, My, Mz);

	//for (int i = 0; i < dim_Wz; i++)
	//	Temp_Wz[i] = 0.5 * (Wz_n[i] + Wz_nplus1[i]);

	//Righthand_step1_wxnplus1(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Temp_Wx, Temp_Wy, Temp_Wz, numSolution, t, Mx, My, Mz);
	//Explicit_Wx(Wx_nplus1, Wx_n, righthandX, xCondition, tau, 1.0, Mx, My, Mz);


	////FILE * file7 = fopen("righthand4.xls","wt");
	////for (int i = 0; i < dim_Wx; i++)
	////	fprintf(file7, "%f \n", righthandX[i]);
	////fclose(file7);
	////FILE * file8 = fopen("wx_1.xls","wt");
	////for (int i = 0; i < dim_Wx; i++)
	////	fprintf(file8, "%f \n", Wx_nplus1[i]);
	////fclose(file8);

	//Righthand_step1_wynplus1(righthandY, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Temp_Wx, Temp_Wy, Temp_Wz, numSolution, t, Mx, My, Mz);
	//Explicit_Wy(Wy_nplus1, Wy_n, righthandY, yCondition, tau, 1.0, Mx, My, Mz);

	////FILE * file9 = fopen("righthand5.xls","wt");
	////for (int i = 0; i < dim_Wy; i++)
	////	fprintf(file9, "%f \n", righthandY[i]);
	////fclose(file9);
	////FILE * file10 = fopen("wy_12.xls","wt");
	////for (int i = 0; i < dim_Wy; i++)
	////	fprintf(file10, "%f \n", Wy_nplus1[i]);
	////fclose(file10);


	////FILE * file11 = fopen("righthand6.xls","wt");
	////for (int i = 0; i < dim_Wz; i++)
	////	fprintf(file11, "%f \n", righthandZ[i]);
	////fclose(file11);
	////FILE * file12 = fopen("wz_1.xls","wt");
	////for (int i = 0; i < dim_Wz; i++)
	////	fprintf(file12, "%f \n", Wz_nplus1[i]);
	////fclose(file12);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, Wz_n, Wz_nplus1, T_n, T_nplus1, Mx, My, Mz);

	Accuracy_calculate(T_n, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);

	return 0.0;
}
int Local1D_3D_Scheme_mainstep(double * T_n, double * T_nplus1, double * F, double * righthandX, double * Vx, double * Wx_n, double * Wx_nplus1, double * righthandY, double * Vy, double * Wy_n, double * Wy_nplus1, double * righthandZ, double * Vz, double * Wz_n, double * Wz_nplus1, int t, double tau, int numSolution, boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition, int Mx, int My, int Mz, int print_step, double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt, double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);
	int dim_T = Mx * My * Mz;
	//...write it down first

	//Righthand_step1_wxnplus1_local(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, numSolution, t, Mx, My, Mz);
	Righthand_step1_wxnplus1(righthandX, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz);
	Implicit_Wx(Wx_nplus1, Wx_n, righthandX, xCondition, tau, 1.0, Mx, My, Mz);

	//steps
	//Righthand_step1_wynplus1_local(righthandY, F, 0.5, 0.5, Vx, 1.0, Wx_n, Wx_nplus1, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, numSolution, t, Mx, My, Mz);
	Righthand_step1_wynplus1(righthandY, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_nplus1, Wy_n, Wz_n, numSolution, t, Mx, My, Mz);
	Implicit_Wy(Wy_nplus1, Wy_n, righthandY, yCondition, tau, 1.0, Mx, My, Mz);

	//Righthand_step1_wznplus1_local(righthandZ, F, 0.5, 0.5, Vx, 1.0, Wx_n, Wx_nplus1, -1.0, xCondition, Vy, 1.0, Wy_n, Wy_nplus1, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, numSolution, t, Mx, My, Mz);
	Righthand_step1_wznplus1(righthandZ, F, 0.5, 0.5, Vx, 1.0, Wx_n, -1.0, xCondition, Vy, 1.0, Wy_n, -1.0, yCondition, Vz, 1.0, Wz_n, -1.0, zCondition, Wx_nplus1, Wy_nplus1, Wz_n, numSolution, t, Mx, My, Mz);
	Implicit_Wz(Wz_nplus1, Wz_n, righthandZ, zCondition, tau, 1.0, Mx, My, Mz);

	//computing new TEMPERATURE T_nplus1
	Temperature_standart_step(T_n, T_nplus1, Vx, 1.0, Wx_n, Wx_nplus1, Vy, 1.0, Wy_n, Wy_nplus1, Vz, 1.0, Wz_n, Wz_nplus1, 0.5, tau, numSolution, t, Mx, My, Mz);

	Temperature_Flux_Update(Wx_n, Wx_nplus1, Wy_n, Wy_nplus1, Wz_n, Wz_nplus1, T_n, T_nplus1, Mx, My, Mz);

	Accuracy_calculate(T_n, Wx_n, Wy_n, Wz_n, numSolution, t, Mx, My, Mz, print_step, eps_max_w_pt, eps_max_T_pt, eps_l2_w_pt, eps_l2_T_pt, eps_relative_max_w_pt, eps_relative_max_T_pt, eps_relative_l2_w_pt, eps_relative_l2_T_pt);

	return 0.0;
}


int Velocity_Init(double *Vx, double *Vy, double *Vz, int Mx, int My, int Mz)
{
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);
	int i1, j1, k1;
	// инициализация Vx, Vy, Vz , причем вектор V=(Vx,Vy,Vz) таков, что V*n(скалярное произведение на нормаль) на границе области = 0
	for ( int i = 0 ; i < dim_Wx ; i++ )
	{
		i1 = i%(Mx+1);
		j1 = ((i - i1)/(Mx+1))%My;
		k1 = (i - i1 - j1*(Mx+1))/((Mx+1)*My);
		//		Vx[i] = v_init(i1*hx, (j1+0.5)*hy, (k1+0.5)*hz);

		Vx[i] = exact_solution1111Vx(xpoints[i1],ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),0,mx, my, mz);
	}
	for ( int i = 0 ; i < dim_Wy ; i++ )
	{
		j1 = i%(My+1);
		k1 = ((i - j1)/(My+1))%Mz;
		i1 = (i - j1 - k1*(My+1))/((My+1)*Mz);

		Vy[i] = exact_solution1111Vy(xpoints[i1] + 0.5*hx(i1),ypoints[j1],zpoints[k1] + 0.5*hz(k1),0,mx, my, mz);
	}

	for ( int i = 0 ; i < dim_Wz ; i++ )
	{
		//zyx
		k1 = i%(Mz+1);
		j1 = ((i - k1)/(Mz+1))%My;
		i1 = (i - k1 - j1*(Mz+1))/((Mz+1)*My);

		Vz[i] = exact_solution1111Vz(xpoints[i1] + 0.5*hx(i1),ypoints[j1] + 0.5 * hy(j1), zpoints[k1],0,mx, my, mz);
	}

	return 0;
}


int fprintdVec ( char * filename, double * vec, int dim )
{
	FILE * f = fopen ( filename, "wt" );
	if ( f == NULL )
	{
		printf ( "Cannot open file %s for wiritng \n", filename );
		return -1;
	}
	else
	{
		for ( int i = 0; i < dim; i++ )
		{
			//printf ( "val[%d] = %f \n", i, vec[i] );
			fprintf (f, "%f \n", vec[i]);
		}

		fclose(f);
	}

	return 0;
}

int	DifferenceProjectors_calc( double * T_av, double * Wx_av, double * Wy_av, double * Wz_av,
							  double time, int Mx, int My, int Mz, int numSolution,
							  boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition,
							  double * diffPr_max_w_pt, double * diffPr_l2_w_pt, double * diffrelPr_max_w_pt, double * diffrelPr_l2_w_pt,
							  double * diffPr_max_T_pt, double * diffPr_l2_T_pt, double * diffrelPr_max_T_pt, double * diffrelPr_l2_T_pt)
{
	double diff_max_w = *diffPr_max_w_pt;
	double diff_l2_w = *diffPr_l2_w_pt;
	double diffrel_max_w = *diffrelPr_max_w_pt;
	double diffrel_l2_w = *diffrelPr_l2_w_pt;
	double norm_l2_w = 0;
	double norm_max_w = 0;

	double diff_max_T = *diffPr_max_T_pt;
	double diff_l2_T = *diffPr_l2_T_pt;
	double diffrel_max_T = *diffrelPr_max_T_pt;
	double diffrel_l2_T = *diffrelPr_l2_T_pt;
	double norm_l2_T = 0;
	double norm_max_T = 0;


	int dim_T = Mx * My * Mz;
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);


	for ( int i = 0 ; i < Mx ; i++ )
	{
		for (int j = 0; j < My; j++)
		{
			for (int k = 0; k < Mz; k++)
			{
				T_av[i + j*Nx + k*Nx*Ny] = exact_solution_averageQuad ( numSolution, time, i, j, k, xpoints, ypoints, zpoints, NQUADPOINTS );
			}
		}
	}
	HeatFlux_Init(T_av, Wx_av, xCondition, Wy_av, yCondition, Wz_av, zCondition, numSolution, Mx, My, Mz);

	double temp, temp2;
	int i1, j1, k1;

	double diff_max_T_tn = 0.0;
	double diff_l2_T_tn = 0.0;
	double norm_max_T_tn = 0.0;
	double norm_l2_T_tn = 0.0;
	double diffrel_max_T_tn = 0.0;
	double diffrel_l2_T_tn = 0.0;

	for (int ind = 0; ind < My*Mz; ind++)
	{
		j1 = ind - (ind/My)*My;
		k1 = ind/My;
		for ( int i = 0 ; i < Mx ; i++ )
		{
			temp = 0;
			temp2 = 0;

			temp = fabs(T_av[ind*Mx + i] - exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), time));
			temp2 = fabs(exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), time));

			//if (ind == 0 && t == -1)
			//printf("temp = %f \n", temp);
			
			if ( fabs(temp) > diff_max_T_tn )
				diff_max_T_tn = fabs(temp);
			if ( fabs(temp2) > norm_max_T_tn )
				norm_max_T_tn = fabs(temp2);

			diff_l2_T_tn += temp*temp; 
			norm_l2_T_tn += temp2*temp2;
		}
	}

	diff_l2_T_tn = sqrt(diff_l2_T_tn/(dim_T));
	norm_l2_T_tn = sqrt(norm_l2_T_tn/dim_T);

	diffrel_max_T_tn = diff_max_T_tn / norm_max_T_tn;
	diffrel_l2_T_tn = diff_l2_T_tn / norm_l2_T_tn;

//	printf ( "diff_max_T_tn = %e diff_l2_T_tn = %e \n", diff_max_T_tn, diff_l2_T_tn );
//	printf ( "norm_max_T_tn = %e norm_l2_T_tn = %e \n", norm_max_T_tn, norm_l2_T_tn );
//	printf ( "diffrel_max_T_tn = %e diffrel_l2_T_tn = %e \n", diffrel_max_T_tn, diffrel_l2_T_tn );

	double diff_max_W_tn = 0.0;
	double diff_l2_W_tn = 0.0;
	double norm_max_W_tn = 0.0;
	double norm_l2_W_tn = 0.0;
	double diffrel_max_W_tn = 0.0;
	double diffrel_l2_W_tn = 0.0;

	double diff_max_Wz_tn = 0.0;
	double diff_l2_Wz_tn = 0.0;
	double norm_max_Wz_tn = 0.0;
	double norm_l2_Wz_tn = 0.0;
	double diffrel_max_Wz_tn = 0.0;
	double diffrel_l2_Wz_tn = 0.0;

	double diffz;
	temp2 = 0;
	for (int i=0; i<My*Mx; i++)
	{
		for (int k=0; k<Mz+1; k++)
		{
			i1 = i/Ny;
			j1 = i - i1*My;

			k1 = k;

			diffz = fabs(Wz_av[i*(Mz+1) + k] - (-heat_conductivity(i1,j1,k)*exact_gradientZ(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1], time) ));
			temp2 = fabs(-heat_conductivity(i1,j1,k)*exact_gradientZ(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1], time) );

			if (diffz > diff_max_Wz_tn)
				diff_max_Wz_tn = diffz;

			if (temp2 > norm_max_Wz_tn)
				norm_max_Wz_tn = temp2;

			diff_l2_Wz_tn += diffz*diffz; 
			norm_l2_Wz_tn += temp2*temp2;
		}
	}
	diff_l2_Wz_tn = diff_l2_Wz_tn/dim_Wz;
	norm_l2_Wz_tn = norm_l2_Wz_tn/dim_Wz;

	if ( diff_max_Wz_tn > diff_max_W_tn)
		diff_max_W_tn = diff_max_Wz_tn;
	if ( norm_max_Wz_tn > norm_max_W_tn)
		norm_max_W_tn = norm_max_Wz_tn;


	/////////////////
	double diff_max_Wx_tn = 0.0;
	double diff_l2_Wx_tn = 0.0;
	double norm_max_Wx_tn = 0.0;
	double norm_l2_Wx_tn = 0.0;
	double diffrel_max_Wx_tn = 0.0;
	double diffrel_l2_Wx_tn = 0.0;

	double diffx;
	temp2 = 0;
	for (int i = 0; i < My*Mz; i++)
	{
		for (int k = 0; k < Mx+1; k++)
		{
			j1 = i - (i/My)*My;
			k1 = i/My;

			diffx = fabs(Wx_av[i*(Mx+1) + k] - (-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), time) ));
			temp2 = fabs(-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), time) );

			//double tempval = (-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );
			//diffx = 1.0;

			if (diffx > diff_max_Wx_tn)
				diff_max_Wx_tn = diffx;
			if (temp2 > norm_max_Wx_tn)
				norm_max_Wx_tn = temp2;

			diff_l2_Wx_tn += diffx*diffx; 
			norm_l2_Wx_tn += temp2*temp2;
		}
	}
	diff_l2_Wx_tn = diff_l2_Wx_tn/dim_Wx;
	norm_l2_Wx_tn = norm_l2_Wx_tn/dim_Wx;

	if ( diff_max_Wx_tn > diff_max_W_tn)
		diff_max_W_tn = diff_max_Wx_tn;
	if ( norm_max_Wx_tn > norm_max_W_tn)
		norm_max_W_tn = norm_max_Wx_tn;

	/////////////////////
	double diff_max_Wy_tn = 0.0;
	double diff_l2_Wy_tn = 0.0;
	double norm_max_Wy_tn = 0.0;
	double norm_l2_Wy_tn = 0.0;
	double diffrel_max_Wy_tn = 0.0;
	double diffrel_l2_Wy_tn = 0.0;

	double diffy;
	temp2 = 0;
	for (int i = 0; i < Mz*Mx; i++)
	{
		for (int k = 0; k < My+1; k++)
		{
			i1 = (i-i%Mz)/Mz;
			j1 = k;
			k1 = i%Mz;

			diffy = fabs(Wy_av[i*(My+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), time) ));
			temp2 = fabs(-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), time) );

			if (diffy > diff_max_Wy_tn)
				diff_max_Wy_tn = diffy;

			if (temp2 > norm_max_Wy_tn)
				norm_max_Wy_tn = temp2;

			diff_l2_Wy_tn += diffy*diffy; 
			norm_l2_Wy_tn += temp2*temp2;
		}
	}
	diff_l2_Wy_tn = diff_l2_Wy_tn/dim_Wy;
	norm_l2_Wy_tn = norm_l2_Wy_tn/dim_Wy;

	if ( diff_max_Wy_tn > diff_max_W_tn)
		diff_max_W_tn = diff_max_Wy_tn;
	if ( norm_max_Wy_tn > norm_max_W_tn)
		norm_max_W_tn = norm_max_Wy_tn;

	diff_l2_W_tn = sqrt(diff_l2_Wx_tn + diff_l2_Wy_tn + diff_l2_Wz_tn);
	norm_l2_W_tn = sqrt(norm_l2_Wx_tn + norm_l2_Wy_tn + norm_l2_Wz_tn);

	//printf("after all eps_wn2 = %f \n", eps_wn2);
	//printf("after all norm_l2_w = %f \n", norm_l2_w);
	//printf("----------------------\n");

	diffrel_max_W_tn = diff_max_W_tn / norm_max_W_tn;
	diffrel_l2_W_tn = diff_l2_W_tn / norm_l2_W_tn;


	if ( diff_max_W_tn > diff_max_w)
		diff_max_w = diff_max_W_tn;

	if ( diff_l2_W_tn > diff_l2_w )
		diff_l2_w = diff_l2_W_tn;

	if ( diff_max_T_tn > diff_max_T)
		diff_max_T = diff_max_T_tn;

	if ( diff_l2_T_tn > diff_l2_T )
		diff_l2_T = diff_l2_T_tn;

	if ( diffrel_max_W_tn > diffrel_max_w )
		diffrel_max_w = diffrel_max_W_tn;

	if ( diffrel_max_T_tn > diffrel_max_T )
		diffrel_max_T = diffrel_max_T_tn;

	if ( diffrel_l2_W_tn > diffrel_l2_w )
		diffrel_l2_w = diffrel_l2_W_tn;
	if ( diffrel_l2_T_tn > diffrel_l2_T )
		diffrel_l2_T = diffrel_l2_T_tn;

	printf("time=%f tstep = %d \n",time, int(time/tau));
	if ( int(time/tau)-( int(time/tau)/print_step)*print_step==0 )
	{
		printf("time=%f tstep = %d \n",time, int(time/tau));
		printf("diff_max_T = %e \n",diff_max_T); 
		printf("diff_l2_T = %e \n",diff_l2_T); 
		printf("diffrel_max_T = %e \n",diffrel_max_T); 
		printf("diffrel_l2_T = %e \n",diffrel_l2_T); 

		printf("diff_max_w = %e \n",diff_max_w); 
		printf("diff_l2_w = %e \n",diff_l2_w); 
		printf("diffrel_max_w = %e \n",diffrel_max_w); 
		printf("diffrel_l2_w = %e \n",diffrel_l2_w); 


		printf("n-th step \n");
		printf("diff_l2_W_tn = %e \n", diff_l2_W_tn);
		printf("diff_l2_T_tn = %e \n", diff_l2_T_tn);

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

	*diffPr_max_w_pt = diff_max_w;
	*diffPr_l2_w_pt = diff_l2_w;
	*diffrelPr_max_w_pt = diffrel_max_w;
	*diffrelPr_l2_w_pt = diffrel_l2_w;
	*diffPr_max_T_pt = diff_max_T;
	*diffPr_l2_T_pt = diff_l2_T;
	*diffrelPr_max_T_pt = diffrel_max_T;
	*diffrelPr_l2_T_pt = diffrel_l2_T;

	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int	Accuracy_calculateL2Pr(double * T_n, double * Wx_n, double * Wy_n, double * Wz_n,
						   double * T_Av, double * Wx_Av, double * Wy_Av, double * Wz_Av,
						   boundaryConditions xCondition, boundaryConditions yCondition, boundaryConditions zCondition,
						   int numSolution, int t, int Mx, int My, int Mz, int print_step, 
						   double * eps_max_w_pt, double * eps_max_T_pt, double * eps_l2_w_pt, double * eps_l2_T_pt,
						   double * eps_relative_max_w_pt, double * eps_relative_max_T_pt, double * eps_relative_l2_w_pt, double * eps_relative_l2_T_pt)
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
	double tz = 0;
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

	for ( int i = 0 ; i < Mx ; i++ )
	{
		for (int j = 0; j < My; j++)
		{
			for (int k = 0; k < Mz; k++)
			{
				T_Av[i + j*Nx + k*Nx*Ny] = exact_solution_averageQuad ( numSolution, (t+1)*tau, i, j, k, xpoints, ypoints, zpoints, NQUADPOINTS );
			}
		}
	}
	HeatFlux_Init(T_Av, Wx_Av, xCondition, Wy_Av, yCondition, Wz_Av, zCondition, numSolution, Mx, My, Mz);


	int i_max_T = 0;
	int j_max_T = 0;
	int k_max_T = 0;
	int i1, j1, k1;
	int dim_T = Mx * My * Mz;
	int dim_Wx = (Mx + 1) * My * Mz;
	int dim_Wy = Mx * (My + 1) * Mz;
	int dim_Wz = Mx * My * (Mz + 1);

	for (int ind = 0; ind < My*Mz; ind++)
	{
		j1 = ind - (ind/My)*My;
		k1 = ind/My;
		for ( int i = 0 ; i < Mx ; i++ )
		{
			temp = 0;
			temp2 = 0;

#ifdef SPECIALTEST
			temp = fabs(T_n[ind*Mx + i]);
			temp2 = temp;
#else
			temp = fabs(T_n[ind*Mx + i] - T_Av[ind*Mx + i]);
			temp2 = fabs(T_Av[ind*Mx + i]);
#endif
			//if (ind == 0 && t == -1)
			//printf("temp = %f \n", temp);
			
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

	//diff_l2z = 0.0;
	//norm_l2_wz = 0.0;
	int imax_z = 0;
	int jmax_z = 0;
	int kmax_z = 0;
	temp2 = 0;

	for (int i=0; i<Ny*Nx; i++)
	{
		for (int k=0; k<Nz+1; k++)
		{
			//i1 = i%Nx;
			//j1 = (i-i%Nx)/Nx;
			i1 = i/Ny;
			j1 = i - i1*Ny;

			k1 = k;

#ifdef SPECIALTEST
			diffz = fabs(Wz_n[i*(Nz+1) + k]);
			temp2 = diffz;
#else
			diffz = fabs(Wz_n[i*(Nz+1) + k] - Wz_Av[i*(Nz+1) + k]);
			temp2 = fabs(Wz_Av[i*(Nz+1) + k]);
#endif
			if (diffz > tz)
			{
				tz = diffz;
				imax_z = i1;
				jmax_z = j1;
				kmax_z = k1;
			}

			if (temp2 > norm_max_w)
				norm_max_w = temp2;

			diff_l2z += diffz*diffz; 
			norm_l2_wz += temp2*temp2;
		}
	}
	diff_l2z = diff_l2z/dim_Wz;
	norm_l2_wz = norm_l2_wz/dim_Wz;

	if ( tz > eps_n1_w)
		eps_n1_w = tz;

	//printf("----------------------\n");
	//printf("after wz_loop eps_n1_w = %f \n", eps_n1_w);
	//printf("after wz_loop diff_l2z = %f \n", diff_l2z);
	//printf("after wz_loop norm_l2_wz = %f \n", norm_l2_wz);


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

#ifdef SPECIALTEST
			diffx = fabs(Wx_n[i*(Mx+1) + k]);
			temp2 = diffx;
#else
			diffx = fabs(Wx_n[i*(Mx+1) + k] - Wx_Av[i*(Mx+1) + k]);
			temp2 = fabs(Wx_Av[i*(Mx+1) + k]);
#endif
			//double tempval = (-heat_conductivity_func(xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), (t+1)*tau) );
			//diffx = 1.0;

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

		//printf("after wx_loop eps_n1_w = %f \n", eps_n1_w);
		//printf("after wx_loop diff_l2x = %f \n", diff_l2x);
		//printf("after wx_loop norm_l2_wx = %f \n", norm_l2_wx);


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

#ifdef SPECIALTEST
				diffy = fabs(Wy_n[i*(My+1) + k]);
				temp2 = diffy;
#else
				diffy = fabs(Wy_n[i*(My+1) + k] - Wy_Av[i*(My+1) + k]);
				temp2 = fabs(Wy_Av[i*(My+1) + k]);
#endif
				if (diffy > ty)
				{
					ty = diffy;
					imax_y = i1;
					jmax_y = j1;
					kmax_y = k1;
				}

//				printf ( "i1 = %d j1 = %d k1 = %d \n", i1,j1,k1);
//				printf ( "val_numer = %3.3e val_ex = %3.3e \n", Wy_n[i*(My+1) + k], (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution, xpoints[i1]+0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1), (t+1)*tau) ) ); 

				if (temp2 > norm_max_w)
					norm_max_w = temp2;

				diff_l2y += diffy*diffy; 
				norm_l2_wy += temp2*temp2;
			}
			diff_l2y = diff_l2y/dim_Wy;
			norm_l2_wy = norm_l2_wy/dim_Wy;

			if ( ty > eps_n1_w)
				eps_n1_w = ty;

			//printf("after wy_loop eps_n1_w = %f \n", eps_n1_w);
			//printf("after wy_loop diff_l2z = %f \n", diff_l2y);
			//printf("after wy_loop norm_l2_wz = %f \n", norm_l2_wy);


			eps_wn2 = sqrt(diff_l2x + diff_l2y + diff_l2z);
			norm_l2_w = sqrt(norm_l2_wx + norm_l2_wy + norm_l2_wz);

			//printf("after all eps_wn2 = %f \n", eps_wn2);
			//printf("after all norm_l2_w = %f \n", norm_l2_w);
			//printf("----------------------\n");

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
			if (t == -1)
			{
				printf("t=%d inside function Accuracy \n",t);
				printf("eps_max_T = %e \n",eps_max_T); 
				printf("eps_l2_T = %e \n",eps_l2_T); 
				printf("eps_max_w = %e \n",eps_max_w); 
				printf("eps_l2_w = %e \n",eps_l2_w); 

				printf("eps_wx_l2 = %f \n", sqrt(diff_l2x));
				printf("eps_wy_l2 = %f \n", sqrt(diff_l2y));
				printf("eps_wz_l2 = %f \n", sqrt(diff_l2z));

			}
			if (t+1-((t+1)/print_step)*print_step==0 && t != -1)
			{
				printf("t=%d \n",t);
				printf("eps_max_T = %e \n",eps_max_T); 
				printf("eps_l2_T = %e \n",eps_l2_T); 
				printf("eps_max_w = %e \n",eps_max_w); 
				printf("eps_l2_w = %e \n",eps_l2_w); 

				printf("eps_wx_l2 = %f \n", sqrt(diff_l2x));
				printf("eps_wy_l2 = %f \n", sqrt(diff_l2y));
				printf("eps_wz_l2 = %f \n", sqrt(diff_l2z));

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
