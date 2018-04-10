#include<stdio.h>
//#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<string.h>
//#include<conio.h>
//#include<windows.h>
//#include<time.h>
//#include<sys/types.h>

// if activated, mesh heat flux W0 is obtained as a pointwise projection of the exact heat flux at t = 0.
#define EXACT_W0_INIT_H
// if activated, mesh temperature T0 is taken NOT as a pointwise projection of the exact temperature at t = 0
// but as a CELL average for each mesh cell.
// Currently works only for numsol = 501
#define AVERAGE_T0_INIT
// if AVERAGE_T0_INIT is active, then defines number of points used for evaluating an integral of T0 over each mesh cell
#define NQUADPOINTS (10)
// if active, then additional comparison between pointwise and L2 projector of the exact analytical flux is carried out
// works only for homogeneous boundary conditions
#define COMPARE_PROJECTORS_H

// if active, the error is measured between numerical flux and L2 projection of the exact flux
// if passive, the error is measured between numerical flux and pointwise projection of the exact flux
#define ERROR_FOR_L2PROJECTION_H

#define SPECIALTEST_H


#include "data.h"
#include "solution.h"
#include "key_functions3d.h"
#include "eigenvectors3D.h"

#define BNDCNDS boundaryConditions

/*
Программа написана для трехмерного случая, для уравнения на тепловой поток реализована схема Дугласа Гана
для сеточной дивергенции, записанная в терминах тепловго потока.
Одномерный вариант: MixedFEM krank nikolson.
Двумерный вариант: MixedFem_krank_nikolson_splitting2D.
Состояние:

09.11.2011
Добавил блоки инциализации по z для краевых условий  Дирихле и Неймана и проверил их.
Добавил все блоки.
На тесте 1266 по x держится 0, считает правдоподобно, но неверно, ищу ошибку.
Такое же поведение для теста 266 (по z поток = 0)
!!! Кое-где только однородные условия Неймана сработают.

Нашел ошибку в самой схеме, исправил. Тем не менее на num_sol = 266 расходится неспешно.
Проверил, что если пересчитывать с 0.5*tau*delta_Wxplus13, то считает step 1 так же, как м
двумерная программа MixedFem_2D_splitting.

Решил проверить на одномерном тесте numsol = 66. Поставил правильную правую часть, f_n и f_n+1-f_n /2 для y
и одномерный случай заработал на постоянных к-тах, по крайней мере. C external = 1 не работает, поток по x становится ненулевым
из-за потока по y. Пока ошибку не нашел. Нашел одну ошибку, не было для numsol66 инициализации правой части
в прогонке по x.

15.11.2011
Теперь работает при tau = 1e0у-6 external 0, 1, 2, 3 на тестах 66 и на тесте 106 (Дирихле и Нейман 1D y).

16.11.2011
Проверил с tau = 1.0e-6 external 0, 1, 2, 3 на тесте 6, external 0, 1, 2, 3 на тесте 366.
Проверил с tau = 1.0e-6 external 0, 1, 2, 3 на тестах 2006 и 2066 (Нейман и Дирихле по z) 
Заодно сделал удобную систему изменения краевых условий, ввел переменные и свитчи.

17.11.2011
Проверил с tau = 1.0e-6 external 0, 1, 2, 3 на тесте 206.(Нейман х + Нейман y)

18.11.2011
Проверил с tau = 1.0e-6 external 0, 1, 2, 3 на тесте 266, исправил одну ошибочку.(Нейман х + Дирихле y)
Проверил с tau = 1.0e-6 external 0, 1, 2, 3 на тесте 1266 (Нейман z + Дирихле y)
Проверил с tau = 1.0e-6 external 0, 1, 1 (c заменой x на z), 2, 3 на тесте 1206 (Нейман z + Нейман y)


22.11.2011
Проверил с tau = 1.0e-6 external 0, 1, 2, 2 (c заменой y на z), 3 на тесте 2206 (Нейман x + Нейман z)
Проверил с tau = 1.0e-6 external 0, 1, 2, 2 (c заменой y на z), 3 на тесте 1366 (Нейман x + Дирихле z)
Проверил с tau = 1.0e-6 external 0, 1, 1 (c заменой x на z), 2, 3 на тесте 1466 (Нейман y + Дирихле z)
Проверил с tau = 1.0e-6 external 0, 1, 2, 3 на тесте 1566 (Нейман y + Дирихле x)
Проверил с tau = 1.0e-6 external 0, 1, 2, 2 (c заменой y на z), 3 на тесте 2366 (Нейман z + Дирихле x)

Тестирование двумерных решений завершено, перехожу к трехмерным.
Проверил с tau = 1.0e-6 external 0 на тесте 9006 (Нейман x + Нейман y + Нейман z)

07.12.2011
Проверил с tau = 1.0e-6 external 0 на тесте 9066 (Дирихле x + Нейман y + Нейман z)
Готовлю программу для тетсирования на кластере.

!!!!!!!!!!!Почему-то не работает нестационарный тест 366 с external = 2, пока ошибки не нашел.
Нашел ошибку, был косяк в номере правой части во второй прогонке по y. Теперь работает.
!!!!!!!!!!!Почему-то с разрывом плтности по x тест 1206 расходится, за счет потока по х, который вообще
должен быть равен 0. Причем с 1266 несмотря на большой поток по x, расходимости нет.
Нашел ошибку, был косяк в номере правой части в прогонке по x
!!!!!!!!!!! тест 1466 не работает, на external 0, уже в deltaWznplus13 большие значения.
Исправил ошибку.
!!!!!!!!!!!Почему-то на тесте 266 с external = 0 погрешность по потоку в какой-то момент после
уменьшения начинает увеличиваться..хотя в итоге вроде начинает обратно уменьшаться, температура
стационируется на привычных цифрах. Аналогично с external = 1, 2. И с другими двумерными тестами.
Непонятно.

23.05.2012
Готовлю программу для передачи Яровенко для создания параллельной реализации. Решили поменять нумерацию для компонент
потока: Wx -- x-y-z; Wy -- y-z-x; Wz -- z-y-x. Сохранил работающую версию кода как 3Dsplitting_douglasgann_cyclic.cpp.
Нумерацию поменял.
С краевыми условиями решил пока не заморачиваться, чистый Нейман и чистый Дирихле по всем компонентам. Потому что 
в реальной задаче пока неясно, какие краевые условия будут.
Убираю явно лишние массивы Wx_nplus13, Wx_nplus23, delta_Wx_nplus23, delta_Wx_nplus1; Wy_nplus23, delta_Wy_nplus1;
Wz_nplus23, delta_Wz_nplus23.

06.12.2014
Сделал конструктор. Работает трехмерная схема расщепления на основе схемы Дугласа-Ганна, есть сходимость со 2 порядком.
Добавил конвекцию (через тепловой поток). Сделаны num_sol = 1166 со скоростью, 9066. Сделал локально-одномерную схему, сравнил на тесте 1166.
Есть сходимость с 1 порядком.

07.12.2014
Сделал схему предиктор-корректор. Есть какая-то ошибка, см описание перед ..._mainstep.

08.12.2014
Исправил ошибку для схемы предиктор-корректор, теперь есть сходимость на 1166 со вторым порядком.
Сделал схему на основе алгоритма Удзавы3D, поток считает так же, как и пред.-корр., а температуру лучше.
Похоже, что устойчива. Есть 2 порядок для теста 1166.

27.09.2016
Для статьи по трехмерному предиктор-корректору нужны результаты численных экспериментов.
Починил пару багов - один в условиях Дирихле по y, в Implicit_Wy(), другой в краевых условиях по z: yCondition стояло вместо zCondition.
Сделал тест с произведением синусов (999), теперь на нем второй порядок.
Вроде у предиктор-корректора тоже.

28.09.2016
Сделал тест с переменной гладкостью (888), исправил еще баги в BzMminus1_F_Add и в самом тестовом решении. Кажется, работает, есть сходимость
со вторым порядком при достаточном уровне гладкости.
Сделал тест где синусы и по одной переменной s * sin (pi s) (через dir_axis, 777) и тест, где по всем переменным так (776).
Точность у Дугласа и Ганна низкая, но с первым порядком все-таки сходится поток.
Зато на тесте 500 с неразделяющимися переменными Дуглас Ганн таки не сходится, этот тест по идее подойдет для статьи.

Исправил ошибку в Invert_Ay, испольующейся в подсчете нормы в H~, сделал осреднение (вкл-выкл по флагу) с исп-м квадратуры,
сделал замеры для различных тестов, но со статьей проблемы.

Сделал тест 889, негладкий по одному направлению и с синусами по остальным, на нем непонятно все.
То ли тест с косяком, то ли что-то еще мешает - нормы ведут себя плохо (и для предиктор-корректора),
сходимость для гладкого случая есть при этом, и цифры для гладкого случая норм, а для негладкого цифры совсем большие.

Создаю тест переменной гладкости с cos^p (444). Работает довольно отвратно во всех смыслах.
Пока его результаты решено не использовать.
*/

double tau;
//int Nx, Ny, Nz;
//int external;
int num_sol;
int scheme;

double* Wx_n;			// массив для значений x-компоненты теплового потока на n-ом временном слое
double* Wy_n;			// массив для значений y-компоненты теплового потока на n-ом временном слое
double* Wz_n;			// массив для значений z-компоненты теплового потока на n-ом временном слое
double* T;				// массив для значений температуры (в процессе вычислений - на n-ом временном слое, но в конце основного цикла хранит результирующие значения температуры

double *xpoints, *ypoints, *zpoints;


//void alfa_calc(int n, double* alfa);		//вычисляет коэфф-ты альфа i-e для применения прогонки к трехдиагональной матрице (см.описание матриц) размера n
//void alfaz_calc(int n, double* alfa);		// --||-- для z

void allocation();
double righthand_my(int k);
double norm_l2(double * massiv, int N, int M, int K);
double scal_l2(double * massiv1, double * massiv2, int N, int M, int K);
double norm_max(double * massiv, int N, int M, int K);
void stability_normFull(double* wx, double* wy, double * wz, int Nx, int Ny, int Nz, double *l2_norm_pt);

int stability_norm3Dxyz ( double * Wx, double * Wy, double * Wz, int Nx, int Ny, int Nz, double tau, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition, double * res_pt);
int stability_norm3Dyzx ( double * Wx, double * Wy, double * Wz, int Nx, int Ny, int Nz, double tau, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition, double * res_pt);
int zero_init ( double * vec, int dim );
int compute_CTxyz ( double * CT, double * T, double tau, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition );
int compute_CTyzx ( double * CT, double * T, double tau, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition );
int	compute_Lambdax ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition);
int	compute_Lambday ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition);
int	compute_Lambdaz ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition);
int	compute_Lambda ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition);
int Invert_Ax(double * sol, double * rhand, boundaryConditions xCondition, int Mx, int My, int Mz);
int Invert_Ay(double * sol, double * rhand, boundaryConditions yCondition, int Mx, int My, int Mz);
int Invert_Az(double * sol, double * rhand, boundaryConditions zCondition, int Mx, int My, int Mz);


FILE* f1 = NULL; //файл для выходных данных


int main(int argc, char **argv)
{
	if(argc == 1) // no parameters
	{
		printf("Enter parameters tau, Nx, Ny, Nz, external and num_sol \n");
		scanf("%lf", &tau);
		scanf("%d", &Nx);
		scanf("%d", &Ny);
		scanf("%d", &Nz);
		scanf("%d", &external);
		scanf("%d", &num_sol);
		scanf("%d", &scheme);
	}else if(argc == 8)
	{
		tau = atof(argv[1]);
		Nx = atoi(argv[2]);
		Ny = atoi(argv[3]);
		Nz = atoi(argv[4]);
		external = atoi(argv[5]);
		num_sol = atoi(argv[6]);
		scheme = atoi(argv[7]);
	} else
	{
		printf("Incorrect parameters \n");
		exit(-1);
	}

	printf ( "tau = %f Nx = %d Ny = %d Nz = %d external = %d num_sol = %d scheme = %d \n", tau, Nx, Ny, Nz, external, num_sol, scheme ); 
	allocation(); //set x(y,z)points 
	setBoundaryConditions(num_sol); //eDirichlet eNeumann eMixed
	printf("xCondition = %d yCondition = %d zCondition = %d \n",xCondition, yCondition, zCondition);


	char filename[60];
	char addname1[15];
	char addname2[15];
	char addname3[15];
	char addname4[15];
	char addname5[15];
	bool neu2D = (  (((num_sol == 6) || (num_sol == 7)) || ((num_sol == 206) || (num_sol == 207))) || (num_sol == 106));
	bool dir2D = (((num_sol == 66) || (num_sol == 67)) || ((num_sol == 266) || (num_sol == 267)));
	bool neu3D = (num_sol == 9006);
	bool dir3D = (num_sol == 9066 || num_sol == 999 );
	if ((neu2D || dir2D) || (neu3D || dir3D))
	{
		if (neu2D)
			strcpy(filename,"resultKrank2DNeu_");
		if (dir2D)
			strcpy(filename,"resultKrank2DDir_");
		if (neu3D)
			strcpy(filename,"resultKrank3DNeu_");
		if (dir3D)
			strcpy(filename,"resultKrank3DDir_");
	}
	else
		strcpy(filename,"resultKrank_wtf_");
	sprintf(addname5,"%d_",external);
	sprintf(addname1,"%f_",tau);
	sprintf(addname2,"%d_",Nx);
	sprintf(addname3,"%d_",Ny);
	sprintf(addname4,"%d",Nz);
	strcat(filename,addname5);
	strcat(filename,addname1);
	strcat(filename,addname2);
	strcat(filename,addname3);
	strcat(filename,addname4);
	strcat(filename,".txt");
	f1 = fopen(filename,"w");

	double* xpoints;
	xpoints = new double [Nx+1];
	xpoints[0] = 0;
	for (int i=1; i<=Nx ; i++)
	{
		xpoints[i] = xpoints[i-1] + hx(i-1);
	}
	double* ypoints;
	ypoints = new double [Ny+1];
	ypoints[0] = 0;
	for (int j=1; j<=Ny ; j++)
	{
		ypoints[j] = ypoints[j-1] + hy(j-1);
	}
	double* zpoints;
	zpoints = new double [Nz+1];
	zpoints[0] = 0;
	for (int k=1; k<=Nz ; k++)
	{
		zpoints[k] = zpoints[k-1] + hz(k-1);
	}

	int N = int(Time/tau);//количество шагов по времени

	//	размерности используемых массивов
	int dim_T = Nx*Ny*Nz;
	int dim_Wx = (Nx+1)*Ny*Nz;
	int dim_Wy = Nx*(Ny+1)*Nz;
	int dim_Wz = Nx*Ny*(Nz+1);
	int dim_W = dim_Wx + dim_Wy + dim_Wz;
	//	printf("dim_T = %d \ndim_Wx = %d \ndim_Wy = %d \ndim_Wz = %d \n",dim_T, dim_Wx, dim_Wy,dim_Wz);

	///////Определяем размеры массивов для компонент скорости и теплового потока
	Wx_n = new double [dim_Wx];
	Wy_n = new double [dim_Wy];
	Wz_n = new double [dim_Wz];
	T = new double [dim_T];

#if defined(COMPARE_PROJECTORS) || defined(ERROR_FOR_L2PROJECTION)
	double * T_averaged = new double [dim_T];
	double * Wx_n_Taver = new double [dim_Wx];
	double * Wy_n_Taver = new double [dim_Wy];
	double * Wz_n_Taver = new double [dim_Wz];
#endif
	///////
	printf("Initializing T...\n");

#ifdef AVERAGE_T0_INIT
	printf ( "ATTENTION: flag AVERAGE_T0_INIT is active, number of quadrat. points = %d \n", NQUADPOINTS );
#endif
	for ( int i = 0 ; i < Nx ; i++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
#ifdef AVERAGE_T0_INIT
				T[i + j*Nx + k*Nx*Ny] = exact_solution_averageQuad ( num_sol, 0, i, j, k, xpoints, ypoints, zpoints, NQUADPOINTS );
#else
				T[i + j*Nx + k*Nx*Ny] = exact_solution(num_sol,xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0);
#endif
#if defined(COMPARE_PROJECTORS) || defined(ERROR_FOR_L2PROJECTION)
				T_averaged[i + j*Nx + k*Nx*Ny] = exact_solution_averageQuad ( num_sol, 0, i, j, k, xpoints, ypoints, zpoints, NQUADPOINTS );
#endif

			}
		}
	}


#ifdef COMPARE_PROJECTORS
	printf ( "ATTENTION: COMPARE_PROJECTORS is activated \n" );
#endif

#ifdef ERROR_FOR_L2PROJECTION
	printf ( "ATTENTION: ERROR_FOR_L2PROJECTION is activated \n" );
#endif

#ifdef SPECIALTEST
	//printf ( "ATTENTION: SPECIALTEST ia activated: No f, T = 0, w = w_exact - w_exact_proj NUMSOLSPECIAL = %d \n", NUMSOLSPECIAL );
	printf ( "ATTENTION: SPECIALTEST ia activated: No f, T = 0, w = w_exact - w_exact_proj \n" );
#endif

	if ( num_sol == 888 )
		printf ( "DEG888 = %f A888 = %f C888 = %f \n", DEG888, A888, C888 );
	if ( num_sol == 889 )
	{
		printf ( "DEG888 = %f A888 = %f C888 = %f \n", DEG888, A888, C888 );
        printf ( "main_axis = %d mx = %d my = %d mz = %d \n", Dir_axis, mx, my, mz );
	}
	if ( num_sol == 500 )
		printf ( "MT500 = %f \n", MT500 );
	if ( num_sol == 999 )
		printf ( "mx = %d my = %d mz = %d \n", mx, my, mz );
	if ( num_sol == 444 )
		printf ( "mx = %d my = %d mz = %d (must be a multiple of 2), DEG444 = %f \n", mx, my, mz, DEG444 );


	//////////////////////////
	printf("Initializing W...\n");
	double *alfax = new double[Nx + 1];
	double *betax = new double[Nx + 1]; 
	double *alfay = new double[Ny + 1];
	double *betay = new double[Ny + 1]; 
	double *alfaz = new double[Nz + 1];
	double *betaz = new double[Nz + 1]; 

	{;}

	// расчет из Aw - BT = G теплового потока w0 ~ инициализация w0
	////////////////////////////////////////***************** Step 0
	//T0 - - > W0 = (Wx_0, Wy_0, Wz_0 )transp
	//( W0 = A(-1) * ( G - B * Tn ) progonka!
	{;}

	double a_0_old, a_i_old, a_iminus1_old, a_n_old = 0;
	double b_0_old, b_i_old, b_n_old = 0;
	double temper_0, temper_k, temper_n = 0;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	//Flux initialization, old
	{;}
/*
	switch(bnd_X)
	{
	case 0:
		//правильное для неоднородных условий Неймана по x
		//инициализация для потока по x - Wx_n
		//Нейман, x
		for (int i = 0; i < Ny*Nz; i++)
		{
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			betax[0] = 0 ;
			betax[1] = wx_0 ;
			alfax[0] = 0 ;
			alfax[1] = 0 ;
			for ( int k = 1 ; k < Nx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = T[k - 1 + j1*Nx + k1*Nx*Ny] - T[k + j1*Nx + k1*Nx*Ny];

				if (k==Nx-1)
				{
					alfax[k + 1] = 0;
					betax[k + 1] = ( ( temper_k - a_i_old*wx_1 ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( ( temper_k - a_iminus1_old*wx_0 ) ) / (  b_i_old );
				}
				else
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( ( temper_k ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
				}
			}
			Wx_n[i*(Nx+1) + Nx ] = wx_1 ;
			for ( int j = 1 ; j < Nx + 1 ; j++ )
			{
				Wx_n[i*(Nx+1) + Nx - j] = Wx_n[i*(Nx+1) + Nx + 1 - j] * alfax[Nx + 1 - j] + betax[Nx + 1 - j];
			}
		}
		break;
	case 1:
		//для неоднородных условий Дирихле
		//Дирихле, x
		{;}

		for (int i = 0; i < Ny*Nz; i++)
		{
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;
			//к-ты матрицы А
			a_0_old = hx(0)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;

			temper_0 = Tx_0 - T[0 + j1*Nx + k1*Nx*Ny];
			betax[0] = 0 ;
			alfax[0] = 0 ;
			alfax[1] = -0.5 ;
			betax[1] =  temper_0 / b_0_old; 
			for ( int k = 1 ; k < Nx ; k++ )
			{
				i1 = k;
				//к-ты матрицы А
				a_i_old = hx(k)/(heat_conductivity(i1,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)/(heat_conductivity(i1-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;
				//правая часть
				temper_k = T[k - 1 + j1*Nx + k1*Nx*Ny] - T[k + j1*Nx + k1*Nx*Ny];

				alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
				betax[k + 1] = ( ( temper_k ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
			}
			a_n_old = hx(Nx-1)/(heat_conductivity(Nx-1,j1,k1)*6.0);
			b_n_old = 2*a_n_old;
			temper_n = T[Nx-1 + j1*Nx + k1*Nx*Ny] - Tx_1;

			Wx_n[i*(Nx+1) + Nx] = ( temper_n  -  betax[Nx]*a_n_old ) / (alfax[Nx]*a_n_old + b_n_old);
			for ( int j = 1 ; j < Nx + 1 ; j++ )
			{
				Wx_n[i*(Nx+1) + Nx - j] = Wx_n [i*(Nx+1) + Nx + 1 - j] * alfax[Nx + 1 - j] + betax[Nx + 1 -j]; 
			}
		}
		break;
	default:
		printf("Bad bnd_X because of wrong Bnd_x\n");
		return -1;
		break;
	}

	switch(bnd_Y)
	{
	case 0:
		//инициализация для потока по y - Wy_n
		//Нейман, y
		{;}
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;

			betay[0] = 0 ;
			betay[1] = wy_0 ;
			alfay[0] = 0 ;
			alfay[1] = 0 ;
			for ( int k = 1 ; k < Ny ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = (T[i1 + (k - 1)*Nx + k1*Nx*Ny] - T[i1 + k*Nx + k1*Nx*Ny]);

				if (k==Ny-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = ( ( temper_k - a_i_old*wy_1 ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k - a_iminus1_old*wy_0 ) ) / (  b_i_old );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
			}
			Wy_n[i*(Ny+1) + Ny ] = wy_1 ;
			for ( int j = 1 ; j < Ny + 1 ; j++ )
			{
				Wy_n[i*(Ny+1) + Ny - j ] = Wy_n [i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
			}
		}
		break;
	case 1:
		//для неоднородных условий Дирихле
		//Дирихле, y
		{;}
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;
			//к-ты матрицы А
			a_0_old = hy(0)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = Ty_0 - T[i1 + 0*Nx + k1*Nx*Ny];

			betay[0] = 0 ;
			alfay[0] = 0 ;
			alfay[1] = -0.5 ;
			betay[1] =  temper_0 / b_0_old; 

			for ( int k = 1 ; k < Ny ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = T[i1 + (k - 1)*Nx + k1*Nx*Ny] - T[i1 + k*Nx + k1*Nx*Ny];

				alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
				betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
			}
			a_n_old = hy(Ny-1)/(heat_conductivity(i1,Ny-1,k1)*6.0);
			b_n_old = 2*a_n_old;
			temper_n = T[i1 + (Ny-1)*Nx + k1*Nx*Ny] - Ty_1;

			Wy_n[i*(Ny+1) + Ny] = ( temper_n  -  betay[Ny]*a_n_old ) / (alfay[Ny]*a_n_old + b_n_old);
			for ( int j = 1 ; j < Ny + 1 ; j++ )
			{
				Wy_n[i*(Ny+1) + Ny - j ] = Wy_n [i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
			}
		}
		break;
	default:
		printf("Bad bnd_Y because of wrong Bnd_y\n");
		return -1;
		break;
	}

	switch(bnd_Z)
	{
	case 0:
		//инициализация для потока по z - Wz_n
		//Нейман, z (новое, проверено)
		{;}
		for (int i = 0; i < Ny*Nx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/Ny;
			j1 = i - i1*Ny;

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
				temper_k = T[i1 + j1*Nx + (k-1)*Nx*Ny] - T[i1 + j1*Nx + k*Nx*Ny];

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
			Wz_n[i*(Nz+1) + Nz ] = wz_1 ;
			for ( int j = 1 ; j < Nz + 1 ; j++ )
			{
				Wz_n[i*(Nz+1) + Nz - j] = Wz_n[i*(Nz+1) + Nz + 1 - j] * alfaz[Nz + 1 - j] + betaz[Nz + 1 - j];
			}
		}
		break;
	case 1:
		//инициализация для потока по z - Wz_n
		//Дирихле, z (новое, проверено)
		{;}

		for (int i = 0; i < Ny*Nx; i++)
		{
			//j1 = i/Nx;
			//i1 = i - (i/Nx)*Nx;
			i1 = i/Ny;
			j1 = i - i1*Ny;

			//к-ты матрицы А
			a_0_old = hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = Tz_0 - T[i1 + j1*Nx + 0*Nx*Ny];

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
				temper_k = T[i1 + j1*Nx + (k-1)*Nx*Ny] - T[i1 + j1*Nx + k*Nx*Ny];

				alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
				betaz[k + 1] =  ( ( temper_k ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				//			printf("temper_k = %f \n",temper_k);
			}
			a_n_old = hz(Nz-1)/(heat_conductivity(i1,j1,Nz-1)*6.0);
			b_n_old = 2*a_n_old;
			temper_n = T[i1 + j1*Nx + (Nz-1)*Nx*Ny] - Tz_1;
			//		printf("temper_n = %f \n",temper_n);

			Wz_n[i*(Nz+1) + Nz] = ( temper_n  -  betaz[Nz]*a_n_old ) / (alfaz[Nz]*a_n_old + b_n_old);
			for ( int j = 1 ; j < Nz + 1 ; j++ )
			{
				Wz_n[i*(Nz+1) + Nz - j] = Wz_n[i*(Nz+1) + Nz + 1 - j] * alfaz[Nz + 1 - j] + betaz[Nz + 1 - j];
			}
		}
		break;
	default:
		printf("Bad bnd_Z because of wrong Bnd_z\n");
		return -1;
		break;
	}
*/
	{;}
	HeatFlux_Init(T, Wx_n, xCondition, Wy_n, yCondition, Wz_n, zCondition, num_sol, Nx, Ny, Nz);

#if defined(COMPARE_PROJECTORS) || defined(ERROR_FOR_L2PROJECTION)
	HeatFlux_Init(T_averaged, Wx_n_Taver, xCondition, Wy_n_Taver, yCondition, Wz_n_Taver, zCondition, num_sol, Nx, Ny, Nz);
#endif


	double *Wx_0_exact = (double*)malloc(dim_Wx*sizeof(double));
	for (int i = 0; i < Ny*Nz; i++)
	{
		for (int k = 0; k < Nx+1; k++)
		{
			i1 = k;
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			Wx_0_exact[i*(Nx+1) + k] = (-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(num_sol,xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) );
		}
	}
	double *Wy_0_exact = (double*)malloc(dim_Wy*sizeof(double));
	for (int i=0; i<Nz*Nx; i++)
	{
		for (int k=0; k<Ny+1; k++)
		{
			i1 = i/Nz;
			j1 = k;
			k1 = i - (i/Nz)*Nz;

			Wy_0_exact[i*(Ny+1) + k] = (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(num_sol,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) );
		}
	}

	double *Wz_0_exact = (double*)malloc(dim_Wz*sizeof(double));
	for (int i=0; i<Ny*Nx; i++)
	{
		for (int k=0; k<Nz+1; k++)
		{
			i1 = i/Ny;
			j1 = i - i1 * Ny;
			k1 = k;

			Wz_0_exact[i*(Nz+1) + k1] = (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[j1]+0.5*hy(j1), zpoints[k1])*exact_gradientZ(num_sol,xpoints[i1]+0.5*hx(i1),ypoints[j1]+0.5*hy(j1), zpoints[k1], 0) );
		}
	}

#ifdef	EXACT_W0_INIT
	printf ( "ATTENTION: flag EXACT_W0_INIT is active \n" );
	for (int i = 0; i < dim_Wx; i++)
		Wx_n[i] = Wx_0_exact[i];
	for (int i = 0; i < dim_Wy; i++)
		Wy_n[i] = Wy_0_exact[i];
	for (int i = 0; i < dim_Wz; i++)
		Wz_n[i] = Wz_0_exact[i];
#endif

	FILE * Wx_0_file_NEW = fopen("Wx_0_NICE_interface_3d.xls","wt");
	for (int i = 0; i < dim_Wx; i++)
		fprintf(Wx_0_file_NEW,"%f \n", Wx_n[i]);
	fclose(Wx_0_file_NEW);
	FILE * Wy_0_file_NEW = fopen("Wy_0_NICE_interface_3d.xls","wt");
	for (int i = 0; i < dim_Wy; i++)
		fprintf(Wy_0_file_NEW,"%f \n", Wy_n[i]);
	fclose(Wy_0_file_NEW);
	FILE * Wz_0_file_NEW = fopen("Wz_0_NICE_interface_3d.xls","wt");
	for (int i = 0; i < dim_Wz; i++)
		fprintf(Wz_0_file_NEW,"%f \n", Wz_n[i]);
	fclose(Wz_0_file_NEW);

	double *Vx = new double [dim_Wx];
	double *Vy = new double [dim_Wy];
	double *Vz = new double [dim_Wz];

	Velocity_Init(Vx, Vy, Vz, Nx, Ny, Nz);

	FILE * vel_x = fopen("velocity_x_3d.xls","wt");
	for (int i = 0; i < dim_Wx; i++)
		fprintf(vel_x, "%f \n", Vx[i]);
	fclose(vel_x);
	FILE * vel_y = fopen("velocity_y_3d.xls","wt");
	for (int i = 0; i < dim_Wy; i++)
		fprintf(vel_y, "%f \n", Vy[i]);
	fclose(vel_y);
	FILE * vel_z = fopen("velocity_z_3d.xls","wt");
	for (int i = 0; i < dim_Wz; i++)
		fprintf(vel_z, "%f \n", Vz[i]);
	fclose(vel_z);


	//for (int i = 0; i < dim_Wx; i++)
	//	Vx[i] = 0.0;

	//for (int i = 0; i < dim_Wy; i++)
	//	Vy[i] = 0.0;

	//for (int i = 0; i < dim_Wz; i++)
	//	Vz[i] = 0.0;


	
//	for (int i = 0; i < Nz+1; i++)
//		printf("Wz_n[%d] = %f \n",i,Wz_n[i]);
//	_getch();

	/*
	//	for (int j = 0; j < Nx*(Ny+1)*Nz; j++)
	for (int j = 0; j < (Ny+1); j++)
	{
	printf("Wy_n[%d] = %f \n",j,Wy_n[j]);
	//		fprintf(filee1,"%f \n",Wy_n[j]);
	//		fprintf(filee2,"%f \n",(-heat_conductivity(0,j,0)*exact_gradient206y(xpoints[0]+0.5*hx(0),ypoints[j], zpoints[0] + 0.5*hz(0), 0, mx, n) ));
	}
	_getch();

	for (int i = 0; i < (Nx+1); i++)
	printf("Wx_n[%d] = %f \n",i,Wx_n[i]);
	_getch();

	for (int i = 0; i < (Nz+1); i++)
	printf("Wz_n[%d] = %f \n",i,Wz_n[i]);
	_getch();
	*/
	////////для остальных компонент потока
	//**********************
	//**********************
	//**********************
	//конец инициализации w0

	//печать массивов Wx_0, Wy_0, Wz_0. первый столбец - посчитанные данные, второй - точные значения для
	//тестового решения num_sol 1111.
	switch(printcase)
	{
	case 0:
		break;
	case 1:
		{
			FILE* wz_check = fopen("Wz0_check_split3D.xls","wt");
			double add1 = 0;
			for (int i = 0; i < dim_Wz; i++)
			{
				k1 = i % (Nz+1);
				j1 = i / ((Nz+1)*Nx);
				i1 = (i - k1 - j1*(Nz+1)*Nx) / (Nz+1);

				add1 = -heat_conductivity(i1,j1,k1)*exact_gradientZ(num_sol, xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0);
				//switch(num_sol)
				//{
				//case 9006:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient9006z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0, mx, my, mz);
				//	break;
				//case 9066:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient9066z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0, mx, my, mz);
				//	break;
				//case 999:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient999z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0, mx, my, mz);
				//	break;
				//case 8:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient8z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0, A888, DEG888, Dir_axis);
				//	break;
				//case 888:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient888z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0, A888, DEG888);
				//	break;
				//case 888:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient888z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0, A888, DEG888);
				//	break;
				//case 500:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient500z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0,mx);
				//	break;
				//case 501:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient501z(xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0);
				//	break;
				//}

				fprintf(wz_check,"%15.15f \t %15.15f \t \t %15.15f \n",Wz_n[i], add1, Wz_n[i] - add1);
			}
			fclose(wz_check);

			FILE* wy_check = fopen("Wy0_check_split3D.xls","wt");
			for (int i = 0; i < dim_Wy; i++)
			{
				j1 = i % (Ny+1);
				i1 = i / ((Ny+1)*Nz);
				k1 = (i - j1 - i1*(Ny+1)*Nz) / (Ny+1);

				add1 = -heat_conductivity(i1,j1,k1)*exact_gradientY(num_sol, xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0);
				//switch(num_sol)
				//{
				//case 9006:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient9006y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0, mx, my, mz);
				//	break;
				//case 9066:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient9066y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0, mx, my, mz);
				//	break;
				//case 999:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient999y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0, mx, my, mz);
				//	break;
				//case 8:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient8y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0, A888, DEG888, Dir_axis);
				//	break;
				//case 888:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient888y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0, A888, DEG888);
				//	break;
				//case 500:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient500y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0,mx);
				//	break;
				//case 501:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient501y(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0);
				//	break;
				//}

				fprintf(wy_check,"%15.15f \t %15.15f \t \t %15.15f \n",Wy_n[i], add1, Wy_n[i] - add1);
			}
			fclose(wy_check);

			FILE* wx_check = fopen("Wx0_check_split3D.xls","wt");
			for (int i = 0; i < dim_Wx; i++)
			{
				i1 = i % (Nx+1);
				k1 = i / ((Nx+1)*Ny);
				j1 = (i - i1 - k1*(Nx+1)*Ny) / (Nx+1);

				add1 = -heat_conductivity(i1,j1,k1)*exact_gradientX(num_sol, xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0);

				//switch(num_sol)
				//{
				//case 9006:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient9006x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0, mx, my, mz);
				//	break;
				//case 9066:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient9066x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0, mx, my, mz);
				//	break;
				//case 999:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient999x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0, mx, my, mz);
				//	break;
				//case 8:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient8x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0, A888, DEG888, Dir_axis);
				//	break;
				//case 888:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient888x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0, A888, DEG888);
				//	break;	
				//case 500:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient500x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0, mx);
				//	break;	
				//case 501:
				//	add1 = -heat_conductivity(i1,j1,k1)*exact_gradient501x(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0);
				//	break;	
				//}
				
				fprintf(wx_check,"%15.15f \t %15.15f \t \t %15.15f \n",Wx_n[i], add1, Wx_n[i] - add1);
			}
			fclose(wx_check);
		}
		break;
	}

#ifdef SPECIALTEST
	for (int i = 0; i < Ny*Nz; i++)
	{
		for (int k = 0; k < Nx+1; k++)
		{
			i1 = k;
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			Wx_n[i*(Nx+1) + k] = Wx_0_exact[i*(Nx+1) + k] - Wx_n_Taver[i*(Nx+1) + k];
		}
	}
	for (int i=0; i<Nz*Nx; i++)
	{
		for (int k=0; k<Ny+1; k++)
		{
			i1 = i/Nz;
			j1 = k;
			k1 = i - (i/Nz)*Nz;

			Wy_n[i*(Ny+1) + k] = Wy_0_exact[i*(Ny+1) + k] - Wy_n_Taver[i*(Ny+1) + k];
		}
	}

	for (int i=0; i<Ny*Nx; i++)
	{
		for (int k=0; k<Nz+1; k++)
		{
			i1 = i/Ny;
			j1 = i - i1 * Ny;
			k1 = k;

			Wz_n[i*(Nz+1) + k1] = Wz_0_exact[i*(Nz+1) + k1] - Wz_n_Taver[i*(Nz+1) + k1];
		}
	}
#endif

	double eps0_max_w = 0.0;
	double eps0_max_T = 0.0;
	double eps0_l2_w = 0.0;
	double eps0_l2_T = 0.0;
	double eps0_relative_max_w = 0.0;
	double eps0_relative_max_T = 0.0;
	double eps0_relative_l2_w = 0.0;
	double eps0_relative_l2_T = 0.0;
	Accuracy_calculate(T, Wx_n, Wy_n, Wz_n, num_sol, -1, Nx, Ny, Nz, print_step, &eps0_max_w, &eps0_max_T, &eps0_l2_w, &eps0_l2_T, &eps0_relative_max_w, &eps0_relative_max_T, &eps0_relative_l2_w, &eps0_relative_l2_T);

	printf("After calling Accuracy_calculate \n");
	printf("eps0_max_w = %e \n",eps0_max_w);
	printf("eps0_l2_w = %e \n",eps0_l2_w);
	printf("eps0relative_max_w = %e \n",eps0_relative_max_w);
	printf("eps0relative_l2_2 = %e \n",eps0_relative_l2_w);

	fprintf(f1,"T = 0: \n");
	fprintf(f1,"eps0_max_T = %e \n", eps0_max_T);
	fprintf(f1,"eps0_l2_T = %e \n", eps0_l2_T);
	fprintf(f1,"eps0_max_w = %e \n", eps0_max_w);
	fprintf(f1,"eps0_l2_w = %e \n", eps0_l2_w);
	fprintf(f1,"releps0_max_T = %e \n", eps0_relative_max_T);
	fprintf(f1,"releps0_l2_T = %e \n", eps0_relative_l2_T);
	fprintf(f1,"releps0_max_w = %e \n", eps0_relative_max_w);
	fprintf(f1,"releps0_l2_w = %e \n", eps0_relative_l2_w);


				printf("Starting main loop...\n");

				///определение размеров массивов, с которыми имеем дело в основном цикле
				//double* delta_Wx_nplus1;
				//double* delta_Wy_nplus1;
				double* delta_Wz_nplus1;
				//double* delta_Wx_nplus23;
				double* delta_Wy_nplus23;
				//double* delta_Wz_nplus23;
				double* delta_Wx_nplus13;
				double* delta_Wy_nplus13;
				double* delta_Wz_nplus13;

				double* Wx_nplus1;
				double* Wy_nplus1;
				double* Wz_nplus1;
				//double* Wx_nplus23;
				//double* Wy_nplus23;
				//double* Wz_nplus23;
				//double* Wx_nplus13;
				double* Wy_nplus13;
				double* Wz_nplus13;
				double* T_nplus1;     //массив, хранящий значения температ. на n-ом слое
				Wx_nplus1 = new double [dim_Wx];
				Wy_nplus1 = new double [dim_Wy];
				Wz_nplus1 = new double [dim_Wz];

				double* Wx_nplus05;
				double* Wy_nplus05;
				double* Wz_nplus05;
				Wx_nplus05 = new double [dim_Wx];
				Wy_nplus05 = new double [dim_Wy];
				Wz_nplus05 = new double [dim_Wz];

				double * F = (double*)malloc(dim_T * sizeof(double));
				double * righthandX = (double*)malloc(dim_Wx * sizeof(double));
				double * righthandY = (double*)malloc(dim_Wy * sizeof(double));
				double * righthandZ = (double*)malloc(dim_Wz * sizeof(double));

				double * Temp_Wx = (double*)malloc(dim_Wx * sizeof(double));
				double * Temp_Wy = (double*)malloc(dim_Wy * sizeof(double));
				double * Temp_Wz = (double*)malloc(dim_Wz * sizeof(double));

				//Wx_nplus23 = new double [dim_Wx];
				//Wy_nplus23 = new double [dim_Wy];
				//Wz_nplus23 = new double [dim_Wz];

				//Wx_nplus13 = new double [dim_Wx];
				Wy_nplus13 = new double [dim_Wy];
				Wz_nplus13 = new double [dim_Wz];

				//delta_Wx_nplus1 = new double [dim_Wx];
				//delta_Wy_nplus1 = new double [dim_Wy];
				delta_Wz_nplus1 = new double [dim_Wz];

				//delta_Wx_nplus23 = new double [dim_Wx];
				delta_Wy_nplus23 = new double [dim_Wy];
				//delta_Wz_nplus23 = new double [dim_Wz];

				delta_Wx_nplus13 = new double [dim_Wx];
				delta_Wy_nplus13 = new double [dim_Wy];
				delta_Wz_nplus13 = new double [dim_Wz];
				T_nplus1 = new double [dim_T];
				//				arr1 = new double[Nx + 1];
				//				arr2 = new double[Nx + 1]; 
				//double eps_max_T = 0;			//(для тестов!)погрешность = sup (по слоям) max модуля разности (точного решения и посчитанного)
				//double eps_l2_T = 0;			//(для тестов!)погрешность = sup (по слоям) l2 нормы разности (точного решения и посчитанного)
				//double eps_max_w = 0;			
				//double eps_l2_w = 0;	
				double special_eps_l2_T = 0;
				double special_eps_l2_w = 0;
				double special_eps_max_T = 0;
				double special_eps_max_w = 0;

				double my_eps_max_T = 0;
				double my_eps_l2_T = 0;	
				double my_eps_max_w = 0;			
				double my_eps_l2_w = 0;	
				double my_releps_max_T = 0;
				double my_releps_l2_T = 0;	
				double my_releps_max_w = 0;			
				double my_releps_l2_w = 0;	

				double eps_max_T = 0;
				double eps_l2_T = 0;
				double eps_max_w = 0;			
				double eps_l2_w = 0;	
				double eps_relative_max_w = 0.0;
				double eps_relative_max_T = 0.0;
				double eps_relative_l2_w = 0.0;
				double eps_relative_l2_T = 0.0;
				double diffPr_max_w = 0.0;
				double diffPr_l2_w = 0.0;
				double diffrelPr_max_w = 0.0;
				double diffrelPr_l2_w = 0.0;
				double diffPr_max_T = 0.0;
				double diffPr_l2_T = 0.0;
				double diffrelPr_max_T = 0.0;
				double diffrelPr_l2_T = 0.0;

				double *Wx_0_error = (double*)malloc(dim_Wx*sizeof(double));
				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_error[i*(Nx+1) + k] = Wx_n[i*(Nx+1) + k] - (-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(num_sol,xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}
				double *Wy_0_error = (double*)malloc(dim_Wy*sizeof(double));
				for (int i=0; i<Nz*Nx; i++)
				{
					for (int k=0; k<Ny+1; k++)
					{
						i1 = i/Nz;
						j1 = k;
						k1 = i - (i/Nz)*Nz;

						Wy_0_error[i*(Ny+1) + k] = Wy_n[i*(Ny+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(num_sol,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}

				double *Wz_0_error = (double*)malloc(dim_Wz*sizeof(double));
				for (int i=0; i<Ny*Nx; i++)
				{
					for (int k=0; k<Nz+1; k++)
					{
						i1 = i/Ny;
						j1 = i - i1 * Ny;
						k1 = k;

						Wz_0_error[i*(Nz+1) + k1] = Wz_n[i*(Nz+1) + k1] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[j1]+0.5*hy(j1), zpoints[k1])*exact_gradientZ(num_sol,xpoints[i1]+0.5*hx(i1),ypoints[j1]+0.5*hy(j1), zpoints[k1], 0) );
					}
				}

#ifdef COMPARE_PROJECTORS
				double *Wx_0_diffPr = new double [dim_Wx];
				double *Wy_0_diffPr = new double [dim_Wy];
				double *Wz_0_diffPr = new double [dim_Wz];

				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_diffPr[i*(Nx+1) + k] = Wx_0_exact[i*(Nx+1) + k] - Wx_n_Taver[i*(Nx+1) + k];
					}
				}
				for (int i=0; i<Nz*Nx; i++)
				{
					for (int k=0; k<Ny+1; k++)
					{
						i1 = i/Nz;
						j1 = k;
						k1 = i - (i/Nz)*Nz;

						Wy_0_diffPr[i*(Ny+1) + k] = Wy_0_exact[i*(Ny+1) + k] - Wy_n_Taver[i*(Ny+1) + k];
					}
				}

				for (int i=0; i<Ny*Nx; i++)
				{
					for (int k=0; k<Nz+1; k++)
					{
						i1 = i/Ny;
						j1 = i - i1 * Ny;
						k1 = k;

						Wz_0_diffPr[i*(Nz+1) + k1] = Wz_0_exact[i*(Nz+1) + k1] - Wz_n_Taver[i*(Nz+1) + k1];
					}
				}
				printf("l2-norm of Wx_0_diff betw. projectors = %f \n", norm_l2(Wx_0_diffPr,Nx+1,Ny, Nz));
				printf("l2-norm of Wy_0_diff betw. projectors = %f \n", norm_l2(Wy_0_diffPr,Ny+1,Nz, Nx));
				printf("l2-norm of Wz_0_diff betw. projectors = %f \n", norm_l2(Wz_0_diffPr,Nz+1,Ny, Nx));
				printf("max-norm of Wx_0_diff betw. projectors = %f \n", norm_max(Wx_0_diffPr,Nx+1,Ny, Nz)); 
				printf("max-norm of Wy_0_diff betw. projectors = %f \n", norm_max(Wy_0_diffPr,Ny+1,Nz, Nx)); 
				printf("max-norm of Wz_0_diff betw. projectors = %f \n", norm_max(Wz_0_diffPr,Nz+1,Ny, Nx));

				printf ( "Computing new stability norm (for W0_error between Projectors) \n" );
				double newstabnormPr = 0.0;
				stability_norm3Dxyz ( Wx_0_diffPr, Wy_0_diffPr, Wz_0_diffPr, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnormPr);
				//stability_norm3Dyzx ( Wx_0_diffPr, Wy_0_diffPr, Wz_0_diffPr, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnormPr);
				printf ( "new stab norm for projectors diff= %f \n", newstabnormPr );

				FILE * xfile = fopen("Wx_0_diffPr_3d.xls","wt");
				for (int i = 0; i < dim_Wx; i++)
					fprintf(xfile,"%f \n", Wx_0_diffPr[i]);
				fclose(xfile);
				FILE * yfile = fopen("Wy_0_diffPr_3d.xls","wt");
				for (int i = 0; i < dim_Wy; i++)
					fprintf(yfile,"%f \n", Wy_0_diffPr[i]);
				fclose(yfile);
				FILE * zfile = fopen("Wz_0_diffPr_3d.xls","wt");
				for (int i = 0; i < dim_Wz; i++)
					fprintf(zfile,"%f \n", Wz_0_diffPr[i]);
				fclose(zfile);

				free (Wx_0_diffPr);
				free (Wy_0_diffPr);
				free (Wz_0_diffPr);
#endif
				//double l2_norm_stab = 0.0;
				//stability_normFull(Wx_0_exact,Wy_0_exact,Wz_0_exact,Nx,Ny,Nz,&l2_norm_stab);
				//printf("stability FULL l2-norm of W_0_exact = %f \n", l2_norm_stab);
				//stability_normFull(Wx_0_error,Wy_0_error,Wz_0_error,Nx,Ny,Nz,&l2_norm_stab);
				//printf("stability FULL l2-norm of W_0_error = %f \n", l2_norm_stab);
				////fprintf(f1,"stability FULL l2_norm (of error) = %e \n", l2_norm_stab);

				printf("l2-norm of Wx_0_error = %f \n", norm_l2(Wx_0_error,Nx+1,Ny, Nz));
				printf("l2-norm of Wy_0_error = %f \n", norm_l2(Wy_0_error,Ny+1,Nz, Nx));
				printf("l2-norm of Wz_0_error = %f \n", norm_l2(Wz_0_error,Nz+1,Ny, Nx));
				printf("max-norm of Wx_0_error = %f \n", norm_max(Wx_0_error,Nx+1,Ny, Nz)); 
				printf("max-norm of Wy_0_error = %f \n", norm_max(Wy_0_error,Ny+1,Nz, Nx)); 
				printf("max-norm of Wz_0_error = %f \n", norm_max(Wz_0_error,Nz+1,Ny, Nx)); 
#ifdef ERROR_FOR_L2PROJECTION
				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_error[i*(Nx+1) + k] = Wx_n[i*(Nx+1) + k] - Wx_n_Taver[i*(Nx+1) + k];
					}
				}
				for (int i=0; i<Nz*Nx; i++)
				{
					for (int k=0; k<Ny+1; k++)
					{
						i1 = i/Nz;
						j1 = k;
						k1 = i - (i/Nz)*Nz;

						Wy_0_error[i*(Ny+1) + k] = Wy_n[i*(Ny+1) + k] - Wy_n_Taver[i*(Ny+1) + k];
					}
				}

				for (int i=0; i<Ny*Nx; i++)
				{
					for (int k=0; k<Nz+1; k++)
					{
						i1 = i/Ny;
						j1 = i - i1 * Ny;
						k1 = k;

						Wz_0_error[i*(Nz+1) + k1] = Wz_n[i*(Nz+1) + k1] - Wz_n_Taver[i*(Nz+1) + k1];
					}
				}
#endif
				printf ( "Checking new stability norm for W0_exact \n" );
				double newstabnormm = 0.0;
				stability_norm3Dxyz ( Wx_0_exact, Wy_0_exact, Wz_0_exact, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnormm);
				//stability_norm3Dyzx ( Wx_0_exact, Wy_0_exact, Wz_0_exact, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnormm);
				printf ( "new stab normm (W0_exact) = %f \n\n", newstabnormm );

				//return 0;
				
#ifdef COMPARE_PROJECTORS
				printf ( "Checking new stability norm for W0_exact_project \n" );
				double newstabnormm_proj = 0.0;
				stability_norm3Dxyz ( Wx_n_Taver, Wy_n_Taver, Wz_n_Taver, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnormm_proj);
				//stability_norm3Dyzx ( Wx_n_Taver, Wy_n_Taver, Wz_n_Taver, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnormm_proj);
				printf ( "new stab normm (W0_exact_projected) = %f \n\n", newstabnormm_proj );
#endif

#ifdef SPECIALTEST
				printf ( "Computing new stability norm (for W0_error which is W0 for the SPECIALTEST) \n" );
				double newstabnorm = 0.0;
				stability_norm3Dxyz ( Wx_n, Wy_n, Wz_n, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnorm);
				//stability_norm3Dyzx ( Wx_n, Wy_n, Wz_n, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnorm);
				printf ( "new stab norm = %f \n", newstabnorm );
#else
				printf ( "Computing new stability norm (for W0_error) \n" );
				double newstabnorm = 0.0;
				stability_norm3Dxyz ( Wx_0_error, Wy_0_error, Wz_0_error, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnorm);
				//stability_norm3Dyzx ( Wx_0_error, Wy_0_error, Wz_0_error, Nx, Ny, Nz, tau, xCondition, yCondition, zCondition, &newstabnorm);
				printf ( "new stab norm = %f \n", newstabnorm );
#endif

				//return 0;

				//int scheme = 1;
				// 1 - DouglasGunn_Scheme_mainstep
				// 2 - Local1D_3D_Scheme_mainstep
				// 3 - PredictorCorrector3D_Scheme_mainstep
				// 4 - Uzawa3D_Scheme_mainstep

				switch(scheme)
				{
				case 1:
					printf ( "Douglas Gunn scheme \n" );
					break;
				case 2:
					printf ( "Local1d scheme \n" );
					break;
				case 3:
					printf ( "Predictor-corrector3d scheme \n" );
					break;
				case 4:
					printf ( "Uzawa3d scheme \n" );
					break;
				default:
					printf ( "Wrong value of scheme = %d - code will not work \n", scheme );
					return 0;
					break;
				}

//				unsigned int StartTime = GetTickCount();
				int count = 0;
				//N = 1;
				for ( int t = 0 ; t < N ; t++ )
				{
//					printf("Time step started t = %d\n",t);
//					_getch();

					/////////////////////////////////////////////////////////////////////////////////////////////
					/*
					//////////////////////////////////////////STEP 1/////////////////////////////////////////////
					*/
					//////////0///////////////////////////////////////////////////////////////////////////////////
					{;}
					/*
					FILE* deleteit1 = fopen("righthand_step1.xls","wt");
					if (deleteit1 == NULL)
					{
					printf("can not open file righthand_step1.txt");
					_getch();
					exit(3);
					}
					*/
					{;}
					switch(scheme)
					{
					case 1:
						DouglasGunn_Scheme_mainstep(T, T_nplus1, F, Temp_Wx, Temp_Wy, Temp_Wz,
							righthandX, Vx, Wx_n, Wx_nplus1, 
							righthandY, Vy, Wy_n, Wy_nplus13, Wy_nplus1, 
							righthandZ, Vz, Wz_n, Wz_nplus13, Wz_nplus1,
							t, tau, num_sol, xCondition, yCondition, zCondition, Nx, Ny, Nz, print_step, 
							&eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, 
							&eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T,
							&diffPr_max_w, &diffPr_l2_w, &diffrelPr_max_w, &diffrelPr_l2_w,
							&diffPr_max_T, &diffPr_l2_T, &diffrelPr_max_T, &diffrelPr_l2_T);
						break;
					case 2:
						Local1D_3D_Scheme_mainstep(T, T_nplus1, F, righthandX, Vx, Wx_n, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus1, righthandZ, Vz, Wz_n, Wz_nplus1, t, tau, num_sol, xCondition, yCondition, zCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
						break;
					case 3:
						PredictorCorrector3D_Scheme_mainstep(T, T_nplus1, F, Temp_Wx, Temp_Wy, Temp_Wz,
							righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, 
							righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 
							righthandZ, Vz, Wz_n, Wz_nplus05, Wz_nplus1,
							t, tau, num_sol, xCondition, yCondition, zCondition, Nx, Ny, Nz, print_step,
							&eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, 
							&eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T,
							&diffPr_max_w, &diffPr_l2_w, &diffrelPr_max_w, &diffrelPr_l2_w,
							&diffPr_max_T, &diffPr_l2_T, &diffrelPr_max_T, &diffrelPr_l2_T);
						break;
					case 4:
						Uzawa3D_Scheme_mainstep(T, T_nplus1, F, Temp_Wx, Temp_Wy, Temp_Wz,
							righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, 
							righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 
							righthandZ, Vz, Wz_n, Wz_nplus1, 
							t, tau, num_sol, xCondition, yCondition, zCondition, Nx, Ny, Nz, print_step, 
							&eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, 
							&eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
						break;
					default:
                        printf ( "Code should not get here in any case. \n");
						return -1;
						break;
					}
					
					{;}
				}

				printf("Final: \n");
				printf("Nx = %d \ntau = %f N = %d\n", Nx, tau, N);
				printf("eps_max_T = %e \n", eps_max_T);
				printf("eps_l2_T = %e \n", eps_l2_T);
				printf("eps_max_w = %e \n", eps_max_w);
				printf("eps_l2_w = %e \n", eps_l2_w);
				printf("special_eps_l2_T = %e \n", eps_relative_l2_T);
				printf("special_eps_l2_w = %e \n", eps_relative_l2_w);
				printf("special_eps_max_T = %e \n", eps_relative_max_T);
				printf("special_eps_max_w = %e \n", eps_relative_max_w);

#ifdef COMPARE_PROJECTORS
				printf("diffPr_max_T = %e \n", diffPr_max_T);
				printf("diffPr_l2_T = %e \n", diffPr_l2_T);
				printf("diffrelPr_max_T = %e \n", diffrelPr_max_T);
				printf("diffrelPr_l2_T = %e \n", diffrelPr_l2_T);

				printf("diffPr_max_w = %e \n", diffPr_max_w);
				printf("diffPr_l2_w = %e \n", diffPr_l2_w);
				printf("diffrelPr_max_w = %e \n", diffrelPr_max_w);
				printf("diffrelPr_l2_w = %e \n", diffrelPr_l2_w);
#endif

Exit:
				//Выдача теплового потока W , температуры T и погрешностей eps_max и eps_l2
				fprintf(f1,"Final: \n");
				fprintf(f1,"Nx = %d, Ny = %d, Nz = %d \ntau = %f N = %d \n", Nx, Ny, Nz, tau, N);
				fprintf(f1,"eps_max_T = %e \n", eps_max_T);
				fprintf(f1,"eps_l2_T = %e \n", eps_l2_T);
				fprintf(f1,"eps_max_w = %e \n", eps_max_w);
				fprintf(f1,"eps_l2_w = %e \n", eps_l2_w);
#ifdef COMPARE_PROJECTORS
				fprintf(f1,"diffPr_max_T = %e \n", diffPr_max_T);
				fprintf(f1,"diffPr_l2_T = %e \n", diffPr_l2_T);
				fprintf(f1,"diffrelPr_max_T = %e \n", diffrelPr_max_T);
				fprintf(f1,"diffrelPr_l2_T = %e \n", diffrelPr_l2_T);

				fprintf(f1,"diffPr_max_w = %e \n", diffPr_max_w);
				fprintf(f1,"diffPr_l2_w = %e \n", diffPr_l2_w);
				fprintf(f1,"diffrelPr_max_w = %e \n", diffrelPr_max_w);
				fprintf(f1,"diffrelPr_l2_w = %e \n", diffrelPr_l2_w);
#endif
				fprintf(f1,"T = 1: \n");
				fprintf(f1,"my_eps_max_T = %e \n", my_eps_max_T);
				fprintf(f1,"my_eps_l2_T = %e \n", my_eps_l2_T);
				fprintf(f1,"my_eps_max_w = %e \n", my_eps_max_w);
				fprintf(f1,"my_eps_l2_w = %e \n", my_eps_l2_w);
				fprintf(f1,"my_releps_max_T = %e \n", my_releps_max_T);
				fprintf(f1,"my_releps_l2_T = %e \n", my_releps_l2_T);
				fprintf(f1,"my_releps_max_w = %e \n", my_releps_max_w);
				fprintf(f1,"my_releps_l2_w = %e \n", my_releps_l2_w);

				fprintf(f1,"num_sol = %d, mx = %d, my = %d mz = %d\nexternal = %d \n", num_sol,mx,my,mz,external);
				fprintf(f1,"Tx0 = %f, Tx1 = %f\nTy0 = %f, Ty_1 = %f \nTz_0 = %f Tz_1 = %f\nwx_0 = %f, wx_1 = %f\nwy_0 = %f, wy_1 = %f \nwz_0 = %f wz_1 = %f \n", Tx_0, Tx_1, Ty_0, Ty_1, Tz_0,Tz_1,wx_0, wx_1, wy_0, wy_1, wz_0, wz_1);
				fprintf(f1,"special_eps_l2_T = %e \n", special_eps_l2_T);
				fprintf(f1,"special_eps_l2_w = %e \n", special_eps_l2_w);
				fprintf(f1,"special_eps_max_T = %e \n", special_eps_max_T);
				fprintf(f1,"special_eps_max_w = %e \n", special_eps_max_w);

//				unsigned int EndTime = GetTickCount();
//				unsigned int Time_ms = EndTime - StartTime;
//				fprintf(f1,"Time total (milliseconds) = %d \n",Time_ms);
//				fprintf(f1,"period of Time = [0,%f] \n",Time);

				////////////
				//Освобождение памяти из-под всех использованных динамических массивов
				free(Wx_n);
				free(Wy_n);
				free(Wz_n);
				free(Wx_nplus1);
				free(Wy_nplus1);
				free(Wz_nplus1);
				//free(Wx_nplus23);
				//free(Wy_nplus23);
				//free(Wz_nplus23);
				//free(Wx_nplus13);
				free(Wy_nplus13);
				free(Wz_nplus13);
				//free(delta_Wx_nplus1);
				//free(delta_Wy_nplus1);
				free(delta_Wz_nplus1);
				//free(delta_Wx_nplus23);
				free(delta_Wy_nplus23);
				//free(delta_Wz_nplus23);
				free(delta_Wx_nplus13);
				free(delta_Wy_nplus13);
				free(delta_Wz_nplus13);

//				free(delta_Wx_nplus05);
//				free(delta_Wy_nplus05);
				free(T);
				free(T_nplus1);
				free(betax);
				free(alfax);
				free(betay);
				free(alfay);
				free(betaz);
				free(alfaz);
				fclose(f1);

#ifdef COMPARE_PROJECTORS
				free (T_averaged);
				free (Wx_n_Taver);
				free (Wy_n_Taver);
				free (Wz_n_Taver);
#endif

				//////////
				return 0;
}


double righthand_my(int k)
{
	double ret_k = 0;
	/*
	//for step 1.1
	//к-ты матрицы А
	double a_0_old = hy(0)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double a_i_old = hy(k)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double a_iminus1_old = hy(k-1)*hx(0)/(heat_conductivity(0,k-1,0)*6.0);
	double a_nminus1_old = hy(Ny-1)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double b_0_old = 2*a_0_old;
	double b_i_old = 2*a_iminus1_old + 2*a_i_old;
	double b_n_old = 2*a_nminus1_old;

	//к-ты матрицы mx
	double m_0 = heat_capacity(0,0,0)*density(0,0,0)*hy(0)*hx(0);
	double m_i = heat_capacity(0,k,0)*density(0,k,0)*hy(k)*hx(0);
	double m_iminus1 = heat_capacity(0,k-1,0)*density(0,k-1,0)*hy(k-1)*hx(0);
	double m_nminus1 = heat_capacity(0,0,0)*density(0,0,0)*hy(0)*hx(0);

	//к-ты матрицы A + tau R(tr) = A + tau/2 * BM(-1)B(tr)
	//	double a_i = a_i_old - tau*(1.0/m_i);
	//	double a_iminus1 = a_iminus1_old - tau*(1.0/m_iminus1);
	//	double b_i = b_i_old + tau*( 1.0/m_i + 1.0/m_iminus1 );

	double a_i = a_i_old ;
	double a_0 = a_0_old ;
	double a_nminus1 = a_nminus1_old ;
	double a_iminus1 = a_iminus1_old ;
	double b_0 = b_0_old ;
	double b_i = b_i_old ;
	double b_n = b_n_old ;
	*/
	/*
	//for step 1.2
	//к-ты матрицы А
	double a_i_old = hx(k)*hy(0)/(heat_conductivity(k,0,0)*6.0);
	double a_iminus1_old = hx(k-1)*hy(0)/(heat_conductivity(k-1,0,0)*6.0);
	double b_i_old = 2*a_iminus1_old + 2*a_i_old;

	//к-ты матрицы mx
	double m_i = heat_capacity(k,0,0)*density(k,0,0)*hx(k)*hy(0);
	double m_iminus1 = heat_capacity(k-1,0,0)*density(k-1,0,0)*hx(k-1)*hy(0);

	//к-ты матрицы Ax + tau/2 * BxM(-1)Bx(tr)
	double a_i = a_i_old - 0.5*tau*hy(0)*hy(0)*(1.0/m_i);
	double a_iminus1 = a_iminus1_old - 0.5*tau*hy(0)*hy(0)*(1.0/m_iminus1);
	double b_i = b_i_old + 0.5*tau*hy(0)*hy(0)*( 1.0/m_i + 1.0/m_iminus1 );
	*/
	/*
	//for step 2.1
	//к-ты матрицы А
	double a_i_old = hy(k)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double a_iminus1_old = hy(k-1)*hx(0)/(heat_conductivity(0,k-1,0)*6.0);
	double b_i_old = 2*a_iminus1_old + 2*a_i_old;

	//к-ты матрицы mx
	double m_i = heat_capacity(0,k,0)*density(0,k,0)*hy(k)*hx(0);
	double m_iminus1 = heat_capacity(0,k-1,0)*density(0,k-1,0)*hy(k-1)*hx(0);

	double a_i = a_i_old ;
	double a_iminus1 = a_iminus1_old ;
	double b_i = b_i_old ;
	*/
	//for step 2.2
	//к-ты матрицы А
	double a_0_old = hy(0)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double a_i_old = hy(k)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double a_iminus1_old = hy(k-1)*hx(0)/(heat_conductivity(0,k-1,0)*6.0);
	double a_nminus1_old = hy(Ny-1)*hx(0)/(heat_conductivity(0,k,0)*6.0);
	double b_0_old = 2*a_0_old;
	double b_i_old = 2*a_iminus1_old + 2*a_i_old;
	double b_n_old = 2*a_nminus1_old;

	//к-ты матрицы mx
	double m_0 = heat_capacity(0,0,0)*density(0,0,0)*hy(0)*hx(0);
	double m_i = heat_capacity(0,k,0)*density(0,k,0)*hy(k)*hx(0);
	double m_iminus1 = heat_capacity(0,k-1,0)*density(0,k-1,0)*hy(k-1)*hx(0);
	double m_nminus1 = heat_capacity(0,0,0)*density(0,0,0)*hy(0)*hx(0);

	//к-ты матрицы A + tau/2 * BM(-1)B(tr)
	double a_0 = a_0_old - 0.5*tau*hx(0)*hx(0)*(1.0/m_0);
	double a_i = a_i_old - 0.5*tau*hx(0)*hx(0)*(1.0/m_i);
	double a_iminus1 = a_iminus1_old - 0.5*tau*hx(0)*hx(0)*(1.0/m_iminus1);
	double a_nminus1 = a_nminus1_old - 0.5*tau*hx(0)*hx(0)*(1.0/m_nminus1);
	double b_0 = b_0_old + 0.5*tau*hx(0)*hx(0)*( 1.0/m_0 );
	double b_i = b_i_old + 0.5*tau*hx(0)*hx(0)*( 1.0/m_i + 1.0/m_iminus1 );
	double b_n = b_n_old + 0.5*tau*hx(0)*hx(0)*( 1.0/m_nminus1 );


	double delta_w_exact_iminus1 = k-1;
	double delta_w_exact_i = k;
	double delta_w_exact_iplus1 = k+1;
	{;}
	if (k==0)
		ret_k = b_0*delta_w_exact_i + a_0*delta_w_exact_iplus1;
	if (k==Nx)
		ret_k = a_nminus1*delta_w_exact_iminus1 + b_n*delta_w_exact_i;
	if (k==1)
		ret_k = a_0*delta_w_exact_iminus1 + b_i*delta_w_exact_i + a_i*delta_w_exact_iplus1;
	if (k==Nx-1)
		ret_k = a_iminus1*delta_w_exact_iminus1 + b_i*delta_w_exact_i + a_nminus1*delta_w_exact_iplus1;
	if ((k>1)&&(k<Nx-1))
		ret_k = a_iminus1*delta_w_exact_iminus1 + b_i*delta_w_exact_i + a_i*delta_w_exact_iplus1;
	return ret_k;
}



void allocation()
{
	int i = 0, j = 0, k = 0;
	xpoints = (double*)malloc(sizeof(double)*(Nx+1));
	xpoints[0] = 0;
	for (i = 1; i <= Nx; i++)
	{
		xpoints[i] = xpoints[i-1] + hx(i-1);
	}
	ypoints = (double*)malloc(sizeof(double)*(Ny+1));
	ypoints[0] = 0;
	for (j = 1; j <=Ny; j++)
	{
		ypoints[j] = ypoints[j-1] + hy(j-1);
	}
	zpoints = (double*)malloc(sizeof(double)*(Nz+1));
	zpoints[0] = 0;
	for (k = 1; k <= Nz; k++)
	{
		zpoints[k] = zpoints[k-1] + hz(k-1);
	}
}



void stability_normFull(double* wx, double* wy, double * wz, int Nx, int Ny, int Nz, double *l2_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny * Nz;
	int dimWy = Nx * (Ny + 1) * Nz;
	int dimWz = Nx * Ny * (Nz + 1);
	int dimT = Nx * Ny * Nz;

	//printf("Calculating ByB_tr begins ... \n");
	double * Btr_w0 = (double*)malloc(dimT * sizeof(double));

	double m_i;
	for (int k = 0; k < Nz; k++ )
	{
		for (int j = 0; j < Ny; j++ )
		{
			for (int i = 0; i < Nx; i++ )
			{
				m_i = heat_capacity(i,j,k)*density(i,j,k)*hy(j)*hx(i)*hz(k);
				Btr_w0[k * Nx * Ny + j * Nx + i] = 0.0;
				Btr_w0[k * Nx * Ny + j * Nx + i] += hy(j) * hz(k) * ( wx[i + 1 + j*(Nx+1) + k*(Nx+1)*Ny] - wx[i + j*(Nx+1) + k*(Nx+1)*Ny] );
				Btr_w0[k * Nx * Ny + j * Nx + i] += hx(i) * hz(k) * ( wy[i*(Ny+1)*Nz + k * (Ny+1) + j + 1] - wy[i*(Ny+1)*Nz + k * (Ny+1) + j] );
				Btr_w0[k * Nx * Ny + j * Nx + i] += hx(i) * hy(j) * ( wz[k + 1 + j*(Nz+1) + i*(Nz+1)*Ny] - wz[k + j*(Nz+1) + i*(Nz+1)*Ny] );
				Btr_w0[k * Nx * Ny + j * Nx + i] *= (1.0 / m_i );
			}
		}
	}

	double temp1, temp2;
	temp1 = norm_l2(Btr_w0, Nx, Ny, Nz);
	temp2 = norm_max(Btr_w0, Nx, Ny, Nz);
	printf("norm_l2 for Btr_w0 = %f \n", temp1);
	printf("norm_max for Btr_w0 = %f \n\n", temp2);

	double *Ayminus1_ByB_tr = (double*)malloc(dimWy * sizeof(double));

	HeatFlux_Wy0_Init( Btr_w0, Ayminus1_ByB_tr, yCondition, num_sol, Nx, Ny, Nz);

	temp1 = norm_l2(Ayminus1_ByB_tr,Ny + 1, Nz, Nx);
	temp2 = norm_max(Ayminus1_ByB_tr,Ny + 1, Nz, Nx);
	printf("norm_l2 for Aminus1_ByB_tr = %f \n", temp1);
	printf("norm_max for Aminus1_ByB_tr = %f \n\n", temp2);

	//printf("Calculating By_tr_Aminus1_ByB_tr w0 begins ... \n");
	//printf("Final right \n");
	double * MhalfLambdaY_Btr_w0 = (double*)malloc(dimT * sizeof(double)); // = By_tr * Aminus1_By_B_tr, но надо M

	for (int k = 0; k < Nz; k++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
				MhalfLambdaY_Btr_w0 [ k * Nx * Ny + j * Nx + i ] = hx(i) * hz(k) * ( Ayminus1_ByB_tr[j + 1 + i*(Ny+1)*Nz + k * (Ny+1)] - Ayminus1_ByB_tr[j + i*(Ny+1)*Nz + k * (Ny+1)] );
			}
		}
	}

	temp1 = norm_l2(MhalfLambdaY_Btr_w0, Nx, Ny, Nz);
	printf("norm_l2 for MhalfLambdaY_Btr_w0 = %f \n", temp1);
	temp2 = norm_max(MhalfLambdaY_Btr_w0, Nx, Ny, Nz);
	printf("norm_max for MhalfLambdaY_Btr_w0 = %f \n\n", temp2);

	double * finalRightpart = (double*)malloc(dimT * sizeof(double));  // = Lambda_y * Btr * w0
	for (int k = 0; k < Nz; k++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
				finalRightpart [ k * Nx * Ny + j * Nx + i ] = MhalfLambdaY_Btr_w0 [ k * Nx * Ny + j * Nx + i ] ; 
			}
		}
	}

	double * MminushalfLambdaY_Btr_w0 = (double*)malloc(dimT * sizeof(double)); // = M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	for (int k = 0; k < Nz; k++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
				m_i = heat_capacity(i,j,k)*density(i,j,k)*hy(j)*hx(i)*hz(k);
				MminushalfLambdaY_Btr_w0 [ k * Nx * Ny + j * Nx + i ] = (1.0 / m_i ) * MhalfLambdaY_Btr_w0 [ k * Nx * Ny + j * Nx + i  ];
			}
		}
	}

	temp1 = norm_l2(MminushalfLambdaY_Btr_w0, Nx, Ny, Nz);
	printf("norm_l2 for MminushalfLambdaY_Btr_w0 = %f \n", temp1);
	temp2 = norm_max(MminushalfLambdaY_Btr_w0, Nx, Ny, Nz);
	printf("norm_max for MminushalfLambdaY_Btr_w0 = %f \n\n", temp2);

	double * Ayminus1_xvost_w0 = (double*)malloc(dimWy * sizeof(double)); // = Ay(-1) * By M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	HeatFlux_Wy0_Init(MminushalfLambdaY_Btr_w0, Ayminus1_xvost_w0, yCondition, num_sol, Nx, Ny, Nz);

	temp1 = norm_l2(Ayminus1_xvost_w0, Ny + 1, Nz, Nx);
	printf("norm_l2 for Ayminus1_xvost_w0 = %f \n", temp1);
	temp2 = norm_max(Ayminus1_xvost_w0, Ny + 1, Nz, Nx);
	printf("norm_max for Ayminus1_xvost_w0 = %f \n\n", temp2);

	double * Mminus1Bytr_Ayminus1_xvost_w0 = (double*)malloc(dimT * sizeof(double)); // = M(-1) By_tr * Ay(-1) * By M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	for (int k = 0; k < Nz; k++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
				m_i = heat_capacity(i,j,k)*density(i,j,k)*hy(j)*hx(i)*hz(k);

				Mminus1Bytr_Ayminus1_xvost_w0 [ k * Nx * Ny + j * Nx + i ] = hx(i) * hz(k) * ( 1.0 / m_i ) * ( Ayminus1_xvost_w0[j + 1 + i*(Ny+1)*Nz + k * (Ny+1)] - Ayminus1_xvost_w0[j + i*(Ny+1)*Nz + k * (Ny+1)] );
			}
		}
	}

	temp1 = norm_l2(Mminus1Bytr_Ayminus1_xvost_w0, Nx, Ny, Nz);
	printf("norm_l2 for Mminus1Bytr_Ayminus1_xvost_w0 = %f \n", temp1);
	temp2 = norm_max(Mminus1Bytr_Ayminus1_xvost_w0, Nx, Ny, Nz);
	printf("norm_max for Mminus1Bytr_Ayminus1_xvost_w0 = %f \n\n", temp2);

	//printf("Final left \n");
	double * finalLeftpart = (double*)malloc(dimT * sizeof(double));  // = Lambda * Lambda_y * Btr * w0
	for (int k = 0; k < Nz; k++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
			{
				finalLeftpart [ k * Nx * Ny + j * Nx + i  ] = Mminus1Bytr_Ayminus1_xvost_w0 [ k * Nx * Ny + j * Nx + i  ]; 
			}
		}
	}

//	temp1 = scal_l2 ( finalLeftpart, MhalfLambdaY_Btr_w0 , Nx, Ny );
//	printf("scalnorm_l2 = %f \n", temp1);

	l2_norm = sqrt ( scal_l2 ( finalLeftpart, finalRightpart, Nx, Ny, Nz ) / (Nx * Ny * Nz) );
	printf ( "full stability norm = %f \n", l2_norm );

	free(Ayminus1_ByB_tr);
	free(Btr_w0);
	//free(Mminus1Bxtr_Axminus1_xvost_w0);
	free(MhalfLambdaY_Btr_w0);
	//free(Axminus1_xvost_w0);
	free(MminushalfLambdaY_Btr_w0);
	free(finalRightpart);
	free(finalLeftpart);
//	free(BxBy_tr);

	*l2_norm_pt = l2_norm;
}

double norm_l2(double * massiv, int N, int M, int K)
{
	double norm = 0.0;
	int index = 0;
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			for (int k = 0; k < K; k++)
			{
				norm += massiv[index]*massiv[index];
				index += 1;
			}
		}
	}
	norm = sqrt(norm / (N * M * K));
	return norm;
}

double scal_l2(double * massiv1, double * massiv2, int N, int M, int K)
{
	double out = 0.0;
	int index = 0;
	for (int k = 0; k < K; k++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int i = 0; i < M; i++)
			{
				out += massiv1[index]*massiv2[index];
				index += 1;
			}
		}
	}
	return out;
}

double norm_max(double * massiv, int N, int M, int K)
{
	double norm_max = 0.0;
	int index = 0;
	int i_max = 0;
	int j_max = 0;
	int k_max = 0;
	for (int k = 0; k < K; k++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int i = 0; i < M; i++)
			{
				if (fabs(massiv[index]) > norm_max)
				{
					norm_max = fabs(massiv[index]);
					i_max = i;
					j_max = j;
					k_max = k;
				}
				index += 1;
			}
		}
	}
	//printf("i_max = %d \n", i_max);
	//printf("j_max = %d \n", j_max);
	//printf("k_max = %d \n", k_max);

	return norm_max;
}

int stability_norm3Dxyz ( double * Wx, double * Wy, double * Wz, int Nx, int Ny, int Nz, double tau, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition, double * res_pt)
{
	double res;
	double term1, term2, term3;

	int dimT = Nx * Ny * Nz;

	double * tempT = (double*) malloc (dimT * sizeof(double));
	double * tempT2 = (double*) malloc (dimT * sizeof(double));

	// 1. compute term1 = || W ||_H = || W ||_Hdivh
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);

	term1 = norm_l2 ( tempT, Nx, Ny, Nz ) * norm_l2 ( tempT, Nx, Ny, Nz );								// term1 = || tempT ||_L2h^2 
	term1 += norm_l2 ( Wx, Nx + 1, Ny, Nz) * norm_l2 ( Wx, Nx + 1, Ny, Nz) +							
		norm_l2 ( Wy, Nx, Ny+1, Nz ) * norm_l2 ( Wy, Nx, Ny+1, Nz ) +
		norm_l2 ( Wz, Nx, Ny, Nz+1 ) * norm_l2 ( Wz, Nx, Ny, Nz+1 );									// term1 = || tempT ||_L2h^2 + || W ||_L2hv^2
	term1 = sqrt (term1);																				// term1 = || W ||_Hdivh

	printf ( "term1 = %f \n", term1 );
	printf ( "norm_BtrW_L2 = %f \n", norm_l2 ( tempT, Nx, Ny, Nz ) );
	printf ( "norm_W_L2 = %f \n", sqrt(norm_l2 ( Wx, Nx + 1, Ny, Nz) * norm_l2 ( Wx, Nx + 1, Ny, Nz) + norm_l2 ( Wy, Nx, Ny+1, Nz ) * norm_l2 ( Wy, Nx, Ny+1, Nz ) +
		norm_l2 ( Wz, Nx, Ny, Nz+1 ) * norm_l2 ( Wz, Nx, Ny, Nz+1 )) );

	// 2. compute term2 = || C Btr W ||_Lambda
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);
	printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	compute_CTxyz ( tempT2, tempT, tau, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT2 = C tempT = C Btr W
	printf ( "norm_CBtrW = %f \n", norm_l2 (tempT2, Nx, Ny, Nz) );
	zero_init (tempT, dimT);
	compute_Lambda (tempT, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT (new) = Lambda * tempT2 = Lambda * C * Btr * W
	printf ( "norm_Lambda_CBtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	term2 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny, Nz ) / dimT );

	printf ( "term2 = %f \n", term2 );

	// 3. compute term3 = || Lambdaz Btr W ||_Lambda
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);
	//fprintdVec ( "interm3.txt", tempT, dimT ):
	printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	compute_Lambdaz ( tempT2, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT2 = Lambdaz tempT = Lambdaz Btr W
	printf ( "norm_Lambdaz_BtrW = %f \n", norm_l2 (tempT2, Nx, Ny, Nz) );

	zero_init (tempT, dimT);
	compute_Lambda (tempT, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT (new) = Lambda * tempT2 = Lambda * Lambdaz * Btr * W
	printf ( "norm_Lambda_Lambdaz_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	term3 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny, Nz ) / dimT );

	//// 2. compute term2 = || C Btr W ||_Lambdax
	//zero_init (tempT, dimT);
	//Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	//for ( int i = 0; i < dimT; i++ )
	//	tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);
	//printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	//compute_CT ( tempT2, tempT, tau, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT2 = C tempT = C Btr W
	//printf ( "norm_CBtrW = %f \n", norm_l2 (tempT2, Nx, Ny, Nz) );
	//zero_init (tempT, dimT);
	//compute_Lambday (tempT, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT (new) = Lambday * tempT2 = Lambday * C * Btr * W
	//printf ( "norm_Lambdax_CBtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	//term2 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny, Nz ) / dimT );

	//printf ( "term2 = %f \n", term2 );

	//// 3. compute term3 = || Lambdaz Btr W ||_Lambda
	//zero_init (tempT, dimT);
	//Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	//for ( int i = 0; i < dimT; i++ )
	//	tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);
	////fprintdVec ( "interm3.txt", tempT, dimT ):
	//printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	//compute_Lambdaz ( tempT2, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT2 = Lambdaz tempT = Lambdaz Btr W
	//printf ( "norm_Lambdaz_BtrW = %f \n", norm_l2 (tempT2, Nx, Ny, Nz) );

	//zero_init (tempT, dimT);
	//compute_Lambday (tempT, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT (new) = Lambday * tempT2 = Lambday * Lambdaz * Btr * W
	//printf ( "norm_Lambday_Lambdaz_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	//term3 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny, Nz ) / dimT );


	double * deleteit1 = (double*)malloc(dimT * sizeof(double));
	compute_Lambdax ( deleteit1, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );				// deleteit1 = Lambdax tempT = Lambdax Btr W
	printf ( "norm_Lambdax_BtrW = %f \n", norm_l2 (deleteit1, Nx, Ny, Nz) );
	free (deleteit1);
	double * deleteit2 = (double*)malloc(dimT * sizeof(double));
	compute_Lambday ( deleteit2, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );				// deleteit2 = Lambday tempT = Lambday Btr W
	printf ( "norm_Lambday_BtrW = %f \n", norm_l2 (deleteit2, Nx, Ny, Nz) );
	free (deleteit2);

	printf ( "term3 = %f \n", term3 );

	printf ( "term1 = %f term2 + term3 with tau ~ %f \n", term1, sqrt (pow(tau,4) * (term2 * term2 + term3 * term3)));

	res = sqrt ( term1 * term1 + pow(tau,4) * (term2 * term2 + term3 * term3) );

	*res_pt = res;

	free ( tempT );
	free ( tempT2 );

	return 0;
}


int stability_norm3Dyzx ( double * Wx, double * Wy, double * Wz, int Nx, int Ny, int Nz, double tau, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition, double * res_pt)
{
	double res;
	double term1, term2, term3;

	int dimT = Nx * Ny * Nz;

	double * tempT = (double*) malloc (dimT * sizeof(double));
	double * tempT2 = (double*) malloc (dimT * sizeof(double));

	// 1. compute term1 = || W ||_H = || W ||_Hdivh
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);

	term1 = norm_l2 ( tempT, Nx, Ny, Nz ) * norm_l2 ( tempT, Nx, Ny, Nz );								// term1 = || tempT ||_L2h^2 
	term1 += norm_l2 ( Wx, Nx + 1, Ny, Nz) * norm_l2 ( Wx, Nx + 1, Ny, Nz) +							
		norm_l2 ( Wy, Nx, Ny+1, Nz ) * norm_l2 ( Wy, Nx, Ny+1, Nz ) +
		norm_l2 ( Wz, Nx, Ny, Nz+1 ) * norm_l2 ( Wz, Nx, Ny, Nz+1 );									// term1 = || tempT ||_L2h^2 + || W ||_L2hv^2
	term1 = sqrt (term1);																				// term1 = || W ||_Hdivh

	printf ( "term1 = %f \n", term1 );
	printf ( "norm_BtrW_L2 = %f \n", norm_l2 ( tempT, Nx, Ny, Nz ) );
	printf ( "norm_W_L2 = %f \n", sqrt(norm_l2 ( Wx, Nx + 1, Ny, Nz) * norm_l2 ( Wx, Nx + 1, Ny, Nz) + norm_l2 ( Wy, Nx, Ny+1, Nz ) * norm_l2 ( Wy, Nx, Ny+1, Nz ) +
		norm_l2 ( Wz, Nx, Ny, Nz+1 ) * norm_l2 ( Wz, Nx, Ny, Nz+1 )) );

	// 2. compute term2 = || C Btr W ||_Lambda
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);
	printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	compute_CTyzx ( tempT2, tempT, tau, Nx, Ny, Nz, xCondition, yCondition, zCondition );				// tempT2 = Cyzx tempT = C Btr W
	printf ( "norm_CBtrW = %f \n", norm_l2 (tempT2, Nx, Ny, Nz) );
	zero_init (tempT, dimT);
	compute_Lambday (tempT, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT (new) = Lambday * tempT2 = Lambday * Cyzx * Btr * W
	printf ( "norm_Lambday_CBtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	term2 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny, Nz ) / dimT );

	printf ( "term2 = %f \n", term2 );

	// 3. compute term3 = || Lambdax Btr W ||_Lambdaz
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Wz, zCondition, Nx, Ny, Nz);		// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0) / hz(0);
	//fprintdVec ( "interm3.txt", tempT, dimT ):
	printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	compute_Lambdax ( tempT2, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT2 = Lambdax tempT = Lambdax Btr W
	printf ( "norm_Lambdax_BtrW = %f \n", norm_l2 (tempT2, Nx, Ny, Nz) );

	zero_init (tempT, dimT);
	compute_Lambdaz (tempT, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition );					// tempT (new) = Lambdaz * tempT2 = Lambdaz * Lambdax * Btr * W
	printf ( "norm_Lambdaz_Lambdaz_BtrW = %f \n", norm_l2 (tempT, Nx, Ny, Nz) );
	term3 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny, Nz ) / dimT );

	double * deleteit1 = (double*)malloc(dimT * sizeof(double));
	compute_Lambday ( deleteit1, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );				// deleteit1 = Lambday tempT = Lambday Btr W
	printf ( "norm_Lambday_BtrW = %f \n", norm_l2 (deleteit1, Nx, Ny, Nz) );
	free (deleteit1);
	double * deleteit2 = (double*)malloc(dimT * sizeof(double));
	compute_Lambdaz ( deleteit2, tempT, Nx, Ny, Nz, xCondition, yCondition, zCondition );				// deleteit2 = Lambdaz tempT = Lambdaz Btr W
	printf ( "norm_Lambdaz_BtrW = %f \n", norm_l2 (deleteit2, Nx, Ny, Nz) );
	free (deleteit2);

	printf ( "term3 = %f \n", term3 );

	printf ( "term1 = %f term2 + term3 with tau ~ %f \n", term1, sqrt (pow(tau,4) * (term2 * term2 + term3 * term3)));

	res = sqrt ( term1 * term1 + pow(tau,4) * (term2 * term2 + term3 * term3) );

	*res_pt = res;

	free ( tempT );
	free ( tempT2 );

	return 0;
}
int zero_init ( double * vec, int dim )
{
	if ( vec == NULL ) 
	{
		printf ( "Zero pointer in zero_init() \n" );
		return -1;
	}
	else
	{
		for ( int i = 0; i < dim; i++ )
			vec[i] = 0.0;
		return 0;
	}
}

// computes C * T, where C = Lambda_y + Lambda_x + tau/2 Lambda_y * Lambda_z
//// computes C * T, where C = Lambda_x + Lambda_z + tau/2 Lambda_z * Lambda_x
int compute_CTxyz ( double * CT, double * T, double tau, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition )
{
	int dimT = Nx * Ny * Nz;

	double * tempT1 = (double*)malloc(dimT * sizeof(double));
	double * tempT2 = (double*)malloc(dimT * sizeof(double));
	double * tempT3 = (double*)malloc(dimT * sizeof(double));

	// 1.
	compute_Lambday ( tempT1, T, Nx, Ny, Nz, xCondition, yCondition, zCondition);			// tempT1 = Lambday * T

	// 2.
	compute_Lambdaz ( tempT2, T, Nx, Ny, Nz, xCondition, yCondition, zCondition);			// tempT2 = Lambdaz * T

	// 3.
	compute_Lambday ( tempT3, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition);		// tempT3 = Lambday * tempT2 = Lambday * Lambdaz * T

	//// 1.
	//compute_Lambdaz ( tempT1, T, Nx, Ny, Nz, xCondition, yCondition, zCondition);			// tempT1 = Lambdaz * T

	//// 2.
	//compute_Lambdax ( tempT2, T, Nx, Ny, Nz, xCondition, yCondition, zCondition);			// tempT2 = Lambdax * T

	//// 3.
	//compute_Lambdaz ( tempT3, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition);		// tempT3 = Lambdaz * tempT2 = Lambdaz * Lambdax * T

	for ( int i = 0; i < dimT; i++ )
		CT[i] = tempT1[i] + tempT2[i] + tau * 0.5 * tempT3[i];

	free (tempT1);
	free (tempT2);
	free (tempT3);

	return 0;
}


// computes C * T, where C = Lambda_x + Lambda_z + tau/2 Lambda_z * Lambda_x
int compute_CTyzx ( double * CT, double * T, double tau, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition )
{
	int dimT = Nx * Ny * Nz;

	double * tempT1 = (double*)malloc(dimT * sizeof(double));
	double * tempT2 = (double*)malloc(dimT * sizeof(double));
	double * tempT3 = (double*)malloc(dimT * sizeof(double));

	// 1.
	compute_Lambdaz ( tempT1, T, Nx, Ny, Nz, xCondition, yCondition, zCondition);			// tempT1 = Lambdaz * T

	// 2.
	compute_Lambdax ( tempT2, T, Nx, Ny, Nz, xCondition, yCondition, zCondition);			// tempT2 = Lambdax * T

	// 3.
	compute_Lambdaz ( tempT3, tempT2, Nx, Ny, Nz, xCondition, yCondition, zCondition);		// tempT3 = Lambdaz * tempT2 = Lambdaz * Lambdax * T

	for ( int i = 0; i < dimT; i++ )
		CT[i] = tempT1[i] + tempT2[i] + tau * 0.5 * tempT3[i];

	free (tempT1);
	free (tempT2);
	free (tempT3);

	return 0;
}

// computes out = Lambdax * in
int	compute_Lambdax ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition)
{
	int dimWx = (Nx + 1) * Ny * Nz;
	int dimWy = Nx * (Ny + 1) * Nz;
	int dimWz = Nx * Ny * (Nz + 1);
	int dimT = Nx * Ny * Nz;

	// 1. tempWx = Bx * in
	double * tempWx = (double*)malloc(dimWx * sizeof(double));
	zero_init (tempWx, dimWx);

	BxMminus1_F_Add( tempWx, in, xCondition, Nx, Ny, Nz);					// tempWx1 = Bx * in
//	printf ( "norm_1x = %e \n", norm_l2 ( tempWx, Nx + 1, Ny, Nz ) );

	// 2. tempWx2 = Ax^(-1) * Bx * in
	double * tempWx2 = (double*)malloc(dimWx * sizeof(double));		
	zero_init (tempWx2, dimWx);

	Invert_Ax( tempWx2, tempWx, xCondition, Nx, Ny, Nz);					// tempWx2 = Ax^(-1) * tempWx1 = Ax^(-1) * Bx * in
	for ( int i = 0; i < dimWx; i++ )
		tempWx2[i] *= hx(0);
//	printf ( "norm_2x = %e \n", norm_l2 ( tempWx2, Nx + 1, Ny, Nz ) );

	// 3. out = Bxtr * Ax^(-1) * Bx * in
	zero_init (out, dimT);
	Bxtr_Wx_Add( out, 1.0, tempWx2, xCondition, Nx, Ny, Nz);

//	printf ( "norm_3x = %e \n", norm_l2 ( out, Nx, Ny, Nz ) );

	for ( int i = 0; i < dimT; i++ )
		out[i] *= 1.0 / hx(0) / hy(0) / hz(0);

	free (tempWx);
	free (tempWx2);

	return 0;
}



// computes out = Lambday * in
int	compute_Lambday ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition)
{
	int dimWx = (Nx + 1) * Ny * Nz;
	int dimWy = Nx * (Ny + 1) * Nz;
	int dimWz = Nx * Ny * (Nz + 1);
	int dimT = Nx * Ny * Nz;

	// 1. tempWy = By * in
	double * tempWy = (double*)malloc(dimWy * sizeof(double));
	zero_init (tempWy, dimWy);

	ByMminus1_F_Add( tempWy, in, yCondition, Nx, Ny, Nz);					// tempWy1 = By * in

//	printf ( "norm_1y = %e \n", norm_l2 ( tempWy, Nx, Ny + 1, Nz ) );
	// 2. tempWy2 = Ay^(-1) * By * in
	double * tempWy2 = (double*)malloc(dimWy * sizeof(double));		
	zero_init (tempWy2, dimWy);

	Invert_Ay( tempWy2, tempWy, yCondition, Nx, Ny, Nz);					// tempWy2 = Ay^(-1) * tempWy1 = Ay^(-1) * By * in
	for ( int i = 0; i < dimWy; i++ )
		tempWy2[i] *= hy(0);
//	printf ( "norm_2y = %e \n", norm_l2 ( tempWy2, Nx, Ny + 1, Nz ) );

	// 3. out = Bytr * Ay^(-1) * By * in
	zero_init (out, dimT);
	Bytr_Wy_Add( out, 1.0, tempWy2, yCondition, Nx, Ny, Nz);
//	printf ( "norm_3y = %e \n", norm_l2 ( out, Nx, Ny, Nz ) );
	for ( int i = 0; i < dimT; i++ )
		out[i] *= 1.0 / hx(0) / hy(0) / hz(0);

	free (tempWy);
	free (tempWy2);

	return 0;
}

// computes out = Lambdaz * in
int	compute_Lambdaz ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition)
{
	int dimWx = (Nx + 1) * Ny * Nz;
	int dimWy = Nx * (Ny + 1) * Nz;
	int dimWz = Nx * Ny * (Nz + 1);
	int dimT = Nx * Ny * Nz;

	// 1. tempWz = Bz * in
	double * tempWz = (double*)malloc(dimWz * sizeof(double));
	zero_init (tempWz, dimWz);

	BzMminus1_F_Add( tempWz, in, zCondition, Nx, Ny, Nz);					// tempWz1 = Bz * in

//	printf ( "norm_1z = %e \n", norm_l2 ( tempWz, Nx, Ny, Nz + 1 ) );

	// 2. tempWz2 = Az^(-1) * Bz * in
	double * tempWz2 = (double*)malloc(dimWz * sizeof(double));		
	zero_init (tempWz2, dimWz);

	Invert_Az( tempWz2, tempWz, zCondition, Nx, Ny, Nz);					// tempWz2 = Az^(-1) * tempWz1 = Az^(-1) * Bz * in
	for ( int i = 0; i < dimWz; i++ )
		tempWz2[i] *= hz(0);
//	printf ( "norm_2z = %e \n", norm_l2 ( tempWz2, Nx, Ny, Nz + 1 ) );

	// 3. out = Bztr * Az^(-1) * Bz * in
	zero_init (out, dimT);
	Bztr_Wz_Add( out, 1.0, tempWz2, zCondition, Nx, Ny, Nz);
//	printf ( "norm_3z = %e \n", norm_l2 ( out, Nx, Ny, Nz ) );
	for ( int i = 0; i < dimT; i++ )
		out[i] *= 1.0 / hx(0) / hy(0) / hz(0);

	free (tempWz);
	free (tempWz2);

	return 0;
}

int	compute_Lambda ( double * out, double * in, int Nx, int Ny, int Nz, BNDCNDS xCondition, BNDCNDS yCondition, BNDCNDS zCondition)
{
	int dimT = Nx * Ny * Nz;

	double * tempT1 = (double*)malloc(dimT * sizeof(double));
	double * tempT2 = (double*)malloc(dimT * sizeof(double));
	double * tempT3 = (double*)malloc(dimT * sizeof(double));

	compute_Lambdax(tempT1, in, Nx, Ny, Nz, xCondition, yCondition, zCondition);
	compute_Lambday(tempT2, in, Nx, Ny, Nz, xCondition, yCondition, zCondition);
	compute_Lambdaz(tempT3, in, Nx, Ny, Nz, xCondition, yCondition, zCondition);

	for ( int i = 0; i < dimT; i++ )
		out[i] = tempT1[i] + tempT2[i] + tempT3[i];

	free (tempT1);
	free (tempT2);
	free (tempT3);

	return 0.0;
}


int Invert_Ax(double * sol, double * rhand, boundaryConditions xCondition, int Mx, int My, int Mz)
{
	double a_0_old, a_i_old, a_iminus1_old, a_nminus1_old;
	double b_0_old, b_i_old, b_n_old = 0;
	double temper_0, temper_k, temper_n = 0;

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

			betax[0] = 0 ;
			betax[1] = rhand[i*(Mx+1) + 0];
			alfax[0] = 0 ;
			alfax[1] = 0 ;
			for ( int k = 1 ; k < Mx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = rhand[i*(Mx+1) + k];

				if (k==Mx-1)
				{
					alfax[k + 1] = 0;
					betax[k + 1] = ( ( temper_k - a_i_old*rhand[i*(Mx+1) + Mx] ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( ( temper_k - a_iminus1_old*rhand[i*(Mx+1) + 0] ) ) / (  b_i_old );
				}
				else
				{
					alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
					betax[k + 1] = ( ( temper_k ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
				}

			}

			sol[i*(Mx+1) + Mx ] = rhand[i*(Mx+1) + Mx] ;
			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				sol[i*(Mx+1) + Mx - j] = sol[i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 - j];
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

			//к-ты матрицы А
			a_0_old = hx(0)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;

			temper_0 = rhand[i*(Mx+1) + 0];
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
				temper_k = rhand[i*(Mx+1) + k];

				alfax[k + 1] = (-1)*a_i_old/(alfax[k]*a_iminus1_old + b_i_old);
				betax[k + 1] = ( ( temper_k ) - betax[k]*a_iminus1_old ) / ( alfax[k]*a_iminus1_old + b_i_old );
					
			}
			a_nminus1_old = hx(Mx-1)/(heat_conductivity(Mx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			temper_n = rhand[i*(Mx+1) + Mx];

			sol[i*(Mx+1) + Mx] = ( temper_n  -  betax[Mx]*a_nminus1_old ) / (alfax[Mx]*a_nminus1_old + b_n_old);

			for ( int j = 1 ; j < Mx + 1 ; j++ )
			{
				sol[i*(Mx+1) + Mx - j] = sol [i*(Mx+1) + Mx + 1 - j] * alfax[Mx + 1 - j] + betax[Mx + 1 -j]; 
			}
		}
		break;
	}

	free(alfax);
	free(betax);

	return 0;
}

int Invert_Ay(double * sol, double * rhand, boundaryConditions yCondition, int Mx, int My, int Mz)
{
	double a_0_old, a_i_old, a_iminus1_old, a_nminus1_old;
	double b_0_old, b_i_old, b_n_old = 0;
	double temper_0, temper_k, temper_n = 0;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	double * alfay = (double *) malloc ((My + 1) * sizeof(double));
	double * betay = (double *) malloc ((My + 1) * sizeof(double));

	switch(yCondition)
	{
	case eNeumann:
		//инициализация для потока по y - Wy_n
		//Нейман, y
		for ( int i = 0 ; i < Mx*Mz ; i++ )
		{
			i1 = i/Mz;
			k1 = i - (i/Mz)*Mz;

			betay[0] = 0 ;
			betay[1] = rhand[i*(My+1) + 0];
			alfay[0] = 0 ;
			alfay[1] = 0 ;
			for ( int k = 1 ; k < My ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = rhand[i*(My+1) + k];

				if (k==My-1)
				{
					alfay[k + 1] = 0;
					betay[k + 1] = ( ( temper_k - a_i_old*rhand[i*(My+1) + My] ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k - a_iminus1_old*rhand[i*(My+1) + 0] ) ) / (  b_i_old );
				}
				else
				{
					alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
					betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				}
			}
			sol[i*(My+1) + My ] = rhand[i*(My+1) + My] ;
			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				sol[i*(My+1) + My - j ] = sol[i*(My+1) + My + 1 - j ] * alfay[My + 1 - j] + betay[My + 1 - j]; 
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


			//к-ты матрицы А
			a_0_old = hy(0)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = rhand[i*(My+1) + 0];

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
				temper_k = rhand[i*(My+1) + k];

				alfay[k + 1] = (-1)*a_i_old/(alfay[k]*a_iminus1_old + b_i_old);
				//betay[k + 1] = ( ( diff ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				//betay[k + 1] = ( ( diff2 ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );
				betay[k + 1] = ( ( temper_k ) - betay[k]*a_iminus1_old ) / ( alfay[k]*a_iminus1_old + b_i_old );

			}
			a_nminus1_old = hy(My-1)/(heat_conductivity(i1,My-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			temper_n = rhand[i*(My+1) + My];

			//Wy_n[i*(My+1) + My] = ( diff  -  betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);
			//Wy_n[i*(My+1) + My] = ( diff2  -  betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);
			sol[i*(My+1) + My] = ( temper_n  -  betay[My]*a_nminus1_old ) / (alfay[My]*a_nminus1_old + b_n_old);

			for ( int j = 1 ; j < My + 1 ; j++ )
			{
				sol[i*(My+1) + My - j ] = sol [i*(My+1) + My + 1 - j ] * alfay[My + 1 - j] + betay[My + 1 - j]; 
			}
		}
		break;

	}

	free(alfay);
	free(betay);

	return 0;
}

int Invert_Az(double * sol, double * rhand, boundaryConditions zCondition, int Mx, int My, int Mz)
{
	double a_i_old, a_iminus1_old, b_i_old, a_0_old, a_n_old, b_0_old, b_n_old;
	double temper_0, temper_k, temper_n;
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

			betaz[0] = 0 ;
			betaz[1] = rhand[i*(Mz+1) + 0];
			alfaz[0] = 0 ;
			alfaz[1] = 0 ;
			for ( int k = 1 ; k < Nz ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hz(k)/(heat_conductivity(i1,j1,k)*6.0);
				a_iminus1_old = hz(k-1)/(heat_conductivity(i1,j1,k-1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				temper_k = rhand[i*(Mz+1) + k];

				if (k==Nz-1)
				{
					alfaz[k + 1] = 0;
					betaz[k + 1] = ( ( temper_k - a_i_old*rhand[i*(Mz+1) + Mz] ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				}
				else if (k==1)
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = ( ( temper_k - a_iminus1_old*rhand[i*(Mz+1) + 0] ) ) / (  b_i_old );
				}
				else
				{
					alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
					betaz[k + 1] = ( ( temper_k ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				}
			}
			sol[i*(Nz+1) + Nz ] = rhand[i*(Mz+1) + Mz] ;
			for ( int j = 1 ; j < Nz + 1 ; j++ )
			{
				sol[i*(Nz+1) + Nz - j] = sol[i*(Nz+1) + Nz + 1 - j] * alfaz[Nz + 1 - j] + betaz[Nz + 1 - j];
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

			//к-ты матрицы А
			a_0_old = hz(0)/(heat_conductivity(i1,j1,0)*6.0);
			b_0_old = 2*a_0_old;
			temper_0 = rhand[i*(Mz+1) + 0];

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

				temper_k = rhand[i*(Mz+1) + k];

				alfaz[k + 1] = (-1)*a_i_old/(alfaz[k]*a_iminus1_old + b_i_old);
				betaz[k + 1] =  ( ( temper_k ) - betaz[k]*a_iminus1_old ) / ( alfaz[k]*a_iminus1_old + b_i_old );
				//			printf("temper_k = %f \n",temper_k);
			}
			a_n_old = hz(Nz-1)/(heat_conductivity(i1,j1,Nz-1)*6.0);
			b_n_old = 2*a_n_old;
			temper_n = rhand[i*(Mz+1) + Mz];
			//		printf("temper_n = %f \n",temper_n);

			sol[i*(Nz+1) + Nz] = ( temper_n  -  betaz[Nz]*a_n_old ) / (alfaz[Nz]*a_n_old + b_n_old);
			for ( int j = 1 ; j < Nz + 1 ; j++ )
			{
				sol[i*(Nz+1) + Nz - j] = sol[i*(Nz+1) + Nz + 1 - j] * alfaz[Nz + 1 - j] + betaz[Nz + 1 - j];
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

	return 0;

}
