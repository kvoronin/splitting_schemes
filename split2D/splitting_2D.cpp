#define GETCH_H
#define DEBUGE_H
#define ACCURACY
#define NEW_BOUNDARY_CONDITION_H
#define EXACT_BOUNDARY_CONDITION_H
#define ZEROING_H
#define FINEMESH_W0_INIT_H
#define FINEMESH_IND 12
#define EXACT_W0_INIT_H
#define SPECIAL_MODE_H //Запуск в схемы в качестве начального потока погрешности иницализации, 1 шаг по времени
#define EIGENVECTOR_WY0_H //использование собственных векторов и чисел для обращения матрицы A = нахождения начального потока через начальную температуру
#define SADDLE_POINT_H //использование седловой системы для нахождения начального потока и начальной температуры по известному лапласу температуры
#define SPECIAL_MODE_T_H  //использование в качестве начальной температуры погрешность между температурой в узлах сетки и температурой как интегральным средним по ячейкам
#define LUMPING_INIT_H
#define STEKLOV_H
#define STEKLOV_FLUX_H

#define HARMONICS_H
#define MIND_GAME_H
#define MINIMIZATION_H
#define MINIMIZATION_W0_FOR_SCHEME
//#define CG_TEST_H
#define CG_ZERO_INIT_H
#define NO_H_WEIGHT2
#define DEBUG_MINIM_H

//#ifdef GETCH
//#include<conio.h>
//#endif

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<string.h>
#include"src/data.h"
#include"src/solution.h"
#include"src/key_functions.h"
#include"src/eigenvectors.h"

#ifdef SPECIAL_MODE
	#ifndef DEBUGE
		#define DEBUGE
	#endif
#endif

#ifdef SADDLE_POINT
	#include<mkl_pardiso.h>
	#include<mkl_spblas.h>
	#include<mkl.h>
	#include"mkl_types.h"
#endif

#define BNDCNDS boundaryConditions

int	myBtr_w0 ( double * Btr_w0, double * wx, double * wy, int Nx, int Ny );

/*
Программа написана для двумерного случая, реализована схема Кранка - Николсон с помощью 
покомпонентного расщепления оператора в уравнении на w_t, есть конвекция (!). Введены тесты Арбогаста
и новый вариант интерфейса для тестирования.

Состояние:

22.10.2012
Скопировал код из Mfem_krnik_split2D_convection, ввожу новый вариант интерфейса и тесты Арбогаста.
Почти переделал, остается ввести g_po_t_0 и g_po_t_1. Вопрос, как их определить. Может, кстати, в этом и проблема
в проекте с алгоритмом Удзавы для переменного к-та теплопроводности.

За исключением g_po_t_0 и g_po_t_1 необходимые изменения внес, есть повторяемость старых результатов, точнее, 
резудьтатов работы программы Mfem_krnik_split2D_convection на тесте numSolution = 105, Dir_axis = 1.

Ввел g_po_t_0 и g_po_t_1, но тест numSolution = 108 нифига не работает, на первом же шаге считает чушь. Крешатся
причем даже квазиодномерные тесты по обоим направлениям.
Временно для поиска ошибки ввел для тестов Арбогаста краевые условия Неймана. Для Неймана тест 108 работает адекватно,
тесты 109 и 110  считают безобразно. Для половинчатых шагов явные косяки с краю.

Разобрался, что проблема в зависимости краевых условий от времени - для теста 108 с Нейманом все отлично, для теста
108 с Дирихле (зависящем от времени как раз) - все плохо. Или нет. Тест 105 с Dir_axis = 0 работает неверно.
Исправил, был косяк с перезатиранием краевых условий. Проблем дл зависящих от времени краевых условий это не решило.

Пришла в голову идея сделать в правой части полусумму и для полушагов по y. Внес соот-ве изменения для краевых условий Дирихле.
Расходимость исчезла, есть второй порядок, но абсолютные цифры не очень - для tau = 0.05, h = 1/20 eps_max_w ~ 4.0e+0
хотя хороший такой второй порядок.
Сделал соответствующие изменения и для краевых условий Неймана.
Если взять краевые условия Неймана по y и Дирихле по x и тест 108, нету сходимости в равномерной норме потока.
Для краевых условий Дирихле для теста 109 сходимость есть только по температуре приличная в л2-норме, по потоку
в л2-норме меньше первого, в равномерной вообще нету сходимости.
Для краевых условий Неймана нету сходимости для теста 110 ни в какой норме.
Для краевых условий Дирихле нету сходимости для теста 110 ни в какой норме.
Кстати, для любых краевых условий для теста 110 есть зависимость краевых условий от времени.

Считаю ключевым тот факт, что для теста sin(2mtPIt) (numsol = 108, mx=my=0) погрешность для потока в равномерной
норме аж порядка 17e+0. Если насильно занулить wy или wx на обоих дробных шагах, то считает нормально, погрешность
для потока порядка e-2, сходимость чуть не на порядок на 0.1-32 и 0.025-64.

18.11.2012
Решил разобраться с расходимостью в случае теста numSolution = 110. Провожу две серии тестов - с относительно маленьким фиксированным h
и с маленьким фиксированным tau.

При малом фиксированном h сходимость по tau какая-то (лучше, чем первый, но хуже, чем второй порядок) есть.
При малом фиксированном tau по h сходимости по потоку нету (а по температуре удивительнейшим образом на некоторых сетках есть), поэтому
написал тест 111 = тест 110 без времени (тупо убрал отовсюду из теста). Написал тест 112 - тест 110, симметричный
относительно середины области. 113 - модуль (x-0.5)^alpha., 114 - модуль (y-0.5)^alpha.
Тест 115 - только зависимость по времени, и не работает. (64 х 0.01) -> (64 x 0.005), (128 x 0.005).
Для теста 115 нашел ошибку в правой части, исправил.

03.12.2012
Решили, что не считает тест 110 из-за того, что есть слагаемое при аппроксимации порядка (tau/h)^2. Т.к. из-за условий на краю есть погрешность порядка tau^2,
от которой потом берется в прогонке по x вторая производная по пространству. Решили сделать тогда на краю точную производную по времени от краевого условия.

04.12.2012
Сделал, но не помогло. Ни тесту 110, ни тесту 109. (у которого едва-заметная сходимость в л2 и нету сходимости в равномерной норме потока).\
Проверку провожу для теста с my = 0. Помогло только тесту 108 c mx = my = 0 - теперь поток тождественно нулевой, а не 19 в равномерной норме на шагах 0.1-32.

05.12.2012
Проверил, что если занулять принудительно и wy_nplus05, и wy_nplus1, то считает хорошо.

17.12.2012
Временно поставил вместо полусуммы в wy_nplus05 в правую часть явный член f_n.

30.11.2013
Вернулся к исследованию поведения двумерной схемы расщепления для неравномерной сетки. ПРосто запустил тест 105 с dir_axis = 1, воспроизвел результаты из файла
Результаты split2DnoLump 01_03_2012.doc на равномерной и неравномерной сетках.

С неравномерной сеткой разобрался. Начал заново разбираться с переменным к-том теплопроводности.
Делаю новый тест 300 на замену 267, потому что у 267 главный член возникающего артефакта зануляется из-за того,
что тест из чистого синуса по y. Но что-то пока тест не работает.

Исправил ошибку в правой части righthand300 - не хватало минуса перед производной по времени.
Работает все равно только с external = 0, mx = 0. 

С external = 4, mx = 0 почти нет сходимости, особенно в L2 норме, и даже при tau/h^2 = const. 
Максимальная погрешность уже после одного шага, с ходом времени уменьшается, но для различных шагов сходимости также нет.
Исправил еще одну ошибку в righthand300, неправильно учитывался переменный к-т теплопроводности. 
Теперь нормально работает с external = 4, mx = 0, сходится как надо.

Но почему-то external = 0, mx = 2 при t/h = const дает сходимость со 2ым порядком в температуре и не сходится для потока.
И погрешность очень большая (коэффициент-то = 1 при external = 0, а погрешность для потока в равномерной норме порядка 1.0e+2)
Погрешность уже после одного шага в потоке по X большая, и потом еще увеличивается заметно.
Погрешность сидит в Wy_nplus05 уже, причем откровенно косяки на обоих краях. Что странно. Пила на краях в Wy0, невидимой амплитуды.
Она дает косяк в Wy_nplus05, а потом за пересчет потока по x еще вырастает, уже в монстра. А в Wx0 пилы нет, все гладенькое.
А на тесте, например, 105 с Dir_axis = 1 погрешность в начальном потоке Wy0 гладкая.
В чем дело?
Если программа работает правильно, то чем тест 300 принципиально хуже теста 105?

Итак, при tau/h константа погрешность потока в равномерной норме слегка расходится или не убывает. При этом в L2 норме температуры есть сходимость со вторым порядком.
Отсюда напрашивается вывод, что косяк в самой схеме, иначе как может сходиться температуры до е-4 и не сходиться поток порядка e+2.
Точная инициализация потока приводит к сходимости и потока, и температуры, и без всяких e+2.

Долго разбирался, но в итоге придумали - дело в том, что если брать в качестве начальной температуры значения в серединах ячеек, то на тесте 300 погрешность аппроксимации
на краю в уравнении Aw = BT имеет первый порядок (а у теста 105 просто зануляется эта производная), а вот если брать в качестве T0 интегральные средине по ячейкам, то и 
на краю второй порядок, что в итоге приводит к тому, что все хорошо и тест ходится со вторым порядком с приемлемыми цифрами при tau/h = const.

Использование интегральных средних решает все проблемы численно, но пока не удалось доказать это теоретически.

18.04.2013
Вернулся к зависящим от времени краевым условиям. Смотрю на тест 108. С условиями Неймана (постоянными) есть сходимость второго порядка.
С Дирихле+Нейман(переменным) одномерный тест - сходимость со 2 порядком. С Дирихле+Нейман двумерный тест в равномерной норме потока дает
константу. А с точной инициализацией потока все работает. Рабочая версия - что косяк в аппроксимации на краю возникает. Но странно, что
на одномерном тесте этого не наблюдается. И интегральные средние не помогают.

Итак
1) Тест 108 двумерный с условиями Нейман+Дирихле дает для температуры и потока в L2 норме сходимость с 2 и 1 порядками, а в равномерной
норме - с 1 и 0 порядком соответственно. При tau/h = const.
2) Тест 108 двумерный с условиями Нейман+Дирихле дает для температуры и потока в L2 норме сходимость с 3 и 2 порядками, а в равномерной
норме - с 2 и 1 соответственно. При tau/h^2 = const.
3) C инициализацией потока точными значениями - во всех нормах сходимость с 2 порядком.
4) Использование интегральных средних не помогло.

ДЛЯ СТАЦИОНАРНОГО ТЕСТА та же фигня типа 1). Сначала надо понять, что на стационарном тесте не так.
более того, проблемы уже для теста 108 с mt = my = 0. Краевые условия - Нейман + Дирихле(у). То есть для одномерного теста.
А для Нейман + Нейман все хорошо. Плохо становится, когда для теста sin(2 pi x) ставить Дирихле по y. Хм....
Точная инициализация потока приводит к сходимости, а интегральные средние температуры - нет.
В правой части для wy_nplus05 на краю всплеск порядка h^2 из-за погрешности в начальном потоке.
Вообще погрешность в углах области (т.е. несмотря на константу в равномерной норме для ошибки в потоке, в L2 есть сходимость
для потока). Т.е. фигня с неверной аппроксимацией краевых условий в углах области.

23.04.2014 
Решил оформить код, как следует - все в виде функций. Черновой вариант готов. Правда, заполнение правой части для прогонок пока как одна функция.
Можно разбить ее на отдельные части - BxF, ByF, BxMBxWx+BxMByWy, ByMBxWx+ByMByWy плюс конвективные члены еще, например.

26.05.2014
Сделал новые функции для правых частей прогонок, за исключением конвективных слагаемых.
Конвекцию тоже сделал, но не сравнивал со старой версией. Разок запустил на тесте 105 с 2 скоростью - разумные цифры вылезли (а в старой версии
почему-то 15 знаков до запятой на тех же шагах  1/32 - 0.1). Не стал разбираться, почему.

02.06.2014
Сделал и прогнал на одном тесте (105, nonstationary, mx = my = 2) алгоритм Удзавы. Поменял местами x и у по сравнеyию с тем кодом, что был, - чтобы больше соответствовало моей схеме расщепления.
Но возможны путаница при сравнении кодов - надо x и y в условиях теста поменять местами, не забыть про краевые условия.

Проверил - если сделать первый шаг по Арбогасту, а потом считать по схеме расщепления, то лучше не становится - "говно" побеждает сразу же.

04.06.2014
Сделал локально одномерную схему в том же ключе. Сравнил на тесте 105 с локально-одномерной 3D, которая была раньше, - работают одинаково.
Пытаюсь реализовать двуциклическую локально-одномерную схему, пока не удается. Погрешность после одного шага большая - поток по форме напоминает точный,
но отличается в странное количество раз, которое сложным образом зависит от шагов сетки. На одномерном тесте 105 считает не так, как должно. 
Как искать ошибку, пока непонятно, тут немного необычно по сравнению с другими схемами включается правая часть.
Погрешность не меняется почти, когда я вообще убираю этап с прибавлением правой части. Погрешность того же порядка (изменения доли процентов)
- что есть этот этап, что нет. Это странно и так быть не должно. При этом Wx_nplus05 (результат этапа с правой частью) меняется сильно. Но несмотря на то,
большое или маленькое Wx_nplus05, Wx_nplus1 все равно маленькое.

05.06.2014
Решил отступить и сделать классический предиктор-корректор сначала. Предиктор-корректор заработал, приемлемые цифры, h^2.
Но для этого пришлось отклониться от того, что в книжке и развернуть потоковую схему так, чтобы не возникало необходимости
отдельно считать поток, где нет вторых производных, а есть только првая часть от исходного дифф. уравнения. Если сделать в лоб (см. step 1 и step 2
вместо alternative step 1-2), то цифры громадные, хотя сходимость со вторым порядком есть.

Переписал схему двуциклического покомпонентного расщепления  в лоб как две локально-одномерные (alernative steps), появился второй порядок.
При этом правые части в t_nplus14 и t_nplus34 аппроксимировал (0.75*f_n + 0.25*f_nplus1) и (0.25*f_n + 0.75*f_nplus1) соответственно.
Есть второй порядок, цифры приемлемые.

15.12.2014
Решил протестировать аппроксимацию конвективных слагаемых на гиперболичесом уравнении переноса.

*/

double tau;
//int Nx, Ny, Nz;
//int external;
int numSolution;


#ifdef GETCH
const int pause = 1;
#endif
#ifdef DEBUGE
const int printcase = 1;
#endif
const int print_step = 10;
//const int nonstationary = 1;
double Time = 1.;		//промежуток времени [0,Time]
double* Wx_n;			// массив для значений x-компоненты теплового потока на n-ом временном слое
double* Wy_n;			// массив для значений x-компоненты теплового потока на n-ом временном слое
double* Wz_n;			// массив для значений x-компоненты теплового потока на n-ом временном слое
double* Vx;
double* Vy;
double* Vz;
double* T;				// массив для значений температуры (в процессе вычислений - на n-ом временном слое, но в конце основного цикла хранит результирующие значения температуры

double *xpoints, *ypoints, *zpoints;

double* alfax;		    // массив переменного размера для хранения коэфиициентов альфа i-ых из прогонки, применяемой к x-компонентам
double* betax;			// массив переменного размера для хранения коэфиициентов бета i-ых из прогонки, применяемой к x-компонентам
double* alfay;
double* betay;

double exact_solution(int number, double x, double y, double z, double t); 
double exact_gradientX(int number, double x, double y, double z, double t);
double exact_gradientY(int number, double x, double y, double z, double t); 
double exact_gradientZ(int number, double x, double y, double z, double t);
double righthand(int number, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
//void setBoundaryConditions(int number);
void Output_filename(char* filename, int numSolution);
void meshGeneration(int version, char **argv);
void allocation();
double norm_l2(double * massiv, int N, int M);
double scal_l2(double * massiv1, double * massiv2, int N, int M);
double norm_max(double * wx, int Nx, int Ny);
void divergence_norm(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt);
void stability_norm(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt);
void stability_normFull(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt);
void stability_norm_sum(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt);

int	myMhalfLambdayBtr_w0 ( double * MhalfLambdaY_Btr_w0, double * wx, double * wy, int Nx, int Ny );
int	myMminushalfLambdaLambdayBtr_w0 ( double * Mminus1Btr_Aminus1_xvost_w0, double * MhalfLambdaY_Btr_w0, double * wx, double * wy, int Nx, int Ny );
int	myAyminus1ByMminus1Btr_w0 ( double * Ayminus1ByMminus1_Btr_w0, double * wx, double * wy, int Nx, int Ny );


double func_gamma_k(int k, int n);
void func_u_k(double *u, int k, int n);
void func_p_k(double *p, int k, int n);
double dot_product(double* a, double* b, int size);
double second_norm(double* a, int size);
void divergence_Wy_norm(double* wy, int size, double *l2_norm_pt, double *max_norm_pt);
double exact_gradient_harm(int harm, double y);
double exact_harm(int harm, double y);
double my_solution(double x,double y,double z,double t);

double saddle_righthand_laplace300(int i, int j, int m, int n);
double saddle_righthand_laplace(int numSolution, int i, int j, int m, int n);
double saddle_righthand_laplace105(int i, int j, int m, int n);
double exact_laplace105(double x, double y, double z, double time, int m, int n);

int myFunctionalGrad ( double * VecX, double * VecY, double * T, double * GradFuncX, double * GradFuncY, int Nx, int Ny, double weight1, double weight2 );
double myFunctional ( double * VecX, double * VecY, double * T, int Nx, int Ny, double weight1, double weight2 );
int	myAvec ( double *tempvecX, double *tempvecY, double *VecX, double *VecY, int Nx, int Ny );
int myNevyazka ( double *nevyazkaX, double * nevyazkaY, double * VecX, double * VecY, double * T, int Nx, int Ny );
double normAQminusF ( double * VecX, double * VecY, double * T, int Nx, int Ny, double * l2_norm_pt );
int myMinimization ( double * Wx_0, double * Wy_0, double * T_0, double ** Wx_0_new_pt, double ** Wy_0_new_pt, double converg_tol, double stagn_tol, double grad_tol, double coeff1, double coeff2, int iter_limit, int stagnation_limit, double h, int Nx, int Ny, int Nz );

int stability_norm2D ( double * Wx, double * Wy, int Nx, int Ny, double tau, BNDCNDS xCondition, BNDCNDS yCondition, double * res_pt);
int	compute_Lambdax2D ( double * out, double * in, int Nx, int Ny, BNDCNDS xCondition, BNDCNDS yCondition);
int	compute_Lambday2D ( double * out, double * in, int Nx, int Ny, BNDCNDS xCondition, BNDCNDS yCondition);
int	compute_Lambda2D ( double * out, double * in, int Nx, int Ny, BNDCNDS xCondition, BNDCNDS yCondition);
int zero_init ( double * vec, int dim );
int Invert_Ax(double * sol, double * rhand, boundaryConditions xCondition, int Mx, int My, int Mz);
int Invert_Ay(double * sol, double * rhand, boundaryConditions yCondition, int Mx, int My, int Mz);


int maximum(int N, int M, int K);
int stepen2(int n);

FILE* f1 = NULL;
FILE* f1_excel = NULL;


int main(int argc, char **argv)
{

	// argc - кол-во параметров командной строки. Т.к. ты запускаешь ехе-шник свой, то само имя программы тоже параметр - потому argc минимум равен 1. argv[0] = ProgrammName.exe
	if(argc == 1) // значит без параметров
	{
		printf("Enter parameters tau, Nx, Ny, Nz, external and numSolution \n");
		scanf("%lf", &tau);
		scanf("%d", &Nx);
		scanf("%d", &Ny);
		scanf("%d", &Nz);
		scanf("%d", &external);
		scanf("%d", &numSolution);
	}else if(argc == 7)
	{
		tau = atof(argv[1]);
		Nx = atoi(argv[2]);
		Ny = atoi(argv[3]);
		Nz = atoi(argv[4]);
		external = atoi(argv[5]);
		numSolution = atoi(argv[6]);
		//там ещё посмотри параметры. Где-то код ошибки смотреть perror() какой-нибудь.. сходу не напишу)
	} else
	{
		printf("Incorrect parameters \n");
		exit(-1);
	}

	Nz = 1;

	//meshGeneration(argc, argv);
	allocation(); //set x(y,z)points 
	setBoundaryConditions(numSolution); //eDirichlet eNeumann eMixed
	printf("xCondition = %d yCondition = %d \n",xCondition, yCondition);

	char filename[100];
	char addname1[15];
	char addname2[15];
	char addname3[15];
	char addname4[15];
	char addname5[15];

	Output_filename(filename, numSolution);

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

	char filename_excel[100];
	char addname1_excel[15];
	char addname2_excel[15];
	char addname3_excel[15];
	char addname4_excel[15];
	char addname5_excel[15];

	Output_filename(filename_excel, numSolution);

	sprintf(addname5_excel,"%d_",external);
	sprintf(addname1_excel,"%f_",tau);
	sprintf(addname2_excel,"%d_",Nx);
	sprintf(addname3_excel,"%d_",Ny);
	sprintf(addname4_excel,"%d",Nz);
	strcat(filename_excel,addname5_excel);
	strcat(filename_excel,addname1_excel);
	strcat(filename_excel,addname2_excel);
	strcat(filename_excel,addname3_excel);
	strcat(filename_excel,addname4_excel);
	strcat(filename_excel,".xls");
	f1_excel = fopen(filename_excel,"w");


	int N = int(Time/tau);
	//	printf("N = %d\nPress any key to start calculations\n",N);
	int dim_T = Nx*Ny*Nz;
	int dim_Wx = (Nx+1)*Ny*Nz;
	int dim_Wy = Nx*(Ny+1)*Nz;
	int dim_Wz = Nx*Ny*(Nz+1);
	int dim_W = dim_Wx + dim_Wy + dim_Wz;
	//	printf("dim_T = %d \ndim_Wx = %d \ndim_Wy = %d \ndim_Wz = %d \n",dim_T, dim_Wx, dim_Wy,dim_Wz);

	///////Определяем размеры массивов для компонент скорости и теплового потока
	Vx = new double [dim_Wx];
	Vy = new double [dim_Wy];
	//Vz = new double [dim_Wz];
	Wx_n = new double [dim_Wx];
	Wy_n = new double [dim_Wy];
	//Wz_n = new double [dim_Wz];
	T = new double [dim_T];

	// инициализация Vx, Vy, Vz , причем вектор V=(Vx,Vy,Vz) таков, что V*n(скалярное произведение на нормаль) на границе области = 0
	for ( int i = 0 ; i < dim_Wx ; i++ )
	{
		int i1 = i%(Nx+1);
		int j1 = ((i - i1)/(Nx+1))%Ny;
		int k1 = (i - i1 - j1*(Nx+1))/((Nx+1)*Ny);
		//		Vx[i] = v_init(i1*hx, (j1+0.5)*hy, (k1+0.5)*hz);

		Vx[i] = exact_solution1111Vx(xpoints[i1],ypoints[j1] + 0.5*hy(j1),zpoints[k1] + 0.5*hz(k1),0,mx, my, mz);
	}
	for ( int i = 0 ; i < dim_Wy ; i++ )
	{
		int j1 = i%(Ny+1);
		int k1 = ((i - j1)/(Ny+1))%Nz;
		int i1 = (i - j1 - k1*(Ny+1))/((Ny+1)*Nz);

		Vy[i] = exact_solution1111Vy(xpoints[i1] + 0.5*hx(i1),ypoints[j1],zpoints[k1] + 0.5*hz(k1),0,mx, my, mz);
	}

	printf("Initializing T...\n");

#ifdef SPECIAL_MODE_T
	double* T_error = (double*)malloc(dim_T*sizeof(double));
#endif
	//////Инициализация T
	for ( int i = 0 ; i < Nx ; i++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				T[i + j*Nx + k*Nx*Ny] = exact_solution(numSolution,xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0);
				//printf ( "value = %f \n", T[i + j*Nx + k*Nx*Ny] );
				//_getch();
#ifdef SPECIAL_MODE_T
				T_error[i + j*Nx + k*Nx*Ny] = T[i + j*Nx + k*Nx*Ny] - exact_average_solution300(i, j, 0, tau, mx, my);
				T[i + j*Nx + k*Nx*Ny] = T_error[i + j*Nx + k*Nx*Ny];
#endif
				//T[i + j*Nx + k*Nx*Ny] = my_solution(xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), 0, 0);
				//T[i + j*Nx + k*Nx*Ny] = exact_average_solution300(i, j, 0, tau, mx, my);
#ifdef STEKLOV			
				if (numSolution == 105)
					T[i + j*Nx + k*Nx*Ny] = exact_average_solution_other105(i, j, 0, tau, mx, my);
				if (numSolution == 108)
					T[i + j*Nx + k*Nx*Ny] = arbSmooth_average_ex1_other_solution(i, j, 0, tau, mx, my, mt);
				if (numSolution == 300)
					T[i + j*Nx + k*Nx*Ny] = exact_average_solution_other300(i, j, 0, tau, mx, my);
				if (numSolution == 400)
					T[i + j*Nx + k*Nx*Ny] = exact_average_solution_other400(i, j, 0, tau, mx, my);
#endif
				
				//double temp = exact_average_solution_other300(i, j, 0, tau, mx, my);
				//if (fabs(temp -  T[i + j*Nx + k*Nx*Ny]) >= 1.0e-13)
				//	printf("something wrong: temp = %e, T[%d] = %e \n", temp, i + j*Nx + k*Nx*Ny, T[i + j*Nx + k*Nx*Ny]);
			}
		}
	}

#ifdef SPECIAL_MODE_T
//	double tem = -(1.0 / (PI * my))*(ypoints[0]*cos(PI * my * ypoints[0]) - (-hy(0))*cos(PI * my * (-hy(0)))) + (1.0 / (PI * my * PI * my))*(sin(PI * my * ypoints[0]) - sin(PI * my * (-hy(0))));
//	printf("T_-1 = %f \n",nonstationary * exp (0 * tau) * tem * hx(0) / (hx(0)*hy(0)));
#endif

	FILE* T0_file = fopen("T0_split2D.xls","wt");
	for (int k = 0; k < Nz; k++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for ( int i = 0 ; i < Nx ; i++ )
			{
#ifdef SPECIAL_MODE_T
				fprintf(T0_file,"%f \t %f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),exact_average_solution_other400(i,j,0,tau,mx,my), T_error[i + j*Nx + k*Nx*Ny]);
#else
				fprintf(T0_file,"%f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),arbSmooth_average_ex1_other_solution(i,j,0,tau,mx,my,mt));
#endif
			}
		}
	}
	fclose(T0_file);

#ifdef DEBUGE
	FILE* T0_file = fopen("T0_split2D.xls","wt");
	for (int k = 0; k < Nz; k++)
	{
		for ( int i = 0 ; i < Nx ; i++ )
//		for (int j = 0; j < Ny; j++)
		{
//			for ( int i = 0 ; i < Nx ; i++ )
			for (int j = 0; j < Ny; j++)
			{
#ifdef SPECIAL_MODE_T
				fprintf(T0_file,"%f \t %f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),exact_average_solution_other400(i,j,0,tau,mx,my), T_error[i + j*Nx + k*Nx*Ny]);
#else
				//fprintf(T0_file,"%f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),exact_average_solution400(i,j,0,tau,mx,my));
				fprintf(T0_file,"%f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),arbSmooth_average_ex1_other_solution(i,j,0,tau,mx,my,mt));
				//fprintf(T0_file,"%f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0));
#endif
			}
		}
	}
	fclose(T0_file);
#ifdef SPECIAL_MODE_T
	printf("l2_norm T_error = %f \n",norm_l2(T_error, Ny, 1));
#endif


	FILE* rhand_file = fopen("righthand_split2D.xls","wt");
	for ( int i = 0 ; i < Nx ; i++ )
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				fprintf(rhand_file,"%f \n",righthand(numSolution,0,i,j,k,xpoints,ypoints,zpoints));
			}
		}
	}
	fclose(rhand_file);
#endif
	//////////////////////////
	printf("Initializing W...\n");

	{;}

	// расчет из Aw - BT = G теплового потока w0 ~ инициализация w0
	////////////////////////////////////////***************** Step 0
	//T0 - - > W0 = (Wx_0, Wy_0, Wz_0 )transp
	//( W0 = A(-1) * ( G - B * Tn ) progonka!
	{;}

	HeatFlux_Init(T, Wx_n, xCondition, Wy_n, yCondition, numSolution, Nx, Ny, Nz);

	///temporary block for checing filtration code
	//printf ( "Wx_0 \n" );
	//_getch();
	//for ( int i = 0; i < dim_Wx; i++ )
	//	printf ( "Wx_0[%d] = %f \n", i, Wx_n[i] );
	//_getch();
	//return 0;

	//HeatFlux_Wx0_Init(T, Wx_n, xCondition, numSolution, Nx, Ny, Nz);
	//HeatFlux_Wy0_Init(T, Wy_n, yCondition, numSolution, Nx, Ny, Nz);

	FILE * Wx_0_file_NEW = fopen("Wx0_special.xls","wt");
	for (int i = 0; i < dim_Wx; i++)
		fprintf(Wx_0_file_NEW,"%f \n", Wx_n[i]);
	fclose(Wx_0_file_NEW);
	FILE * Wy_0_file_NEW = fopen("Wy0_special.xls","wt");
	for (int i = 0; i < dim_Wy; i++)
		fprintf(Wy_0_file_NEW,"%f \n", Wy_n[i]);
	fclose(Wy_0_file_NEW);

#ifdef FINEMESH_W0_INIT
	int fine_ind = FINEMESH_IND;
	int fineNx = fine_ind * Nx;
	int fineNy = fine_ind * Ny;

	double * fine_xpoints, * fine_ypoints;
	fine_xpoints = (double*)malloc(sizeof(double)*(fineNx+1));
	fine_xpoints[0] = 0;
	for (int i = 1; i <= fineNx; i++)
	{
		fine_xpoints[i] = fine_xpoints[i-1] + hx( (i/fine_ind)-1 ) / fine_ind;
		printf ( "fine_xpoints[%d] = %f \n", i, fine_xpoints[i]);
		if (i % fine_ind == 0 )
			printf ( "xpoints[%d] = %f \n", i/fine_ind, xpoints[i/fine_ind]);
	}
	fine_ypoints = (double*)malloc(sizeof(double)*(fineNy+1));
	fine_ypoints[0] = 0;
	for (int j = 1; j <=fineNy; j++)
	{
		fine_ypoints[j] = fine_ypoints[j-1] + hy( (j/fine_ind)-1 ) / fine_ind;
		printf ( "fine_ypoints[%d] = %f \n", j, fine_ypoints[j]);
		if (j % fine_ind == 0 )
			printf ( "ypoints[%d] = %f \n", j/fine_ind, ypoints[j/fine_ind]);
	}



	double * fineT0 = (double *)malloc ( fineNx * fineNy * sizeof ( double ) );
	
	for ( int i = 0 ; i < fineNx ; i++ )
	{
		for (int j = 0; j < fineNy; j++)
		{
			fineT0[i + j*fineNx] = exact_solution(numSolution,fine_xpoints[i]+0.5*hx(i/fine_ind)/fine_ind, fine_ypoints[j]+0.5*hy(j/fine_ind)/fine_ind, 0.5, 0);
//			printf ( "fineT0[%d] = %f \n", i + j*fineNx, fineT0[i + j*fineNx],  );
//			if ( i % fine_ind == 0 && j % fine_ind == 0 )
//				printf ( "T[%d] = %f \n",i/fine_ind + (j / fine_ind) * Nx, T[i/fine_ind + (j / fine_ind) * Nx] );

//			_getch();
		}
	}

	//T0_file = fopen("T0_check_FINEMESH.xls","wt");
	//for (int k = 0; k < Nz; k++)
	//{
	//	for (int j = 0; j < fineNy; j++)
	//	{
	//		for ( int i = 0 ; i < fineNx ; i++ )
	//		{
	//			fprintf(T0_file,"%f \t %f \n",fineT0[i + j*fineNx],exact_solution(numSolution,fine_xpoints[i]+0.5*hx(i/fine_ind)/fine_ind, fine_ypoints[j]+0.5*hy(j/fine_ind)/fine_ind, 0.5, 0) );
	//		}
	//	}
	//}
	//fclose(T0_file);

	double * fineWx0 = (double *) malloc ( ( fineNx + 1 ) * fineNy * sizeof ( double ) );
	double * fineWy0 = (double *) malloc ( fineNx * ( fineNy + 1 ) * sizeof ( double ) );

	//works incorrectly because of usage of global hx, xpoints everywhere inside
	//HeatFlux_Init(fineT0, fineWx0, xCondition, fineWy0, yCondition, numSolution, fineNx, fineNy, Nz);

	FILE * f_spec;
//	char filename[100];
	char addname[10];

	strcpy ( filename, "Wx0_special" );
	sprintf ( addname, "_%d.xls", fine_ind * Nx );
	strcat ( filename, addname );

	f_spec = fopen ( filename, "rt" );
	if ( f_spec == NULL )
	{
		printf ( "Cannot open file %s for reading Wx0 special \n", filename );
		return -1;
	}

	for ( int k = 0; k < ( fineNx + 1 ) * fineNy; k++ )
		fscanf ( f_spec, "%lf \n", & (fineWx0[k]) );

	fclose (f_spec);

	strcpy ( filename, "Wy0_special" );
	sprintf ( addname, "_%d.xls", fine_ind * Ny );
	strcat ( filename, addname );

	f_spec = fopen ( filename, "rt" );
	if ( f_spec == NULL )
	{
		printf ( "Cannot open file %s for reading Wy0 special \n", filename );
		return -1;
	}

	for ( int k = 0; k < ( fineNy + 1 ) * fineNx; k++ )
		fscanf ( f_spec, "%lf \n", & (fineWy0[k]) );

	fclose (f_spec);


	//works incorrectly because of hx, hy usage inside
	//{
	//double eps0_max_w = 0.0;
	//double eps0_max_T = 0.0;
	//double eps0_l2_w = 0.0;
	//double eps0_l2_T = 0.0;
	//double eps0_relative_max_w = 0.0;
	//double eps0_relative_max_T = 0.0;
	//double eps0_relative_l2_w = 0.0;
	//double eps0_relative_l2_T = 0.0;
	//Accuracy_calculate(fineT0, fineWx0, fineWy0, numSolution, -1, fineNx, fineNy, Nz, print_step, &eps0_max_w, &eps0_max_T, &eps0_l2_w, &eps0_l2_T, &eps0_relative_max_w, &eps0_relative_max_T, &eps0_relative_l2_w, &eps0_relative_l2_T);
	//}

	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i <= Nx; i++ )
		{
//			printf ( "i = %d j = %d \n", i, j);
//			if ( i == 0 && j == 0 )
//				printf ( "indx = %d indy = %d index_total = %d \n", fine_ind * i, ( (fine_ind/2) + fine_ind * j ), fine_ind * i + ( (fine_ind/2) + fine_ind * j ) * ( fineNx + 1 ));
			Wx_n[i + j * ( Nx + 1 )] = fineWx0[fine_ind * i + ( (fine_ind/2) + fine_ind * j ) * ( fineNx + 1 ) ];
//			_getch();
		}
	}

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j <= Ny; j++ )
		{
//			if ( i == 0 && j == 0 )
//				printf ( "indx = %d indy = %d index_total = %d \n", fine_ind * j, ( (fine_ind/2) + fine_ind * i ), fine_ind * j + ( (fine_ind/2) + fine_ind * i ) * ( fineNy + 1 ) );
			Wy_n[j + i * ( Ny + 1 )] = fineWy0 [fine_ind * j + ( (fine_ind/2) + fine_ind * i ) * ( fineNy + 1 ) ];
		}
	}

	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
//			if ( i == 0 && j == 0 )
//				printf ( "indx = %d indy = %d index_sum = %d \n", (fine_ind/2 ) + fine_ind * i, ( (fine_ind/2 + 1) + fine_ind * j ), (fine_ind/2) + fine_ind * i + ( (fine_ind/2 + 1) + fine_ind * j ) * fineNx );
			T[i + j * Nx ] = fineT0[ (fine_ind/2 ) + fine_ind * i + ( (fine_ind/2) + fine_ind * j ) * fineNx ];
		}
	}

	Wx_0_file_NEW = fopen("Wx_0_NICE_interface_FINEMESH.xls","wt");
	for (int i = 0; i < dim_Wx; i++)
		fprintf(Wx_0_file_NEW,"%f \n", Wx_n[i]);
	fclose(Wx_0_file_NEW);
	Wy_0_file_NEW = fopen("Wy_0_NICE_interface_FINEMESH.xls","wt");
	for (int i = 0; i < dim_Wy; i++)
		fprintf(Wy_0_file_NEW,"%f \n", Wy_n[i]);
	fclose(Wy_0_file_NEW);

	T0_file = fopen("T0_split2D_FINEMESH.xls","wt");
	for (int k = 0; k < Nz; k++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for ( int i = 0 ; i < Nx ; i++ )
			{
#ifdef SPECIAL_MODE_T
				fprintf(T0_file,"%f \t %f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),exact_average_solution_other400(i,j,0,tau,mx,my), T_error[i + j*Nx + k*Nx*Ny]);
#else
				fprintf(T0_file,"%f \t %f \t %f \n",T[i + j*Nx + k*Nx*Ny],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), zpoints[k]+0.5*hz(k), 0),arbSmooth_average_ex1_other_solution(i,j,0,tau,mx,my,mt));
#endif
			}
		}
	}
	fclose(T0_file);

#endif
	

	double a_0_old, a_i_old, a_iminus1_old, a_nminus1_old = 0;
	double b_0_old, b_i_old, b_n_old = 0;
	double temper_0, temper_k, temper_n = 0;
	double Tx_0_temp, Tx_1_temp, Ty_0_temp, Ty_1_temp, wx_0_temp, wx_1_temp, wy_0_temp, wy_1_temp;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	
	alfax = new double[Nx + 1];
	betax = new double[Nx + 1]; 
	alfay = new double[Ny + 1];
	betay = new double[Ny + 1]; 

#ifdef LUMPING_INIT
	double* Aydiag = (double*)malloc(dim_Wy * sizeof(double));
	switch(yCondition)
	{
	case eDirichlet:
		//Дирихле, y
		{;}
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;

			a_0_old = hy(0)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;

			Aydiag[i*(Ny+1) + 0] = b_0_old + a_0_old;
			for ( int k = 1 ; k < Ny ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				Aydiag[i*(Ny+1) + k] = a_i_old + a_iminus1_old + b_i_old;
			}

			a_nminus1_old = hy(Ny-1)/(heat_conductivity(i1,Ny-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;
			Aydiag[i*(Ny+1) + Ny] = a_nminus1_old + b_n_old;
		}
		break;
	default:
		break;
	}

	double * Wy_0_lumping = (double*)malloc(dim_Wy * sizeof(double));
	double * AlumpWy0_exact = (double*)malloc(dim_Wy * sizeof(double));
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny + 1; j++)
			AlumpWy0_exact[i*(Ny + 1) + j] = Aydiag[i*(Ny + 1) + j] * (-exact_gradientY(numSolution, xpoints[i] + 0.5*hx(i), ypoints[j], 0, 0)); 
	double ByT_0, ByT_k, ByT_n;
	FILE* rhand_lumping = fopen("righthand_lumping.xls","wt");
	double diff, diff2;

	switch(yCondition)
	{
	case eDirichlet:
		//Дирихле, y
		{;}
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;

			ByT_0 = (Ty_0 - T[i1 + 0*Nx + k1*Nx*Ny]);
			diff = ByT_0 - (-hy(0)*exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[0], 0, 0));
			diff2 = ByT_0 - AlumpWy0_exact[i*(Ny + 1) + 0];
			fprintf(rhand_lumping,"%f \t %f \t %f \t %f \n",ByT_0, -hy(0)*exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[0], 0, 0), diff, diff2);
			//Wy_0_lumping[i*(Ny+1) + 0] = diff / Aydiag[i*(Ny+1) + 0];
			//Wy_0_lumping[i*(Ny+1) + 0] = diff2 / Aydiag[i*(Ny+1) + 0];
			Wy_0_lumping[i*(Ny+1) + 0] = ByT_0 / Aydiag[i*(Ny+1) + 0];
			for ( int k = 1 ; k < Ny ; k++ )
			{
				ByT_k = (T[i1 + (k - 1)*Nx + k1*Nx*Ny] - T[i1 + k*Nx + k1*Nx*Ny]);
				diff = ByT_k - (-hy(0)*exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[k], 0, 0));
				diff2 = ByT_k - AlumpWy0_exact[i*(Ny + 1) + k];
				fprintf(rhand_lumping,"%f \t %f \t %f \t %f \n",ByT_k,-hy(0)*exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[k], 0, 0), diff, diff2);
				//Wy_0_lumping[i*(Ny+1) + k] = diff / Aydiag[i*(Ny+1) + k] ;
				//Wy_0_lumping[i*(Ny+1) + k] = diff2 / Aydiag[i*(Ny+1) + k] ;
				Wy_0_lumping[i*(Ny+1) + k] = ByT_k / Aydiag[i*(Ny+1) + k] ;
			}
 
			ByT_n = (T[i1 + (Ny-1)*Nx + k1*Nx*Ny] - Ty_1);
			diff = ByT_n - (-hy(0)*exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[Ny], 0, 0));
			diff2 = ByT_n - AlumpWy0_exact[i*(Ny + 1) + Ny];
			fprintf(rhand_lumping,"%f \t %f \t %f \t %f \n",ByT_n,-hy(0)*exact_gradientY(numSolution, xpoints[i1] + 0.5*hx(i1), ypoints[Ny], 0, 0), diff, diff2);
			//Wy_0_lumping[i*(Ny+1) + Ny] = diff / Aydiag[i*(Ny+1) + Ny];
			//Wy_0_lumping[i*(Ny+1) + Ny] = diff2 / Aydiag[i*(Ny+1) + Ny];
			Wy_0_lumping[i*(Ny+1) + Ny] = ByT_n / Aydiag[i*(Ny+1) + Ny];
		}
		break;
	default:
		break;
	}
	fclose(rhand_lumping);

	double * Wy_0_lumping_error = (double*)malloc(dim_Wy * sizeof(double));
	for (int i=0; i<Nz*Nx; i++)
	{
		for (int k=0; k<Ny+1; k++)
		{
			i1 = i/Nz;
			j1 = k;
			k1 = i - (i/Nz)*Nz;

			Wy_0_lumping_error[i*(Ny+1) + k] = Wy_0_lumping[i*(Ny+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) );
		}
	}
#ifdef DEBUGE
	{
		char* ext = ".xls";
		char nx_string[5];
		sprintf(nx_string,"_%d",Nx);
		char filename[30];
		strcpy(filename,"Wy_lumping");
		strcat(filename,nx_string);
		strcat(filename,ext);
		FILE * Wy_file = fopen(filename,"wt");
		for (int ind = 0; ind < dim_Wy; ind++)
			//fprintf(Wy_file,"%f \n",Wy_0_lumping_error[ind]);
			fprintf(Wy_file,"%f \n",Wy_0_lumping[ind]);
		FilePrint2D("Wy_lumping_error_2D.xls",Wy_0_lumping_error, Ny + 1, Nx);
		fclose(Wy_file);

		FILE* T0_y_file = fopen("T0_y_split2D.xls","wt");
		for ( int i = 0 ; i < Nx ; i++ )
		{
			for (int j = 0; j < Ny; j++)
			{
				fprintf(T0_y_file,"%f \t %f \t %f \n",T[i + j*Nx],exact_solution(numSolution, xpoints[i]+0.5*hx(i), ypoints[j]+0.5*hy(j), 0, 0), exact_average_solution300(i,j,0,tau,mx,my));
			}
		}
		fclose(T0_y_file);
		printf("Ty_0 = %f Ty_1 = %f \n", Ty_0, Ty_1);

	}
#endif

	double div_W1D_lump_l2 = 0.0, div_W1D_lump_max = 0.0;
	divergence_Wy_norm(Wy_0_lumping, Ny, &div_W1D_lump_l2, &div_W1D_lump_max);
	printf("div_W1D_lump_l2 = %f \n",div_W1D_lump_l2);
	printf("div_W1D_lump_max = %f \n",div_W1D_lump_max);
#endif //LUMPING_INIT

#ifdef EIGENVECTOR_WY0
	{
		double * u_k = (double*)malloc(Ny * sizeof(double));
		double * p_k = (double*)malloc((Ny + 1) * sizeof(double));
		double * c_k = (double*)malloc(Ny * sizeof(double));
		double * c_k_new = (double*)malloc(Ny * sizeof(double));
		double * c_k_div = (double*)malloc(Ny * sizeof(double));
		double * Wy_0_eigenv = (double*)malloc((Ny + 1) * sizeof(double));
		double * Wy_0_error_eigenv = (double*)malloc((Ny + 1) * sizeof(double));
		double * T0 = (double*)malloc(Ny * sizeof(double));
		double * Btr_Wy_0_eigenv = (double*)malloc(Ny * sizeof(double));

		double * diver_harmonics_l2 = (double*)malloc(Ny * sizeof(double));
		double * diver_harmonics_max = (double*)malloc(Ny * sizeof(double));

		for (int harm = 1; harm <= Ny; harm++)
		{
			func_u_k(T0, harm, Ny);

			//decomposing T_0 onto span of u_k = calculating c_k's;
			for (int k = 1; k <= Ny; k++)
			{
				func_u_k(u_k, k, Ny);
				c_k[k-1] = dot_product(T0, u_k, Ny) / dot_product(u_k, u_k, Ny);
			}

			//if (harm == 4)
			//	for (int i = 0; i < Ny; i++)
			//		printf("c_k[%d] = %f \n",i,c_k[i]);


			//calculating koeff's of decomposition of Wy0 onto p_k
			double gamma_k = 0.0;
			for (int k = 1; k <= Ny; k++)
			{
				gamma_k = func_gamma_k(k,Ny);
				//if (harm == 1)
				//	printf("gamma_%d = %f \n",k,gamma_k);
				c_k_new[k-1] = - c_k[k-1] * gamma_k * (1.0 / hy(0)) * (6.0 / (6.0 - gamma_k * gamma_k)); //minus sign because of B = -grad for me and B = grad for Popov
			}
			//if (harm == 4)
			//	for (int i = 0; i < Ny; i++)
			//		printf("c_k_new[%d] = %f \n",i,c_k_new[i]);

			//calculating Wy_0
			for (int i = 0; i < Ny + 1; i++)
				Wy_0_eigenv[i] = 0.0;
			for (int k = 1; k <= Ny; k++)
			{
				func_p_k(p_k, k, Ny);
				for (int i = 0; i < Ny + 1; i++)
					Wy_0_eigenv[i] += c_k_new[k-1] * p_k[i];
			}

			for (int i = 0; i < Ny + 1; i++)
				Wy_0_error_eigenv[i] = Wy_0_eigenv[i] - (-exact_gradient_harm(harm, ypoints[i]));

#ifdef DEBUGE
//			if (harm == Nx)
			{
				char filename[30];
				char* ext = ".xls";
				char nx_string[5];
				char harm_string[5];
				sprintf(nx_string,"_%d",Nx);
				sprintf(harm_string,"_%d",harm);
				strcpy(filename,"Harm_error");
				strcat(filename,nx_string);
				strcat(filename,harm_string);
				strcat(filename,ext);
				FILE * Harm_error_file = fopen(filename,"wt");
				for (int ind = 0; ind < Ny + 1; ind++)
					fprintf(Harm_error_file,"%f \t %f \t %f \n",Wy_0_eigenv[ind], -exact_gradient_harm(harm, ypoints[ind]), Wy_0_error_eigenv[ind]);
				fclose(Harm_error_file);

			}
#endif
			//if (harm == 4)
			//	for (int i = 0; i < Ny + 1; i++)
			//		printf("Wy_0_error_eigenv[%d] = %f \n",i,Wy_0_error_eigenv[i]);

			double div_harm_l2 = 0.0, div_harm_max = 0.0;
			divergence_Wy_norm(Wy_0_error_eigenv, Ny + 1, &div_harm_l2, &div_harm_max);
			diver_harmonics_l2[harm - 1] = div_harm_l2;
			diver_harmonics_max[harm - 1] = div_harm_max;

			//printf("diver_l2[%d] = %f \n",harm - 1, diver_harmonics_l2[harm - 1]);
		}
#ifdef DEBUGE
		{
			char* ext = ".xls";
			char nx_string[5];
			sprintf(nx_string,"_%d",Nx);
			char filename[30];
			strcpy(filename,"Harm_l2");
			strcat(filename,nx_string);
			strcat(filename,ext);
			FILE * Harm_l2_file = fopen(filename,"wt");
			for (int ind = 0; ind < Ny + 1; ind++)
				fprintf(Harm_l2_file,"%f \n",diver_harmonics_l2[ind]);
			fclose(Harm_l2_file);

			char filename2[30];
			strcpy(filename2,"Harm_max");
			strcat(filename2,nx_string);
			strcat(filename2,ext);
			FILE * Harm_max_file = fopen(filename2,"wt");
			for (int ind = 0; ind < Ny + 1; ind++)
				fprintf(Harm_max_file,"%f \n",diver_harmonics_max[ind]);
			fclose(Harm_max_file);

		}
#endif

		double * Psi_h = (double*)malloc((Ny + 1) * sizeof(double));
		double * Eps_h = (double*)malloc((Ny + 1) * sizeof(double));
		double * B_tr_Eps_h = (double*)malloc(Ny* sizeof(double));
		double * c_k_Psi = (double*)malloc((Ny + 1) * sizeof(double));
		double * d_k_Eps = (double*)malloc((Ny + 1) * sizeof(double));
		double * beta_k_divEps = (double*)malloc(Ny * sizeof(double));

		//filling in T0 and Psi_h
		double temper_i, awy0_exact_i, diff2;
		for (int i = 0; i <=  Ny; i++)
		{
			T0[i] = T[i*Nx];
			if (i == 0)
			{
				a_0_old = hy(0)/(heat_conductivity(0,0,0)*6.0);
				b_0_old = 2*a_0_old;

				temper_i = Ty_0_temp - T0[0];
				awy0_exact_i = b_0_old * Wy0_exact[0] + a_0_old * Wy0_exact[1];

				printf("i = %d \n",i);
				printf("temper_0 = %f \n",temper_i);
				printf("awy0_exact_0 = %f \n",awy0_exact_i);

			}
			else
				if (i == Ny)
				{
					a_nminus1_old = hy(Ny-1)/(heat_conductivity(0,0,0)*6.0);
					b_n_old = 2*a_nminus1_old;

					temper_i = T0[Ny-1] - Ty_1_temp;
					awy0_exact_i = a_nminus1_old * Wy0_exact[Ny-1] + b_n_old * Wy0_exact[Ny];
				}
				else
				{
					a_i_old = hy(i)/(heat_conductivity(0,0,0)*6.0);
					a_iminus1_old = hy(i-1)/(heat_conductivity(0,0,0)*6.0);
					b_i_old = 2*a_iminus1_old + 2*a_i_old;

					temper_i = T0[i - 1] - T0[i];
					awy0_exact_i = a_iminus1_old * Wy0_exact[i-1] + b_i_old * Wy0_exact[i] + a_i_old * Wy0_exact[i + 1];
				}
			diff2 = temper_i - awy0_exact_i;
			if (i == 0)
				printf("diff2 = %f \n", diff2);
			Psi_h[i] = diff2;
		}
		func_p_k(Psi_h, 0, Ny);
		func_p_k(p_k, 2, Ny);
		double result = dot_product(Psi_h, p_k, Ny + 1);
		printf("RESULT = %f \n", result);
		FILE* Psi_init_file = fopen("Psi_init.xls","wt");
		for (int i = 0; i <= Ny; i++)
			fprintf(Psi_init_file,"%f \n",Psi_h[i]);
		fclose(Psi_init_file);


		//func_u_k(T0, 4, Ny);

		//for (int i = 0; i < Ny; i++)
		//	printf("T0[%d] = %f \n",i,T0[i]);

		//decomposing T_0 onto span of u_k = calculating c_k's;
		for (int k = 1; k <= Ny; k++)
		{
			func_u_k(u_k, k, Ny);
			c_k[k-1] = dot_product(T0, u_k, Ny) / dot_product(u_k, u_k, Ny);
		}
		FILE* ck_coeffs_file = fopen("ck_koeffs.xls","wt");
		for (int i = 0; i < Ny; i++)
			fprintf(ck_coeffs_file,"%f \n",c_k[i]);
		fclose(ck_coeffs_file);
		//for (int i = 0; i < Ny; i++)
		//	printf("c_k[%d] = %f \n",i,c_k[i]);

		//calculating koeff's of decomposition of Wy0 onto p_k
		double gamma_k = 0.0;
		for (int k = 1; k <= Ny; k++)
		{
			gamma_k = func_gamma_k(k,Ny);
			c_k_new[k-1] = - c_k[k-1] * gamma_k * (1.0 / hy(0)) * (6.0 / (6.0 - gamma_k * gamma_k)); //minus sign because of B = -grad for me and B = grad for Popov
			c_k_div[k-1] = c_k_new[k-1] * gamma_k / hy(0);
		}

		FILE* ck_div_coeffs_file = fopen("ck_div_koeffs.xls","wt");
		for (int i = 0; i < Ny; i++)
			fprintf(ck_div_coeffs_file,"%f \n",c_k_div[i]);
		fclose(ck_div_coeffs_file);

		//for (int i = 0; i < Ny; i++)
		//	printf("c_k_new[%d] = %f \n",i,c_k_new[i]);

		//calculating Wy_0
		for (int i = 0; i < Ny + 1; i++)
			Wy_0_eigenv[i] = 0.0;
		for (int k = 1; k <= Ny; k++)
		{
			func_p_k(p_k, k, Ny);
			for (int i = 0; i < Ny + 1; i++)
				Wy_0_eigenv[i] += c_k_new[k-1] * p_k[i];
		}

		//calculating B(tr) Wy_0
		for (int i = 0; i < Ny; i++)
			Btr_Wy_0_eigenv[i] = 0.0;
		for (int k = 1; k <= Ny; k++)
		{
			func_u_k(u_k, k, Ny);
			for (int i = 0; i < Ny; i++)
				Btr_Wy_0_eigenv[i] += c_k_div[k-1] * u_k[i];
		}


		//decomposing Psi_h onto span of p_k = calculating c_k_Psi's; 
		for (int k = 0; k <= Ny; k++)
		{
			func_p_k(p_k, k, Ny);
			c_k_Psi[k] = dot_product(Psi_h, p_k, Ny + 1) / dot_product(p_k, p_k, Ny + 1);
		}
		//НЕ РАБОТАЕТ!!! СЧИТАЕТ НЕВЕРНО
		FILE* ck_Psi_coeffs_file = fopen("ck_Psi_koeffs.xls","wt");
		for (int i = 0; i <= Ny; i++)
			fprintf(ck_Psi_coeffs_file,"%f \n",c_k_Psi[i]);
		fclose(ck_Psi_coeffs_file);
		//calculating Psi_h
		for (int i = 0; i <= Ny; i++)
			Psi_h[i] = 0.0;
		for (int k = 0; k <= Ny; k++)
		{
			func_p_k(p_k, k, Ny);
			for (int i = 0; i <= Ny; i++)
				Psi_h[i] += c_k_Psi[k] * p_k[i];
		}
		FILE* Psi_file = fopen("Psi.xls","wt");
		for (int i = 0; i <= Ny; i++)
			fprintf(Psi_file,"%f \n",Psi_h[i]);
		fclose(Psi_file);

		//double gamma_k = 0.0;
		d_k_Eps[0] = 0;
		for (int k = 1; k <= Ny; k++)
		{
			gamma_k = func_gamma_k(k,Ny);
			d_k_Eps[k] = - c_k_Psi[k] * gamma_k * (1.0 / hy(0)) * (6.0 / (6.0 - gamma_k * gamma_k)); //minus sign because of B = -grad for me and B = grad for Popov
		}
		FILE* dk_Eps_coeffs_file = fopen("dk_Eps_koeffs.xls","wt");
		for (int i = 0; i <= Ny; i++)
			fprintf(dk_Eps_coeffs_file,"%f \n",d_k_Eps[i]);
		fclose(dk_Eps_coeffs_file);
		//calculating Eps_h
		for (int i = 0; i < Ny; i++)
			Eps_h[i] = 0.0;
		for (int k = 0; k <= Ny; k++)
		{
			func_p_k(p_k, k, Ny);
			for (int i = 0; i <= Ny; i++)
				Eps_h[i] += d_k_Eps[k] * p_k[i];
		}
		FILE* Eps_file = fopen("Eps.xls","wt");
		for (int i = 0; i <= Ny; i++)
			fprintf(Eps_file,"%f \n",Eps_h[i]);
		fclose(Eps_file);


		for (int k = 1; k <= Ny; k++)
		{
			gamma_k = func_gamma_k(k,Ny);
			beta_k_divEps[k - 1] = d_k_Eps[k - 1] * gamma_k * (1.0 / hy(0));
		}
		FILE* betak_divEps_coeffs_file = fopen("betak_divEps_koeffs.xls","wt");
		for (int i = 0; i < Ny; i++)
			fprintf(betak_divEps_coeffs_file,"%f \n",beta_k_divEps[i]);
		fclose(betak_divEps_coeffs_file);
		
		//calculating B_tr_Eps_h
		for (int i = 0; i < Ny; i++)
			B_tr_Eps_h[i] = 0.0;
		for (int k = 1; k <= Ny; k++)
		{
			func_u_k(u_k, k, Ny);
			for (int i = 0; i < Ny; i++)
				B_tr_Eps_h[i] += beta_k_divEps[k - 1] * u_k[i];
		}
		FILE* B_tr_Eps_file = fopen("B_tr_Eps.xls","wt");
		for (int i = 0; i < Ny; i++)
			fprintf(B_tr_Eps_file,"%f \n",B_tr_Eps_h[i]);
		fclose(B_tr_Eps_file);

		free(Psi_h);
		free(Eps_h);
		free(B_tr_Eps_h);
		free(c_k_Psi);
		free(d_k_Eps);
		free(beta_k_divEps);


#ifdef DEBUGE
		char* ext = ".xls";
		char nx_string[5];
		sprintf(nx_string,"_%d",Nx);
		char filename[30];
		strcpy(filename,"Wy_0_eigenv");
		strcat(filename,nx_string);
		strcat(filename,ext);
		FILE * Wy_eigenv_file = fopen(filename,"wt");
		for (int ind = 0; ind < Ny + 1; ind++)
			fprintf(Wy_eigenv_file,"%f \t %f \n",Wy_0_eigenv[ind], Wy_n[ind]);
		//FilePrint2D("Wy_lumping_error_2D.xls",Wy_0_lumping_error, Ny + 1, Nx);
		fclose(Wy_eigenv_file);

		char filename2[30];
		strcpy(filename2,"Btr_Wy_0_eigenv");
		strcat(filename2,nx_string);
		strcat(filename2,ext);
		FILE * Btr_Wy_eigenv_file = fopen(filename2,"wt");
		double temp1, temp2;
		for (int ind = 0; ind < Ny; ind++)
		{
			temp1 = -exact_gradientY(numSolution,0.5*hx(0),ypoints[ind],0,0);
			temp2 = -exact_gradientY(numSolution,0.5*hx(0),ypoints[ind + 1],0,0);
			fprintf(Btr_Wy_eigenv_file,"%f \t %f \t %f \t \t %f \n",Btr_Wy_0_eigenv[ind], (Wy_n[ind + 1] - Wy_n[ind])/hy(0), (temp2 - temp1) / hy(0), Btr_Wy_0_eigenv[ind] + ((temp2 - temp1) / hy(0)) );
		}
		fclose(Btr_Wy_eigenv_file);

#endif

	}
#endif //#endif EIGENVECTOR_WY0

#ifdef SADDLE_POINT
	MKL_INT nnz_dir = ( (2 + 1) + (2 + 1) + (Ny+1 - 2)*5) + (2 * Ny);
	MKL_INT size = (Ny + 1) + Ny;
	int *row = (int*)malloc(nnz_dir*sizeof(int));
	int *col = (int*)malloc(nnz_dir*sizeof(int));
	double *acoo = (double*)malloc(nnz_dir*sizeof(double));
	int count = 0;
	for (int row_ind = 0; row_ind < Ny + 1; row_ind++)
	{
		if (row_ind > 0)
		{
			row[count] = row_ind;
			col[count] = row_ind - 1;
			acoo[count] = hx(0)*hy(0)*1.0/6.0;
			count++;
		}

		row[count] = row_ind;
		col[count] = row_ind;
		if (row_ind == 0 || row_ind == Ny)
			acoo[count] = 2.0*hx(0)*hy(0)*1.0/6.0;
		else
			acoo[count] = 4.0*hx(0)*hy(0)*1.0/6.0;
		count++;

		if (row_ind < Ny)
		{
			row[count] = row_ind;
			col[count] = row_ind + 1;
			acoo[count] = hx(0)*hy(0)*1.0/6.0;
			count++;
		}

		if (row_ind > 0)
		{
			row[count] = row_ind;
			col[count] = (Ny + 1) + row_ind - 1;
			acoo[count] = - hx(0);  //Aw - BT = 0, B = -gradient -> Aw + gradient T = 0
			count++;
		}

		if (row_ind < Ny)
		{
			row[count] = row_ind;
			col[count] = (Ny + 1) + row_ind;
			acoo[count] = hx(0);  //Aw - BT = 0, B = -gradient -> Aw + gradient T = 0
			count++;
		}
	}

	for (int row_ind = 0; row_ind < Ny; row_ind++)
	{
		row[count] = Ny + 1 + row_ind;
		col[count] = row_ind;
		acoo[count] = -hx(0);  //B(tr) w = 0
		count++;

		row[count] = Ny + 1 + row_ind;
		col[count] = row_ind + 1;
		acoo[count] =  hx(0);  //B(tr) w = 0
		count++;
	}

	if (count != nnz_dir){
		printf("FUCK YOU: count=%d != nnz_dir=%d \n", count, nnz_dir);return -1;}

	MKL_INT* job = (MKL_INT*)malloc(6*sizeof(MKL_INT));
	job[0] = 1;
	job[1] = 1;
	job[2] = 0;
	job[4] = nnz_dir;
	job[5] = 0;
	MKL_INT info = 0;
	MKL_INT * ia = (MKL_INT*)malloc((size)*sizeof(MKL_INT));
	MKL_INT * ja = (MKL_INT*)malloc(nnz_dir*sizeof(MKL_INT));
	double * acsr = (double*)malloc(nnz_dir*sizeof(double));
	mkl_dcsrcoo (job, &size, acsr, ja, ia, &nnz_dir, acoo, row, col, &info);

	if (info != 0){
		printf("Error in calling mkl_dcsrcoo, info = %d \n", info); return -1;}
	else
		printf("Success in mkl_dcsrcoo, info = %d \n", info);

	double *rhand = (double*)malloc(size * sizeof(double));
	double *solution = (double*)malloc(size * sizeof(double));
	double * saddle_w_solution = (double*)malloc((Ny + 1) * sizeof(double));
	double * saddle_T_solution = (double*)malloc(Ny * sizeof(double));
	for (int j = 0; j < Ny + 1; j++)
		rhand[j] = 0.0;
	int index = 0;

	//double *B_tr = (double*)malloc(Ny * sizeof(double));
	//for (int j = 0; j < Ny; j++)
	//{
	//	double m_i = heat_capacity(0,j,0)*density(0,j,0)*hy(j);

	//	//B_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i) * (wy[i*(Ny+1) + j + 1] - wy[i*(Ny+1) + j]) + hy(j) * (wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]) );
	//	B_tr[j] = hx(0) * (Wy_n[j + 1] - Wy_n[j]) ;
	//	//B_tr[j] = (1.0/m_i)* (Wy_n[j + 1] - Wy_n[j]) ;
	//}

	for (int j = 0; j < Ny; j++)
	{
		index = (Ny + 1)  + j;
		rhand[index] = -saddle_righthand_laplace(numSolution, 0, j, mx, my);
		//rhand[index] = B_tr[j];
		//rhand[index] = hx(0) * B_tr[j];
		//rhand[index] = hx(0) * hy(0);
	}

	void* pt[64] ;
	MKL_INT maxfct ;
	MKL_INT mnum ;
	MKL_INT mtype ; 
	MKL_INT phase ;
	MKL_INT* perm;
	MKL_INT nrhs ;
	MKL_INT error ;
	MKL_INT msglvl ; 
	MKL_INT iparm[64];
	MKL_INT idum;
	//MKL_INT n1;

	for (MKL_INT i=0; i< 64; i++)
	{
		pt[i] = 0;
		iparm[i] = 0;
	}
	maxfct = 1;
	mnum = 1;
	mtype = 11; //-2 = real and symmetric indefinite matrix, 11 = real and unsymmetric matrix
	phase = 11;
	//n1 = dim_Matrix;
	nrhs = 1;
	error = 100;
	msglvl = 0; 
	//	iparm[0]= 0;
	iparm[0] = 1; // No solver default 
	iparm[1] = 2; // Fill-in reordering from METIS 
	iparm[2] = 1; // Not in use
	//iparm[3] = 80; // Preconditioned CGS
	iparm[3] = 0; // Preconditioned CGS
	iparm[4] = 0; // User permutation 
	iparm[5] = 0; // Write solution on x
	iparm[6] = 0; // Not in use 
	//iparm[7] = 150; // Max numbers of iterative refinement steps 
	iparm[7] = 0; // Max numbers of iterative refinement steps 
	iparm[8] = 0; // Not in use 
	iparm[9] = 13; // Perturb the pivot elements with 1E-13 
	iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS = scaling vectors. 
	iparm[11] = 0; // Not in use 
	iparm[12] = 1; // Input = improved accuracy using (non-)symmetric weighted matchings. 
	iparm[13] = 0; // Output = number of perturbed pivots
	iparm[14] = 0; // Output = peak memory symbolic factorization
	iparm[15] = 0; // Output = permanent memory symbolic	factorization
	iparm[16] = 0; // Output = memory numerical factorization and solution 
	iparm[17] = -1; // Output: Number of nonzeros in the factor LU /
	iparm[18] = -1; // Output: Mflops for LU factorization 
	iparm[19] = 0; // Output: Numbers of CG Iterations 
	iparm[26] = 1; // проверяет структуру матрицы
	iparm[27] = 0; //0 = двойная точность

	printf("\nBefore solving...\n");
	pardiso (pt, &maxfct, &mnum, &mtype, &phase, &size, acsr, ia, ja, &idum, &nrhs, iparm, &msglvl, rhand, solution, &error);
	printf("Stage 1 completed. \nerror = %d \n",error);
	phase = 22;
	pardiso (pt, &maxfct, &mnum, &mtype, &phase, &size, acsr, ia, ja, &idum, &nrhs, iparm, &msglvl, rhand, solution, &error);
	printf("Stage 2 completed. \nerror = %d \n",error);
	phase = 33;
	pardiso (pt, &maxfct, &mnum, &mtype, &phase, &size, acsr, ia, ja, &idum, &nrhs, iparm, &msglvl, rhand, solution, &error);
	printf("Stage 3 completed. \nerror = %d \n",error);

	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	else
	{
		printf("Success in pardiso, info = %d \n", error);
		double *residual =(double*)malloc(size * sizeof(double));
		char *transa = "N";
		mkl_dcsrgemv (transa, &size, acsr, ia, ja, solution, residual);
		for (int i = 0; i < size; i++)
			residual[i] -= rhand[i];
		printf("residual_max \t %f \n",norm_max(residual,size, 1));
		printf("residual_l2 \t %f \n",norm_l2(residual,size, 1));


		for (int i = 0; i < Ny + 1; i++)
			saddle_w_solution[i] = solution[i];
		for (int i = 0; i < Ny; i++)
			saddle_T_solution[i] = solution[(Ny + 1) + i];

		double * saddle_w_exact = (double*)malloc((Ny + 1) * sizeof(double));
		double * saddle_T_exact = (double*)malloc(Ny * sizeof(double));
		for (int i = 0; i < Ny + 1; i++)
		{
			saddle_w_exact[i] = (-exact_gradientY(numSolution, 0, ypoints[i], 0, 0));
			//saddle_w_exact[i] = (1.0 - 2*ypoints[i]) / 2.0;
			//saddle_w_exact[i] = 0.0;
		}
		for (int i = 0; i < Ny; i++)
		{
			saddle_T_exact[i] = exact_solution(numSolution, 0, ypoints[i] + 0.5*hy(i), 0, 0);
			//saddle_T_exact[i] = (ypoints[i] + 0.5*hy(i))*(1.0 - (ypoints[i] + 0.5*hy(i))) / 2.0;
			//saddle_T_exact[i] = (1.0 - 2*(ypoints[i] + 0.5*hy(i))) / 2.0;
		}


		double * saddle_exactsolution = (double*)malloc(size * sizeof(double));
		for (int i = 0; i < Ny + 1; i++)
		{
			saddle_exactsolution[i] = saddle_w_exact[i];
			residual[i] = 0;
		}
		for (int i = 0; i < Ny; i++)
		{
			saddle_exactsolution[(Ny + 1) + i] = saddle_T_exact[i];
			residual[(Ny + 1) + i] = 0;
		}
		FILE * Sol_exact = fopen("Sol_exact.xls","wt");
		for (int i = 0; i < size; i++)
			fprintf(Sol_exact,"%f \n", saddle_exactsolution[i]);
		fclose(Sol_exact);

		mkl_dcsrgemv (transa, &size, acsr, ia, ja, saddle_exactsolution, residual);
		//for (int i = 0; i < nnz_dir; i++)
		//{
		//	if (row[i] == Ny || row[i] == 0)
		//	{
		//		printf("i = %d row[i] = %d col[i] = %d acoo[i] = %f x[col[i]] = %f \n",i,row[i],col[i],acoo[i],saddle_exactsolution[col[i]]);

		//	}
		//}

		//mkl_dcoogemv (transa, &size, acoo, row, col, &nnz_dir, saddle_exactsolution, residual);
		FILE * Ax_exact = fopen("Ax_exact.xls","wt");
		for (int i = 0; i < size; i++)
			fprintf(Ax_exact,"%f \n", residual[i]);
		fclose(Ax_exact);
		FILE * Rhand_exact = fopen("Rhand_exact.xls","wt");
		for (int i = 0; i < size; i++)
			fprintf(Rhand_exact,"%f \n", rhand[i]);
		fclose(Rhand_exact);
		for (int i = 0; i < size; i++)
			residual[i] -= rhand[i];

		char nx_string[10];
		sprintf(nx_string,"_%d.xls",Nx);
		char filename0[30];
		strcpy(filename0,"saddle_approx_residual");
		strcat(filename0, nx_string);
		FILE * saddle_approx_file = fopen(filename0,"wt");
		for (int i = 0; i < size; i++)
			fprintf(saddle_approx_file,"%f \n",residual[i]);
		fclose(saddle_approx_file);


		double * saddle_w_error = (double*)malloc((Ny + 1) * sizeof(double));
		double * saddle_T_error = (double*)malloc(Ny * sizeof(double));
		for (int i = 0; i < Ny + 1; i++)
			saddle_w_error[i] = saddle_w_solution[i] - saddle_w_exact[i];
		for (int i = 0; i < Ny; i++)
			saddle_T_error[i] = saddle_T_solution[i] - saddle_T_exact[i];

		char filename1[30], filename2[30];
		strcpy(filename1,"saddle_w");
		strcat(filename1, nx_string);
		FILE * file1 = fopen(filename1,"wt");
		for (int i = 0; i < Ny + 1; i++)
			fprintf(file1,"%f \t %f \t %f \n",saddle_w_solution[i], saddle_w_exact[i], saddle_w_error[i]);
		fprintf(file1,"eps_max \t %f \n",norm_max(saddle_w_error,Ny + 1, 1));
		fprintf(file1,"eps_l2 \t %f \n",norm_l2(saddle_w_error,Ny + 1, 1));
		fclose(file1);

		strcpy(filename2,"saddle_T");
		strcat(filename2, nx_string);
		FILE * file2 = fopen(filename2,"wt");
		for (int i = 0; i < Ny; i++)
			fprintf(file2,"%f \t %f \t %f \n",saddle_T_solution[i], saddle_T_exact[i], saddle_T_error[i]);
		fprintf(file2,"eps_max \t %f \n",norm_max(saddle_T_error,Ny, 1));
		fprintf(file2,"eps_l2 \t %f \n",norm_l2(saddle_T_error,Ny, 1));
		fclose(file2);

		double div_saddleW_l2 = 0.0, div_saddleW_max = 0.0;
		divergence_Wy_norm(saddle_w_error, Ny, &div_saddleW_l2, &div_saddleW_max);
		printf("div_saddleW_l2 = %f \n",div_saddleW_l2);
		printf("div_saddleW_max = %f \n",div_saddleW_max);

	}

#endif //#ifdef SADDLE_POINT

#ifdef DEBUGE
	//печать массивов Wx_0, Wy_0, Wz_0. 
	//печать массивов Wx_0, Wy_0, Wz_0. первый столбец - посчитанные данные, второй - точные значения для
	//тестового решения numSolution 1111.
	switch(printcase)
	{
	case 0:
		break;
	case 1:
		{
			double add1 = 0;

			//FILE* wz_check = fopen("Wz0_check_split3D.xls","wt");
			//for (int i = 0; i < dim_Wz; i++)
			//{
			//	k1 = i % (Nz+1);
			//	j1 = i / ((Nz+1)*Nx);
			//	i1 = (i - k1 - j1*(Nz+1)*Nx) / (Nz+1);

			//     add1 = -heat_conductivity(i1,j1,k1)*exact_gradientZ(numSolution,xpoints[i1] + 0.5*hx(i1), ypoints[j1] + 0.5*hy(j1), zpoints[k1],0);

			//	fprintf(wz_check,"%15.15f \t %15.15f \n",Wz_n[i], add1);
			//}
			//fclose(wz_check);

			FILE* wy_check = fopen("Wy0_check_split3D.xls","wt");
			for (int i = 0; i < dim_Wy; i++)
			{
				j1 = i % (Ny+1);
				i1 = i / ((Ny+1)*Nz);
				k1 = (i - j1 - i1*(Ny+1)*Nz) / (Ny+1);

				add1 = -heat_conductivity_func(xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1] + 0.5*hx(i1), ypoints[j1], zpoints[k1] + 0.5*hz(k1),0);

				fprintf(wy_check,"%15.15f \t %15.15f \n",Wy_n[i], add1);
			}
			fclose(wy_check);

			FILE* wx_check = fopen("Wx0_check_split3D.xls","wt");
			for (int i = 0; i < dim_Wx; i++)
			{
				i1 = i % (Nx+1);
				k1 = i / ((Nx+1)*Ny);
				j1 = (i - i1 - k1*(Nx+1)*Ny) / (Nx+1);

				add1 = -heat_conductivity_func(xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[i1], ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1),0);

				fprintf(wx_check,"%15.15f \t %15.15f \n",Wx_n[i], add1);
			}
			fclose(wx_check);
		}
		break;
	}
#endif


#ifdef EXACT_W0_INIT
	for (int i = 0; i <= Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
#ifdef STEKLOV_FLUX
			if ( numSolution != 300 )
			{
				printf ( "ERROR: numsol must be 300 for STEKLOV_FLUX \n" );
				return -1;
			}
			Wx_n[j*(Nx+1) + i] = - exact_gradient_average_300x(i, j, 0, tau, mx, my);
			//printf ( "val1 = %f val2 = %f \n", - exact_gradient_average_300x(i, j, 0, tau, mx, my), - exact_gradientX(numSolution,xpoints[i], ypoints[j] + 0.5*hy(j), 0, 0) );
			//_getch();
#else
			Wx_n[j*(Nx+1) + i] = - exact_gradientX(numSolution,xpoints[i], ypoints[j] + 0.5*hy(j), 0, 0);
#endif
		}

	FILE* wx_check = fopen("Wx0_to4noe.xls","wt");
	for (int i = 0; i < dim_Wx; i++)
		fprintf(wx_check,"%15.15f \n",Wx_n[i]);
	fclose(wx_check);

	for (int i = 0; i < Nx; i++)
		for (int j = 0; j <= Ny; j++)
		{
#ifdef STEKLOV_FLUX
			if ( numSolution != 300 )
			{
				printf ( "ERROR: numsol must be 300 for STEKLOV_FLUX \n" );
				return -1;
			}
			Wy_n[i*(Ny+1) + j] = - exact_gradient_average_300y(i, j, 0, tau, mx, my);
#else
			Wy_n[i*(Ny+1) + j] = - exact_gradientY(numSolution,xpoints[i]+0.5*hx(i), ypoints[j], 0, 0);
#endif
		}

	FILE* wy_check = fopen("Wy0_to4noe.xls","wt");
	for (int i = 0; i < dim_Wy; i++)
		fprintf(wy_check,"%15.15f \n",Wy_n[i]);
	fclose(wy_check);

#endif

	{;}

	{;}
	double eps0_max_w = 0.0;
	double eps0_max_T = 0.0;
	double eps0_l2_w = 0.0;
	double eps0_l2_T = 0.0;
	double eps0_relative_max_w = 0.0;
	double eps0_relative_max_T = 0.0;
	double eps0_relative_l2_w = 0.0;
	double eps0_relative_l2_T = 0.0;
	Accuracy_calculate(T, Wx_n, Wy_n, numSolution, -1, Nx, Ny, Nz, print_step, &eps0_max_w, &eps0_max_T, &eps0_l2_w, &eps0_l2_T, &eps0_relative_max_w, &eps0_relative_max_T, &eps0_relative_l2_w, &eps0_relative_l2_T);
/*
	double norm0_l2_wx = 0;
	double norm0_l2_wy = 0;
	double norm0_l2_wz = 0;
	double norm0_l2_w = 0;
	double norm0_max_w = 0;

	double diffx0 = 0;
	double epsx0 = 0;
	double epsx0_l2 = 0;
	int imax0_x = 0;
	int jmax0_x = 0;
	int kmax0_x = 0;
	double temprr = 0;
	for (int i = 0; i < Ny*Nz; i++)
		for (int k = 0; k < Nx+1; k++)
		{
			i1 = k;
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			diffx0 = fabs(Wx_n[i*(Nx+1) + k] - (-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution,xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) ));
			temprr = fabs(-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution, xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) );

			if (diffx0 > epsx0)
			{
				epsx0 = diffx0;
				imax0_x = k;
				jmax0_x = j1;
				kmax0_x = k1;
			}
			if (temprr > norm0_max_w)
			{
				norm0_max_w = temprr;
			}
			norm0_l2_wx += temprr*temprr;
			epsx0_l2 += diffx0*diffx0;
		}
		epsx0_l2 = sqrt(epsx0_l2 /(dim_Wx));
		norm0_l2_wx = sqrt(norm0_l2_wx / (dim_Wx));

		////////////////
		double diffy0 = 0;
		double epsy0 = 0;
		double epsy0_l2 = 0;
		int imax0_y = 0;
		int jmax0_y = 0;
		int kmax0_y = 0;

		for (int i=0; i<Nz*Nx; i++)
			for (int k=0; k<Ny+1; k++)
			{
				i1 = i/Nz;
				j1 = k;
				k1 = i - (i/Nz)*Nz;

				diffy0 = fabs(Wy_n[i*(Ny+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) ));
				temprr = fabs(-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0));

				if (diffy0 > epsy0)
				{
					epsy0 = diffy0;
					imax0_y = i1;
					jmax0_y = j1;
					kmax0_y = k1;
				}
				if (temprr > norm0_max_w)
				{
					norm0_max_w = temprr;
				}
				norm0_l2_wy += temprr*temprr;
				epsy0_l2 += diffy0*diffy0;
			}
			epsy0_l2 = sqrt(epsy0_l2 /(dim_Wy));
			norm0_l2_wy = sqrt(norm0_l2_wy / (dim_Wy));

			double diffz0 = 0;
			double epsz0 = 0;
			double epsz0_l2 = 0;
			int imax0_z = 0;
			int jmax0_z = 0;
			int kmax0_z = 0;
			for (int i=0; i<Ny*Nx; i++)
				for (int k=0; k<Nz+1; k++)
				{
					i1 = i%Nx;
					j1 = (i-i%Nx)/Nx;
					k1 = k;
					switch(numSolution)
					{
					default:
						diffz0 = 0;
						temprr = 0;
						break;
					}
					if (diffz0 > epsz0)
					{
						epsz0 = diffz0;
						imax0_z = i1;
						jmax0_z = j1;
						kmax0_z = k1;
					}
					if (temprr > norm0_max_w)
						norm0_max_w = temprr;

					norm0_l2_wz += temprr*temprr;
					epsz0_l2 += diffz0*diffz0;
				}
				epsz0_l2 = 0;
				norm0_l2_wz = 0;
				norm0_l2_w = sqrt(norm0_l2_wx*norm0_l2_wx + norm0_l2_wy*norm0_l2_wy);

				eps0_l2_w = sqrt(epsx0_l2*epsx0_l2 + epsy0_l2*epsy0_l2);
				eps0_relative_l2_w = eps0_l2_w / norm0_l2_w;
				eps0_max_w = 0;
				if (epsx0_l2 > eps0_max_w)
					eps0_max_w = epsx0_l2;
				if (epsy0_l2 > eps0_max_w)
					eps0_max_w = epsy0_l2;
				if (epsz0_l2 > eps0_max_w)
					eps0_max_w = epsz0_l2;
				eps0_relative_max_w = eps0_max_w / norm0_max_w;

				printf("epsx0 = %f \nimax0_x = %d, jmax0_x = %d, kmax0_x = %d \n",epsx0,imax0_x,jmax0_x,kmax0_x);
				printf("epsx0 = %e \n",epsx0);
				printf("epsx0_l2 = %e \n",epsx0_l2);
				printf("epsy0 = %f \nimax0_y = %d, jmax0_y = %d, kmax0_y = %d \n",epsy0,imax0_y,jmax0_y,kmax0_y);
				printf("epsy0 = %e \n",epsy0);
				printf("epsy0_l2 = %e \n",epsy0_l2);
				printf("epsz0 = %f \nimax0_z = %d, jmax0_z = %d, kmax0_z = %d \n",epsz0,imax0_z,jmax0_z,kmax0_z);
				printf("epsz0 = %e \n",epsz0);
				printf("epsz0_l2 = %e \n",epsz0_l2);
				printf("norm0_max_w = %e \n",norm0_max_w);
				printf("norm0_l2_w = %e \n",norm0_l2_w);
*/
				printf("eps0_l2_w = %e \n",eps0_max_w);
				printf("eps0_l2_w = %e \n",eps0_l2_w);
				printf("eps0relative_max_w = %e \n",eps0_relative_max_w);
				printf("eps0relative_l2_2 = %e \n",eps0_relative_l2_w);

				double *Wx_0_error = (double*)malloc(dim_Wx*sizeof(double));
				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_error[i*(Nx+1) + k] = Wx_n[i*(Nx+1) + k] - (-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution,xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) );
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

						Wy_0_error[i*(Ny+1) + k] = Wy_n[i*(Ny+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}

				double *Wx_0_ex = (double*)malloc(dim_Wx*sizeof(double));
				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_ex[i*(Nx+1) + k] = - (-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution,xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}
				double *Wy_0_ex = (double*)malloc(dim_Wy*sizeof(double));
				for (int i=0; i<Nz*Nx; i++)
				{
					for (int k=0; k<Ny+1; k++)
					{
						i1 = i/Nz;
						j1 = k;
						k1 = i - (i/Nz)*Nz;

						Wy_0_ex[i*(Ny+1) + k] = - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}
				printf("l2-norm of Wx_0_error = %f \n", norm_l2(Wx_0_error,Nx+1,Ny));
				printf("l2-norm of Wy_0_error = %f \n", norm_l2(Wy_0_error,Ny+1,Nx));
				printf("max-norm of Wx_0_error = %f \n", norm_max(Wx_0_error,Nx+1,Ny));
				printf("max-norm of Wy_0_error = %f \n", norm_max(Wy_0_error,Ny+1,Nx));

				double l2_norm_div, max_norm_div;

				//for (int i = 0; i < dim_Wx; i++)
				//	Wx_0_error[i] = 1.0;
				//for (int i = 0; i < dim_Wy; i++)
				//	Wy_0_error[i] = 1.0;

				double div_W1D_l2 = 0.0, div_W1D_max = 0.0;
				//divergence_Wy_norm(Wy_n, Ny, &div_W1D_l2, &div_W1D_max);
				//printf("div_W1D_l2 = %f \n",div_W1D_l2);
				//printf("div_W1D_max = %f \n",div_W1D_max);
				//divergence_norm(Wx_n,Wy_n,Nx,Ny,&l2_norm_div, &max_norm_div);
				//printf("div l2-norm of W_0 = %f \n", l2_norm_div);
				//printf("div max-norm of W_0 = %f \n", max_norm_div);
				//fprintf(f1,"div l2_norm = %e \n", l2_norm_div);
				//fprintf(f1,"div max_norm = %e \n", max_norm_div);

				//double l2_norm_stab, max_norm_stab;
				//stability_norm(Wx_n,Wy_n,Nx,Ny,&l2_norm_stab, &max_norm_stab);
				//printf("stability l2-norm of W_0 = %f \n", l2_norm_stab);
				//printf("stability max-norm of W_0 = %f \n", max_norm_stab);
				//fprintf(f1,"stability l2_norm = %e \n", l2_norm_stab);
				//fprintf(f1,"stability max_norm = %e \n", max_norm_stab);

				divergence_Wy_norm(Wy_0_error, Ny, &div_W1D_l2, &div_W1D_max);
				//printf("div_W1D_error_l2 = %f \n",div_W1D_l2);
				//printf("div_W1D_error_max = %f \n",div_W1D_max);
				divergence_norm(Wx_0_error,Wy_0_error,Nx,Ny,&l2_norm_div, &max_norm_div);
/*				printf("div l2-norm of W_0_error = %f \n", l2_norm_div);
				printf("div max-norm of W_0_error = %f \n", max_norm_div);
				fprintf(f1,"div l2_norm (of error) = %e \n", l2_norm_div);
				fprintf(f1,"div max_norm (of error) = %e \n", max_norm_div)*/;


				double l2_norm_stab, max_norm_stab;
				//stability_norm(Wx_0_error,Wy_0_error,Nx,Ny,&l2_norm_stab, &max_norm_stab);
				//printf("stability l2-norm of W_0_error = %f \n", l2_norm_stab);
				//printf("stability max-norm of W_0_error = %f \n", max_norm_stab);
				//fprintf(f1,"stability l2_norm (of error) = %e \n", l2_norm_stab);
				//fprintf(f1,"stability max_norm (of error) = %e \n", max_norm_stab);

				stability_normFull(Wx_n,Wy_n,Nx,Ny,&l2_norm_stab);
				printf("stability FULL l2-norm of W_0 = %f \n", l2_norm_stab);
				fprintf(f1,"stability FULL l2_norm W_0 = %e \n", l2_norm_stab);

				printf ( " \nCHECK \n" );
				stability_normFull(Wx_0_ex,Wy_0_ex,Nx,Ny,&l2_norm_stab);
				printf("stability FULL l2-norm of W_0_ex = %f \n", l2_norm_stab);
				fprintf(f1,"stability FULL l2_norm W_0_ex = %e \n", l2_norm_stab);

				stability_normFull(Wx_0_error,Wy_0_error,Nx,Ny,&l2_norm_stab);
				printf("stability FULL l2-norm of W_0_error = %f \n", l2_norm_stab);
				fprintf(f1,"stability FULL l2_norm (of error) = %e \n", l2_norm_stab);


				printf ( "\nChecking new stability norm2D for W0_exact \n" );
				double newstabnormm = 0.0;
				stability_norm2D ( Wx_0_ex, Wy_0_ex, Nx, Ny, tau, xCondition, yCondition, &newstabnormm);
				printf ( "new stab normm (W0_exact) = %f \n\n", newstabnormm );
				printf ( "Computing new stability norm2D (for W0_error) \n" );
				double newstabnorm = 0.0;
				stability_norm2D ( Wx_0_error, Wy_0_error, Nx, Ny, tau, xCondition, yCondition, &newstabnorm);
				printf ( "new stab norm = %f \n", newstabnorm );

                                //_getch();

#ifdef HARMONICS
				double * Vector_T0;
				Vector_T0 = ( double * ) malloc ( dim_T * sizeof ( double ) );

				//int ind_Tx0 = 5;
				//int ind_Ty0 = 9;
				//harmonics2D_T ( Vector_T0, ind_Tx0, ind_Ty0, Nx, Ny, xCondition, yCondition );

				for ( int i = 0; i < dim_T; i++ )
					Vector_T0[i] = T[i];

				myprint ( Vector_T0, dim_T, "Vector_T0.xls" );

				double * Coeffs_T0 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				coeffs_T ( Coeffs_T0, Vector_T0, Nx, Ny, xCondition, yCondition );
				myprint ( Coeffs_T0, dim_T, "Coeffs_T0.xls" );

				double * Vector_T0back = ( double * ) malloc ( dim_T * sizeof ( double ) );
				vector_T ( Vector_T0back, Coeffs_T0, Nx, Ny, xCondition, yCondition );
				myprint ( Vector_T0back, dim_T, "Vector_T0back.xls" );


				double * Coeffs_Wx0 = ( double * ) malloc ( dim_Wx * sizeof ( double ) );
				double * Coeffs_Wy0 = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				coeffs_Aminus1B ( Coeffs_Wx0, Coeffs_Wy0, Coeffs_T0, Nx, Ny, xCondition, yCondition );

				myprint ( Coeffs_Wx0, dim_Wx, "Coeff_Wx0.xls" );
				myprint ( Coeffs_Wy0, dim_Wy, "Coeff_Wy0.xls" );

				double * Vector_x, * Vector_y;
				Vector_x = ( double * ) malloc ( dim_Wx * sizeof ( double ) );
				Vector_y = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				double gamma_DD, gamma_NN;

				//FILE * f = fopen ( "Harmonics_Wx.xls","wt");
				//if ( f == NULL )
				//{
				//	printf ( "cannot open file for writing 1 \n" );
				//	return -1;
				//}
				////char filename[100], addname[10];
				//for ( int j = 0; j < Ny; j++ )
				//{
				//	for ( int i = 0; i <= Nx; i++ )
				//	{
				//		harmonics2D_Wx ( Vector_x, i, j, Nx, Ny, xCondition, yCondition );
				//		//strcpy ( filename, "Harmonics_Wx");
				//		//sprintf (addname, "_%d_%d.txt" );
				//		//strcat ( filename, addname );
				//		for ( int k = 0; k < dim_Wx; k++ )
				//			fprintf (f, "%f \t", Vector_x[k] );
				//		fprintf (f, "\n" );
				//	}
				//}
				//fclose (f);
				//f = fopen ( "Harmonics_Wy.xls","wt");
				//if ( f == NULL )
				//{
				//	printf ( "cannot open file for writing 2 \n" );
				//	return -1;
				//}
				////char filename[100], addname[10];
				//for ( int i = 0; i < Nx; i++ )
				//{
				//	for ( int j = 0; j <= Ny; j++ )
				//	{
				//		harmonics2D_Wy ( Vector_y, j, i, Nx, Ny, xCondition, yCondition );
				//		//strcpy ( filename, "Harmonics_Wx");
				//		//sprintf (addname, "_%d_%d.txt" );
				//		//strcat ( filename, addname );
				//		for ( int k = 0; k < dim_Wy; k++ )
				//			fprintf (f, "%f \t", Vector_y[k] );
				//		fprintf (f, "\n" );
				//	}
				//}
				//fclose (f);

				//double * Vector_y1, * Vector_y2;
				//Vector_y1 = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				//Vector_y2 = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				//harmonics2D_Wy ( Vector_y1, 1, 0, Nx, Ny, xCondition, yCondition );
				//harmonics2D_Wy ( Vector_y2, 3, 0, Nx, Ny, xCondition, yCondition );

				//double tempr = scal_l2 ( Vector_y1, Vector_y2, Nx, Ny + 1 );
				//printf ( "tempr = %f \n", tempr );

				//double * Vector_y1, * Vector_y2;
				//Vector_y1 = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				//Vector_y2 = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				//harmonics2D_Wy ( Vector_y1, 1, 0, Nx, Ny, xCondition, yCondition );
				//harmonics2D_Wy ( Vector_y2, 5, 0, Nx, Ny, xCondition, yCondition );

				//double tempr = scal_l2 ( Vector_y1, Vector_y2, Nx, Ny + 1 );
				//printf ( "tempr = %f \n", tempr );

				int ind_Wx = 0, ind_Ty = 0;
				int ind_Wy = 1, ind_Tx = 0;
				harmonics2D_Wx ( Vector_x, ind_Wx, ind_Ty, Nx, Ny, xCondition, yCondition );
				harmonics2D_Wy ( Vector_y, ind_Wy, ind_Tx, Nx, Ny, xCondition, yCondition );

				for ( int k = 0; k < dim_Wx; k++ )
					Vector_x[k] = Wx_n[k];
				for ( int k = 0; k < dim_Wy; k++ )
					Vector_y[k] = Wy_n[k];

				myprint ( Vector_x, dim_Wx, "Vector_x.xls" );
				myprint ( Vector_y, dim_Wy, "Vector_y.xls" );


				double * Coeff_Vec_x, * Coeff_Vec_y;
				Coeff_Vec_x = ( double * ) malloc ( dim_Wx * sizeof ( double ) );
				Coeff_Vec_y = ( double * ) malloc ( dim_Wy * sizeof ( double ) );

				coeffs_Wx ( Coeff_Vec_x, Vector_x, Nx, Ny, xCondition, yCondition );
				coeffs_Wy ( Coeff_Vec_y, Vector_y, Nx, Ny, xCondition, yCondition );

				//return 0;

				myprint ( Coeff_Vec_x, dim_Wx, "Coeff_Vec_x.xls" );
				myprint ( Coeff_Vec_y, dim_Wy, "Coeff_Vec_y.xls" );

				for ( int k = 0; k < dim_Wx; k++ )
					Coeff_Vec_x[k] = Coeffs_Wx0[k];
				for ( int k = 0; k < dim_Wy; k++ )
					Coeff_Vec_y[k] = Coeffs_Wy0[k];


				vector_Wx ( Vector_x, Coeff_Vec_x, Nx, Ny, xCondition, yCondition );
				vector_Wy ( Vector_y, Coeff_Vec_y, Nx, Ny, xCondition, yCondition );

				myprint ( Vector_x, dim_Wx, "Vector_x_after.xls" );
				myprint ( Vector_y, dim_Wy, "Vector_y_after.xls" );

				//for ( int k = 0; k < dim_Wx; k++ )
				//	Coeff_Vec_x[k] = 0.0;
				//for ( int k = 0; k < dim_Wy; k++ )
				//	Coeff_Vec_y[k] = 0.0;
				//Coeff_Vec_y[ind_Wy] = 1;

				//return 0;

				double * Btr_w0 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				myBtr_w0 ( Btr_w0, Vector_x, Vector_y, Nx, Ny );

				double * Coeffs_0 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				coeffs_T ( Coeffs_0, Btr_w0, Nx, Ny, xCondition, yCondition );
				myprint ( Coeffs_0, dim_T, "Coeffs_0.xls" );

				double wtf1 = scal_l2 ( Btr_w0, Btr_w0, Nx, Ny );
				wtf1 = sqrt ( wtf1 / dim_T );
				printf ( "wtf1 = %e \n", wtf1 );

				myprint ( Btr_w0, dim_T, "Btr_Vec.xls" );

				double * Coeffs_Btr;
				Coeffs_Btr = ( double * ) malloc ( dim_T * sizeof ( double ) );

				coeffs_Btr ( Coeffs_Btr, Coeff_Vec_x, Coeff_Vec_y, Nx, Ny, xCondition, yCondition );
//				for ( int i = 0; i < dim_T; i++ )
//					Coeffs_Btr [i] *= hx(0) * hy (0);					

				//double wtf2 = scal_l2 ( Coeffs_Btr, Coeffs_Btr, Nx, Ny );
				double * Harmonic_ij = (double *) malloc ( dim_T * sizeof ( double ) );
				double wtf2 = 0.0;
				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
						double temp = scal_l2 ( Harmonic_ij, Harmonic_ij, Nx, Ny ) / dim_T;
						wtf2 += Coeffs_Btr [ j * Nx + i ] * Coeffs_Btr [ j * Nx + i ] * temp;
					}
				}
				wtf2 = sqrt ( wtf2 );
				printf ( "wtf2 = %e \n", wtf2 );
				printf ( "wtf1/wtf2 = %f \n", wtf1/wtf2);

				double wtf3 = 0.0;
				double lambda_j, lambda_i, gamma_i, gamma_j, h_i, h_j;
				for ( int j = 0; j < Ny; j++ )
				{
					h_j = 1.0 / Ny;
					gamma_j = gamma ( j, Ny, yCondition );
					lambda_j = 6.0 * gamma_j * gamma_j / ( ( 6.0 - gamma_j * gamma_j ) * h_j * h_j ) ;
					for ( int i = 0; i < Nx; i++ )
					{
						h_i = 1.0 / Nx;
						gamma_i = gamma ( i, Nx, xCondition );
						lambda_i = 6.0 * gamma_i * gamma_i / ( ( 6.0 - gamma_i * gamma_i ) * h_i * h_i ) ;

						harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
						double temp = scal_l2 ( Harmonic_ij, Harmonic_ij, Nx, Ny ) / dim_T;
						wtf3 += ( lambda_j + lambda_i ) * Coeffs_T0[j * Nx + i] * ( lambda_j + lambda_i ) * Coeffs_T0[j * Nx + i] * temp;
					}
				}
				wtf3 = sqrt ( wtf3 );
				wtf3 *= hx(0) * hy (0);
				printf ( "wtf3 = %e \n", wtf3 );
				printf ( "wtf1/wtf3 = %f \n", wtf1/wtf3);


				myprint ( Coeffs_Btr, dim_T, "Coeffs_Btr.xls" );
				
				double * Vector_T = ( double * ) malloc ( dim_T * sizeof ( double ) );
				vector_T ( Vector_T, Coeffs_Btr, Nx, Ny, xCondition, yCondition );

				myprint ( Vector_T, dim_T, "Btr_Vec_after.xls" );

//				return 0;

//!!!				ДЛЯ (Wx, 0) дивергенция правильная, а для (0, Wy) - отличается знаком, почему?

				double * MhalfLambdaY_Btr_w0 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				myMhalfLambdayBtr_w0 ( MhalfLambdaY_Btr_w0, Vector_x, Vector_y, Nx, Ny );
				myprint ( MhalfLambdaY_Btr_w0, dim_T, "MhalfLambdayBtr_Vec.xls" );

				double * Coeffs_1 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				coeffs_T ( Coeffs_1, MhalfLambdaY_Btr_w0, Nx, Ny, xCondition, yCondition );
				myprint ( Coeffs_1, dim_T, "Coeffs_1.xls" );

				double * Coeffs_LambdayBtr = ( double * ) malloc ( dim_T * sizeof ( double ) );
				coeff_applyLambday ( Coeffs_LambdayBtr, Coeffs_Btr, Nx, Ny, xCondition, yCondition );
				//for ( int i = 0; i < dim_T; i++ )
				//	Coeffs_LambdayBtr [i] *= hx(0) * hy (0);					
				myprint ( Coeffs_LambdayBtr, dim_T, "Coeffs_MhalfLambdaBtr.xls" );

				vector_T ( Vector_T, Coeffs_LambdayBtr, Nx, Ny, xCondition, yCondition );
				//for ( int i = 0; i < dim_T; i++ )
				//	Vector_T [i] *= hy(0) * hx(0) ;					
				myprint ( Vector_T, dim_T, "MhalfLambdayBtr_Vec_after.xls" );


				double ddummy;
				double * Mminus1Btr_Aminus1_xvost_w0 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				myMminushalfLambdaLambdayBtr_w0 ( Mminus1Btr_Aminus1_xvost_w0, MhalfLambdaY_Btr_w0, &ddummy, &ddummy, Nx, Ny );
				myprint ( Mminus1Btr_Aminus1_xvost_w0, dim_T, "Mminus1LambdaLambdayBtr_Vec.xls" );

				double * Coeffs_2 = ( double * ) malloc ( dim_T * sizeof ( double ) );
				coeffs_T ( Coeffs_2, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny, xCondition, yCondition );
				myprint ( Coeffs_2, dim_T, "Coeffs_2.xls" );


				double * Coeffs_LambdaLambdayBtr = ( double * ) malloc ( dim_T * sizeof ( double ) );
				coeff_applyLambda ( Coeffs_LambdaLambdayBtr, Coeffs_LambdayBtr, Nx, Ny, xCondition, yCondition );
//				for ( int i = 0; i < dim_T; i++ )
//					Coeffs_LambdaLambdayBtr [i] /= hx(0) * hy (0);					
				myprint ( Coeffs_LambdaLambdayBtr, dim_T, "Coeffs_Mminus1LambdaLambdayBtr.xls" );

				vector_T ( Vector_T, Coeffs_LambdaLambdayBtr, Nx, Ny, xCondition, yCondition );
//				for ( int i = 0; i < dim_T; i++ )
//					Vector_T [i] *= hy(0) * hx(0) ;					
				myprint ( Vector_T, dim_T, "Mminus1LambdaLambdayBtr_Vec_after.xls" );

				double l2_normm;
				l2_normm = 0.0;
				stability_normFull( Vector_x, Vector_y, Nx, Ny, &l2_normm);
				printf ( "stab_norm = %e \n", l2_normm );

				//l2_normm = sqrt ( scal_l2 ( Coeffs_LambdaLambdayBtr, Coeffs_LambdayBtr, Nx, Ny ) / ( dim_T ) );
				double norm_sq;
				l2_normm = 0.0;
				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
						norm_sq = scal_l2 ( Harmonic_ij, Harmonic_ij, Nx, Ny ) / dim_T;
						l2_normm += Coeffs_LambdaLambdayBtr [ j * Nx + i ] * Coeffs_LambdayBtr [ j * Nx + i ] * norm_sq;
					}
				}
				l2_normm = sqrt ( l2_normm );
				printf ( "stab_norm_coeffs = %e \n", l2_normm );

				l2_normm = 0.0;
				//double gamma_i, gamma_j;
				//double lambda_i, lambda_j;
				h_i = 1.0 / Nx;
				h_j = 1.0 / Ny;
				for ( int j = 0; j < Ny; j++ )
				{
					gamma_j = gamma ( j, Ny, yCondition );
					lambda_j = 6.0 * gamma_j * gamma_j / ( ( 6.0 - gamma_j * gamma_j ) * h_j * h_j ) ;
					for ( int i = 0; i < Nx; i++ )
					{
						gamma_i = gamma ( i, Nx, xCondition );
						lambda_i = 6.0 * gamma_i * gamma_i / ( ( 6.0 - gamma_i * gamma_i ) * h_i * h_i ) ;

						harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
						norm_sq = scal_l2 ( Harmonic_ij, Harmonic_ij, Nx, Ny ) / dim_T;

						l2_normm += lambda_j * lambda_j * ( lambda_j + lambda_i ) * ( lambda_j + lambda_i ) * ( lambda_j + lambda_i ) * Coeffs_T0 [ j * Nx + i ] * Coeffs_T0 [ j * Nx + i ] * norm_sq;

					}
				}
				l2_normm = sqrt ( l2_normm );
				//l2_normm *= h_i; 
				printf ( "stab_norm_coeffsfromT0 = %e \n", l2_normm );

				//for ( int j = 0; j < Ny; j++ )
				//{
				//	printf ( "j = %d \n", j );
				//	gamma_j = gamma ( j, Ny, yCondition );
				//	lambda_j = 6.0 * gamma_j * gamma_j / ( ( 6.0 - gamma_j * gamma_j ) * h_j * h_j ) ;
				//	printf ( "gamma_j = %f \n", gamma_j );
				//	printf ( "lambda_j = %f \n", lambda_j );
				//	for ( int i = 0; i < Nx; i++ )
				//	{
				//		gamma_i = gamma ( i, Nx, xCondition );
				//		lambda_i = 6.0 * gamma_i * gamma_i / ( ( 6.0 - gamma_i * gamma_i ) * h_i * h_i ) ;

				//		harmonics2D_T ( Harmonic_ij, i, j, Nx, Ny, xCondition, yCondition );
				//		norm_sq = scal_l2 ( Harmonic_ij, Harmonic_ij, Nx, Ny ) / dim_T;

				//		printf ( "i = %d \n", i );
				//		printf ( "gamma_i = %f \n", gamma_i );
				//		printf ( "lambda_i = %f \n", lambda_i );
				//		printf ( "norm_sq = %f \n", norm_sq );
				//	}
				//	_getch();
				//}
				//return 0;
				//_getch();
				//
				//harmonics2D_T ( Harmonic_ij, Nx/2, Ny/2, Nx, Ny, xCondition, yCondition );
				//double tempp = scal_l2 ( Harmonic_ij, Harmonic_ij, Nx, Ny );
				//printf ( "tempp = %f \n", tempp );
				//l2_normm *= tempp;
				//printf ( "stab_norm_coeffs after = %e \n", l2_normm );
				
//				stability_normFull_coeffs( Vector_x, Vector_y, Nx, Ny, &l2_normm, &max_normm);

#ifdef MIND_GAME
				printf ( "Mind game begins \n" );

				FILE * f_h;
				f_h = fopen ( "h_dependance_Wx_scheme_300_32.xls", "wt" );
				if ( f_h == NULL )
				{
					printf ( "Cannot open file f_h for Wx \n" );
					return -1;
				}
				double res;
				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = 0; i <= Nx; i++ )
					{
						res = 0.0;
						//if ( fabs ( Coeff_Vec_x [ j * ( Nx + 1 ) + i ] ) > 1.0e-12 )
						{
							int ind_Wx = i, ind_Ty = j;
							int ind_Wy = 0, ind_Tx = 0;
							harmonics2D_Wx ( Vector_x, ind_Wx, ind_Ty, Nx, Ny, xCondition, yCondition );
							harmonics2D_Wy ( Vector_y, ind_Wy, ind_Tx, Nx, Ny, xCondition, yCondition );

							stability_normFull( Vector_x, Vector_y, Nx, Ny, &res);


							//coeffs_Wx ( Coeff_Vec_x, Vector_x, Nx, Ny, xCondition, yCondition );
							//coeffs_Wy ( Coeff_Vec_y, Vector_y, Nx, Ny, xCondition, yCondition );

							//coeffs_Btr ( Coeffs_Btr, Coeff_Vec_x, Coeff_Vec_y, Nx, Ny, xCondition, yCondition );
							//for ( int i = 0; i < dim_T; i++ )
							//	Coeffs_Btr [i] *= hx(0) * hy (0);					
			
							//coeff_applyLambday ( Coeffs_LambdayBtr, Coeffs_Btr, Nx, Ny, xCondition, yCondition );

							//coeff_applyLambda ( Coeffs_LambdaLambdayBtr, Coeffs_LambdayBtr, Nx, Ny, xCondition, yCondition );
							//for ( int i = 0; i < dim_T; i++ )
							//	Coeffs_LambdaLambdayBtr [i] /= hx(0) * hy (0);					

						}
						//else
						//	res = 0.0;
						fprintf ( f_h, "%f \n", res );

					}
				}
				fclose ( f_h );

				f_h = fopen ( "h_dependance_Wy_scheme_300_32.xls", "wt" );
				if ( f_h == NULL )
				{
					printf ( "Cannot open file f_h for Wy \n" );
					return -1;
				}
				for ( int i = 0; i < Nx; i++ )
				{
					for ( int j = 0; j <= Ny; j++ )
					{
						res = 0.0;
						//if ( fabs ( Coeff_Vec_y [ i * ( Ny + 1 ) + j ] ) > 1.0e-12 )
						{
							int ind_Wx = 0, ind_Ty = 0;
							int ind_Wy = j, ind_Tx = i;
							harmonics2D_Wx ( Vector_x, ind_Wx, ind_Ty, Nx, Ny, xCondition, yCondition );
							harmonics2D_Wy ( Vector_y, ind_Wy, ind_Tx, Nx, Ny, xCondition, yCondition );

							stability_normFull( Vector_x, Vector_y, Nx, Ny, &res);
						}
						//else
						//	res = 0.0;
						fprintf ( f_h, "%f \n", res );

					}
				}
				fclose ( f_h );
				//coeffs_Wx ( Coeff_Vec_x, Wx_n, Nx, Ny, xCondition, yCondition );
				//coeffs_Wy ( Coeff_Vec_y, Wy_n, Nx, Ny, xCondition, yCondition );

				//coeffs_Btr ( Coeffs_Btr, Coeff_Vec_x, Coeff_Vec_y, Nx, Ny, xCondition, yCondition );
				//for ( int i = 0; i < dim_T; i++ )
				//	Coeffs_Btr [i] *= hx(0) * hy (0);					

				//coeff_applyLambday ( Coeffs_LambdayBtr, Coeffs_Btr, Nx, Ny, xCondition, yCondition );

				printf ( "Mind games ends \n" );
#endif

				free ( Btr_w0 );
				free ( Vector_x );
				free ( Vector_y );
				free ( Vector_T );
#endif
                                //_getch();
//				return 0;


#ifdef MINIMIZATION
				printf ( "\nStarting minimization procedure for W_0 \n" );
				_getch();
				double * Wx_0_new, * Wy_0_new;
//#ifdef CG_TEST
//				double coeff1 = 1.0;
//				double coeff2 = 0.0;
//#else
				double h = 1.0 / Nx;
//				double coeff1 = 1.0;
//				double coeff2 = 0.0;
				double coeff1 = 1.0 / ( h * h );
				double coeff2 = 1.0 * h * h * h * h;
//#endif
				double converge_tol = 1.0e-2 * h * h * h *h;
				double stagn_tol = 1.0e-4;
				double grad_tol = 1.0e-8;
				int stag_limit = 200;
				int iter_limit = 50000;

				//int kcount = 0;
				//double tempk = 0;
				//while ( true )
				//{
				//	//myByMminus1Btr_w0 ( Wy_n, Wx_n, Wy_n, Nx, Ny );
				//	//myAyminus1ByMminus1Btr_w0 ( Wy_n, Wx_n, Wy_n, Nx, Ny );
				//	//myMhalfLambdayBtr_w0 ( T, Wx_n, Wy_n, Nx, Ny );
				//	stability_normFull( Wx_n, Wy_n, Nx, Ny, &tempk );
				//	printf ( "kcount = %d \n", kcount );
				//	kcount++;
				//	_getch();
				//}

#ifdef CG_ZERO_INIT
				double * Wx_start = ( double * ) malloc ( dim_Wx * sizeof ( double ) );
				double * Wy_start = ( double * ) malloc ( dim_Wy * sizeof ( double ) );
				for ( int i = 0; i < dim_Wx; i++ )
					Wx_start [ i ] = 0.0;
				for ( int i = 0; i < dim_Wy; i++ )
					Wy_start [ i ] = 0.0;
				myMinimization ( Wx_start, Wy_start, T, &Wx_0_new, &Wy_0_new, converge_tol, stagn_tol, grad_tol, coeff1, coeff2, iter_limit, stag_limit, h, Nx, Ny, Nz );

				free ( Wx_start );
				free ( Wy_start );
#else
				myMinimization ( Wx_n, Wy_n, T, &Wx_0_new, &Wy_0_new, converge_tol, stagn_tol, grad_tol, coeff1, coeff2, iter_limit, stag_limit, h, Nx, Ny, Nz );
#endif

				printf ( "Minimization core finished \n" );
				_getch();

				double temp1 = (1.0 /  (h * h * h * h) ) * myFunctional ( Wx_n, Wy_n, T, Nx, Ny, 1.0, 0.0 );
				//printf ( "check 1 \n" );
				//_getch();
				//return 0;
				double temp2 = (1.0 /  (h * h * h * h) ) * myFunctional ( Wx_0_new, Wy_0_new, T, Nx, Ny, 1.0, 0.0 );
				printf ( "Nevyazka^2 (W_0) = %e Nevyazka^2 (W_0_optimal) = %e \n", temp1, temp2 );
				temp1 = 0.0;
				normAQminusF ( Wx_n, Wy_n, T, Nx, Ny, &temp1 );
				temp1 *= temp1 /  (h * h * h * h);
				temp2 = 0.0;
				normAQminusF ( Wx_0_new, Wy_0_new, T, Nx, Ny, &temp2 );
				temp2 *= temp2 /  (h * h * h * h);
				printf ( "Nevyazka^2 (W_0) = %e Nevyazka^2 (W_0_optimal) = %e \n", temp1, temp2 );
				//return 0;

				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_error[i*(Nx+1) + k] = Wx_0_new[i*(Nx+1) + k] - Wx_n[i*(Nx+1) + k];
					}
				}
				for (int i=0; i<Nz*Nx; i++)
				{
					for (int k=0; k<Ny+1; k++)
					{
						i1 = i/Nz;
						j1 = k;
						k1 = i - (i/Nz)*Nz;

						Wy_0_error[i*(Ny+1) + k] = Wy_0_new[i*(Ny+1) + k] - Wy_n[i*(Ny+1) + k];
					}
				}

				double temprrr = 0.0;
				temprrr += norm_l2 (Wx_0_error, Nx + 1, Ny ) * norm_l2 (Wx_0_error, Nx + 1, Ny );
				temprrr += norm_l2 (Wy_0_error, Ny + 1, Nx ) * norm_l2 (Wy_0_error, Ny + 1, Nx );
				temprrr = sqrt (temprrr);

				printf ( "difference between W_0 and 'optimal' W_0 in L2 = %e \n", temprrr );
				_getch();

				for (int i = 0; i < Ny*Nz; i++)
				{
					for (int k = 0; k < Nx+1; k++)
					{
						i1 = k;
						j1 = i - (i/Ny)*Ny;
						k1 = i/Ny;

						Wx_0_error[i*(Nx+1) + k] = Wx_0_new[i*(Nx+1) + k] - (-heat_conductivity_func(xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1))*exact_gradientX(numSolution,xpoints[k],ypoints[j1] + 0.5*hy(j1), zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}
				for (int i=0; i<Nz*Nx; i++)
				{
					for (int k=0; k<Ny+1; k++)
					{
						i1 = i/Nz;
						j1 = k;
						k1 = i - (i/Nz)*Nz;

						Wy_0_error[i*(Ny+1) + k] = Wy_0_new[i*(Ny+1) + k] - (-heat_conductivity_func(xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1))*exact_gradientY(numSolution,xpoints[i1]+0.5*hx(i1),ypoints[k], zpoints[k1] + 0.5*hz(k1), 0) );
					}
				}
				temprrr = 0.0;
				temprrr += norm_l2 (Wx_0_error, Nx + 1, Ny ) * norm_l2 (Wx_0_error, Nx + 1, Ny );
				temprrr += norm_l2 (Wy_0_error, Ny + 1, Nx ) * norm_l2 (Wy_0_error, Ny + 1, Nx );
				temprrr = sqrt (temprrr);
				printf ( "difference between W_0_exact and 'optimal' W_0 in L2 = %e \n", temprrr );
				_getch();

				stability_normFull(Wx_0_new,Wy_0_new,Nx,Ny,&l2_norm_stab);
				printf("stability FULL l2-norm of W_0_new after minim = %f \n", l2_norm_stab);
				stability_normFull(Wx_0_error,Wy_0_error,Nx,Ny,&l2_norm_stab);
				printf("stability FULL l2-norm of W_0_error after minim = %f \n", l2_norm_stab);
				printf ( "\nEnding the entire minimization procedure for W_0 \n" );
#ifdef MINIMIZATION_W0_FOR_SCHEME
				for ( int i = 0; i < dim_Wx; i++ )
					Wx_n [i] = Wx_0_new [i];
				for ( int i = 0; i < dim_Wy; i++ )
					Wy_n [i] = Wy_0_new [i];
#endif
				{;}
				//stability_normFull(Wx_n,Wy_n,Nx,Ny,&l2_norm_stab);
				//printf("stability FULL l2-norm of W_0 = %f \n", l2_norm_stab);
				//fprintf(f1,"stability FULL l2_norm = %e \n", l2_norm_stab);

				//double l2_norm_stab2, max_norm_stab2;
				//printf("FINAL STABILITY NORM ||Lamb2 * B^tr w0 ||_(Lamb2 + Lamb1) \n");
				//stability_norm_sum(Wx_0_error,Wy_0_error,Nx,Ny,&l2_norm_stab2, &max_norm_stab2);
				//printf("stability l2-norm sum of W_0_error = %f \n", l2_norm_stab2);
				//printf("stability max-norm sum of W_0_error = %f \n", max_norm_stab2);
				//fprintf(f1,"stability l2_norm sum (of error) = %e \n", l2_norm_stab2);
				//fprintf(f1,"stability max_norm sum (of error) = %e \n", max_norm_stab2);

#ifdef SPECIAL_MODE
				for (int i = 0; i < dim_Wx; i++)
					Wx_n[i] = Wx_0_error[i];
				for (int i = 0; i < dim_Wy; i++)
					Wy_n[i] = Wy_0_error[i];
#endif

#endif //for #ifdef MINIMIZATION
				//return 0;
				//				_getch();

				//double eps0_max_T = 0;
				//double eps0_l2_T = 0;
				//double eps0_relative_max_T = 0;
				//double eps0_relative_l2_T = 0;

				//printf("FUCK YOU \n");
				//if (f1 == NULL)
				//  printf("FUCK YOU twice \n");

				fprintf(f1,"T = 0: \n");
				fprintf(f1,"eps0_max_T = %e \n", eps0_max_T);
				fprintf(f1,"eps0_l2_T = %e \n", eps0_l2_T);
				fprintf(f1,"eps0_max_w = %e \n", eps0_max_w);
				fprintf(f1,"eps0_l2_w = %e \n", eps0_l2_w);
				fprintf(f1,"releps0_max_T = %e \n", eps0_relative_max_T);
				fprintf(f1,"releps0_l2_T = %e \n", eps0_relative_l2_T);
				fprintf(f1,"releps0_max_w = %e \n", eps0_relative_max_w);
				fprintf(f1,"releps0_l2_w = %e \n", eps0_relative_l2_w);

				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				printf("Starting main loop...\n");
				double* Wx_nplus1 = (double*)malloc(dim_Wx * sizeof(double));
				double* Wy_nplus1 = (double*)malloc(dim_Wy * sizeof(double));
				double* Wx_nplus05 = (double*)malloc(dim_Wx * sizeof(double));
				double* Wy_nplus05 = (double*)malloc(dim_Wy * sizeof(double));
				double * Temp_Wx = (double*)malloc(dim_Wx * sizeof(double));
				double * Temp_Wy = (double*)malloc(dim_Wy * sizeof(double));
				double * T_plus1 = (double*)malloc(dim_T * sizeof(double));
				double * F = (double*)malloc(dim_T * sizeof(double));
				double * righthandX = (double*)malloc(dim_Wx * sizeof(double));
				double * righthandY = (double*)malloc(dim_Wy * sizeof(double));

				double eps_max_T = 0;
				double eps_l2_T = 0;
				double eps_max_w = 0;			
				double eps_l2_w = 0;	
				double eps_relative_max_w = 0.0;
				double eps_relative_max_T = 0.0;
				double eps_relative_l2_w = 0.0;
				double eps_relative_l2_T = 0.0;

#ifdef DEBUGE
				N = 1;
#endif
//				N = 1;
				for ( int t = 0 ; t < N ; t++ )
				{
					Splitting2D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//Uzawa2D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//Local1D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//PredictorCorrector_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//PredictorPaper2D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//TwoCyclic_Local1D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);

					//if (t == 0)
					//	Uzawa2D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//else
					//	Splitting2D_Scheme_mainstep(T, T_plus1, F, Temp_Wx, Temp_Wy, righthandX, Vx, Wx_n, Wx_nplus05, Wx_nplus1, righthandY, Vy, Wy_n, Wy_nplus05, Wy_nplus1, 0.5, t, tau, numSolution, xCondition, yCondition, Nx, Ny, Nz, print_step, &eps_max_w, &eps_max_T, &eps_l2_w, &eps_l2_T, &eps_relative_max_w, &eps_relative_max_T, &eps_relative_l2_w, &eps_relative_l2_T);
					//_getch();
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

Exit:
				fprintf(f1_excel,"Final: \n");
#ifdef EXACT_W0_INIT
				fprintf(f1_excel,"tags: EXACT_W0_INIT \n");
#endif
				fprintf(f1_excel,"Nx = \t %d \t Ny = %d \t Nz = %d \ntau = \t %8.8f \t N = \t %d \n\n", Nx, Ny, Nz, tau, N);
				//fprintf(f1_excel,"eps_max_T = \t%.6g\t", eps_max_T);
				//fprintf(f1_excel,"eps_l2_T = \t%.6g\t\n", eps_l2_T);
				//fprintf(f1_excel,"eps_max_w = \t%.6g\t\n", eps_max_w);
				//fprintf(f1_excel,"eps_l2_w = \t%.6g\t\n", eps_l2_w);

				fprintf(f1_excel,"eps_max_T\teps_l2_T\teps_max_w\teps_l2_w\t\n");
				fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n",eps_max_T,eps_l2_T,eps_max_w,eps_l2_w);

				fprintf(f1_excel,"eps_max_T\teps_l2_T\teps_max_w\teps_l2_w\t\n");
				fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n",eps_relative_max_T,eps_relative_l2_T,eps_relative_max_w,eps_relative_l2_w);

				fprintf(f1,"Nx = %d, Ny = %d, Nz = %d \ntau = %f N = %d \n", Nx, Ny, Nz, tau, N);
				fprintf(f1,"eps_max_T = %e \n", eps_max_T);
				fprintf(f1,"eps_l2_T = %e \n", eps_l2_T);
				fprintf(f1,"eps_max_w = %e \n", eps_max_w);
				fprintf(f1,"eps_l2_w = %e \n", eps_l2_w);
				fprintf(f1,"eps_relative_max_T = %e \n", eps_relative_max_T);
				fprintf(f1,"eps_relative_l2_T = %e \n", eps_relative_l2_T);
				fprintf(f1,"eps_relative_max_w = %e \n", eps_relative_max_w);
				fprintf(f1,"eps_relative_l2_w = %e \n", eps_relative_l2_w);

				////for calculating accruacy for the t = N - 1
				//eps_max_w = 0;
				//eps_max_T = 0;
				//eps_l2_w = 0;
				//eps_l2_T = 0;
				//eps_relative_max_w = 0;
				//eps_relative_max_T = 0;
				//eps_relative_l2_w = 0;
				//eps_relative_l2_T = 0;

				eps0_max_w = 0;
				eps0_max_T = 0;
				eps0_l2_w = 0;
				eps0_l2_T = 0;
				eps0_relative_max_w = 0;
				eps0_relative_max_T = 0;
				eps0_relative_l2_w = 0;
				eps0_relative_l2_T = 0;
				
				Accuracy_calculate(T, Wx_n, Wy_n, numSolution, N-1, Nx, Ny, Nz, print_step, &eps0_max_w, &eps0_max_T, &eps0_l2_w, &eps0_l2_T, &eps0_relative_max_w, &eps0_relative_max_T, &eps0_relative_l2_w, &eps0_relative_l2_T);

				//printf("WTTTTTTTTTTTTFFFFFFF \n");
				//printf("eps0_max_w = %f \n", eps0_max_w);

				fprintf(f1_excel,"T = 1: \n");
				fprintf(f1_excel,"eps_max_T\teps_l2_T\teps_max_w\teps_l2_w\t\n");
				fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n",eps0_max_T,eps0_l2_T,eps0_max_w,eps0_l2_w);

				fprintf(f1_excel,"eps_max_T\teps_l2_T\teps_max_w\teps_l2_w\t\n");
				fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n",eps0_relative_max_T,eps0_relative_l2_T,eps0_relative_max_w,eps0_relative_l2_w);

				fprintf(f1,"T = 1: \n");
				fprintf(f1,"eps_max_T = %e \n", eps0_max_T);
				fprintf(f1,"eps_l2_T = %e \n", eps0_l2_T);
				fprintf(f1,"eps_max_w = %e \n", eps0_max_w);
				fprintf(f1,"eps_l2_w = %e \n", eps0_l2_w);
				fprintf(f1,"eps_max_T = %e \n", eps0_relative_max_T);
				fprintf(f1,"eps_l2_T = %e \n", eps0_relative_l2_T);
				fprintf(f1,"eps_max_w = %e \n", eps0_relative_max_w);
				fprintf(f1,"eps_l2_w = %e \n", eps0_relative_l2_w);


				//fprintf(f1_excel,"special_eps_max_T\tspecial_eps_l2_T\tspecial_eps_max_w\tspecial_eps_l2_w\t\n");
				//fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n\n",special_eps_max_T,special_eps_l2_T,special_eps_max_w,special_eps_l2_w);

				fprintf(f1_excel,"T = 1: \n");

				//fprintf(f1_excel,"my_releps_max_T\tmy_releps_l2_T\tmy_releps_max_w\tmy_releps_l2_w\t\n");
				//fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n",my_releps_max_T,my_releps_l2_T,my_releps_max_w,my_releps_l2_w);

				//fprintf(f1_excel,"my_eps_max_T\tmy_eps_l2_T\tmy_eps_max_w\tmy_eps_l2_w\t\n");
				//fprintf(f1_excel,"%.6e\t%.6e\t%.6e\t%.6e\t\n\n",my_eps_max_T,my_eps_l2_T,my_eps_max_w,my_eps_l2_w);

				fprintf(f1_excel,"numSolution = \t %d \t mx = \t %f \t my = \t %f\nexternal = \t %d \n", numSolution,mx,my,external);
				fprintf(f1_excel,"Tx0 = \t %8.8f \t Tx1 = \t  %8.8f\nTy0 = \t %8.8f \t Ty_1 = \t %8.8f \nwx_0 = \t %8.8f \t wx_1 = \t  %8.8f\nwy_0 = \t %8.8f \t wy_1 = \t %8.8f \n", Tx_0, Tx_1, Ty_0, Ty_1, wx_0, wx_1, wy_0, wy_1);

				//Выдача теплового потока W , температуры T и погрешностей eps_max и eps_l2
				fprintf(f1,"Final: \n");
#ifdef GETCH
				fprintf(f1,"tags: GETCH \n");
#endif
#ifdef ACCURACY
				fprintf(f1,"tags: ACCURACY \n");
#endif
#ifdef NEW_BOUNDARY_CONDITION
				fprintf(f1,"tags: NEW_BOUNDARY_CONDITION \n");
#endif
#ifdef EXACT_BOUNDARY_CONDITION
				fprintf(f1,"tags: EXACT_BOUNDARY_CONDITION \n");
#endif
#ifdef ZEROING
				fprintf(f1,"tags: ZEROING \n");
#endif
#ifdef EXACT_W0_INIT
				fprintf(f1,"tags: EXACT_W0_INIT \n");
#endif
#ifdef SPECIAL_MODE
				fprintf(f1,"tags: SPECIAL_MODE \n");
#endif
#ifdef EIGENVECTOR_WY0
				fprintf(f1,"tags: EIGENVECTOR_WY0 \n");
#endif
#ifdef SADDLE_POINT
				fprintf(f1,"tags: SADDLE_POINT \n");
#endif
#ifdef SPECIAL_MODE_T
				fprintf(f1,"tags: SPECIAL_MODE_T \n");
#endif
#ifdef LUMPING_INIT
				fprintf(f1,"tags: LUMPING_INIT \n");
#endif
#ifdef STEKLOV
				fprintf(f1,"tags: STEKLOV \n");
#endif

				fprintf(f1,"numSolution = %d, mx = %f, my = %f mz = %f\nexternal = %d \n", numSolution,mx,my,mz,external);
				fprintf(f1,"velocityX_num = %d, velocityY_num = %d velocityZ_num = %d\n", velocityX_num,velocityY_num, velocityZ_num);
				fprintf(f1,"Tx0 = %f, Tx1 = %f\nTy0 = %f, Ty_1 = %f \nwx_0 = %f, wx_1 = %f\nwy_0 = %f, wy_1 = %f \n", Tx_0, Tx_1, Ty_0, Ty_1, wx_0, wx_1, wy_0, wy_1);
				fprintf(f1,"vx_const = %f, vy_const = %f, vz_const = %f \n", vx_const, vy_const, vz_const);
				//fprintf(f1,"special_eps_l2_T = %e \n", special_eps_l2_T);
				//fprintf(f1,"special_eps_l2_w = %e \n", special_eps_l2_w);
				//fprintf(f1,"special_eps_max_T = %e \n", special_eps_max_T);
				//fprintf(f1,"special_eps_max_w = %e \n", special_eps_max_w);

				//unsigned int EndTime = GetTickCount();
				//unsigned int Time_ms = EndTime - StartTime;
				//fprintf(f1,"Time total (milliseconds) = %d \n",Time_ms);
				//fprintf(f1,"period of Time = [0,%f] \n",Time);

				////////////
				//Освобождение памяти из-под всех использованных динамических массивов
				free(Wx_n);
				free(Wy_n);
				//free(Wz_n);
				free(Vx);
				free(Vy);
				//free(Vz);
				free(Wx_nplus1);
				free(Wy_nplus1);
				free(Wx_nplus05);
				free(Wy_nplus05);
				//free(Wx_temp);
				//free(delta_Wx_nplus1);
				//free(delta_Wy_nplus1);
				//free(delta_Wx_nplus05);
				//free(delta_Wy_nplus05);
#ifdef NEW_BOUNDARY_CONDITION
				free(Gx_temp);
				free(Gy_temp);
#endif

				free(T);
				free(T_plus1);
				free(betax);
				free(alfax);
				free(betay);
				free(alfay);
				fclose(f1);
				fclose(f1_excel);
				//////////
				return 0;
}



double v_init(double x)
{
	return x*(1-x);
}


int maximum(int N, int M, int K)
{
	if ( N >=M )
	{
		if ( N > K )
			return N;
		else
			return K;
	}

	else
	{
		if ( M > K )
			return M;
		else 
			return K;
	}
}

int stepen2(int n)
{
	if (n==0)
		return 1;
	else
		return 2*stepen2(n-1);
}




void Output_filename(char* filename, int numSolution)
{
	char addname_alpha[10];
	switch(numSolution)
	{
	case 101:
		strcpy(filename,"output/resultSplit_2D_DirZ_");
		break;

	case 102:
		strcpy(filename,"output/resultSplit_2D_DirY_");
		break;
	case 103:
		strcpy(filename,"output/resultSplit_2D_DirX_");
		break;
	case 104:
		if (mx != 0)
			strcpy(filename,"output/resultSplit_2D_NeuX_");
		else
			if (my != 0)
				strcpy(filename,"output/resultSplit_2D_NeuY_");
			else
				strcpy(filename,"output/resultSplit_2D_NeuZ_");
		break;
	case 105:
		if (Dir_axis == 0)
			strcpy(filename,"output/resultSplit_2D_DirX_");
		else
			if (Dir_axis == 1)
				strcpy(filename,"output/resultSplit_2D_DirY_");
			else
				strcpy(filename,"output/resultSplit_2D_DirZ_");
		break;
	case 106:
		strcpy(filename,"output/resultSplit_2D_MixedZ_");
		break;
	case 107:
		strcpy(filename,"output/resultSplit_2D_NeuY_");
		break;
	case 108:
		strcpy(filename,"output/resultSplit_2D_ArbEx1_");
		break;
	case 109:
		strcpy(filename,"output/resultSplit_2D_ArbEx2_");
		break;
	case 110:
		strcpy(filename,"output/resultSplit_2D_ArbEx3_");
		sprintf(addname_alpha,"%3.1f",alpha);
		strcat(filename,addname_alpha);
		break;
	case 111:
		strcpy(filename,"output/resultSplit_2D_ArbEx4noTime_");
		sprintf(addname_alpha,"%3.1f",alpha);
		strcat(filename,addname_alpha);
		break;
	case 112:
		strcpy(filename,"output/resultSplit_2D_ArbEx5noTimeSymmetric_");
		sprintf(addname_alpha,"%3.1f",alpha);
		strcat(filename,addname_alpha);
		break;
	case 113:
		strcpy(filename,"output/resultSplit_2D_ArbEx5noTimeSymmetricX_");
		sprintf(addname_alpha,"%3.1f",alpha);
		strcat(filename,addname_alpha);
		break;
	case 114:
		strcpy(filename,"output/resultSplit_2D_ArbEx5noTimeSymmetricY_");
		sprintf(addname_alpha,"%3.1f",alpha);
		strcat(filename,addname_alpha);
		break;
	case 115:
		strcpy(filename,"output/resultSplit_2D_ArbEx5onlyTime_");
		sprintf(addname_alpha,"%3.1f",alpha);
		strcat(filename,addname_alpha);
		break;
	default:
		strcpy(filename,"output/resultSplit_2D_wtf_");
		break;
	}
}


void meshGeneration(int version, char **argv)
{
	if (version == 1)
	{
		FILE *input;
		input = fopen("data.txt", "r");
		if (input == NULL)
		{
			printf("FUCK YOU! No file data.txt is found and no parameters are given \n");
		}
		else
		{
			fscanf(input, "%lf", &tau);
			fscanf(input, "%d", &Nx);
			fscanf(input, "%d", &Ny);
			fscanf(input, "%d", &Nz);
			fscanf(input, "%d", &external);
			fscanf(input, "%d", &numSolution);
		}
	}
	else
	{
		printf("Enter parameters tau, Nx, Ny, Nz, external and numSolution \n");
		tau = atof(argv[1]);
		Nx = atoi(argv[2]);
		Ny = atoi(argv[3]);
		Nz = atoi(argv[4]);
		external = atoi(argv[5]);
		numSolution = atoi(argv[6]);
	}
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

void FilePrint2D(char *filename, double *Array, int dim1_array, int dim2_array)
{
	FILE *file = fopen(filename,"wt");
	for (int j = 0;j < dim2_array; j++)
	{
		for (int i = 0;i < dim1_array; i++)
			fprintf(file,"%f \t",Array[j*dim1_array + i]);
		fprintf(file,"\n");
	}
	fclose(file);
}

void stability_norm(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
#ifdef DEBUGE
	char* ext = ".xls";
	char nx_string[5];
	sprintf(nx_string,"_%d",Nx);
	char filename[30];
	strcpy(filename,"Wy");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * Wy_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(Wy_file,"%f \n",wy[ind]);
	fclose(Wy_file);
	FilePrint2D("Wy_2D.xls",wy,Ny + 1, Nx);

	strcpy(filename,"Wx");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * Wx_file = fopen("Wx.xls","wt");
	FILE * Wx_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWx; ind++)
		fprintf(Wx_file,"%f \n",wx[ind]);
	fclose(Wx_file);
	FilePrint2D("Wx_2D.xls",wx,Nx + 1, Ny);
#endif
	//printf("Calculating ByB_tr begins ... \n");

	double *ByB_tr = (double*)malloc(dimWy * sizeof(double));
	//computing ByB_tr W
	double m_iminus1, m_i;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			ByB_tr[i*(Ny + 1) + j] = 0.0;
			if (j > 0)
				m_iminus1 = heat_capacity(i,j-1,0)*density(i,j-1,0)*hy(j-1)*hx(i);
			if (j < Ny)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (j > 0 && j < Ny)
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1])  - (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)])  - hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else if (j == 0)
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*(- (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(- hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)]) );
			}
		}
	}

	//printf("Calculating ByB_tr finishes ... \n");

#ifdef DEBUGE
	strcpy(filename,"ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * ByB_tr_file = fopen("ByB_tr.xls","wt");
	FILE * ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(ByB_tr_file,"%f \n",ByB_tr[ind]);
	fclose(ByB_tr_file);
	FilePrint2D("ByB_tr_2D.xls",ByB_tr,Ny + 1, Nx);
#endif
	double temp1, temp2;
	temp1 = norm_l2(ByB_tr,Ny + 1, Nx);
	temp2 = norm_max(ByB_tr,Ny + 1, Nx);
	printf("norm_l2 for By_B_tr = %f \n", temp1);
	printf("norm_max for By_B_tr = %f \n", temp2);


	//return 0.0;

	double * Aminus1_By_B_tr = (double*)malloc(dimWy * sizeof(double));
	//computing Ay(-1) * BxBy_tr W 
	//...

	//printf("Calculating Aminus1_By_B_tr begins ... \n");

	double *alfay = (double*)malloc((Ny + 1)* sizeof(double));
	double *betay = (double*)malloc((Ny + 1)* sizeof(double));
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	double temper_0, temper_j, temper_n;
	for ( int i = 0 ; i < Nx ; i++ )
	{

		//к-ты матрицы А
		a_0_old = hy(0)*hx(i)/(heat_conductivity(i,0,0)*6.0);
		b_0_old = 2*a_0_old;
		temper_0 = ByB_tr[i*(Ny + 1) + 0];

		betay[0] = 0 ;
		alfay[0] = 0 ;
		alfay[1] = -0.5 ;
		betay[1] =  temper_0 / b_0_old; 

		for ( int j = 1 ; j < Ny ; j++ )
		{
			//к-ты матрицы А
			a_i_old = hy(j)*hx(i)/(heat_conductivity(i,j,0)*6.0);
			a_iminus1_old = hy(j-1)*hx(i)/(heat_conductivity(i,j-1,0)*6.0);
			b_i_old = 2*a_iminus1_old + 2*a_i_old;

			temper_j = ByB_tr[i*(Ny + 1) + j];

			alfay[j + 1] = (-1)*a_i_old/(alfay[j]*a_iminus1_old + b_i_old);
			betay[j + 1] = ( ( temper_j ) - betay[j]*a_iminus1_old ) / ( alfay[j]*a_iminus1_old + b_i_old );

		}
		a_nminus1_old = hy(Ny-1)*hx(i)/(heat_conductivity(i,Ny-1,0)*6.0);
		b_n_old = 2*a_nminus1_old;

		temper_n = ByB_tr[i*(Ny + 1) + Ny];

		Aminus1_By_B_tr[i*(Ny+1) + Ny] = ( temper_n  -  betay[Ny]*a_nminus1_old ) / (alfay[Ny]*a_nminus1_old + b_n_old);

		for ( int j = 1 ; j < Ny + 1 ; j++ )
		{
			Aminus1_By_B_tr[i*(Ny+1) + Ny - j ] = Aminus1_By_B_tr[i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
		}
	}
	free(alfay);
	free(betay);

	//printf("Calculating Aminus1_By_B_tr finishes ... \n");

	//return 0.0;


#ifdef DEBUGE
	strcpy(filename,"Aminus1_ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * Aminus1_ByB_tr_file = fopen("Aminus1_ByB_tr.xls","wt");
	FILE * Aminus1_ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(Aminus1_ByB_tr_file,"%f \n",Aminus1_By_B_tr[ind]);
	fclose(Aminus1_ByB_tr_file);
	FilePrint2D("Aminus1_ByB_tr_2D.xls",Aminus1_By_B_tr,Ny + 1, Nx);
#endif
	//double temp1, temp2;
	temp1 = norm_l2(Aminus1_By_B_tr,Ny + 1, Nx);
	temp2 = norm_max(Aminus1_By_B_tr,Ny + 1, Nx);
	printf("norm_l2 for Aminus1_By_B_tr = %f \n", temp1);
	printf("norm_max for Aminus1_By_B_tr = %f \n", temp2);

	double *BxBy_tr = (double*)malloc(dimWx * sizeof(double));

	int dimT = Nx * Ny;
	double *By_tr_Aminus1_ByB_tr = (double*)malloc(dimT * sizeof(double));

	//printf("Calculating By_tr_Aminus1_ByB_tr and BxBy_tr begins ... \n");

	//computing BxBy_tr * Ay(-1) * ByB_tr W
	//...

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx + 1; i++)
		{
			BxBy_tr[j*(Nx + 1) + i] = 0.0;
			if (i > 0)
				m_iminus1 = heat_capacity(i-1,j,0)*density(i-1,j,0)*hy(j)*hx(i-1);
			if (i < Nx)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (i > 0 && i < Nx)
			{
				BxBy_tr[j*(Nx + 1) + i] = hy(j)*(hx(i-1)*(1.0/m_iminus1)*(Aminus1_By_B_tr[j + 1 + (i-1)*(Ny+1)] - Aminus1_By_B_tr[j + (i-1)*(Ny+1)]) - hx(i)*(1.0/m_i)*(Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)]));
			}
			else if (i == 0)
			{
				BxBy_tr[j*(Nx + 1) + i] = 0.0;
			}
			else
			{
				BxBy_tr[j*(Nx + 1) + i] = 0.0;
			}

			if (i < Nx)
				By_tr_Aminus1_ByB_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i)*(Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)]));
		}
	}

	//return 0.0;


	//printf("Calculating BxBy_tr finishes ... \n");

#ifdef DEBUGE
	strcpy(filename,"BxBy_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * BxBy_tr_file = fopen("BxBy_tr.xls","wt");
	FILE * BxBy_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWx; ind++)
		fprintf(BxBy_tr_file,"%f \n",BxBy_tr[ind]);
	fclose(BxBy_tr_file);
	FilePrint2D("BxBy_tr_2D.xls",BxBy_tr,Nx + 1, Ny);

	strcpy(filename,"By_tr_Aminus1_ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * By_tr_Aminus1_ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimT; ind++)
		fprintf(By_tr_Aminus1_ByB_tr_file,"%f \n",By_tr_Aminus1_ByB_tr[ind]);
	fclose(By_tr_Aminus1_ByB_tr_file);
	FilePrint2D("By_tr_Aminus1_ByB_tr_2D.xls",By_tr_Aminus1_ByB_tr, Nx, Ny);

#endif
	temp1 = norm_l2(By_tr_Aminus1_ByB_tr, Ny, Nx);
	temp2 = norm_max(By_tr_Aminus1_ByB_tr, Ny, Nx);
	printf("norm_l2 for By_tr_Aminus1_By_B_tr = %f \n", temp1);
	printf("norm_max for By_tr_Aminus1_By_B_tr = %f \n", temp2);

	//printf("Calling norm_l2 begins ... \n");
	l2_norm = norm_l2(BxBy_tr, Nx + 1, Ny);
	//return 0.0;
	//printf("Calling norm_l2 finishes ... \n");

	//printf("Calling norm_max begins ... \n");
	max_norm = norm_max(BxBy_tr, Nx + 1, Ny);
	//return 0.0;
	//printf("Calling norm_max finishes ... \n");

	free(Aminus1_By_B_tr);
	free(ByB_tr);
	free(BxBy_tr);

	*l2_norm_pt = l2_norm;
	*max_norm_pt = max_norm;
}

int	myBtr_w0 ( double * Btr_w0, double * wx, double * wy, int Nx, int Ny )
{
	int dimWx = ( Nx + 1 ) *  Ny;
	int dimWy = Nx * ( Ny + 1 );
	int dimT = Nx * Ny;

	//computing B_tr W
	double m_iminus1, m_i;
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			Btr_w0 [ j * Nx + i ] = hy(j) * ( wx[j*(Nx+1) + i + 1] - wx[j*(Nx+1) + i] ) + hx(i) * ( wy[i*(Ny+1) + j + 1] - wy[i*(Ny+1) + j]);
		}
	}
#ifdef DEBUG_MINIM
	double temp1, temp2;
	temp1 = norm_l2(Btr_w0, Ny, Nx);
	temp2 = norm_max(Btr_w0, Ny, Nx);
	printf("norm_l2 for B_tr = %f \n", temp1);
	printf("norm_max for B_tr = %f \n", temp2);
#endif
	return 0;
}

int	myByMminus1Btr_w0 ( double * ByMminus1_Btr_w0, double * wx, double * wy, int Nx, int Ny )
{
	int dimWx = ( Nx + 1 ) *  Ny;
	int dimWy = Nx * ( Ny + 1 );
	int dimT = Nx * Ny;


	double * Btr_w0 = ( double * ) malloc ( dimT * sizeof ( double ) );
	myBtr_w0 ( Btr_w0, wx, wy, Nx, Ny );
//
////	double *ByB_tr = (double*)malloc(dimWy * sizeof(double));
//	//computing ByB_tr W
//	double m_iminus1, m_i;
//	for (int i = 0; i < Nx; i++)
//	{
//		for (int j = 0; j < Ny + 1; j++)
//		{
//			ByMminus1_Btr_w0[i*(Ny + 1) + j] = 0.0;
//			if (j > 0)
//				m_iminus1 = heat_capacity(i,j-1,0)*density(i,j-1,0)*hy(j-1)*hx(i);
//			if (j < Ny)
//				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);
//
//			if (j > 0 && j < Ny)
//			{
//				ByMminus1_Btr_w0[i*(Ny + 1) + j] = hx(i)*((1.0/m_iminus1)* Btr_w0 [ (j - 1) * Nx + i ]  - (1.0/m_i)* Btr_w0 [ j * Nx + i ] );
//			}
//			else if (j == 0)
//			{
//				ByMminus1_Btr_w0[i*(Ny + 1) + j] = hx(i)*( - (1.0/m_i)* Btr_w0 [ j * Nx + i ] );
//			}
//			else
//			{
//				ByMminus1_Btr_w0[i*(Ny + 1) + j] = hx(i)*((1.0/m_iminus1)* Btr_w0 [ (j - 1) * Nx + i ] );
//			}
//		}
//	}

//	double *ByB_tr = (double*)malloc(dimWy * sizeof(double));
	//computing ByB_tr W
	double m_iminus1, m_i;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			ByMminus1_Btr_w0[i*(Ny + 1) + j] = 0.0;
			if (j > 0)
				m_iminus1 = heat_capacity(i,j-1,0)*density(i,j-1,0)*hy(j-1)*hx(i);
			if (j < Ny)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (j > 0 && j < Ny)
			{
				ByMminus1_Btr_w0[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1])  - (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByMminus1_Btr_w0[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)])  - hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else if (j == 0)
			{
				ByMminus1_Btr_w0[i*(Ny + 1) + j] = hx(i)*hx(i)*(- (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByMminus1_Btr_w0[i*(Ny + 1) + j] += hx(i)*(- hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else
			{
				ByMminus1_Btr_w0[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1]));
				ByMminus1_Btr_w0[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)]) );
			}
		}
	}

	free ( Btr_w0 );
	return 0;
}

int	myAyminus1ByMminus1Btr_w0 ( double * Ayminus1ByMminus1_Btr_w0, double * wx, double * wy, int Nx, int Ny )
{
	int dimWx = ( Nx + 1 ) *  Ny;
	int dimWy = Nx * ( Ny + 1 );
	int dimT = Nx * Ny;

	double *ByB_tr = (double*)malloc(dimWy * sizeof(double));
	myByMminus1Btr_w0 ( ByB_tr, wx, wy, Nx, Ny );

//	double * Aminus1_By_B_tr = (double*)malloc(dimWy * sizeof(double));

	//printf("Calculating Aminus1_By_B_tr begins ... \n");

	double *alfay = (double*)malloc((Ny + 1)* sizeof(double));
	double *betay = (double*)malloc((Ny + 1)* sizeof(double));
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	double temper_0, temper_j, temper_n;
	for ( int i = 0 ; i < Nx ; i++ )
	{

		//к-ты матрицы А
		a_0_old = hy(0)*hx(i)/(heat_conductivity(i,0,0)*6.0);
		b_0_old = 2*a_0_old;
		temper_0 = ByB_tr[i*(Ny + 1) + 0];

		betay[0] = 0 ;
		alfay[0] = 0 ;
		alfay[1] = -0.5 ;
		betay[1] =  temper_0 / b_0_old; 

		for ( int j = 1 ; j < Ny ; j++ )
		{
			//к-ты матрицы А
			a_i_old = hy(j)*hx(i)/(heat_conductivity(i,j,0)*6.0);
			a_iminus1_old = hy(j-1)*hx(i)/(heat_conductivity(i,j-1,0)*6.0);
			b_i_old = 2*a_iminus1_old + 2*a_i_old;

			temper_j = ByB_tr[i*(Ny + 1) + j];

			alfay[j + 1] = (-1)*a_i_old/(alfay[j]*a_iminus1_old + b_i_old);
			betay[j + 1] = ( ( temper_j ) - betay[j]*a_iminus1_old ) / ( alfay[j]*a_iminus1_old + b_i_old );

		}
		a_nminus1_old = hy(Ny-1)*hx(i)/(heat_conductivity(i,Ny-1,0)*6.0);
		b_n_old = 2*a_nminus1_old;

		temper_n = ByB_tr[i*(Ny + 1) + Ny];

		Ayminus1ByMminus1_Btr_w0[i*(Ny+1) + Ny] = ( temper_n  -  betay[Ny]*a_nminus1_old ) / (alfay[Ny]*a_nminus1_old + b_n_old);

		for ( int j = 1 ; j < Ny + 1 ; j++ )
		{
			Ayminus1ByMminus1_Btr_w0[i*(Ny+1) + Ny - j ] = Ayminus1ByMminus1_Btr_w0[i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
		}
	}
	free(alfay);
	free(betay);

	free (ByB_tr);
//	free (Aminus1_By_B_tr);

	return 0;
}

int	myMhalfLambdayBtr_w0 ( double * MhalfLambdaY_Btr_w0, double * wx, double * wy, int Nx, int Ny )
{
	int dimWx = ( Nx + 1 ) *  Ny;
	int dimWy = Nx * ( Ny + 1 );
	int dimT = Nx * Ny;

	//double *ByB_tr = (double*)malloc(dimWy * sizeof(double));

	//myByMminus1Btr_w0 ( ByB_tr, wx, wy, Nx, Ny );

	////computing ByB_tr W
	//double m_iminus1, m_i;
	//for (int i = 0; i < Nx; i++)
	//{
	//	for (int j = 0; j < Ny + 1; j++)
	//	{
	//		ByB_tr[i*(Ny + 1) + j] = 0.0;
	//		if (j > 0)
	//			m_iminus1 = heat_capacity(i,j-1,0)*density(i,j-1,0)*hy(j-1)*hx(i);
	//		if (j < Ny)
	//			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

	//		if (j > 0 && j < Ny)
	//		{
	//			ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1])  - (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
	//			ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)])  - hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
	//		}
	//		else if (j == 0)
	//		{
	//			ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*(- (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
	//			ByB_tr[i*(Ny + 1) + j] += hx(i)*(- hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
	//		}
	//		else
	//		{
	//			ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1]));
	//			ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)]) );
	//		}
	//	}
	//}

//	//printf("Calculating ByB_tr finishes ... \n");
//
//#ifdef DEBUGE
//	strcpy(filename,"ByB_tr");
//	strcat(filename,nx_string);
//	strcat(filename,ext);
//	//FILE * ByB_tr_file = fopen("ByB_tr.xls","wt");
//	FILE * ByB_tr_file = fopen(filename,"wt");
//	for (int ind = 0; ind < dimWy; ind++)
//		fprintf(ByB_tr_file,"%f \n",ByB_tr[ind]);
//	fclose(ByB_tr_file);
//	FilePrint2D("ByB_tr_2D.xls",ByB_tr,Ny + 1, Nx);
//#endif
//
//#ifdef DEBUG_MINIM
//	double temp1, temp2;
//	temp1 = norm_l2(ByB_tr,Ny + 1, Nx);
//	temp2 = norm_max(ByB_tr,Ny + 1, Nx);
//	printf("norm_l2 for By_B_tr = %f \n", temp1);
//	printf("norm_max for By_B_tr = %f \n", temp2);
//#endif
	double * Aminus1_By_B_tr = (double*)malloc(dimWy * sizeof(double));

	myAyminus1ByMminus1Btr_w0 ( Aminus1_By_B_tr, wx, wy, Nx, Ny );

	////printf("Calculating Aminus1_By_B_tr begins ... \n");

	//double *alfay = (double*)malloc((Ny + 1)* sizeof(double));
	//double *betay = (double*)malloc((Ny + 1)* sizeof(double));
	//double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	//double temper_0, temper_j, temper_n;
	//for ( int i = 0 ; i < Nx ; i++ )
	//{

	//	//к-ты матрицы А
	//	a_0_old = hy(0)*hx(i)/(heat_conductivity(i,0,0)*6.0);
	//	b_0_old = 2*a_0_old;
	//	temper_0 = ByB_tr[i*(Ny + 1) + 0];

	//	betay[0] = 0 ;
	//	alfay[0] = 0 ;
	//	alfay[1] = -0.5 ;
	//	betay[1] =  temper_0 / b_0_old; 

	//	for ( int j = 1 ; j < Ny ; j++ )
	//	{
	//		//к-ты матрицы А
	//		a_i_old = hy(j)*hx(i)/(heat_conductivity(i,j,0)*6.0);
	//		a_iminus1_old = hy(j-1)*hx(i)/(heat_conductivity(i,j-1,0)*6.0);
	//		b_i_old = 2*a_iminus1_old + 2*a_i_old;

	//		temper_j = ByB_tr[i*(Ny + 1) + j];

	//		alfay[j + 1] = (-1)*a_i_old/(alfay[j]*a_iminus1_old + b_i_old);
	//		betay[j + 1] = ( ( temper_j ) - betay[j]*a_iminus1_old ) / ( alfay[j]*a_iminus1_old + b_i_old );

	//	}
	//	a_nminus1_old = hy(Ny-1)*hx(i)/(heat_conductivity(i,Ny-1,0)*6.0);
	//	b_n_old = 2*a_nminus1_old;

	//	temper_n = ByB_tr[i*(Ny + 1) + Ny];

	//	Aminus1_By_B_tr[i*(Ny+1) + Ny] = ( temper_n  -  betay[Ny]*a_nminus1_old ) / (alfay[Ny]*a_nminus1_old + b_n_old);

	//	for ( int j = 1 ; j < Ny + 1 ; j++ )
	//	{
	//		Aminus1_By_B_tr[i*(Ny+1) + Ny - j ] = Aminus1_By_B_tr[i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
	//	}
	//}
	//free(alfay);
	//free(betay);


	//double temp1, temp2;
//#ifdef DEBUG_MINIM
//	temp1 = norm_l2(Aminus1_By_B_tr,Ny + 1, Nx);
//	temp2 = norm_max(Aminus1_By_B_tr,Ny + 1, Nx);
//	printf("norm_l2 for Aminus1_By_B_tr = %f \n", temp1);
//	printf("norm_max for Aminus1_By_B_tr = %f \n", temp2);
//#endif

	//printf("Calculating By_tr_Aminus1_ByB_tr w0 begins ... \n");
	//printf("Final right \n");

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			MhalfLambdaY_Btr_w0 [ j * Nx + i ] = hx(i) * ( Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)] );
		}
	}

#ifdef DEBUG_MINIM
	double temp1, temp2;
	temp1 = norm_l2(MhalfLambdaY_Btr_w0, Ny, Nx);
	printf("norm_l2 for MhalfLambdaY_Btr_w0 = %f \n", temp1);
	temp2 = norm_max(MhalfLambdaY_Btr_w0, Ny, Nx);
	printf("norm_max for MhalfLambdaY_Btr_w0 = %f \n", temp2);
#endif
	//free (ByB_tr);
	free (Aminus1_By_B_tr);

	return 0;
}

int	myMminushalfLambdaLambdayBtr_w0 ( double * Mminus1Btr_Aminus1_xvost_w0, double * MhalfLambdaY_Btr_w0, double * wx, double * wy, int Nx, int Ny )
{
	int dimWx = ( Nx + 1 ) * Ny;
	int dimWy = Nx * ( Ny + 1 );
	int dimT = Nx * Ny;

	double temp1, temp2;
	double m_i;

	double * MminushalfLambdaY_Btr_w0 = (double*)malloc(dimT * sizeof(double)); // = M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			MminushalfLambdaY_Btr_w0 [ j * Nx + i ] = (1.0 / m_i ) * MhalfLambdaY_Btr_w0 [ j * Nx + i ];
		}
	}
#ifdef DEBUG_MINIM
	temp1 = norm_l2(MminushalfLambdaY_Btr_w0, Ny, Nx);
	printf("norm_l2 for MminushalfLambdaY_Btr_w0 = %f \n", temp1);
	temp2 = norm_max(MminushalfLambdaY_Btr_w0, Ny, Nx);
	printf("norm_max for MminushalfLambdaY_Btr_w0 = %f \n", temp2);
#endif
	//printf("Calculating Axminus1_Bx_Mminus1_By_tr_Aminus1_By_B_tr w0 begins ... \n");
	double * Axminus1_xvost_w0 = (double*)malloc(dimWx * sizeof(double)); // = Ax(-1) * Bx M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	HeatFlux_Wx0_Init(MminushalfLambdaY_Btr_w0, Axminus1_xvost_w0, xCondition, numSolution, Nx, Ny, 1);
#ifdef DEBUG_MINIM
	temp1 = norm_l2(Axminus1_xvost_w0, Nx + 1, Ny);
	printf("norm_l2 for Axminus1_xvost_w0 = %f \n", temp1);
	temp2 = norm_max(Axminus1_xvost_w0, Nx+1, Ny);
	printf("norm_max for Axminus1_xvost_w0 = %f \n", temp2);
#endif
	double * Mminus1Bxtr_Axminus1_xvost_w0 = (double*)malloc(dimT * sizeof(double)); // = M(-1) Bx_tr * Ax(-1) * Bx M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			Mminus1Bxtr_Axminus1_xvost_w0 [ j * Nx + i ] = hy(j) * ( 1.0 / m_i ) * ( Axminus1_xvost_w0[i + 1 + j*(Nx+1)] - Axminus1_xvost_w0[i + j*(Nx+1)] );
		}
	}
#ifdef DEBUG_MINIM
	temp1 = norm_l2(Mminus1Bxtr_Axminus1_xvost_w0, Ny, Nx);
	printf("norm_l2 for Mminus1Bxtr_Axminus1_xvost_w0 = %f \n", temp1);
	temp2 = norm_max(Mminus1Bxtr_Axminus1_xvost_w0, Ny, Nx);
	printf("norm_max for Mminus1Bxtr_Axminus1_xvost_w0 = %f \n", temp2);
#endif

	double * Ayminus1_xvost_w0 = (double*)malloc(dimWy * sizeof(double)); // = Ay(-1) * By M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	HeatFlux_Wy0_Init(MminushalfLambdaY_Btr_w0, Ayminus1_xvost_w0, yCondition, numSolution, Nx, Ny, 1);

//#ifdef DEBUG_MINIM
//	temp1 = norm_l2(Ayminus1_xvost_w0, Nx, Ny + 1);
//	printf("norm_l2 for Ayminus1_xvost_w0 = %f \n", temp1);
//	temp2 = norm_max(Ayminus1_xvost_w0, Nx, Ny+1);
//	printf("norm_max for Ayminus1_xvost_w0 = %f \n", temp2);
//#endif
	double * Mminus1Bytr_Ayminus1_xvost_w0 = (double*)malloc(dimT * sizeof(double)); // = M(-1) By_tr * Ay(-1) * By M(-1) * By_tr * Aminus1_By_B_tr, но надо M

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] = hx(i) * ( 1.0 / m_i ) * ( Ayminus1_xvost_w0[j + 1 + i*(Ny+1)] - Ayminus1_xvost_w0[j + i*(Ny+1)] );
		}
	}
//#ifdef DEBUG_MINIM
//	temp1 = norm_l2(Mminus1Bytr_Ayminus1_xvost_w0, Ny, Nx);
//	printf("norm_l2 for Mminus1Bytr_Ayminus1_xvost_w0 = %f \n", temp1);
//	temp2 = norm_max(Mminus1Bytr_Ayminus1_xvost_w0, Ny, Nx);
//	printf("norm_max for Mminus1Bytr_Ayminus1_xvost_w0 = %f \n", temp2);
//#endif
	//printf("Final left \n");
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ] = Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] + Mminus1Bxtr_Axminus1_xvost_w0 [ j * Nx + i ]; 
		}
	}
#ifdef DEBUG_MINIM
	temp1 = norm_l2(Mminus1Btr_Aminus1_xvost_w0, Ny, Nx);
	printf("norm_l2 for Mminus1Bytr_Aminus1_xvost_w0 = %f \n", temp1);
	temp2 = norm_max(Mminus1Btr_Aminus1_xvost_w0, Ny, Nx);
	printf("norm_max for Mminus1Bytr_Aminus1_xvost_w0 = %f \n", temp2);
#endif

	free ( MminushalfLambdaY_Btr_w0 );
	free ( Axminus1_xvost_w0 );
	free ( Ayminus1_xvost_w0 );
	free ( Mminus1Bytr_Ayminus1_xvost_w0 );
	free ( Mminus1Bxtr_Axminus1_xvost_w0 );

	return 0;
}


int	myByMminus1LambdaLambdayBtr_w0 ( double* ByMminus1LambdaLambday_Btr_w0, double* Mminus1Btr_Aminus1_xvost_w0, int Nx, int Ny )
{
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	double temp1, temp2;
	double m_i;

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j < Ny + 1; j++ )
		{
			if ( j > 0 && j < Ny )
			{
				ByMminus1LambdaLambday_Btr_w0 [ i * ( Ny + 1 ) + j ] = hx(i) * ( - Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ] + Mminus1Btr_Aminus1_xvost_w0 [ (j-1) * Nx + i ] );
			}

			if ( j == 0 )
			{
				ByMminus1LambdaLambday_Btr_w0 [ i * ( Ny + 1 ) + j ] = hx(i) * ( - Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ] );
			}

			if ( j == Ny )
			{
				ByMminus1LambdaLambday_Btr_w0 [ i * ( Ny + 1 ) + j ] = hx(i) * ( Mminus1Btr_Aminus1_xvost_w0 [ (j-1) * Nx + i ] );
			}

		}
	}

	return 0;
}


int	myAyminus1ByMminus1LambdaLambdayBtr_w0 (double* Ayminus1ByMminus1LambdaLambday_Btr_w0, double* Mminus1Btr_Aminus1_xvost_w0, int Nx, int Ny )
{
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	double temp1, temp2;
	double m_i;

	//for (int i = 0; i < Nx; i++)
	//{
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);
	//		Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ] = m_i  * Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ];
	//		//Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] = hx(i) * ( Ayminus1_xvost_w0[j + 1 + i*(Ny+1)] - Ayminus1_xvost_w0[j + i*(Ny+1)] );
	//	}
	//}

	//double * Ayminus1_xvost_w0 = (double*)malloc(dimWy * sizeof(double)); 

	HeatFlux_Wy0_Init(Mminus1Btr_Aminus1_xvost_w0, Ayminus1ByMminus1LambdaLambday_Btr_w0, yCondition, numSolution, Nx, Ny, 1);

	return 0;
}


int	myMminus1BytrAyminus1ByMminus1LambdaLambdayBtr_w0 ( double * Y3, double * Mminus1Btr_Aminus1_xvost_w0, int Nx, int Ny )
{
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;
	double m_i;

	double * Ayminus1ByMminus1LambdaLambday_Btr_w0 = ( double * ) malloc ( dimWy * sizeof ( double ) );
	myAyminus1ByMminus1LambdaLambdayBtr_w0 (Ayminus1ByMminus1LambdaLambday_Btr_w0, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny );

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);
			Y3 [ j * Nx + i ] = hx(i) * ( 1.0 / m_i ) * ( Ayminus1ByMminus1LambdaLambday_Btr_w0[j + 1 + i*(Ny+1)] - Ayminus1ByMminus1LambdaLambday_Btr_w0[j + i*(Ny+1)] );
		}
	}



	return 0;
}


int	myBLambdayLambdaLambdayBtr_w0 ( double* tempvec2X, double* tempvec2Y, double* Mminus1Btr_Aminus1_xvost_w0, int Nx, int Ny )
{
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	double temp1, temp2;
	double m_i;

//	double * Y3 = ( double * ) malloc ( dimWy * sizeof ( double ) );
	double * Ayminus1_xvost_w0 = (double*)malloc(dimWy * sizeof(double)); 
	myAyminus1ByMminus1LambdaLambdayBtr_w0 ( Ayminus1_xvost_w0, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny );

	////for (int i = 0; i < Nx; i++)
	////{
	////	for (int j = 0; j < Ny; j++)
	////	{
	////		m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);
	////		Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ] = m_i  * Mminus1Btr_Aminus1_xvost_w0 [ j * Nx + i ];
	////		//Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] = hx(i) * ( Ayminus1_xvost_w0[j + 1 + i*(Ny+1)] - Ayminus1_xvost_w0[j + i*(Ny+1)] );
	////	}
	////}

	//double * Ayminus1_xvost_w0 = (double*)malloc(dimWy * sizeof(double)); 

	//HeatFlux_Wy0_Init(Mminus1Btr_Aminus1_xvost_w0, Ayminus1_xvost_w0, yCondition, numSolution, Nx, Ny, 1);

	double * Mminus1Bytr_Ayminus1_xvost_w0 = (double*)malloc(dimT * sizeof(double)); 
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);
			Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] = hx(i) * ( 1.0 / m_i ) * ( Ayminus1_xvost_w0[j + 1 + i*(Ny+1)] - Ayminus1_xvost_w0[j + i*(Ny+1)] );
			//Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] = hx(i) * ( Ayminus1_xvost_w0[j + 1 + i*(Ny+1)] - Ayminus1_xvost_w0[j + i*(Ny+1)] );
		}
	}

	// calculating tempvec = B * Mminus1Bytr_Ayminus1_xvost_w0
	for ( int j = 0; j < Ny; j++ )
	{
		for ( int i = 0; i < Nx + 1; i++ )
		{
			if ( i > 0 && i < Nx )
			{
				tempvec2X [ j * ( Nx + 1 ) + i ] = hy(j) * (  - Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] + Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i - 1 ] );
			}

			if ( i == 0 )
			{
				tempvec2X [ j * ( Nx + 1 ) + i ] = hy(j) * ( - Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] );
			}

			if ( i == Nx )
			{
				tempvec2X [ j * ( Nx + 1 ) + i ] = hy(j) * ( Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i - 1 ] );
			}

		}
	}

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j < Ny + 1; j++ )
		{
			if ( j > 0 && j < Ny )
			{
				tempvec2Y [ i * ( Ny + 1 ) + j ] = hx(i) * ( - Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] + Mminus1Bytr_Ayminus1_xvost_w0 [ (j-1) * Nx + i ] );
			}

			if ( j == 0 )
			{
				tempvec2Y [ i * ( Ny + 1 ) + j ] = hx(i) * ( - Mminus1Bytr_Ayminus1_xvost_w0 [ j * Nx + i ] );
			}

			if ( j == Ny )
			{
				tempvec2Y [ i * ( Ny + 1 ) + j ] = hx(i) * ( Mminus1Bytr_Ayminus1_xvost_w0 [ (j-1) * Nx + i ] );
			}

		}
	}


	free ( Ayminus1_xvost_w0 );
	free ( Mminus1Bytr_Ayminus1_xvost_w0 );

	return 0;
}



void stability_normFull(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	double temp1, temp2;

	double * MhalfLambdaY_Btr_w0 = (double*)malloc(dimT * sizeof(double));
	myMhalfLambdayBtr_w0 ( MhalfLambdaY_Btr_w0, wx, wy, Nx, Ny );

#ifdef DEBUG_MINIM
	temp1 = norm_l2(MhalfLambdaY_Btr_w0, Ny, Nx);
	printf("norm_l2 for MhalfLambdaY_Btr_w0 = %f \n", temp1);
	temp2 = norm_max(MhalfLambdaY_Btr_w0, Ny, Nx);
	printf("norm_max for MhalfLambdaY_Btr_w0 = %f \n", temp2);
#endif

	double * Mminus1Btr_Aminus1_xvost_w0 = (double*)malloc(dimT * sizeof(double));  // = M(-1/2) * Lambda * Lambda_y * Btr * w0
	myMminushalfLambdaLambdayBtr_w0 ( Mminus1Btr_Aminus1_xvost_w0, MhalfLambdaY_Btr_w0, wx, wy, Nx, Ny );

//	temp1 = scal_l2 ( Mminus1Btr_Aminus1_xvost_w0, MhalfLambdaY_Btr_w0 , Nx, Ny ) / dimT;
//	printf ( "temp1 = %f \n", temp1 );

	//l2_norm = sqrt ( scal_l2 ( Mminus1Btr_Aminus1_xvost_w0, MhalfLambdaY_Btr_w0 , Nx, Ny ) / (Nx * Ny) );
	l2_norm = sqrt ( scal_l2 ( Mminus1Btr_Aminus1_xvost_w0, MhalfLambdaY_Btr_w0 , Nx, Ny ) );
#ifdef DEBUG_MINIM
	printf ( "full stability norm = %f \n", l2_norm );
#endif
	free ( MhalfLambdaY_Btr_w0 );
	free ( Mminus1Btr_Aminus1_xvost_w0 );

	*l2_norm_pt = l2_norm;
}

void stability_norm_Matlab(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
#ifdef DEBUGE
	char* ext = ".xls";
	char nx_string[5];
	sprintf(nx_string,"_%d",Nx);
	char filename[30];
	strcpy(filename,"Wy");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * Wy_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(Wy_file,"%f \n",wy[ind]);
	fclose(Wy_file);
	FilePrint2D("Wy_2D.xls",wy,Ny + 1, Nx);

	strcpy(filename,"Wx");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * Wx_file = fopen("Wx.xls","wt");
	FILE * Wx_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWx; ind++)
		fprintf(Wx_file,"%f \n",wx[ind]);
	fclose(Wx_file);
	FilePrint2D("Wx_2D.xls",wx,Nx + 1, Ny);
#endif
	//printf("Calculating ByB_tr begins ... \n");

	double *ByB_tr = (double*)malloc(dimWy * sizeof(double));
	//computing ByB_tr W
	double m_iminus1, m_i;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			ByB_tr[i*(Ny + 1) + j] = 0.0;
			if (j > 0)
				m_iminus1 = heat_capacity(i,j-1,0)*density(i,j-1,0)*hy(j-1)*hx(i);
			if (j < Ny)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (j > 0 && j < Ny)
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1])  - (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)])  - hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else if (j == 0)
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*(- (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(- hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)]) );
			}
		}
	}

	//printf("Calculating ByB_tr finishes ... \n");

#ifdef DEBUGE
	strcpy(filename,"ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * ByB_tr_file = fopen("ByB_tr.xls","wt");
	FILE * ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(ByB_tr_file,"%f \n",ByB_tr[ind]);
	fclose(ByB_tr_file);
	FilePrint2D("ByB_tr_2D.xls",ByB_tr,Ny + 1, Nx);
#endif
	double temp1, temp2;
	temp1 = norm_l2(ByB_tr,Ny + 1, Nx);
	temp2 = norm_max(ByB_tr,Ny + 1, Nx);
	printf("norm_l2 for By_B_tr = %f \n", temp1);
	printf("norm_max for By_B_tr = %f \n", temp2);


	//return 0.0;

	double * Aminus1_By_B_tr = (double*)malloc(dimWy * sizeof(double));
	//computing Ay(-1) * BxBy_tr W 
	//...

	//printf("Calculating Aminus1_By_B_tr begins ... \n");

	double *alfay = (double*)malloc((Ny + 1)* sizeof(double));
	double *betay = (double*)malloc((Ny + 1)* sizeof(double));
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	double temper_0, temper_j, temper_n;
	for ( int i = 0 ; i < Nx ; i++ )
	{

		//к-ты матрицы А
		a_0_old = hy(0)*hx(i)/(heat_conductivity(i,0,0)*6.0);
		b_0_old = 2*a_0_old;
		temper_0 = ByB_tr[i*(Ny + 1) + 0];

		betay[0] = 0 ;
		alfay[0] = 0 ;
		alfay[1] = -0.5 ;
		betay[1] =  temper_0 / b_0_old; 

		for ( int j = 1 ; j < Ny ; j++ )
		{
			//к-ты матрицы А
			a_i_old = hy(j)*hx(i)/(heat_conductivity(i,j,0)*6.0);
			a_iminus1_old = hy(j-1)*hx(i)/(heat_conductivity(i,j-1,0)*6.0);
			b_i_old = 2*a_iminus1_old + 2*a_i_old;

			temper_j = ByB_tr[i*(Ny + 1) + j];

			alfay[j + 1] = (-1)*a_i_old/(alfay[j]*a_iminus1_old + b_i_old);
			betay[j + 1] = ( ( temper_j ) - betay[j]*a_iminus1_old ) / ( alfay[j]*a_iminus1_old + b_i_old );

		}
		a_nminus1_old = hy(Ny-1)*hx(i)/(heat_conductivity(i,Ny-1,0)*6.0);
		b_n_old = 2*a_nminus1_old;

		temper_n = ByB_tr[i*(Ny + 1) + Ny];

		Aminus1_By_B_tr[i*(Ny+1) + Ny] = ( temper_n  -  betay[Ny]*a_nminus1_old ) / (alfay[Ny]*a_nminus1_old + b_n_old);

		for ( int j = 1 ; j < Ny + 1 ; j++ )
		{
			Aminus1_By_B_tr[i*(Ny+1) + Ny - j ] = Aminus1_By_B_tr[i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
		}
	}
	free(alfay);
	free(betay);

	//printf("Calculating Aminus1_By_B_tr finishes ... \n");

	//return 0.0;


#ifdef DEBUGE
	strcpy(filename,"Aminus1_ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * Aminus1_ByB_tr_file = fopen("Aminus1_ByB_tr.xls","wt");
	FILE * Aminus1_ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(Aminus1_ByB_tr_file,"%f \n",Aminus1_By_B_tr[ind]);
	fclose(Aminus1_ByB_tr_file);
	FilePrint2D("Aminus1_ByB_tr_2D.xls",Aminus1_By_B_tr,Ny + 1, Nx);
#endif
	//double temp1, temp2;
	temp1 = norm_l2(Aminus1_By_B_tr,Ny + 1, Nx);
	temp2 = norm_max(Aminus1_By_B_tr,Ny + 1, Nx);
	printf("norm_l2 for Aminus1_By_B_tr = %f \n", temp1);
	printf("norm_max for Aminus1_By_B_tr = %f \n", temp2);

	double *BxBy_tr = (double*)malloc(dimWx * sizeof(double));

	int dimT = Nx * Ny;
	double *By_tr_Aminus1_ByB_tr = (double*)malloc(dimT * sizeof(double));

	//printf("Calculating By_tr_Aminus1_ByB_tr and BxBy_tr begins ... \n");

	//computing BxBy_tr * Ay(-1) * ByB_tr W
	//...

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx + 1; i++)
		{
			BxBy_tr[j*(Nx + 1) + i] = 0.0;
			if (i > 0)
				m_iminus1 = heat_capacity(i-1,j,0)*density(i-1,j,0)*hy(j)*hx(i-1);
			if (i < Nx)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (i > 0 && i < Nx)
			{
				BxBy_tr[j*(Nx + 1) + i] = hy(j)*(hx(i-1)*(1.0/m_iminus1)*(Aminus1_By_B_tr[j + 1 + (i-1)*(Ny+1)] - Aminus1_By_B_tr[j + (i-1)*(Ny+1)]) - hx(i)*(1.0/m_i)*(Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)]));
			}
			else if (i == 0)
			{
				BxBy_tr[j*(Nx + 1) + i] = 0.0;
			}
			else
			{
				BxBy_tr[j*(Nx + 1) + i] = 0.0;
			}

			if (i < Nx)
				By_tr_Aminus1_ByB_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i)*(Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)]));
		}
	}

	//return 0.0;


	//printf("Calculating BxBy_tr finishes ... \n");

#ifdef DEBUGE
	strcpy(filename,"BxBy_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * BxBy_tr_file = fopen("BxBy_tr.xls","wt");
	FILE * BxBy_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWx; ind++)
		fprintf(BxBy_tr_file,"%f \n",BxBy_tr[ind]);
	fclose(BxBy_tr_file);
	FilePrint2D("BxBy_tr_2D.xls",BxBy_tr,Nx + 1, Ny);

	strcpy(filename,"By_tr_Aminus1_ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * By_tr_Aminus1_ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimT; ind++)
		fprintf(By_tr_Aminus1_ByB_tr_file,"%f \n",By_tr_Aminus1_ByB_tr[ind]);
	fclose(By_tr_Aminus1_ByB_tr_file);
	FilePrint2D("By_tr_Aminus1_ByB_tr_2D.xls",By_tr_Aminus1_ByB_tr, Nx, Ny);

#endif
	temp1 = norm_l2(By_tr_Aminus1_ByB_tr, Ny, Nx);
	temp2 = norm_max(By_tr_Aminus1_ByB_tr, Ny, Nx);
	printf("norm_l2 for By_tr_Aminus1_By_B_tr = %f \n", temp1);
	printf("norm_max for By_tr_Aminus1_By_B_tr = %f \n", temp2);

	//printf("Calling norm_l2 begins ... \n");
	l2_norm = norm_l2(BxBy_tr, Nx + 1, Ny);
	//return 0.0;
	//printf("Calling norm_l2 finishes ... \n");

	//printf("Calling norm_max begins ... \n");
	max_norm = norm_max(BxBy_tr, Nx + 1, Ny);
	//return 0.0;
	//printf("Calling norm_max finishes ... \n");

	free(Aminus1_By_B_tr);
	free(ByB_tr);
	free(BxBy_tr);

	*l2_norm_pt = l2_norm;
	*max_norm_pt = max_norm;
}

void stability_norm_sum(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt)
{
	double l2_norm_x = 0.0, max_norm_x = 0.0, l2_norm_y = 0.0, max_norm_y = 0.0, l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
#ifdef DEBUGE
	char* ext = ".xls";
	char nx_string[5];
	sprintf(nx_string,"_%d",Nx);
	char filename[30];
	strcpy(filename,"Wy");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * Wy_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(Wy_file,"%f \n",wy[ind]);
	fclose(Wy_file);
	FilePrint2D("Wy_2D.xls",wy,Ny + 1, Nx);

	strcpy(filename,"Wx");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * Wx_file = fopen("Wx.xls","wt");
	FILE * Wx_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWx; ind++)
		fprintf(Wx_file,"%f \n",wx[ind]);
	fclose(Wx_file);
	FilePrint2D("Wx_2D.xls",wx,Nx + 1, Ny);
#endif
	//printf("Calculating ByB_tr begins ... \n");

	double *ByB_tr = (double*)malloc(dimWy * sizeof(double));
	//computing ByB_tr W
	double m_iminus1, m_i;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			ByB_tr[i*(Ny + 1) + j] = 0.0;
			if (j > 0)
				m_iminus1 = heat_capacity(i,j-1,0)*density(i,j-1,0)*hy(j-1)*hx(i);
			if (j < Ny)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (j > 0 && j < Ny)
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1])  - (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)])  - hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else if (j == 0)
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*(- (1.0/m_i)*(wy[i*(Ny+1) + j+1] - wy[i*(Ny+1) + j]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(- hy(j)*(1.0/m_i)*(wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]));
			}
			else
			{
				ByB_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(wy[i*(Ny+1) + j] - wy[i*(Ny+1) + j-1]));
				ByB_tr[i*(Ny + 1) + j] += hx(i)*(hy(j-1)*(1.0/m_iminus1)*(wx[i + 1 + (j-1)*(Nx+1)] - wx[i + (j-1)*(Nx+1)]) );
			}
		}
	}

	//printf("Calculating ByB_tr finishes ... \n");

#ifdef DEBUGE
	strcpy(filename,"ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * ByB_tr_file = fopen("ByB_tr.xls","wt");
	FILE * ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(ByB_tr_file,"%f \n",ByB_tr[ind]);
	fclose(ByB_tr_file);
	FilePrint2D("ByB_tr_2D.xls",ByB_tr,Ny + 1, Nx);
#endif
	double temp1, temp2;
	temp1 = norm_l2(ByB_tr,Ny + 1, Nx);
	temp2 = norm_max(ByB_tr,Ny + 1, Nx);
	printf("norm_l2 for By_B_tr = %f \n", temp1);
	printf("norm_max for By_B_tr = %f \n", temp2);


	//return 0.0;

	double * Aminus1_By_B_tr = (double*)malloc(dimWy * sizeof(double));
	//computing Ay(-1) * BxBy_tr W 
	//...

	//printf("Calculating Aminus1_By_B_tr begins ... \n");

	double *alfay = (double*)malloc((Ny + 1)* sizeof(double));
	double *betay = (double*)malloc((Ny + 1)* sizeof(double));
	double a_0_old, a_iminus1_old, a_i_old, a_nminus1_old, b_0_old, b_i_old, b_n_old;
	double temper_0, temper_j, temper_n;
	for ( int i = 0 ; i < Nx ; i++ )
	{

		//к-ты матрицы А
		a_0_old = hy(0)*hx(i)/(heat_conductivity(i,0,0)*6.0);
		b_0_old = 2*a_0_old;
		temper_0 = ByB_tr[i*(Ny + 1) + 0];

		betay[0] = 0 ;
		alfay[0] = 0 ;
		alfay[1] = -0.5 ;
		betay[1] =  temper_0 / b_0_old; 

		for ( int j = 1 ; j < Ny ; j++ )
		{
			//к-ты матрицы А
			a_i_old = hy(j)*hx(i)/(heat_conductivity(i,j,0)*6.0);
			a_iminus1_old = hy(j-1)*hx(i)/(heat_conductivity(i,j-1,0)*6.0);
			b_i_old = 2*a_iminus1_old + 2*a_i_old;

			temper_j = ByB_tr[i*(Ny + 1) + j];

			alfay[j + 1] = (-1)*a_i_old/(alfay[j]*a_iminus1_old + b_i_old);
			betay[j + 1] = ( ( temper_j ) - betay[j]*a_iminus1_old ) / ( alfay[j]*a_iminus1_old + b_i_old );

		}
		a_nminus1_old = hy(Ny-1)*hx(i)/(heat_conductivity(i,Ny-1,0)*6.0);
		b_n_old = 2*a_nminus1_old;

		temper_n = ByB_tr[i*(Ny + 1) + Ny];

		Aminus1_By_B_tr[i*(Ny+1) + Ny] = ( temper_n  -  betay[Ny]*a_nminus1_old ) / (alfay[Ny]*a_nminus1_old + b_n_old);

		for ( int j = 1 ; j < Ny + 1 ; j++ )
		{
			Aminus1_By_B_tr[i*(Ny+1) + Ny - j ] = Aminus1_By_B_tr[i*(Ny+1) + Ny + 1 - j ] * alfay[Ny + 1 - j] + betay[Ny + 1 - j]; 
		}
	}
	free(alfay);
	free(betay);

	//printf("Calculating Aminus1_By_B_tr finishes ... \n");

	//return 0.0;


#ifdef DEBUGE
	strcpy(filename,"Aminus1_ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * Aminus1_ByB_tr_file = fopen("Aminus1_ByB_tr.xls","wt");
	FILE * Aminus1_ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(Aminus1_ByB_tr_file,"%f \n",Aminus1_By_B_tr[ind]);
	fclose(Aminus1_ByB_tr_file);
	FilePrint2D("Aminus1_ByB_tr_2D.xls",Aminus1_By_B_tr,Ny + 1, Nx);
#endif
	//double temp1, temp2;
	temp1 = norm_l2(Aminus1_By_B_tr,Ny + 1, Nx);
	temp2 = norm_max(Aminus1_By_B_tr,Ny + 1, Nx);
	printf("norm_l2 for Aminus1_By_B_tr = %f \n", temp1);
	printf("norm_max for Aminus1_By_B_tr = %f \n", temp2);

	double *BxBy_tr = (double*)malloc(dimWx * sizeof(double));
	double *ByBy_tr = (double*)malloc(dimWy * sizeof(double));

	int dimT = Nx * Ny;
	double *By_tr_Aminus1_ByB_tr = (double*)malloc(dimT * sizeof(double));

	//printf("Calculating By_tr_Aminus1_ByB_tr and BxBy_tr begins ... \n");

	//computing BxBy_tr * Ay(-1) * ByB_tr W
	//...

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx + 1; i++)
		{
			BxBy_tr[j*(Nx + 1) + i] = 0.0;
			if (i > 0)
				m_iminus1 = heat_capacity(i-1,j,0)*density(i-1,j,0)*hy(j)*hx(i-1);
			if (i < Nx)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (i > 0 && i < Nx)
			{
				BxBy_tr[j*(Nx + 1) + i] = hy(j)*(hx(i-1)*(1.0/m_iminus1)*(Aminus1_By_B_tr[j + 1 + (i-1)*(Ny+1)] - Aminus1_By_B_tr[j + (i-1)*(Ny+1)]) - hx(i)*(1.0/m_i)*(Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)]));
			}
			else if (i == 0)
			{
				BxBy_tr[j*(Nx + 1) + i] = 0.0;
			}
			else
			{
				BxBy_tr[j*(Nx + 1) + i] = 0.0;
			}

			if (i < Nx)
				By_tr_Aminus1_ByB_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i)*(Aminus1_By_B_tr[j + 1 + i*(Ny+1)] - Aminus1_By_B_tr[j + i*(Ny+1)]));
		}
	}


	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			ByBy_tr[i*(Ny + 1) + j] = 0.0;
			if (j > 0)
				m_iminus1 = heat_capacity(i-1,j,0)*density(i-1,j,0)*hy(j - 1)*hx(i);
			if (j < Ny)
				m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			if (j == 0)
			{
				ByBy_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*( - (1.0/m_i)*(Aminus1_By_B_tr[i*(Ny+1) + j+1] - Aminus1_By_B_tr[i*(Ny+1) + j]));
			}
			else
				if (j == Ny)
				{
					ByBy_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(Aminus1_By_B_tr[i*(Ny+1) + j] - Aminus1_By_B_tr[i*(Ny+1) + j-1]) );
				}
				else
				{
					ByBy_tr[i*(Ny + 1) + j] = hx(i)*hx(i)*((1.0/m_iminus1)*(Aminus1_By_B_tr[i*(Ny+1) + j] - Aminus1_By_B_tr[i*(Ny+1) + j-1])  - (1.0/m_i)*(Aminus1_By_B_tr[i*(Ny+1) + j+1] - Aminus1_By_B_tr[i*(Ny+1) + j]));
				}
		}
	}
	//return 0.0;


	//printf("Calculating BxBy_tr and ByBy_tr finishes ... \n");
	//exit;

#ifdef DEBUGE
	strcpy(filename,"BxBy_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	//FILE * BxBy_tr_file = fopen("BxBy_tr.xls","wt");
	FILE * BxBy_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWx; ind++)
		fprintf(BxBy_tr_file,"%f \n",BxBy_tr[ind]);
	fclose(BxBy_tr_file);
	FilePrint2D("BxBy_tr_2D.xls",BxBy_tr,Nx + 1, Ny);

	strcpy(filename,"By_tr_Aminus1_ByB_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * By_tr_Aminus1_ByB_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimT; ind++)
		fprintf(By_tr_Aminus1_ByB_tr_file,"%f \n",By_tr_Aminus1_ByB_tr[ind]);
	fclose(By_tr_Aminus1_ByB_tr_file);
	FilePrint2D("By_tr_Aminus1_ByB_tr_2D.xls",By_tr_Aminus1_ByB_tr, Nx, Ny);

	strcpy(filename,"ByBy_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * ByBy_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimWy; ind++)
		fprintf(ByBy_tr_file,"%f \n",ByBy_tr[ind]);
	fclose(ByBy_tr_file);
	FilePrint2D("ByBy_tr_2D.xls",ByBy_tr,Ny + 1, Nx);

#endif
	temp1 = norm_l2(By_tr_Aminus1_ByB_tr, Ny, Nx);
	temp2 = norm_max(By_tr_Aminus1_ByB_tr, Ny, Nx);
	printf("norm_l2 for By_tr_Aminus1_By_B_tr = %f \n", temp1);
	printf("norm_max for By_tr_Aminus1_By_B_tr = %f \n", temp2);

	//printf("Calling norm_l2 begins ... \n");
	//exit;
	l2_norm_x = norm_l2(BxBy_tr, Nx + 1, Ny);
	l2_norm_y = norm_l2(ByBy_tr, Ny + 1, Nx);
	l2_norm = sqrt(l2_norm_x * l2_norm_x + l2_norm_y * l2_norm_y);
	printf("l2_norm Bx_By_tr_A_minus1_By_B_tr = %f \n",l2_norm_x);
	printf("l2_norm By_By_tr_A_minus1_By_B_tr = %f \n",l2_norm_y);
	//l2_norm = l2_norm_y;
	//return 0.0;
	//printf("Calling norm_l2 finishes ... \n");

	//printf("Calling norm_max begins ... \n");
	max_norm_x = norm_max(BxBy_tr, Nx + 1, Ny);
	max_norm_y = norm_max(ByBy_tr, Ny + 1, Nx);
	printf("max_norm Bx_By_tr_A_minus1_By_B_tr = %f \n",max_norm_x);
	printf("max_norm By_By_tr_A_minus1_By_B_tr = %f \n",max_norm_y);
	max_norm = max_norm_x;
	if (max_norm_y >= max_norm)
		max_norm = max_norm_y;

	//max_norm = max_norm_y;
	//return 0.0;
	//printf("Calling norm_max finishes ... \n");

	free(Aminus1_By_B_tr);
	free(ByB_tr);
	free(BxBy_tr);
	free(ByBy_tr);

	*l2_norm_pt = l2_norm;
	*max_norm_pt = max_norm;
}


double exact_harm(int harm, double y)
{
	return sin(harm * PI * y);
}

double exact_gradient_harm(int harm, double y)
{
	return harm * PI * cos(harm * PI * y);
}
void divergence_Wy_norm(double* wy, int Ny, double *l2_norm_pt, double *max_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	//printf("Calculating B_tr begins ... \n");
	double *B_tr = (double*)malloc(Ny * sizeof(double));
	//computing B_tr W
	double m_i;
	for (int j = 0; j < Ny; j++)
	{
		m_i = heat_capacity(0,j,0)*density(0,j,0)*hy(j);

		//B_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i) * (wy[i*(Ny+1) + j + 1] - wy[i*(Ny+1) + j]) + hy(j) * (wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]) );
		B_tr[j] = (1.0/m_i)* (wy[j + 1] - wy[j]) ;
	}
	//printf("Calculating B_tr finishes ... \n");

	l2_norm = norm_l2(B_tr, Ny, 1);
	max_norm = norm_max(B_tr, Ny, 1);

#ifdef DEBUGE
	char filename[50];
	strcpy(filename,"diverg_Wy_SPECIAL_MODE_T");
	char nx_string[5];
	sprintf(nx_string,"_%d.xls",Nx);
	strcat(filename,nx_string);
	FILE *result_div = fopen(filename,"wt");
	for (int i = 0; i < Ny; i++)
		fprintf(result_div,"%f \n",B_tr[i]);
	fclose(result_div);
#endif
	//printf("l2_norm = %f \n", l2_norm);
	//printf("max_norm = %f \n", max_norm);

	free(B_tr);

	*l2_norm_pt = l2_norm;
	*max_norm_pt = max_norm;

}

void divergence_norm(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	//printf("Calculating B_tr begins ... \n");
	double *B_tr = (double*)malloc(dimT * sizeof(double));
	//computing B_tr W
	double m_i;
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			//B_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i) * (wy[i*(Ny+1) + j + 1] - wy[i*(Ny+1) + j]) + hy(j) * (wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]) );
			B_tr[j*Nx + i] = (1.0/m_i)*( hx(i) * (wy[i*(Ny+1) + j + 1] - wy[i*(Ny+1) + j]) + hy(j) * (wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]) );
		}
	}
	//printf("Calculating B_tr finishes ... \n");

#ifdef DEBUGE
	char* ext = ".xls";
	char nx_string[5];
	sprintf(nx_string,"_%d",Nx);
	char filename[30];
	strcpy(filename,"B_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * B_tr_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimT; ind++)
		fprintf(B_tr_file,"%f \n",B_tr[ind]);
	fclose(B_tr_file);

	char filename2[30];
	strcpy(filename2,"B_tr_yx");
	strcat(filename2,nx_string);
	strcat(filename2,ext);
	B_tr_file = fopen(filename2,"wt");
	for (int ind = 0; ind < Ny; ind++)
		fprintf(B_tr_file,"%f \n",B_tr[ind*Nx]);
	fclose(B_tr_file);

	FilePrint2D("B_tr_2D.xls",B_tr, Nx , Ny);

#endif

	//printf("Calling norm_l2 begins ... \n");
	l2_norm = norm_l2(B_tr, Nx, Ny);
	//return 0.0;
	//printf("Calling norm_l2 finishes ... \n");

	//printf("Calling norm_max begins ... \n");
	max_norm = norm_max(B_tr, Nx, Ny);
	//return 0.0;
	//printf("Calling norm_max finishes ... \n");

	free(B_tr);

	*l2_norm_pt = l2_norm;
	*max_norm_pt = max_norm;
}

void Bytr_By_divergence_norm(double* wx, double* wy, int Nx, int Ny, double *l2_norm_pt, double *max_norm_pt)
{
	double l2_norm = 0.0, max_norm = 0.0;
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	//printf("Calculating B_tr begins ... \n");
	double *B_tr = (double*)malloc(dimT * sizeof(double));
	//computing B_tr W
	double m_i;
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			m_i = heat_capacity(i,j,0)*density(i,j,0)*hy(j)*hx(i);

			B_tr[j*Nx + i] = (1.0/sqrt(m_i))*( hx(i) * (wy[i*(Ny+1) + j + 1] - wy[i*(Ny+1) + j]) + hy(j) * (wx[i + 1 + j*(Nx+1)] - wx[i + j*(Nx+1)]) );
		}
	}
	//printf("Calculating B_tr finishes ... \n");

	double *Bytr_By_B_tr = (double*)malloc(dimT * sizeof(double));



#ifdef DEBUGE
	char* ext = ".xls";
	char nx_string[5];
	sprintf(nx_string,"_%d",Nx);
	char filename[30];
	strcpy(filename,"Bytr_By_B_tr");
	strcat(filename,nx_string);
	strcat(filename,ext);
	FILE * Bytr_By_file = fopen(filename,"wt");
	for (int ind = 0; ind < dimT; ind++)
		fprintf(Bytr_By_file,"%f \n",Bytr_By_B_tr[ind]);
	fclose(Bytr_By_file);
	FilePrint2D("Bytr_By_B_tr_2D.xls",B_tr, Nx , Ny);

#endif

	//printf("Calling norm_l2 begins ... \n");
	l2_norm = norm_l2(B_tr, Nx, Ny);
	//return 0.0;
	//printf("Calling norm_l2 finishes ... \n");

	//printf("Calling norm_max begins ... \n");
	max_norm = norm_max(B_tr, Nx, Ny);
	//return 0.0;
	//printf("Calling norm_max finishes ... \n");

	free(B_tr);

	*l2_norm_pt = l2_norm;
	*max_norm_pt = max_norm;
}

double norm_l2(double * massiv, int N, int M)
{
	double norm = 0.0;
	int index = 0;
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			norm += massiv[index]*massiv[index];
			index += 1;
		}
	}
	norm = sqrt(norm / (N * M));
	return norm;
}


double scal_l2_bnd ( char * TorW, char * XorY, double * Vec1, double * Vec2, int Nx, int Ny, boundaryConditions xCondition, boundaryConditions yCondition  )
{
	double res = 0.0;

	if ( XorY == "X" && TorW == "W" )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			for ( int i = 0; i <= Nx; i++ )
			{
				if ( i == 0 || i == Nx )
					res += 0.5 * Vec1[j*(Nx+1)+i] * Vec2[j*(Nx+1)+i];
				else
					res += Vec1[j*(Nx+1)+i] * Vec2[j*(Nx+1)+i];
			}
		}
	}

	if ( XorY == "Y" && TorW == "W" )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			for ( int j = 0; j <= Ny; j++ )
			{
				if ( j == 0 || j == Ny )
					res += 0.5 * Vec1[i*(Ny+1)+j] * Vec2[i*(Ny+1)+j];
				else
					res += Vec1[i*(Ny+1)+j] * Vec2[i*(Ny+1)+j];
			}
		}
	}

	if ( TorW == "T" )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			for ( int i = 0; i < Nx; i++ )
			{
				res += Vec1[j*(Nx)+i] * Vec2[j*(Nx)+i];
			}
		}
	}
	return res;
}

double scal_l2(double * massiv1, double * massiv2, int N, int M)
{
	double out = 0.0;
	int index = 0;
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			out += massiv1[index]*massiv2[index];
			index += 1;
		}
	}
	return out;
}


double norm_max(double * massiv, int N, int M)
{
	double norm_max = 0.0;
	int index = 0;
	int i_max = 0;
	int j_max = 0;
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			if (fabs(massiv[index]) > norm_max)
			{
				norm_max = fabs(massiv[index]);
				i_max = i;
				j_max = j;
			}
			index += 1;
		}
	}
	//printf("i_max = %d \n", i_max);
	//printf("j_max = %d \n", j_max);
	return norm_max;
}

void func_p_k(double *p, int k, int n)//k = 0..n
{
	for (int i = 0; i < n + 1; i++)
		p[i] = cos(k * PI * i * 1.0/n);
}
void func_u_k(double *u, int k, int n) //k = 1..n
{
	for (int i = 1; i <= n; i++)
		u[i-1] = sin(k * PI * (2.0 * i - 1.0) * 1.0 / (2.0 * n));
}

double func_gamma_k(int k, int n)
{
	return 2.0 * sin(k * PI / (2.0 * n));
}

double dot_product(double* a, double* b, int size)
{
	double result = 0.0;
	for (int i = 0; i < size; i++)
		result += a[i] * b[i];
	return result;
}

double second_norm(double* a, int size)
{
	double result = 0.0;
	for (int i = 0; i < size; i++)
		result += a[i] * a[i];
	result = sqrt(result);
	return result;
}

double my_solution(double x,double y,double z,double t)
{
	return cos(32*PI*y);
}


double saddle_righthand_laplace300(int i, int j, int m, int n)
{
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	return hx(i)*hy(j)*exact_laplace300(x_mid,y_mid,0,0,m,n);
}

double exact_laplace105(double x, double y, double z, double time, int m, int n)
{
	return (-my*PI*my*PI)*sin(my*PI*y);
}

double saddle_righthand_laplace105(int i, int j, int m, int n)
{
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	return hx(i)*hy(j)*exact_laplace105(x_mid,y_mid,0,0,m,n);
}

double saddle_righthand_laplace(int numSolution, int i, int j, int m, int n)
{
	switch(numSolution)
	{
	case 105:
		return saddle_righthand_laplace105(i, j, m, n);
	case 300:
		return saddle_righthand_laplace300(i, j, m, n);
	default:
		break;
	}
}

int	myLambdaCalc ( double * VecX, double * VecY, double * GradFuncX, double * GradFuncY, double * tempvecX, double * tempvecY, double * nevyazkaX, double * nevyazkaY, int Nx, int Ny, double weight1, double weight2, double * lambda_pt )
{
	int dimWx = ( Nx + 1 ) * Ny;
	int dimWy = ( Ny + 1 ) * Nx;

	double lambda;

	double temp11, temp12, temp1, temp21, temp22, temp2;

	//double old_norm_sq;
	//old_norm_sq = 0.0;
	//old_norm_sq += norm_l2 ( GradFuncX, Nx + 1, Ny ) * norm_l2 ( GradFuncX, Nx + 1, Ny );
	//old_norm_sq += norm_l2 ( GradFuncY, Ny + 1, Nx ) * norm_l2 ( GradFuncY, Ny + 1, Nx );

	temp11 = 0.0;
//#ifndef CG_TEST
	if ( weight1 != 0 )
	{
		myFunctionalGrad ( VecX, VecY, T, tempvecX, tempvecY, Nx, Ny, 1, 0 );
		temp11 += scal_l2 ( tempvecX, GradFuncX, Nx + 1, Ny ) / dimWx;
		temp11 += scal_l2 ( tempvecY, GradFuncY, Ny + 1, Nx ) / dimWy;
		temp11 /= 2.0;

		//myNevyazka ( nevyazkaX, nevyazkaY, VecX, VecY, T, Nx, Ny );
		//myAvec ( tempvecX, tempvecY, nevyazkaX, nevyazkaY, Nx, Ny );
		//temp11 += scal_l2 ( tempvecX, GradFuncX, Nx + 1, Ny ) / dimWx;
		//temp11 += scal_l2 ( tempvecY, GradFuncY, Ny + 1, Nx ) / dimWy;
	}
//#endif

	temp12 = 0.0;
	if ( weight2 != 0 )
	{
		myFunctionalGrad ( VecX, VecY, T, tempvecX, tempvecY, Nx, Ny, 0, 1 );
		temp12 += scal_l2 ( tempvecX, GradFuncX, Nx + 1, Ny ) / dimWx;
		temp12 += scal_l2 ( tempvecY, GradFuncY, Ny + 1, Nx ) / dimWy;
		temp12 /= 2.0;
	}
	//temp12 += old_norm_sq;

	temp1 = weight1 * temp11 + weight2 * temp12;
	//temp1 = temp12;

	temp21 = 0.0;
	if ( weight1 != 0 )
	{
		myAvec ( tempvecX, tempvecY, GradFuncX, GradFuncY, Nx, Ny );
		temp21 += norm_l2 ( tempvecX, Nx + 1, Ny ) * norm_l2 ( tempvecX, Nx + 1, Ny );
		temp21 += norm_l2 ( tempvecY, Ny + 1, Nx ) * norm_l2 ( tempvecY, Ny + 1, Nx );
	}

	temp22 = 0.0;
//#ifndef CG_TEST
	if ( weight2 != 0 )
	{
		stability_normFull ( GradFuncX, GradFuncY, Nx, Ny, &temp22 );
		temp22 *= temp22;
	}
//#endif

	temp2 = weight1 * temp21 + weight2 * temp22;

	lambda = temp1 / temp2;

	lambda = lambda;

	//lambda = - lambda;

	*lambda_pt = lambda;
	return 0;
}


int myMinimization ( double * Wx_0, double * Wy_0, double * T_0, double ** Wx_0_new_pt, double ** Wy_0_new_pt, double converge_tol, double stagn_tol, double grad_tol, double coeff1, double coeff2, int iter_limit, int stagnation_limit, double h, int Nx, int Ny, int Nz )
{
	double weight1 = coeff1;
//#ifdef NO_H_WEIGHT2
	double weight2 = coeff2;
//#else
//	double weight2 = coeff2 * h * h * h * h;
//#endif
	double lambda;
	double omega;
	double old_value;
//#ifdef NO_H_WEIGHT2
	double new_value;
//#else
//	double new_value = 2 * param * h * h * h * h;
//#endif
	int stagnation_flag = 0;
	int stagnation_count = 0;
	int dimWx = (Nx+1) * Ny * Nz;
	int dimWy = Nx * (Ny+1) * Nz;
	int dimT = Nx * Ny * Nz;
	double rel_change = 0;

	int exit_reason = 0; // = 1, if functional is small; = 2, if stagnation; = 3, if functional increases; = 4, if gradient is small; = 5, if iteration limit is reached.
	bool converge_check = false; 
	bool stagnation_check = false; 
	bool increase_check = false; 
	bool grad_check = false; 
	bool iter_check = false;

	double old_norm_sq, new_norm_sq;
	double temp1, temp2, temp11, temp12, temp21, temp22;
	double temper;
	int iter = 0;

//	double * Wx_0_new = (double *) malloc ( dimWx * sizeof(double));
//	double * Wy_0_new = (double *) malloc ( dimWy * sizeof(double));

	double *tempvecX = (double *) malloc ( dimWx * sizeof(double));
	double *tempvecY = (double *) malloc ( dimWy * sizeof(double));

	double *nevyazkaX = (double *) malloc ( dimWx * sizeof(double));
	double *nevyazkaY = (double *) malloc ( dimWy * sizeof(double));

	double *DirectionX = (double *) malloc ( dimWx * sizeof(double));
	double *DirectionY = (double *) malloc ( dimWy * sizeof(double));
	double *newDirectionX = (double *) malloc ( dimWx * sizeof(double));
	double *newDirectionY = (double *) malloc ( dimWy * sizeof(double));

	double *GradFuncX = (double *) malloc ( dimWx * sizeof(double));
	double *GradFuncY = (double *) malloc ( dimWy * sizeof(double));

//	double *newGradFuncX = (double *) malloc ( dimWx * sizeof(double));
//	double *newGradFuncY = (double *) malloc ( dimWy * sizeof(double));


	// Wx_0, Wy_0, Wz_0 = initial vector w0 components for the minimization procedure
	double *VecX = (double *) malloc ( dimWx * sizeof(double));
	double *VecY = (double *) malloc ( dimWy * sizeof(double));

	// setting Vec = W_0
	for (int i = 0; i < dimWx; i++ )
		VecX [i] = Wx_0[i];
	for (int i = 0; i < dimWy; i++ )
		VecY [i] = Wy_0[i];

	double *newVecX = (double *) malloc ( dimWx * sizeof(double));
	double *newVecY = (double *) malloc ( dimWy * sizeof(double));

	old_value = myFunctional ( VecX, VecY, T, Nx, Ny, weight1, weight2 );
	new_value = 2 * old_value;
	rel_change = 1;

	printf ( "old_value = %e \n", old_value );
        //_getch();

	myFunctionalGrad ( VecX, VecY, T, GradFuncX, GradFuncY, Nx, Ny, weight1, weight2 );

	old_norm_sq = 0.0;
	old_norm_sq += norm_l2 ( GradFuncX, Nx + 1, Ny ) * norm_l2 ( GradFuncX, Nx + 1, Ny );
	old_norm_sq += norm_l2 ( GradFuncY, Ny + 1, Nx ) * norm_l2 ( GradFuncY, Ny + 1, Nx );

	//checking the FunctionalGrad
	printf ( "checking the myFunctionalGrad implementation, weight1 = %f should be 0, weight2= %f \n", weight1, weight2 );
	if ( weight1 == 0 )
	{
		temper = 0.0;
		temper += scal_l2 ( GradFuncX, VecX, Nx + 1, Ny ) / dimT ;
		temper += scal_l2 ( GradFuncY, VecY, Ny + 1, Nx ) / dimT ;
		temper /= 2.0;
		printf ( "old_value = %e temper = %e \n", old_value, temper );
                //_getch();
	}

	for (int i = 0; i < dimWx; i++ )
		DirectionX [i] = - GradFuncX [i];
	for (int i = 0; i < dimWy; i++ )
		DirectionY [i] = - GradFuncY [i];


	myLambdaCalc ( VecX, VecY, GradFuncX, GradFuncY, tempvecX, tempvecY, nevyazkaX, nevyazkaY, Nx, Ny, weight1, weight2, &lambda );

	printf ( "initial minimization parameters: \n" );
	printf ( "weight1 = %e weight2 = %e \n", weight1, weight2 );
	printf ( "converge_tol = %e, stagn_tol = %e, grad_tol = %e, coeff2 = %f, iter_limit = %d, stagnation_limit = %d \nold_value = %e new_value = %e \n", converge_tol, stagn_tol, grad_tol, coeff2, iter_limit, stagnation_limit, old_value, new_value );
	printf ( "old_norm_sq = %e temp2 = %e old_lambda = %e \n", old_norm_sq, temp2, lambda );
        //_getch();
	
//#ifdef NO_H_WEIGHT2
	while ( ! ( iter_check || converge_check || stagnation_check || increase_check || grad_check ) )
//#else
//	while ( new_value > param * h * h * h * h || ( fabs (new_value - old_value) < 1.0e-12 && stagnation_count == stagnation_limit ) )
//#endif
	{
		// newVec = Vec + lambda * Direction
		for (int i = 0; i < dimWx; i++ )
			newVecX [i] = VecX[i] + lambda * DirectionX[i];
		for (int i = 0; i < dimWy; i++ )
			newVecY [i] = VecY[i] + lambda * DirectionY[i];

		// calculating new_value
		new_value = myFunctional ( newVecX, newVecY, T, Nx, Ny, weight1, weight2 );
		rel_change = ( old_value - new_value ) / old_value;
		printf ( "new val = %e old_val = %e \n", new_value, old_value );
		//_getch();

		// newDirection = - GradFunctional ( newVec ) + omega * Direction

		myFunctionalGrad ( newVecX, newVecY, T, GradFuncX, GradFuncY, Nx, Ny, weight1, weight2 );

		// omega = norm^2 ( - grad Functional ( newVec ) ) / norm^2 ( - grad Functional ( Vec ) );
		new_norm_sq = 0.0;
		new_norm_sq += norm_l2 ( GradFuncX, Nx + 1, Ny ) * norm_l2 ( GradFuncX, Nx + 1, Ny );
		new_norm_sq += norm_l2 ( GradFuncY, Ny + 1, Nx ) * norm_l2 ( GradFuncY, Ny + 1, Nx );

		omega = new_norm_sq / old_norm_sq ;

		for (int i = 0; i < dimWx; i++ )
			newDirectionX [i] = - GradFuncX [i] + omega * DirectionX[i];
		for (int i = 0; i < dimWy; i++ )
			newDirectionY [i] = - GradFuncY [i] + omega * DirectionY[i];

		// checking for stagnation of iterations
		if ( rel_change > 0 && rel_change < stagn_tol )
		{
			if ( stagnation_flag == 0 )
			{
				stagnation_count = 1;
				stagnation_flag = 1;
			}
			else
			{
				printf ( "new stagnation_count = %d \n", stagnation_count );
				stagnation_count ++;
			}
		}
		else
		{
			stagnation_flag = 0;
			stagnation_count = 0;
		}

		//finishing iteration
		for (int i = 0; i < dimWx; i++ )
			VecX [i] = newVecX[i];
		for (int i = 0; i < dimWy; i++ )
			VecY [i] = newVecY[i];

		for (int i = 0; i < dimWx; i++ ) 
			DirectionX [i] = newDirectionX [i];
		for (int i = 0; i < dimWy; i++ )
			DirectionY [i] = newDirectionY [i];

		old_norm_sq = new_norm_sq;
		old_value = new_value;

		// calculating lambda = ...
		myLambdaCalc ( VecX, VecY, GradFuncX, GradFuncY, tempvecX, tempvecY, nevyazkaX, nevyazkaY, Nx, Ny, weight1, weight2, &lambda );
	
		iter++;
		printf ( "iteration %d: finished, new value = %e \n", iter, new_value );
		printf ( "iteration %d: rel_change = %e, stagn_count = %d, iter_gradtol = %e \n", iter, rel_change, stagnation_count, old_norm_sq );
		
		//old
		//converge_check = ( new_value < converge_tol * converge_tol * h * h * h * h );
		//new
		converge_check = ( new_value < converge_tol * converge_tol );
		stagnation_check = ( rel_change > 0 && rel_change < stagn_tol && stagnation_count == stagnation_limit );
		increase_check = ( rel_change < 0 ) ;
		grad_check = ( old_norm_sq < grad_tol * grad_tol );
		iter_check = ( iter + 1 > iter_limit );

		if ( converge_check == true )
			exit_reason = 1;
		if ( stagnation_check == true )
			exit_reason = 2;
		if ( increase_check == true )
			exit_reason = 3;
		if ( grad_check == true )
			exit_reason = 4;
		if ( iter_check == true )
			exit_reason = 5;
		//_getch();
	}

	printf ( "\nexit_reason: %d \n", exit_reason );
	switch ( exit_reason )
	{
	case 0:
		printf ( "no iterations performed \n" );
		break;
	case 1:
		printf ( "iterations converged (functional is small) \n" );
		break;
	case 2:
		printf ( "stagnation occured \n" );
		break;
	case 3:
		printf ( "functional suddenly increases \n" );
		break;
	case 4:
		printf ( "gradient is small enough \n" );
		break;
	case 5:
		printf ( "iterations limit is reached \n" );
		break;
	default:
		printf ( "bad value of exit_reason = %d (should be 0, 1, 2, 3, 4 or 5) \n", exit_reason );
		return -1;
		break;
	}
	printf ( "final value = %e \n", new_value );

	//setting output W_0_new vector's components to newVecX, newVecY, newVecZ
	*Wx_0_new_pt = VecX;
	*Wy_0_new_pt = VecY;

	//for ( int i = 0; i < Nx; i++ )
	//	printf ( "Wx_0_new [%d] = %f \n", i, VecX[i] );
	//for ( int i = 0; i < Ny; i++ )
	//	printf ( "Wy_0_new [%d] = %f \n", i, VecWy_0_new[i] );

	free ( tempvecX );
	free ( tempvecY );

	free ( nevyazkaX );
	free ( nevyazkaY );

	free ( newVecX );
	free ( newVecY );

	free ( GradFuncX );
	free ( GradFuncY );

	free ( DirectionX );
	free ( DirectionY );
	free ( newDirectionX );
	free ( newDirectionY );

	return 0;
}


int myFunctionalGrad ( double * VecX, double * VecY, double * T, double * GradFuncX, double * GradFuncY, int Nx, int Ny, double weight1, double weight2 )
{
	int dimWx = ( Nx + 1 ) * Ny;
	int dimWy = Nx * ( Ny + 1 );
	int dimT = Nx * Ny;
	double temp1, temp2;

	double * nevyazkaX, * nevyazkaY, * tempvec1X, * tempvec1Y;
	double * MhalfLambdaY_Btr_w0, * Mminus1Btr_Aminus1_xvost_w0, * tempvec2X, * tempvec2Y;

	// 1) calculating A* * (A * vec - f ) = A * ( A * vec - f )
	if ( weight1 != 0 )
	{
		// 1.1) calculating A * vec - f = nevyazka
		nevyazkaX = ( double * ) malloc ( dimWx * sizeof ( double ) );
		nevyazkaY = ( double * ) malloc ( dimWy * sizeof ( double ) );

		myNevyazka ( nevyazkaX, nevyazkaY, VecX, VecY, T, Nx, Ny );

		tempvec1X = ( double * ) malloc ( dimWx * sizeof ( double ) );
		tempvec1Y = ( double * ) malloc ( dimWy * sizeof ( double ) );

		// 1.2) calculating A * ( A * vec - f ) = A * nevyazka
		myAvec ( tempvec1X, tempvec1Y, nevyazkaX, nevyazkaY, Nx, Ny );
	}
	// 2) calculating Lambda * Lamda_y * Btr * vec
//#ifndef CG_TEST
	if ( weight2 != 0 )
	{
		// 2.1) calculating Lambda_y * Btr * vec
		MhalfLambdaY_Btr_w0 = (double*)malloc(dimT * sizeof(double));
		myMhalfLambdayBtr_w0 ( MhalfLambdaY_Btr_w0, VecX, VecY, Nx, Ny );
//			temp1 = norm_l2(MhalfLambdaY_Btr_w0, Ny, Nx);
//			printf("norm_l2 for MhalfLambdaY_Btr_w0 = %f \n", temp1);

		// 2.2) calculating Mminushalf * Lambda * Lambda_y * Btr * vec based on results of 1)
		Mminus1Btr_Aminus1_xvost_w0 = (double*)malloc(dimT * sizeof(double));  // = Lambda * Lambda_y * Btr * w0
		myMminushalfLambdaLambdayBtr_w0 ( Mminus1Btr_Aminus1_xvost_w0, MhalfLambdaY_Btr_w0, VecX, VecY, Nx, Ny );
//			temp1 = norm_l2(Mminus1Btr_Aminus1_xvost_w0, Ny, Nx);
//			printf("norm_l2 for Mminus1Btr_Aminus1_xvost_w0 = %f \n", temp1);

//			temp1 = scal_l2(MhalfLambdaY_Btr_w0, Mminus1Btr_Aminus1_xvost_w0, Ny, Nx) / dimT ;
//			printf("right value = %f \n", temp1);

		// 2.3) calculating B Lambday * Lambda * Lambday * Btr * vec based on results of 2)
		tempvec2X = ( double * ) malloc ( dimWx * sizeof ( double ) );
		tempvec2Y = ( double * ) malloc ( dimWy * sizeof ( double ) );
		myBLambdayLambdaLambdayBtr_w0 ( tempvec2X, tempvec2Y, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny );
		// 2.4) checking calculating B Lambday * Lambda * Lambday * Btr * vec based on results of 2)
//			temp1 = norm_l2(tempvec2X, Nx + 1, Ny);
//			printf("norm_l2 for tempvec2X = %f \n", temp1);
//			temp1 = norm_l2(tempvec2Y, Ny + 1, Nx);
//			printf("norm_l2 for tempvec2Y = %f \n", temp1);
//			temp1 = scal_l2(tempvec2X, tempvec2X, Nx + 1, Ny) + scal_l2(tempvec2Y, tempvec2Y, Nx, Ny + 1 ) / dimWy;
//			temp1 = sqrt ( temp1 );
//			printf("norm_l2 for tempvec2 = %f \n", temp1);

		//	double * Y1 = ( double * ) malloc ( dimWy * sizeof ( double ) );
		//myByMminus1LambdaLambdayBtr_w0 ( Y1, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny );
		//	double * Y2 = ( double * ) malloc ( dimWy * sizeof ( double ) );
		//myAyminus1ByMminus1Btr_w0 ( Y2, VecX, VecY, Nx, Ny );
		//	temp1 = scal_l2( Y1, Y2, Ny + 1, Nx) / dimT;
		//	//temp1 = sqrt ( temp1 );
		//	printf( "(Y1, Y2) = %f \n", temp1);

		//	double * Y3 = ( double * ) malloc ( dimWy * sizeof ( double ) );
		//myAyminus1ByMminus1LambdaLambdayBtr_w0 ( Y3, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny );
		//	double * Y4 = ( double * ) malloc ( dimWy * sizeof ( double ) );
		//myByMminus1Btr_w0 ( Y4, VecX, VecY, Nx, Ny );
		//	temp1 = scal_l2( Y3, Y4, Ny + 1, Nx) / dimT;
		//	//temp1 = sqrt ( temp1 );
		//	printf( "(Y3, Y4) = %f \n", temp1);

		//	double * Y5 = ( double * ) malloc ( dimT * sizeof ( double ) );
		//myMminus1BytrAyminus1ByMminus1LambdaLambdayBtr_w0 ( Y5, Mminus1Btr_Aminus1_xvost_w0, Nx, Ny );
		//	double * Y6 = ( double * ) malloc ( dimT * sizeof ( double ) );
		//myBtr_w0 ( Y6, VecX, VecY, Nx, Ny );
		//	temp1 = scal_l2( Y5, Y6, Nx, Ny) / dimT;
		//	//temp1 = sqrt ( temp1 );
		//	printf( "(Y5, Y6) = %f \n", temp1);

	}

	// 3) calculating GradFunc = tempvec1 + weight * tempvec2
	for ( int i = 0; i < dimWx; i++ )
		if ( weight1 != 0 && weight2 != 0 )
			//GradFuncX [i] = weight1 * tempvec1X [i] + weight2 * tempvec2X [i];
			GradFuncX [i] = weight1 * 2.0 * tempvec1X [i] + weight2 * 2.0 * tempvec2X [i];
		else
		{
			if ( weight1 != 0 )
				//GradFuncX [i] = weight1 * tempvec1X [i];
				GradFuncX [i] = weight1 * 2.0 * tempvec1X [i];
			else //if ( weight2 != 0 )
				//GradFuncX [i] = weight2 * tempvec2X [i];
				GradFuncX [i] = weight2 * 2.0 * tempvec2X [i];
		}
	for ( int i = 0; i < dimWy; i++ )
		if ( weight1 != 0 && weight2 != 0 )
			//GradFuncY [i] = weight1 * tempvec1Y [i] + weight2 *  tempvec2Y [i];
			GradFuncY [i] = weight1 * 2.0 * tempvec1Y [i] + weight2 * 2.0 * tempvec2Y [i];
		else
		{
			if ( weight1 != 0 )
				//GradFuncY [i] = weight1 * tempvec1Y [i];
				GradFuncY [i] = weight1 * 2.0 * tempvec1Y [i];
			else
				//GradFuncY [i] = weight2 * tempvec2Y [i];
				GradFuncY [i] = weight2 * 2.0 * tempvec2Y [i];
		}

//#else
	//// 3) calculating GradFunc = tempvec1 
	//for ( int i = 0; i < dimWx; i++ )
	//	GradFuncX [i] = tempvec1X [i] ;
	//for ( int i = 0; i < dimWy; i++ )
	//	GradFuncY [i] = tempvec1Y [i] ;
//#endif

	if ( weight1 != 0 )
	{
		free (nevyazkaX);
		free (nevyazkaY);
		free (tempvec1X);
		free (tempvec1Y);
	}
	if ( weight2 != 0 )
	{
		//printf ( "checking myFunctionalGrad calculation \n" );
		//double temp1 = 0.0;
		//temp1 += scal_l2 ( tempvec2X, VecX, Nx + 1, Ny ) / dimT;
		//temp1 += scal_l2 ( tempvec2Y, VecY, Ny + 1, Nx ) / dimT;
		//double temp2 = myFunctional ( VecX, VecY, T, Nx, Ny, weight1, weight2 );
		//printf ( "temp1 = %e temp2 = %e must be equal temp1/temp2 = %f \n", temp1, temp2, temp1 / temp2 );

		free (tempvec2X);
		free (tempvec2Y);
		free ( MhalfLambdaY_Btr_w0 );
		free ( Mminus1Btr_Aminus1_xvost_w0 );
	}

	return 0;
}

double myFunctional ( double * VecX, double * VecY, double * T, int Nx, int Ny, double weight1, double weight2 )
{
	double out = 0.0;
	double temp1 = 0.0, temp2 = 0.0;

//#ifndef CG_TEST
	// temp1 = stability_normFull ^2
	if ( weight2 != 0 )
	{
		stability_normFull( VecX, VecY, Nx, Ny, &temp2 );
		//printf ( "stab_full_norm = %f \n", temp2 );
		temp2 *= temp2;
	}
//#endif

	// temp2 = norm (Aq - f) ^2
	if ( weight1 != 0 )
	{
		normAQminusF ( VecX, VecY, T, Nx, Ny, &temp1 );
		//printf ( "nevyazka_norm = %e wegith1 = %e weight2 = %e \n", temp1, weight1, weight2 );
		temp1 *= temp1;
	}

	printf ( "weight1 = %e temp1 = %e weight2 =%e temp2 = %e \n", weight1, temp1, weight2, temp2 );
	out = weight1 * temp1 + weight2 * temp2;

	return out;
}


int	myAvec ( double *tempvecX, double *tempvecY, double *VecX, double *VecY, int Nx, int Ny )
{

	double a_0_old, a_i_old, a_iminus1_old, a_nminus1_old;
	double b_0_old, b_i_old, b_n_old = 0;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	switch(xCondition)
	{
	case eNeumann:
		//Нейман, x
		for (int i = 0; i < Ny*Nz; i++)
		{
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			tempvecX [ i * ( Nx  + 1 ) + 0 ] = 0.0;

			for ( int k = 1 ; k < Nx ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)/(heat_conductivity(k,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)/(heat_conductivity(k-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==Nx-1)
				{
					a_i_old = 0.0;
				}
				
				if (k==1)
				{
					a_iminus1_old = 0.0;
				}

				tempvecX [ i * (Nx + 1) + k ] = a_iminus1_old * VecX [ i * (Nx + 1 ) + k - 1 ] + b_i_old * VecX [ i * (Nx + 1 ) + k ] + a_i_old * VecX [ i * (Nx + 1 ) + k + 1 ];
			}

			tempvecX [ i * ( Nx  + 1 ) + Nx ] = 0.0;

		}
		break;
	case eDirichlet:
		//для неоднородных условий Дирихле
		//Дирихле, x
		{;}

		for (int i = 0; i < Ny*Nz; i++)
		{
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			//к-ты матрицы А
			a_0_old = hx(0)*hy(j1)/(heat_conductivity(0,j1,k1)*6.0);
			b_0_old = 2*a_0_old;

			tempvecX [ i * (Nx + 1) + 0 ] = b_0_old * VecX [ i * (Nx + 1) + 0 ] + a_0_old * VecX [ i * (Nx + 1) + 1 ];

			for ( int k = 1 ; k < Nx ; k++ )
			{
				i1 = k;
				//к-ты матрицы А
				a_i_old = hx(k)*hy(j1)/(heat_conductivity(i1,j1,k1)*6.0);
				a_iminus1_old = hx(k-1)*hy(j1)/(heat_conductivity(i1-1,j1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				tempvecX [ i * (Nx + 1) + k ] = a_iminus1_old * VecX [ i * (Nx + 1 ) + k - 1 ] + b_i_old * VecX [ i * (Nx + 1 ) + k ] + a_i_old * VecX [ i * (Nx + 1 ) + k + 1 ];
			}
			a_nminus1_old = hx(Nx-1)*hy(j1)/(heat_conductivity(Nx-1,j1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			tempvecX [ i * (Nx + 1) + Nx ] = a_nminus1_old * VecX [ i * (Nx + 1) + Nx - 1 ] + b_n_old * VecX [ i * (Nx + 1) + Nx ] ;

		}
		break;
	}

	switch(yCondition)
	{
	case eNeumann:
		//Нейман, y
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;

			tempvecY [ i * (Ny + 1) + 0 ] = 0.0;

			for ( int k = 1 ; k < Ny ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				if (k==Ny-1)
				{
					a_i_old = 0.0;
				}
				
				if (k==1)
				{
					a_iminus1_old = 0.0;
				}
				tempvecY [ i * (Ny + 1) + k ] = a_iminus1_old * VecY [ i * (Ny + 1 ) + k - 1 ] + b_i_old * VecY [ i * (Ny + 1 ) + k ] + a_i_old * VecY [ i * (Ny + 1 ) + k + 1 ];
			}

			tempvecY [ i * (Ny + 1) + Ny ] = 0.0;
		}
		break;
	case eDirichlet:
		//Дирихле, y
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;


			//к-ты матрицы А
			a_0_old = hy(0)*hx(i1)/(heat_conductivity(i1,0,k1)*6.0);
			b_0_old = 2*a_0_old;

			tempvecY [ i * (Ny + 1) + 0 ] = b_0_old * VecY [ i * (Ny + 1 ) + 0 ] + a_0_old * VecY [ i * (Ny + 1 ) + 1 ];

			for ( int k = 1 ; k < Ny ; k++ )
			{
				//к-ты матрицы А
				a_i_old = hy(k)*hx(i1)/(heat_conductivity(i1,k,k1)*6.0);
				a_iminus1_old = hy(k-1)*hx(i1)/(heat_conductivity(i1,k-1,k1)*6.0);
				b_i_old = 2*a_iminus1_old + 2*a_i_old;

				tempvecY [ i * (Ny + 1) + k ] = a_iminus1_old * VecY [ i * (Ny + 1 ) + k - 1 ] + b_i_old * VecY [ i * (Ny + 1 ) + k ] + a_i_old * VecY [ i * (Ny + 1 ) + k + 1 ];
			}
			a_nminus1_old = hy(Ny-1)*hx(i1)/(heat_conductivity(i1,Ny-1,k1)*6.0);
			b_n_old = 2*a_nminus1_old;

			tempvecY [ i * (Ny + 1) + Ny ] = a_nminus1_old * VecY [ i * (Ny + 1 ) + Ny - 1 ] + b_n_old * VecY [ i * (Ny + 1 ) + Ny ];

		}
		break;
	}

	return 0;
}

int myNevyazka ( double *nevyazkaX, double * nevyazkaY, double * VecX, double * VecY, double * T, int Nx, int Ny )
{

	myAvec ( nevyazkaX, nevyazkaY, VecX, VecY, Nx, Ny ); // A * Vec = nevyazka

	double temper_0, temper_k, temper_n = 0;
	double Tx_0_temp, Tx_1_temp, Ty_0_temp, Ty_1_temp, wx_0_temp, wx_1_temp, wy_0_temp, wy_1_temp;

	int i1 = 0; 
	int j1 = 0;
	int k1 = 0;

	switch(xCondition)
	{
	case eNeumann:
		//Нейман, x
		for (int i = 0; i < Ny*Nz; i++)
		{
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			for ( int k = 1 ; k < Nx ; k++ )
			{
				temper_k = hy(j1) * ( T[k - 1 + j1*Nx + k1*Nx*Ny] - T[k + j1*Nx + k1*Nx*Ny] );

				nevyazkaX [ i * (Nx + 1) + k ] -= temper_k;
			}
		}
		break;
	case eDirichlet:
		//Дирихле, x

		for (int i = 0; i < Ny*Nz; i++)
		{
			j1 = i - (i/Ny)*Ny;
			k1 = i/Ny;

			Tx_0_temp = Tx_0_bound(numSolution,0,j1,0*tau);
			Tx_1_temp = Tx_1_bound(numSolution,Nx-1,j1,0*tau);

			temper_0 = hy(j1) * ( Tx_0_temp - T[0 + j1*Nx + k1*Nx*Ny] );
			
			nevyazkaX [ i * (Nx + 1) + 0 ] -= temper_0;

			for ( int k = 1 ; k < Nx ; k++ )
			{
				i1 = k;

				//правая часть
				temper_k = hy(j1) * ( T[k - 1 + j1*Nx + k1*Nx*Ny] - T[k + j1*Nx + k1*Nx*Ny] );

				nevyazkaX [ i * (Nx + 1) + k ] -= temper_k;
			}
			temper_n = hy (j1) * ( T[Nx-1 + j1*Nx + k1*Nx*Ny] - Tx_1_temp );

			nevyazkaX [ i * (Nx + 1) + Nx ] -= temper_n;
		}
		break;
	}

	switch(yCondition)
	{
	case eNeumann:
		//Нейман, y
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;

			for ( int k = 1 ; k < Ny ; k++ )
			{
				temper_k = hx(i1)*(T[i1 + (k - 1)*Nx + k1*Nx*Ny] - T[i1 + k*Nx + k1*Nx*Ny]);

				nevyazkaY [ i * (Ny + 1) + k ] -= temper_k;
			}
		}
		break;
	case eDirichlet:
		//Дирихле, y
		for ( int i = 0 ; i < Nx*Nz ; i++ )
		{
			i1 = i/Nz;
			k1 = i - (i/Nz)*Nz;

			Ty_0_temp = Ty_0_bound(numSolution,i1,0,0*tau);
			Ty_1_temp = Ty_1_bound(numSolution,i1,Ny-1,0*tau);

			temper_0 = hx(i1) * ( Ty_0_temp - T[i1 + 0*Nx + k1*Nx*Ny] );

			nevyazkaY [ i * (Ny + 1) + 0 ] -= temper_0;

			for ( int k = 1 ; k < Ny ; k++ )
			{
				temper_k = hx(i1) * ( T[i1 + (k - 1)*Nx + k1*Nx*Ny] - T[i1 + k*Nx + k1*Nx*Ny] );

				nevyazkaY [ i * (Ny + 1) + k ] -= temper_k;
			}
			temper_n = hx(i1) * ( T[i1 + (Ny-1)*Nx + k1*Nx*Ny] - Ty_1_temp );

			nevyazkaY [ i * (Ny + 1) + Ny ] -= temper_n;
		}
		break;
	}

	return 0;
}
double normAQminusF ( double * VecX, double * VecY, double * T, int Nx, int Ny, double * l2_norm_pt )
{
	double out = 0.0;
	int dimWx = (Nx + 1 ) * Ny;
	int dimWy = Nx * ( Ny + 1 );
	int Nz = 1;

	double * nevyazkaX = ( double * ) malloc ( dimWx * sizeof ( double ) );
	double * nevyazkaY = ( double * ) malloc ( dimWy * sizeof ( double ) );

	myNevyazka ( nevyazkaX, nevyazkaY, VecX, VecY, T, Nx, Ny );

	double temp = 0.0;
	temp += norm_l2 (nevyazkaX, Nx + 1, Ny ) * norm_l2 (nevyazkaX, Nx + 1, Ny );
	temp += norm_l2 (nevyazkaY, Ny + 1, Nx ) * norm_l2 (nevyazkaY, Ny + 1, Nx );

	////NEW
	//temp *= Nx * Ny;

	out = sqrt(temp);

	free( nevyazkaX );
	free( nevyazkaY );

	*l2_norm_pt = out;
	return out;
}



int stability_norm2D ( double * Wx, double * Wy, int Nx, int Ny, double tau, BNDCNDS xCondition, BNDCNDS yCondition, double * res_pt)
{
	double res;
	double term1, term2;

	int dimT = Nx * Ny;

	double * tempT = (double*) malloc (dimT * sizeof(double));
	double * tempT2 = (double*) malloc (dimT * sizeof(double));

	// 1. compute term1 = || W ||_H = || W ||_Hdivh
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Nx, Ny);							// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0);

	term1 = norm_l2 ( tempT, Nx, Ny ) * norm_l2 ( tempT, Nx, Ny );								// term1 = || tempT ||_L2h^2 
	term1 += norm_l2 ( Wx, Nx + 1, Ny) * norm_l2 ( Wx, Nx + 1, Ny) +							
		norm_l2 ( Wy, Nx, Ny+1 ) * norm_l2 ( Wy, Nx, Ny+1 );									// term1 = || tempT ||_L2h^2 + || W ||_L2hv^2
	term1 = sqrt (term1);																		// term1 = || W ||_Hdivh

	printf ( "term1 = %f \n", term1 );
	printf ( "norm_BtrW_L2 = %f \n", norm_l2 ( tempT, Nx, Ny ) );
	printf ( "norm_W_L2 = %f \n", sqrt(norm_l2 ( Wx, Nx + 1, Ny) * norm_l2 ( Wx, Nx + 1, Ny) + norm_l2 ( Wy, Nx, Ny+1 ) * norm_l2 ( Wy, Nx, Ny+1 ) ) );

	// 2. compute term2 = || C Btr W ||_Lambda
	zero_init (tempT, dimT);
	Btr_W_Add( tempT, 1.0, Wx, xCondition, 1.0, Wy, yCondition, 1.0, Nx, Ny);					// tempT = Btr W
	for ( int i = 0; i < dimT; i++ )
		tempT[i] *= 1.0 / hx(0) / hy(0);
	printf ( "norm_BtrW = %f \n", norm_l2 (tempT, Nx, Ny) );

	compute_Lambday2D ( tempT2, tempT, Nx, Ny, xCondition, yCondition );						// tempT2 = Lambday tempT = Lambday Btr W
	printf ( "norm_Lambday_BtrW = %f \n", norm_l2 (tempT2, Nx, Ny) );
	zero_init (tempT, dimT);
	compute_Lambda2D (tempT, tempT2, Nx, Ny, xCondition, yCondition );							// tempT (new) = Lambda * tempT2 = Lambda * Lambday * Btr * W
	printf ( "norm_Lambda_CBtrW = %f \n", norm_l2 (tempT, Nx, Ny) );
	term2 = sqrt ( scal_l2 (tempT, tempT2, Nx, Ny ) / dimT );

	printf ( "term2 = %f \n", term2 );

	res = sqrt ( term1 * term1 + pow(tau,4) * term2 * term2 );

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

// computes out = Lambdax * in
int	compute_Lambdax2D( double * out, double * in, int Nx, int Ny, BNDCNDS xCondition, BNDCNDS yCondition)
{
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	// 1. tempWx = Bx * in
	double * tempWx = (double*)malloc(dimWx * sizeof(double));
	zero_init (tempWx, dimWx);

	BxMminus1_F_Add( tempWx, in, xCondition, Nx, Ny, 1);					// tempWx1 = Bx * in

	// 2. tempWx2 = Ax^(-1) * Bx * in
	double * tempWx2 = (double*)malloc(dimWx * sizeof(double));		
	zero_init (tempWx2, dimWx);

	Invert_Ax( tempWx2, tempWx, xCondition, Nx, Ny, 1);						// tempWx2 = Ax^(-1) * tempWx1 = Ax^(-1) * Bx * in
	for ( int i = 0; i < dimWx; i++ )
		tempWx2[i] *= hx(0);

	// 3. out = Bxtr * Ax^(-1) * Bx * in
	zero_init (out, dimT);
	Bxtr_Wx_Add( out, 1.0, tempWx2, xCondition, Nx, Ny, 1);
	for ( int i = 0; i < dimT; i++ )
		out[i] *= 1.0 / hx(0) / hy(0);

	free (tempWx);
	free (tempWx2);

	return 0;
}



// computes out = Lambday * in
int	compute_Lambday2D ( double * out, double * in, int Nx, int Ny, BNDCNDS xCondition, BNDCNDS yCondition)
{
	int dimWx = (Nx + 1) * Ny;
	int dimWy = Nx * (Ny + 1);
	int dimT = Nx * Ny;

	// 1. tempWy = By * in
	double * tempWy = (double*)malloc(dimWy * sizeof(double));
	zero_init (tempWy, dimWy);

	ByMminus1_F_Add( tempWy, in, yCondition, Nx, Ny, 1);					// tempWy1 = By * in

	// 2. tempWy2 = Ay^(-1) * By * in
	double * tempWy2 = (double*)malloc(dimWy * sizeof(double));		
	zero_init (tempWy2, dimWy);

	Invert_Ay( tempWy2, tempWy, yCondition, Nx, Ny, 1);						// tempWy2 = Ay^(-1) * tempWy1 = Ay^(-1) * By * in
	for ( int i = 0; i < dimWy; i++ )
		tempWy2[i] *= hy(0);

	// 3. out = Bytr * Ay^(-1) * By * in
	zero_init (out, dimT);
	Bytr_Wy_Add( out, 1.0, tempWy2, yCondition, Nx, Ny, 1);
	for ( int i = 0; i < dimT; i++ )
		out[i] *= 1.0 / hx(0) / hy(0);

	free (tempWy);
	free (tempWy2);

	return 0;
}


int	compute_Lambda2D ( double * out, double * in, int Nx, int Ny, BNDCNDS xCondition, BNDCNDS yCondition)
{
	int dimT = Nx * Ny;

	double * tempT1 = (double*)malloc(dimT * sizeof(double));
	double * tempT2 = (double*)malloc(dimT * sizeof(double));

	compute_Lambdax2D(tempT1, in, Nx, Ny, xCondition, yCondition);
	compute_Lambday2D(tempT2, in, Nx, Ny, xCondition, yCondition);

	for ( int i = 0; i < dimT; i++ )
		out[i] = tempT1[i] + tempT2[i];

	free (tempT1);
	free (tempT2);

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
				temper_k = rhand[i*(My+1) + 0];

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
