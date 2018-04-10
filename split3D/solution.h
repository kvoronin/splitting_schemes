#include "math.h"

typedef enum{
  eDirichlet = 0,
  eNeumann,
  eMixed
} boundaryConditions;


//! the initial data
boundaryConditions xCondition, yCondition, zCondition;

extern double tau;
extern double *xpoints, *ypoints, *zpoints;

const int velocityX_num = 0;//номер скорости Vx для exact_solution1111
const int velocityY_num = 0;//номер скорости Vy для exact_solution1111
const int velocityZ_num = 0;//номер скорости Vz для exact_solution1111

//const int num_sol = 206;//номер тестового решения
const int mx = 2; //номер гармоники по x для тестового решения 
//const int n = 2; //номер гармоники по y для тестового решения
const int my = 2; //номер гармоники по y для тестового решения
const int mz = 2;//номер гармоники по z для тестового решения
const int Dir_axis = 2;

#define MT500 (1.0)
//also used in 889
#define DEG888 (3.1)
#define A888 (0.1)
//#define A888 (0.1875)
#define C888 (64.0)
//#define C888 (1.0)
// must be an integer value
#define DEG444 (4.0)

//#ifdef SPECIALTEST
//	#define NUMSOLSPECIAL (889)
//#endif


const int pause_var = 0;     //делать ли паузу, не работает на кластере, поэтому в нужном месте закомментировано.
const int printcase = 1;     //наличие выдачи в файл промежуточных выкладок
const int print_step = 20; //шаг выдачи промежуточных значений погрешности на экран и в файлы.
const int nonstationary = 1; //нестационарность теста
//const int special_print = 10000;
//double tau = 0.000001;	//временной шаг
double Time = 1;		//промежуток времени [0,Time]
const double PI = 3.14159265358979323;
double Tx_0 = 0;		
double Tx_1 = 0;		
double Ty_0 = 0;		
double Ty_1 = 0;		
double Tz_0 = 0;		
double Tz_1 = 0;		
double wx_0 = 0;		
double wx_1 = 0;		
double wy_0 = 0;		
double wy_1 = 0;		
double wz_0 = 0;		
double wz_1 = 0;

double vx_const = 1.0;  //значение для постоянной скорости Vx
double vy_const = 1.0;  //значение для постоянной скорости Vy
double vz_const = 1.0;  //значение для постоянной скорости Vz


//char* Bnd_x = "D";
//char* Bnd_y = "N";
//char* Bnd_z = "N";
//int bnd_X, bnd_Y, bnd_Z;

void setBoundaryConditions(int number);
double exact_solution(int number, double x, double y, double z, double t); 
double exact_gradientX(int number, double x, double y, double z, double t);
double exact_gradientY(int number, double x, double y, double z, double t); 
double exact_gradientZ(int number, double x, double y, double z, double t); 
double righthand(int number, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double exact_solution_averageQuad ( int number, double time, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints,
									   int NQuadpoints );

double exact_solution1(double x, double y, double z, double t); //первый тест
double exact_solution2(double x, double y, double z, double t); //второй тест
double exact_solution3(double x, double y, double z, double t); //третий тест
double exact_solution4(double x, double y, double z, double t); //четвертый тест
double exact_solution5(double x, double y, double z, double t); 

double exact_solution6(double x, double y, double z, double t, int mx); //неоднородные условия Неймана, x
double exact_gradient6(double x, double y, double z, double t, int mx); 
double righthand6(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution2006(double x, double y, double z, double t, int mz); //неоднородные условия Неймана, x
double exact_gradient2006(double x, double y, double z, double t, int mz); 
double righthand2006(int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution2066(double x, double y, double z, double t, int mz); //неоднородные условия Неймана, x
double exact_gradient2066(double x, double y, double z, double t, int mz); 
double righthand2066(int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution366(double x, double y, double z, double t, int mx); //неоднородные условия Дирихле
double exact_gradient366x(double x, double y, double z, double t, int mx); 
double righthand366(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);


double exact_solution66(double x, double y, double z, double t, int my); //неоднородные условия Дирихле
double exact_gradient66y(double x, double y, double z, double t, int my); 
double righthand66(int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution106(double x, double y, double z, double t, int mx); //неоднородные условия Неймана, y
double exact_gradient106(double x, double y, double z, double t, int mx); 
double righthand106(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution206(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient206x(double x, double y, double z, double t, int mx, int my); 
double exact_gradient206y(double x, double y, double z, double t, int mx, int my); 
double righthand206(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution207(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient207x(double x, double y, double z, double t, int mx, int my); 
double exact_gradient207y(double x, double y, double z, double t, int mx, int my); 
double righthand207(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution266(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient266x(double x, double y, double z, double t, int mx, int my); 
double exact_gradient266y(double x, double y, double z, double t, int mx, int my); 
double righthand266(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1566(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient1566x(double x, double y, double z, double t, int mx, int my); 
double exact_gradient1566y(double x, double y, double z, double t, int mx, int my); 
double righthand1566(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution267(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient267x(double x, double y, double z, double t, int mx, int my); 
double exact_gradient267y(double x, double y, double z, double t, int mx, int my); 
double righthand267(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution66(double x, double y, double z, double t, int mx); //неоднородные условия Дирихле
double exact_gradient66(double x, double y, double z, double t, int mx); 
double righthand66(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution10(double x, double y, double z, double t, int mx); 
double exact_gradient10(double x, double y, double z, double t, int mx); 
double righthand10(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution11(double x, double y, double z, double t, int mx); 
double exact_gradient11(double x, double y, double z, double t, int mx); 
double righthand11(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution12(double x, double y, double z, double t, int mx); 
double exact_gradient12(double x, double y, double z, double t, int mx); 
double righthand12(int mx, int t_discr, int k, double* xpoints);

double exact_solution1266(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient1266z(double x, double y, double z, double t, int mx, int my); 
double exact_gradient1266y(double x, double y, double z, double t, int mx, int my); 
double righthand1266(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1466(double x, double y, double z, double t, int mz, int my); //неоднородные условия Неймана y, xy
double exact_gradient1466z(double x, double y, double z, double t, int mz, int my); 
double exact_gradient1466y(double x, double y, double z, double t, int mz, int my); 
double righthand1466(int mz, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);


double exact_solution1206(double x, double y, double z, double t, int mz, int my); //неоднородные условия Неймана y, xy
double exact_gradient1206z(double x, double y, double z, double t, int mz, int my); 
double exact_gradient1206y(double x, double y, double z, double t, int mz, int my); 
double righthand1206(int mz, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);


double exact_solution2206(double x, double y, double z, double t, int mz, int mx); //неоднородные условия Неймана y, xy
double exact_gradient2206z(double x, double y, double z, double t, int mz, int mx); 
double exact_gradient2206x(double x, double y, double z, double t, int mz, int mx); 
double righthand2206(int mz, int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1366(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient1366z(double x, double y, double z, double t, int mx, int my); 
double exact_gradient1366x(double x, double y, double z, double t, int mx, int my); 
double righthand1366(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution2366(double x, double y, double z, double t, int mx, int my); //неоднородные условия Неймана y, xy
double exact_gradient2366z(double x, double y, double z, double t, int mx, int my); 
double exact_gradient2366x(double x, double y, double z, double t, int mx, int my); 
double righthand2366(int mx, int my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution9006(double x, double y, double z, double t, int mx, int my, int mz); //неоднородные условия Неймана y, xy
double exact_gradient9006x(double x, double y, double z, double t, int mx, int my, int mz); 
double exact_gradient9006y(double x, double y, double z, double t, int mx, int my, int mz); 
double exact_gradient9006z(double x, double y, double z, double t, int mx, int my, int mz); 
double righthand9006(int mx, int my, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution9066(double x, double y, double z, double t, int mx, int my, int mz); //неоднородные условия Неймана y, xy
double exact_gradient9066x(double x, double y, double z, double t, int mx, int my, int mz); 
double exact_gradient9066y(double x, double y, double z, double t, int mx, int my, int mz); 
double exact_gradient9066z(double x, double y, double z, double t, int mx, int my, int mz); 
double righthand9066(int mx, int my, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1111(double x, double y, double z, double t, int mx, int my, int mz); //неоднородные условия Неймана y, xy
double exact_solution1111Vx(double x, double y, double z, double t, int mx, int my, int mz); //неоднородные условия Неймана y, xy
double exact_solution1111Vy(double x, double y, double z, double t, int mx, int my, int mz); //неоднородные условия Неймана y, xy
double exact_solution1111Vz(double x, double y, double z, double t, int mx, int my, int mz); //неоднородные условия Неймана y, xy
double exact_gradient1111x(double x, double y, double z, double t, int mx, int my, int mz); 
double exact_gradient1111y(double x, double y, double z, double t, int mx, int my, int mz); 
double exact_gradient1111z(double x, double y, double z, double t, int mx, int my, int mz); 
double righthand1111(int mx, int my, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

/*
double exact_solution1166(double x, double y, double z, double t, int mx, int my, int mz, int dir_axis); //неоднородные условия Неймана y, xy
double exact_gradient1166x(double x, double y, double z, double t, int mx, int my, int mz, int dir_axis); 
double exact_gradient1166y(double x, double y, double z, double t, int mx, int my, int mz, int dir_axis); 
double exact_gradient1166z(double x, double y, double z, double t, int mx, int my, int mz, int dir_axis); 
double righthand1166(int dir_axis, int mx, int my, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
*/
double exact_solution1166(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); //неоднородные условия Неймана y, xy
double exact_gradient1166x(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double exact_gradient1166y(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double exact_gradient1166z(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double righthand1166(int dir_axis, double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

//////////////////////////////////////////////////////////////////////////////////
double exact_solution999(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient999x(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient999y(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient999z(double x, double y, double z, double t, double mx, double my, double mz);
double righthand999(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double suppl_xyz_function888 ( double x, double y, double z );
double suppl_xyzt_function888 ( double x, double y, double z, double t, double a888 );
double suppl_xyzt_function888_dx ( double x, double y, double z, double t, double a888 );
double suppl_xyzt_function888_dy ( double x, double y, double z, double t, double a888 );
double suppl_xyzt_function888_dz ( double x, double y, double z, double t, double a888 );
double suppl_function888 (double x, double y, double z, double t, double a888);
double exact_solution888(double x, double y, double z, double t, double a888, double deg888);
double exact_gradient888x(double x, double y, double z, double t, double a888, double deg888);
double exact_gradient888y(double x, double y, double z, double t, double a888, double deg888);
double exact_gradient888z(double x, double y, double z, double t, double a888, double deg888);
double suppl_xyzt_function888_dxdx ( double x, double y, double z, double t, double a888 );
double suppl_xyzt_function888_dydy ( double x, double y, double z, double t, double a888 );
double suppl_xyzt_function888_dzdz ( double x, double y, double z, double t, double a888 );
double suppl_xyzt_function888_dt( double x, double y, double z, double t, double a888 );
double exact_gradient888xx(double x, double y, double z, double t, double a888, double deg888);
double exact_gradient888yy(double x, double y, double z, double t, double a888, double deg888);
double exact_gradient888zz(double x, double y, double z, double t, double a888, double deg888);
double exact_laplace888(double x, double y, double z, double t, double a888, double deg888);
double exact_derivdt888(double x, double y, double z, double t, double a888, double deg888);
double righthand888(double a888, double deg888, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);


double suppl_x_function8 ( double x, double y, double z, int main_axis );
double suppl_x_function8_d ( double x, double y, double z, int main_axis );
double suppl_x_function8_dd ( double x, double y, double z, int main_axis );
double suppl_xyzt_function8 ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dx ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dy ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dz ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dxdx ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dydy ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dzdz ( double x, double y, double z, double t, double a888, int main_axis );
double suppl_xyzt_function8_dt( double x, double y, double z, double t, double a888, int main_axis );
double suppl_function8 (double x, double y, double z, double t, double a888, int main_axis);
double exact_solution8(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_gradient8x(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_gradient8y(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_gradient8z(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_gradient8xx(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_gradient8yy(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_gradient8zz(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_laplace8(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double exact_derivdt8(double x, double y, double z, double t, double a888, double deg888, int main_axis);
double righthand8(int main_axis, double a888, double deg888, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double function_g(double s, int n);
double function_g_deriv(double s, int n);
double function_g_deriv2(double s, int n);
double exact_solution777(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_gradient777x(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_gradient777y(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_gradient777z(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_gradient777xx(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_gradient777yy(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_gradient777zz(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double righthand777(int dir_axis, double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double exact_solution777_onlytime(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);
double exact_laplace777(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis);

double exact_solution776(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient776x(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient776y(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient776z(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient776xx(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient776yy(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient776zz(double x, double y, double z, double t, double mx, double my, double mz);
double righthand776(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double exact_solution776_onlytime(double x, double y, double z, double t, double mx, double my, double mz);
double exact_laplace776(double x, double y, double z, double t, double mx, double my, double mz);

double exact_solution500(double x, double y, double z, double t, int m);
double exact_gradient500x(double x, double y, double z, double t, int m);
double exact_gradient500y(double x, double y, double z, double t, int m);
double exact_gradient500z(double x, double y, double z, double t, int m);
double exact_solution500_po_t(double x, double y, double z, double t, int m);
double exact_laplace500(double x, double y, double z, double t, int m);
double righthand500(int m, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution501(double x, double y, double z, double t);
double exact_gradient501x(double x, double y, double z, double t);
double exact_gradient501y(double x, double y, double z, double t);
double exact_gradient501z(double x, double y, double z, double t);
double exact_solution501_po_t(double x, double y, double z, double t);
double exact_laplace501(double x, double y, double z, double t);
double righthand501(int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double exact_solution_501_averageQuad ( int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints,
									   int NQuadpoints );


double exact_solution889(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_gradient889x(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_gradient889y(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_gradient889z(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_gradient889xx(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_gradient889yy(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_gradient889zz(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_laplace889(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double exact_derivdt889(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis);
double righthand889(int main_axis, int mx, int my, int mz, double a888, double deg888, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution444(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_gradient444x(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_gradient444y(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_gradient444z(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_gradient444xx(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_gradient444yy(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_laplace444(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double exact_derivdt444(double x, double y, double z, double t, double mx, double my, double mz, double deg444);
double righthand444(int mx, int my, int mz, double deg444, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);



double Tx_0_bound(int numSolution, int i, int j, int k, double time);
double Tx_1_bound(int numSolution, int i, int j, int k, double time);
double Ty_0_bound(int numSolution, int i, int j, int k, double time);
double Ty_1_bound(int numSolution, int i, int j, int k, double time);
double Tz_0_bound(int numSolution, int i, int j, int k, double time);
double Tz_1_bound(int numSolution, int i, int j, int k, double time);
double wx_0_bound(int numSolution, int i, int j, int k, double time);
double wx_1_bound(int numSolution, int i, int j, int k, double time);
double wy_0_bound(int numSolution, int i, int j, int k, double time);
double wy_1_bound(int numSolution, int i, int j, int k, double time);
double wz_0_bound(int numSolution, int i, int j, int k, double time);
double wz_1_bound(int numSolution, int i, int j, int k, double time);


void setBoundaryConditions(int number)
{
//  eDirichlet = 0,
//  eNeumann,
//  eMixed

 switch(number)
  {
  case 1:
      xCondition = eNeumann;
      yCondition = eNeumann;
      zCondition = eNeumann;
    break;
  case 9066:
      xCondition = eDirichlet;
      yCondition = eNeumann;
      zCondition = eNeumann;
	break;
  case 9006:
	  xCondition = eNeumann;
	  yCondition = eNeumann;
	  zCondition = eNeumann;
	  break;
  case 1111:
	  if (mx != 0)
	  {
		xCondition = eDirichlet;
		yCondition = eNeumann;
		zCondition = eNeumann;
	  }
	  else if (my != 0)
	  {
		xCondition = eNeumann;
		yCondition = eDirichlet;
		zCondition = eNeumann;
	  }
	  else
	  {
		xCondition = eNeumann;
		yCondition = eNeumann;
		zCondition = eDirichlet;
	  }
	  break;
  case 1166:
	  switch(Dir_axis)
	  {
	  case 0:
		  {
			  xCondition = eDirichlet;
			  yCondition = eNeumann;
			  zCondition = eNeumann;
		  }
		  break;
	  case 1:
		  {
			  xCondition = eNeumann;
			  yCondition = eDirichlet;
			  zCondition = eNeumann;
		  }
		  break;
	  case 2:
		  {
			  xCondition = eNeumann;
			  yCondition = eNeumann;
			  zCondition = eDirichlet;
		  }
		  break;
	  }
	  break;
  case 999:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 888:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 889:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 444:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 8:
	  switch(Dir_axis)
	  {
	  case 0:
		  {
			  xCondition = eDirichlet;
			  yCondition = eNeumann;
			  zCondition = eNeumann;
		  }
		  break;
	  case 1:
		  {
			  xCondition = eNeumann;
			  yCondition = eDirichlet;
			  zCondition = eNeumann;
		  }
		  break;
	  case 2:
		  {
			  xCondition = eNeumann;
			  yCondition = eNeumann;
			  zCondition = eDirichlet;
		  }
		  break;
	  }
	  break;
  case 777:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 776:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 500:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  case 501:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eDirichlet;
	  break;
  default:
    printf("FUCK YOU ALL (hello from default case of setBoundaryConditions \n");
    break;
  }
}

double exact_solution(int number, double x, double y, double z, double t) 
{
	switch(number)
	{
	case 1:
		return 1;
		break;
	case 9066:
		return exact_solution9066(x, y, z, t, mx, my, mz);
		break;
	case 9006:
		return exact_solution9006(x, y, z, t, mx, my, mz);
		break;
	case 1166:
		return exact_solution1166(x, y, z, t, mx, my, mz, Dir_axis);
		break;
	case 1111:
		return exact_solution1111(x, y, z, t, mx, my, mz);
		break;
	case 999:
		return exact_solution999(x, y, z, t, mx, my, mz);
		break;
	case 888:
		return exact_solution888(x, y, z, t, A888, DEG888  );
		break;
	case 889:
		return exact_solution889(x, y, z, t, A888, DEG888, mx, my, mz, Dir_axis  );
		break;
	case 8:
		return exact_solution8(x, y, z, t, A888, DEG888, Dir_axis  );
		break;
	case 777:
		return exact_solution777(x, y, z, t, mx, my, mz, Dir_axis);
		break;
	case 776:
		return exact_solution776(x, y, z, t, mx, my, mz);
		break;
	case 500:
		return exact_solution500(x, y, z, t, mx);
		break;
	case 501:
		return exact_solution501(x, y, z, t);
		break;
	case 444:
		return exact_solution444(x, y, z, t, mx, my, mz, DEG444);
		break;
/*
  case 2:
    return 1 - z;
    break;

  case 101:
	  return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t) + Tz_0 + (Tz_1 - Tz_0)*z;
	  break;
  case 102:
	  return nonstationary*cos(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t) + Ty_0 + (Ty_1 - Ty_0)*y;
	  break;
  case 103:
	  return nonstationary*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x;
	  break;
  case 104:
	  //if mx > 0 - inhomog. Neumann for x, homog. Neumann for x,y; 
	  //if mx = 0, my > 0 - inhomog. Neumann for y, Neumann for z (no x, quasi 2D test);
	  //if mx = 0, my = 0 - inhomog. Neumann for z (no x,y, quasi 1D test);
	  return exact_solution1111(x, y, z, t, mx, my, mz);
	  break;
  case 105:
	  //Dir_axis (0,1,2 -> x,y,z) defines which direction has Dirichlet boundary condition, others have Neumann bnd cnds.
	  //e.g. Dir_axis = 1 means that y-axis has Dirichlet conditions, x- and z-axis - Neumann.
	  return exact_solution1166(x,y,z,t,mx,my,mz, Dir_axis);
	  break;
  case 305:
	  //Dir_axis (0,1,2 -> x,y,z) defines which direction has Dirichlet boundary condition, others have Neumann bnd cnds.
	  //e.g. Dir_axis = 1 means that y-axis has Dirichlet conditions, x- and z-axis - Neumann.
	  if (Dir_axis == 1 && external == 4)
		  return exact_solution267(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 300:
	  //Dir_axis (0,1,2 -> x,y,z) defines which direction has Dirichlet boundary condition, others have Neumann bnd cnds.
	  //e.g. Dir_axis = 1 means that y-axis has Dirichlet conditions, x- and z-axis - Neumann.
	  if (Dir_axis == 1)
		  return exact_solution300(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 400:
	  //Dir_axis (0,1,2 -> x,y,z) defines which direction has Dirichlet boundary condition, others have Neumann bnd cnds.
	  //e.g. Dir_axis = 1 means that y-axis has Dirichlet conditions, x- and z-axis - Neumann.
	  if (Dir_axis == 1)
		  return exact_solution400(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 106:
	  //works for certain mx,my,mz, not for all values (special test).
	  //i suppose it works for mx = my = 2, mz = 1.
	  return exact_solution1199(x,y,z,t,mx,my,mz);
	  break;
  case 107:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return exact_solution206(x,y,z,t,mx,my);
    break;
  case 108:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbSmooth_ex1_solution(x,y,z,t,mx,my,mt);
    break;
  case 109:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbSmooth_ex2_solution(x,y,z,t,mx,my,mt);
    break;
  case 110:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbNSmooth_ex3_solution(x,y,z,t,mx,my,mt,alpha);
    break;
  case 111:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbNSmooth_ex4_solution(x,y,z,t,mx,my,mt,alpha);
    break;
  case 112:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbNSmooth_ex5_solution(x,y,z,t,mx,my,mt,alpha);
    break;
  case 113:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbNSmooth_ex6_solution(x,y,z,t,mx,my,mt,alpha);
    break;
  case 114:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbNSmooth_ex7_solution(x,y,z,t,mx,my,mt,alpha);
    break;
  case 115:
    //Neumann test with inhomogeneous Neumann for y (impossible for 1111)
    return arbNSmooth_ex8_solution(x,y,z,t,mx,my,mt,alpha);
    break;

  case 10:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
    break;
  case 3:
    return sin(PI*x)*exp(-t);
    break;
  case 4:
    return cos(PI*x)*exp(-t);
    break;
  case 5:
    return cos(2*PI*x)*exp(-t);
    break;
  case 6:
    return nonstationary*cos(mx*PI*x)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
    break;
*/
  default:
    return 0;
    break;
  }
  return 0;
}

double exact_solution_averageQuad ( int number, double time, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints,
									   int NQuadpoints )
{
	double res = 0.0;

	double x1 = xpoints[i];
	double x2 = xpoints[i+1];
	double y1 = ypoints[j];
	double y2 = ypoints[j+1];
	double z1 = zpoints[k];
	double z2 = zpoints[k+1];

	double volume = (x2 - x1) * (y2 - y1) * (z2 - z1);

	double hxQuad = (x2 - x1) / NQuadpoints;
	double hyQuad = (y2 - y1) / NQuadpoints;
	double hzQuad = (z2 - z1) / NQuadpoints;

//	printf ( "x1 = %3.3f x2 = %3.3f y1 = %3.3f y2 = %3.3f z1 = %3.3f z2 = %3.3f \n", x1, x2, y1, y2, z1, z2 );
//	printf ( "hxQuad = %3.3f, hyQuad = %3.3f hzQuad = %3.3f \n", hxQuad, hyQuad, hzQuad );
//	printf ( "volume = %f \n", volume );

	double xmidQuad, ymidQuad, zmidQuad, volumeQuad;

	for ( int kQuad = 0; kQuad < NQuadpoints; kQuad++ )
	{
		zmidQuad = z1 +  (kQuad + 0.5) * hzQuad;
		for ( int jQuad = 0; jQuad < NQuadpoints; jQuad++ )
		{
			ymidQuad = y1 + (jQuad + 0.5) * hyQuad;
			for ( int iQuad = 0; iQuad < NQuadpoints; iQuad++ )
			{
//				printf ( "iQuad = %d jQuad = %d kQuad = %d \n", iQuad, jQuad, kQuad );
				xmidQuad = x1 + (iQuad + 0.5) * hxQuad;
				volumeQuad = hxQuad * hyQuad * hzQuad;
				res += volumeQuad * exact_solution ( number, xmidQuad, ymidQuad, zmidQuad, time );
//				printf ( "xmidQuad = %3.3f ymidQuad = %3.3f zmidQuad = %3.3f \n", xmidQuad, ymidQuad, zmidQuad );
//				printf ( "res = %f \n", res );
			}
		}
	}

	res /= volume;
	return res;
}

//double exact_solutionT0(int number, double x, double y, double z, double t) 
//{
//	switch(number)
//	{
//	case 501:
//		return exact_solution501(x, y, z, t);
//		break;
//	default:
//		printf ( "exact_solutionAverage() is not implemented for numsol = %d \n", number );
//		return 0;
//		break;
//  }
//  return 0;
//}
double exact_gradientX(int number, double x, double y, double z, double t) 
{
  switch(number)
  {
  case 1:
    return 0;
    break;
  case 9006:
	  return exact_gradient9006x(x,y,z,t, mx, my, mz);
	  break;
  case 9066:
	  return exact_gradient9066x(x,y,z,t, mx, my, mz);
	  break;
  case 1166:
	  return exact_gradient1166x(x,y,z,t, mx, my, mz, Dir_axis);
	  break;
  case 1111:
	  return exact_gradient1111x(x,y,z,t, mx, my, mz);
	  break;
  case 999:
	  return exact_gradient999x(x,y,z,t, mx, my, mz);
	  break;
  case 888:
	  return exact_gradient888x(x,y,z,t, A888, DEG888);
	  break;
  case 889:
	  return exact_gradient889x(x,y,z,t, A888, DEG888, mx, my, mz, Dir_axis);
	  break;
  case 8:
	  return exact_gradient8x(x,y,z,t, A888, DEG888, Dir_axis);
	  break;
  case 777:
	  return exact_gradient777x(x,y,z,t, mx, my, mz, Dir_axis);
	  break;
  case 776:
	  return exact_gradient776x(x,y,z,t, mx, my, mz);
	  break;
  case 500:
	  return exact_gradient500x(x,y,z,t, mx);
	  break;
  case 501:
	  return exact_gradient501x(x,y,z,t);
	  break;
  case 444:
	  return exact_gradient444x(x,y,z,t, mx, my, mz, DEG444);
	  break;
/*
  case 2:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
    break;

  case 101:
	  return nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
	  break;
  case 102:
	  return nonstationary*(-mx*PI)*sin(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
	  break;
  case 103:
	  return nonstationary*(mx*PI)*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (Tx_1 - Tx_0) ;
	  break;
  case 104:
	  return exact_gradient1111x(x,y,z,t,mx,my,mz);
	  break;
  case 105:
	  return exact_gradient1166x(x,y,z,t,mx,my,mz, Dir_axis);
	  break;
  case 305:
	  if (Dir_axis == 1 && external == 4)
		  return exact_gradient267x(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 300:
	  if (Dir_axis == 1)
		  return exact_gradient300x(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 400:
	  if (Dir_axis == 1)
		  return exact_gradient400x(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 106:
	  return exact_gradient1199x(x,y,z,t,mx,my,mz);
	  break;
  case 107:
    return exact_gradient206x(x,y,z,t,mx,my);
    break;
  case 108:
    return arbSmooth_ex1_gradientX(x,y,z,t,mx,my, mt);
    break;
  case 109:
    return arbSmooth_ex2_gradientX(x,y,z,t,mx,my, mt);
    break;
  case 110:
    return arbNSmooth_ex3_gradientX(x,y,z,t,mx,my,mt,alpha);
    break;
  case 111:
    return arbNSmooth_ex4_gradientX(x,y,z,t,mx,my,mt,alpha);
    break;
  case 112:
    return arbNSmooth_ex5_gradientX(x,y,z,t,mx,my,mt,alpha);
    break;
  case 113:
    return arbNSmooth_ex6_gradientX(x,y,z,t,mx,my,mt,alpha);
    break;
  case 114:
    return arbNSmooth_ex7_gradientX(x,y,z,t,mx,my,mt,alpha);
    break;
  case 115:
    return arbNSmooth_ex8_gradientX(x,y,z,t,mx,my,mt,alpha);
    break;

  case 3:
    return sin(PI*x)*exp(-t);
    break;
  case 4:
    return cos(PI*x)*exp(-t);
    break;
  case 5:
    return cos(2*PI*x)*exp(-t);
    break;
  case 6:
    return nonstationary*cos(mx*PI*x)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
    break;
*/
  default:
	return 0;
    break;
  }
  return 0;
}

double exact_gradientY(int number, double x, double y, double z, double t) 
{
  switch(number)
  {
  case 1:
    return 0;
    break;
  case 9006:
	  return exact_gradient9006y(x,y,z,t, mx, my, mz);
	  break;
  case 9066:
	  return exact_gradient9066y(x,y,z,t, mx, my, mz);
	  break;
  case 1166:
	  return exact_gradient1166y(x,y,z,t, mx, my, mz, Dir_axis);
	  break;
  case 1111:
	  return exact_gradient1111y(x,y,z,t, mx, my, mz);
	  break;
  case 999:
	  return exact_gradient999y(x,y,z,t, mx, my, mz);
	  break;
  case 888:
	  return exact_gradient888y(x,y,z,t, A888, DEG888);
	  break;
  case 889:
	  return exact_gradient889y(x,y,z,t, A888, DEG888, mx, my, mz, Dir_axis);
	  break;
  case 8:
	  return exact_gradient8y(x,y,z,t, A888, DEG888, Dir_axis);
	  break;
  case 777:
	  return exact_gradient777y(x,y,z,t, mx, my, mz, Dir_axis);
	  break;
  case 776:
	  return exact_gradient776y(x,y,z,t, mx, my, mz);
	  break;
  case 500:
	  return exact_gradient500y(x,y,z,t, mx);
	  break;
  case 501:
	  return exact_gradient501y(x,y,z,t);
	  break;
  case 444:
	  return exact_gradient444y(x,y,z,t, mx, my, mz, DEG444);
	  break;
/*
  case 2:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
	break;

  case 101:
	  return nonstationary*cos(mx*PI*x)*(-my*PI)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
	  break;
  case 102:
	  return nonstationary*cos(mx*PI*x)*(my*PI)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (Ty_1 - Ty_0);
	  break;
  case 103:
	  return nonstationary*sin(mx*PI*x)*(-my*PI)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
	  break;
  case 104:
	  return exact_gradient1111y(x,y,z,t,mx,my,mz);
	  break;
  case 105:
	  return exact_gradient1166y(x,y,z,t,mx,my,mz,Dir_axis);
	  break;
  case 305:
	  if (Dir_axis == 1 && external == 4)
		  return exact_gradient267y(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 300:
	  if (Dir_axis == 1)
		  return exact_gradient300y(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 400:
	  if (Dir_axis == 1)
		  return exact_gradient400y(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 106:
	  return exact_gradient1199y(x,y,z,t,mx,my,mz);
	  break;
  case 107:
    return exact_gradient206y(x,y,z,t,mx,my);
    break;
  case 108:
    return arbSmooth_ex1_gradientY(x,y,z,t,mx,my, mt);
    break;
  case 109:
    return arbSmooth_ex2_gradientY(x,y,z,t,mx,my, mt);
    break;
  case 110:
    return arbNSmooth_ex3_gradientY(x,y,z,t,mx,my,mt,alpha);
    break;
  case 111:
    return arbNSmooth_ex4_gradientY(x,y,z,t,mx,my,mt,alpha);
    break;
  case 112:
    return arbNSmooth_ex5_gradientY(x,y,z,t,mx,my,mt,alpha);
    break;
  case 113:
    return arbNSmooth_ex6_gradientY(x,y,z,t,mx,my,mt,alpha);
    break;
  case 114:
    return arbNSmooth_ex7_gradientY(x,y,z,t,mx,my,mt,alpha);
    break;
  case 115:
    return arbNSmooth_ex8_gradientY(x,y,z,t,mx,my,mt,alpha);
    break;


  case 3:
    return sin(PI*x)*exp(-t);
    break;
  case 4:
    return cos(PI*x)*exp(-t);
    break;
  case 5:
    return cos(2*PI*x)*exp(-t);
    break;
  case 6:
    return nonstationary*cos(mx*PI*x)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
    break;
*/
  default:
    return 0;
    break;
  }
  return 0;
}

double exact_gradientZ(int number, double x, double y, double z, double t) 
{
 switch(number)
  {
  case 1:
    return 0;
    break;
  case 9006:
	  return exact_gradient9006z(x,y,z,t, mx, my, mz);
	  break;
  case 9066:
	  return exact_gradient9066z(x,y,z,t, mx, my, mz);
	  break;
  case 1166:
	  return exact_gradient1166z(x,y,z,t, mx, my, mz, Dir_axis);
	  break;
  case 1111:
	  return exact_gradient1111z(x,y,z,t, mx, my, mz);
	  break;
  case 999:
	  return exact_gradient999z(x,y,z,t, mx, my, mz);
	  break;
  case 888:
	  return exact_gradient888z(x,y,z,t,A888,DEG888);
	  break;
  case 889:
	  return exact_gradient889z(x,y,z,t,A888,DEG888,mx,my,mz,Dir_axis);
	  break;
  case 8:
	  return exact_gradient8z(x,y,z,t,A888,DEG888, Dir_axis);
	  break;
  case 777:
	  return exact_gradient777z(x,y,z,t, mx, my, mz, Dir_axis);
	  break;
  case 776:
	  return exact_gradient776z(x,y,z,t, mx, my, mz);
	  break;
  case 500:
	  return exact_gradient500z(x,y,z,t, mx);
	  break;
  case 501:
	  return exact_gradient501z(x,y,z,t);
	  break;
  case 444:
	  return exact_gradient444z(x,y,z,t, mx, my, mz, DEG444);
	  break;
/*
  case 2:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
    break;

  case 101:
	  return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*(mz*PI)*cos(mz*PI*z)*exp(-t) + (Tz_1 - Tz_0);
	  break;
  case 102:
	  return nonstationary*cos(mx*PI*x)*sin(my*PI*y)*(-mz*PI)*sin(mz*PI*z)*exp(-t) ;
	  break;
  case 103:
	  return nonstationary*sin(mx*PI*x)*cos(my*PI*y)*(-mz*PI)*sin(mz*PI*z)*exp(-t);
	  break;
  case 104:
	  return exact_gradient1111z(x,y,z,t,mx,my,mz);
	  break;
  case 105:
	  return exact_gradient1166z(x,y,z,t,mx,my,mz,Dir_axis);
	  break;
  case 106:
	  return exact_gradient1199z(x,y,z,t,mx,my,mz);
	  break;
  case 107:
    return 0.0;
    break;
  case 108:
    return 0.0;
    break;
  case 109:
    return 0.0;
    break;

  case 3:
    return sin(PI*x)*exp(-t);
    break;
  case 4:
    return cos(PI*x)*exp(-t);
    break;
  case 5:
    return cos(2*PI*x)*exp(-t);
    break;
  case 6:
    return nonstationary*cos(mx*PI*x)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
    break;
*/
  default:
	return 0;
    break;
  }
  return 0;
}


double righthand(int number, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
#ifdef SPECIAL_MODE
	return 0;
#endif
#ifdef SPECIALTEST
	return 0.0;
#else
  switch(number)
  {
  case 1:
	  return 0;
	  break;
  case 9066:
	  return righthand9066(mx, my, mz, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 1166:
	  return righthand1166(Dir_axis, mx, my, mz, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 1111:
	  return righthand1111(mx, my, mz, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 999:
	  return righthand999(mx, my, mz, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 888:
	  return righthand888 (A888, DEG888, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 889:
	  return righthand889 (Dir_axis,mx,my,mz, A888, DEG888, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 8:
	  return righthand8 (Dir_axis, A888, DEG888, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 777:
	  return righthand777(Dir_axis, mx, my, mz, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 776:
	  return righthand776(mx, my, mz, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 500:
	  return righthand500(mx, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 501:
	  return righthand501(t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
  case 444:
	  return righthand444(mx, my, mz, DEG444, t_discr, i, j, k, xpoints, ypoints, zpoints);
	  break;
/*
  case 101:
	  return nonstationary*((-1)*heat_capacity(i,j,k)*heat_capacity(i,j,k) - (-1)*heat_conductivity(i,j,k)*(mx*PI*mx*PI + my*PI*my*PI + mz*PI*mz*PI))*(1.0/(mx*PI*my*PI*(-1)*mz*PI))*(sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]))*(sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]))*(cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]))*exp(-t_discr*tau) ;
	  break;
  case 102:
	  return nonstationary*((-1)*heat_capacity(i,j,k)*heat_capacity(i,j,k) - (-1)*heat_conductivity(i,j,k)*(mx*PI*mx*PI + my*PI*my*PI + mz*PI*mz*PI))*(1.0/(mx*PI*(-1)*my*PI*mz*PI))*(sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]))*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]))*(sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]))*exp(-t_discr*tau) ;
	  //return nonstationary*()*cos(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t) ;
	  break;
  case 103:
	  return nonstationary*((-1)*heat_capacity(i,j,k)*heat_capacity(i,j,k) - (-1)*heat_conductivity(i,j,k)*(mx*PI*mx*PI + my*PI*my*PI + mz*PI*mz*PI))*(1.0/((-1)*mx*PI*my*PI*mz*PI))*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]))*(sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]))*(sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]))*exp(-t_discr*tau) ;
	  //return nonstationary*()*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) ;
	  break;
  case 104:
	  return righthand1111(mx,my,mz,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;
  case 105:
	  return righthand1166(Dir_axis,mx,my,mz,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;
  case 305:
	  if (Dir_axis == 1 && external == 4)
		  return righthand267(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  else
		  return 0.0;
	  break;
  case 300:
	  if (Dir_axis == 1)
		  return righthand300(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  else
		  return 0.0;
	  break;
  case 400:
	  if (Dir_axis == 1)
		  return righthand400(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  else
		  return 0.0;
	  break;
  case 106:
	  return righthand1199(mx,my,mz,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;
  case 107:
    return righthand206(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 108:
    return arbSmooth_ex1_righthand(mx,my,mt,t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 109:
    return arbSmooth_ex2_righthand(mx,my,mt,t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 110:
    return arbNSmooth_ex3_righthand(mx,my,mt,alpha, t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 111:
    return arbNSmooth_ex4_righthand(mx,my,mt,alpha, t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 112:
    return arbNSmooth_ex5_righthand(mx,my,mt,alpha, t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 113:
    return arbNSmooth_ex6_righthand(mx,my,mt,alpha, t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 114:
    return arbNSmooth_ex7_righthand(mx,my,mt,alpha, t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
  case 115:
    return arbNSmooth_ex8_righthand(mx,my,mt,alpha, t_discr,i,j,k,xpoints,ypoints, zpoints);
    break;
*/
  default:
	  break;
  }
#endif
  return 0;
}



double exact_solution1(double x, double y, double z, double t)
{
	return 10;
}

double exact_solution2(double x, double y, double z, double t)
{
	if (x<=0.5)
		return 0.25*cos(PI*x)*exp(-t);
	else
		return cos(PI*x)*exp(-t);
}

double exact_solution3(double x, double y, double z, double t)
{
	return sin(PI*x)*exp(-t);
}

double exact_solution4(double x, double y, double z, double t)
{
	return cos(PI*x)*exp(-t);
}

double exact_solution5(double x, double y, double z, double t)
{
	return cos(2*PI*x)*exp(-t);
}

///////////////////////////////////////////////////////////////////////
double exact_solution6(double x, double y, double z, double t, int mx)
{
	return nonstationary*cos(mx*PI*x)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
}
double exact_gradient6(double x, double y, double z, double t, int mx)
{
	return nonstationary*(-mx*PI)*sin(mx*PI*x)*exp(-t) + x*wx_1/(-heat_conductivity(Nx-1,0,0)) + (1-x)*wx_0/(-heat_conductivity(0,0,0));
}
double righthand6(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,j,0)));
	return hy(j)*hz(k)*temp;
	//	return temp;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double exact_solution2006(double x, double y, double z, double t, int mz)
{
	return nonstationary*cos(mz*PI*z)*exp(-t) + (z*z/2)*wz_1/(-heat_conductivity(0,0,Nz-1)) + (z-z*z/2)*wz_0/(-heat_conductivity(0,0,0));
}
double exact_gradient2006(double x, double y, double z, double t, int mz)
{
	return nonstationary*(-mz*PI)*sin(mz*PI*z)*exp(-t) + z*wz_1/(-heat_conductivity(0,0,Nz-1)) + (1-z)*wz_0/(-heat_conductivity(0,0,0));
}
double righthand2006(int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*mz*PI - density(i,j,k)*heat_capacity(i,j,k)/(mz*PI) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hz(k)*(wz_1/(-heat_conductivity(i,j,Nz-1))  - wz_0/(-heat_conductivity(i,j,0)));
	return hy(j)*hx(i)*temp;
	//	return temp;
}
////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////
double exact_solution2066(double x, double y, double z, double t, int mz)
{
	return nonstationary*sin(mz*PI*z)*exp(-t) + Tz_0 + (Tz_1 - Tz_0)*z;
}
double exact_gradient2066(double x, double y, double z, double t, int mz)
{
	return nonstationary*(mz*PI)*cos(mz*PI*z)*exp(-t) + (Tz_1 - Tz_0);
}
double righthand2066(int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*mz*PI - density(i,j,k)*heat_capacity(i,j,k)/(mz*PI) ) * (-1)* ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) ;
	return hx(i)*hy(j)*temp;
}
///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
double exact_solution106(double x, double y, double z, double t, int mx)
{
	return nonstationary*cos(mx*PI*y)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
}
double exact_gradient106(double x, double y, double z, double t, int mx)
{
	return nonstationary*(-mx*PI)*sin(mx*PI*y)*exp(-t) + y*wy_1/(-heat_conductivity(0,Ny-1,0)) + (1-y)*wy_0/(-heat_conductivity(0,0,0));
}
double righthand106(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*ypoints[j+1]) - sin(mx*PI*ypoints[j]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hy(j)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
	return hx(i)*hz(k)*temp;
	//	return temp;
}

///////////////////////////////////////////////////////////////////////
double exact_solution206(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*cos(mx*PI*x)*cos(n*PI*y)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
}
double exact_gradient206x(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(n*PI*y)*exp(-t); 
}
double exact_gradient206y(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(-n*PI)*sin(n*PI*y)*cos(mx*PI*x)*exp(-t) + y*wy_1/(-heat_conductivity(0,Ny-1,0)) + (1-y)*wy_0/(-heat_conductivity(0,0,0));
}

double righthand206(int mx, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
	return hz(k)*temp;
	//	return temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
double exact_solution207(double x, double y, double z, double t, int mx, int n)//для разрывного посередине к-та теплопроводности (lambda_0->lambda_1)
{
	if (y<0.5)
	{
		return nonstationary*cos(mx*PI*x)*cos(n*PI*y)*exp(-t) + (y*y)*(-wy_0-wy_0)/(-heat_conductivity(0,0,0)) + y*wy_0/(-heat_conductivity(0,0,0));
		printf("zna4enie1 = %f \n",nonstationary*cos(mx*PI*x)*cos(n*PI*y)*exp(-t) + (y*y)*(-wy_0-wy_0)/(-heat_conductivity(0,0,0)) + y*wy_0/(-heat_conductivity(0,0,0)));
	}
	else
	{
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*cos(mx*PI*x)*cos(n*PI*y)*exp(-t) + (y-0.5)*(y-0.5)*(wy_1+wy_0)/(-heat_conductivity(0,Ny-1,0)) + (y-0.5)*(-wy_0)/(-heat_conductivity(0,Ny-1,0));
		printf("zna4enie2 = %f \n",(heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*cos(mx*PI*x)*cos(n*PI*y)*exp(-t) + (y-0.5)*(y-0.5)*(wy_1+wy_0)/(-heat_conductivity(0,Ny-1,0)) + (y-0.5)*(-wy_0)/(-heat_conductivity(0,Ny-1,0)));
	}
}
double exact_gradient207y(double x, double y, double z, double t, int mx, int n)
{
	if (y<0.5)
		return nonstationary*cos(mx*PI*x)*(-n*PI)*sin(n*PI*y)*exp(-t) + (2*y)*(-wy_0-wy_0)/(-heat_conductivity(0,0,0)) + wy_0/(-heat_conductivity(0,0,0));
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*cos(mx*PI*x)*(-n*PI)*nonstationary*sin(n*PI*y)*exp(-t) + 2*(y-0.5)*(wy_1+wy_0)/(-heat_conductivity(0,Ny-1,0)) + (-wy_0)/(-heat_conductivity(0,Ny-1,0));
	//	return nonstationary*(-mx*PI)*sin(mx*PI*x)*exp(-t) + x*w_1/(-heat_conductivity(Nx-1)) + (1-x)*w_0/(-heat_conductivity(0));
}
double exact_gradient207x(double x, double y, double z, double t, int mx, int n)
{
	if (y<0.5)
		return nonstationary*cos(n*PI*y)*(-mx*PI)*sin(mx*PI*x)*exp(-t); 
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*cos(n*PI*y)*(-mx*PI)*nonstationary*sin(mx*PI*x)*exp(-t) ;
	//	return nonstationary*(-mx*PI)*sin(mx*PI*x)*exp(-t) + x*w_1/(-heat_conductivity(Nx-1)) + (1-x)*w_0/(-heat_conductivity(0));
}
double righthand207(int mx, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = 0;
	if (ypoints[j]<0.5)
		temp = nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*2*(-wy_0-wy_0)/(-heat_conductivity(0,0,0));  
	else
		temp = (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*2*(wy_1+wy_0)/(-heat_conductivity(0,Ny-1,0));
	//	double temp = nonstationary*(heat_conductivity(k)*mx*PI - density(k)*heat_capacity(k)/(mx*PI) ) * ( sin(mx*PI*xpoints[k+1]) - sin(mx*PI*xpoints[k]) ) * exp(-t_discr*tau) - heat_conductivity(k)*hx(k)*(w_1/(-heat_conductivity(Nx-1))  - w_0/(-heat_conductivity(0)));
	return temp;
}
//////////////////////////////////////////////
double exact_solution266(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*sin(n*PI*y)*cos(mx*PI*x)*exp(-t) + Ty_0 + (Ty_1 - Ty_0)*y;
}
double exact_gradient266y(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(n*PI)*cos(n*PI*y)*cos(mx*PI*x)*exp(-t) + (Ty_1 - Ty_0);
}
double exact_gradient266x(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(-mx*PI)*sin(n*PI*y)*sin(mx*PI*x)*exp(-t);
}
double righthand266(int mx, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) ;
	return hz(k)*temp;
}
//////////////////////////////////////////////
double exact_solution267(double x, double y, double z, double t, int mx, int n)
{
	if (y<0.5)
		return nonstationary*sin(n*PI*y)*cos(mx*PI*x)*exp(-t) + Ty_0 + 2*y*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*sin(n*PI*y)*cos(mx*PI*x)*exp(-t) + (heat_conductivity(0,0,0)*Ty_0 + heat_conductivity(0,Ny-1,0)*Ty_1)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0)) + 2*(y-0.5)*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient267y(double x, double y, double z, double t, int mx, int n)
{
	if (y<0.5)
	{
		return (n*PI)*nonstationary*cos(n*PI*y)*cos(mx*PI*x)*exp(-t) + 2*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
		//		printf("heat_conductivity(0,Ny-1,0) = %f \n",heat_conductivity(0,Ny-1,0));
	}
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*(n*PI)*nonstationary*cos(n*PI*y)*cos(mx*PI*x)*exp(-t) + 2*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient267x(double x, double y, double z, double t, int mx, int n)
{
	if (y<0.5)
		return nonstationary*(-mx*PI)*sin(n*PI*y)*sin(mx*PI*x)*exp(-t); 
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-mx*PI)*sin(n*PI*y)*sin(mx*PI*x)*exp(-t); 
}
double righthand267(int mx, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = 0;
	if (ypoints[j]<0.5)
		temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) ;
	else
		temp = (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau);
	return hz(k)*temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double exact_solution10(double x, double y, double z, double t, int mx)
{
	return exp(-t);
}
double exact_gradient10(double x, double y, double z, double t, int mx)
{
	return 0;
}
double righthand10(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	return -hx(i)*density(i,j,k)*heat_capacity(i,j,k) * exp(-t_discr*tau); 
}

double exact_solution11(double x, double y, double z, double t, int mx)
{
	return cos(mx*PI*x);
}
double exact_gradient11(double x, double y, double z, double t, int mx)
{
	return -mx*PI*sin(mx*PI*x);
}
double righthand11(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	return heat_conductivity(i,j,k)*mx*PI * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) ;
}

//////////////////////////////////////////////
double exact_solution66(double x, double y, double z, double t, int n)
{
	return nonstationary*sin(n*PI*y)*exp(-t) + Ty_0 + (Ty_1 - Ty_0)*y;
}
double exact_gradient66y(double x, double y, double z, double t, int n)
{
	return nonstationary*(n*PI)*cos(n*PI*y)*exp(-t) + (Ty_1 - Ty_0);
}
double righthand66(int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*n*PI - density(i,j,k)*heat_capacity(i,j,k)/(n*PI) ) * (-1)* ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * exp(-t_discr*tau) ;
	return hx(i)*hz(k)*temp;
}
//////////////////////////////////////////////
//////////////////////////////////////////////
double exact_solution366(double x, double y, double z, double t, int mx)
{
	return nonstationary*sin(mx*PI*x)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x;
}
double exact_gradient366x(double x, double y, double z, double t, int mx)
{
	return nonstationary*(mx*PI)*cos(mx*PI*x)*exp(-t) + (Tx_1 - Tx_0);
}
double righthand366(int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * (-1)* ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) ;
	return hy(j)*hz(k)*temp;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
double exact_solution1266(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*sin(n*PI*y)*cos(mx*PI*z)*exp(-t) + Ty_0 + (Ty_1 - Ty_0)*y;
}
double exact_gradient1266y(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(n*PI)*cos(n*PI*y)*cos(mx*PI*z)*exp(-t) + (Ty_1 - Ty_0);
}
double exact_gradient1266z(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(-mx*PI)*sin(n*PI*y)*sin(mx*PI*z)*exp(-t);
}
double righthand1266(int mx, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * ( sin(mx*PI*zpoints[k+1]) - sin(mx*PI*zpoints[k]) ) * exp(-t_discr*tau) ;
	return hx(i)*temp;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
double exact_solution1466(double x, double y, double z, double t, int mz, int n)
{
	return nonstationary*cos(n*PI*y)*sin(mz*PI*z)*exp(-t) + Tz_0 + (Tz_1 - Tz_0)*z;
}
double exact_gradient1466z(double x, double y, double z, double t, int mz, int n)
{
	return nonstationary*(mz*PI)*cos(mz*PI*z)*cos(n*PI*y)*exp(-t) + (Tz_1 - Tz_0);
}
double exact_gradient1466y(double x, double y, double z, double t, int mz, int n)
{
	return nonstationary*(-n*PI)*sin(n*PI*y)*sin(mz*PI*z)*exp(-t);
}
double righthand1466(int mz, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mz*mz + n*n)/(mz*n) - density(i,j,k)*heat_capacity(i,j,k)/((mz*PI)*(n*PI)) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) ;
	return hx(i)*temp;
}
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////
double exact_solution1566(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*sin(mx*PI*x)*cos(n*PI*y)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x;
}
double exact_gradient1566x(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(mx*PI)*cos(mx*PI*x)*cos(n*PI*y)*exp(-t) + (Tx_1 - Tx_0);
}
double exact_gradient1566y(double x, double y, double z, double t, int mx, int n)
{
	return nonstationary*(-n*PI)*sin(mx*PI*x)*sin(n*PI*y)*exp(-t);
}
double righthand1566(int mx, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + n*n)/(mx*n) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * exp(-t_discr*tau) ;
	return hz(k)*temp;
}
//////////////////////////////////////////////

double exact_solution1206(double x, double y, double z, double t, int mz, int n)
{
	return nonstationary*cos(mz*PI*z)*cos(n*PI*y)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
}
double exact_gradient1206z(double x, double y, double z, double t, int mz, int n)
{
	return nonstationary*(-mz*PI)*sin(mz*PI*z)*cos(n*PI*y)*exp(-t); 
}
double exact_gradient1206y(double x, double y, double z, double t, int mz, int n)
{
	return nonstationary*(-n*PI)*sin(n*PI*y)*cos(mz*PI*z)*exp(-t) + y*wy_1/(-heat_conductivity(0,Ny-1,0)) + (1-y)*wy_0/(-heat_conductivity(0,0,0));
}

double righthand1206(int mz, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*(mz*mz + n*n)/(mz*n) - density(i,j,k)*heat_capacity(i,j,k)/((mz*PI)*(n*PI)) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hy(j)*hz(k)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
	return hx(i)*temp;
	//	return temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double exact_solution2206(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*cos(mz*PI*z)*cos(mx*PI*x)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
}
double exact_gradient2206z(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*(-mz*PI)*sin(mz*PI*z)*cos(mx*PI*x)*exp(-t); 
}
double exact_gradient2206x(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(mz*PI*z)*exp(-t) + x*wx_1/(-heat_conductivity(Nx-1,0,0)) + (1-x)*wx_0/(-heat_conductivity(0,0,0));
}

double righthand2206(int mz, int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*(mz*mz + mx*mx)/(mz*mx) - density(i,j,k)*heat_capacity(i,j,k)/((mz*PI)*(mx*PI)) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hz(k)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(0,j,k)));
	return hy(j)*temp;
	//	return temp;
}
///////////////////////////////////////////////////////

double exact_solution1366(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*sin(mz*PI*z)*cos(mx*PI*x)*exp(-t) + Tz_0 + (Tz_1 - Tz_0)*z;
}
double exact_gradient1366z(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*(mz*PI)*cos(mz*PI*z)*cos(mx*PI*x)*exp(-t) + (Tz_1 - Tz_0);
}
double exact_gradient1366x(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*(-mz*PI)*sin(mz*PI*z)*sin(mx*PI*x)*exp(-t);
}
double righthand1366(int mz, int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mz*mz + mx*mx)/(mz*mx) - density(i,j,k)*heat_capacity(i,j,k)/((mz*PI)*(mx*PI)) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) ;
	return hy(j)*temp;
}

///////////////////////////////////////////////////////
double exact_solution2366(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*sin(mx*PI*x)*cos(mz*PI*z)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x;
}
double exact_gradient2366x(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*(mx*PI)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t) + (Tx_1 - Tx_0);
}
double exact_gradient2366z(double x, double y, double z, double t, int mz, int mx)
{
	return nonstationary*(-mz*PI)*sin(mx*PI*x)*sin(mz*PI*z)*exp(-t);
}
double righthand2366(int mz, int mx, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mz*mz + mx*mx)/(mz*mx) - density(i,j,k)*heat_capacity(i,j,k)/((mz*PI)*(mx*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) ;
	return hy(j)*temp;
}


///////////////////////////////////////////////////////////////////////
double exact_solution9006(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*cos(mx*PI*x)*cos(n*PI*y)*cos(mz*PI*z)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
}
double exact_gradient9006x(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(n*PI*y)*cos(mz*PI*z)*exp(-t); 
}
double exact_gradient9006y(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*(-n*PI)*sin(n*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t) + y*wy_1/(-heat_conductivity(0,Ny-1,0)) + (1-y)*wy_0/(-heat_conductivity(0,0,0));
}

double exact_gradient9006z(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*cos(mx*PI*x)*cos(n*PI*y)*(-mz*PI)*sin(mz*PI*z)*exp(-t); 
}
double righthand9006(int mx, int n, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n + mz*mz)/(mx*n*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)*(mz*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
	return temp;
	//	return temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double exact_solution9066(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*sin(mx*PI*x)*cos(n*PI*y)*cos(mz*PI*z)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x; 
}
double exact_gradient9066x(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*(mx*PI)*cos(mx*PI*x)*cos(n*PI*y)*cos(mz*PI*z)*exp(-t) + (Tx_1 - Tx_0); 
}
double exact_gradient9066y(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*(-n*PI)*sin(n*PI*y)*sin(mx*PI*x)*cos(mz*PI*z)*exp(-t); 
}

double exact_gradient9066z(double x, double y, double z, double t, int mx, int n, int mz)
{
	return nonstationary*sin(mx*PI*x)*cos(n*PI*y)*(-mz*PI)*sin(mz*PI*z)*exp(-t); 
}
double righthand9066(int mx, int n, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = (-1)*nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n + mz*mz)/(mx*n*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)*(mz*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) ;
	return temp;
	//	return temp;
}
///////////////////////////////////////////////////////

double Tx_0_bound(int numSolution, int i, int j, int k, double time)
{
  switch (numSolution)
  {
  /*case 108:
    return arbSmooth_ex1_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 109:
    return arbSmooth_ex2_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 110:
    return arbNSmooth_ex3_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 111:
    return arbNSmooth_ex4_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 112:
    return arbNSmooth_ex5_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 113:
    return arbNSmooth_ex6_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 114:
    return arbNSmooth_ex7_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 115:
    return arbNSmooth_ex8_solution(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;*/
  default:
    return exact_solution(numSolution,0,ypoints[j] + 0.5*hy(j),zpoints[k] + 0.5 * hz(k),time);
    break;
  }
}

double Tx_1_bound(int numSolution, int i, int j, int k, double time)
{
  switch (numSolution)
  {
 /* case 108:
    return arbSmooth_ex1_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 109:
    return arbSmooth_ex2_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 110:
    return arbNSmooth_ex3_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 111:
    return arbNSmooth_ex4_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 112:
    return arbNSmooth_ex5_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 113:
    return arbNSmooth_ex6_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 114:
    return arbNSmooth_ex7_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 115:
    return arbNSmooth_ex8_solution(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;*/
  default:
    return exact_solution(numSolution,1,ypoints[j] + 0.5*hy(j),zpoints[k] + 0.5 * hz(k),time);
    break;
  }
}


double Ty_0_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
 /* case 108:
#ifndef STEKLOV
    return arbSmooth_ex1_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt);
#else
	return arbSmooth_ex1_solution_other_boundTy0(i,j,time);  
#endif
    break;
  case 109:
    return arbSmooth_ex2_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt);
    break;
  case 110:
    return arbNSmooth_ex3_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt, alpha);
    break;
  case 111:
    return arbNSmooth_ex4_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt, alpha);
    break;
  case 112:
    return arbNSmooth_ex5_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt, alpha);
    break;
  case 113:
    return arbNSmooth_ex6_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt, alpha);
    break;
  case 114:
    return arbNSmooth_ex7_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt, alpha);
    break;
  case 115:
    return arbNSmooth_ex8_solution(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt, alpha);
    break;*/
  default:
    return exact_solution(numSolution,xpoints[i] + 0.5*hx(i),0,zpoints[k] + 0.5 * hz(k),time);
    break;
  }
}

double Ty_1_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
 /* case 108:
#ifndef STEKLOV
    return arbSmooth_ex1_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt);
#else
	return arbSmooth_ex1_solution_other_boundTy1(i,j,time);
#endif
    break;
  case 109:
    return arbSmooth_ex2_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt);
    break;
  case 110:
    return arbNSmooth_ex3_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt, alpha);
    break;
  case 111:
    return arbNSmooth_ex4_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt, alpha);
    break;
  case 112:
    return arbNSmooth_ex5_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt, alpha);
    break;
  case 113:
    return arbNSmooth_ex6_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt, alpha);
    break;
  case 114:
    return arbNSmooth_ex7_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt, alpha);
    break;
  case 115:
    return arbNSmooth_ex8_solution(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt, alpha);
    break;*/
  default:
    return exact_solution(numSolution,xpoints[i] + 0.5*hx(i),1,zpoints[k] + 0.5 * hz(k),time);
    break;
  }
}


double Tz_0_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
  default:
    return exact_solution(numSolution,xpoints[i] + 0.5*hx(i),ypoints[j] + 0.5*hy(j),0,time);
    break;
  }
}

double Tz_1_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
  default:
    return exact_solution(numSolution,xpoints[i] + 0.5*hx(i),ypoints[j] + 0.5*hy(j),1,time);
    break;
  }
}
double wx_0_bound(int numSolution, int i, int j, int k, double time)
{
  switch (numSolution)
  {
 // case 108:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
	////return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX_integ(i, j, time,mx,my,mt);
 //   break;
 // case 109:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex2_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
 //   break;
 // case 110:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex3_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
 //   break;
 // case 111:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex4_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
 //   break;
 // case 112:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex5_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
 //   break;
 // case 113:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex6_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
 //   break;
 // case 114:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex7_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
 //   break;
 // case 115:
 //   return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex8_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
 //   break;
  default:
    return wx_0;
    break;
  }
}

double wx_1_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
 // case 108:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
	////return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX_integ(i, j, time,mx,my,mt);
 //   break;
 // case 109:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex2_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
 //   break;
 // case 110:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex3_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
 //   break;
 // case 111:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex4_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
 //   break;
 // case 112:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex5_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
 //   break;
 // case 113:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex6_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
 //   break;
 // case 114:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex7_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
 //   break;
 // case 115:
 //   return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex8_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
 //   break;
  default:
    //return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*exact_gradientX(numSolution,1,ypoints[j] + 0.5*hy(j),0,time);
    return wx_1;
    break;
  }
}

double wy_0_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
  //case 108:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbSmooth_ex1_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt);
  //  break;
  //case 109:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbSmooth_ex2_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt);
  //  break;
  //case 110:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex3_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
  //  break;
  //case 111:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex4_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
  //  break;
  //case 112:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex5_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
  //  break;
  //case 113:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex6_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
  //  break;
  //case 114:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex7_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
  //  break;
  //case 115:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex8_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
  //  break;
  default:
    //return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*exact_gradientX(numSolution,xpoints[i] + 0.5*hx(i),0,0,time);
    return wy_0;
    break;
  }
}

double wy_1_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
  //case 108:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbSmooth_ex1_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt);
  //  break;
  //case 109:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbSmooth_ex2_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt);
  //  break;
  //case 110:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex3_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
  //  break;
  //case 111:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex4_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
  //  break;
  //case 112:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex5_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
  //  break;
  //case 113:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex6_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
  //  break;
  //case 114:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex7_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
  //  break;
  //case 115:
  //  return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex8_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
  //  break;
  default:
    //return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*exact_gradientX(numSolution,xpoints[i] + 0.5*hx(i),1,0,time);
    return wy_1;
    break;
  }
}



double wz_0_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
  default:
    return wz_0;
    break;
  }
}

double wz_1_bound(int numSolution, int i, int j, int k,  double time)
{
  switch (numSolution)
  {
  default:
    return wz_1;
    break;
  }
}


double v_init(double x, double y, double z)
{
	return x*(1-x)*y*(1-y)*z*(1-z);
}

double Vx_1(double x, double y, double z, double t)
{
	return vx_const;
}

double Vy_1(double x, double y, double z, double t)
{
	return vy_const;
}

double Vz_1(double x, double y, double z, double t)
{
	return vz_const;
}

double Vx_2(double x, double y, double z, double t)
{
	return x*(1 - x);
}

double Vy_2(double x, double y, double z, double t)
{
	return y*(1 - y);
}

double Vz_2(double x, double y, double z, double t)
{
	return 2.0*(x + y - 1)*z;
}

double Vx_3(double x, double y, double z, double t)
{
	return x;
}

double Vy_3(double x, double y, double z, double t)
{
	return y;
}

double Vz_3(double x, double y, double z, double t)
{
	return -2.0*z;
}


///////////////////////////////////////////////////////////////////////
double exact_solution1111(double x, double y, double z, double t, int mx, int my, int mz)
{
	if (mx != 0)
		return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
	else if (my != 0)
		return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
	else
		return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (z*z/2)*wz_1/(-heat_conductivity(0,0,Nz-1)) + (z-z*z/2)*wz_0/(-heat_conductivity(0,0,0));
}
double exact_solution1111Vx(double x, double y, double z, double t, int mx, int my, int mz)
{
	//	return vx_const;
	switch(velocityX_num)
	{
	case 1:
		return Vx_1(x,y,z,t);
		break;
	case 2:
		return Vx_2(x,y,z,t);
		break;
	case 3:
		return Vx_3(x,y,z,t);
		break;
	default:
		return 0;
	}
	//	return Vx_1(x,y,z,t);
	//	return Vx_2(x,y,z,t);
	//	return 0;
	//	return x*(1-x);
	//	if ((x > 0) && (x < 1))
	//		return 1;
	//	else
	//		return 0;
}
double exact_solution1111Vy(double x, double y, double z, double t, int mx, int my, int mz)
{
	//	return vy_const;
	switch(velocityY_num)
	{
	case 1:
		return Vy_1(x,y,z,t);
		break;
	case 2:
		return Vy_2(x,y,z,t);
		break;
	case 3:
		return Vy_3(x,y,z,t);
		break;
	default:
		return 0;
	}
	//	return y*(1-y);
	//	if ((y > 0) && (y < 1))
	//		return 1;
	//	else
	//		return 0;
}
double exact_solution1111Vz(double x, double y, double z, double t, int mx, int my, int mz)
{
	//	return vz_const;
	switch(velocityZ_num)
	{
	case 1:
		return Vz_1(x,y,z,t);
		break;
	case 2:
		return Vz_2(x,y,z,t);
		break;
	case 3:
		return Vz_3(x,y,z,t);
		break;
	default:
		return 0;
	}

	//	return Vz_1(x,y,z,t);
	//	return vz_const + z*(1-z);
	//	if ((z > 0) && (z < 1))
	//		return 1;
	//	else
	//		return 0;
}
double exact_gradient1111x(double x, double y, double z, double t, int mx, int my, int mz)
{
	if (mx > 0)
	{
		double lambdagrad = nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + x*wx_1/(-heat_conductivity(Nx-1,0,0)) + (1-x)*wx_0/(-heat_conductivity(0,0,0)); 
		return lambdagrad;
	}
	else
	{
		//double lambdagrad = + x*wx_1/(-heat_conductivity(Nx-1,0,0)) + (1-x)*wx_0/(-heat_conductivity(0,0,0)); 
		//return lambdagrad;
		return 0;
	}
	//	double convect = heat_capacity(x,y,z)*density(x,y,z)*exact_solution1111Vx(x,y,z,t,mx, my, mz)*exact_solution1111(x,y,z,t,mx, my, mz);
	//	return lambdagrad + convect;
}
double exact_gradient1111y(double x, double y, double z, double t, int mx, int my, int mz)
{
	if (my > 0)
	{
		double lambdagrad = nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t);
		if (mx > 0)
			return lambdagrad;
		else
		{
			lambdagrad += + y*wy_1/(-heat_conductivity(0,Ny-1,0)) + (1-y)*wy_0/(-heat_conductivity(0,0,0));
			return lambdagrad;
		}
	}
	else
	{
		return 0;
	}
	//	double convect = heat_capacity(x,y,z)*density(x,y,z)*exact_solution1111Vy(x,y,z,t,mx, my, mz)*exact_solution1111(x,y,z,t,mx, my, mz);
	//	return lambdagrad + convect;
}
double exact_gradient1111z(double x, double y, double z, double t, int mx, int my, int mz)
{
	if (mz > 0)
	{
		double lambdagrad = nonstationary*(-mz*PI)*sin(mz*PI*z)*cos(mx*PI*x)*cos(my*PI*y)*exp(-t);
		if ((mx > 0) || (my > 0))
			return lambdagrad;
		else
		{
			lambdagrad += + z*wz_1/(-heat_conductivity(0,0,Nz-1)) + (1-z)*wz_0/(-heat_conductivity(0,0,0));
			return lambdagrad;
		}
	} 
	else
	{
		return 0;
	}
	//	double convect = heat_capacity(x,y,z)*density(x,y,z)*exact_solution1111Vz(x,y,z,t,mx, my, mz)*exact_solution1111(x,y,z,t,mx, my, mz);
	//	return lambdagrad + convect;
}

double righthand1111(int mx, int my, int mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempbound, tempv, tempTw, tempdivvx, tempdivvy, tempdivvz;

	tempTw = nonstationary*(heat_conductivity(i,j,k)*(mx*mx + my*my + mz*mz)/(mx*my*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)*(mz*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,0,k)));

	if (mz == 0)
	{
		tempTw = nonstationary*(heat_conductivity(i,j,k)*(mx*mx + my*my)/(mx*my) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * hz(k) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,0,k)));
	} 

	if ((my == 0) && (mz == 0))
	{
		tempTw = nonstationary*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,0,k)));
		//tempTw = nonstationary*( - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,0,k)));
	}
	if ((mx == 0) && (mz == 0))
	{
		tempTw = nonstationary*(heat_conductivity(i,j,k)*my*PI - density(i,j,k)*heat_capacity(i,j,k)/(my*PI) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * hx(i) * hz(k) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
	} 
	if ((mx == 0) && (my == 0))
	{
		tempTw = nonstationary*(heat_conductivity(i,j,k)*mz*PI - density(i,j,k)*heat_capacity(i,j,k)/(mz*PI) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * hy(j) * hx(i) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wz_1/(-heat_conductivity(i,j,Nz-1))  - wz_0/(-heat_conductivity(i,j,0)));
		//tempTw = nonstationary*( - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*hz(k)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,0,k)));
	}

	//double temp = nonstationary*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,j,0)));
	//return hy(j)*hz(k)*temp;


	switch(velocityX_num)
	{
	case 0:
		tempvx = 0;
		break;
	case 1:
		//tempbound = vx_const*(0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
		//if ((my == 0) && (mz == 0))
		//	tempvx = vx_const*nonstationary*density(i,j,k)*heat_capacity(i,j,k)* ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau) + tempbound;
		//else
		//	tempvx = vx_const*nonstationary*(mx*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) + tempbound;
		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		break;
	case 2:
		//tempvnonstat = cell_volume * z1minusz * nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)* sin(mx*PI*(xpoints[i] + 0.5*hx(i))) * sin(my*PI*(ypoints[j] + 0.5*hy(j))) * cos(mz*PI*(zpoints[k] + 0.5*hz(k))) * exp(-t_discr*tau);
		//tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		if ((my == 0) && (mz == 0))
			tempdivvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1111(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) ;
		else
			tempdivvx = 0;
		tempvx += tempdivvx;

		break;
	case 3:
		//tempvnonstat = cell_volume * z1minusz * nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)* sin(mx*PI*(xpoints[i] + 0.5*hx(i))) * sin(my*PI*(ypoints[j] + 0.5*hy(j))) * cos(mz*PI*(zpoints[k] + 0.5*hz(k))) * exp(-t_discr*tau);
		//tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		break;
	default:
		tempvx = 0;
		break;
	}
	switch(velocityY_num)
	{
	case 0:
		tempvy = 0;
		break;
	case 1:
		/*
		if ((my == 0) && (mz == 0))
		tempvy = 0;
		else if ()
		{
		} 
		else
		{
		tempvy = vy_const*nonstationary*(my*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
		}
		*/
		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		break;
	case 2:
		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		if ((mx == 0) && (mz == 0))
			tempdivvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*y_mid) * exact_solution1111(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) ;
		else
			tempdivvy = 0;
		tempvy += tempdivvy;

		break;
	case 3:
		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		break;
	default:
		tempvy = 0;
		break;
	}
	switch(velocityZ_num)
	{
	case 0:
		tempvz = 0;
		break;
	case 1:
		/*
		if (((my == 0) && (mz == 0)) || ((mx == 0) && (mz == 0))) 
		tempvz = 0;
		else
		tempvz = vz_const*nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
		*/
		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		break;
	case 2:
		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);

		if ((mx == 0) && (my == 0))
			tempdivvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * 2.0*(x_mid + y_mid - 1) * exact_solution1111(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) ;
		else
			tempdivvz = 0;
		tempvz += tempdivvz;

		break;
	case 3:
		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1111z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
		break;
	default:
		tempvz = 0;
		break;
	}

	//double tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
	//double tempvx = nonstationary*(mx*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) + tempbound;
	//double tempvy = nonstationary*(my*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
	//double tempvz = nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
	//double tempv = vx_const*tempvx + vy_const*tempvy + vz_const*tempvz;
	tempv = tempvx + tempvy + tempvz;
	//if (((t_discr > 950) && (t_discr % 20 == 0))&&(((i==4)&&(j==4))&&(k==4)))
	//{
	//	printf("t_discr = %d tempv = %f \n",t_discr, tempv);
	//	_getch();
	//}
	//printf("tempvx = %f tempvy = %f tempvz = %f \n",tempvx,tempvy,tempvz);
	//printf("my = %d mz = %d \n", my, mz);
	//_getch();
	return tempTw + tempv;
	//double z1minusz = (zpoints[k] + 0.5*hz(k))*(1 - (zpoints[k] + 0.5*hz(k)));
	//double cell_volume = hx(i)*hy(j)*hz(k);
	//double tempvnonstat = cell_volume * z1minusz * nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)* sin(mx*PI*(xpoints[i] + 0.5*hx(i))) * sin(my*PI*(ypoints[j] + 0.5*hy(j))) * cos(mz*PI*(zpoints[k] + 0.5*hz(k))) * exp(-t_discr*tau);
	//return tempTw + tempv + tempvnonstat;
	//return tempTw;
}

//////////////////////////////////////////////////////////////////////////////////
double exact_solution1166(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x;
		break;
	case 1:
		return nonstationary*sin(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t) + Ty_0 + (Ty_1 - Ty_0)*y;
		break;
	case 2:
		return nonstationary*sin(mz*PI*z)*cos(mx*PI*x)*cos(my*PI*y)*exp(-t) + Tz_0 + (Tz_1 - Tz_0)*z;
		break;
	default:
		printf("FUCK YOU \n");
		return 0;
		break;
	}
	//if (mx != 0)
	//	return nonstationary*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + Tx_0 + (Tx_1 - Tx_0)*x;
	//else if (my != 0)
	//	return nonstationary*sin(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t) + Ty_0 + (Ty_1 - Ty_0)*y;
	//else
	//	return nonstationary*sin(mz*PI*z)*cos(mx*PI*x)*cos(my*PI*y)*exp(-t) + Tz_0 + (Tz_1 - Tz_0)*y;
}
///////////////////////////////////////////////////////////////////////
double exact_gradient1166x(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	double lambdagrad = 0;
	switch(dir_axis)
	{
	case 0:
		lambdagrad = nonstationary*(mx*PI)*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (Tx_1 - Tx_0);
		return lambdagrad;
		break;
	case 1:
		lambdagrad = nonstationary*(-mx*PI)*sin(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 2:
		lambdagrad = nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	default:
		printf("FUCK YOU \n");
		return 0;
		break;
	}

	//	double convect = heat_capacity(x,y,z)*density(x,y,z)*exact_solution1111Vx(x,y,z,t,mx, my, mz)*exact_solution1111(x,y,z,t,mx, my, mz);
	//	return lambdagrad + convect;
}

double exact_gradient1166y(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	double lambdagrad = 0;
	switch(dir_axis)
	{
	case 0:
		lambdagrad = nonstationary*(-my*PI)*sin(my*PI*y)*sin(mx*PI*x)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 1:
		lambdagrad = nonstationary*(my*PI)*cos(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t) + (Ty_1 - Ty_0);
		return lambdagrad;
		break;
	case 2:
		lambdagrad = nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*sin(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	default:
		printf("FUCK YOU \n");
		return 0;
		break;
	}

	//	double convect = heat_capacity(x,y,z)*density(x,y,z)*exact_solution1111Vy(x,y,z,t,mx, my, mz)*exact_solution1111(x,y,z,t,mx, my, mz);
	//	return lambdagrad + convect;
}

double exact_gradient1166z(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	double lambdagrad = 0;
	switch(dir_axis)
	{
	case 0:
		lambdagrad = nonstationary*(-mz*PI)*cos(my*PI*y)*sin(mx*PI*x)*sin(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 1:
		lambdagrad = nonstationary*(-mz*PI)*sin(my*PI*y)*cos(mx*PI*x)*sin(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 2:
		lambdagrad = nonstationary*(mz*PI)*cos(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t) + (Tz_1 - Tz_0);
		return lambdagrad;
		break;
	default:
		printf("FUCK YOU \n");
		return 0;
		break;
	}


	//if (mz > 0)
	//{
	//	double lambdagrad = nonstationary*(mz*PI)*cos(mz*PI*z)*cos(mx*PI*x)*cos(my*PI*y)*exp(-t);
	//	if ((mx > 0) || (my > 0))
	//		return lambdagrad;
	//	else
	//	{
	//		lambdagrad += + z*wz_1/(-heat_conductivity(0,0,Nz-1)) + (1-z)*wz_0/(-heat_conductivity(0,0,0));
	//		return lambdagrad;
	//	}
	//} 
	//else
	//{
	//	return 0;
	//}
	//	double convect = heat_capacity(x,y,z)*density(x,y,z)*exact_solution1111Vz(x,y,z,t,mx, my, mz)*exact_solution1111(x,y,z,t,mx, my, mz);
	//	return lambdagrad + convect;
}

double righthand1166(int dir_axis, double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempv, tempTw;//, tempbound, tempdivvx, tempdivvy, tempdivvz;

	switch(dir_axis)
	{
	case 0:
		if (mx == 0)
			tempTw = 0;
		else
		{
			tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my + mz*mz)/(mx*my*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)*(mz*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);

			if (mz == 0)
			{
				tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my)/(mx*my) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * hz(k) * exp(-t_discr*tau) ;
			} 

			if ((my == 0) && (mz == 0))
			{
				tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau);
				//tempTw = nonstationary*(-1)*(- density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau);
			}

		}
		break;
	case 1:
		if (my == 0)
			tempTw = 0;
		else
		{
			tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my + mz*mz)/(mx*my*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)*(mz*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);

			if (mz == 0)
			{
				tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my)/(mx*my) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * hz(k) * exp(-t_discr*tau) ;
			} 

			if ((mx == 0) && (mz == 0))
			{
				tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*my*PI - density(i,j,k)*heat_capacity(i,j,k)/(my*PI) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * hx(i) * hz(k) * exp(-t_discr*tau);
			}

		}
		break;
	case 2:
		if (mz == 0)
			tempTw = 0;
		else
		{
			tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my + mz*mz)/(mx*my*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)*(mz*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);

			if ((mx == 0) && (my == 0))
			{
				tempTw = nonstationary*(-1)*(heat_conductivity(i,j,k)*mz*PI - density(i,j,k)*heat_capacity(i,j,k)/(mz*PI) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * hx(i) * hy(j) * exp(-t_discr*tau);
			}

		}

		break;
	default:
		printf("FUCK YOU \n");
		return 0;
		break;
	}

	//double temp = nonstationary*(heat_conductivity(i,j,k)*mx*PI - density(i,j,k)*heat_capacity(i,j,k)/(mx*PI) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*(wx_1/(-heat_conductivity(Nx-1,j,k))  - wx_0/(-heat_conductivity(i,j,0)));
	//return hy(j)*hz(k)*temp;


	switch(velocityX_num)
	{
	case 0:
		tempvx = 0;
		break;
	case 1:
		//tempbound = vx_const*(0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
		//if ((my == 0) && (mz == 0))
		//	tempvx = vx_const*nonstationary*density(i,j,k)*heat_capacity(i,j,k)* ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * hy(j) * hz(k) * exp(-t_discr*tau) + tempbound;
		//else
		//	tempvx = vx_const*nonstationary*(mx*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) + tempbound;
		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//if (external == 1)
		//{
		//	//tempvx = cell_volume * spec_density(x_mid)*heat_capacity(i,j,k) * (Vx[i + j*(Nx+1) + k*(Nx+1)*Ny] + Vx[i + 1 + j*(Nx+1) + k*(Nx+1)*Ny]) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//	//v_exact * ro_exact = 1.
		//	tempvx = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vx(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//}
		break;
	case 2:
		//tempvnonstat = cell_volume * z1minusz * nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)* sin(mx*PI*(xpoints[i] + 0.5*hx(i))) * sin(my*PI*(ypoints[j] + 0.5*hy(j))) * cos(mz*PI*(zpoints[k] + 0.5*hz(k))) * exp(-t_discr*tau);
		//tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//if ((my == 0) && (mz == 0))
		//{
		//	if (external == 1)
		//	{
		//		tempvx = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//		tempdivvx = cell_volume * heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//	else
		//	{
		//		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//		tempdivvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//}
		//else if (mz == 0)
		//{
		//	if (external == 1)
		//	{
		//		tempvx = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//		tempdivvx = cell_volume * heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//	else
		//	{
		//		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//		tempdivvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//}
		//else
		//{
		//	tempdivvx = 0;
		//}
		//tempvx += tempdivvx;
		break;
	case 3:
		//tempvnonstat = cell_volume * z1minusz * nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)* sin(mx*PI*(xpoints[i] + 0.5*hx(i))) * sin(my*PI*(ypoints[j] + 0.5*hy(j))) * cos(mz*PI*(zpoints[k] + 0.5*hz(k))) * exp(-t_discr*tau);
		//tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
		tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		break;
	default:
		tempvx = 0;
		break;
	}
	switch(velocityY_num)
	{
	case 0:
		tempvy = 0;
		break;
	case 1:
		/*
		if ((my == 0) && (mz == 0))
		tempvy = 0;
		else if ()
		{
		} 
		else
		{
		tempvy = vy_const*nonstationary*(my*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
		}
		*/
		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);

		//if (external == 2)
		//{
		//	//v_exact * c_p_exact = 1.
		//	tempvy = cell_volume * density(i,j,k) * exact_solution1111Vy(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//}

		break;
	case 2:
		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//if ((mx == 0) && (mz == 0))
		//{
		//	if (external == 2)
		//	{
		//		//tempvx = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//		//tempdivvx = cell_volume * heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//		tempvy = cell_volume * density(i,j,k)* exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//		tempdivvy = cell_volume * density(i,j,k)* (1 - 2*y_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//	else
		//	{
		//		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//		tempdivvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*y_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//}
		//else if (mz == 0)
		//{
		//	if (external == 2)
		//	{
		//		//tempvx = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		//		//tempdivvx = cell_volume * heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//		tempvy = cell_volume * density(i,j,k)* exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//		tempdivvy = cell_volume * density(i,j,k)* (1 - 2*y_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//	else
		//	{
		//		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//		tempdivvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*y_mid) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//}
		//else
		//{
		//	tempdivvy = 0;
		//}
		//tempvy += tempdivvy;

		break;
	case 3:
		tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		break;
	default:
		tempvy = 0;
		break;
	}
	switch(velocityZ_num)
	{
	case 0:
		tempvz = 0;
		break;
	case 1:
		/*
		if (((my == 0) && (mz == 0)) || ((mx == 0) && (mz == 0))) 
		tempvz = 0;
		else
		tempvz = vz_const*nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
		*/
		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//if (external == 11)
		//{
		//	tempvz = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//}
		break;
	case 2:
		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);

		//if ((mx == 0) && (my == 0))
		//{
		//	if (external == 11)
		//	{
		//		tempvz = cell_volume * heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//		tempdivvz = cell_volume * heat_capacity(i,j,k) * 2.0*(x_mid + y_mid - 1) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//	else
		//	{
		//		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis);
		//		tempdivvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * 2.0*(x_mid + y_mid - 1) * exact_solution1166(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz, dir_axis) ;
		//	}
		//}
		//else
		//{
		//	tempdivvz = 0;
		//}
		//tempvz += tempdivvz;
		break;
	case 3:
		tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1166z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz,dir_axis);
		break;
	default:
		tempvz = 0;
		break;
	}

	//double tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
	//double tempvx = nonstationary*(mx*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) + tempbound;
	//double tempvy = nonstationary*(my*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
	//double tempvz = nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
	//double tempv = vx_const*tempvx + vy_const*tempvy + vz_const*tempvz;
	tempv = tempvx + tempvy + tempvz;
	//if (((t_discr > 950) && (t_discr % 20 == 0))&&(((i==4)&&(j==4))&&(k==4)))
	//{
	//	printf("t_discr = %d tempv = %f \n",t_discr, tempv);
	//	_getch();
	//}
	//printf("tempvx = %f tempvy = %f tempvz = %f \n",tempvx,tempvy,tempvz);
	//printf("my = %d mz = %d \n", my, mz);
	//_getch();
	return tempTw + tempv;
	//double z1minusz = (zpoints[k] + 0.5*hz(k))*(1 - (zpoints[k] + 0.5*hz(k)));
	//double cell_volume = hx(i)*hy(j)*hz(k);
	//double tempvnonstat = cell_volume * z1minusz * nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)* sin(mx*PI*(xpoints[i] + 0.5*hx(i))) * sin(my*PI*(ypoints[j] + 0.5*hy(j))) * cos(mz*PI*(zpoints[k] + 0.5*hz(k))) * exp(-t_discr*tau);
	//return tempTw + tempv + tempvnonstat;
	//return tempTw;
}

///////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
double exact_solution999(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*sin(mx*PI*x)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
}
///////////////////////////////////////////////////////////////////////
double exact_gradient999x(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*(mx*PI)*cos(mx*PI*x)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
}

double exact_gradient999y(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*sin(mx*PI*x)*(my*PI)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
}

double exact_gradient999z(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*sin(mx*PI*x)*sin(my*PI*y)*(mz*PI)*cos(mz*PI*z)*exp(-t);
}

double righthand999(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempv, tempTw;//, tempbound, tempdivvx, tempdivvy, tempdivvz;

	return nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my + mz*mz)/(mx*my*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)*(mz*PI)) ) *
		( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * 
		( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) *
		( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);

//	double temp = (-1)*nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n + mz*mz)/(mx*n*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)*(mz*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) ;

}

//////////////////////////////////////////////////////////////////////////////////
double function_g(double s, int n)
{
	return s * sin(PI*n*s);
}
double function_g_deriv(double s, int n)
{
	return sin(PI*n*s) + s * (PI * n) * cos(PI*n*s);
}
double function_g_deriv2(double s, int n)
{
	return 2*PI*n*cos(PI*n*s) - s * (PI * n) * (PI * n) * sin(PI*n*s);
}

double exact_solution777(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g(x,mx)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 1:
		return nonstationary*sin(mx*PI*x)*function_g(y,my)*sin(mz*PI*z)*exp(-t);
		break;
	case 2:
		return nonstationary*sin(mx*PI*x)*sin(my*PI*y)*function_g(z,mz)*exp(-t);
		break;
	default:
		printf ( "Wrong dir_axis = %d in solution 777 \n", dir_axis );
	}
	return 0.0;
}
///////////////////////////////////////////////////////////////////////
double exact_gradient777x(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g_deriv(x,mx)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 1:
		return nonstationary*(mx*PI)*cos(mx*PI*x)*function_g(y,my)*sin(mz*PI*z)*exp(-t);
		break;
	case 2:
		return nonstationary*(mx*PI)*cos(mx*PI*x)*sin(my*PI*y)*function_g(z,mz)*exp(-t);
		break;
	default:
		printf ( "Wrong dir_axis = %d in solution 777 \n", dir_axis );
	}
	return 0.0;
}

double exact_gradient777y(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g(x,mx)*(my*PI)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 1:
		return nonstationary*sin(mx*PI*x)*function_g_deriv(y,my)*sin(mz*PI*z)*exp(-t);
		break;
	case 2:
		return nonstationary*sin(mx*PI*x)*(my*PI)*cos(my*PI*y)*function_g(z,mz)*exp(-t);
		break;
	default:
		printf ( "Wrong dir_axis = %d in solution 777 \n", dir_axis );
	}
	return 0.0;
//	return nonstationary*sin(mx*PI*x)*(my*PI)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
}

double exact_gradient777z(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g(x,mx)*sin(my*PI*y)*(mz*PI)*cos(mz*PI*z)*exp(-t);
		break;
	case 1:
		return nonstationary*sin(mx*PI*x)*function_g(y,my)*(mz*PI)*cos(mz*PI*z)*exp(-t);
		break;
	case 2:
		return nonstationary*sin(mx*PI*x)*sin(my*PI*y)*function_g_deriv(z,mz)*exp(-t);
		break;
	default:
		printf ( "Wrong dir_axis = %d in solution 777 \n", dir_axis );
	}
	return 0.0;
//	return nonstationary*sin(mx*PI*x)*sin(my*PI*y)*(mz*PI)*cos(mz*PI*z)*exp(-t);
}

//double integral777_simple(double x, int m)
//{
//	if (m == 0)
//		return x;
//	else
//		return (1.0 / (PI * m)) * sin (PI * m * x);
//}
//double integral777_notsimple(double y, int n)
//{
//	if (n == 0)
//		return y;
//	else
//		return - (1.0 / (PI * n)) * y * cos (PI * n * y) + (1.0 / (PI * n * PI * n)) * sin (PI * n * y);
//}
double exact_gradient777xx(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g_deriv2(x,mx)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 1:
		return nonstationary*((-PI*mx*PI*mx) * sin(mx*PI*x))*function_g(y,my)*sin(mz*PI*z)*exp(-t);
		break;
	case 2:
		return nonstationary*((-PI*mx*PI*mx) * sin(mx*PI*x))*sin(my*PI*y)*function_g(z,mz)*exp(-t);
		break;
	}
	return 0.0;
}
double exact_gradient777yy(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g(x,mx)*((-PI*my*PI*my) * sin(my*PI*y))*sin(mz*PI*z)*exp(-t);
		break;
	case 1:
		return nonstationary*sin(mx*PI*x)*function_g_deriv2(y,my)*sin(mz*PI*z)*exp(-t);
		break;
	case 2:
		return nonstationary*sin(mx*PI*x)*((-PI*my*PI*my) * sin(my*PI*y))*function_g(z,mz)*exp(-t);
		break;
	}
	return 0.0;
}
double exact_gradient777zz(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return nonstationary*function_g(x,mx)*sin(my*PI*y)*((-PI*mz*PI*mz) * sin(mz*PI*z))*exp(-t);
		break;
	case 1:
		return nonstationary*sin(mx*PI*x)*function_g(y,my)*((-PI*mz*PI*mz) * sin(mz*PI*z))*exp(-t);
		break;
	case 2:
		return nonstationary*sin(mx*PI*x)*sin(my*PI*y)*function_g_deriv2(z,mz)*exp(-t);
		break;
	}
	return 0.0;
}
double exact_laplace777(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	return exact_gradient777xx(x,y,z,t,mx,my,mz,dir_axis) + exact_gradient777yy(x,y,z,t,mx,my,mz,dir_axis) + exact_gradient777zz(x,y,z,t,mx,my,mz,dir_axis);
}
double exact_solution777_onlytime(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	return exact_solution777(x,y,z,t,mx,my,mz,dir_axis);
}
double righthand777(int dir_axis, double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);

	//double integ_x, integ_y, integ_z;
	//switch (dir_axis)
	//{
	//case 0:
	//	integ_x = integral_notsimple(xpoints[i+1],mx) - integral_notsimple(xpoints[i],mx);
	//	integ_y = integral_simple(ypoints[j+1],my) - integral_simple(ypoints[j],my);
	//	integ_z = integral_simple(zpoints[k+1],mz) - integral_simple(zpoints[k],mz);
	//	break;
	//case 1:
	//	integ_x = integral_simple(xpoints[i+1],mx) - integral_simple(xpoints[i],mx);
	//	integ_y = integral_notsimple(ypoints[j+1],my) - integral_notsimple(ypoints[j],my);
	//	integ_z = integral_simple(zpoints[k+1],mz) - integral_simple(zpoints[k],mz);
	//	break;
	//case 2:
	//	integ_x = integral_simple(xpoints[i+1],mx) - integral_simple(xpoints[i],mx);
	//	integ_y = integral_simple(ypoints[j+1],my) - integral_simple(ypoints[j],my);
	//	integ_z = integral_notsimple(zpoints[k+1],mz) - integral_notsimple(zpoints[k],mz);
	//	break;
	//}
	//return nonstationary * exp (-t_discr * tau) * integ_x * integ_y * integ_z / volume;

	//return nonstationary*(-1)*(heat_conductivity(i,j,k)*(mx*mx + my*my + mz*mz)/(mx*my*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)*(mz*PI)) ) *
	//	( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * 
	//	( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) *
	//	( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);

//	double temp = (-1)*nonstationary*(heat_conductivity(i,j,k)*(mx*mx + n*n + mz*mz)/(mx*n*mz*PI) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(n*PI)*(mz*PI)) ) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(n*PI*ypoints[j+1]) - sin(n*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) ;
	double temp = cell_volume * (-density(i,j,k)*heat_capacity(i,j,k)*exact_solution777_onlytime(x_mid, y_mid, z_mid, t_discr * tau, mx, my, mz, dir_axis) -
		heat_conductivity(i,j,k)*exact_laplace777(x_mid, y_mid, z_mid, t_discr * tau, mx, my, mz, dir_axis));
	return temp;

}
///////////////////////////////////////////////////////////////////////////////////////////////
double exact_solution776(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g(x,mx)*function_g(y,my)*function_g(z,mz)*exp(-t);
}
///////////////////////////////////////////////////////////////////////
double exact_gradient776x(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g_deriv(x,mx)*function_g(y,my)*function_g(z,mz)*exp(-t);
}

double exact_gradient776y(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g(x,mx)*function_g_deriv(y,my)*function_g(z,mz)*exp(-t);
}

double exact_gradient776z(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g(x,mx)*function_g(y,my)*function_g_deriv(z,mz)*exp(-t);
}
double exact_gradient776xx(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g_deriv2(x,mx)*function_g(y,my)*function_g(z,mz)*exp(-t);
}
double exact_gradient776yy(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g(x,mx)*function_g_deriv2(y,my)*function_g(z,mz)*exp(-t);
}
double exact_gradient776zz(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*function_g(x,mx)*function_g(y,my)*function_g_deriv2(z,mz)*exp(-t);
}
double exact_laplace776(double x, double y, double z, double t, double mx, double my, double mz)
{
	return exact_gradient776xx(x,y,z,t,mx,my,mz) + exact_gradient776yy(x,y,z,t,mx,my,mz) + exact_gradient776zz(x,y,z,t,mx,my,mz);
}
double exact_solution776_onlytime(double x, double y, double z, double t, double mx, double my, double mz)
{
	return exact_solution776(x,y,z,t,mx,my,mz);
}
double righthand776(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);

	double temp = cell_volume * (-density(i,j,k)*heat_capacity(i,j,k)*exact_solution776_onlytime(x_mid, y_mid, z_mid, t_discr * tau, mx, my, mz) -
		heat_conductivity(i,j,k)*exact_laplace776(x_mid, y_mid, z_mid, t_discr * tau, mx, my, mz));
	return temp;

}
///////////////////////////////////////////////////////////////////////////////////////////////
double exact_solution500(double x, double y, double z, double t, int m)
{
	return nonstationary * cos ( MT500 * PI * t ) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) );
}
double exact_gradient500x(double x, double y, double z, double t, int m)
{
	return nonstationary * cos ( MT500 * PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * ( 1 - 2 * x ) * y * ( 1 - y ) * z * (1 - z);
}
double exact_gradient500y(double x, double y, double z, double t, int m)
{
	return nonstationary * cos ( MT500 * PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * ( 1 - 2 * y ) * z * (1 - z);
}
double exact_gradient500z(double x, double y, double z, double t, int m)
{
	return nonstationary * cos ( MT500 * PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * y * (1 - y) * ( 1 - 2 * z );
}
double exact_solution500_po_t(double x, double y, double z, double t, int m)
{
	return nonstationary * (- MT500 * PI) * sin ( PI * t ) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) );
}
double exact_laplace500(double x, double y, double z, double t, int m)
{
	double out = 0.0 ;
	out += nonstationary * cos ( MT500 * PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * ( - 2 ) * y * ( 1 - y ) * z * (1 - z);
	out += nonstationary * cos ( MT500 * PI * t ) * (-1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * ( 1 - 2 * x ) * (1 - 2 * x ) * y * ( 1 - y ) * y * ( 1 - y ) * z * (1 - z) * z * (1 - z);
	out += nonstationary * cos ( MT500 * PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * ( - 2 ) * z * (1 - z);
	out += nonstationary * cos ( MT500 * PI * t ) * ( -1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * x * ( 1 - x ) * ( 1 - 2 * y ) * ( 1 - 2 * y ) * z * (1 - z) * z * (1 - z); 
	out += nonstationary * cos ( MT500 * PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * ( - 2 ) * y * (1 - y);
	out += nonstationary * cos ( MT500 * PI * t ) * ( -1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * x * ( 1 - x ) * ( 1 - 2 * z ) * ( 1 - 2 * z ) * y * (1 - y) * y * (1 - y); 

	return out;
}
double righthand500(int m, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * ( density(i,j,k)*heat_capacity(i,j,k)*exact_solution500_po_t(x_mid, y_mid, z_mid, t_discr * tau, m) - heat_conductivity(i,j,k)*exact_laplace500(x_mid, y_mid, z_mid, t_discr * tau, m));
	return temp;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
double exact_solution501(double x, double y, double z, double t)
{
	return nonstationary * exp(-t) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) );
}
double exact_gradient501x(double x, double y, double z, double t)
{
	return nonstationary * exp(-t) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * ( 1 - 2 * x ) * y * ( 1 - y ) * z * (1 - z);
}
double exact_gradient501y(double x, double y, double z, double t)
{
	return nonstationary * exp(-t) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * ( 1 - 2 * y ) * z * (1 - z);
}
double exact_gradient501z(double x, double y, double z, double t)
{
	return nonstationary * exp(-t) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * x * ( 1 - x ) * y * (1 - y) * ( 1 - 2 * z );
}
double exact_solution501_po_t(double x, double y, double z, double t)
{
	return nonstationary * (-1) * exp(-t) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) );
}
double exact_laplace501(double x, double y, double z, double t)
{
	double out = 0.0 ;
	out += nonstationary * exp(-t) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * ( - 2 ) * y * ( 1 - y ) * z * (1 - z);
	out += nonstationary * exp(-t) * (-1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) *
		( 1 - 2 * x ) * (1 - 2 * x ) * y * ( 1 - y ) * y * ( 1 - y ) * z * (1 - z) * z * (1 - z);
	out += nonstationary * exp(-t) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) * 
		x * ( 1 - x ) * ( - 2 ) * z * (1 - z);
	out += nonstationary * exp(-t) * (-1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) *
		x * ( 1 - x ) * x * ( 1 - x ) * ( 1 - 2 * y ) * ( 1 - 2 * y ) * z * (1 - z) * z * (1 - z); 
	out += nonstationary * exp(-t) * cos ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) *
		x * ( 1 - x ) * ( - 2 ) * y * (1 - y);
	out += nonstationary * exp(-t) * (-1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) * z * (1 - z) ) *
		x * ( 1 - x ) * x * ( 1 - x ) * ( 1 - 2 * z ) * ( 1 - 2 * z ) * y * (1 - y) * y * (1 - y); 

	return out;
}
double righthand501(int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * ( density(i,j,k)*heat_capacity(i,j,k)*exact_solution501_po_t(x_mid, y_mid, z_mid, t_discr * tau) - heat_conductivity(i,j,k)*exact_laplace501(x_mid, y_mid, z_mid, t_discr * tau));
	return temp;
}

//// computes Frenel C function up to the tolerance of 3.0e-6 for |x| < 0.63
//double Frenel_C ( double x )
//{
//	return x - pow(PI,2) * 0.025 * pow(x,5) + pow(PI,4) * pow(x,9) / 3456.0;
//}
//
//// computes Frenel S function up to the tolerance of 3.0e-6 for |x| < 0.73
//double Frenel_S ( double x )
//{
//	return PI/5.0 - pow(PI,3) * pow(x,7) / 336.0 + pow(PI,5) * pow(x,11) / 42240.0;
//}
//
//// computes integral of sin(x*(1-x)*y*(1-y)*z*(1-z)) with undefined constant taken as 0
//// expression taken from Mathematica
//double exact_integralmean5hundred( double x, double y, dou)
//{
//	double argument_forC = 
//}

double exact_solution_501_averageQuad ( int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints,
									   int NQuadpoints )
{
	double res = 0.0;

	double x1 = xpoints[i];
	double x2 = xpoints[i+1];
	double y1 = ypoints[j];
	double y2 = ypoints[j+1];
	double z1 = zpoints[k];
	double z2 = zpoints[k+1];

	double volume = (x2 - x1) * (y2 - y1) * (z2 - z1);

	double hxQuad = (x2 - x1) / NQuadpoints;
	double hyQuad = (y2 - y1) / NQuadpoints;
	double hzQuad = (z2 - z1) / NQuadpoints;

//	printf ( "x1 = %3.3f x2 = %3.3f y1 = %3.3f y2 = %3.3f z1 = %3.3f z2 = %3.3f \n", x1, x2, y1, y2, z1, z2 );
//	printf ( "hxQuad = %3.3f, hyQuad = %3.3f hzQuad = %3.3f \n", hxQuad, hyQuad, hzQuad );
//	printf ( "volume = %f \n", volume );

	double xmidQuad, ymidQuad, zmidQuad, volumeQuad;

	for ( int kQuad = 0; kQuad < NQuadpoints; kQuad++ )
	{
		zmidQuad = z1 +  (kQuad + 0.5) * hzQuad;
		for ( int jQuad = 0; jQuad < NQuadpoints; jQuad++ )
		{
			ymidQuad = y1 + (jQuad + 0.5) * hyQuad;
			for ( int iQuad = 0; iQuad < NQuadpoints; iQuad++ )
			{
//				printf ( "iQuad = %d jQuad = %d kQuad = %d \n", iQuad, jQuad, kQuad );
				xmidQuad = x1 + (iQuad + 0.5) * hxQuad;
				volumeQuad = hxQuad * hyQuad * hzQuad;
				res += volumeQuad * exact_solution501 ( xmidQuad, ymidQuad, zmidQuad, t_discr * tau );
//				printf ( "xmidQuad = %3.3f ymidQuad = %3.3f zmidQuad = %3.3f \n", xmidQuad, ymidQuad, zmidQuad );
//				printf ( "res = %f \n", res );
			}
		}
	}

	res /= volume;
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////
double suppl_xyz_function888 ( double x, double y, double z )
{
	return x * ( 1 - x ) * y * ( 1 - y ) * z * ( 1 - z );
}

double suppl_xyz_function888_dx ( double x, double y, double z )
{
	return (1 - 2.0 * x) * y * ( 1 - y ) * z * ( 1 - z );
}
double suppl_xyz_function888_dy ( double x, double y, double z )
{
	return (1 - 2.0 * y) * x * ( 1 - x ) * z * ( 1 - z );
}
double suppl_xyz_function888_dz ( double x, double y, double z )
{
	return (1 - 2.0 * z) * y * ( 1 - y ) * x * ( 1 - x );
}
double suppl_xyz_function888_dxdx ( double x, double y, double z )
{
	return (- 2.0) * y * ( 1 - y ) * z * ( 1 - z );
}
double suppl_xyz_function888_dydy ( double x, double y, double z )
{
	return (- 2.0) * x * ( 1 - x ) * z * ( 1 - z );
}
double suppl_xyz_function888_dzdz ( double x, double y, double z )
{
	return (- 2.0) * y * ( 1 - y ) * x * ( 1 - x );
}

double suppl_xyzt_function888 ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888(x,y,z) * exp(-t) - a888;
}
double suppl_xyzt_function888_dx ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888_dx(x,y,z) * exp(-t);
}
double suppl_xyzt_function888_dy ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888_dy(x,y,z) * exp(-t);
}
double suppl_xyzt_function888_dz ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888_dz(x,y,z) * exp(-t);
}
double suppl_xyzt_function888_dxdx ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888_dxdx(x,y,z) * exp(-t);
}
double suppl_xyzt_function888_dydy ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888_dydy(x,y,z) * exp(-t);
}
double suppl_xyzt_function888_dzdz ( double x, double y, double z, double t, double a888 )
{
	return C888 * suppl_xyz_function888_dzdz(x,y,z) * exp(-t);
}

double suppl_xyzt_function888_dt( double x, double y, double z, double t, double a888 )
{
	return (-1.0) * C888 * suppl_xyz_function888(x,y,z) * exp(-t);
}
double suppl_function888 (double x, double y, double z, double t, double a888)
{
	return fabs(suppl_xyzt_function888(x,y,z,t,a888));
}
double exact_solution888(double x, double y, double z, double t, double a888, double deg888)
{
	return pow(suppl_function888(x,y,z,t,a888), deg888) - pow(a888, deg888);
}
double exact_solution888noconst(double x, double y, double z, double t, double a888, double deg888)
{
	return pow(suppl_function888(x,y,z,t,a888), deg888);
}
///////////////////////////////////////////////////////////////////////
double exact_gradient888x(double x, double y, double z, double t, double a888, double deg888)
{
	double res = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function888_dx (x,y,z,t,a888);
	return res;
}
double exact_gradient888y(double x, double y, double z, double t, double a888, double deg888)
{
	double res = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function888_dy (x,y,z,t,a888);
	return res;
}
double exact_gradient888z(double x, double y, double z, double t, double a888, double deg888)
{
	double res = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function888_dz (x,y,z,t,a888);
	return res;
}
///////////////////////////////////////////////////////////////////////
double exact_gradient888xx(double x, double y, double z, double t, double a888, double deg888)
{
	double res1 = deg888 * (deg888 - 1) * exact_solution888noconst(x, y, z, t, a888, deg888 - 2);
	res1 *= pow(suppl_xyzt_function888_dx (x,y,z,t,a888), 2);

	double res2 = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res2 *= -1.0;
	res2 *= suppl_xyzt_function888_dxdx (x,y,z,t,a888);

	return res1 + res2;
}
double exact_gradient888yy(double x, double y, double z, double t, double a888, double deg888)
{
	double res1 = deg888 * (deg888 - 1) * exact_solution888noconst(x, y, z, t, a888, deg888 - 2);
	res1 *= pow(suppl_xyzt_function888_dy (x,y,z,t,a888), 2);

	double res2 = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res2 *= -1.0;
	res2 *= suppl_xyzt_function888_dydy (x,y,z,t,a888);

	return res1 + res2;
}
double exact_gradient888zz(double x, double y, double z, double t, double a888, double deg888)
{
	double res1 = deg888 * (deg888 - 1) * exact_solution888noconst(x, y, z, t, a888, deg888 - 2);
	res1 *= pow(suppl_xyzt_function888_dz (x,y,z,t,a888), 2);

	double res2 = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res2 *= -1.0;
	res2 *= suppl_xyzt_function888_dzdz (x,y,z,t,a888);

	return res1 + res2;
}

///////////////////////////////////////////////////////////////////////
double exact_laplace888(double x, double y, double z, double t, double a888, double deg888)
{
	return exact_gradient888xx (x, y, z, t, a888, deg888 ) + exact_gradient888yy (x, y, z, t, a888, deg888 ) + exact_gradient888zz (x, y, z, t, a888, deg888 );
}
double exact_derivdt888(double x, double y, double z, double t, double a888, double deg888)
{
	double res = deg888 * exact_solution888noconst(x, y, z, t, a888, deg888 - 1);
	if ( suppl_xyzt_function888(x, y, z, t, a888) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function888_dt (x,y,z,t,a888);

	return res;
}
double righthand888(double a888, double deg888, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempv, tempTw;//, tempbound, tempdivvx, tempdivvy, tempdivvz;

	return cell_volume *
		( density(i,j,k) * heat_capacity(i,j,k) * exact_derivdt888(x_mid, y_mid, z_mid, t_discr*tau, a888, deg888) +
		(-1) * heat_conductivity(i,j,k) * exact_laplace888(x_mid, y_mid, z_mid, t_discr*tau, a888, deg888) ); 
}

//--------------------------------------------------------
double suppl_x_function8 ( double x, double y, double z, int main_axis )
{
	switch (main_axis)
	{
	case 0:
		return x * (1 - x);
		break;
	case 1:
		return y * (1 - y);
		break;
	case 2:
		return z * (1 - z);
		break;
	default:
		printf ( "Wrong main axis = %d \n", main_axis );
		return 0;
		break;
	}
	return 0;
}

double suppl_x_function8_d ( double x, double y, double z, int main_axis )
{
	switch (main_axis)
	{
	case 0:
		return 1.0 - 2.0 * x;
		break;
	case 1:
		return 1.0 - 2.0 * y;
		break;
	case 2:
		return 1.0 - 2.0 * z;
		break;
	default:
		printf ( "Wrong main axis = %d \n", main_axis );
		return 0;
		break;
	}
	return 0;
}
double suppl_x_function8_dd ( double x, double y, double z, int main_axis )
{
	return -2.0;
}

double suppl_xyzt_function8 ( double x, double y, double z, double t, double a888, int main_axis )
{
	return C888 * suppl_x_function8 ( x, y, z, main_axis ) * exp(-t) - a888;
}
double suppl_xyzt_function8_dx ( double x, double y, double z, double t, double a888, int main_axis )
{
	if ( main_axis == 0 )
		return C888 * suppl_x_function8_d ( x, y, z, main_axis ) * exp(-t);
	else
		return 0.0;
}
double suppl_xyzt_function8_dy ( double x, double y, double z, double t, double a888, int main_axis )
{
	if ( main_axis == 1 )
		return C888 * suppl_x_function8_d ( x, y, z, main_axis ) * exp(-t);
	else
		return 0.0;
}
double suppl_xyzt_function8_dz ( double x, double y, double z, double t, double a888, int main_axis )
{
	if ( main_axis == 2 )
		return C888 * suppl_x_function8_d ( x, y, z, main_axis ) * exp(-t);
	else
		return 0.0;
}

double suppl_xyzt_function8_dxdx ( double x, double y, double z, double t, double a888, int main_axis )
{
	if ( main_axis == 0 )
		return C888 * suppl_x_function8_dd ( x, y, z, main_axis ) * exp(-t);
	else
		return 0.0;
}
double suppl_xyzt_function8_dydy ( double x, double y, double z, double t, double a888, int main_axis )
{
	if ( main_axis == 1 )
		return C888 * suppl_x_function8_dd ( x, y, z, main_axis ) * exp(-t);
	else
		return 0.0;
}
double suppl_xyzt_function8_dzdz ( double x, double y, double z, double t, double a888, int main_axis )
{
	if ( main_axis == 2 )
		return C888 * suppl_x_function8_dd ( x, y, z, main_axis ) * exp(-t);
	else
		return 0.0;
}
double suppl_xyzt_function8_dt( double x, double y, double z, double t, double a888, int main_axis )
{
	return (-1.0) * C888 * suppl_x_function8 ( x, y, z, main_axis ) * exp(-t);
}

double suppl_function8 (double x, double y, double z, double t, double a888, int main_axis)
{
	return fabs(suppl_xyzt_function8(x,y,z,t,a888,main_axis));
}
double exact_solution8(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	return pow(suppl_function8(x,y,z,t,a888,main_axis), deg888);
}
///////////////////////////////////////////////////////////////////////
double exact_gradient8x(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function8_dx (x,y,z,t,a888,main_axis);
	return res;
}
double exact_gradient8y(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function8_dy (x,y,z,t,a888,main_axis);
	return res;
}
double exact_gradient8z(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
	{
		//printf ( "val = %e \n", suppl_xyzt_function8(x, y, z, t, a888, main_axis)  );
		res *= -1.0;
	}
	res *= suppl_xyzt_function8_dz (x,y,z,t,a888,main_axis);
	return res;
}

///////////////////////////////////////////////////////////////////////
double exact_gradient8xx(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res1 = deg888 * (deg888 - 1) * exact_solution8(x, y, z, t, a888, deg888 - 2, main_axis);
	res1 *= pow(suppl_xyzt_function8_dx (x,y,z,t,a888, main_axis), 2);

	double res2 = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
		res2 *= -1.0;
	res2 *= suppl_xyzt_function8_dxdx (x,y,z,t,a888, main_axis);

	return res1 + res2;
}
double exact_gradient8yy(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res1 = deg888 * (deg888 - 1) * exact_solution8(x, y, z, t, a888, deg888 - 2, main_axis);
	res1 *= pow(suppl_xyzt_function8_dy (x,y,z,t,a888, main_axis), 2);

	double res2 = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
		res2 *= -1.0;
	res2 *= suppl_xyzt_function8_dydy (x,y,z,t,a888, main_axis);

	return res1 + res2;
}
double exact_gradient8zz(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res1 = deg888 * (deg888 - 1) * exact_solution8(x, y, z, t, a888, deg888 - 2, main_axis);
	res1 *= pow(suppl_xyzt_function8_dz (x,y,z,t,a888, main_axis), 2);

	double res2 = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
		res2 *= -1.0;
	res2 *= suppl_xyzt_function8_dzdz (x,y,z,t,a888, main_axis);

	return res1 + res2;
}

///////////////////////////////////////////////////////////////////////
double exact_laplace8(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	return exact_gradient8xx (x, y, z, t, a888, deg888, main_axis ) + exact_gradient8yy (x, y, z, t, a888, deg888, main_axis ) + exact_gradient8zz (x, y, z, t, a888, deg888, main_axis );
}
double exact_derivdt8(double x, double y, double z, double t, double a888, double deg888, int main_axis)
{
	double res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
	if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
		res *= -1.0;
	res *= suppl_xyzt_function8_dt (x,y,z,t,a888, main_axis);

	return res;
}
double righthand8(int main_axis, double a888, double deg888, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempv, tempTw;//, tempbound, tempdivvx, tempdivvy, tempdivvz;

	return cell_volume *
		( density(i,j,k) * heat_capacity(i,j,k) * exact_derivdt8(x_mid, y_mid, z_mid, t_discr*tau, a888, deg888, main_axis) +
		(-1) * heat_conductivity(i,j,k) * exact_laplace8(x_mid, y_mid, z_mid, t_discr*tau, a888, deg888, main_axis) );
}
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
double exact_solution889(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	switch(main_axis)
	{
	case 0:
		return pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * sin ( my * PI * y ) * sin ( mz * PI * z );
		break;
	case 1:
		return sin ( mx * PI * x ) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * sin ( mz * PI * z );
		break;
	case 2:
		return sin ( mx * PI * x ) * sin ( my * PI * y ) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888);
		break;
	}
	return 0.0;
}
///////////////////////////////////////////////////////////////////////
double exact_gradient889x(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res = 0.0;

	switch (main_axis)
	{
	case 0:
		res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res *= -1.0;
		res *= suppl_xyzt_function8_dx (x,y,z,t,a888, main_axis);
		res *= sin ( my * PI * y ) * sin ( mz * PI * z );
		break;
	case 1:
		res = (mx * PI) * cos (mx * PI * x) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * sin ( mz * PI * z );
		break;
	case 2:
		res = (mx * PI) * cos (mx * PI * x) * sin ( my * PI * y ) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888);
		break;
	}
	return res;
}
double exact_gradient889y(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res = 0.0;

	switch (main_axis)
	{
	case 0:
		res = pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * (my * PI) * cos (my * PI * y) * sin ( mz * PI * z );
		break;
	case 1:
		res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res *= -1.0;
		res *= suppl_xyzt_function8_dy (x,y,z,t,a888, main_axis);
		res *= sin ( mx * PI * x ) * sin ( mz * PI * z );
		break;
	case 2:
		res = sin ( mx * PI * x ) * (my * PI) * cos (my * PI * y) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888);
		break;
	}

	return res;
}
double exact_gradient889z(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res = 0.0;

	switch (main_axis)
	{
	case 0:
		res =  pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * sin ( my * PI * y ) * (mz * PI) * cos (mz * PI * z);
		break;
	case 1:
		res = sin ( mx * PI * x ) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * (mz * PI) * cos (mz * PI * z);
		break;
	case 2:
		res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res *= -1.0;
		res *= suppl_xyzt_function8_dz (x,y,z,t,a888, main_axis);
		res *= sin ( mx * PI * x ) * sin ( my * PI * y );
		break;
	}

	return res;
}
/////////////////////////////////////////////////////////////////////
double exact_gradient889xx(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res1, res2;
	switch (main_axis)
	{
	case 0:
		res1 = deg888 * (deg888 - 1) * exact_solution8(x, y, z, t, a888, deg888 - 2, main_axis);
		res1 *= pow(suppl_xyzt_function8_dx (x,y,z,t,a888, main_axis), 2);

		res2 = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res2 *= -1.0;
		res2 *= suppl_xyzt_function8_dxdx (x,y,z,t,a888, main_axis);
		return (res1 + res2) * sin ( my * PI * y ) * sin ( mz * PI * z );
		break;
	case 1:
		return (-1) * (mx * PI) * (mx * PI) * sin (mx * PI * x) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * sin ( mz * PI * z );
		break;
	case 2:
		return (-1) * (mx * PI) * (mx * PI) * sin (mx * PI * x) * sin ( my * PI * y ) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888);
		break;
	}

	return 0.0;
}
double exact_gradient889yy(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res1, res2;

	switch(main_axis)
	{
	case 0:
		return pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * (-1) * (my * PI) * (my * PI) * sin (my * PI * y) * sin ( mz * PI * z );
		break;
	case 1:
		res1 = deg888 * (deg888 - 1) * exact_solution8(x, y, z, t, a888, deg888 - 2, main_axis);
		res1 *= pow(suppl_xyzt_function8_dy (x,y,z,t,a888, main_axis), 2);

		res2 = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res2 *= -1.0;
		res2 *= suppl_xyzt_function8_dydy (x,y,z,t,a888, main_axis);
		return (res1 + res2) * sin ( mx * PI * x ) * sin ( mz * PI * z );
		break;
	case 2:
		return sin ( mx * PI * x ) * (-1)* (my * PI) * (my * PI) * sin (my * PI * y) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888);
		break;
	}

	return 0.0;
}
double exact_gradient889zz(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res1, res2;

	switch(main_axis)
	{
	case 0:
		return pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * sin ( my * PI * y ) * (-1) * (mz * PI) * (mz * PI) * sin (mz * PI * z);
		break;
	case 1:
		return sin ( mx * PI * x ) * pow(suppl_function8(x,y,z,t,a888,main_axis), deg888) * (-1) * (mz * PI) * (mz * PI) * sin (mz * PI * z);
		break;
	case 2:
		res1 = deg888 * (deg888 - 1) * exact_solution8(x, y, z, t, a888, deg888 - 2, main_axis);
		res1 *= pow(suppl_xyzt_function8_dz (x,y,z,t,a888, main_axis), 2);

		res2 = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res2 *= -1.0;
		res2 *= suppl_xyzt_function8_dzdz (x,y,z,t,a888, main_axis);
		return (res1 + res2) * sin ( mx * PI * x ) * sin (my * PI * y);
		break;
	}

	return 0.0;
}

///////////////////////////////////////////////////////////////////////
double exact_laplace889(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	return exact_gradient889xx (x, y, z, t, a888, deg888, mx, my, mz, main_axis ) 
		+ exact_gradient889yy (x, y, z, t, a888, deg888, mx, my, mz, main_axis )
		+ exact_gradient889zz (x, y, z, t, a888, deg888, mx, my, mz, main_axis );
}
double exact_derivdt889(double x, double y, double z, double t, double a888, double deg888, int mx, int my, int mz, int main_axis)
{
	double res;
	
	switch (main_axis)
	{
	case 0:
		res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res *= -1.0;
		res *= suppl_xyzt_function8_dt (x,y,z,t,a888, main_axis);
		res *= sin ( my * PI * y ) * sin ( mz * PI * z );
		break;
	case 1:
		res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res *= -1.0;
		res *= suppl_xyzt_function8_dt (x,y,z,t,a888, main_axis);
		res *= sin ( mx * PI * x ) * sin ( mz * PI * z );
		break;
	case 2:
		res = deg888 * exact_solution8(x, y, z, t, a888, deg888 - 1, main_axis);
		if ( suppl_xyzt_function8(x, y, z, t, a888, main_axis) < 0.0 ) 
			res *= -1.0;
		res *= suppl_xyzt_function8_dt (x,y,z,t,a888, main_axis);
		res *= sin ( mx * PI * x ) * sin ( my * PI * y );
		break;
	}

	return res;
}
double righthand889(int main_axis, int mx, int my, int mz, double a888, double deg888, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempv, tempTw;//, tempbound, tempdivvx, tempdivvy, tempdivvz;

	return cell_volume *
		( density(i,j,k) * heat_capacity(i,j,k) * exact_derivdt889(x_mid, y_mid, z_mid, t_discr*tau, a888, deg888, mx, my, mz, main_axis) +
		(-1) * heat_conductivity(i,j,k) * exact_laplace889(x_mid, y_mid, z_mid, t_discr*tau, a888, deg888, mx, my, mz, main_axis) ); 
}

///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
double suppl_func4_s (double s, double ms, double deg4 )
{
	if ( s < 0.25 )
		return 0;
	else if ( s < 0.75 )
		return pow(cos(ms * PI * s), deg4 );
	else
		return 0.0;
}
double suppl_func4_s1deriv (double s, double ms, double deg4 )
{
	if ( s < 0.25 )
		return 0;
	else if ( s < 0.75 )
		return deg4 * pow(cos(ms * PI * s), deg4 - 1 ) * (-ms * PI) * sin (ms * PI * s);
	else
		return 0.0;
}
double suppl_func4_s2deriv (double s, double ms, double deg4 )
{
	if ( s < 0.25 )
		return 0;
	else if ( s < 0.75 )
	{
		double res = deg4 * (deg4 - 1) * pow(cos(ms * PI * s), deg4 - 2 ) * (-ms * PI) * sin (ms * PI * s) * (-ms * PI) * sin (ms * PI * s);
		res += deg4 * pow(cos(ms * PI * s), deg4 - 1 ) * (-ms * PI) * (ms * PI) * cos (ms * PI * s);
		return res;
	}
	else
		return 0.0;
}
// mx,my,mz must be equal to 2k, k>=1
double exact_solution444(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s(x,mx,deg444)*suppl_func4_s(y,my,deg444)*suppl_func4_s(z,mz,deg444)*exp(-t);
}
double exact_gradient444x(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s1deriv(x,mx,deg444)*suppl_func4_s(y,my,deg444)*suppl_func4_s(z,mz,deg444)*exp(-t);
}

double exact_gradient444y(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s(x,mx,deg444)*suppl_func4_s1deriv(y,my,deg444)*suppl_func4_s(z,mz,deg444)*exp(-t);
}

double exact_gradient444z(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s(x,mx,deg444)*suppl_func4_s(y,my,deg444)*suppl_func4_s1deriv(z,mz,deg444)*exp(-t);
}

double exact_gradient444xx(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s2deriv(x,mx,deg444)*suppl_func4_s(y,my,deg444)*suppl_func4_s(z,mz,deg444)*exp(-t);
}
double exact_gradient444yy(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s(x,mx,deg444)*suppl_func4_s2deriv(y,my,deg444)*suppl_func4_s(z,mz,deg444)*exp(-t);
}
double exact_gradient444zz(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return nonstationary*suppl_func4_s(x,mx,deg444)*suppl_func4_s(y,my,deg444)*suppl_func4_s2deriv(z,mz,deg444)*exp(-t);
}
double exact_laplace444(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return exact_gradient444xx(x,y,z,t,mx,my,mz,deg444) + exact_gradient444yy(x,y,z,t,mx,my,mz,deg444) + exact_gradient444zz(x,y,z,t,mx,my,mz,deg444);
}
double exact_derivdt444(double x, double y, double z, double t, double mx, double my, double mz, double deg444)
{
	return -exact_solution444(x,y,z,t,mx,my,mz,deg444);
}
double righthand444(int mx, int my, int mz, double deg444, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	double tempvx, tempvy, tempvz, tempv, tempTw;//, tempbound, tempdivvx, tempdivvy, tempdivvz;

	return cell_volume *
		( density(i,j,k) * heat_capacity(i,j,k) * exact_derivdt444(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mz, deg444) +
		(-1) * heat_conductivity(i,j,k) * exact_laplace444(x_mid, y_mid, z_mid, t_discr*tau,  mx, my, mz, deg444) ); 
}