#include "math.h"

#define STEKLOV_H

#define A (1.0)
#define B (1.0)

typedef enum{
  eDirichlet = 0,
  eNeumann,
  eMixed
} boundaryConditions;


//! the initial data
boundaryConditions xCondition, yCondition, zCondition;

extern double tau;
extern double *xpoints, *ypoints, *zpoints;

const double PI = 3.14159265358979323;
const int nonstationary = 1;

const int velocityX_num = 0;//номер скорости Vx для exact_solution1111
const int velocityY_num = 0;//номер скорости Vy для exact_solution1111
const int velocityZ_num = 0;//номер скорости Vz для exact_solution1111
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
double vx_const = 0.0;  //значение для постоянной скорости Vx
double vy_const = 0.0;  //значение для постоянной скорости Vy
double vz_const = 0.0;  //значение для постоянной скорости Vz

const int Dir_axis = 1;
//const int Neu_axis = 1;
const double mx = 2; //nomer garmoniki x
const double my = 2; //nomer garmoniki y
const double mz = 0; //nomer garmoniki z
const double mt = 0;
const double alpha = 4.0;

double Tx_0_bound(int numSolution, int i, int j, double time);
double Tx_1_bound(int numSolution, int i, int j, double time);
double Ty_0_bound(int numSolution, int i, int j, double time);
double Ty_1_bound(int numSolution, int i, int j, double time);
double wx_0_bound(int numSolution, int i, int j, double time);
double wx_1_bound(int numSolution, int i, int j, double time);
double wy_0_bound(int numSolution, int i, int j, double time);
double wy_1_bound(int numSolution, int i, int j, double time);

double exact_solution9006(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_gradient9006x(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient9006y(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient9006z(double x, double y, double z, double t, double mx, double my, double mz); 
double righthand9006(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution9066(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_gradient9066x(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient9066y(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient9066z(double x, double y, double z, double t, double mx, double my, double mz); 
double righthand9066(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1111(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_solution1111Vx(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_solution1111Vy(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_solution1111Vz(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_gradient1111x(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1111y(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1111z(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_solution1111_po_t(double x, double y, double z, double t, double mx, double my, double mz); //неоднородные условия Неймана y, xy
double exact_gradient1111x_po_t(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1111y_po_t(double x, double y, double z, double t, double mx, double my, double mz); 
double righthand1111(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1166(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); //неоднородные условия Неймана y, xy
double exact_gradient1166x(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double exact_gradient1166y(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double exact_gradient1166z(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double exact_solution1166_po_t(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); //неоднородные условия Неймана y, xy
double exact_gradient1166x_po_t(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double exact_gradient1166y_po_t(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis); 
double righthand1166(int dir_axis, double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution1199(double x, double y, double z, double t, double mx, double my, double mz);
double exact_laplace1199(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_solution1199_po_t(double x, double y, double z, double t, double mx, double my, double mz);
double exact_gradient1199x(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1199y(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1199z(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1199x_po_t(double x, double y, double z, double t, double mx, double my, double mz); 
double exact_gradient1199y_po_t(double x, double y, double z, double t, double mx, double my, double mz); 
double righthand1199(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution206(double x, double y, double z, double t, double mx, double my); //неоднородные условия Неймана y, xy
double exact_gradient206x(double x, double y, double z, double t, double mx, double my); 
double exact_gradient206y(double x, double y, double z, double t, double mx, double my); 
double exact_solution206_po_t(double x, double y, double z, double t, double mx, double my); //неоднородные условия Неймана y, xy
double exact_gradient206x_po_t(double x, double y, double z, double t, double mx, double my); 
double exact_gradient206y_po_t(double x, double y, double z, double t, double mx, double my); 
double righthand206(double mx, double my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution267(double x, double y, double z, double t, int m, int n); //неоднородные условия Неймана y, xy
double exact_gradient267x(double x, double y, double z, double t, int m, int n); 
double exact_gradient267y(double x, double y, double z, double t, int m, int n); 
double righthand267(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution500(double x, double y, double z, double t, int m, int n);
double exact_gradient500x(double x, double y, double z, double t, int m, int n);
double exact_gradient500y(double x, double y, double z, double t, int m, int n);
double exact_solution500_po_t(double x, double y, double z, double t, int m, int n);
double exact_laplace500(double x, double y, double z, double t, int m, int n);
double righthand500(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution700(double x, double y, double z, double t, int m, int n);
double exact_gradient700x(double x, double y, double z, double t, int m, int n);
double exact_gradient700y(double x, double y, double z, double t, int m, int n);
double exact_solution700_po_t(double x, double y, double z, double t, int m, int n);
double exact_laplace700(double x, double y, double z, double t, int m, int n);
double righthand700(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution701(double x, double y, double z, double t, int m, int n);
double exact_gradient701x(double x, double y, double z, double t, int m, int n);
double exact_gradient701y(double x, double y, double z, double t, int m, int n);
double exact_solution701_po_t(double x, double y, double z, double t, int m, int n);
double exact_laplace701 (double x, double y, double z, double t, int m, int n);
double righthand701(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double exact_solution702(double x, double y, double z, double t, int m, int n);
double exact_gradient702x(double x, double y, double z, double t, int m, int n);
double exact_gradient702y(double x, double y, double z, double t, int m, int n);
double exact_solution702_po_t(double x, double y, double z, double t, int m, int n);
double exact_laplace702 (double x, double y, double z, double t, int m, int n);
double righthand702(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);


double exact_solution300(double x, double y, double z, double t, int m, int n); 
double exact_solution300_onlytime(double x, double y, double z, double t, int m, int n);
double exact_gradient300x(double x, double y, double z, double t, int m, int n); 
double exact_gradient300y(double x, double y, double z, double t, int m, int n); 
double exact_gradient_average_300y(int i, int j, int t_discr, double tau, int m, int n);
double exact_gradient_average_300x(int i, int j, int t_discr, double tau, int m, int n);
double exact_laplace300(double x, double y, double z, double t, int m, int n); 
double righthand300(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double function_g(double y, int n);
double function_g_gradient(double y, int n);
double function_g_2derivative(double y, int n);
double exact_average_solution300(int i, int j, int t_discr, double tau, int m, int n);
double exact_average_solution_other300(int i, int j, int t_discr, double tau, int m, int n);

double exact_solution600(double x, double y, double z, double t, int m, int n); 
double exact_solution600_onlytime(double x, double y, double z, double t, int m, int n);
double exact_gradient600x(double x, double y, double z, double t, int m, int n); 
double exact_gradient600y(double x, double y, double z, double t, int m, int n); 
//double exact_gradient_average_300y(int i, int j, int t_discr, double tau, int m, int n);
//double exact_gradient_average_300x(int i, int j, int t_discr, double tau, int m, int n);
double exact_laplace600(double x, double y, double z, double t, int m, int n); 
double righthand600(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double function_gdirneu(double y, int n);
double function_gdirneu_gradient(double y, int n);
double function_gdirneu_2derivative(double y, int n);
double exact_average_solution600(int i, int j, int t_discr, double tau, int m, int n);
double exact_average_solution_other600(int i, int j, int t_discr, double tau, int m, int n);

double exact_solution400(double x, double y, double z, double t, int m, int n); 
double exact_solution400_onlytime(double x, double y, double z, double t, int m, int n);
double exact_gradient400x(double x, double y, double z, double t, int m, int n); 
double exact_gradient400y(double x, double y, double z, double t, int m, int n); 
double exact_laplace400(double x, double y, double z, double t, int m, int n); 
double righthand400(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double function_glambda(double y, int n);
double function_glambda_gradient(double y, int n);
double function_glambda_2derivative(double y, int n);
double exact_average_solution_other400(int i, int j, int t_discr, double tau, int m, int n);


double exact_average_solution_other105(int i, int j, int t_discr, double tau, int m, int n);

double arbSmooth_ex1_solution(double x, double y, double z, double t, double mx, double my, double mt); //неоднородные условия Неймана y, xy
double arbSmooth_ex1_gradientX(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex1_gradientY(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex1_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt); //неоднородные условия Неймана y, xy
double arbSmooth_ex1_gradientX_po_t(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex1_gradientY_po_t(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex1_righthand(double mx, double my, double mt, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);
double arbSmooth_average_ex1_other_solution(int i, int j, int t_discr, double tau, int mx, int my, int mt);
double arbSmooth_ex1_gradientX_integ(int i, int j, double t, double mx, double my, double mt);

double arbSmooth_ex2_solution(double x, double y, double z, double t, double mx, double my, double mt); //неоднородные условия Неймана y, xy
double arbSmooth_ex2_gradientX(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex2_gradientY(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex2_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt); //неоднородные условия Неймана y, xy
double arbSmooth_ex2_gradientX_po_t(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex2_gradientY_po_t(double x, double y, double z, double t,  double mx, double my, double mt); 
double arbSmooth_ex2_righthand(double mx, double my, double mt, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double arbNSmooth_ex3_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex3_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex3_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex3_gradientX(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex3_gradientY(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex3_gradientX_po_t(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex3_gradientY_po_t(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex3_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double arbNSmooth_ex4_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex4_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex4_gradientX(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex4_gradientY(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex4_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double arbNSmooth_ex5_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex5_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex5_gradientX(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex5_gradientY(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex5_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double arbNSmooth_ex6_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex6_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex6_gradientX(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex6_gradientY(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex6_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double arbNSmooth_ex7_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex7_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex7_gradientX(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex7_gradientY(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex7_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);

double arbNSmooth_ex8_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex8_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt, double alpha); //неоднородные условия Неймана y, xy
double arbNSmooth_ex8_gradientX(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex8_gradientY(double x, double y, double z, double t,  double mx, double my, double mt, double alpha); 
double arbNSmooth_ex8_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints);


double a_gradientX(double x, double y, double z);
double a_gradientY(double x, double y, double z);


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
  case 2:
      xCondition = eNeumann;
      yCondition = eNeumann;
      zCondition = eDirichlet;
    break;
  case 7:
      xCondition = eDirichlet;
      yCondition = eNeumann;
      zCondition = eNeumann;
    break;
  case 101:
      xCondition = eNeumann;
      yCondition = eNeumann;
      zCondition = eDirichlet;
    break;
  case 102:
      xCondition = eNeumann;
      yCondition = eDirichlet;
      zCondition = eNeumann;
    break;
  case 103:
      xCondition = eDirichlet;
      yCondition = eNeumann;
      zCondition = eNeumann;
     break;
  case 104:
    xCondition = eNeumann;
    yCondition = eNeumann;
    zCondition = eNeumann;
    break;
  case 105:
    switch(Dir_axis)
    {
    case 0:
      xCondition = eDirichlet;
      yCondition = eNeumann;
      zCondition = eNeumann;
      break;
    case 1:
      xCondition = eNeumann;
      yCondition = eDirichlet;
      zCondition = eNeumann;
      break;
    case 2:
      xCondition = eNeumann;
      yCondition = eNeumann;
      zCondition = eDirichlet;
      break;
    }
    break;
  case 305:
	  if (Dir_axis == 1)
	  {
		  xCondition = eNeumann;
		  yCondition = eDirichlet;
		  zCondition = eNeumann;
	  }
	  else
		  printf("wrong Dir_axis for solutioin 305 (267) (setBoundaryConditions) \n");
	  break;
  case 106:
    xCondition = eNeumann;
    yCondition = eNeumann;
    zCondition = eMixed;
    break;
  case 107:
    xCondition = eNeumann;
    yCondition = eNeumann;
    zCondition = eNeumann;
    break;
  case 108:
    xCondition = eNeumann;
    //yCondition = eNeumann;
    //xCondition = eDirichlet;
    yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 109:
    //xCondition = eNeumann;
    //yCondition = eNeumann;
    xCondition = eDirichlet;
    yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 110:
    //xCondition = eNeumann;
    //yCondition = eNeumann;
    xCondition = eDirichlet;
    yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 111:
    //xCondition = eNeumann;
    //yCondition = eNeumann;
    xCondition = eDirichlet;
    yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 112:
    xCondition = eNeumann;
    yCondition = eNeumann;
    //xCondition = eDirichlet;
    //yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 113:
    //xCondition = eNeumann;
    //yCondition = eNeumann;
    xCondition = eDirichlet;
    yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 114:
    //xCondition = eNeumann;
    //yCondition = eNeumann;
    xCondition = eDirichlet;
    yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 115:
    xCondition = eNeumann;
    yCondition = eNeumann;
    //xCondition = eDirichlet;
    //yCondition = eDirichlet;
    zCondition = eNeumann;
    break;
  case 300:
	  xCondition = eNeumann;
	  //yCondition = eNeumann; works badly
	  yCondition = eDirichlet;
	  zCondition = eNeumann;
	  break;
  case 600:
	  xCondition = eNeumann;
	  //yCondition = eNeumann;
	  yCondition = eDirichlet;
	  zCondition = eNeumann;
	  break;
  case 400:
	  xCondition = eNeumann;
	  yCondition = eDirichlet;
	  zCondition = eNeumann;
	  break;
  case 500:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eNeumann;
	  break;
  case 700:
	  xCondition = eDirichlet;
	  yCondition = eDirichlet;
	  zCondition = eNeumann;
	  break;
  case 701:
	  xCondition = eDirichlet;
	  yCondition = eNeumann;
	  zCondition = eNeumann;
	  break;
  case 702:
	  xCondition = eNeumann;
	  yCondition = eDirichlet;
	  zCondition = eNeumann;
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
  case 2:
    return 1 - z;
    break;
  case 7:
    return sin ( 2 * PI * x );
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
  case 500:
	  return exact_solution500(x,y,z,t,mx,my);
	  break;
  case 600:
	  return exact_solution600(x,y,z,t,mx,my);
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
  //case 116:
	 // //Dir_axis (0,1,2 -> x,y,z) defines which direction has Dirichlet boundary condition, others have Neumann bnd cnds.
	 // //e.g. Dir_axis = 1 means that y-axis has Dirichlet conditions, x- and z-axis - Neumann.
	 // return exact_solution1170(x,y,z,t,mx,my,mz, Dir_axis);
	 // break;

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
  case 700:
	  return exact_solution700(x,y,z,t,mx,my);
	  break;
  case 701:
	  return exact_solution701(x,y,z,t,mx,my);
	  break;
  case 702:
	  return exact_solution702(x,y,z,t,mx,my);
	  break;
  default:
    break;
  }
  return 0;
}

double exact_gradientX(int number, double x, double y, double z, double t) 
{
  switch(number)
  {
  case 1:
    return 0;
    break;
  case 2:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
    break;
  case 7:
    return 2 * PI * cos ( 2 * PI * x );
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
  case 600:
	  return exact_gradient600x(x,y,z,t,mx,my);
	  break;
  case 400:
	  if (Dir_axis == 1)
		  return exact_gradient400x(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 500:
	  return exact_gradient500x(x,y,z,t,mx,my);
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
  case 700:
	  return exact_gradient700x(x,y,z,t,mx,my);
	  break;
  case 701:
	  return exact_gradient701x(x,y,z,t,mx,my);
	  break;
  case 702:
	  return exact_gradient702x(x,y,z,t,mx,my);
	  break;
  default:
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
  case 2:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
	break;
  case 7:
    return 0.0;
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
  case 600:
	  return exact_gradient600y(x,y,z,t,mx,my);
	  break;
  case 400:
	  if (Dir_axis == 1)
		  return exact_gradient400y(x,y,z,t,mx,my);
	  else
		  return 0.0;
	  break;
  case 500:
	  return exact_gradient500y(x,y,z,t,mx,my);
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
  case 700:
	  return exact_gradient700y(x,y,z,t,mx,my);
	  break;
  case 701:
	  return exact_gradient701y(x,y,z,t,mx,my);
	  break;
  case 702:
	  return exact_gradient702y(x,y,z,t,mx,my);
	  break;
  default:
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
  case 2:
    if (x<=0.5)
      return 0.25*cos(PI*x)*exp(-t);
    else
      return cos(PI*x)*exp(-t);
    break;
  case 7:
    return 0.0;
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
  case 500:
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
  case 700:
    return 0.0;
    break;
  case 701:
    return 0.0;
    break;
  case 702:
    return 0.0;
    break;
  default:
    break;
  }
  return 0;
}

double exact_solution_po_t(int number, double x, double y, double z, double t) 
{
	switch(number)
	{
	case 1:
		return 0;
		break;
	case 2:
		return 0;
		break;
	case 7:
		return 0.0;
		break;

	case 101:
		return -nonstationary*cos(mx*PI*x)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 102:
		return -nonstationary*cos(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 103:
		return -nonstationary*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 104:
		//if mx > 0 - inhomog. Neumann for x, homog. Neumann for x,y; 
		//if mx = 0, my > 0 - inhomog. Neumann for y, Neumann for z (no x, quasi 2D test);
		//if mx = 0, my = 0 - inhomog. Neumann for z (no x,y, quasi 1D test);
		return exact_solution1111_po_t(x, y, z, t, mx, my, mz);
		break;
	case 105:
		//Dir_axis (0,1,2 -> x,y,z) defines which direction has Dirichlet boundary condition, others have Neumann bnd cnds.
		//e.g. Dir_axis = 1 means that y-axis has Dirichlet conditions, x- and z-axis - Neumann.
		return exact_solution1166_po_t(x,y,z,t,mx,my,mz, Dir_axis);
		break;
	case 106:
		//works for certain mx,my,mz, not for all values (special test).
		//i suppose it works for mx = my = 2, mz = 1.
		return exact_solution1199_po_t(x,y,z,t,mx,my,mz);
		break;
	case 107:
		//Neumann test with inhomogeneous Neumann for y (impossible for 1111)
		return exact_solution206_po_t(x,y,z,t,mx,my);
		break;
	case 108:
		//Neumann test with inhomogeneous Neumann for y (impossible for 1111)
		return arbSmooth_ex1_solution_po_t(x,y,z,t,mx,my,mt);
		break;
	case 109:
		//Neumann test with inhomogeneous Neumann for y (impossible for 1111)
		return arbSmooth_ex2_solution_po_t(x,y,z,t,mx,my,mt);
		break;
	case 110:
		//Neumann test with inhomogeneous Neumann for y (impossible for 1111)
		return arbNSmooth_ex3_solution_po_t(x,y,z,t,mx,my,mt,alpha);
		break;
	case 111:
		return 0.0;
		break;
	case 112:
		return 0.0;
		break;
	case 113:
		return 0.0;
		break;
	case 114:
		return 0.0;
		break;
	case 115:
		//Neumann test with inhomogeneous Neumann for y (impossible for 1111)
		return arbNSmooth_ex8_solution_po_t(x,y,z,t,mx,my,mt,alpha);
		break;
	case 500:
		return exact_solution500_po_t(x,y,z,t,mx,my);
		break;
	case 700:
		return 0.0;
		break;
	case 701:
		return 0.0;
		break;
	case 702:
		return 0.0;
		break;
	default:
		break;
	}
	return 0;
}

double exact_gradientX_po_t(int number, double x, double y, double z, double t) 
{
	switch(number)
	{
	case 1:
		return 0;
		break;
	case 2:
		if (x<=0.5)
			return -0.25*cos(PI*x)*exp(-t);
		else
			return -cos(PI*x)*exp(-t);
		break;
	case 7:
		return 0.0;
		break;

	case 101:
		return -nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 102:
		return -nonstationary*(-mx*PI)*sin(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 103:
		return -nonstationary*(mx*PI)*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 104:
		return exact_gradient1111x_po_t(x,y,z,t,mx,my,mz);
		break;
	case 105:
		return exact_gradient1166x_po_t(x,y,z,t,mx,my,mz, Dir_axis);
		break;
	case 106:
		return exact_gradient1199x_po_t(x,y,z,t,mx,my,mz);
		break;
	case 107:
		return exact_gradient206x_po_t(x,y,z,t,mx,my);
		break;
	case 108:
		return arbSmooth_ex1_gradientX_po_t(x,y,z,t,mx,my, mt);
		break;
	case 109:
		return arbSmooth_ex2_gradientX_po_t(x,y,z,t,mx,my, mt);
		break;
	case 110:
		return arbNSmooth_ex3_gradientX_po_t(x,y,z,t,mx,my,mt,alpha);
		break;
	case 111:
		return 0.0;
		break;
	case 112:
		return 0.0;
		break;
	case 113:
		return 0.0;
		break;
	case 114:
		return 0.0;
		break;
	case 115:
		return 0.0;
		break;
	case 700:
		return 0.0;
		break;
	case 701:
		return 0.0;
		break;
	case 702:
		return 0.0;
		break;
	default:
		break;
	}
	return 0;
}

double exact_gradientY_po_t(int number, double x, double y, double z, double t) 
{
	switch(number)
	{
	case 1:
		return 0;
		break;
	case 2:
		if (x<=0.5)
			return -0.25*cos(PI*x)*exp(-t);
		else
			return -cos(PI*x)*exp(-t);
		break;
	case 7:
		return 0.0;
		break;

	case 101:
		return -nonstationary*cos(mx*PI*x)*(-my*PI)*sin(my*PI*y)*sin(mz*PI*z)*exp(-t);
		break;
	case 102:
		return -nonstationary*cos(mx*PI*x)*(my*PI)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 103:
		return -nonstationary*sin(mx*PI*x)*(-my*PI)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 104:
		return exact_gradient1111y_po_t(x,y,z,t,mx,my,mz);
		break;
	case 105:
		return exact_gradient1166y_po_t(x,y,z,t,mx,my,mz,Dir_axis);
		break;
	case 106:
		return exact_gradient1199y_po_t(x,y,z,t,mx,my,mz);
		break;
	case 107:
		return exact_gradient206y_po_t(x,y,z,t,mx,my);
		break;
	case 108:
		return arbSmooth_ex1_gradientY_po_t(x,y,z,t,mx,my, mt);
		break;
	case 109:
		return arbSmooth_ex2_gradientY_po_t(x,y,z,t,mx,my, mt);
		break;
	case 110:
		return 0.0;
		break;
	case 111:
		return 0.0;
		break;
	case 112:
		return 0.0;
		break;
	case 113:
		return 0.0;
		break;
	case 114:
		return 0.0;
		break;
	case 115:
		return 0.0;
		break;
	case 700:
		return 0.0;
		break;
	case 701:
		return 0.0;
		break;
	case 702:
		return 0.0;
		break;

	default:
		break;
	}
	return 0;
}
double righthand(int number, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
#ifdef SPECIAL_MODE
	return 0;
#endif
  switch(number)
  {
  case 1:
    return 0;
    break;
  case 7:
    printf ( "Should not be called for this numsol \n" );
    return 0.0;
    break;

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
  case 600:
	  return righthand600(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;
  case 400:
	  if (Dir_axis == 1)
		  return righthand400(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  else
		  return 0.0;
	  break;
  case 500:
	  return righthand500(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
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
  case 700:
	  return righthand700(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;
  case 701:
	  return righthand701(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;
  case 702:
	  return righthand702(mx,my,t_discr,i,j,k,xpoints,ypoints, zpoints);
	  break;

  default:
	  break;
  }
  return 0;
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
	//return z*(1.0 - z);
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


double exact_solution1111Vx(double x, double y, double z, double t, double mx, double my, double mz)
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
double exact_solution1111Vy(double x, double y, double z, double t, double mx, double my, double mz)
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
double exact_solution1111Vz(double x, double y, double z, double t, double mx, double my, double mz)
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

double v_init(double x, double y, double z)
{
	return x*(1-x)*y*(1-y)*z*(1-z);
}


///////////////////////////////////////////////////////////////////////
double exact_solution1111(double x, double y, double z, double t, double mx, double my, double mz)
{
	if (mx != 0)
		return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (x*x/2)*wx_1/(-heat_conductivity(Nx-1,0,0)) + (x-x*x/2)*wx_0/(-heat_conductivity(0,0,0));
	else if (my != 0)
		return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
	else
		return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t) + (z*z/2)*wz_1/(-heat_conductivity(0,0,Nz-1)) + (z-z*z/2)*wz_0/(-heat_conductivity(0,0,0));
}

double exact_solution1111_po_t(double x, double y, double z, double t, double mx, double my, double mz)
{
	if (mx != 0)
		return -nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
	else if (my != 0)
		return -nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
	else
		return -nonstationary*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
}

double exact_gradient1111x(double x, double y, double z, double t, double mx, double my, double mz)
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
double exact_gradient1111y(double x, double y, double z, double t, double mx, double my, double mz)
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
double exact_gradient1111x_po_t(double x, double y, double z, double t, double mx, double my, double mz)
{
	if (mx > 0)
	{
		double lambdagrad = -nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
	}
	else
	{
		return 0;
	}
}
double exact_gradient1111y_po_t(double x, double y, double z, double t, double mx, double my, double mz)
{
	if (my > 0)
	{
		double lambdagrad = -nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t);
		if (mx > 0)
			return lambdagrad;
		else
		{
			lambdagrad += 0;
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
double exact_gradient1111z(double x, double y, double z, double t, double mx, double my, double mz)
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

double righthand1111(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
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
		//if ((my == 0) && (mz == 0))
		//	tempdivvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*x_mid) * exact_solution1111(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) ;
		//else
		//	tempdivvx = 0;
		//tempvx += tempdivvx;

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
		//if ((mx == 0) && (mz == 0))
		//	tempdivvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * (1 - 2*y_mid) * exact_solution1111(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) ;
		//else
		//	tempdivvy = 0;
		//tempvy += tempdivvy;

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

		//if ((mx == 0) && (my == 0))
		//	tempdivvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * 2.0*(x_mid + y_mid - 1) * exact_solution1111(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) ;
		//else
		//	tempdivvz = 0;
		//tempvz += tempdivvz;

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
	//	printf("t_discr = %d tempv = %f \my",t_discr, tempv);
	//	_getch();
	//}
	//printf("tempvx = %f tempvy = %f tempvz = %f \my",tempvx,tempvy,tempvz);
	//printf("my = %d mz = %d \my", my, mz);
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
double exact_solution1166_po_t(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	switch(dir_axis)
	{
	case 0:
		return -nonstationary*sin(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
		break;
	case 1:
		return -nonstationary*sin(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t);
		break;
	case 2:
		return -nonstationary*sin(mz*PI*z)*cos(mx*PI*x)*cos(my*PI*y)*exp(-t);
		break;
	default:
		printf("FUCK YOU \n");
		return 0;
		break;
	}
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

///////////////////////////////////////////////////////////////////////
double exact_gradient1166x_po_t(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	double lambdagrad = 0;
	switch(dir_axis)
	{
	case 0:
		lambdagrad = -nonstationary*(mx*PI)*cos(mx*PI*x)*cos(my*PI*y)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 1:
		lambdagrad = -nonstationary*(-mx*PI)*sin(mx*PI*x)*sin(my*PI*y)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 2:
		lambdagrad = -nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*sin(mz*PI*z)*exp(-t);
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

double exact_gradient1166y_po_t(double x, double y, double z, double t, double mx, double my, double mz, int dir_axis)
{
	double lambdagrad = 0;
	switch(dir_axis)
	{
	case 0:
		lambdagrad = -nonstationary*(-my*PI)*sin(my*PI*y)*sin(mx*PI*x)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 1:
		lambdagrad = -nonstationary*(my*PI)*cos(my*PI*y)*cos(mx*PI*x)*cos(mz*PI*z)*exp(-t);
		return lambdagrad;
		break;
	case 2:
		lambdagrad = -nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*sin(mz*PI*z)*exp(-t);
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
	double tempvx, tempvy, tempvz, tempbound, tempv, tempTw, tempdivvx, tempdivvy, tempdivvz;

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
double exact_solution1199(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*sin(0.5*PI*z)*exp(-t) + Tz_0 + z*wz_1/(-heat_conductivity(0,0,Nz-1));
}
double exact_solution1199_po_t(double x, double y, double z, double t, double mx, double my, double mz)
{
	return -nonstationary*cos(mx*PI*x)*cos(my*PI*y)*sin(0.5*PI*z)*exp(-t);
}
double exact_gradient1199x(double x, double y, double z, double t, double Harm_xx, double Harm_xy, double mz)
{
	return nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*sin(0.5*PI*z)*exp(-t);
}
double exact_gradient1199y(double x, double y, double z, double t, double Harm_xx, double Harm_xy, double mz)
{
	return nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*sin(0.5*PI*z)*exp(-t);
}
double exact_gradient1199x_po_t(double x, double y, double z, double t, double Harm_xx, double Harm_xy, double mz)
{
	return -nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*sin(0.5*PI*z)*exp(-t);
}
double exact_gradient1199y_po_t(double x, double y, double z, double t, double Harm_xx, double Harm_xy, double mz)
{
	return -nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*sin(0.5*PI*z)*exp(-t);
}
double exact_gradient1199z(double x, double y, double z, double t, double Harm_xx, double Harm_xy, double mz)
{
	return nonstationary*(0.5*PI)*cos(0.5*PI*z)*cos(mx*PI*x)*cos(my*PI*y)*exp(-t) + wz_1/(-heat_conductivity(0,0,Nz-1));
}

double exact_laplace1199(double x, double y, double z, double t, double mx, double my, double mz)
{
	return nonstationary*(-1)*(mx*mx + my*my + 0.5*0.5)*(PI*PI)*cos(mx*PI*x)*cos(my*PI*y)*sin(0.5*PI*z)*exp(-t);
}

double exact_solution_po_t_1199(double x, double y, double z, double t, double mx, double my, double mz)
{
	return (-1)*nonstationary*cos(mx*PI*x)*cos(my*PI*y)*sin(0.5*PI*z)*exp(-t);
}

double righthand1199(double mx, double my, double mz, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double cell_volume = hx(i)*hy(j)*hz(k);
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);

	double tempvx, tempvy, tempvz, tempbound, tempv, tempTw, tempdivvx, tempdivvy, tempdivvz;

	tempTw = cell_volume*nonstationary*(density(i,j,k)*heat_capacity(i,j,k)*exact_solution_po_t_1199(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mz));
	tempTw += cell_volume*nonstationary*(-heat_conductivity(i,j,k)*exact_laplace1199(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mz));

	tempvx = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vx(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1199x(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
	tempvy = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vy(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1199y(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);
	tempvz = cell_volume * density(i,j,k)*heat_capacity(i,j,k) * exact_solution1111Vz(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz) * exact_gradient1199z(x_mid,y_mid,z_mid,t_discr*tau,mx,my,mz);

	//double tempbound = (0.5*(xpoints[i+1]*xpoints[i+1] - xpoints[i]*xpoints[i])*( (wx_1/(-heat_conductivity(Nx-1,j,k)))  - (wx_0/(-heat_conductivity(i,0,k))) ) + (xpoints[i+1] - xpoints[i])*wx_0/(-heat_conductivity(i,0,k)))*hy(j)*hz(k);
	//double tempvx = nonstationary*(mx*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau) + tempbound;
	//double tempvy = nonstationary*(my*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) ) * ( sin(mz*PI*zpoints[k+1]) - sin(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
	//double tempvz = nonstationary*(mz*PI)*density(i,j,k)*heat_capacity(i,j,k)*1.0/((mx*PI)*(my*PI)*(mz*PI)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * ( cos(mz*PI*zpoints[k+1]) - cos(mz*PI*zpoints[k]) ) * exp(-t_discr*tau);
	//double tempv = vx_const*tempvx + vy_const*tempvy + vz_const*tempvz;
	tempv = tempvx + tempvy + tempvz;
	//if (((t_discr > 950) && (t_discr % 20 == 0))&&(((i==4)&&(j==4))&&(k==4)))
	//{
	//	printf("t_discr = %d tempv = %f \my",t_discr, tempv);
	//	_getch();
	//}
	//printf("tempvx = %f tempvy = %f tempvz = %f \my",tempvx,tempvy,tempvz);
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

///////////////////////////////////////////////////////////////////////
double exact_solution206(double x, double y, double z, double t, double mx, double my)
{
  return nonstationary*cos(mx*PI*x)*cos(my*PI*y)*exp(-t) + (y*y/2)*wy_1/(-heat_conductivity(0,Ny-1,0)) + (y-y*y/2)*wy_0/(-heat_conductivity(0,0,0));
}
double exact_solution206_po_t(double x, double y, double z, double t, double mx, double my)
{
	return -nonstationary*cos(mx*PI*x)*cos(my*PI*y)*exp(-t);
}
double exact_gradient206x(double x, double y, double z, double t, double mx, double my)
{
  return nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*exp(-t); 
}
double exact_gradient206y(double x, double y, double z, double t, double mx, double my)
{
  return nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*exp(-t) + y*wy_1/(-heat_conductivity(0,Ny-1,0)) + (1-y)*wy_0/(-heat_conductivity(0,0,0));
}
double exact_gradient206x_po_t(double x, double y, double z, double t, double mx, double my)
{
	return -nonstationary*(-mx*PI)*sin(mx*PI*x)*cos(my*PI*y)*exp(-t); 
}
double exact_gradient206y_po_t(double x, double y, double z, double t, double mx, double my)
{
	return -nonstationary*(-my*PI)*sin(my*PI*y)*cos(mx*PI*x)*exp(-t);
}

double righthand206(double mx, double my, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double temp = 0;
  if (mx != 0)
    temp = nonstationary*(heat_conductivity(i,j,k)*(mx*mx + my*my)/(mx*my) - density(i,j,k)*heat_capacity(i,j,k)/((mx*PI)*(my*PI)) ) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
  else
    temp = nonstationary*(heat_conductivity(i,j,k)*my*PI - density(i,j,k)*heat_capacity(i,j,k)/(my*PI) ) * ( sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]) ) * hx(i) * exp(-t_discr*tau) - heat_conductivity(i,j,k)*hx(i)*hy(j)*(wy_1/(-heat_conductivity(i,Ny-1,k))  - wy_0/(-heat_conductivity(i,0,k)));
  return hz(k)*temp;
  //	return temp;
}
///////////////////////////////////////////////////////
double exact_solution267(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary*sin(n*PI*y)*cos(m*PI*x)*exp(-t) + Ty_0 + 2*y*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*sin(n*PI*y)*cos(m*PI*x)*exp(-t) + (heat_conductivity(0,0,0)*Ty_0 + heat_conductivity(0,Ny-1,0)*Ty_1)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0)) + 2*(y-0.5)*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient267y(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
	{
		return (n*PI)*nonstationary*cos(n*PI*y)*cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
		//		printf("heat_conductivity(0,Ny-1,0) = %f \n",heat_conductivity(0,Ny-1,0));
	}
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*(n*PI)*nonstationary*cos(n*PI*y)*cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient267x(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary*(-m*PI)*sin(n*PI*y)*sin(m*PI*x)*exp(-t); 
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-m*PI)*sin(n*PI*y)*sin(m*PI*x)*exp(-t); 
}
double righthand267(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp = 0;
	if (m != 0)
	{
		if (ypoints[j]<0.5)
			temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(m*m + n*n)/(m*n) - density(i,j,k)*heat_capacity(i,j,k)/((m*PI)*(n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * ( sin(m*PI*xpoints[i+1]) - sin(m*PI*xpoints[i]) ) * exp(-t_discr*tau) ;
		else
			temp = (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-1)*(heat_conductivity(i,j,k)*(m*m + n*n)/(m*n) - density(i,j,k)*heat_capacity(i,j,k)/((m*PI)*(n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * ( sin(m*PI*xpoints[i+1]) - sin(m*PI*xpoints[i]) ) * exp(-t_discr*tau);
	}
	else
	{
		if (ypoints[j]<0.5)
			temp = nonstationary*(-1)*(heat_conductivity(i,j,k)*(n*n)/(n) - density(i,j,k)*heat_capacity(i,j,k)/((n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * hx(i) * exp(-t_discr*tau) ;
		else
			temp = (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-1)*(heat_conductivity(i,j,k)*(n*n)/(n) - density(i,j,k)*heat_capacity(i,j,k)/((n*PI)) ) * ( cos(n*PI*ypoints[j+1]) - cos(n*PI*ypoints[j]) ) * hx(i) * exp(-t_discr*tau);
	}
	return temp;
}
///////////////////////////////////////////////////////

double exact_solution500(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * cos ( PI * t ) * sin ( x * ( 1 - x ) * y * ( 1 - y ) );
}
double exact_gradient500x(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * cos ( PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) ) * ( 1 - 2 * x ) * y * ( 1 - y );
}
double exact_gradient500y(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * cos ( PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) ) * x * ( 1 - x ) * ( 1 - 2 * y );
}
double exact_solution500_po_t(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * (-PI) * sin ( PI * t ) * sin ( x * ( 1 - x ) * y * ( 1 - y ) );
}
double exact_laplace500(double x, double y, double z, double t, int m, int n)
{
	double out = 0.0 ;
	out += nonstationary * cos ( PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) ) * ( - 2 ) * y * ( 1 - y );
	out += nonstationary * cos ( PI * t ) * (-1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) ) * ( 1 - 2 * x ) * (1 - 2 * x ) * y * ( 1 - y ) * y * ( 1 - y );
	out += nonstationary * cos ( PI * t ) * cos ( x * ( 1 - x ) * y * ( 1 - y ) ) * x * ( 1 - x ) * ( - 2 );
	out += nonstationary * cos ( PI * t ) * ( -1) * sin ( x * ( 1 - x ) * y * ( 1 - y ) ) * x * ( 1 - x ) * x * ( 1 - x ) * ( 1 - 2 * y ) * ( 1 - 2 * y ); 

	return out;
}
double righthand500(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * ( density(i,j,k)*heat_capacity(i,j,k)*exact_solution500_po_t(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace500(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}

double exact_solution700(double x, double y, double z, double t, int m, int n)
{
	return x * (1 - x) * y * (1 - y);
}
double exact_gradient700x(double x, double y, double z, double t, int m, int n)
{
	return (1 - 2.0 * x) * y * (1 - y);
}
double exact_gradient700y(double x, double y, double z, double t, int m, int n)
{
	return (1 - 2.0 * y) * x * (1 - x);
}
double exact_solution700_po_t(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double exact_laplace700(double x, double y, double z, double t, int m, int n)
{
	return (-2.0) * y * (1 - y) + (-2.0) * x * (1 - x);
}
double righthand700(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * ( density(i,j,k)*heat_capacity(i,j,k)*exact_solution700_po_t(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace700(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}

double exact_solution701(double x, double y, double z, double t, int m, int n)
{
	return A + B * x;
}
double exact_gradient701x(double x, double y, double z, double t, int m, int n)
{
	return B;
}
double exact_gradient701y(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double exact_solution701_po_t(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double exact_laplace701(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double righthand701(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * ( density(i,j,k)*heat_capacity(i,j,k)*exact_solution701_po_t(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace701(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}

double exact_solution702(double x, double y, double z, double t, int m, int n)
{
	return A + B * y;
}
double exact_gradient702x(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double exact_gradient702y(double x, double y, double z, double t, int m, int n)
{
	return B;
}
double exact_solution702_po_t(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double exact_laplace702(double x, double y, double z, double t, int m, int n)
{
	return 0.0;
}
double righthand702(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * ( density(i,j,k)*heat_capacity(i,j,k)*exact_solution702_po_t(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace702(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}
///////////////////////////////////////////////////////
double function_g(double y, int n)
{
	return y * sin(PI*n*y);
}

double exact_solution300_onlytime(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary * function_g(y,n) * cos(m*PI*x)*exp(-t);
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * function_g(y,n) * cos(m*PI*x)*exp(-t);
}
double exact_solution300(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary * function_g(y,n) * cos(m*PI*x)*exp(-t) + Ty_0 + 2*y*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * function_g(y,n) * cos(m*PI*x)*exp(-t) + (heat_conductivity(0,0,0)*Ty_0 + heat_conductivity(0,Ny-1,0)*Ty_1)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0)) + 2*(y-0.5)*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_average_solution300(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	if (n == 0)
		integ_y = hy(j);
	else
		integ_y = -(1.0 / (PI * n))*(ypoints[j+1]*cos(PI * n * ypoints[j+1]) - ypoints[j]*cos(PI * n * ypoints[j])) + (1.0 / (PI * n * PI * n))*(sin(PI * n * ypoints[j+1]) - sin(PI * n * ypoints[j]));
	if (m == 0)
		integ_x = hx(i);
	else
		integ_x = (1.0 / (PI * m)) * (sin(PI * m * xpoints[i+1]) - sin(PI * m * xpoints[i]));

	if (ypoints[j] + 0.5 * hy(j) <= 0.5)
		return nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
}

double integral_x(double x, int m)
{
	if (m == 0)
		return x;
	else
		return (1.0 / (PI * m)) * sin (PI * m * x);
}
double integral_y(double y, int n)
{
	if (n == 0)
		return y;
	else
		return - (1.0 / (PI * n)) * y * cos (PI * n * y) + (1.0 / (PI * n * PI * n)) * sin (PI * n * y);
}
double exact_average_solution_other300(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	integ_x = integral_x(xpoints[i+1],m) - integral_x(xpoints[i],m);
	integ_y = integral_y(ypoints[j+1],n) - integral_y(ypoints[j],n);

	if (ypoints[j] + 0.5 * hy(j) <= 0.5)
		return nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
}

double function_g_gradient(double y, int n)
{
	return sin(PI*n*y) + y * (PI * n) * cos(PI*n*y);
}

double exact_gradient300y(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
	{
		return nonstationary * function_g_gradient(y,n) *cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
	}
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0)) * nonstationary * function_g_gradient(y,n) *cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient300x(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary*(-m*PI)*function_g(y,n)*sin(m*PI*x)*exp(-t); 
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-m*PI)*function_g(y,n)*sin(m*PI*x)*exp(-t); 
}


double exact_gradient_average_300y(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	integ_x = integral_x(xpoints[i+1],m) - integral_x(xpoints[i],m);

	if (ypoints[j] <= 0.5)
		return nonstationary * exp (-t_discr * tau) * function_g_gradient(ypoints[j],n) * integ_x / hx(i);
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * exp (-t_discr * tau) * function_g_gradient(ypoints[j],n) * integ_x / hx(i);
}

double exact_gradient_average_300x(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	integ_y = integral_y(ypoints[j+1],n) - integral_y(ypoints[j],n);

	if (ypoints[j] <= 0.5)
		return nonstationary * (-m*PI) * exp (-t_discr * tau) * sin(m*PI*(xpoints[i])) * integ_y / hy(j);
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0)) * (-m*PI) * nonstationary * exp (-t_discr * tau) * sin(m*PI*(xpoints[i])) * integ_y / hy(j);
}

double function_g_2derivative(double y, int n)
{
	return 2*PI*n*cos(PI*n*y) - y * (PI * n) * (PI * n) * sin(PI*n*y);
}

double exact_laplace300(double x, double y, double z, double t, int m, int n)
{
	if (y < 0.5)
		return nonstationary * (cos(m*PI*x) * function_g_2derivative(y,n) + (-PI*m*PI*m) * cos(m*PI*x) * function_g(y,n)) * exp(-t);
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))* nonstationary * (cos(m*PI*x) * function_g_2derivative(y,n) + (-PI*m*PI*m) * cos(m*PI*x) * function_g(y,n)) * exp(-t);
}
double righthand300(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * (-density(i,j,k)*heat_capacity(i,j,k)*exact_solution300_onlytime(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace300(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
double function_gdirneu(double y, int n)
{
	if ( n == 0 )
		return 1.0;
	else
		return y * y * (1 - y) * (1 - y);
}
double integral_600x(double x, int m)
{
	if (m == 0)
		return x;
	else
		return (1.0 / (PI * m)) * sin (PI * m * x);
}
double integral_600y(double y, int n)
{
	if (n == 0)
		return y;
	else
		return (1.0/3.0) * y * y * y - (1.0/2.0) * y * y * y * y + (1.0/5.0) * y * y * y * y * y;
}
double exact_solution600_onlytime(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * function_gdirneu(y,n) * cos(m*PI*x)*exp(-t);
}
double exact_solution600(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * function_gdirneu(y,n) * cos(m*PI*x)*exp(-t) + Ty_0 + 2*y*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_average_solution600(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);

	integ_x = integral_600x(xpoints[i+1],m) - integral_600x(xpoints[i],m);
	integ_y = integral_600y(ypoints[j+1],m) - integral_600y(ypoints[j],m);

	return nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
}

double exact_average_solution_other600(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	integ_x = integral_600x(xpoints[i+1],m) - integral_600x(xpoints[i],m);
	integ_y = integral_600y(ypoints[j+1],n) - integral_600y(ypoints[j],n);

	return nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
}

double function_gdirneu_gradient(double y, int n)
{
	if ( n == 0 )
		return 0;
	else
		return 2 * y * (1 - y) * (1 - y) - y * y * 2 * (1 - y) ;
}
double exact_gradient600y(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * function_gdirneu_gradient(y,n) *cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient600x(double x, double y, double z, double t, int m, int n)
{
	return nonstationary*(-m*PI)*function_gdirneu(y,n)*sin(m*PI*x)*exp(-t); 
}

//double exact_gradient_average_300y(int i, int j, int t_discr, double tau, int m, int n)
//{
//	double integ_x, integ_y, volume = hx(i)*hy(j);
//	integ_x = integral_x(xpoints[i+1],m) - integral_x(xpoints[i],m);
//
//	if (ypoints[j] <= 0.5)
//		return nonstationary * exp (-t_discr * tau) * function_g_gradient(ypoints[j],n) * integ_x / hx(i);
//	else
//		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * exp (-t_discr * tau) * function_g_gradient(ypoints[j],n) * integ_x / hx(i);
//}
//
//double exact_gradient_average_300x(int i, int j, int t_discr, double tau, int m, int n)
//{
//	double integ_x, integ_y, volume = hx(i)*hy(j);
//	integ_y = integral_y(ypoints[j+1],n) - integral_y(ypoints[j],n);
//
//	if (ypoints[j] <= 0.5)
//		return nonstationary * (-m*PI) * exp (-t_discr * tau) * sin(m*PI*(xpoints[i])) * integ_y / hy(j);
//	else
//		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0)) * (-m*PI) * nonstationary * exp (-t_discr * tau) * sin(m*PI*(xpoints[i])) * integ_y / hy(j);
//}

double function_gdirneu_2derivative(double y, int n)
{
	if ( n == 0 )
		return 0.0;
	else
		return 2.0 - 12.0 * y + 12.0 * y * y;
}

double exact_laplace600(double x, double y, double z, double t, int m, int n)
{
	return nonstationary * (cos(m*PI*x) * function_gdirneu_2derivative(y,n) + (-PI*m*PI*m) * cos(m*PI*x) * function_gdirneu(y,n)) * exp(-t);
}
double righthand600(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * (-density(i,j,k)*heat_capacity(i,j,k)*exact_solution600_onlytime(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace600(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}
///////////////////////////////////////////////////////
double function_glambda(double y, int n)
{
	return (2.0 * y * y * y - 3.0 * y * y) * sin(PI*n*y);
}
double exact_solution400_onlytime(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary * function_glambda(y,n) * cos(m*PI*x)*exp(-t);
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * function_glambda(y,n) * cos(m*PI*x)*exp(-t);
}
double exact_solution400(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary * function_glambda(y,n) * cos(m*PI*x)*exp(-t) + Ty_0 + 2*y*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * function_glambda(y,n) * cos(m*PI*x)*exp(-t) + (heat_conductivity(0,0,0)*Ty_0 + heat_conductivity(0,Ny-1,0)*Ty_1)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0)) + 2*(y-0.5)*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double integral_lambda_x(double x, int m)
{
	if (m == 0)
		return x;
	else
		return (1.0 / (PI * m)) * sin (PI * m * x);
}
double integral_lambda_y(double y, int n)
{
	double nPi = n * PI;
	if (n == 0)
		return y;
	else
		return - (2.0 * y * y * y - 3.0 * y * y) * cos(PI*n*y) / nPi + (6.0 * y * (y - 1)) * sin(PI*n*y) / (nPi * nPi) + (6.0 * (2.0*y - 1)) * cos(PI*n*y) / (nPi * nPi * nPi) - 12.0 * sin(PI*n*y) / (nPi * nPi * nPi * nPi) ;
}
double exact_average_solution_other400(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	integ_x = integral_lambda_x(xpoints[i+1],m) - integral_lambda_x(xpoints[i],m);
	integ_y = integral_lambda_y(ypoints[j+1],n) - integral_lambda_y(ypoints[j],n);

	if (ypoints[j] + 0.5 * hy(j) <= 0.5)
		return nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
}

double function_glambda_gradient(double y, int n)
{
	return 6.0 * y * (y - 1) *sin(PI*n*y) + (2.0 * y * y * y - 3.0 * y * y) * (PI * n) * cos(PI*n*y);
}
double exact_gradient400y(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
	{
		return nonstationary * function_glambda_gradient(y,n) *cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,Ny-1,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
	}
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0)) * nonstationary * function_glambda_gradient(y,n) *cos(m*PI*x)*exp(-t) + 2*heat_conductivity(0,0,0)*(Ty_1 - Ty_0)/(heat_conductivity(0,0,0)+heat_conductivity(0,Ny-1,0));
}
double exact_gradient400x(double x, double y, double z, double t, int m, int n)
{
	if (y<0.5)
		return nonstationary*(-m*PI)*function_glambda(y,n)*sin(m*PI*x)*exp(-t); 
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))*nonstationary*(-m*PI)*function_glambda(y,n)*sin(m*PI*x)*exp(-t); 
}
double function_glambda_2derivative(double y, int n)
{
	return (6.0 * (2.0 * y - 1)) * sin(PI * n * y) + (6.0 * y * (y - 1)) * 2.0 * n * PI * cos(PI * n * y) + (2.0 * y * y * y - 3.0 * y * y) * (- n * PI) * (n * PI) * sin(PI * n * y);
}

double exact_laplace400(double x, double y, double z, double t, int m, int n)
{
	if (y < 0.5)
		return nonstationary * (cos(m*PI*x) * function_glambda_2derivative(y,n) + (-PI*m*PI*m) * cos(m*PI*x) * function_glambda(y,n)) * exp(-t);
	else
		return (heat_conductivity(0,0,0)/heat_conductivity(0,Ny-1,0))* nonstationary * (cos(m*PI*x) * function_glambda_2derivative(y,n) + (-PI*m*PI*m) * cos(m*PI*x) * function_glambda(y,n)) * exp(-t);
}
double righthand400(int m, int n, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
	double temp;
	double x_mid = xpoints[i] + 0.5*hx(i);
	double y_mid = ypoints[j] + 0.5*hy(j);
	double z_mid = zpoints[k] + 0.5*hz(k);
	temp = hx(i) * hy(j) * hz(k) * (-density(i,j,k)*heat_capacity(i,j,k)*exact_solution400_onlytime(x_mid, y_mid, z_mid, t_discr * tau, m, n) - heat_conductivity(i,j,k)*exact_laplace400(x_mid, y_mid, z_mid, t_discr * tau, m, n));
	return temp;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

double arbSmooth_ex1_solution(double x, double y, double z, double t, double mx, double my, double mt)
{
	return nonstationary*sin(mt*PI*t) + sin(mx*PI*x) + sin(my*PI*y);
}
double arbSmooth_ex1_solution_other_boundTy0(int i, int j, double time)
{
	return nonstationary * sin(mt * PI * time) + (-1.0/ (mx * PI)) * (cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i])) / hx(i) ;
}
double arbSmooth_ex1_solution_other_boundTy1(int i, int j, double time)
{
	return nonstationary * sin(mt * PI * time) + (-1.0/ (mx * PI)) * (cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i])) / hx(i) ;
}
double arbSmooth_average_ex1_other_solution(int i, int j, int t_discr, double tau, int mx, int my, int mt)
{
	double integ_x, integ_y, integ_t, volume = hx(i)*hy(j);
	if (mx == 0)
		integ_x = 0.0;
	else
		integ_x = hy(j) * (-1.0/ (mx * PI)) * (cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]));
	if (my == 0)
		integ_y = 0.0;
	else
		integ_y = hx(i) * (-1.0/ (my * PI)) * (cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]));
	integ_t = nonstationary * sin(mt * PI * t_discr * tau) * volume;

	return (1.0/volume) * (integ_x + integ_y + integ_t) ;
}

double arbSmooth_ex1_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt)
{
	return nonstationary*mt*PI*cos(mt*PI*t);
}
double arbSmooth_ex1_gradientX(double x, double y, double z, double t, double mx, double my, double mt)
{
  return (mx*PI)*cos(mx*PI*x); 
}
//double arbSmooth_ex1_gradientX_integ_overX(int i, int j, double t, double mx, double my, double mt)
//{
//	if (mx == 0)
//		return 0.0;
//	else
//		return (1.0 / hx(i)) * ( sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i])); 
//}
double arbSmooth_ex1_gradientY(double x, double y, double z, double t, double mx, double my, double mt)
{
  return (my*PI)*cos(my*PI*y);
}

double arbSmooth_ex1_gradientX_po_t(double x, double y, double z, double t, double mx, double my, double mt)
{
	return 0.0; 
}
double arbSmooth_ex1_gradientY_po_t(double x, double y, double z, double t, double mx, double my, double mt)
{
	return 0.0;
}

double arbSmooth_ex1_righthand(double mx, double my, double mt, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double temp_t, temp_xx, temp_yy, temp, temp_lambdagrad, temp_lambd_x, temp_lambd_y;
  if (external <= 9)
  {
    //temp_t = nonstationary*density(i,j,k)*heat_capacity(i,j,k)*((sin(mt*PI*(t_discr+1)*tau) - sin(mt*PI*t_discr*tau))/tau)*hx(i)*hy(j);
    temp_t = nonstationary*density(i,j,k)*heat_capacity(i,j,k)*(mt*PI)*cos(mt*PI*t_discr*tau)*hx(i)*hy(j);
    //if (mx != 0)
      temp_xx = -heat_conductivity(i,j,k)*(mx*PI)*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]))*hy(j);
    //if (my != 0)
      temp_yy = -heat_conductivity(i,j,k)*hx(i)*(my*PI)*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]));

    temp = temp_t + temp_xx + temp_yy;
  }
  else
  {
    double x_mid = xpoints[i] + 0.5*hx(i);
    double y_mid = ypoints[j] + 0.5*hy(j);
    double z_mid = zpoints[k] + 0.5*hz(k);

    temp_t = nonstationary*density(i,j,k)*heat_capacity(i,j,k)*(mt*PI)*cos(mt*PI*t_discr*tau)*hx(i)*hy(j);
    //if (mx != 0)
      temp_xx = -heat_conductivity(i,j,k)*(mx*PI)*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]))*hy(j);
   // if (my != 0)
      temp_yy = -heat_conductivity(i,j,k)*hx(i)*(my*PI)*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]));

    temp_lambd_x = -hx(i)*hy(j)*a_gradientX(x_mid, y_mid, z_mid)*arbSmooth_ex1_gradientX(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mt);
    temp_lambd_y = -hx(i)*hy(j)*a_gradientY(x_mid, y_mid, z_mid)*arbSmooth_ex1_gradientY(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mt);
    //temp_lambd_x = -a_gradientX(x_mid, y_mid, z_mid)*(sin(mx*PI*xpoints[i+1]) - sin(mx*PI*xpoints[i]))*hy(j);
    //temp_lambd_y = -a_gradientY(x_mid, y_mid, z_mid)*hx(i)*(sin(my*PI*ypoints[j+1]) - sin(my*PI*ypoints[j]));
    temp_lambdagrad = temp_lambd_x + temp_lambd_y;
    temp = temp_t + temp_xx + temp_yy + temp_lambdagrad;

  }

  return hz(k)*temp;
  //	return temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbSmooth_ex2_solution(double x, double y, double z, double t, double mx, double my, double mt)
{
  return (nonstationary*sin(mt*PI*t) + 1) * (sin(mx*PI*x) + 1) * (sin(my*PI*y) + 1) ;//+ x*Tx_1 + (1-x)*Tx_0;
}
double arbSmooth_ex2_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt)
{
	return nonstationary*mt*PI*cos(mt*PI*t) * (sin(mx*PI*x) + 1) * (sin(my*PI*y) + 1) ;//+ x*Tx_1 + (1-x)*Tx_0;
}
double arbSmooth_ex2_gradientX(double x, double y, double z, double t, double mx, double my, double mt)
{
  return (nonstationary*sin(mt*PI*t) + 1)*(mx*PI)*cos(mx*PI*x) * (sin(my*PI*y) + 1) ;//+ Tx_1 - Tx_0; 
}
double arbSmooth_ex2_gradientY(double x, double y, double z, double t, double mx, double my, double mt)
{
  return (nonstationary*sin(mt*PI*t) + 1) * (sin(mx*PI*x) + 1) * (my*PI)*cos(my*PI*y);
}
double arbSmooth_ex2_gradientX_po_t(double x, double y, double z, double t, double mx, double my, double mt)
{
	return nonstationary*mt*cos(mt*PI*t)*(mx*PI)*cos(mx*PI*x) * (sin(my*PI*y) + 1) ;//+ Tx_1 - Tx_0; 
}
double arbSmooth_ex2_gradientY_po_t(double x, double y, double z, double t, double mx, double my, double mt)
{
	return nonstationary*mt*cos(mt*PI*t) * (sin(mx*PI*x) + 1) * (my*PI)*cos(my*PI*y);
}

double arbSmooth_ex2_righthand(double mx, double my, double mt, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double temp_t, temp_xx, temp_yy, temp, temp_lambdagrad, temp_lambd_x, temp_lambd_y, mnogitel_x, mnogitel_y;
  if (external <= 9)
  {
    if (my != 0)
      temp_xx = heat_conductivity(i,j,k)*(nonstationary*sin(mt*PI*t_discr*tau) + 1)*(mx*PI)*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]))*(1.0/(my*PI))*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) - my*PI*hy(j));
    else
      temp_xx = -heat_conductivity(i,j,k)*(nonstationary*sin(mt*PI*t_discr*tau) + 1)*(mx*PI)*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]))*hy(j);

    if (mx != 0)
      temp_yy = heat_conductivity(i,j,k)*(nonstationary*sin(mt*PI*t_discr*tau) + 1)*(1.0/(mx*PI))*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) - mx*PI*hx(i))*(my*PI)*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]));
    else
      temp_yy = -heat_conductivity(i,j,k)*(nonstationary*sin(mt*PI*t_discr*tau) + 1)*hx(i)*(my*PI)*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]));


    if (mx == 0)
      mnogitel_x = hx(i);
    else
      mnogitel_x = -(1.0/(mx*PI))*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) - mx*PI*hx(i));

    if (my == 0)
      mnogitel_y = hy(j);
    else
      mnogitel_y = -(1.0/(my*PI))*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) - my*PI*hy(j));

    //temp_t = nonstationary*density(i,j,k)*heat_capacity(i,j,k)*(mt*PI)*cos(mt*PI*t_discr*tau)*(1.0/(mx*PI))*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) - mx*PI*hx(i))*(1.0/(my*PI))*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) - my*PI*hy(j));
    temp_t = nonstationary*density(i,j,k)*heat_capacity(i,j,k)*(mt*PI)*cos(mt*PI*t_discr*tau)*mnogitel_x*mnogitel_y;


    //printf("temp_xx = %f temp_yy = %f temp_t = %f temp = %f \n",temp_xx, temp_yy, temp_t, temp);
    //_getch();
    temp = temp_t + temp_xx + temp_yy;
  }
  else
  {
    double x_mid = xpoints[i] + 0.5*hx(i);
    double y_mid = ypoints[j] + 0.5*hy(j);
    double z_mid = zpoints[k] + 0.5*hz(k);

    temp_t = nonstationary*density(i,j,k)*heat_capacity(i,j,k)*(mt*PI)*cos(mt*PI*t_discr*tau)*(1.0/(mx*PI))*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) - mx*PI*hx(i))*(1.0/(my*PI))*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) - my*PI*hy(j));
    if (my != 0)
      temp_xx = heat_conductivity(i,j,k)*(nonstationary*sin(mt*PI*t_discr*tau) + 1)*(mx*PI)*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]))*(1.0/(my*PI))*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]) - my*PI*hy(j));
    if (mx != 0)
      temp_yy = heat_conductivity(i,j,k)*(nonstationary*sin(mt*PI*t_discr*tau) + 1)*(1.0/(mx*PI))*(cos(mx*PI*xpoints[i+1]) - cos(mx*PI*xpoints[i]) - mx*PI*hx(i))*(my*PI)*(cos(my*PI*ypoints[j+1]) - cos(my*PI*ypoints[j]));

    temp_lambd_x = -hx(i)*hy(j)*a_gradientX(x_mid, y_mid, z_mid)*arbSmooth_ex2_gradientX(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mt);
    temp_lambd_y = -hx(i)*hy(j)*a_gradientY(x_mid, y_mid, z_mid)*arbSmooth_ex2_gradientY(x_mid, y_mid, z_mid, t_discr*tau, mx, my, mt);
    temp_lambdagrad = temp_lambd_x + temp_lambd_y;
    temp = temp_t + temp_xx + temp_yy + temp_lambdagrad;

  }

  return hz(k)*temp;
  //	return temp;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbNSmooth_ex3_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  //printf("x = %f y = %f t = %f alpha = %f \n",x,y,t,alpha);
  //printf("|x*y*t - 0.25| = %f log |x*y*t - 0.25| = %f \n", fabs(x*y*t - 0.25), log(fabs(x*y*t - 0.25)));
  //printf("arg = %f exp = %f \n",alpha*log(fabs(x*y*t - 0.25)),exp(alpha*log(fabs(x*y*t - 0.25))));
  return exp(alpha*log(fabs(x*y*t - 0.25)));
}
double arbNSmooth_ex3_gradientX(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if (x*y*t - 0.25 >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(y*t);
  else
    return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(-y*t);
}
double arbNSmooth_ex3_gradientY(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if (x*y*t - 0.25 >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(x*t);
  else
    return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(-x*t);
}
double arbNSmooth_ex3_gradientX_po_t(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
	if (x*y*t - 0.25 >= 0)
		return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(y) + alpha*(alpha - 1.0)*exp((alpha - 2.0)*log(fabs(x*y*t - 0.25)))*(y*t)*(x*y);
	else
		return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(-y) + alpha*(alpha - 1.0)*exp((alpha - 2.0)*log(fabs(x*y*t - 0.25)))*(y*t)*(x*y);
}
double arbNSmooth_ex3_gradientY_po_t(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
	if (x*y*t - 0.25 >= 0)
		return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(x) + alpha*(alpha - 1.0)*exp((alpha - 2.0)*log(fabs(x*y*t - 0.25)))*(x*t)*(x*y);
	else
		return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(-x) + alpha*(alpha - 1.0)*exp((alpha - 2.0)*log(fabs(x*y*t - 0.25)))*(x*t)*(x*y);
}

double arbNSmooth_ex3_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double volume = hx(i)*hy(j)*hz(k);
  double x_mid = xpoints[i] + 0.5*hx(i);
  double y_mid = ypoints[j] + 0.5*hy(j);
  double z_mid = zpoints[k] + 0.5*hz(k);
  double temp = volume*(density(i,j,k)*heat_capacity(i,j,k)*arbNSmooth_ex3_solution_po_t(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha) - heat_conductivity(i,j,k)*arbNSmooth_ex3_solution_laplas(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  //temp += - volume*a_gradientX(x_mid, y_mid, z_mid)*arbNSmooth_ex3_gradientX(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha);
  //temp += - volume*a_gradientY(x_mid, y_mid, z_mid)*arbNSmooth_ex3_gradientY(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha);
  return temp;
}
double arbNSmooth_ex3_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if (x*y*t - 0.25 >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(x*y);
  else
    return alpha*exp((alpha - 1.0)*log(fabs(x*y*t - 0.25)))*(-x*y);

}
double arbNSmooth_ex3_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  double t_xx, t_yy;
  t_xx =  alpha*(alpha - 1.0)*(y*t)*(y*t)*exp((alpha - 2.0)*log(fabs(x*y*t - 0.25)));
  t_yy =  alpha*(alpha - 1.0)*(x*t)*(x*t)*exp((alpha - 2.0)*log(fabs(x*y*t - 0.25)));
  return t_xx + t_yy;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbNSmooth_ex4_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  //printf("x = %f y = %f t = %f alpha = %f \n",x,y,t,alpha);
  //printf("|x*y*t - 0.25| = %f log |x*y*t - 0.25| = %f \n", fabs(x*y*t - 0.25), log(fabs(x*y*t - 0.25)));
  //printf("arg = %f exp = %f \n",alpha*log(fabs(x*y*t - 0.25)),exp(alpha*log(fabs(x*y*t - 0.25))));
  return exp(alpha*log(fabs(x*y - 0.25)));
}
double arbNSmooth_ex4_gradientX(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if (x*y - 0.25 >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs(x*y - 0.25)))*(y);
  else
    return alpha*exp((alpha - 1.0)*log(fabs(x*y - 0.25)))*(-y);
}
double arbNSmooth_ex4_gradientY(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if (x*y - 0.25 >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs(x*y - 0.25)))*(x);
  else
    return alpha*exp((alpha - 1.0)*log(fabs(x*y - 0.25)))*(-x);
}

double arbNSmooth_ex4_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double volume = hx(i)*hy(j)*hz(k);
  double x_mid = xpoints[i] + 0.5*hx(i);
  double y_mid = ypoints[j] + 0.5*hy(j);
  double z_mid = zpoints[k] + 0.5*hz(k);
  double temp = volume*(- heat_conductivity(i,j,k)*arbNSmooth_ex4_solution_laplas(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  return temp;
}
double arbNSmooth_ex4_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  double t_xx, t_yy;
  t_xx =  alpha*(alpha - 1.0)*(y)*(y)*exp((alpha - 2.0)*log(fabs(x*y - 0.25)));
  t_yy =  alpha*(alpha - 1.0)*(x)*(x)*exp((alpha - 2.0)*log(fabs(x*y - 0.25)));
  return t_xx + t_yy;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbNSmooth_ex5_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  //printf("x = %f y = %f t = %f alpha = %f \n",x,y,t,alpha);
  //printf("|x*y*t - 0.25| = %f log |x*y*t - 0.25| = %f \n", fabs(x*y*t - 0.25), log(fabs(x*y*t - 0.25)));
  //printf("arg = %f exp = %f \n",alpha*log(fabs(x*y*t - 0.25)),exp(alpha*log(fabs(x*y*t - 0.25))));
  return exp(alpha*log(fabs((x - 0.5)*(y - 0.5))));
}
double arbNSmooth_ex5_gradientX(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if ((x-0.5)*(y-0.5) >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs((x-0.5)*(y-0.5))))*(y-0.5);
  else
    return alpha*exp((alpha - 1.0)*log(fabs((x-0.5)*(y-0.5))))*(-1)*(y-0.5);
}
double arbNSmooth_ex5_gradientY(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if ((x-0.5)*(y-0.5) >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs((x-0.5)*(y-0.5))))*(x-0.5);
  else
    return alpha*exp((alpha - 1.0)*log(fabs((x-0.5)*(y-0.5))))*(-1)*(x-0.5);
}
double arbNSmooth_ex5_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double volume = hx(i)*hy(j)*hz(k);
  double x_mid = xpoints[i] + 0.5*hx(i);
  double y_mid = ypoints[j] + 0.5*hy(j);
  double z_mid = zpoints[k] + 0.5*hz(k);
  double temp = volume*(- heat_conductivity(i,j,k)*arbNSmooth_ex5_solution_laplas(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  return temp;
}
double arbNSmooth_ex5_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  double t_xx, t_yy;
  t_xx =  alpha*(alpha - 1.0)*(y-0.5)*(y-0.5)*exp((alpha - 2.0)*log(fabs((x-0.5)*(y-0.5))));
  t_yy =  alpha*(alpha - 1.0)*(x-0.5)*(x-0.5)*exp((alpha - 2.0)*log(fabs((x-0.5)*(y-0.5))));
  return t_xx + t_yy;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbNSmooth_ex6_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  //printf("x = %f y = %f t = %f alpha = %f \n",x,y,t,alpha);
  //printf("|x*y*t - 0.25| = %f log |x*y*t - 0.25| = %f \n", fabs(x*y*t - 0.25), log(fabs(x*y*t - 0.25)));
  //printf("arg = %f exp = %f \n",alpha*log(fabs(x*y*t - 0.25)),exp(alpha*log(fabs(x*y*t - 0.25))));
  return exp(alpha*log(fabs((x - 0.5))));
}
double arbNSmooth_ex6_gradientX(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if ((x-0.5) >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs((x-0.5))))*1;
  else
    return alpha*exp((alpha - 1.0)*log(fabs((x-0.5))))*(-1);
}
double arbNSmooth_ex6_gradientY(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  return 0.0;
}
double arbNSmooth_ex6_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double volume = hx(i)*hy(j)*hz(k);
  double x_mid = xpoints[i] + 0.5*hx(i);
  double y_mid = ypoints[j] + 0.5*hy(j);
  double z_mid = zpoints[k] + 0.5*hz(k);
  double temp = volume*(- heat_conductivity(i,j,k)*arbNSmooth_ex6_solution_laplas(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  return temp;
}
double arbNSmooth_ex6_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  double t_xx, t_yy;
  t_xx =  alpha*(alpha - 1.0)*exp((alpha - 2.0)*log(fabs((x-0.5))));
  t_yy = 0.0;
  return t_xx + t_yy;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbNSmooth_ex7_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  //printf("x = %f y = %f t = %f alpha = %f \n",x,y,t,alpha);
  //printf("|x*y*t - 0.25| = %f log |x*y*t - 0.25| = %f \n", fabs(x*y*t - 0.25), log(fabs(x*y*t - 0.25)));
  //printf("arg = %f exp = %f \n",alpha*log(fabs(x*y*t - 0.25)),exp(alpha*log(fabs(x*y*t - 0.25))));
  return exp(alpha*log(fabs((y - 0.5))));
}
double arbNSmooth_ex7_gradientY(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if ((y-0.5) >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs((y-0.5))))*1;
  else
    return alpha*exp((alpha - 1.0)*log(fabs((y-0.5))))*(-1);
}
double arbNSmooth_ex7_gradientX(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  return 0.0;
}
double arbNSmooth_ex7_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double volume = hx(i)*hy(j)*hz(k);
  double x_mid = xpoints[i] + 0.5*hx(i);
  double y_mid = ypoints[j] + 0.5*hy(j);
  double z_mid = zpoints[k] + 0.5*hz(k);
  double temp = volume*(- heat_conductivity(i,j,k)*arbNSmooth_ex7_solution_laplas(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  return temp;
}
double arbNSmooth_ex7_solution_laplas(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  double t_xx, t_yy;
  t_yy =  alpha*(alpha - 1.0)*exp((alpha - 2.0)*log(fabs((y-0.5))));
  t_xx = 0.0;
  return t_xx + t_yy;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
double arbNSmooth_ex8_solution(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  //printf("x = %f y = %f t = %f alpha = %f \n",x,y,t,alpha);
  //printf("|x*y*t - 0.25| = %f log |x*y*t - 0.25| = %f \n", fabs(x*y*t - 0.25), log(fabs(x*y*t - 0.25)));
  //printf("arg = %f exp = %f \n",alpha*log(fabs(x*y*t - 0.25)),exp(alpha*log(fabs(x*y*t - 0.25))));
  return exp(alpha*log(fabs((t - 0.25))));
}
double arbNSmooth_ex8_gradientY(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  return 0.0;
}
double arbNSmooth_ex8_gradientX(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  return 0.0;
}

double arbNSmooth_ex8_righthand(double mx, double my, double mt, double alpha, int t_discr, int i, int j, int k, double* xpoints, double* ypoints, double* zpoints)
{
  double volume = hx(i)*hy(j)*hz(k);
  double x_mid = xpoints[i] + 0.5*hx(i);
  double y_mid = ypoints[j] + 0.5*hy(j);
  double z_mid = zpoints[k] + 0.5*hz(k);
  double temp = volume*(density(i,j,k)*heat_capacity(i,j,k)*arbNSmooth_ex8_solution_po_t(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  //printf("i = %d j = %d k = %d \n",i,j,k);
  //printf("x_mid = %f, y_mid = %f, z_mid = %f, mx = %f, my = %f, mz = %f, alpha = %alpha\n",x_mid, y_mid, z_mid, mx, my, mz, alpha);
  //printf("temp = %15.15f volume = %15.15f solution_po_t = %15.15f \n",temp, volume, arbNSmooth_ex8_solution_pot(x_mid, y_mid, z_mid, t_discr*tau, mx,my,mt,alpha));
  //_getch();
  return temp;
}
double arbNSmooth_ex8_solution_po_t(double x, double y, double z, double t, double mx, double my, double mt, double alpha)
{
  if (t - 0.25 >= 0)
    return alpha*exp((alpha - 1.0)*log(fabs(t - 0.25)))*(1);
  else
    return alpha*exp((alpha - 1.0)*log(fabs(t - 0.25)))*(-1);

}

///////////////////////////////////////////////////////
double Tx_0_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
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
    break;
  default:
    return exact_solution(numSolution,0,ypoints[j] + 0.5*hy(j),0,time);
    break;
  }
}

double Tx_1_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
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
    break;
  default:
    return exact_solution(numSolution,1,ypoints[j] + 0.5*hy(j),0,time);
    break;
  }
}


double Ty_0_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
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
    break;
  default:
    return exact_solution(numSolution,xpoints[i] + 0.5*hx(i),0,0,time);
    break;
  }
}

double Ty_1_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
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
    break;
  default:
    return exact_solution(numSolution,xpoints[i] + 0.5*hx(i),1,0,time);
    break;
  }
}

double wx_0_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
	//return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX_integ(i, j, time,mx,my,mt);
    break;
  case 109:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex2_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 110:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex3_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 111:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex4_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 112:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex5_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 113:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex6_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 114:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex7_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 115:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex8_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  default:
    return wx_0;
    break;
  }
}

double wx_0_bound_noglobal(int numSolution, int i, int j, double time, double * xpoints_loc, double * ypoints_loc)
{
  switch (numSolution)
  {
  case 108:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
	//return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX_integ(i, j, time,mx,my,mt);
    break;
  case 109:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex2_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 110:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex3_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 111:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex4_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 112:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex5_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 113:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex6_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 114:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex7_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  case 115:
    return -heat_conductivity_func(0,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex8_gradientX(0,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt, alpha);
    break;
  default:
    return wx_0;
    break;
  }
}

double wx_1_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
	//return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex1_gradientX_integ(i, j, time,mx,my,mt);
    break;
  case 109:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbSmooth_ex2_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt);
    break;
  case 110:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex3_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
    break;
  case 111:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex4_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
    break;
  case 112:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex5_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
    break;
  case 113:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex6_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
    break;
  case 114:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex7_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
    break;
  case 115:
    return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*arbNSmooth_ex8_gradientX(1,ypoints[j] + 0.5*hy(j),0,time,mx,my,mt,alpha);
    break;
  default:
    //return -heat_conductivity_func(1,ypoints[j] + 0.5*hy(j),0)*exact_gradientX(numSolution,1,ypoints[j] + 0.5*hy(j),0,time);
    return wx_1;
    break;
  }
}

double wy_0_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbSmooth_ex1_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt);
    break;
  case 109:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbSmooth_ex2_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt);
    break;
  case 110:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex3_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
    break;
  case 111:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex4_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
    break;
  case 112:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex5_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
    break;
  case 113:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex6_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
    break;
  case 114:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex7_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
    break;
  case 115:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*arbNSmooth_ex8_gradientY(xpoints[i] + 0.5*hx(i),0,0,time,mx,my,mt,alpha);
    break;
  case 300:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*exact_gradient300y(xpoints[i] + 0.5*hx(i),0,0,time,mx,my);
    break;
  default:
    //return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),0,0)*exact_gradientX(numSolution,xpoints[i] + 0.5*hx(i),0,0,time);
    return wy_0;
    break;
  }
}

double wy_1_bound(int numSolution, int i, int j, double time)
{
  switch (numSolution)
  {
  case 108:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbSmooth_ex1_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt);
    break;
  case 109:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbSmooth_ex2_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt);
    break;
  case 110:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex3_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
    break;
  case 111:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex4_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
    break;
  case 112:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex5_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
    break;
  case 113:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex6_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
    break;
  case 114:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex7_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
    break;
  case 115:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*arbNSmooth_ex8_gradientY(xpoints[i] + 0.5*hx(i),1,0,time,mx,my,mt,alpha);
    break;
  case 300:
    return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*exact_gradient300y(xpoints[i] + 0.5*hx(i),1,0,time,mx,my);
    break;
  default:
    //return -heat_conductivity_func(xpoints[i] + 0.5*hx(i),1,0)*exact_gradientX(numSolution,xpoints[i] + 0.5*hx(i),1,0,time);
    return wy_1;
    break;
  }
}

double a_gradientX(double x, double y, double z)
{
  switch(external)
  {
  case 10:
    return arb_a2_coeff_gradX(x,y,z);
    break;
  case 11:
    return arb_a3_coeff_gradX(x,y,z);
    break;
  case 12:
    return arb_a4_coeff_gradX(x,y,z);
    break;
  default:
    return 0.0;
    break;
  }
}

double a_gradientY(double x, double y, double z)
{
  switch(external)
  {
  case 10:
    return arb_a2_coeff_gradY(x,y,z);
    break;
  case 11:
    return arb_a3_coeff_gradY(x,y,z);
    break;
  case 12:
    return arb_a4_coeff_gradY(x,y,z);
    break;
  default:
    return 0.0;
    break;
  }
}

double exact_average_solution_other105(int i, int j, int t_discr, double tau, int m, int n)
{
	double integ_x, integ_y, volume = hx(i)*hy(j);
	if (m == 0)
		integ_x = hx(i);
	else
		integ_x =  (1.0 / (PI * m)) * (sin (PI * m * xpoints[i+1]) - sin (PI * m * xpoints[i]));

	if (n == 0)
		integ_y = hy(j);
	else
		integ_y =  - (1.0 / (PI * n)) * (cos (PI * n * ypoints[j+1]) - cos (PI * n * ypoints[j]));

	return nonstationary * exp (-t_discr * tau) * integ_y * integ_x / volume;
}
