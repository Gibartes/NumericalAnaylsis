#pragma once

#ifndef __MATRIX_H__
#define __MATRIX_H__
#define QUADRIC 3
#define S_CUBIC 4
#define _MT_SIN 	((unsigned int)0)
#define _MT_COS 	((unsigned int)1)
#define _MT_TAN 	((unsigned int)2)
#define _MT_ASIN 	((unsigned int)3)
#define _MT_ACOS 	((unsigned int)4)
#define _MT_ATAN 	((unsigned int)5)
#define _MT_EXP 	((unsigned int)6)
#define _MT_LN 		((unsigned int)7)
#define _MT_LOG10	((unsigned int)8)
#define _MT_INV		((unsigned int)9)
#define _MT_SQRT	((unsigned int)10)
#define _MT_RADIAN	((unsigned int)11)
#define _MT_DEGREE	((unsigned int)12)
#ifdef _MSC_VER
#define dllpublic extern __declspec(dllimport)
#else
#define dllpublic extern
#endif
/******************************************************************************* 
* This file is a matrix library / header for assginment                        *
*                                                                              *
* Made by Gibartes                                    *
*                                                                              *
*******************************************************************************/

/* -------------------------------------------------------------------------- */

typedef struct __vector{
	unsigned int dim;		/* the number of valid ingredients										*/
	unsigned int max;		/* reseved. max capacity of vector  (dim <= max. default dim=max) 	*/
	double *element;		/* element array in vector											*/
	int flag;				/* reseved. (do not change this)									*/
}Vector, *pVector;

typedef struct __result{
	unsigned int flg;		/* error flag or index flag */
	double value;			/* return value */
	double error;			/* error value */
	void *ptr;				/* return pointer */
}result;

typedef struct __solution{
	unsigned int dim;		/* dimension of solution */
	Vector *index;			/* ingredient initial position of solution */
	Vector *solution;		/* solution vector  */
}Solution;

typedef struct _Poly {
	double  c;
	int 	e;
	char	symbol;
	int 	overwrite;
	struct _Poly *link;
}Poly;
/* -------------------------------------------------------------------------- */

/* construct result object */

dllpublic result construct(int flg, double value, double error, void *ptr);

/* -------------------------------------------------------------------------- */
dllpublic result get_module_version(void);
/* create vector with double array */

/* -------------------------------------------------------------------------- */
/* Vector ?????? ?????? */

/* create vector from static allocated double array */
dllpublic Vector *vector_create(unsigned int dim, double *elem);
/* create vector from dynamic allocated double array */
dllpublic Vector *vector_array(unsigned int dim, double *elem);
/* create vector from string expression */
dllpublic Vector *vector_string(unsigned int dim, char *expr, size_t size);
/* create O-vector */
dllpublic Vector *zero_vector(unsigned int dim);
/* create vector which has all ingredients value is 1 */
dllpublic Vector *one_vector(unsigned int dim);
/* create index vector */
dllpublic Vector *index_vector(unsigned int dim);
/* release vector */
dllpublic void erase_vector(Vector *src);
/* release matrix */
dllpublic void erase_matrix(Vector **M, const unsigned int row);

/* -------------------------------------------------------------------------- */
/* Vector & Solution ?????? ?????? */ 

/* print vector by rows */
dllpublic void print_vector(Vector *src);
/* print vector by columns */
dllpublic void print_col(Vector *src);
/* print solution if mode == 0 then, print row vector type */
dllpublic void print_solution(Solution *x, int mode);
/* print matrix solution system */
dllpublic void print_sys(Vector *src, Vector *base, unsigned int line_no);
/* print matrix */
dllpublic void print_matrix(Vector **matrix, const unsigned int row);

/* -------------------------------------------------------------------------- */
/* vector base arithmetic (????????????) */

/* set vector dimension : dim <= max */
dllpublic result set_dim(Vector *x, unsigned int dim);
/* restore vector dimension */
dllpublic result restore_max(Vector *x);
/* get current ingredient of vector[index] */
dllpublic result get_value(Vector *vec,const unsigned int index);
/* change one ingredient in vector */
dllpublic void change(Vector *src, const unsigned int index, const double value);
/* swap two ingreients in one vector */
dllpublic result swap(Vector *src, const unsigned int index_o, const unsigned int index_t);
/* interchange position between two vectors */
dllpublic result interchange(Vector *src, Vector *dst);
/* swap vector and original vector is freed */
dllpublic result vector_move(Vector *src, Vector *dst);
/* copy vector */
dllpublic result vector_copy(Vector *src, Vector *dst);
/* summmation of ingredient of vector */
dllpublic double ingredient_sum(Vector *src);
/* vector size */
dllpublic double vector_length(Vector *src);
/* get vector dimension */
dllpublic  unsigned int dims(Vector *vec);
/* vector_head_normalize is equivalent with elimination + mult */
/* ????????? 0??? ?????? head ????????? ????????? ????????? vector head ??? 1??? ?????? */
dllpublic void vector_head_normalize(Vector *ori);
/* ?????? ?????? ?????????(?????????) */
dllpublic result find_max(Vector *src, unsigned int start);
/* ?????? ?????? ????????? ?????? (bound ??????) */
dllpublic result vector_cmp(Vector *src, Vector *dst, const double bound, unsigned int start, unsigned int end);
/* vector addition and subtraction */
dllpublic result add(Vector *src, Vector *dst, Vector *new, const char op);
/* inner product between two vectors */
dllpublic result dot(Vector *src, Vector *dst);
/* ?????? ?????? ????????? */
dllpublic result mult(Vector*src, Vector *new, const double c, unsigned int start);
/* ????????? ??? ?????? ?????? ?????? */
dllpublic result point_add(Vector *src, unsigned int inx, double pt);
/* ?????? ????????? ?????? ?????? ????????? ?????? ?????? op<0?????? ?????? */
dllpublic result vector_add(Vector *src, Vector *dst, int op, unsigned int start, unsigned int end);
/* ?????? ????????? ?????? ?????? ????????? ?????? ????????? */
dllpublic result vector_mult(Vector *src, Vector *dst, unsigned int start, unsigned int end);
/* ?????? ????????? ?????? ?????? ????????? ?????? ?????????  ????????? */
dllpublic result vector_divide(Vector *src, Vector *dst, unsigned int start, unsigned int end,unsigned int mode);
/* ????????? ?????? ????????? ?????????????????? ????????? ??????????????? ????????????. */
dllpublic result vector_conversion(Vector *src, unsigned int func);


/* -------------------------------------------------------------------------- */

/* ?????? ?????? ?????? */

/* make column vector from row vector organzied matrix n x n */
dllpublic Vector *column_vector(Vector **M, unsigned int col);
/* make column vector from row vector organzied matrix m x n */
dllpublic Vector *column_vector_mn(Vector **M, unsigned int row, unsigned int col);
/* make column vector from row vector organzied matrix m x n */
dllpublic result fast_column_vector(Vector **M, Vector *v, unsigned int row, unsigned int col);
/* excavate diagonal ingredients in matrix */
dllpublic result trace(Vector **Matrix, Vector *diagonal);
/* transpose matrix */
dllpublic result transpose(Vector **original, Vector **tmatrix, unsigned int row);
/* ?????? ????????? ?????????(?????????) */
dllpublic result find_max_col(Vector **M, unsigned int row, unsigned int start,unsigned int col);
/* ?????? ????????? 0??? ?????? */
dllpublic void matrix_zeros(Vector **M, const unsigned int row, const unsigned int col);
/* ?????? ?????? */
dllpublic void matrix_eigen(Vector **M, const unsigned int row);
/* Ax = b ?????? b?????? ?????? */
dllpublic result matrix_single_mult(Vector **M, Vector *x, Vector *b, unsigned int row);
/* AB = M ?????? ?????? */
dllpublic result matrix_mult(Vector **A, Vector **B, Vector **M, unsigned int rowA, unsigned int rowB);
/* A[j][i] = A[i][j] */ 
dllpublic result symmestry(Vector **A, unsigned int row);
/* ?????? ?????? ?????? |aii| > |sum(aij)-aii| 
*  flag = 1 : diagonally dominant 		*/
dllpublic int Diagonally_Dominant(Vector **M,unsigned int row);

/*-------------------------------------------------------------------------- */

/* Solution ?????? ????????? ?????? */ 

/* consider "x" vector in Ax = b*/
/* Solution ?????? ?????? */
dllpublic Solution *solution_set(unsigned int dim);
/* x?????? ????????? */
dllpublic void solution_reset(Solution *x, const unsigned int dim);
/* ???????????? ????????? x?????? data ????????? ?????? */
dllpublic void solution_value_set(Solution *x, const unsigned int index, const double data);
/* ?????? src ?????? ????????? ????????? ??? ?????? ?????? ????????? ?????? */
dllpublic void solution_vector_set(Vector *src, Solution *x);
/* ???????????? ?????? ?????? ???????????? */
dllpublic void swap_solution_all(Solution *x,const unsigned int index_o, const unsigned int index_t);
/* ???????????? ???????????? */
dllpublic void swap_solution(Solution *x,const unsigned int index_o, const unsigned int index_t);
/* ????????? ?????? */
dllpublic void erase_solution(Solution *x);
/* ????????? ????????? ????????????. */
dllpublic void solution_vector(Solution *x, Vector *target);
/* ????????? ????????? ???????????? ????????? ?????????. */
dllpublic void solution_handle(Solution *x, Vector *target);
/* copy solution vector to target vector */
dllpublic void get_solution_vector(Solution *x, Vector *target);

/*-------------------------------------------------------------------------- */
/* ?????? ?????? */

/* ????????? ????????? */
/* ????????? ?????? */
dllpublic result elimination(Vector *ori, Vector *tgt, unsigned int col);
/* b ?????? ?????? : row0 - ?????????, row1 - ????????? */
dllpublic result elimination_base(Vector *base, double m, unsigned int row0, unsigned int row1);
/* Pivot : elimination ready process */
/* ??? ?????? PIVOTING */
dllpublic result row_pivoting(Vector **M, Solution *x, Vector *base, unsigned int row);
/* ??? ?????? PIVOTING */
dllpublic result col_pivoting(Vector **M, Solution *x, Vector *base, unsigned int row);
/* ?????? PIVOTING */
dllpublic result pivoting(Vector **M, Solution *x, Vector *base, unsigned int row);
/* ???????????? */
dllpublic result calculate_one_line(Vector *row, Vector *base, Solution *x, unsigned int row_no);
/* ???????????? */
dllpublic result pre_calculate_one_line(Vector *row, Vector *base, Solution *x, unsigned int row_no);
/* ????????? ????????? ????????? ????????? ???????????? */
dllpublic result arrange_matrix_by_index(Vector **M, Vector *base, Solution *x);
/* ?????? ????????? ????????? */
dllpublic result gauss_elimination(Vector **M, Vector *base, Solution *x);


/* LU ?????? */
/* U ?????? ????????? - 1??? ?????? ??? ?????? ??? ??? 0?????? ????????? */
dllpublic result U_init(Vector **ori, Vector **L, Vector **A);
/* L ?????? ????????? - 1??? ?????? ??? ???????????? ?????? ??? 0?????? ????????? */
dllpublic result L_init(Vector **ori, Vector **A);
/* LU ?????? ?????? */
dllpublic result LU_decompose(Vector **L, Vector **U, Vector **A);



/* ????????? ????????? */
dllpublic result Jacobi_single_step(Vector **M, Vector *base, Solution *x);
/* ????????? ????????? ????????? */
dllpublic result Gauss_Seidel_single_step(Vector **M, Vector *base, Solution *x);



/* ????????? ????????? ????????? (?????????) */
dllpublic result newton_interpolation(Vector **Matrix, Vector *data_x, Vector *data_y);
/* ????????? ?????? ????????? (????????????)*/
dllpublic result post_newton_interpolation(Vector *Newtons, Vector *x, Vector *y);
/* ????????? ?????? */
dllpublic result newton_interpolation_calculate(Vector *data_x, Vector *result, double value);



/* 2??? ???????????? ?????? ????????? ?????? ?????? */
dllpublic result quadric_spline(Vector **q, Vector *x, Vector *y);
/* 2??? ???????????? ?????? ???????????? ??? ?????? ?????? ?????? */
dllpublic result quadric_spot(Vector **q, Vector *x, double value);
/* 2??? ???????????? ?????? ???????????? ??????????????? ????????? time??????????????? ??????????????? ?????? mode 1 : reset */
dllpublic result quadric_graph(Vector **q, Vector *x, Vector *time, Vector *output,int mode);




dllpublic result cubic_spot(Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Vector *s, double value, int mode);
#ifdef _MSC_VER
/* declare cubic matrix yourself -> ROW&COL must be dp-2 */
dllpublic result cubic_spline(Vector **cubic, Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Solution *s, unsigned int dp, int mode);
#else
/* automatic declare & free cubic matrix */
dllpublic result cubic_spline(Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Solution *s, unsigned int dp, int mode);
#endif
/* You have to declare more infomation than basical function or pre-build system */
dllpublic result cublic_spline_pro(Vector **cubic, Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Solution *s, unsigned int dp, int mode);




/* ?????? ?????? ?????? ?????? */
dllpublic result least_square_matrix(Vector **A, Vector *x, Vector *y, Solution *base, unsigned int e);
/* ?????? ?????? ?????? ?????? ?????? */
#ifdef _MSC_VER
dllpublic result polyfit(Vector **A, Vector *x, Vector *y, Solution *base, unsigned int e);
#else
dllpublic result polyfit(Vector *x, Vector *y, Solution *base, unsigned int e);
#endif
/* ?????? ?????? ?????? ?????? ??? */
dllpublic result polyval(Vector *coeff, double x);
/* RSQUARE ??? */
dllpublic result rsquare(Vector *coeff, Vector *x, Vector *y);
/* SSE(ERROR ???) */
dllpublic result least_error(Vector *coeff, Vector *x, Vector *y);




/*-------------------------------------------------------------------------- */
/* ????????? 																	 */

dllpublic Poly *poly_set(char symbol);								// create polynomial object at memory 
dllpublic Poly *poly_concat(Poly *head1, Poly *head2);				// concat' between two polynomial objects 
dllpublic void poly_insert(Poly **head, Poly *p, int e, double c);	// insert poly's
dllpublic void poly_pop(Poly **head, Poly *p, Poly *removed);		// pop term of poly'
dllpublic void poly_clean(Poly *head);								// release polynomial object
dllpublic void poly_dump(Poly *head);								// reset polynomial object ( start with 0 )
dllpublic void poly_show(Poly *head);								// print out current poly'
dllpublic void poly_differentiate(Poly *dst, Poly *ori);			// defferentiate poly' 
dllpublic void poly_integrate(Poly *dst, Poly *ori, double C);		// integrate poly' ( C : integral constant )
dllpublic void poly_expansion(Poly *p, unsigned int n);				// build finite expansion with poly' to n
dllpublic void set_overwrite(Poly *head);							// set overwrite flag (add, integrate, differentiate)
dllpublic void poly_const_mult(Poly *derv, Poly *ori, double k);	// multiply constant
dllpublic result poly_add(Poly *new, Poly *head1, Poly *head2, int op);	// add(+)/sub(-) between two poly's
dllpublic result poly_calculate(Poly *head,double v);					// get a value of poly at V
dllpublic result poly_swap(Poly *head, int e, int new_e, double c,int op);// change ingredient of poly'
dllpublic result poly_extract(Poly *head, Vector *contain[2]);		// extract from to Vector[2] (exp, coeff)
dllpublic result poly_from_vector(Poly *head, Vector *contain[2]);	// construct poly' from Vector[2] (exp, coeff)

/*-------------------------------------------------------------------------- */


/* ?????? ????????? 																 */

/* 1??? ????????? For static coefficient & exponent */
dllpublic result diff_static(double (*func)(double),double x, double h);
/* 1???????????? For dynamic coefficient & exponent */
dllpublic result diff_dynamic(double (*func)(Vector *, Vector *, double, int),Vector *c, Vector *e, int i, double x, double h);
/* 1???????????? For numerical data points ( except for first  point & last point ) */
dllpublic result diff_numeric(Vector *x, Vector *y, unsigned int index);
/* Newton optimization  */
dllpublic result newton_minimalization(double (*func)(double),double ini, double fin, double precision);
/* simple trapezoidal method */
dllpublic result trapezoidal(Vector *y, double h);
/* simple simpson method */
dllpublic result simpson_quad(Vector *y, double h);
/* generate Romberg's quadnature matrix */
dllpublic void romberg_matrix(Vector **Romberg, double (*func)(double), unsigned int cols, double info[3]);
/* get value of Romberg's quadnature from Romberg Matrix */
dllpublic result romberg_integrate(Vector **Romberg, unsigned int row, unsigned int col);
/* get value of Romberg's quadnature */
dllpublic result integrate(double (*func)(double), unsigned int iter, double a, double b);
/* solve 1-st diff equation with euler method */
dllpublic result euler_method(double (*func)(double,double), double bot, double top, double init,unsigned int term);
/* solve 1-st system of differential equations using Runge-Kutta method order 4 */
dllpublic result RK4(double (*func)(double,double), double bot, double top, double init, unsigned int term);
/* solve 2-nd system of differential equation using Runge-Kutta method order 4 */	
dllpublic result RK4_2(double (*func)(double,double,double),double info[4], unsigned int N);
/* solve 2-nd system of differential equation using Runge-Kutta method order 4 */
dllpublic result RK4_22(double (*f)(double,double,double), double (*g)(double,double,double),double info[4], unsigned int term);
/* solve 1-st system of differential equations using Adamas Bashforth method order 4 */
dllpublic result AB4(double (*func)(double,double), Vector *ytable, double bot, double top, double init, unsigned int term);
/* solve 1-st system of differential equations using Adamas Moulton method order 4 */
dllpublic result AM4(double (*func)(double,double), Vector *ytable, double bot, double top, double init, unsigned int term);
/* 2-order finite difference with matrix, Solve y" + p(x)y' + q(x)y = r(x) */
/* info : [start,end,start_value,end_value] return : solution vector (not include two initial values) */
#ifdef _MSC_VER
dllpublic result fdf2_M(Vector **M, Solution *x,double (*p)(double),double (*q)(double),double (*r)(double),double info[4], unsigned int term);
#else
dllpublic result fdf2_M(Solution *x,double (*p)(double),double (*q)(double),double (*r)(double),double info[4], unsigned int term);
#endif
#endif
