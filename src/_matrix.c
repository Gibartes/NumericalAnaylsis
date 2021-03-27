/******************************************************************************* 
* This file is a matrix library / header for assginment                        *
*                                                                              *
* Made by Gibartes                                                             *
*                                                                              *
*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define etha 	1e-15
#ifndef _MATRIX_VERSION
#define _MATRIX_VERSION "1.1.0"
#endif
#ifndef __MATRIX_H__
#define QUADRIC 3
#define S_CUBIC 4
#define private static inline
#ifdef _MSC_VER
#define public __declspec(dllexport)
#else
#define public
#endif
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
#define LCHILD(x) 2 * (x) + 1
#define RCHILD(x) 2 * (x) + 2
#define PARENT(x) ((x) - 1) / 2
#endif

/* -------------------------------------------------------------------------- */

private double _vinverse(double x){
	if(x!=0){return 1/x;}
	else{return strtod("Inf", NULL);}}
private double _rad(double x){return (x)*180.0*M_PI;}
private double _deg(double x){return (x)*180.0/M_PI;}
private double _rksys2(double t, double x, double y){return x;}
	
public typedef struct __vector{
	unsigned int dim;		/* current the number of ingredients 					*/
	unsigned int max;		/* max capacity of vector (dim <= max. default dim=max) */
	double *element;		/* element array in vector   							*/
	int flag;				/* reseved */
}Vector, *pVector;

public typedef struct __result{
	unsigned int flg;		/* flg flag */
	double value;			/* return value */
	double error;
	void *ptr;				/* return pointer */
}result;

public typedef struct __solution{
	unsigned int dim;		/* dimension of solution */
	Vector *index;			/* ingredient initial position of solution */
	Vector *solution;		/* solution vector  */
}Solution;

public typedef struct _Poly{
	double  c;
	int 	e;
	char	symbol;
	int 	overwrite;
	struct _Poly *link;
}Poly;

private struct mftables{
	unsigned int index;
	double (*eval)(double x);
}mFT[13]={
	{0, sin},
	{1, cos},
	{2, tan},
	{3, asin},
	{4, acos},
	{5, atan},
	{6, exp},
	{7, log},
	{8, log10},
	{9, _vinverse},
	{10, sqrt},
	{11,_rad},
	{12,_deg},
};

private struct mftables *mftfind(unsigned int index){
	if(index<13){return mFT+index;}
	else{return mFT;}
}

private struct mftables mfts;

/* -------------------------------------------------------------------------- */

/* construct result object */

public result construct(int flg, double value, double error, void *ptr){
	result res;
	res.flg = flg;
	res.value = value;
	res.error = error;
	res.ptr   = ptr;
	return res;}

/* -------------------------------------------------------------------------- */

public result get_module_version(void){
	printf("Ttile 	: _matrix_\n");
	printf("Made by : Hongkyun\n");
	printf("version : %s/n",_MATRIX_VERSION);
	return construct(0,0,0,_MATRIX_VERSION);}

public void print_vector(Vector *src){
	printf("[");
	for(unsigned int i=0;i<src->dim;i++){
		printf("%*.6lf",12,*(src->element+i));}
	printf("\t]\n");}
public void print_col(Vector *src){
	for(unsigned int i=0;i<src->dim;i++){
		printf("\t| %*.6lf |\n",12,*(src->element+i));}
	printf("\n");}
public void print_solution(Solution *x, int mode){
	if(mode){print_col(x->solution);}
	else{print_vector(x->solution);}}
public void print_matrix(Vector **matrix, const unsigned int row){
	for(unsigned int i=0;i<row;i++){print_vector(matrix[i]);}}
public void print_sys(Vector *src, Vector *base,unsigned int line_no){
	printf("[");
	for(unsigned int i=0;i<src->dim;i++){
		printf("%*.6lf",12,*(src->element+i));}
	printf("\t] [x_%d] ",line_no);
	printf("[%*.6lf\t]\n",12,*(base->element+line_no));}
	
/* set vector dimension */	
public result set_dim(Vector *x, unsigned int dim){
	if(dim>x->max){return construct(1,0,0,"dim <= max");}
	x->dim = dim;return construct(0,0,0,NULL);}
public result restore_max(Vector *x){
	x->dim = x->max;return construct(0,0,0,NULL);}	
	
/* release memory */
public void erase_vector(Vector *src){
	if(src->flag){free(src->element);}
	free(src);}
public void erase_matrix(Vector **M, const unsigned int row){
	for(unsigned int i=0;i<row;i++){erase_vector(*(M+i));}}
/* change one ingredient in vector */
public void change(Vector *src, const unsigned int index, const double value){
	if(index < src->dim){*(src->element+index) = value;}}


/* create vector from static allocated double array */
public Vector *vector_create(unsigned int dim, double *elem){
	if(dim==0){return NULL;}
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= dim;
	vec->dim		= dim;
	vec->element	= elem;
	vec->flag		= 0;
	return vec;}
/* create vector from dynamic allocated double array */	
public Vector *vector_array(unsigned int dim, double *elem){
	if(dim==0){return NULL;}
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= dim;
	vec->dim		= dim;
	vec->element	= elem;
	vec->flag		= 1;
	return vec;}
/* create O-vector */
public Vector *zero_vector(unsigned int dim){
	if(dim==0){return NULL;}
	double *elem	= (double *)calloc(dim,sizeof(double));
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= dim;
	vec->dim		= dim;
	vec->element	= elem;
	vec->flag	 	= 1;
	return vec;}
/* create 1-vector */
public Vector *one_vector(unsigned int dim){
	if(dim==0){return NULL;}
	double *elem	= (double *)calloc(dim,sizeof(double));
	for(unsigned int i=0;i<dim;i++){*(elem+i) = 1;}
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= dim;	
	vec->dim		= dim;
	vec->element	= elem;
	vec->flag	 	= 1;
	return vec;}
/* create index vector */
public Vector *index_vector(unsigned int dim){
	if(dim==0){return NULL;}
	double *elem	= (double *)calloc(dim,sizeof(double));
	for(unsigned int i=0;i<dim;i++){*(elem+i) = i;}
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= dim;
	vec->dim		= dim;
	vec->element	= elem;
	vec->flag	 	= 1;
	return vec;}
/* make column vector from row vector organzied matrix nxn */
public Vector *column_vector(Vector **_Matrix, unsigned int col){
	unsigned int dim= _Matrix[0]->dim;
	double *elem	= (double *)calloc(dim,sizeof(double));
	for(unsigned int i=0;i<dim;i++){*(elem+i) = (_Matrix[i]->element)[col];}
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= dim;
	vec->dim		= dim;
	vec->element	= elem;
	vec->flag	 	= 1;
	return vec;}
/* make column vector from row vector organzied matrix mxn */
public Vector *column_vector_mn(Vector **_Matrix, unsigned int row, unsigned int col){
	double *elem	= (double *)calloc(row,sizeof(double));
	for(unsigned int i=0;i<row;i++){*(elem+i) = (_Matrix[i]->element)[col];}
	Vector *vec 	= (Vector *)malloc(sizeof(Vector));
	vec->max		= row;
	vec->dim		= row;
	vec->element	= elem;
	vec->flag	 	= 1;
	return vec;}
/* believe programmer */
public result fast_column_vector(Vector **M, Vector *v, unsigned int row, unsigned int col){
	if(v->dim!=row){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<row;i++){*(v->element+i) = (M[i]->element)[col];}
	return construct(0,0,0,NULL);}
public result get_value(Vector *src,const unsigned int index){
	if(index < src->dim){return construct(0,*(src->element+index),0,NULL);}
	return construct(1,0,0,"Dimension match error");}
/* swap two ingredients in one vector --> for exterior using */
public result swap(Vector *src, const unsigned int index_o, const unsigned int index_t){
	if(src->dim <= index_o || src->dim <= index_t){return construct(1,0,0,"Dimension match error");}
	double tmp = *(src->element+index_o);
	*(src->element+index_o) = *(src->element+index_t);
	*(src->element+index_t) = tmp;
	return construct(0,0,0,NULL);}
/* swap two ingredients in one vector --> for interior function */
private void __v_swap(Vector *src, const unsigned int index_o, const unsigned int index_t){
	double tmp = *(src->element+index_o);
	*(src->element+index_o) = *(src->element+index_t);
	*(src->element+index_t) = tmp;}
/* interchange position between two vectors */
public result interchange(Vector *src, Vector *dst){
	if(src->dim != dst->dim){return construct(1,0,0,"Dimension match error");}
	double *tmp  = src->element;
	int flag	 = dst->flag;
	src->element = dst->element;
	dst->element = tmp;
	dst->flag	 = src->flag;
	src->flag	 = flag;
	return construct(0,0,0,NULL);}
private void __interchange(Vector *src, Vector *dst){
	double *tmp  = src->element;
	int flag	 = dst->flag;
	src->element = dst->element;
	dst->element = tmp;
	dst->flag	 = src->flag;
	src->flag	 = flag;}
/* vector pointer swap and original vector release btw same size vectors */
public result vector_move(Vector *src, Vector *dst){
	if(src->dim == dst->dim){
		__interchange(src,dst);
		erase_vector(src);
		src = NULL;
		return construct(0,0,0,NULL);}
	return construct(1,0,0,"Dimension match error");}
/* vector copy --> for exterior using */
public result vector_copy(Vector *src, Vector *dst){
	if(src->dim == dst->dim){
		for(unsigned int i=0;i<src->dim;i++){
			*(dst->element+i) = *(src->element+i);}
		return construct(0,0,0,NULL);}
	return construct(1,0,0,"Dimension match error");}
/* vector copy --> for interior function */
private void __vector_copy(Vector *src, Vector *dst){
	for(unsigned int i=0;i<src->dim;i++){
		*(dst->element+i) = *(src->element+i);}}
/* vector ingredient shift : if dir < 0 : shift left else shift right */
public result vector_shift(Vector *src, int direction){
	Vector *temp = zero_vector(src->dim);
	if(direction<0){*(temp+(src->dim-1))=*(src+(src->dim-1));}
	else{*(temp)=*(src);}
	memmove(src, src+direction, sizeof(src) -sizeof(*src));	
	__vector_copy(src,temp);
	erase_vector(src);
	src = temp;
	return construct(0,0,0,NULL);	}
public result trace(Vector **Matrix, Vector *diagonal){
	Vector *temp = zero_vector(Matrix[0]->dim);
	for(unsigned int i=0;i<Matrix[0]->dim;i++){change(temp,i,*(Matrix[i]->element+i));}
	interchange(temp,diagonal);
	erase_vector(temp);
	return construct(0,0,0,NULL);}
public result transpose(Vector **original, Vector **tmatrix, unsigned int row){
	if(tmatrix[0]->dim!=row){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<original[0]->dim;i++){
		fast_column_vector(original,tmatrix[i], row, i);
	}return construct(0,0,0,NULL);}
/* vector base arithmetic */
/* vector addition/subtraction */
public result add(Vector *src, Vector *dst, Vector *__new__, const char op){
	if(src->dim != dst->dim){return construct(1,0,0,"Dimension match error");}
	else if(src->dim != __new__->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<src->dim;i++){
		if(op=='+'){
			*(__new__->element+i) = *(src->element+i) + *(dst->element+i);}
		else{*(__new__->element+i) = *(src->element+i) - *(dst->element+i);}
	}return construct(0,0,0,NULL);}
/* vector inner product --> for exterior using  */
public result dot(Vector *src, Vector *dst){
	if(src->dim != dst->dim){return construct(1,0,0,"Dimension match error");}
	register double summation = 0;
	for(unsigned int i=0;i<src->dim;i++){
		summation = summation + (*(src->element+i))*(*(dst->element+i));
	}return construct(0,summation,0,NULL);
}
/* vector inner product --> for interior using  */
private result __v_dot(Vector *src, Vector *dst){
	register double summation = 0;
	for(unsigned int i=0;i<src->dim;i++){
		summation = summation + (*(src->element+i))*(*(dst->element+i));
	}return construct(0,summation,0,NULL);}
/* vector constant product  */
public result mult(Vector *src, Vector *__new__, const double c, unsigned int start){
	if(src->dim != __new__->dim){return construct(1,0,0,"Dimension match error");}
	else if(start >= src->dim){return construct(1,0,0,"Dimension match error");}
	for(;start<src->dim;start++){
		*(__new__->element+start) = c * (*(src->element+start));
	}return construct(0,0,0,NULL);}
/* summation of ingredient in src vector */
public double ingredient_sum(Vector *src){
	register double temp=0;
	for(unsigned int i=0;i<src->dim;i++){temp+=*(src->element+i);}
	return temp;}
/* vector addition and subtraction */
public result vector_add(Vector *src, Vector *dst, int op, unsigned int start, unsigned int end){
	if(src->dim!=dst->dim){return construct(1,0,0,"Dimension match error");}
	else if(start >= src->dim){return construct(1,0,0,"Dimension match error");}
	else if(end   > src->dim){return construct(1,0,0,"Dimension match error");}
	else if(start>=end){return construct(1,0,0,"Boundary error");}
	op = op < 0 ? -1 : 1;	// set arithmetic sign
	for(;start<end;start++){
		*(dst->element+start) = *(dst->element+start)+op*(*(src->element+start));
	}return construct(0,0,0,NULL);}
/*  multiplicate ingredients between two vectors */
public result vector_mult(Vector *src, Vector *dst, unsigned int start, unsigned int end){
	if(src->dim!=dst->dim){return construct(1,0,0,"Dimension match error");}
	else if(start >= src->dim){return construct(1,0,0,"Dimension match error");}
	else if(end   > src->dim){return construct(1,0,0,"Dimension match error");}
	else if(start>=end){return construct(1,0,0,"Boundary error");}
	for(;start<end;start++){
		*(dst->element+start) = *(dst->element+start) * (*(src->element+start));
	}return construct(0,0,0,NULL);}
/* divide ingredients between two vectors */
public result vector_divide(Vector *src, Vector *dst, unsigned int start, unsigned int end,unsigned int mode){
	if(src->dim!=dst->dim){return construct(1,0,0,"Dimension match error");}
	else if(start >= src->dim){return construct(1,0,0,"Dimension match error");}
	else if(end   > src->dim){return construct(1,0,0,"Dimension match error");}
	else if(start>=end){return construct(1,0,0,"Boundary error");}
	for(;start<end;start++){
		if(mode){
			*(dst->element+start) = *(src->element+start)/(*(dst->element+start));
			if(*(dst->element+start)==0){return construct(1,0,0,"Division by zero");}}
		else{
			*(dst->element+start) = *(dst->element+start)/(*(src->element+start));
			if(*(src->element+start)==0){return construct(1,0,0,"Division by zero");}}
	}return construct(0,0,0,NULL);}
public result point_add(Vector *src, unsigned int inx, double pt){
	if(src->dim>inx){*(src->element+inx)+=pt;}
	else{return construct(1,0,0,"Dimension match error");}
	return construct(0,0,0,NULL);}
/* vector = func(vector) */
public result vector_conversion(Vector *src, unsigned int func){
	mfts				= *mftfind(func);
	for(unsigned int i=0;i<src->dim;i++){
		*(src->element+i)	= mfts.eval(*(src->element+i));}
	return construct(0,0,0,NULL);}
/* return absolutly max ingredient of vector */
public result find_max(Vector *src, unsigned int start){
	if(start >= src->dim){return construct(1,0,0,"Dimension match error");}
	double max = fabs(*(src->element+start));
	unsigned int index = start;
	for(;start < src->dim;start++){
		if(max<fabs(*(src->element+start))){
			max = fabs(*(src->element+start));
			index = start;}
	}return construct(index,max,0,NULL);}
/* return absolutly max row ingredient of matrix */
public result find_max_col(Vector **src, unsigned int row, unsigned int start, unsigned int col){
	if(start >= src[0]->dim){return construct(-1,0,0,"Dimension match error");}
	Vector *temp = column_vector_mn(src,row,col);
	double max = fabs(temp->element[start]);int index = start;
	for(;start < temp->dim;start++){
		if(max<fabs(*(temp->element+start))){
			max = fabs(*(temp->element+start));
			index = start;}
	}erase_vector(temp);
	return construct(index,max,0,NULL);}
/* get vector dimension */
public unsigned int dims(Vector *vec){return vec->dim;}
/* compare vector ingredient whether the difference of each elements is in the boundary */
public result vector_cmp(Vector *src, Vector *dst, const double bound, unsigned int start, unsigned int end){
	if(src->dim!=dst->dim){return construct(-1,0,0,"Dimension match error");}
	else if(start>=end){return construct(-1,0,0,"Start point must be smaller than end");}
	else if(start < src->dim && end <= src->dim){
		int flag = 0;
		for(unsigned int i=start;i<end;i++){
			flag = (fabs(*(src->element+i)-*(dst->element+i)) <= bound) ? 1:0 ;}
	return construct(0,flag,0,NULL);}
	return construct(-1,0,0,"Dimension match error");}
/* size of the vector */
public double vector_length(Vector *src){
	register double summation = 0;
	for(unsigned int i=0;i<src->dim;i++){
		summation += pow(*(src->element+i),2);}
	return sqrt(summation);}
/* zero matrix */
public void matrix_zeros(Vector **M, const unsigned int row, const unsigned int col){
	for(unsigned int i=0;i<row;i++){M[i] = zero_vector(col);}}
/* eigen matrix */
public void matrix_eigen(Vector **M, const unsigned int row){
	for(unsigned int i=0;i<row;i++){
		M[i] = zero_vector(row);
		*(M[i]->element+i)=1;}}
/* from Ax = b , get b*/
public result matrix_single_mult(Vector **M, Vector *x, Vector *b, unsigned int row){
	if(M[0]->dim!=x->dim){return construct(1,0,0,"Dimension match error");}
	else if(b->dim!=row){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<row;i++){*(b->element+i) = __v_dot(*(M+i),x).value;}
	return construct(0,0,0,NULL);}
/* AB = M */
public result matrix_mult(Vector **A, Vector **B, Vector **M, unsigned int rowA, unsigned int rowB){
	if(A[0]->dim!=rowB){return construct(1,0,0,"Dimension match error");}
	else if(M[0]->dim!=rowA){return construct(1,0,0,"Dimension match error");}
	Vector *temp = zero_vector(rowB);
	for(unsigned int i=0;i<B[0]->dim;i++){
		fast_column_vector(B,temp,rowB,i);
		for(unsigned int j=0;j<rowA;j++){*(M[j]->element+i) = __v_dot(*(A+j),temp).value;}
	}erase_vector(temp);
	return construct(0,0,0,NULL); }
/* A[j][i] = A[i][j] */ 
public result symmestry(Vector **A, unsigned int row){
	if(A[0]->dim!=row){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<row;i++){
		for(unsigned int j=i+1;j<row;j++){*(A[j]->element+i)=*(A[i]->element+j);}}
	return construct(0,0,0,NULL);}
/*-------------------------------------------------------------------------- */
/* Gauss elimination */

/* elimination rows */
public result elimination(Vector *ori, Vector *tgt, unsigned int col){
	if(ori->dim < col || tgt->dim < col){return construct(1,0,0,"Dimension match error");}
	else if(fabs(ori->element[col])< etha){return construct(2,0,0,NULL);}
	double m = tgt->element[col] / ori->element[col];
	for(;col<ori->dim;col++){
		*(tgt->element+col) = *(tgt->element+col) - m*(*(ori->element+col));
	}return construct(0,m,0,NULL);}
/* elimination base vector */
public result elimination_base(Vector *base, double m, unsigned int row0, unsigned int row1){
	if(base->dim < row0 || base->dim < row1 ){return construct(1,0,0,"Dimension match error");}
	*(base->element+row1) = *(base->element+row1) - m*(*(base->element+row0));
	return construct(0,m,0,NULL);}
/* make the value 1 which are posed on the first non-zero point in the vector */
public void vector_head_normalize(Vector *ori){
	double norm;int col = 0;
	while(*(ori->element+col)==0){col++;}
	norm = *(ori->element+col);
	for(unsigned int i = col;i<ori->dim;i++){
		*(ori->element+i) = *(ori->element+i)/norm;}}
public Vector *vector_string(unsigned int dim, char *expr, size_t size){
	Vector *init_vec= zero_vector(dim);
	char *splited 	= NULL;
	char *ptr 		= NULL;
	char *cpyd		= (char*)malloc(size*sizeof(char));
	double data;unsigned int i=0;
	memcpy(cpyd,expr,size);
	#ifdef _MSC_VER								/* For visual studio */
	char *ptrs		= NULL;
	splited = strtok_s(cpyd," ,",&ptr);
	while(splited){
		if(i==dim){break;}
		data 	= strtod(splited,&ptrs);
		change(init_vec,i++,data);
		splited = strtok_s(NULL," ,",&ptr);}	
	#else 										/* For other compiler */
	splited = strtok(cpyd," ,");
	while(splited){
		if(i==dim){break;}
		data 	= strtod(splited,&ptr);
		change(init_vec,i++,data);
		splited = strtok(NULL," ,");}
	#endif
	free(cpyd);
	return init_vec;}
/*-------------------------------------------------------------------------- */

/* consider "x" vector in Ax = b*/
public Solution *solution_set(unsigned int dim){
	if(dim==0){return NULL;}
	Solution *x = (Solution *)malloc(sizeof(Solution));
	Vector *inx = index_vector(dim);
	Vector *sol = zero_vector(dim);
	x->dim		= dim;
	x->index 	= inx;
	x->solution = sol;
	return x;}
public void solution_reset(Solution *x, const unsigned int dim){
	Vector *inx = index_vector(dim);
	Vector *sol = zero_vector(dim);
	__vector_copy(inx,x->index);
	__vector_copy(sol,x->solution);
	erase_vector(inx);
	erase_vector(sol);}
public void solution_value_set(Solution *x, const unsigned int index, const double data){
	if(x->dim > index){*(x->solution->element+index) = data;}
	else{return;}	}
/* initialize elements as src vector */
public void solution_vector_set(Vector *src, Solution *x){
	if(x->dim == src->dim){__vector_copy(src,x->solution);}}
/* exchange indices and elements */
public void swap_solution_all(Solution *x,const unsigned int index_o, const unsigned int index_t){
	__v_swap(x->index,index_o,index_t);
	__v_swap(x->solution,index_o,index_t);}
/* exchange indices only */
public void swap_solution(Solution *x,const unsigned int index_o, const unsigned int index_t){
	__v_swap(x->index,index_o,index_t);}
public void erase_solution(Solution *s){
	erase_vector(s->index);
	erase_vector(s->solution);
	free(s);}
/* copy solution vector to another vector. */
public void solution_vector(Solution *x, Vector *target){
	__vector_copy(x->solution,target);}
/* create a vector which manipulates solution vector */
public void solution_handle(Solution *x, Vector *target){
	erase_vector(target);
	target = x->solution;}
/*-------------------------------------------------------------------------- */

/* Pivot : elimination ready process */
/* row pivoting */
public result row_pivoting(Vector **M, Solution *x, Vector *base, unsigned int row){
	double tmp; double max; unsigned int indicator = 0;
	if(row != M[0]->dim && row != x->solution->dim && row != base->dim){
		return construct(3,0,0,NULL);}
	for(unsigned int i=0;i<row;i++){		
		tmp 	  = find_max(M[i],i).value;	/* find maximun ingredient in a row for scale */
		max 	  = tmp;					/* initialized max */
		indicator = i;						/* swap dest */
		for(unsigned int j=i+1;j<row;j++){	/* upper triangle matrix */
			tmp = find_max(M[j],i).value;
			if(tmp>max){
				indicator = j;
				max = tmp;}}
		if(i!=indicator){
			__interchange(M[i],M[indicator]);	/* row exchange */
			swap_solution(x,i,indicator);
			__v_swap(base,i,indicator);}		
	}return construct(0,0,0,NULL);}
/* column pivoting */
public result col_pivoting(Vector **M, Solution *x, Vector *base, unsigned int row){
	if(row != M[0]->dim && row != x->solution->dim && row != base->dim){
		return construct(3,0,0,NULL);}
	unsigned int index;
	for(unsigned int i=0;i<row && i < M[0]->dim;i++){
		index  = find_max_col(M,row,i,i).flg;	/* find maximun ingredient in a row for scale */
		if(i!=index){
			__interchange(M[i],M[index]);		/* row exchange */
			swap_solution(x,i,index);
			__v_swap(base,i,index);}		
	}
	if(row != M[0]->dim){return construct(1,0,0,"Linearly dependent");}
	return construct(0,0,0,NULL);}
/* scaled pivoting */
public result pivoting(Vector **M, Solution *x, Vector *base, unsigned int row){
	if(row != M[0]->dim && row != x->solution->dim && row != base->dim){
		return construct(3,0,0,NULL);}
	result res;double rows;double tmp;unsigned int index = 0;
	Vector *scale_vec = zero_vector(M[0]->dim);				 /* scale vector */
	for(unsigned int i=0;i<row;i++){
		for(unsigned int j=i;j<row;j++){
			rows  					= find_max(M[j],j).value;/* row_max */
			tmp  					= *(M[j]->element+j);	 /* get current col value */
			*(scale_vec->element+j)	= tmp/rows;				 /* save at scale vector */
		}res = find_max(scale_vec,i);						 /* find max ing' in scale vector */
		index = res.flg;									 /* get index */
		if(i!=index){
			__interchange(M[i],M[index]);						 /* row exchange */
			swap_solution(x,i,index);						 /* solution index change */
			__v_swap(base,i,index);}							 /* base exchange */
	}erase_vector(scale_vec);
	return construct(0,0,0,NULL);}
/* 후치대입 */
public result calculate_one_line(Vector *row, Vector *base, Solution *x, unsigned int row_no){
	if(row_no == row->dim-1){
		if(fabs(*(row->element+row_no))<=etha){return construct(3,0,0,"Division by near zero");}	/* underflow */
		*(x->solution->element+row_no) = *(base->element+row_no)/(*(row->element+row_no));
		return construct(0,0,0,NULL);}
	else{
		/*  used vector's inner product ; Mij = sum(Ai*Bj)*
			Solution vector is initialized 0 such that above ingredients have no effects on current term 
		*/
		result res = __v_dot(row,x->solution);
		if(res.flg){return construct(3,0,0,res.ptr);}
		else if(fabs(*(row->element+row_no)) < etha){return construct(10,0,0,"Linearly dependent");}
		*(x->solution->element+row_no) = (*(base->element+row_no)-res.value)/(*(row->element+row_no));
		return res;}}
/* 전치대입 */
public result pre_calculate_one_line(Vector *row, Vector *base, Solution *x, unsigned int row_no){
	if(row_no == 0){
		if(fabs(*(row->element+row_no))<=etha){return construct(3,0,0,"Division by near zero");}	/* overflow */
		*(x->solution->element+row_no) = *(base->element+row_no)/(*(row->element+row_no));
		return construct(0,0,0,NULL);}
	else{
		/*  used vector's inner product ; Mij = sum(Ai*Bj)*
			Solution vector is initialized 0 such that above ingredients have no effects on current term 
		*/
		result res = __v_dot(row,x->solution);
		if(res.flg){return construct(3,0,0,res.ptr);}
		else if(fabs(*(row->element+row_no)) < etha){return construct(10,0,0,"Linearly dependent");}
		*(x->solution->element+row_no) = (*(base->element+row_no)-res.value)/(*(row->element+row_no));
		return res;}}
/* 순수 가우스 소거법 */
public result gauss_elimination(Vector **M, Vector *base, Solution *x){
	result res;
	register const unsigned int dim = M[0]->dim;
	if(dim != x->solution->dim && dim != base->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=i+1;j<dim;j++){
			res 	  = elimination(M[i],M[j],i);
			elimination_base(base, res.value, i, j);}
	}for (int i=dim-1;i>=0;i--){calculate_one_line(M[i],base,x,i);}
	return construct(0,0,0,NULL);}
/* 인덱스 벡터를 이용해 행렬을 재정렬함 */
public result arrange_matrix_by_index(Vector **M, Vector *base, Solution *x){
	register const unsigned int dim = M[0]->dim;
	if(dim != x->solution->dim && dim != base->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=i+1;j<dim;j++){
			if(*(x->index->element+j)==i){
				__interchange(M[i],M[j]);
				__v_swap(x->index,i,j);
				__v_swap(base,i,j);}}
	}return construct(0,0,0,NULL);
}
/*-------------------------------------------------------------------------- */
/* LU decomposition */
/* U 행렬 초기화 - 1행 계산 대각 성분 1로 초기화 U행렬은 반드시 O행렬로 초기화 되어있어야 함  */
public result U_init(Vector **ori, Vector **L, Vector **A){
	if(ori[0]->dim!=L[0]->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<ori[0]->dim;i++){
		*(ori[0]->element+i) = *(A[0]->element+i)/(*(L[0]->element));
		*(ori[i]->element+i) = 1;}
	return construct(0,0,0,NULL);}
/* L 행렬 초기화 - 1열 계산 L행렬은 반드시 0행렬로 초기화되어있어야 함 */
public result L_init(Vector **ori, Vector **A){
	if(ori[0]->dim!=A[0]->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<ori[0]->dim;i++){*(ori[i]->element) = *(A[i]->element);}
	return construct(0,0,0,NULL);}
/* U 행렬 성분 계산 -- 코드의 계산원리는 전치대입/후치대입 원리와 동일(계산식은 U분해 알고리즘) */
private void __calculate_U_component(Vector **L, Vector **U, Vector **A, const unsigned int row, const unsigned int col){
	Vector *UCol = column_vector(U,col);
	double res = __v_dot(L[row],UCol).value;
	*(U[row]->element+col) = (*(A[row]->element+col) - res) / (*(L[row]->element+row));
	erase_vector(UCol);}
/* L 행렬 성분 계산 -- 코드의 계산원리는 전치대입/후치대입 원리와 동일(계산식은 L분해 알고리즘) */
private void __calculate_L_component(Vector **L, Vector **U, Vector **A, const unsigned int row, const unsigned int col){
	Vector *UCol = column_vector(U,col);
	double res = __v_dot(L[row],UCol).value;
	*(L[row]->element+col) = (*(A[row]->element+col)) - res;
	erase_vector(UCol);}
/* LU decomposition function */
public result LU_decompose(Vector **L, Vector **U, Vector **A){
	L_init(L,A);
	U_init(U,L,A);
	unsigned int ROW = L[0]->dim;
	if(L[0]->dim != U[0]->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=1;i<ROW;i++){
		for(unsigned int j=1;j<ROW;j++){
			__calculate_L_component(L, U, A, j, i);
			if(j>i){__calculate_U_component(L, U, A, i, j);}	/* upper triangle matrix */
	}}return construct(0,0,0,NULL);}
/*-------------------------------------------------------------------------- */
/* Jacobi repeat method */
public result Jacobi_single_step(Vector **M, Vector *base, Solution *x){
	if(M[0]->dim != x->dim){return construct(1,0,0,"Dimension match error");}
	else if(base->dim != x->dim){return construct(1,0,0,"Dimension match error");}
	double *tmp = (double *)calloc(x->dim,sizeof(double));
	for(unsigned int i=0;i<x->dim;i++){
		tmp[i] = __v_dot(M[i],x->solution).value - *(x->solution->element+i) * (*(M[i]->element+i));}
	for(unsigned int i=0;i<x->dim;i++){ /* update next "terms" */
		if(*(M[i]->element+i)==0){return construct(1,0,0,"Linearly dependent");};
		*(x->solution->element+i) = (*(base->element+i)-tmp[i])/(*(M[i]->element+i));}
	free(tmp);
	return construct(0,0,0,NULL);}
/* Gauss Seidel repeat method */
public result Gauss_Seidel_single_step(Vector **M, Vector *base, Solution *x){
	if(M[0]->dim != x->dim){return construct(1,0,0,"Dimension match error");}
	else if(base->dim != x->dim){return construct(1,0,0,"Dimension match error");}
	double res;
	for(unsigned int i=0;i<x->dim;i++){
		if(*(M[i]->element+i)==0){return construct(1,0,0,"Linearly dependent");};
		res = __v_dot(M[i],x->solution).value;
		res = res - *(x->solution->element+i) * (*(M[i]->element+i));
		*(x->solution->element+i) = (*(base->element+i)-res)/(*(M[i]->element+i));/* update next "term"*/
	}return construct(0,0,0,NULL);}
/* discrite diagonally dominant |aii| > |sum(aij)-aii| 
*  flag = 1 : diagonally dominant 		*/
public int Diagonally_Dominant(Vector **M,unsigned int row){
	int flag=0;double sum=0;
	for(unsigned int i=0;i<row;i++){
		sum = 0;
		for(unsigned int j=0;j<M[0]->dim;j++){sum += *(M[i]->element+j);}
		flag = fabs(*(M[i]->element+i)) > fabs(sum - *(M[i]->element+i)) ? 1:0;
	}return flag;}
/*-------------------------------------------------------------------------- */
/* newton interpolation single step */
private void __newton_interpolation_single_step(Vector **M, Vector *data_x, unsigned int col, unsigned int lv){
	if(lv == 0){return;}
	Vector *cols = column_vector(M,lv-1);
	double *d_x  = data_x->element;
	for(unsigned int st_pt=0,i=lv;st_pt+lv<col;i++){
		change(M[i],lv,(get_value(cols,i).value-get_value(cols,i-1).value)/(d_x[st_pt+lv]-d_x[st_pt]));st_pt++;
	}}
/* newton interpolation (pre-) */
public result newton_interpolation(Vector **Matrix, Vector *data_x, Vector *data_y){
	if(data_x->dim!=data_y->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<data_x->dim;i++){change(Matrix[i],0,*(data_y->element+i));}
	for(unsigned int level = 0; level<data_x->dim; level++){
	__newton_interpolation_single_step(Matrix, data_x, data_x->dim, level);}
	return construct(0,0,0,NULL);}
/* newton interpolation (post-) */
public result post_newton_interpolation(Vector *Newtons, Vector *x, Vector *y){
	if(x->dim!=y->dim || Newtons->dim!=x->dim){return construct(1,0,0,"Dimension match error");}
	Vector *operand = zero_vector(x->dim);vector_copy(y,operand);
	for(unsigned int i=0, st=1;i<x->dim;st++,i++){
		for(unsigned int j=0,k=0;j<x->dim-1;k++,j++){
			*(operand->element+j) = (*(operand->element+j+1)-*(operand->element+j)) /
			(*(x->element+k+st)-*(x->element+k));}
		*(Newtons->element+i+1) = *(operand->element+x->dim-i-2);}
	*(Newtons->element)=*(operand->element+x->dim-1);
	erase_vector(operand);return construct(0,0,0,NULL);}
/* value of newton interpolation function at x */
public result newton_interpolation_calculate(Vector *data_x, Vector *result, double value){
	double sum = 0;double mul = 1;
	for(unsigned int i = 0;i<data_x->dim;i++){
		if(i!=0){mul = mul * (value-*(data_x->element+i));}	// i == 0 : constant, continuos multiplication
		sum = sum + mul * *(result->element+i);
	}return construct(0,sum,0,NULL);}
/*-------------------------------------------------------------------------- */
/* Quadric spline */
private result __next_d(Vector *d, Vector *x, Vector *y, unsigned int current){
	double curd = *(d->element+current);
	double curx = *(x->element+current);
	double cury = *(y->element+current);
	double nxtx = *(x->element+current+1);
	double nxty = *(y->element+current+1);
	if(nxtx==curx){return construct(1,0,0,"Data point cannot be not overwritten");}
	*(d->element+(current+1)) = -curd+2*(nxty-cury)/(nxtx-curx);
	return construct(0,0,0,NULL);}
private void __quadric_spline_coeff(Vector **q, Vector *d, Vector *x, Vector *y){
	for(unsigned int i=0;i<x->dim;i++){
		*(q[0]->element+i) = *(y->element+i);
		*(q[1]->element+i) = *(d->element+i);
		*(q[2]->element+i) = (*(d->element+i+1)-*(d->element+i)) / 
							 (2*(*(x->element+i+1)-*(x->element+i)));}}
private result __calculate_quadric_spline(Vector **q, Vector *x, double value, int mode){
	static int last_point = 0;double temp;	// persist last_point for drawing graph
	result res = construct(0,0,0,NULL);
	if(mode){last_point=0;}					// if mode == 1 : reinit
	while(value>=*(x->element+last_point)){
		if(last_point==x->dim){break;}
		last_point++;}
	if(last_point){last_point--;}
	temp = value-*(x->element+last_point);
	res.value = *(q[0]->element+last_point) + \
				*(q[1]->element+last_point) * temp + \
				*(q[2]->element+last_point) * temp * temp;
	return res;}
public result quadric_spline(Vector **q, Vector *x, Vector *y){	// O(2n)
	Vector *d = zero_vector(x->dim);
	result res;
	for(unsigned int i=0;i<x->dim-1;i++){
		res = __next_d(d,x,y,i);
		if(res.ptr!= NULL){return res;}}
	__quadric_spline_coeff(q, d, x, y);
	erase_vector(d);
	return construct(0,0,0,NULL);}							// return coefficient matrix
public result quadric_spot(Vector **q, Vector *x, double value){
	return __calculate_quadric_spline(q,x,value,1);}		// return value of function point
public result quadric_graph(Vector **q, Vector *x, Vector *time, Vector *output, int mode){
	if(mode){quadric_spot(q, x, *(time->element));}			// initialized
	result res = construct(0,0,0,NULL);						// quadric_spline -> quadric_graph
	if(time->dim!=output->dim){return construct(1,0,0,"Dimension match error");}
	for(unsigned int i=0;i<time->dim;i++){
		res = __calculate_quadric_spline(q,x,*(time->element+i),0);
		if(res.ptr==NULL){*(output->element+i) = res.value;}
		else{return res;}}
	return construct(0,0,0,NULL);}
/*-------------------------------------------------------------------------- */
/* Cubic spline */	
private void __cubic_coeff_matrix(Vector *cubic[S_CUBIC], Vector *s, Vector *y, Vector *h){
	for(unsigned int i=0;i<s->dim;i++){
		double h1 	= *(h->element+i);
		double c 	= h1*(2*(*(s->element+i))+*(s->element+i+1))/6;
		*(cubic[0]->element+i) = (*(s->element+i+1)-*(s->element+i))/(6*h1);	// ax^3
		*(cubic[1]->element+i) = *(s->element+i)*0.5;							// bx^2
		*(cubic[2]->element+i) = (*(y->element+i+2)-*(y->element+i+1)/h1-c);	// cx
		*(cubic[3]->element+i) = *(y->element+i+1);}}
private double __cubic_cal_y(double y1, double y2, double y3, double h1, double h2){
	return 6*((y3-y2)/h2-(y2-y1)/h1);}
private void __cubic_h_vector(Vector *x, Vector *h){		// get interval of x
	for(unsigned int i=0;i<dims(x)-1;i++){
		*(h->element+i) = *(x->element+i+1) - *(x->element+i);}
	*(h->element+x->dim-1) = *(x->element+x->dim-1) - *(x->element+x->dim-2);}
private double __cubic_initiate(Vector *s, Vector *h,int mode){
	if(mode==1){return 0;}
	else if(mode==2){return *(s->element);}
	else{
		double h1 = *(h->element);
		double h2 = *(h->element+1);
		return ((h1+h2)*(*(s->element))-h1*(*(s->element+1)))/h2;}}
private void __cubic_evaluate(Vector **cubic, Vector **coeff, Vector *base, Solution *s){
	unsigned int limit = base->dim-2;
	Vector *temp = zero_vector(limit);
	for(unsigned int i=0;i<limit;i++){*(temp->element+i)=*(base->element+i);}
	for(unsigned int i=0;i<limit;i++){
		*(cubic[i]->element+i) = *(coeff[1]->element+i);
		if(i!=0){*(cubic[i]->element+i-1)=*(coeff[0]->element+i);}
		*(cubic[i]->element+i+1) = *(coeff[2]->element+i);}
	gauss_elimination(cubic,temp,s);
	erase_vector(temp);}
private result __cubic_CAL_coeffs_matrix(Vector **c, Vector *base, Vector *h, Vector *y,int mode){
	double temp = 0;
	for(unsigned int i=0;i<h->dim-1;i++){
		temp = __cubic_cal_y(*(y->element+i),*(y->element+i+1), \
							*(y->element+i+2),*(h->element+i),*(h->element+i+1));
		*(c[0]->element+i) = *(h->element+i+1);
		*(c[1]->element+i) = 2*(*(h->element+i)+*(h->element+i+1));
		*(c[2]->element+i) = *(h->element+i+1);
		*(base->element+i) = temp;}
	*(c[0]->element) = 0;*(c[2]->element+h->dim-2) = 0;
	if(mode==1){return construct(0,0,0,"Natural spline");}
	else if(mode==2){
		*(c[1]->element) = 3*(*(h->element))+2*(*(h->element+1));
		*(c[1]->element+h->dim-2) = 3*(*(h->element+h->dim-2))+2*(*(h->element+h->dim-1));
		return construct(0,0,0,"Clamped spline");}
	else{
		double h1 = *(h->element);
		double h2 = *(h->element+1);
		*(c[1]->element) = (h1+h2)*(h1+2*h2)/h2;
		*(c[2]->element) = (h2*h2-h1*h1)/h2;
		h1 = *(h->element+h->dim-2);
		h2 = *(h->element+h->dim-1);
		*(c[1]->element+h->dim-2) = (h1+h2)*(2*h1+h2)/h1;
		*(c[0]->element+h->dim-2) = (h1*h1-h2*h2)/h1;}
	return construct(0,0,0,"Not-a-knot spline.");}
public result cubic_spot(Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Vector *s, double value, int mode){
	static int last_point = 0;double temp = 0;
	if(mode){last_point = 0;}				// mode unset -> continuous mode set -> search
	result res = construct(0,0,0,NULL);
	while(value>=*(x->element+last_point)){
		if(last_point++==x->dim-2){break;}} 
	if(last_point>1){
		last_point--;
		temp = value-*(x->element+last_point);
		last_point--;
		res.value =  *(coeff[3]->element+last_point) + \
					(*(coeff[2]->element+last_point)) * temp + \
					(*(coeff[1]->element+last_point)) * temp * temp + \
					(*(coeff[0]->element+last_point)) * temp * temp * temp;
		return res;}
	else{									// first interval
		double s_0 = __cubic_initiate(s,h,mode);
		double h1  = *(h->element);
		double c1  = h1*(2*s_0+*(s->element))/6;
		temp = value-*(x->element);	
		res.value = *(y->element) + ((*(y->element+1)-*(y->element))/h1-c1) * temp + \
					s_0/2 * temp * temp + (*(s->element)-s_0)/(6*h1) * temp * temp * temp;
		return res;}}	
#ifdef _MSC_VER
public result cubic_spline(Vector **cubic, Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Solution *s, unsigned int dp, int mode) {
	result res;
	erase_vector(s->index);
	erase_vector(s->solution);
	Vector *index = zero_vector(dp-2);
	Vector *newsl = zero_vector(dp-2);
	s->index = index; s->solution = newsl;
	matrix_zeros(coeff,S_CUBIC,dp);
	matrix_zeros(cubic,dp-2,dp-2);
	Vector *base = zero_vector(dp);
	__cubic_h_vector(x, h);
	res = __cubic_CAL_coeffs_matrix(coeff, base, h, y, mode);
	res.flg = dp-2;
	if (res.ptr != NULL) {
		__cubic_evaluate(cubic, coeff, base, s);
		__cubic_coeff_matrix(coeff, s->solution, y, h);
	}erase_vector(base);
	return res;}
#else
public result cubic_spline(Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Solution *s, unsigned int dp, int mode){
	result res;
	erase_vector(s->index);
	erase_vector(s->solution);
	Vector *cubic[dp-2];	
	Vector *index = zero_vector(dp-2);
	Vector *newsl = zero_vector(dp-2);
	s->index = index; s->solution = newsl;
	matrix_zeros(coeff,S_CUBIC,dp);
	matrix_zeros(cubic,dp-2,dp-2);
	Vector *base = zero_vector(dp);
	__cubic_h_vector(x, h);
	res = __cubic_CAL_coeffs_matrix(coeff, base, h, y, mode);
	res.flg = dp-2;
	if(res.ptr!=NULL){
		__cubic_evaluate(cubic,coeff,base,s);
		__cubic_coeff_matrix(coeff, s->solution, y, h);}
	erase_vector(base);erase_matrix(cubic,dp-2);
	return res;}
#endif
public result cublic_spline_pro(Vector **cubic, Vector *coeff[S_CUBIC], Vector *h, Vector *x, Vector *y, Solution *s, unsigned int dp, int mode) {
	result res;
	Vector *base = zero_vector(dp);
	solution_reset(s,dp);
	matrix_zeros(cubic, dp - 2, dp - 2);
	__cubic_h_vector(x, h);
	res = __cubic_CAL_coeffs_matrix(coeff, base, h, y, mode);
	res.flg = dp - 2;
	if (res.ptr != NULL) {
		__cubic_evaluate(cubic, coeff, base, s);
		__cubic_coeff_matrix(coeff, s->solution, y, h);
	}erase_vector(base);
	return res;
}
/*-------------------------------------------------------------------------- */
/* least square matrix */
private result __least_square(Vector **A, Vector *x, Vector *y, Vector *base, unsigned int e){
	Vector *operand 	= zero_vector(x->dim);
	register double tmp;
	register unsigned int inx=0;
	for(unsigned int i=0;i<x->dim;i++){tmp=1;				/* calculate first row */
		for(unsigned int j=1;j<e;j++){						/* save multiplication counts using post multiplication result */
			tmp *= *(x->element+i);
			point_add(A[0],j,tmp);							
			point_add(base,j,tmp * *(y->element+i));	
			if(j<e-1){continue;}
			else{*(operand->element+(inx++))=tmp;}}}		/* Save last point for next */ 
	*(A[0]->element) = x->dim;								/* A[0][0] = n(DATA) */
	*(base->element) = ingredient_sum(y);					/* B[0] = sum(y) 	 */
	for(unsigned int i=1;i<e;i++){
		for(unsigned int j=0;j<e;j++){*(A[i]->element+j) = *(A[i-1]->element+j+1);} // shift 1 left
		vector_mult(x,operand,0,x->dim);					/* save last point for next */
		*(A[i]->element+e-1) = ingredient_sum(operand);		/* set B[i] */
	}erase_vector(operand);
	return construct(0,0,0,NULL);}
/* get least square matrix */
public result least_square_matrix(Vector **A, Vector *x, Vector *y, Solution *base, unsigned int e){
	result res  = construct(0,0,0,NULL);
	res.ptr 	= __least_square(A,x,y,base->solution,e+1).ptr;
	__vector_copy(base->solution,base->index);
	return res;}
/* calculate polynomial least square coefficient (oneshot) */
#ifdef _MSC_VER															/* For visual studio */
public result polyfit(Vector **A, Vector *x, Vector *y, Solution *base, unsigned int e){
	result res = construct(0,0,0,NULL);
	matrix_zeros(A,e+1,e+1);
	Vector *temp = zero_vector(e+1);								/* for reset */
	res = least_square_matrix(A,x,y,base,e);						/* get Least Squares Matrix */
	__vector_copy(temp, base->solution);							/* Solution vector reset for elimination */
	if (res.ptr == NULL){gauss_elimination(A,base->index,base);}	/* Symmestry matrix */
	erase_vector(temp);												/* release heap */
	return res;}
#else
public result polyfit(Vector *x, Vector *y, Solution *base, unsigned int e){
	result res  = construct(0,0,0,NULL);
	Vector *A[e+1];
	matrix_zeros(A,e+1,e+1);
	Vector *temp = zero_vector(e+1);							/* for reset */
	res = least_square_matrix(A, x, y, base,e);					/* get Least Squares Matrix */
	__vector_copy(temp,base->solution);							/* Solution vector reset for elimination */
	if(res.ptr==NULL){gauss_elimination(A, base->index, base);}	/* Symmestry matrix */
	erase_matrix(A,e+1);erase_vector(temp);						/* release heap */
	return res;}
#endif
/* the value of polynomial least square function at point x */
public result polyval(Vector *coeff, double x){
	double sum = 0;
	for(unsigned int i=coeff->dim-1;i>0;i--){sum = sum + *(coeff->element+i)*pow(x,i);}
	sum = sum + *(coeff->element);
	return construct(0,sum,0,NULL);}
/* Goodness” of fitness about the method of least-squares */
public result rsquare(Vector *coeff, Vector *x, Vector *y){
	double SSE = 0;double SST = 0;double value = 0;double mean = 0;
	for(unsigned int i=0;i<y->dim;i++){mean += *(y->element+i);}
	mean = mean / ((double)y->dim);
	for(unsigned int i=0;i<y->dim;i++){
		SSE += pow(*(y->element+i)-polyval(coeff,*(x->element+i)).value,2);
		SST += pow(*(y->element+i)-mean,2);}
	value = SST == 0 ? 0 : fabs(1-SSE/SST);						/* SST = SSE + SSR */
	value = SSE == 0 && SST !=0 ? 1 : value;
	return construct(0,value,SSE,NULL);}
/* get SSE from least-squares */	
public result least_error(Vector *coeff, Vector *x, Vector *y){
	double SSE = 0;
	for(unsigned int i=0;i<y->dim;i++){							/* SSE = sum((y - LM(x))^2) */
		SSE += pow(*(y->element+i)-polyval(coeff,*(x->element+i)).value,2);}
	return construct(0,SSE,0,NULL);}
/*-------------------------------------------------------------------------- */
private double __middle_approximation(double (*func)(double),double x, double h){
	return (func(x+h)-func(x-h))/(2*h);}
private double __middle_approximation_dy(double (*func)(Vector *, Vector *,double, int), Vector *c, Vector *e, int i, double x, double h){
	return (func(c,e, x+h,i)-func(c,e, x-h,i))/(2*h);}	
/* For static coefficient & exponent -- Richardson Extrapolaration */
public result diff_static(double (*func)(double),double x, double h){
	if(h<0){h = etha;}
	double a_2 	= __middle_approximation(func, x, h/2.0);
	double a 	= __middle_approximation(func, x, h);
	return construct(0,a_2 + (1.0/3.0)*(a_2-a),0,NULL);}
/* For dynamic coefficient & exponent */
public result diff_dynamic(double (*func)(Vector *, Vector *, double, int),Vector *c, Vector *e, int i, double x, double h){
	if(h<0){h = etha;}
	double a_2 	= __middle_approximation_dy(func, c, e, i, x, h/2.0);
	double a 	= __middle_approximation_dy(func, c, e, i, x, h);
	return construct(0,a_2 + (1.0/3.0)*(a_2-a),0,NULL);}
/* For numerical data points (value returns) */
public result diff_numeric(Vector *x, Vector *y, unsigned int index){
	if(x->dim!=y->dim){return construct(1,0,0,"Dimension match error");}
	else if(index>1||index+2<x->dim){
		double h	= *(x->element+index+1)-*(x->element+index);
		double a_2 	= *(y->element+index+1)-*(y->element+index-1)/(0.5*h);
		double a 	= *(y->element+index+2)-*(y->element+index-2)/h;
		return construct(0,a_2 + (1.0/3.0)*(a_2-a),0,NULL);}
	return construct(1,0,0,"Condition is not enough");}
/* Newton optimization  */		
public result newton_minimalization(double (*func)(double),double ini, double fin, double h){
	double range[5]={ini,fin,0,0,0};
	if(h<=0){h = etha;}
	do{
		if(fabs(range[0]-range[1])<h){break;}
		range[2] = diff_static(func,range[0],etha).value;	// f'(x_n-1)
		range[3] = diff_static(func,range[1],etha).value;	// f'(x_n)
		range[4] = (range[3]-range[2])/(range[1]-range[0]);	// f"(x)
		range[0] = range[1];								// shift
		range[1] = range[0]-range[3]/range[4];
	}while(1);return construct(0,range[1],0,NULL);}
/*-------------------------------------------------------------------------- */
/* Numerical integration */
/* simple trapezoidal method */
public result trapezoidal(Vector *y, double h){
	double s = 0.5*h*(*(y->element)+*(y->element+y->dim-1));
	for(unsigned int i=1;i<dims(y)-1;i++){s += *(y->element+i)*h;}
	return construct(0,s,0,NULL);}
/* simple simpson method */
public result simpson_quad(Vector *y, double h){
	double s = 0;
	for(unsigned int i=1;i<dims(y)-1;i*=2){
		s += h/3*(*(y->element+i-1)+4*(*(y->element+i))+*(y->element+i+1));}
	return construct(0,s,0,NULL);}	
private void __rombergtrap(Vector *first, double (*func)(double), int level, double set[3]){
	double sum = 0;double store = 1;
	set[0] = 0.5*(set[2]-set[1])*(func(set[2])+func(set[1]));
	for(unsigned int i=0; i<dims(first);i++){
		sum = set[0];
		for(unsigned int j=1;j<=(unsigned int)(store*0.5);j++){
			sum += (set[2]-set[1])/store*(func(set[1]+((set[2]-set[1])*(j-0.5)/(pow(2,i-1)))));}
		change(first,i,sum);
		store*=2;set[0]= sum*0.5;}}
private double __romberginternal(Vector *R, unsigned int rows, unsigned long lv){
	return (*(R->element+rows)*lv - *(R->element+rows-1))/(lv-1);}
/* Romberg's quadnature (Mem usage : n^2) */
public void romberg_matrix(Vector **Romberg, double (*func)(double), unsigned int cols, double set[3]){
	__rombergtrap(Romberg[0],func,0,set);unsigned long lv = 4;
	for(unsigned int i=1; i<cols;i++){
		for(unsigned int j=i; j<dims(Romberg[0]);j++){
			*(Romberg[i]->element+j)=__romberginternal(Romberg[i-1],j,lv);}
		lv*=4;}}
/* get value of Romberg's quadnature */
public result romberg_integrate(Vector **Romberg, unsigned int row, unsigned int col){
	if(Romberg[0]->dim<=row){return construct(1,0,0,"Dimension match error");}
	return construct(0,*(Romberg[col]->element+row),0,NULL);}
/* Romberg's quadnature (Mem usage : n) */
public result integrate(double (*func)(double), unsigned int iter, double a, double b){
	Vector *Romberg = zero_vector(iter);
	double set[3] = {0,a,b};
	__rombergtrap(Romberg,func,0,set);unsigned long lv = 4;
	for(unsigned int i=1;i<iter;i++){
		for(unsigned int j=1;j<iter;j++){
			*(Romberg->element+j-1) = __romberginternal(Romberg,j,lv);}lv*=4;}
	set[0] = *(Romberg->element);
	erase_vector(Romberg);return construct(0,set[0],0,NULL);}
/* solve 1-st diff equation with euler method */	
public result euler_method(double (*func)(double,double),double bot, double top, double init,unsigned int N){
	double h = fabs(top-bot)/(double)N;double y = init; double t = bot;
	for(unsigned int i=0;i<N;t+=h,i++){y = y+h*func(t,y);}
	return construct(0,y,h,NULL);}
/* solve 1-st diff equation with Runge-Kutta lv 4 */	
public result RK4(double (*func)(double,double),double bot, double top, double init, unsigned int N){
	double h 	= fabs(top-bot)/(double)N;double y = init;double t = bot;
	double k[4] = {0,0,0,0};
	for(unsigned int i=0;i<N;t+=h,i++){
		k[0] = func(t,y);
		k[1] = func(t+0.5*h,y+0.5*h*k[0]);
		k[2] = func(t+0.5*h,y+0.5*h*k[1]);
		k[3] = func(t+h,y+h*k[2]);
		y += h/6*(k[0]+2*k[1]+2*k[2]+k[3]);
	}return construct(0,y,h,NULL);}
/* solve 2-nd system of differential equation using Runge-Kutta method order 4 */	
public result RK4_2(double (*func)(double,double,double),double set[4], unsigned int N){
	double h = fabs(set[1]-set[0])/(double)N;
	double y = set[2];double x = set[3];double t = set[0];
	double k[8] = {0,0,0,0,0,0,0,0};
	for(unsigned int i=0;i<N;t+=h,i++){
		k[0] = func(t,x,y);
		k[4] = _rksys2(t,x,y);
		k[1] = func(t+0.5*h,x+0.5*h*k[0],y+0.5*h*k[4]);
		k[5] = _rksys2(t+0.5*h,x+0.5*h*k[0],y+0.5*h*k[4]);
		k[2] = func(t+0.5*h,x+0.5*h*k[1],y+0.5*h*k[5]);
		k[6] = _rksys2(t+0.5*h,x+0.5*h*k[1],y+0.5*h*k[5]);
		k[3] = func(t+h,x+h*k[2],y+h*k[6]);
		k[7] = _rksys2(t+h,x+h*k[2],y+h*k[6]);
		x += h/6*(k[0]+2*k[1]+2*k[2]+k[3]);
		y += h/6*(k[4]+2*k[5]+2*k[6]+k[7]);
	}return construct(0,y,x,NULL);}
/* solve 2-nd system of differential equations using Runge-Kutta method order 4 */
public result RK4_22(double (*f)(double,double,double),double (*g)(double,double,double),double set[4], unsigned int N){
	double h = fabs(set[1]-set[0])/(double)N;
	double y = set[2];double x = set[3];double t = set[0];
	double k[8] = {0,0,0,0,0,0,0,0};
	for(unsigned int i=0;i<N;t+=h,i++){
		k[0] = f(t,x,y);
		k[4] = g(t,x,y);
		k[1] = f(t+0.5*h,x+0.5*h*k[0],y+0.5*h*k[4]);
		k[5] = g(t+0.5*h,x+0.5*h*k[0],y+0.5*h*k[4]);
		k[2] = f(t+0.5*h,x+0.5*h*k[1],y+0.5*h*k[5]);
		k[6] = g(t+0.5*h,x+0.5*h*k[1],y+0.5*h*k[5]);
		k[3] = f(t+h,x+h*k[2],y+h*k[6]);
		k[7] = g(t+h,x+h*k[2],y+h*k[6]);
		x += h/6*(k[0]+2*k[1]+2*k[2]+k[3]);
		y += h/6*(k[4]+2*k[5]+2*k[6]+k[7]);		
	}return construct(0,x,y,NULL);}
public result AB4(double (*func)(double,double), Vector *tbl, double bot, double top, double init, unsigned int N){
	if(tbl->dim!=3){return construct(1,0,0,"Dimension match error");}
	double h = fabs(top-bot)/(double)N;double t = bot+3*h;
	double y = *(tbl->element+3);
	for(unsigned int i=3;i<N;t+=h,i++){
		y = y + h/24.0*(55*func(t,*(tbl->element+3))-59*func(t-h,*(tbl->element+2))
				  +37*func(t-2*h,*(tbl->element+1))-9*func(t-3*h,*(tbl->element)));
		for(unsigned int j=0;j<3;j++){*(tbl->element+j) = *(tbl->element+j+1);}
		*(tbl->element+3) = y;
	}return construct(0,y,h,NULL);}
public result AM4(double (*func)(double,double), Vector *tbl, double bot, double top, double init, unsigned int N){
	double h = fabs(top-bot)/(double)N;double t = bot+3*h;
	double y = *(tbl->element+3);double temp = 0;
	for(unsigned int i=3;i<N;t+=h,i++){
		temp = y + h/24.0*(55*func(t,*(tbl->element+3))-59*func(t-h,*(tbl->element+2))
				  +37*func(t-2*h,*(tbl->element+1))-9*func(t-3*h,*(tbl->element+0)));
		y = y + h/24.0*(9*func(t+h,temp)+19*func(t,*(tbl->element+3))
			-5*func(t-h,*(tbl->element+2))+func(t-2*h,*(tbl->element+1)));
		for(unsigned int j=0;j<3;j++){*(tbl->element+j) = *(tbl->element+j+1);}
		*(tbl->element+3) = y;
	}return construct(0,y,h,NULL);}
/* 2-order finite difference with matrix, Solve y" + p(x)y' + q(x)y = r(x) */
#ifdef _MSC_VER	
public result fdf2_M(Vector **M, Solution *x,double (*p)(double),double (*q)(double),double (*r)(double),double bs[4], unsigned int N){
	double h = (bs[1]-bs[0])/(double)N;double t = bs[0]+h;
	Vector *Ds = zero_vector(N);
	matrix_zeros(M,N,N);
	for(unsigned int i=0;i<N;i++,t+=h){
		if(i!=0){*(M[i]->element+i-1)=1-0.5*h*p(t);}
		if(i!=N-1){*(M[i]->element+i+1)=1+0.5*h*p(t);}
		*(M[i]->element+i)=-2+h*h*q(t);
		*(Ds->element+i)=h*h*r(t);}
	*(Ds->element) = h*h*r(bs[0]+h)-(1-0.5*h*p(bs[0]+h))*bs[2];
	*(M[N-1]->element+N-2) = 0;
	*(M[N-1]->element+N-1) = 1;
	gauss_elimination(M,Ds,x);
	erase_vector(Ds);
	return construct(0,0,h*h,NULL);}
#else
public result fdf2_M(Solution *x,double (*p)(double),double (*q)(double),double (*r)(double),double bs[4], unsigned int N){
	double h = (bs[1]-bs[0])/(double)N;double t = bs[0]+h;
	Vector *M[N-1];
	Vector *Ds = zero_vector(N);
	matrix_zeros(M,N,N);
	for(unsigned int i=0;i<N;i++,t+=h){
		if(i!=0){*(M[i]->element+i-1)=1-0.5*h*p(t);}
		if(i!=N-1){*(M[i]->element+i+1)=1+0.5*h*p(t);}
		*(M[i]->element+i)=-2+h*h*q(t);
		*(Ds->element+i)=h*h*r(t);}
	*(Ds->element) = h*h*r(bs[0]+h)-(1-0.5*h*p(bs[0]+h))*bs[2];
	*(M[N-1]->element+N-2) = 0;
	*(M[N-1]->element+N-1) = 1;	
	gauss_elimination(M,Ds,x);
	erase_vector(Ds);erase_matrix(M,N);
	return construct(0,0,h*h,NULL);}
#endif
/*-------------------------------------------------------------------------- */
/* Polynomial Operator */
public Poly *poly_set(char symbol){
	Poly *_poly;
	_poly 		  	= (Poly *)malloc(sizeof(Poly));
	_poly->e 	  	= -1;
	_poly->c 		= 0;
	_poly->overwrite= 0;
	_poly->symbol 	= symbol;
	_poly->link 	= NULL;
	return _poly;}
public Poly *poly_concat(Poly *head1, Poly *head2){
	Poly *p;
	if(!head1) return head2;
	else if(!head2) return head1;
	else {
		p = head1;
		while(p->link) p = p->link;
		p->link = head2;
		return head1;}}	
public void poly_insert(Poly **head, Poly *p, int e, double c){
	Poly *new_node;
	new_node = (Poly *)malloc(sizeof(Poly));
	new_node->e = e;
	new_node->c = c;
	new_node->symbol = p->symbol;
	if( *head == NULL ){ 	// blank list
		new_node->link = NULL;
		*head = new_node; }
	else if( p == NULL ){ 	// if p is NULL --> Insert as first
		new_node->link = *head;
		*head = new_node; }
	else { 					// insert next
		new_node->link = p->link;
		p->link = new_node;}}
public void poly_pop(Poly **head, Poly *p, Poly *removed){
	if(p==NULL && *head!=NULL){*head = (*head)->link;}
	else if(p!=NULL){
		p->link = removed->link;
		free(removed);}}
public void poly_clean(Poly *s){
	while(s->link){poly_pop(&s,s,s->link);}
	free(s);}
public void poly_dump(Poly *s){
	Poly *q = s->link;
	while(q->link){poly_pop(&q,q,q->link);}
	s->link=NULL;}
public void poly_show(Poly *head){
	Poly *p=head;int init=0;
	if(p->link){p = p->link;}
	while(p){
		if(p->e==-1){printf("%d",0);}
		else if(p->c>=0 && init){printf("+%lf*%c^(%d)",p->c,p->symbol,p->e);}
		else{printf("%lf*%c^(%d)", p->c,p->symbol,p->e);init=1;}
		p = p->link;}
	printf("\n");}
public void poly_differentiate(Poly *derv, Poly *ori){
	Poly *q;q = ori;Poly *tmp = q;
	if(ori->overwrite){while(q){
		if (q->e!=0 && q->e!=-1){
			q->c=(q->e*q->c);q->e=q->e-1;
			if(q->c==0){poly_pop(&ori,tmp,q);}}
		else if(q->e==0){poly_pop(&ori,tmp,q);}
		tmp	= q;q = q->link;}}
	else{while(q){
		if (q->e!=0 && q->e!=-1){poly_insert(&derv,derv,q->e-1,q->e*q->c);}
		q = q->link;}}}
public void poly_integrate(Poly *derv, Poly *ori, double C){
	Poly *q;q = ori;
	if(ori->overwrite){Poly *tmp=q;
		while(q){
			if (q->e!=-1){q->e=q->e+1;q->c=q->c/(q->e);}
			tmp=q;q = q->link;}
			poly_insert(&ori,tmp,0,C);}
	else{while(q){
		if (q->e!=-1){poly_insert(&derv,derv,q->e+1,q->c/(q->e+1));}
		q = q->link;}poly_insert(&derv,derv,0,C);}}
public result poly_add(Poly *new, Poly *head1, Poly *head2, int op){
	double pbit = 0,qbit = 0, ex = 0, co = 0;
	Poly *p;Poly *q;p = head1, q = head2;
	op = op<0 ? -1 : 1;
	if	(p->symbol!=q->symbol){return construct(1,0,0,"Symbol dismatch");}
	while (p!=NULL && q!=NULL){
		pbit = p->e,qbit = q->e;
		if (pbit == qbit){
			co = p->c + op*(q->c);
			ex = p->e;
			p = p->link;q = q->link;
			if(co==0){continue;}}
		else if(pbit > qbit){
			co = p->c;ex = p->e;
			p  = p->link;}
		else if(pbit < qbit){
			co = q->c;ex = q->e;
			q = q->link;}
		poly_insert(&new,new,(int)ex,co);}
		if(p==NULL){
			while(q){poly_insert(&new,new,q->e,q->c);q=q->link;}}
		else if(q==NULL){
			while(p){poly_insert(&new,new,p->e,p->c);p=p->link;}}
		if(head1->overwrite){
			poly_dump(head1);head1->link=new->link;free(new);}
		return construct(0,0,0,NULL);}
public void poly_const_mult(Poly *derv, Poly *ori, double k){
	Poly *q;q = ori;
	if(ori->overwrite){Poly *tmp=q;
		while(q){
			if (q->e!=-1){q->c=k*q->c;}
			tmp=q;q = q->link;}
			poly_insert(&ori,tmp,0,k);}
	else{while(q){
		if (q->e!=-1){poly_insert(&derv,derv,q->e,k*q->c);}
		q = q->link;}poly_insert(&derv,derv,0,k);}}
public result poly_swap(Poly *s, int e, int new_e, double c,int op){
	while(s->link){
		if(s->e!=e){s=s->link;}
		else{
			if(op==0){s->c = c;}
			else if(op>0){s->e = new_e;}
			else{s->c=c;s->e=new_e;}
			return construct(0,0,0,NULL);}}
	if(s==NULL){return construct(1,0,0,"It doesn't have such exponent");}
	return construct(0,0,0,NULL);}
public result poly_calculate(Poly *s,double v){
	Poly *p = s;double sum = 0;
	while(p){
		sum += (p->c)*pow(v,(double)p->e);
		p = p->link;}
	return construct(0,sum,0,NULL);}
public result poly_extract(Poly *s, Vector *contain[2]){
	unsigned int safe_cnt = 0;
	while(s->link!=NULL){
		*(contain[0]->element+safe_cnt) = (double)s->e;
		*(contain[1]->element+(safe_cnt++)) = s->c;
		if(safe_cnt == contain[0]->dim){break;}
		else{s = s->link;}}
	return construct(0,0,0,NULL);}
public result poly_from_vector(Poly *s, Vector *contain[2]){
	Poly *tmp = s;
	for(unsigned int i=0;i<contain[0]->dim;i++){
		poly_insert(&s,tmp,(int)*(contain[0]->element+i),*(contain[1]->element+i));
		tmp = tmp->link;}
	return construct(0,0,0,NULL);}
public void poly_expansion(Poly *s, unsigned int n){
	Poly *tmp = s;
	for(unsigned int i=0;i<=n;i++){
		poly_insert(&s,tmp,(signed)i,0);
		tmp = tmp->link;}}
public void set_overwrite(Poly *s){s->overwrite =1;}
/*-------------------------------------------------------------------------- */
/*
private typedef struct _node{
	unsigned int priority;
	void *element;
}Node;

public typedef struct _Heap{
	unsigned char type;
	unsigned int current;
	Node *elem;
}Heap;

public Heap *heap_init(char type, unsigned int capacity, unsigned int struct_size){
	Heap *_heap; 
	_heap->type 	= type;
	_heap->elem		= (Heap *)malloc((sizeof(Node));
	_heap->current 	= 0;
	_heap->capacity = capacity;
	return _heap;}

private int _perloacte(Heap *heap, int i){
	int root = (LCHILD(i) < heap->size && heap->elem[LCHILD(i)].data < heap->elem[i].data) ? LCHILD(i) : i ;
	if(RCHILD(i) < hp->size && hp->elem[RCHILD(i)].data < hp->elem[root].data) {root = RCHILD(i) ;}
	if(root != i){
		Node temp 	= &(hp->elem[i]);
		Node *n1	= temp;
		Node *n2 	= &(hp->elem[root]);
    	*n1 = *n2;*n2 = temp;	
		return root;}
	return -1;}
private void perlocate(Heap *heap, int cur){
	while(cur!=-1){cur = _perloacte(heap, cur);}}
public result heap_push(Heap *heap, void *new, unsigned int priority){
    if(heap->size){heap->elem = realloc(heap->elem,(heap->size+1)*sizeof(Node));}
    else{heap->elem = malloc(sizeof(node));}
	Node nd;
	nd.priority = priority; 
    nd.element 	= new;
	unsigned int i = (heap->size)++ ;
	while(i && nd.priority < heap->elem[PARENT(i)].priority){
        heap->elem[i] = heap->elem[PARENT(i)];
		i = PARENT(i);}
	heap->elem[i] = nd;
	return construct(0,0,0,NULL);}
public result heap_out(Heap *heap){
	if(heap->size) {
		heap->elem[0] = heap->elem[--(heap->size)] ;
        heap->elem = realloc(heap->elem, heap->size * sizeof(Node)) ;
        perlocate(heap,0);return construct(0,0,0,NULL);
	}else {
		free(heap->elem);
		return construct(1,0,0,"EMPTY_HEAP");}}
public result heap_scheduled(Heap *heap){
	for(int i = PARENT(heap->size);i>=0;i--){perloacte(heap, i);}
}
public void heap_remove(Heap *heap){free(hp->elem);}

public void levelorderTraversal(Heap *heap){
	for(int i=0;i<heap->size;i++){printf("%d ",heap->elem[i].data);}}*/
