#include "pzelgpr3d.h"
#include "pzrefprism.h"
#include "pzshapeprism.h"
#include "pzgeoel.h"

static int nsubeldata[21] = {1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,7,9,9,9,7,21};

static int subeldata[21][21][2] = {
/*00*/{{0,0}},
/*01*/{{1,1}},
/*02*/{{2,2}},
/*03*/{{4,3}},
/*04*/{{5,4}},
/*05*/{{6,5}},
/*06*///{{0,6},{0,1},{1,6}},
        {{0,1},{0,6},{1,6}},
/*07*///{{1,7},{1,2},{2,7}},
        {{1,2},{1,7},{2,7}},
/*08*///{{0,8},{0,2},{2,8}},
        {{0,2},{0,8},{2,8}},
/*09*///{{0,9},{0,3},{4,9}},
        {{0,3},{0,9},{4,9}},
/*10*///{{1,10},{1,4},{5,10}},
        {{1,4},{1,10},{5,10}},
/*11*///{{2,11},{2,5},{6,11}},
        {{2,5},{2,11},{6,11}},
/*12*///{{4,12},{4,4},{5,12}},
        {{4,4},{4,12},{5,12}},
/*13*///{{5,13},{5,5},{6,13}},
        {{5,5},{5,13},{6,13}},
/*14*///{{4,14},{4,5},{6,14}},
        {{4,5},{4,14},{6,14}},
/*15*///{{0,15},{1,15},{2,15},{3,19},{0,7},{1,8},{2,6}},
        {{0,7},{1,8},{2,6},{0,15},{1,15},{2,15},{3,19}},
/*16*///{{0,16},{1,16},{4,16},{5,16},{0,10},{0,12},{5,6},{5,9},{0,4}},
        {{0,10},{0,12},{5,6},{5,9},{0,4},{0,16},{1,16},{4,16},{5,16}},
/*17*///{{1,17},{2,17},{5,17},{6,17},{1,11},{1,13},{6,7},{6,10},{1,5}},
        {{1,11},{1,13},{6,7},{6,10},{1,5},{1,17},{2,17},{5,17},{6,17}},
/*18*///{{0,18},{2,18},{4,18},{6,18},{0,11},{0,14},{6,8},{6,9},{2,3}},
        {{0,11},{0,14},{6,8},{6,9},{2,3},{0,18},{2,18},{4,18},{6,18}},
/*19*///{{4,19},{5,19},{6,19},{7,15},{4,13},{5,14},{6,12}},
        {{4,13},{5,14},{6,12},{4,19},{5,19},{6,19},{7,15}},
/*20*///{{0,20},{1,20},{2,20},{3,20},{4,20},{5,20},{6,20},{7,20},{0,17},
      // {0,19},{4,17},{1,18},{1,19},{5,18},{2,16},{2,19},{6,16},{3,15},
      // {0,13},{1,14},{2,12}}
        {{0,13},{1,14},{2,12},{0,20},{1,20},{2,20},{3,20},{4,20},{5,20},{6,20},{7,20},{0,17},
	 {0,19},{4,17},{1,18},{1,19},{5,18},{2,16},{2,19},{6,16},{3,15}}
};

//static int MidSideNodes[12][2]  = { 
//	{0,1},{1,2},{0,2},{0,3},
//	{1,4},{2,5},{4,4},{5,5}, //original
//	{4,5},{0,4},{1,5},{2,3} 
//};
//Esta estrutura define os novos n�s introduzidos  pela divis�o
//nos  medios dos lados do pai. Os n�s locais s�o dados pela or-
//dena��o  crescente dos lados do pai.  O par {a,b}  na K-�sima 
//linha da matriz define o n� medio no K-�simo lado do pai como
//o canto local b do filho a.
//os cantos e interior do pai n�o entram nesta defini��o, 
//primeiros 6 lados mais o �ltimo (total de 21 lados).
static int MidSideNodes[15][2]  = {
	{0,1},{1,2},{0,2},
	{0,3},{1,4},{2,5},                                      //CORRIGIDO
	{4,4},{5,5},{4,5},
	{0,-10},//a face triangular inferior n�o tem n� medio
	{0,4},{1,5},{2,3},//faces quadrilaterais
	{0,-20},//a face triangular superior n�o tem n� medio
	{0,-30}//o interior n�o tem n� medio
};
//coordadas dos n�s nos medios dos lados do pai,
//dadas em ordem crescente dos lados do pai
static REAL MidCoord[15][3] = { 
	{0.5,0.,-1.},{0.5,0.5,-1.},{0.,0.5,-1.},//arestas
	{0.0,0., 0.},{1.0,0.0, 0.},{0.,1.0, 0.},//arestas
	{0.5,0., 1.},{0.5,0.5, 1.},{0.,0.5, 1.},//arestas                 //CORRIGIDO
	{-99,-99,-99},//face triangular inferior: n�o existe
	{0.5,0., 0.},{0.5,0.5, 0.},{0.,0.5, 0.},//faces quadril�terais
	{-99,-99,-99},//face triangular superior: n�o existe
    {-99,-99,-99}//interior n�o existe
};

/**
 * define as conectividades entre sub-elementos
 * linha i � filho i, {a,b,c} = {lado do filho atual,
 * irm�o vizinho,lado do vizinho}
 */
const int NumInNeigh = 19;
static int InNeigh[8][NumInNeigh][3] =  { 
	{{1,1,0},{2,2,0},{3,4,0},{4,4,1},{5,4,2},{7,3,14},{10,1,9},{11,2,9},{12,4,6},{13,4,7},{14,4,8},{17,3,18},{19,4,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
	{{0,3,5},{2,2,1},{3,7,5},{4,5,1},{5,2,4},{8,3,13},{9,3,11},{11,2,10},{12,5,6},{13,5,7},{14,3,7},{18,3,17},{19,5,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,3,3},{1,3,4},{3,7,3},{4,3,1},{5,6,2},{6,3,12},{9,3,9},{10,3,10},{12,3,6},{13,6,7},{14,6,8},{16,3,16},{19,6,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,0,5},{1,5,2},{2,0,4},{3,0,2},{4,1,2},{5,0,1},{6,6,6},{7,5,8},{8,0,13},{9,0,11},{10,1,11},{11,0,10},{12,2,6},{13,1,8},{14,0,7},{15,7,19},{16,2,16},{17,1,18},{18,0,17}},
    {{0,0,3},{1,5,0},{2,6,0},{4,5,3},{5,6,3},{6,0,12},{7,7,14},{8,0,14},{10,5,9},{11,6,9},{13,7,8},{15,0,19},{17,7,18},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,1,3},{1,1,4},{2,6,1},{3,7,2},{5,6,4},{6,1,12},{7,1,13},{8,7,13},{9,7,11},{11,6,10},{14,7,7},{15,1,19},{18,7,17},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,2,3},{1,7,4},{2,2,5},{3,7,0},{4,7,1},{6,7,12},{7,2,13},{8,2,14},{9,7,9},{10,7,10},{12,7,6},{15,2,19},{16,7,16},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,4,5},{1,5,5},{2,4,4},{3,3,0},{4,1,5},{5,3,2},{6,6,12},{7,5,14},{8,4,13},{9,4,11},{10,5,11},{11,4,10},{12,2,12},{13,1,14},{14,3,8},{16,6,16},{17,5,18},{18,4,17},{19,3,15}} 
};

/**
 * define os cantos locais dos fihos
 */
static int CornerSons[8][6] = { 
	{0,6,8,9,16,18},
	{6,1,7,16,10,17},
	{8,7,2,18,17,11},
	{18,17,16,8,7,6},          //CORRIGIDO
	{9,16,18,3,12,14},
	{16,10,17,12,4,13},
	{18,17,11,14,13,5},
	{14,13,12,18,17,16}
};

/* static int CornerSons[8][6] = {  */
/* 	{0,6,8,9,15,17}, */
/* 	{6,1,7,15,10,16}, */
/* 	{8,7,2,17,16,11}, //CORRIGIR */
/* 	{17,16,15,8,7,6}, */
/* 	{9,15,17,3,12,14}, */
/* 	{15,10,16,12,4,13}, */
/* 	{17,16,11,14,13,5}, */
/* 	{14,13,12,17,16,15} */
/* }; */

static REAL buildt[8][21][4][3] = {//por colunas
/*S0*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,-1}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*14*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*17*/{{-0.25,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.25,-0.5}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,-.5}}},
/*S1*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1.,0.,-1.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*14*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*18*/{{0.,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.25,-0.5}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,-.5}}},
/*S2*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,1.,-1.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*14*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.5,-0.5}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,.5,-.5}}},
/*S3*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*10*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*14*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,-0.5},{0.,0.,0.},{0.25,0.5,-0.5}},
      /*17*/{{0.,-0.25,0.},{0.,0.,-0.5},{0.,0.,0.},{0.5,0.25,-0.5}},
	  /*18*/{{0.25,-0.25,0},{0,0,-0.5},{0.,0.,0.},{0.25,0.25,-0.5}},      
	  /*19*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{.5,-.5,0.},{0.,0.,-.5},{0.,.5,-.5}}},
/*S4*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,1.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*08*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*14*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*17*/{{-0.25,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.25,0.5}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,.5}}},
/*S5*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1.,0.,1.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*14*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*18*/{{0.,0.25,0.},{0.,0.,0.5},{0.,0.,0.},{0.5,0.25,0.5}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,.5}}},
/*S6*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,1.,1.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*14*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.25,0.5,0.5}},
      /*17*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*19*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,.5,.5}}},
/*S7*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*08*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*10*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*11*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*12*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*13*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*14*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*15*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.25,0.,0.},{0.,0.,-0.5},{0.,0.,0.},{0.25,0.5,0.5}},
      /*17*/{{0.,-0.25,0.},{0.,0.,-0.5},{0.,0.,0.},{0.5,0.25,0.5}},
      /*18*/{{0.25,-0.25,0.},{0.,0.,-0.5},{0.,0.,0.},{0.25,0.25,0.5}},
      /*19*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*20*/{{.5,0.,0.},{.5,-.5,0.},{0.,0.,-.5},{0.,.5,.5}}}
};

static int fatherside[8][21] = {
/*00*/{0,6,8,9,16,18,6,15,8,9,16,18,16,20,18,15,16,20,18,20,20},
/*01*/{6,1,7,16,10,17,6,7,15,16,10,17,16,17,20,15,16,17,20,20,20},
/*02*/{8,7,2,18,17,11,15,7,8,18,17,11,20,17,18,15,20,17,18,20,20},
/*03*/{18,17,16,8,7,6,20,20,20,18,17,16,15,15,15,20,20,20,20,15,20},
/*04*/{9,16,18,3,12,14,16,20,18,9,16,18,12,19,14,20,16,20,18,19,20},
/*05*/{16,10,17,12,4,13,16,17,20,16,10,17,12,13,19,20,16,17,20,19,20},
/*06*/{18,17,11,14,13,5,20,17,18,18,17,11,19,13,14,20,20,17,18,19,20},
/*07*/{14,13,12,18,17,16,19,19,19,18,17,16,20,20,20,19,20,20,20,20,20},
};


//Into Divides is necesary to consider the connectivity with the all neighboards
void TPZRefPrism::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
	int i;
	SubElVec.Resize(NSubEl);
	if(geo->HasSubElement()) {
		for(i=0;i<NSubEl;i++) SubElVec[i] = geo->SubElement(i);
		return;//If exist fSubEl return this sons
	}
	int j,sub,matid=geo->MaterialId(),index;
	int np[TPZShapePrism::NSides];//guarda conectividades dos 8 subelementos

	for(j=0;j<TPZShapePrism::NNodes;j++) np[j] = geo->NodeIndex(j);
	for(j=TPZShapePrism::NNodes;j<TPZShapePrism::NSides;j++) {
		NewMidSideNode(geo,j,index);
		np[j] = index;
	}
	// creating new subelements
	TPZGeoElPr3d *pr3d = (TPZGeoElPr3d *) geo;
	TPZGeoElPr3d *sub0;//TESTE
	for(i=0;i<TPZRefPrism::NSubEl;i++) {
		TPZManVector<int> cornerindexes(TPZShapePrism::NNodes);
		for(int j=0;j<TPZShapePrism::NNodes;j++) 
			cornerindexes[j] = np[CornerSons[i][j]];
		TPZGeoElPr3d *subel = new TPZGeoElPr3d(cornerindexes,matid,*geo->Mesh());
		if(i == 0) sub0 = subel;//TESTE
		pr3d->SetSubElement(i , subel);
		SubElVec[i] = subel;
		subel->SetFather(geo);
	}

	for(sub=0;sub<NSubEl;sub++) {
		SubElVec[sub] = geo->SubElement(sub);
		SubElVec[sub]->SetFather(geo);
	}
	for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
		for(j=0;j<NumInNeigh;j++) {        //lado do subel                    numero do filho viz.    lado do viz.
			if(InNeigh[i][j][0] == -1) continue;		
			geo->SubElement(i)->SetNeighbour(InNeigh[i][j][0],TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));			
		}
	}
	geo->SetSubElementConnectivities();
}

void TPZRefPrism::NewMidSideNode(TPZGeoEl *gel,int side,int &index) {

	MidSideNodeIndex(gel,side,index);
	if(side == 15 || side > 18){
		return;//o n� geom�trico n�o pode ser criado
	}
	if(index < 0) {
		TPZGeoElSide gelside = gel->Neighbour(side);
		if(gelside.Element()) {
			while(gelside.Element() != gel) {
				gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
				if(index!=-1) return;
				gelside = gelside.Neighbour();
			}	
		}
		TPZVec<REAL> par(3,0.);
		TPZVec<REAL> coord(3,0.);
		if(side < TPZShapePrism::NNodes) {
			index = gel->NodeIndex(side); 
			return;
		}
		//aqui side = 6 a 20
		side-=TPZShapePrism::NNodes;//0,1,..,13
		par[0] = MidCoord[side][0];
		par[1] = MidCoord[side][1];
		par[2] = MidCoord[side][2];
		gel->X(par,coord);
		index = gel->Mesh()->NodeVec().AllocateNewElement();
		gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
	}
}

void TPZRefPrism::MidSideNodeIndex(TPZGeoEl *gel,int side,int &index) {
	index = -1;
	if(side == 15 || side > 18) return;
	if(side<0 || side>TPZShapePrism::NSides-1) {
		PZError << "TPZRefPrism::MidSideNodeIndex. Bad parameter side = " << side << endl;
		return;
	}
	//sides 0 a 7
	if(side<TPZShapePrism::NNodes) {//o n� medio do lado 0 � o 0 etc.
		index = (gel)->NodeIndex(side);
		return; 
	}
	//o n� medio da face � o centro e o n� medio do centro � o centro
	//como n� de algum filho se este existir
	//caso tenha filhos � o canto de algum filho, se n�o tiver filhos retorna -1
	if(gel->HasSubElement()) {
		side-=TPZShapePrism::NNodes;
		index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
	}
}

void TPZRefPrism::GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel){

	subel.Resize(0);
	if(side<0 || side>TPZShapePrism::NSides || !father->HasSubElement()){
		PZError << "TPZRefPrism::GetSubelements2 called with error arguments\n";
		return;
	}
	int nsub = NSideSubElements(side);//nsubeldata[side];
	for(int i=0;i<nsub;i++)
		subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),
												subeldata[side][i][1]));
}

int TPZRefPrism::NSideSubElements(int side) {  
	if(side<0 || side>TPZShapePrism::NSides-1){
		PZError << "TPZRefPrism::NSideSubelements2 called with error arguments\n";
		return -1;
	}
	return nsubeldata[side];
}

//int TPZRefPrism::NSideSubElements(int side) {
//  if(side < 0 || side > 26) {
//    PZError << "TPZRefPrism::NSideSubElements called for side " << side << endl;
//    return 0;
//  }
//  if(side==26) return 8;//centro
//  if(side>19 && side<26) return 4;//faces
//  if(side>7) return 2;//lados
//  return 1;//cantos
//}


TPZTransform TPZRefPrism::GetTransform(int side,int whichsubel){
  
	if(side<0 || side>TPZShapePrism::NSides-1){
		PZError << "TPZRefPrism::GetTransform side out of range or father null\n";
		return TPZTransform(0,0);
	}
	int smalldim = TPZShapePrism::SideDimension(side);
	int fatherside = FatherSide(side,whichsubel);
	int largedim = TPZShapePrism::SideDimension(fatherside);
	TPZTransform trans(largedim,smalldim);
	int i,j;
	for(i=0; i<largedim; i++) {
		for(j=0; j<smalldim; j++) {
			trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
		}
		trans.Sum() (i,0) = buildt[whichsubel][side][3][i];
	}
	return trans;
}

int TPZRefPrism::FatherSide(int side,int whichsubel){

	if(side<0 || side>TPZShapePrism::NSides-1){
		PZError << "TPZRefPrism::Father2 called error" << endl;
		return -1;
	}
	return fatherside[whichsubel][side];
}

