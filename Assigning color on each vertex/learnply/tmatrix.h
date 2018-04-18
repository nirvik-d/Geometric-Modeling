/*

Matrix and Vector header.

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#ifndef __GL_GL_H__
typedef float Matrix[4][4];
#endif

typedef float Vector[3];
typedef float Vector4[4];
typedef float Plane[4];
typedef float Quaternion[4];

#define X 0
#define Y 1
#define Z 2
#define W 3

float vnorm(Vector);
float vdot(Vector,Vector);
void vadd(Vector,Vector,Vector);
void vsub(Vector,Vector,Vector);
float vlen(Vector);
void vset(float *, float, float, float);
void vec_scale(float *, float);
void vapply(Matrix, Vector, Vector);
void vcopy(Vector, Vector);
void push();
void pop();
void mat_ident(Matrix);
void mat_mult (Matrix prod, Matrix a, Matrix b);
void translate (float x, float y, float z);
void mat_copy (Matrix dest, Matrix source);
void mat_transpose (Matrix m);
void mat_translate (Matrix m, float x, float y, float z);
void mat_rotate (Matrix mat, float angle, char axis);
void mat_apply (Matrix m, Vector v);
void mat_print (Matrix m);

void quat_to_mat(Quaternion q, Matrix mat);
void mat_to_quat(Matrix mat, Quaternion q);

float mat_invert (Matrix mat);

#endif /* __MATRIX_H__ */

