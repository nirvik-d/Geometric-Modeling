/*
 *  Implementation of a virtual trackball.
 *  Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *  the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 *
 *  This file now includes code from vect.C and quaternion.C, so that
 *  the trackball functionality is contained in one .C file.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "tmatrix.h"
#include "trackball.h"

/*

Manipulation of quaternions.

Some of these routines were taken from:

  Ken Shoemake, "Quaternion Calculus for Animation", published where?

*/


/******************************************************************************
Convert a quaternion into a rotation matrix.  Does *not* assume a unit
quaternion.  From Ken Shoemake.

Entry:
  q - quaternion

Exit:
  mat - rotation matrix
******************************************************************************/

void quat_to_mat(Quaternion q, Matrix mat)
{
 float s;
 float xs,ys,zs;
 float wx,wy,wz;
 float xx,xy,xz;
 float yy,yz,zz;

 /* for unit q, just set s = 2 or set xs = q[X] + q[X], etc. */

 s = 2 / (q[X]*q[X] + q[Y]*q[Y] + q[Z]*q[Z] + q[W]*q[W]);

 xs = q[X] * s;
 ys = q[Y] * s;
 zs = q[Z] * s;

 wx = q[W] * xs;
 wy = q[W] * ys;
 wz = q[W] * zs;

 xx = q[X] * xs;
 xy = q[X] * ys;
 xz = q[X] * zs;

 yy = q[Y] * ys;
 yz = q[Y] * zs;
 zz = q[Z] * zs;

 mat[X][X] = 1 - (yy + zz);
 mat[X][Y] = xy - wz;
 mat[X][Z] = xz + wy;
 mat[X][W] = 0;

 mat[Y][X] = xy + wz;
 mat[Y][Y] = 1 - (xx + zz);
 mat[Y][Z] = yz - wx;
 mat[Y][W] = 0;

 mat[Z][X] = xz - wy;
 mat[Z][Y] = yz + wx;
 mat[Z][Z] = 1 - (xx + yy);
 mat[Z][W] = 0;

 mat[W][X] = 0;
 mat[W][Y] = 0;
 mat[W][Z] = 0;
 mat[W][W] = 1;
}


/******************************************************************************
Convert a rotation matrix into a unit quaternion.  From Ken Shoemake.

Entry:
  mat - rotation matrix

Exit:
  q - quaternion
******************************************************************************/

void mat_to_quat(Matrix mat, Quaternion q)
{
  int i,j,k;
  float tr,s;
  static int nxt[3] = {Y, Z, X};  /* hey, this is a new one on me! */

  tr = mat[X][X] + mat[Y][Y] + mat[Z][Z];

  if (tr > 0) {
    s = sqrt (tr + 1);
    q[W] = s * 0.5;
    s = 0.5 / s;
    q[X] = (mat[Z][Y] - mat[Y][Z]) * s;
    q[Y] = (mat[X][Z] - mat[Z][X]) * s;
    q[Z] = (mat[Y][X] - mat[X][Y]) * s;
  }
  else {
    i = X;
    if (mat[Y][Y] > mat[X][X])
      i = Y;
    if (mat[Z][Z] > mat[i][i])
      i = Z;
    j = nxt[i];
    k = nxt[j];
    s = sqrt (1 + (mat[i][i] - (mat[j][j] + mat[k][k])));
    q[i] = s * 0.5;
    s = 0.5 / s;
    q[W] = (mat[k][j] - mat[j][k]) * s;
    q[j] = (mat[j][i] + mat[i][j]) * s;
    q[k] = (mat[k][i] + mat[i][k]) * s;
  }
}

/*
 * vect:
 *	Functions to support operations on vectors and matrices.
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 */

static float *this_vnew()
{
	register float *v;

	v = (float *) malloc(sizeof(float)*3);
	return v;
}

static void this_vcopy(const float *v1, float *v2)
{
	register int i;
	for (i = 0 ; i < 3 ; i++)
		v2[i] = v1[i];
}

static void this_vset(float *v, float x, float y, float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

static void this_vzero(float *v)
{
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
}

static float this_vlength(const float *v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void this_vscale(float *v, float div)
{
	v[0] *= div;
	v[1] *= div;
	v[2] *= div;
}

static void this_vnormal(float *v)
{
	this_vscale(v,1.0/this_vlength(v));
}

static void this_vmult(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] * src2[0];
	dst[1] = src1[1] * src2[1];
	dst[2] = src1[2] * src2[2];
}

static void this_vadd(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] + src2[0];
	dst[1] = src1[1] + src2[1];
	dst[2] = src1[2] + src2[2];
}

static void this_vsub(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] - src2[0];
	dst[1] = src1[1] - src2[1];
	dst[2] = src1[2] - src2[2];
}

static float this_vdot(const float *v1, const float *v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

static void this_vcross(const float *v1, const float *v2, float *cross)
{
	float temp[3];

	temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
	this_vcopy(temp, cross);
}

static void vhalf(const float *v1, const float *v2, float *half)
{
	float len;

	this_vadd(v2,v1,half);
	len = this_vlength(half);
	if(len>0.0001)
		this_vscale(half,1.0/len);
	else
		this_vcopy(v1, half);
}

static void vdirection(const float *v1, float *dir)
{
	this_vcopy(v1, dir);
	this_vnormal(dir);
}

static void vreflect(const float *in, const float *mirror, float *out)
{
	float temp[3];

	this_vcopy(mirror, temp);
	this_vscale(temp,this_vdot(mirror,in));
	this_vsub(temp,in,out);
	this_vadd(temp,out,out);
}

static void vmultmatrix(const Matrix m1, const Matrix m2, Matrix prod)
{
	register int row, col;
	Matrix temp;

	for(row=0 ; row<4 ; row++) 
		for(col=0 ; col<4 ; col++)
			temp[row][col] = m1[row][0] * m2[0][col]
						   + m1[row][1] * m2[1][col]
						   + m1[row][2] * m2[2][col]
						   + m1[row][3] * m2[3][col];
	for(row=0 ; row<4 ; row++) 
		for(col=0 ; col<4 ; col++)
		prod[row][col] = temp[row][col];
}

static void vtransform(const float *v, const Matrix mat, float *vt)
{
	float t[3];

	t[0] = v[0]*mat[0][0] + v[1]*mat[1][0] + v[2]*mat[2][0] + mat[3][0];
	t[1] = v[0]*mat[0][1] + v[1]*mat[1][1] + v[2]*mat[2][1] + mat[3][1];
	t[2] = v[0]*mat[0][2] + v[1]*mat[1][2] + v[2]*mat[2][2] + mat[3][2];
	this_vcopy(t, vt);
}
static void vtransform4(const float *v, const Matrix mat, float *vt)
{
	float t[4];

	t[0] = v[0]*mat[0][0] + v[1]*mat[1][0] + v[2]*mat[2][0] + mat[3][0];
	t[1] = v[0]*mat[0][1] + v[1]*mat[1][1] + v[2]*mat[2][1] + mat[3][1];
	t[2] = v[0]*mat[0][2] + v[1]*mat[1][2] + v[2]*mat[2][2] + mat[3][2];
	this_vcopy(t, vt);
	t[3] = v[0]*mat[0][3] + v[1]*mat[1][3] + v[2]*mat[2][3] + mat[3][3];
	vt[3] = t[3];
}

Matrix idmatrix =
{
	{ 1.0, 0.0, 0.0, 0.0,},
	{ 0.0, 1.0, 0.0, 0.0,},
	{ 0.0, 0.0, 1.0, 0.0,},
	{ 0.0, 0.0, 0.0, 1.0,},
};

static void mcopy(const Matrix m1, Matrix m2)
{
	int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			m2[i][j] = m1[i][j];
}

static void minvert(const Matrix mat, Matrix result)
{
	int i, j, k;
	double temp;
	double m[8][4];
	/*   Declare identity matrix   */

	mcopy(idmatrix, result);
	for (i = 0;  i < 4;  i++) {
		for (j = 0;  j < 4;  j++) {
			m[i][j] = mat[i][j];
			m[i+4][j] = result[i][j];
		}
	}

	/*   Work across by columns   */

	for (i = 0;  i < 4;  i++) {
		for (j = i;  (m[i][j] == 0.0) && (j < 4);  j++)
				;
		if (j == 4) {
			fprintf (stderr, "error:  cannot do inverse matrix\n");
			exit (2);
		} 
		else if (i != j) {
			for (k = 0;  k < 8;  k++) {
				temp = m[k][i];   
				m[k][i] = m[k][j];   
				m[k][j] = temp;
			}
		}

		/*   Divide original row   */

		for (j = 7;  j >= i;  j--)
			m[j][i] /= m[i][i];

		/*   Subtract other rows   */

		for (j = 0;  j < 4;  j++)
			if (i != j)
				for (k = 7;  k >= i;  k--)
					m[k][j] -= m[k][i] * m[i][j];
	}

	for (i = 0;  i < 4;  i++)
		for (j = 0;  j < 4;  j++)
				result[i][j] = m[i+4][j];
}


/* Trackball stuff below */


/*
 * This size should really be based on the distance from the center of
 * rotation to the point on the object underneath the mouse.  That
 * point would then track the mouse as closely as possible.  This is a
 * simple example, though, so that is left as an Exercise for the
 * Programmer.
 */
#define TRACKBALLSIZE  (0.8)

/*
 * Local function prototypes (not defined in trackball.h)
 */
float tb_project_to_sphere(float, float, float);
void normalize_quat(float [4]);

/*
 * Ok, simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center.  This particular function was chosen after trying out
 * several variations.
 * 
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0)
 */
void
trackball(float q[4], float p1x, float p1y, float p2x, float p2y)
{
	float a[3];	/* Axis of rotation */
	float phi;	/* how much to rotate about axis */
	float p1[3], p2[3], d[3];
	float t;

	if (p1x == p2x && p1y == p2y)
	{
		this_vzero(q); q[3] = 1.0; /* Zero rotation */
		return;
	}

/*
 * First, figure out z-coordinates for projection of P1 and P2 to
 * deformed sphere
 */
	this_vset(p1,p1x,p1y,tb_project_to_sphere(TRACKBALLSIZE,p1x,p1y));
	this_vset(p2,p2x,p2y,tb_project_to_sphere(TRACKBALLSIZE,p2x,p2y));

/*
 *	Now, we want the cross product of P1 and P2
 */
	this_vcross(p2,p1,a);

/*
 *	Figure out how much to rotate around that axis.
 */
	this_vsub(p1,p2,d);
	t = this_vlength(d) / (2.0*TRACKBALLSIZE);
	/*
	 * Avoid problems with out-of-control values...
	 */
	if (t > 1.0) t = 1.0;
	if (t < -1.0) t = -1.0;
	phi = 2.0 * asin(t);

	axis_to_quat(a,phi,q);
}

/*
 *	Given an axis and angle, compute quaternion.
 */
void
axis_to_quat(float a[3], float phi, float q[4])
{
	this_vnormal(a);
	this_vcopy(a,q);
	this_vscale(q,sin(phi/2.0));
	q[3] = cos(phi/2.0);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
static float
tb_project_to_sphere(float r, float x, float y)
{
	float d, t, z;

	d = sqrt(x*x + y*y);
	if (d < r*sqrt(0.5))  	/* Inside sphere */
	z = sqrt(r*r - d*d);
	else
	{ 			/* On hyperbola */
		t = r / sqrt(2.0);
		z = t*t / d;
	}
	return z;
}

/*
 * Given two rotations, e1 and e2, expressed as quaternion rotations,
 * figure out the equivalent single rotation and stuff it into dest.
 * 
 * This routine also normalizes the result every RENORMCOUNT times it is
 * called, to keep error from creeping in.
 *
 * NOTE: This routine is written so that q1 or q2 may be the same
 * as dest (or each other).
 */

#define RENORMCOUNT 97

void
add_quats(float q1[4], float q2[4], float dest[4])
{
	static int count=0;
	int i;
	float t1[4], t2[4], t3[4];
	float tf[4];

	this_vcopy(q1,t1); 
	this_vscale(t1,q2[3]);

	this_vcopy(q2,t2); 
	this_vscale(t2,q1[3]);

	this_vcross(q2,q1,t3);
	this_vadd(t1,t2,tf);
	this_vadd(t3,tf,tf);
	tf[3] = q1[3] * q2[3] - this_vdot(q1,q2);

	dest[0] = tf[0];
	dest[1] = tf[1];
	dest[2] = tf[2];
	dest[3] = tf[3];

	if (++count > RENORMCOUNT)
	{
		count = 0;
		normalize_quat(dest);
	}
}

/*
 * Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * If they don't add up to 1.0, dividing by their magnitued will
 * renormalize them.
 *
 * Note: See the following for more information on quaternions:
 * 
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
static void
normalize_quat(float q[4])
{
	int i;
	float mag;

	mag = (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	for (i = 0; i < 4; i++) q[i] /= mag;
}

/*
 * Build a rotation matrix, given a quaternion rotation.
 *
 * Assumes unit quaternion.
 *
 */
void
build_rotmatrix(float m[4][4], float q[4])
{
	m[0][0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
	m[0][1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
	m[0][2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);
	m[0][3] = 0.0;

	m[1][0] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
	m[1][1] = 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
	m[1][2] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
	m[1][3] = 0.0;

	m[2][0] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
	m[2][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
	m[2][2] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);
	m[2][3] = 0.0;

	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
}

