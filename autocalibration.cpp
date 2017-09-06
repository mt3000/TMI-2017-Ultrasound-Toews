
/*

Copyright (c) 2017, Matthew Toews
All rights reserved.

By downloading, copying, installing or using the software you agree to this license.
If you do not agree to this license, do not download, install, copy or use the software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
 3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

*/




#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <algorithm>
#include <vector>
#include <map>
using namespace std;

#include "neldermead.h"
//#include "nlopt.h"

#define PARAMETER_DIM_ROT 3
#define PARAMETER_DIM_RIG 7

void
vec3D_cross_3d(
			 float *pv1,
			 float *pv2,
			 float *pvCross
			 )
{
	pvCross[0] =  pv1[1]*pv2[2] - pv1[2]*pv2[1];
	pvCross[1] = -pv1[0]*pv2[2] + pv1[2]*pv2[0];
	pvCross[2] =  pv1[0]*pv2[1] - pv1[1]*pv2[0];
}

void
vec3D_diff_3d(
			 float *pv1,
			 float *pv2,
			 float *pv12
			 )
{
	pv12[0] = pv2[0]-pv1[0];
	pv12[1] = pv2[1]-pv1[1];
	pv12[2] = pv2[2]-pv1[2];
}

float
vec3D_dot_3d(
			 float *pv1,
			 float *pv2
			 )
{
	return pv2[0]*pv1[0] + pv2[1]*pv1[1] + pv2[2]*pv1[2];
}


float
vec3D_dist_3d(
			 float *pv1,
			 float *pv2
			 )
{

	float dx = pv2[0]-pv1[0];
	float dy = pv2[1]-pv1[1];
	float dz = pv2[2]-pv1[2];
	return sqrt( dx*dx+dy*dy+dz*dz );
}

float
vec3D_distsqr_3d(
			 float *pv1,
			 float *pv2
			 )
{

	float dx = pv2[0]-pv1[0];
	float dy = pv2[1]-pv1[1];
	float dz = pv2[2]-pv1[2];
	return dx*dx+dy*dy+dz*dz;
}

void
vec3D_mult_scalar(
				  float *pf1,
				  float mult
				  )
{
	pf1[0] *= mult;
	pf1[1] *= mult;
	pf1[2] *= mult;
}

void
vec3D_norm_3d(
		 float *pf1
		 )
{
	float fSumSqr = pf1[0]*pf1[0] + pf1[1]*pf1[1] + pf1[2]*pf1[2];
	if( fSumSqr > 0 )
	{
		float fDiv = 1.0/sqrt( fSumSqr );
		pf1[0] *= fDiv;
		pf1[1] *= fDiv;
		pf1[2] *= fDiv;
	}
	else
	{
		pf1[0] = 1;
		pf1[1] = 0;
		pf1[2] = 0;
	}
}

float
vec3D_mag(
		 float *pf1
		 )
{
	float fSumSqr = pf1[0]*pf1[0] + pf1[1]*pf1[1] + pf1[2]*pf1[2];
	if( fSumSqr > 0 )
	{
		return sqrt( fSumSqr );
	}
	else
	{
		return 0;
	}
}

void
mult_3x3_matrix(
			float *mat1,
			float *mat2,
			float *mat_out
			)
{
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			mat_out[i*3+j]=0;
			for( int ii = 0; ii < 3; ii++ )
			{
				mat_out[i*3+j] += mat1[i*3+ii]*mat2[ii*3+j];
			}
		}
	}
}

void
mult_3x3_scalar(
			float *mat1,
			float multiple,
			float *mat_out
			)
{
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			mat_out[i*3+j] = mat1[i*3+j]*multiple;
		}
	}
}

void
mult_1x3_scalar(
			float *vec,
			float multiple,
			float *vec_out
			)
{
	for( int i = 0; i < 3; i++ )
	{
		vec_out[i] = vec[i]*multiple;
	}
}

void
sum_3x3_matrix(
			float *mat1,
			float *mat2,
			float *mat_out
			)
{
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			mat_out[i*3+j] = mat1[i*3+j] + mat2[i*3+j];
		}
	}
}

void
mult_3x3_vector(
			float *mat,
			float *vec_in,
			float *vec_out
			)
{
	for( int i = 0; i < 3; i++ )
	{
		vec_out[i] = 0;
		for( int j = 0; j < 3; j++ )
		{
			vec_out[i] += mat[i*3+j]*vec_in[j];
		}
	}
}

void
mult_4x4_matrix(
			float *mat1,
			float *mat2,
			float *mat_out
			)
{
	for( int i = 0; i < 4; i++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			mat_out[i*4+j]=0;
			for( int ii = 0; ii < 4; ii++ )
			{
				mat_out[i*4+j] += mat1[i*4+ii]*mat2[ii*4+j];
			}
		}
	}
}

void
mult_4x4_vector(
			float *mat,
			float *vec_in,
			float *vec_out
			)
{
	for( int i = 0; i < 4; i++ )
	{
		vec_out[i] = 0;
		for( int j = 0; j < 4; j++ )
		{
			vec_out[i] += mat[i*4+j]*vec_in[j];
		}
	}
}

//
// Assume last row of mat is [0,0,0,1] and that
// vec_in[3] = 1
// vec_out[3] = 1
//
// Here, mat is a 4x4 matrix, vec_in is 1x3, vec_out is 1x3
//
void
mult_4x4_vector_homogenous(
			float *mat,
			float *vec_in,
			float *vec_out
			)
{
	for( int i = 0; i < 3; i++ )
	{
		vec_out[i] = 0;
		for( int j = 0; j < 3; j++ )
		{
			vec_out[i] += mat[i*4+j]*vec_in[j];
		}
		vec_out[i] += mat[i*4+3];
	}
	//vec_out[3] = 1;
}

void
mult_4x4_scalar(
			float *mat1,
			float multiple,
			float *mat_out
			)
{
	for( int i = 0; i < 4; i++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			mat_out[i*4+j] = mat1[i*4+j]*multiple;
		}
	}
}

int
invert_4x4(
		   float *m,
		   float *invOut
		   )
{
    float inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}



nlopt_result 
nldrmd_minimize(
							 int n,
							 nlopt_func f,
							 void *f_data,
			    
							 const double *lb,
							 const double *ub, /* bounds */
			     
							 double *x, /* in: initial guess, out: minimizer */
			     
							 double *minf,
			     
							 const double *xstep, /* initial step sizes */
			     
							 nlopt_stopping *stop
							 );

//
// Declare constraint data
//
typedef struct _my_constraint_data
{
	int iPoints;
	float *pfPoints1;
	float *pfPoints2;

	float fPercentToKeep; // usuall 0.75 of lowest error points, for robustness

	int iErrorCalib; // Boolean, if true, algorithm evaluates error of pfPoints1/2 using ground truth calib
	float pfCalib[PARAMETER_DIM_RIG]; // Ground truth calibration matrix: 16 floats

	int bPerformAlign; // If true, generate and apply the alignment matrix
	float pfAlign[PARAMETER_DIM_RIG]; // Alignment matrix: 16 floats

	// *** 2014 3D reconstruction
	// Keep track of frame indices 
	int *piFrames1;
	int *piFrames2;

} my_constraint_data;

int
convert_vector_to_rotation_matrix(
								  float *vec,
								  float *mat
								  )
{
	float pfOmega[3];

	memcpy( pfOmega, vec, sizeof(pfOmega) );

	float fTheta = vec3D_mag( pfOmega );

	vec3D_norm_3d( pfOmega );

	float fThetaSin = sin( fTheta );
	float fThetaCos = cos( fTheta );
	float fThetaCosInv = 1.0f-fThetaCos;

	float pfOmegaHat[9];
	memset( pfOmegaHat, 0, sizeof(pfOmegaHat) );
	pfOmegaHat[1] = -pfOmega[2];
	pfOmegaHat[2] =  pfOmega[1];
	pfOmegaHat[3] =  pfOmega[2];
	pfOmegaHat[5] = -pfOmega[0];
	pfOmegaHat[6] = -pfOmega[1];
	pfOmegaHat[7] =  pfOmega[0];

	float pfOmegaHatSqr[9];
	mult_3x3_matrix( pfOmegaHat, pfOmegaHat, pfOmegaHatSqr );

	memset( mat, 0, sizeof(float)*9 );
	mat[0]=1;
	mat[4]=1;
	mat[8]=1;

	mult_3x3_scalar( pfOmegaHat, fThetaSin, pfOmegaHat );
	mult_3x3_scalar( pfOmegaHatSqr, fThetaCosInv, pfOmegaHatSqr );
	
	sum_3x3_matrix( mat, pfOmegaHat, mat );
	sum_3x3_matrix( mat, pfOmegaHatSqr, mat );

	return 1;
}

int
convert_rotation_matrix_to_vector(
								  float *vec,
								  float *mat
								  )
{

	float fTrace = mat[0] + mat[4] + mat[8];
	float fTheta = acos( (fTrace - 1.0) / 2.0 );

	vec[0] = mat[7]-mat[5];
	vec[1] = mat[2]-mat[6];
	vec[2] = mat[3]-mat[1];
	vec3D_norm_3d( vec );

	vec[0] *= fTheta;
	vec[1] *= fTheta;
	vec[2] *= fTheta;

	return 1;
}

int
set_vector_to_identity(
					   float *vec
					   )
{
	for( int m = 0; m < 7; m++ )
	{
		// zero rotation, 
		vec[m] = 0;
	}
	return 1;
}

int
set_vector_to_identity(
					   double *vec
					   )
{
	for( int m = 0; m < 7; m++ )
	{
		// zero rotation, 
		vec[m] = 0;
	}
	return 1;
}


int
convert_vector_to_matrix4x4(
								  float *vec,
								  float *mat44 // length 16 vector
								  )
{
	convert_vector_to_rotation_matrix( vec, mat44 );

	//mult_3x3_scalar( mat44, exp(vec[6]), mat44 );

	mat44[10] = mat44[8];
	mat44[9] = mat44[7];
	mat44[8] = mat44[6];

	mat44[6] = mat44[5];
	mat44[5] = mat44[4];
	mat44[4] = mat44[3];

	mat44[3]  = vec[3];
	mat44[7]  = vec[4];
	mat44[11] = vec[5];

	mat44[12] = 0;
	mat44[13] = 0;
	mat44[14] = 0;
	mat44[15] = 1;



	return 1;
}

//
// This code has been verified.
// The plane parameters (a,b,c) are normalized to a unit vector.
//
int
determine_plane_3point(
					   float *pfP01, float *pfP02, float *pfP03, // points in image 1 (3d)
					   float *pfPlane // ax + by + cz + d = 0;
					   )
{
	float pfV0_12[3];
	float pfV0_13[3];
	float pfV0_nm[3];

	// Subtract point 1 to convert to vectors
	vec3D_diff_3d( pfP01, pfP02, pfV0_12 );
	vec3D_diff_3d( pfP01, pfP03, pfV0_13 );

	// Normalize vectors
	//vec3D_norm_3d( pfV0_12 );
	//vec3D_norm_3d( pfV0_13 );

	// Cross product between 2 vectors to get normal 
	vec3D_cross_3d( pfV0_12, pfV0_13, pfV0_nm );
	vec3D_norm_3d(  pfV0_nm );

	pfPlane[0] = pfV0_nm[0]; // a
	pfPlane[1] = pfV0_nm[1]; // b
	pfPlane[2] = pfV0_nm[2]; // c
	pfPlane[3] = -(pfV0_nm[0]*pfP01[0] + pfV0_nm[1]*pfP01[1] + pfV0_nm[2]*pfP01[2]); // d
	pfPlane[3] = -(pfV0_nm[0]*pfP02[0] + pfV0_nm[1]*pfP02[1] + pfV0_nm[2]*pfP02[2]); // d
	pfPlane[3] = -(pfV0_nm[0]*pfP03[0] + pfV0_nm[1]*pfP03[1] + pfV0_nm[2]*pfP03[2]); // d

	return 1;
}

//
// ransac_plane_estimation()
//
// Determine robust plane estimate from 3 or more points.
//
int
ransac_plane_estimation(
						int iPoints,
						float *pts,
						int iIterations,
						float fDistThreshold,
						float *pfPlane // 4 parameters: a b c d, abc is normalized
						)
{
	if( iPoints < 3 )
	{
		// 3 unique points are required
		return -1;
	}

	if( iPoints == 3 )
	{
		// Determine exact plane
		determine_plane_3point( pts + 3*0, pts + 3*1, pts + 3*2, pfPlane );
		return 3;
	}

	int iMaxInliers = 0;
	for( int i = 0; i < iIterations; i++ )
	{
		// Find a set of unique points
		int i1 = (rand()*iPoints)/(RAND_MAX+1.0f);
		int i2 = (rand()*iPoints)/(RAND_MAX+1.0f);
		int i3 = (rand()*iPoints)/(RAND_MAX+1.0f);
		while( i1 == i2 || i1 == i3 || i2 == i3 )
		{
			i1 = (rand()*iPoints)/(RAND_MAX+1.0f);
			i2 = (rand()*iPoints)/(RAND_MAX+1.0f);
			i3 = (rand()*iPoints)/(RAND_MAX+1.0f);
		}

		float *pfP01 = pts + 3*i1;
		float *pfP02 = pts + 3*i2;
		float *pfP03 = pts + 3*i3;
		float pfPTest[4];

		determine_plane_3point( pfP01, pfP02, pfP03, pfPTest );
		if( pfPTest[0]*pfPTest[1]*pfPTest[2]*pfPTest[3] == 0 )
		{
			// Error - probably duplicate points
			continue;
		}

		int iInliers = 0;
		for( int j = 0; j < iPoints; j++ )
		{
			// Compute distance from each point to plane
			float *pfPt = pts + 3*j;
			float fDist = pfPTest[0]*pfPt[0] + pfPTest[1]*pfPt[1] + pfPTest[2]*pfPt[2] + pfPTest[3];
			if( fDist < 0 )
			{
				fDist = -fDist;
			}
			if( fDist < fDistThreshold )
			{
				iInliers++;
			}
		}

		if( iInliers > iMaxInliers )
		{
			iMaxInliers = iInliers;
			for( int j = 0; j < 4; j++ )
			{
				pfPlane[j] = pfPTest[j];
			}

			// Go through and compute error again
			for( int j = 0; j < iPoints; j++ )
			{
				// Compute distance from each point to plane
				float *pfPt = pts + 3*j;
				float fDist = pfPTest[0]*pfPt[0] + pfPTest[1]*pfPt[1] + pfPTest[2]*pfPt[2] + pfPTest[3];
				if( fDist < 0 )
				{
					fDist = -fDist;
				}
				if( fDist < fDistThreshold )
				{
					iInliers++;
				}
			}
		}

	}

	return iMaxInliers;
}



int
convert_vector_to_rotation_matrix_4X4(
								  float *vec,
								  float *mat
								  )
{
	convert_vector_to_rotation_matrix( vec, mat );

	mat[10] = mat[8];
	mat[9] = mat[7];
	mat[8] = mat[6];
	mat[6] = mat[5];
	mat[5] = mat[4];
	mat[4] = mat[3];
	mat[3] = 0;
	mat[7] = 0;
	mat[11] = 0;
	mat[12] = mat[13] = mat[14] = 0;
	mat[15] = 1;

	return 1;
}

int
generate_random_rigid_transformation_matrix_4x4(
								  float *mat
	)
{
	float vec[3];
	vec[0] = rand()/(RAND_MAX+1.0f) - 0.5;
	vec[1] = rand()/(RAND_MAX+1.0f) - 0.5;
	vec[2] = rand()/(RAND_MAX+1.0f) - 0.5;

	convert_vector_to_rotation_matrix_4X4( vec, mat );

	// Put in translation
	mat[3] = 100*(rand()/(RAND_MAX+1.0f) - 0.5);
	mat[7] = 100*(rand()/(RAND_MAX+1.0f) - 0.5);
	mat[11] = 100*(rand()/(RAND_MAX+1.0f) - 0.5);

	return 1;
}

int
test_convert_vector_to_rotation_matrix(
								  )
{
	//float vec[3] = {-0.406417548, -0.247176698, -0.334497981 };
	//float vec[3] = {-2.19563893, -1.335352724, -1.80709911};
	float vec[3] = {-6.586916789, -4.006058171, -5.421297331};

	float vec2[3];
	float mat[9];

	convert_vector_to_rotation_matrix( vec, mat );

	float vec_out[3];
	mult_3x3_vector( mat, vec, vec_out );

	mult_1x3_scalar( vec, 102, vec2 );
	convert_vector_to_rotation_matrix( vec2, mat );

	mult_3x3_vector( mat, vec, vec_out );
	return 1;
}

//
// select_half_data()
//
// Select a random subset of half the data.
//
int
select_half_data(
				 my_constraint_data &dat,
								int iSubset
								)
{
	int iNewPoints =  dat.iPoints / 2;
	if( iSubset > 0 )
	{
		dat.pfPoints1 += 15*iNewPoints;
		dat.pfPoints2 += 15*iNewPoints;
	}
	dat.iPoints = iNewPoints;

	return 1;
}

int
select_half_data_interleaved(
				 my_constraint_data &dat,
								int iSubset
								)
{
	int iBit = iSubset != 0 ? 1 : 0;

	int iNewPoints =  dat.iPoints / 2;
	for( int i = 0; i < iNewPoints; i++ )
	{
		memcpy( dat.pfPoints1 + 15*i, dat.pfPoints1 + 15*(2*i+iBit), 15*sizeof(float) );
		memcpy( dat.pfPoints2 + 15*i, dat.pfPoints2 + 15*(2*i+iBit), 15*sizeof(float) );
	}
	dat.iPoints = iNewPoints;

	return 1;
}


int
generate_minimization_test_data(
								my_constraint_data &dat,
								int iPoints
								)
{
	dat.iPoints = iPoints;

	dat.pfPoints1 = new float[3*iPoints];
	dat.pfPoints2 = new float[3*iPoints];

	float vec_original[3] = {-6.586916789, -4.006058171, -5.421297331};

	for( int i = 0; i < iPoints; i++ )
	{
		// Create random vector & rotation matrix
		float vec2[3];
		memcpy( vec2, vec_original, sizeof(vec2) );
		vec2[0] += rand()/(RAND_MAX+1.0) - 0.5f;
		vec2[1] += rand()/(RAND_MAX+1.0) - 0.5f;
		vec2[2] += rand()/(RAND_MAX+1.0) - 0.5f;
		float mat[9];
		convert_vector_to_rotation_matrix( vec2, mat );

		float *point_in = dat.pfPoints1 + 3*i;
		float *point_out = dat.pfPoints2 + 3*i;

		// Create random input point (-10 to 10)
		point_in[0] = (rand()/(RAND_MAX+1.0) - 0.5f)*10;
		point_in[1] = (rand()/(RAND_MAX+1.0) - 0.5f)*10;
		point_in[2] = (rand()/(RAND_MAX+1.0) - 0.5f)*10;

		// Create random output point: multiply input point by random rotation matrix
		mult_3x3_vector( mat, point_in, point_out );

	}

	return 1;
}

//
// Format: | rot-vec 3 | trs-vec 3 | log scale 1 |
//
int
generate_minimization_test_data_rigid(
								my_constraint_data &dat,
								int iPoints
								)
{
	dat.iPoints = iPoints;

	dat.pfPoints1 = new float[3*iPoints];
	dat.pfPoints2 = new float[3*iPoints];

	//srand(50);

#define LN_2 0.69314718

	float vec_original[PARAMETER_DIM_RIG] = {-6.586916789, -4.006058171, -5.421297331, 30, -25, -5, LN_2 };

	for( int i = 0; i < iPoints; i++ )
	{
		// Create random vector

		float vec2[PARAMETER_DIM_RIG];
		memcpy( vec2, vec_original, sizeof(vec2) );

		// Theta +- 0.5
		vec2[0] += rand()/(RAND_MAX+1.0) - 0.5f;
		vec2[1] += rand()/(RAND_MAX+1.0) - 0.5f;
		vec2[2] += rand()/(RAND_MAX+1.0) - 0.5f;

		// Trans +- 10
		vec2[3] += (rand()/(RAND_MAX+1.0) - 0.5f)*10;
		vec2[4] += (rand()/(RAND_MAX+1.0) - 0.5f)*10;
		vec2[5] += (rand()/(RAND_MAX+1.0) - 0.5f)*10;

		// Scale +- 0.25
		vec2[6] += (rand()/(RAND_MAX+1.0) - 0.5f)*.5;

		float mat[9];

		// Create rotation matrix
		convert_vector_to_rotation_matrix( vec2, mat );
		// Multiply scale
		mult_3x3_scalar( mat, exp(vec2[6]), mat );

		float *point_in = dat.pfPoints1 + 3*i;
		float *point_out = dat.pfPoints2 + 3*i;

		// Create random input point (-50 to 50)
		point_in[0] = (rand()/(RAND_MAX+1.0) - 0.5f)*100;
		point_in[1] = (rand()/(RAND_MAX+1.0) - 0.5f)*100;
		point_in[2] = (rand()/(RAND_MAX+1.0) - 0.5f)*100;

		// Create random output point: multiply input point by random rotation matrix
		mult_3x3_vector( mat, point_in, point_out );

		// Add translation and we're done
		point_out[0] += vec2[3];
		point_out[1] += vec2[4];
		point_out[2] += vec2[5];
	}

	return 1;
}

//
// Format: | rot-vec 3 | trs-vec 3 | log scale 1 |
//
int
generate_minimization_test_data_rigid_US(
								my_constraint_data &dat,
								int iPoints
								)
{
	dat.iPoints = iPoints;

	dat.pfPoints1 = new float[15*iPoints];
	dat.pfPoints2 = new float[15*iPoints];

#define LN_2 0.69314718

	float vec_original[PARAMETER_DIM_RIG] = {0.0001, 0.0001, 0.0001, 30, -25, -5, LN_2 };

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];

	// Create rotation matrix
	convert_vector_to_rotation_matrix( vec_original, mat );
	// Multiply scale
	mult_3x3_scalar( mat, exp(vec_original[6]), mat );

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	memset( &mat44[0], 0, sizeof(mat44) );
	memcpy( &mat44[0], &mat[0], 3*sizeof(float) );
	memcpy( &mat44[4], &mat[3], 3*sizeof(float) );
	memcpy( &mat44[8], &mat[6], 3*sizeof(float) );

	// Copy translations to 4x4 homegenous matrix
	mat44[3] = vec_original[3];
	mat44[7] = vec_original[4];
	mat44[11] = vec_original[5];
	mat44[15] = 1;

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	for( int i = 0; i < iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		float *p_in = point_in+12;
		float *p_out = point_out+12;

		// Generate left hand side
		p_in[0] = (rand()*100)/(RAND_MAX+1.0);
		p_in[1] = (rand()*100)/(RAND_MAX+1.0);
		p_in[2] = 0;

		generate_random_rigid_transformation_matrix_4x4( mat44_1 );

		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );

		float point_test_in[3];
		mult_4x4_vector_homogenous( mat44_1_p,p_in, point_test_in );


		// Generate right hand side
		p_out[0] = (rand()*100)/(RAND_MAX+1.0);
		p_out[1] = (rand()*100)/(RAND_MAX+1.0);
		p_out[2] = 0;

		float point_test_out[3];
		mult_4x4_vector_homogenous( mat44,p_out, point_test_out );

		mat44_2[0] = mat44_2[5] = mat44_2[10] = 1; mat44_2[15] = 1;

		// Add Translation
		mat44_2[ 3] = point_test_in[0] - point_test_out[0];
		mat44_2[ 7] = point_test_in[1] - point_test_out[1];
		mat44_2[11] = point_test_in[2] - point_test_out[2];

		// Save matrices (points already saved)
		memcpy( point_in, &mat44_1[0], sizeof(float)*12 );
		memcpy( point_out, &mat44_2[0], sizeof(float)*12 );


		// Now perform self test
		//
		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Multiply current calibration matrix
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_matrix( mat44_2, mat44, mat44_2_p );

		// Generate output 3D test points

		mult_4x4_vector_homogenous( mat44_1_p, point_in+12, point_test_in );
		mult_4x4_vector_homogenous( mat44_2_p, point_out+12, point_test_out );
	

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		assert( fErrorSqr*fErrorSqr < 0.1*0.1 );

	}

	return 1;
}

int
allocateArray(
		  float ***pppdX,
		  int iCols,
		  int iRows
		  )
{
	float *pdDataX = new float[iRows*iCols];
	if( !pdDataX )
	{
		return -1;
	}

	float **pdX = new float*[iRows];
	for( int i = 0; i < iRows; i++ )
	{
		pdX[i] = &pdDataX[i*iCols];
	}

	*pppdX = pdX;
	return 1;
}

int
deallocateArray(
		  float ***pppdX
				)
{
	float **pdX = *pppdX;
	delete [] pdX[0];
	delete *pppdX;
	return 1;
}


int
readArray(
		  char *pcFileNameBase,
		  float ***pppdX,
		  int &iCols,
		  int &iRows,
		  char *pcFormat
		  )
{
	FILE *infile;
	char pcFileName[300];
	int iReturn;
	
	*pppdX = 0;

	// Open header
	sprintf( pcFileName, "%s.txt", pcFileNameBase );
	infile = fopen( pcFileName, "rt" );
	if( !infile )
	{
		return -1;
	}
	iReturn = fscanf( infile, "cols:\t%d\n", &iCols );
	iReturn = fscanf( infile, "rows:\t%d\n", &iRows );
	iReturn = fscanf( infile, "format:\t%s\n", pcFormat );
	fclose( infile );

	// Open data
	sprintf( pcFileName, "%s.bin", pcFileNameBase );
	infile = fopen( pcFileName, "rb" );
	if( !infile )
	{
		return -1;
	}
	float *pdDataX = new float[iRows*iCols];
	if( !pdDataX )
	{
		fclose( infile );
		return -1;
	}
	if( pcFormat[0] == 'i' )
	{
		for( int i = 0; i < iRows*iCols; i++ )
		{
			int iData;
			fread( &iData, sizeof(int), 1, infile );
			pdDataX[i] = (float)iData;
		}
	}
	else if( pcFormat[0] == 'f' )
	{
		for( int i = 0; i < iRows*iCols; i++ )
		{
			double dData;
			fread( &dData, sizeof(double), 1, infile );
			pdDataX[i] = (float)dData;
		}
	}
	fclose( infile );

	float **pdX = new float*[iRows];
	for( int i = 0; i < iRows; i++ )
	{
		pdX[i] = &pdDataX[i*iCols];
	}

	*pppdX = pdX;

	return 0;
}


//
// readArrayNoAllocation()
//
// Read the array from file, allocation of pdX has taken place elsewhere
//
int
readArrayNoAllocation(
		  char *pcFileNameBase,
		  float *pdX,
		  int &iCols,
		  int &iRows,
		  char *pcFormat
		  )
{
	FILE *infile;
	char pcFileName[300];
	int iReturn;
	// Open header
	sprintf( pcFileName, "%s.txt", pcFileNameBase );
	infile = fopen( pcFileName, "rt" );
	if( !infile )
	{
		return -1;
	}
	iReturn = fscanf( infile, "cols:\t%d\n", &iCols );
	iReturn = fscanf( infile, "rows:\t%d\n", &iRows );
	iReturn = fscanf( infile, "format:\t%s\n", pcFormat );
	fclose( infile );

	// Open data
	sprintf( pcFileName, "%s.bin", pcFileNameBase );
	infile = fopen( pcFileName, "rb" );
	if( !infile )
	{
		return -1;
	}

	if( pcFormat[0] == 'i' )
	{
		for( int i = 0; i < iRows*iCols; i++ )
		{
			int iData;
			fread( &iData, sizeof(int), 1, infile );
			pdX[i] = (float)iData;
		}
	}
	else if( pcFormat[0] == 'f' )
	{
		for( int i = 0; i < iRows*iCols; i++ )
		{
			double dData;
			fread( &dData, sizeof(double), 1, infile );
			pdX[i] = (float)dData;
		}
	}
	fclose( infile );

	return 0;
}

//
// readArraySize()
//
// Just to size up the arrays.
//
int
readArraySize(
		  char *pcFileNameBase,
		  int &iCols,
		  int &iRows,
		  char *pcFormat
		  )
{
	FILE *infile;
	char pcFileName[300];
	int iReturn;

	// Open header
	sprintf( pcFileName, "%s.txt", pcFileNameBase );
	infile = fopen( pcFileName, "rt" );
	if( !infile )
	{
		return -1;
	}
	iReturn = fscanf( infile, "cols:\t%d\n", &iCols );
	iReturn = fscanf( infile, "rows:\t%d\n", &iRows );
	iReturn = fscanf( infile, "format:\t%s\n", pcFormat );
	fclose( infile );
	return 0;
}

int
read_point_data_US_multiple(
							int argc,
							char **argv,
				   my_constraint_data &dat
)
{
	char pcFileName[400];
	int iCols;
	int iRows;
	int iDataCount = 0;
	char pcFormat[10];

	// Pass1: figure out how much space is necessary
	dat.iPoints = 0;
	for( int i = 0; i < argc; i++ )
	{
		sprintf( pcFileName, "%s.non-linear.matches.X", argv[i] );
		if( readArraySize( pcFileName, iCols, iRows, pcFormat ) < 0 )
		{
			return -1;
		}
		dat.iPoints += iRows;
		iDataCount += iCols*iRows;
	}

	// Allocate arrays
	dat.pfPoints1 = new float[iDataCount];
	dat.pfPoints2 = new float[iDataCount];
	if( !dat.pfPoints1 || !dat.pfPoints2 )
	{
		return  -1;
	}

	// Pass2: read in
	int iCurrCount = 0;
	for( int i = 0; i < argc; i++ )
	{
		sprintf( pcFileName, "%s.non-linear.matches.X", argv[i] );
		readArrayNoAllocation( pcFileName, dat.pfPoints1+iCurrCount, iCols, iRows, pcFormat );

		sprintf( pcFileName, "%s.non-linear.matches.Y", argv[i] );
		readArrayNoAllocation( pcFileName, dat.pfPoints2+iCurrCount, iCols, iRows, pcFormat );

		iCurrCount += iCols*iRows;
	}

	//sprintf( pcFileName, "%s.non-linear.matches.X", pcFNameBase );readArray( pcFileName, &p1, iCols, iRows, pf1 );
	//sprintf( pcFileName, "%s.non-linear.matches.Y", pcFNameBase );readArray( pcFileName, &p2, iCols, iRows, pf2 );
	return 1;
}


// Experimental - add offset to points in the image plane to 
// change center of rotation / reduce optimizer bias?
// It seems the optimizer is choosing a translation vs a rotation.
int
add_point_offset(
				   my_constraint_data &dat,
				   int iColOffset=315, // Center of ultrasound fan in MNI data
				   int iRowOffset=315
				   )
{
	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		float *p_in = point_in+12;
		float *p_out = point_out+12;

		// Center around 
		p_in[0] -= iColOffset;
		p_in[1] -= iRowOffset;
		p_out[0] -= iColOffset;
		p_out[1] -= iRowOffset;
	}
	return 1;
}


int
read_point_data_US(
				   char *pcFNameBase,
				   my_constraint_data &dat
)
{
	// Read data from files
	int iRows, iCols;
	float ** p1, ** p2;
	char pf1[2], pf2[2];

	char pcFileName[400];
	sprintf( pcFileName, "%s.non-linear.matches.X", pcFNameBase );readArray( pcFileName, &p1, iCols, iRows, pf1 );
	sprintf( pcFileName, "%s.non-linear.matches.Y", pcFNameBase );readArray( pcFileName, &p2, iCols, iRows, pf2 );

	//readArray( "C:\\writting\\cvs_stuff\\papers\\miccai-2013\\pics\\reconstruction\\non-linear.matches.X", &p1, iCols, iRows, pf1 );
	//readArray( "C:\\writting\\cvs_stuff\\papers\\miccai-2013\\pics\\reconstruction\\non-linear.matches.Y", &p2, iCols, iRows, pf2 );

	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\non-linear.matches.X", &p1, iCols, iRows, pf1 );
	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\non-linear.matches.Y", &p2, iCols, iRows, pf2 );

	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\mr-non-linear.matches.X", &p1, iCols, iRows, pf1 );
	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\mr-non-linear.matches.Y", &p2, iCols, iRows, pf2 );

	dat.iPoints = iRows;
	dat.pfPoints1 = p1[0];
	dat.pfPoints2 = p2[0];
	return 1;
}

//
// read_point_data_US_frame_indices()
//
// *** 2014 3D reconstruction
//
int
read_point_data_US_frame_indices(
				   char *pcFNameBase,
				   my_constraint_data &dat
)
{
	// Read data from files
	int iRows, iCols;
	float ** p1, ** p2;
	char pf1[2], pf2[2];

	char pcFileName[400];
	sprintf( pcFileName, "%s.non-linear.matches.X", pcFNameBase );readArray( pcFileName, &p1, iCols, iRows, pf1 );
	sprintf( pcFileName, "%s.non-linear.matches.Y", pcFNameBase );readArray( pcFileName, &p2, iCols, iRows, pf2 );

	//readArray( "C:\\writting\\cvs_stuff\\papers\\miccai-2013\\pics\\reconstruction\\non-linear.matches.X", &p1, iCols, iRows, pf1 );
	//readArray( "C:\\writting\\cvs_stuff\\papers\\miccai-2013\\pics\\reconstruction\\non-linear.matches.Y", &p2, iCols, iRows, pf2 );

	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\non-linear.matches.X", &p1, iCols, iRows, pf1 );
	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\non-linear.matches.Y", &p2, iCols, iRows, pf2 );

	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\mr-non-linear.matches.X", &p1, iCols, iRows, pf1 );
	//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\mr-non-linear.matches.Y", &p2, iCols, iRows, pf2 );

	dat.iPoints = iRows;
	dat.pfPoints1 = p1[0];
	dat.pfPoints2 = p2[0];
	return 1;
}

int
read_point_data_US(
				   char *pcFNameBase1,
				   char *pcFNameBase2,
				   my_constraint_data &dat
)
{
	my_constraint_data dat1;
	my_constraint_data dat2;
	read_point_data_US( pcFNameBase1, dat1 );
	read_point_data_US( pcFNameBase2, dat2 );

	dat.iPoints = dat1.iPoints + dat2.iPoints;
	dat.pfPoints1 = new float[15*dat.iPoints];
	dat.pfPoints2 = new float[15*dat.iPoints];

	for( int i = 0; i < dat1.iPoints*15; i++ )
	{
		dat.pfPoints1[i] = dat1.pfPoints1[i];
		dat.pfPoints2[i] = dat1.pfPoints2[i];
	}
	for( int i = 0; i < dat2.iPoints*15; i++ )
	{
		dat.pfPoints1[dat1.iPoints*15+i] = dat2.pfPoints1[i];
		dat.pfPoints2[dat1.iPoints*15+i] = dat2.pfPoints2[i];
	}

	delete [] dat1.pfPoints1;
	delete [] dat1.pfPoints2;
	delete [] dat2.pfPoints1;
	delete [] dat2.pfPoints2;

	return 1;
}


int
read_point_data_US_intersection(
				   char *pcFNameBase1,
				   char *pcFNameBase2,
				   my_constraint_data &dat
)
{
	read_point_data_US( pcFNameBase1, dat );

	my_constraint_data dat2;
	read_point_data_US( pcFNameBase2, dat2 );

	// Hash everything into a code to identify duplicates to be kept
	// The idea is to keep the intersection of points in dat & dat2, which 
	// should be the most reliable correspondences (bipartitle nearest neighbors)
	map<float,int> mapMatchDuplicateFinder;
	for( int i = 0; i < dat2.iPoints; i++ )
	{
		float *points21 = dat2.pfPoints1 + 15*i;
		float *points22 = dat2.pfPoints2 + 15*i;
		float fKey = 1;
		for( int j = 0; j < 15; j++ )
		{
			float fNum = points21[j]*points22[j]+1;
			if( fNum != 0 )
			{
				fKey *= fNum;
			}
		}

		mapMatchDuplicateFinder.insert( pair<float,int>(fKey,i) );
	}

	int iSavePoint = 0;
	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *points21 = dat.pfPoints1 + 15*i;
		float *points22 = dat.pfPoints2 + 15*i;
		float fKey = 1;
		for( int j = 0; j < 15; j++ )
		{
			float fNum = points21[j]*points22[j]+1;
			if( fNum != 0 )
			{
				fKey *= fNum;
			}
		}

		// Test this hash code - to avoid multiple matches between precisely the same locations
		map<float,int>::iterator itAmazingDupeFinder = mapMatchDuplicateFinder.find(fKey);
		if( itAmazingDupeFinder != mapMatchDuplicateFinder.end() )
		{
			// Save
			for( int j = 0; j < 15; j++ )
			{
				dat.pfPoints1[15*iSavePoint+j] = dat.pfPoints1[15*i+j];
				dat.pfPoints2[15*iSavePoint+j] = dat.pfPoints2[15*i+j];
			} 
			iSavePoint++;
		}
	}

	dat.iPoints = iSavePoint;

	// Delete secondary array 
	delete [] dat2.pfPoints1;
	delete [] dat2.pfPoints2;

	return 1;
}



//
// Fitness function for rotation data
//
double 
nlopt_func_Fitness(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector
	float vec2[3];
	vec2[0] = x[0];
	vec2[1] = x[1];
	vec2[2] = x[2];
	float mat[9];
	convert_vector_to_rotation_matrix( vec2, mat );

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 3*i;
		float *point_out = dat.pfPoints2 + 3*i;

		float point_out_test[3];

		// Multiply input point by current rotation matrix
		mult_3x3_vector( mat, point_in, point_out_test );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_out, point_out_test );

		// Sum up 
		dSumSqrError += fErrorSqr;
	}

	return dSumSqrError;
}

//
// Fitness function for scaled rigid transform
//
//double 
//nlopt_func_Fitness_rigid(
//		   unsigned n,
//		   const double *x,
//			double *gradient, /* NULL if not needed */
//			void *func_data
//			)
//{
//	double dSumSqrError = 0;
//
//    my_constraint_data &dat = *(my_constraint_data *) func_data;
//
//	// Create rotation matrix from current vector
//	float vec2[PARAMETER_DIM_RIG];
//	float mat[9];
//
//	for( int i = 0; i < 7; i++ )
//		vec2[i] = x[i];
//
//	// Create rotation matrix
//	convert_vector_to_rotation_matrix( vec2, mat );
//	// Multiply scale
//	mult_3x3_scalar( mat, exp(vec2[6]), mat );
//
//	for( int i = 0; i < dat.iPoints; i++ )
//	{
//		float *point_in = dat.pfPoints1 + 3*i;
//		float *point_out = dat.pfPoints2 + 3*i;
//
//		float point_out_test[3];
//
//		// Multiply input point by current rotation matrix
//		mult_3x3_vector( mat, point_in, point_out_test );
//
//		// Add translation
//		point_out_test[0] += vec2[3];
//		point_out_test[1] += vec2[4];
//		point_out_test[2] += vec2[5];
//
//		// Compute squared error
//		float fErrorSqr = vec3D_distsqr_3d( point_out, point_out_test );
//
//		// Sum up 
//		dSumSqrError += fErrorSqr;
//	}
//
//	return dSumSqrError;
//}

int 
_compareFloat(const void *v1, const void *v2) 
{
	float *pf1 = (float*)v1;
	float *pf2 = (float*)v2;
	if( *pf1 < *pf2 )
	{
		return -1;
	}
	else
	{
		return *pf1 > *pf2;
	}
}

//
// Fitness function for scaled rigid transform of ultrasound data
//
// Here, points in my_constraint_data are 15 dimensional vectors.
// storing the first 3 rows of a 3D transform matrix (12), then input point (3)
//
//double 
//nlopt_func_Fitness_rigid_US(
//		   unsigned n,
//		   const double *x,
//			double *gradient, /* NULL if not needed */
//			void *func_data
//			)
//{
//	double dSumSqrError = 0;
//
//    my_constraint_data &dat = *(my_constraint_data *) func_data;
//
//	// Create rotation matrix from current vector (length 7)
//	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
//	float vec2[PARAMETER_DIM_RIG];
//	float mat[9];
//
//	for( int i = 0; i < 7; i++ )
//		vec2[i] = x[i];
//
//	// Create rotation matrix
//	convert_vector_to_rotation_matrix( vec2, mat );
//	// Multiply scale
//	//mult_3x3_scalar( mat, exp(vec2[6]), mat );
//
//	// Copy scaled rotation matrix to 4x4 homogenous matrix
//	float mat44[16];
//	memset( &mat44[0], 0, sizeof(mat44) );
//	memcpy( &mat44[0], &mat[0], 3*sizeof(float) );
//	memcpy( &mat44[4], &mat[3], 3*sizeof(float) );
//	memcpy( &mat44[8], &mat[6], 3*sizeof(float) );
//
//	// Copy translations to 4x4 homegenous matrix
//	mat44[3] = vec2[3];
//	mat44[7] = vec2[4];
//	mat44[11] = vec2[5];
//	mat44[15] = 1;
//
//	// Two US tracker transform matrices
//	float mat44_1[16];
//	float mat44_2[16];
//	// Set last row to: [0,0,0,1]
//	memset( &mat44_1[0], 0, sizeof(mat44_1) );
//	memset( &mat44_2[0], 0, sizeof(mat44_2) );
//	mat44_1[15] = 1;
//	mat44_2[15] = 1;
//
//	// Two product transform matrices
//	float mat44_1_p[16];
//	float mat44_2_p[16];
//
//	float *pfDistSqr = new float[dat.iPoints];
//
//	for( int i = 0; i < dat.iPoints; i++ )
//	{
//		float *point_in = dat.pfPoints1 + 15*i;
//		float *point_out = dat.pfPoints2 + 15*i;
//
//		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
//		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
//		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );
//
//		// Multiply current calibration matrix
//		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
//		mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
//
//		// Generate output 3D test points
//		float point_test_in[3];
//		float point_test_out[3];
//		mult_4x4_vector_homogenous( mat44_1_p, point_in+12, point_test_in );
//		mult_4x4_vector_homogenous( mat44_2_p, point_out+12, point_test_out );
//
//		// Compute squared error
//		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );
//
//		// Sum up 
//		dSumSqrError += fErrorSqr;
//		pfDistSqr[i] = fErrorSqr;
//	}
//
//	// Robust: keep top half of distance 
//	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
//	dSumSqrError = 0;
//	float fPercentToKeep = 3.0/4.0;
//	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
//	{
//		dSumSqrError += pfDistSqr[i];
//	}
//
//	delete [] pfDistSqr;
//
//	return dSumSqrError;
//}

//
// nlopt_func_Fitness_rigid_US_scaled()
//
// Original US autocalibration function (e.g. MICCAI 2013)
//
double 
nlopt_func_Fitness_rigid_US_scaled(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];

	for( int i = 0; i < 7; i++ )
		vec2[i] = x[i];

	float fScale = 1;
	if( n == 7 )
	{
		fScale = exp(vec2[6]);
	}

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	convert_vector_to_matrix4x4( vec2, mat44 );

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

		// Point 2: multiply current calibration matrix
		float point_test_out[3];
		if( dat.iErrorCalib )
		{
			printf( "Error: this code has been changed & untested\n" );
			//// Pass Point 1 through ground truth calibration, for evaluating calibratin error
			//mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			//mult_4x4_vector_homogenous( mat44_2_p, p_in, point_test_out );

			//// Alternative 1 - pass point 2 though ground truth calibration
			//mult_4x4_matrix( mat44_2, &dat.pfCalib[0], mat44_2_p );
			//mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

			//// Alternative 2 - for shifted image points, reverse shift,
			//// pass through ground truth. 
			//float point_out_shifted[3];
			//point_out_shifted[0] = p_in[0] + 315;
			//point_out_shifted[1] = p_in[1] + 315;
			//point_out_shifted[2] = p_in[2];
			//mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			//mult_4x4_vector_homogenous( mat44_2_p, point_out_shifted, point_test_out );
		}
		else
		{
			// Pass Point 2 through current calibration matrix, for autocalibration optimization
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );
		}

		// Rescale to pixel units??
		//vec3D_mult_scalar( point_test_in, 1.0f/fScale );
		//vec3D_mult_scalar( point_test_out, 1.0f/fScale );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		for( int k = 0; k < 3; k++ )
			pdAbsErrorXYZ[k] = fabs(point_test_in[k] - point_test_out[k]);

		// Sum up 
		dSumSqrError += fErrorSqr;
		pfDistSqr[i] = fErrorSqr;
	}

	// Robust: keep top half of distance 
	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
	dSumSqrError = 0;
	//float fPercentToKeep = 3.0/4.0;
	float fPercentToKeep = dat.fPercentToKeep;
	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
	{
		dSumSqrError += sqrt(pfDistSqr[i]);
	}

	// Original return value
	float fReturnValueMM = dSumSqrError / (fPercentToKeep*dat.iPoints);

	return fReturnValueMM;



	// This here is all experimental stuff
	//    didn't really pan out ...  yet
	dSumSqrError = 1;
	for( int i = 0; i < dat.iPoints; i++ )
	{
		float fDist = sqrt(pfDistSqr[i]);
		if( fDist < 4.0f )
		{
			dSumSqrError += (4.0f-fDist);
		}
		else
		{
			dSumSqrError += 10;
		}
	}

	delete [] pfDistSqr;

	float fReturnValueRANSAC = dSumSqrError;//1.0f / dSumSqrError;

	return fReturnValueRANSAC;
}

//
// nlopt_func_Fitness_rigid_US_scaled_test_alignment()
//
// Modified nlopt_func_Fitness_rigid_US_scaled_test() to 
// include case where data are misaligned.
//
double 
nlopt_func_Fitness_rigid_US_scaled_test_alignment(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];

	for( int i = 0; i < 7; i++ )
 		vec2[i] = x[i];

	float fScale = 1;
	if( n == 7 )
	{
		fScale = exp(vec2[6]);
	}

	//vec2[0] = 0;	vec2[1] = 0;	vec2[2] = 0;

	// Copy scaled rotation matrix to 4x4 homogenous matrix

	float mat44[16];
	float mat44_align[16];
	convert_vector_to_matrix4x4( dat.pfCalib, mat44 );
	convert_vector_to_matrix4x4( dat.pfAlign, mat44_align );

	// Initialze
	if( dat.bPerformAlign )
	{
		convert_vector_to_matrix4x4( vec2, mat44_align );
	}
	else
	{
		convert_vector_to_matrix4x4( vec2, mat44 );
	}

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];
	float mat44_2_p_align[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

		// Point 2: multiply current calibration matrix, left multiply 3D alignment matrix.
		float point_test_out[3];
		mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
		mult_4x4_matrix( mat44_align, mat44_2_p, mat44_2_p_align );
		mult_4x4_vector_homogenous( mat44_2_p_align, p_out, point_test_out );

		//mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
		//mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		for( int k = 0; k < 3; k++ )
			pdAbsErrorXYZ[k] = fabs(point_test_in[k] - point_test_out[k]);

		// Sum up 
		dSumSqrError += fErrorSqr;
		pfDistSqr[i] = fErrorSqr;
	}

	// Robust: keep top half of distance 
	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
	dSumSqrError = 0;
	//float fPercentToKeep = 3.0/4.0;
	float fPercentToKeep = dat.fPercentToKeep;
	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
	{
		dSumSqrError += sqrt(pfDistSqr[i]);
	}

	// Original return value
	float fReturnValueMM = dSumSqrError / (fPercentToKeep*dat.iPoints);

	return fReturnValueMM;
}






//
// nlopt_func_Fitness_rigid_US_scaled_centered()
//
// This function was derived from nlopt_func_Fitness_rigid_US_scaled(), the goal 
// is to compute 
//
double 
nlopt_func_Fitness_rigid_US_scaled_centered(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];

	for( int i = 0; i < 7; i++ )
		vec2[i] = x[i];

	float fScale = 1;
	if( n == 7 )
	{
		fScale = exp(vec2[6]);
	}

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	convert_vector_to_matrix4x4( vec2, mat44 );

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

		// Point 2: multiply current calibration matrix
		float point_test_out[3];
		if( dat.iErrorCalib )
		{
			// Pass Point 1 through ground truth calibration, for evaluating calibratin error
			mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_in, point_test_out );

			// Alternative 1 - pass point 2 though ground truth calibration
			mult_4x4_matrix( mat44_2, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

			// Alternative 2 - for shifted image points, reverse shift,
			// pass through ground truth. 
			float point_out_shifted[3];
			point_out_shifted[0] = p_in[0] + 315;
			point_out_shifted[1] = p_in[1] + 315;
			point_out_shifted[2] = p_in[2];
			mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, point_out_shifted, point_test_out );
		}
		else
		{
			// Pass Point 2 through current calibration matrix, for autocalibration optimization
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );
		}

		// Rescale to pixel units??
		//vec3D_mult_scalar( point_test_in, 1.0f/fScale );
		//vec3D_mult_scalar( point_test_out, 1.0f/fScale );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		for( int k = 0; k < 3; k++ )
			pdAbsErrorXYZ[k] = fabs(point_test_in[k] - point_test_out[k]);

		// Sum up 
		dSumSqrError += fErrorSqr;
		pfDistSqr[i] = fErrorSqr;
	}

	// Robust: keep top half of distance 
	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
	dSumSqrError = 0;
	//float fPercentToKeep = 3.0/4.0;
	float fPercentToKeep = dat.fPercentToKeep;
	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
	{
		dSumSqrError += sqrt(pfDistSqr[i]);
	}

	delete [] pfDistSqr;

	return dSumSqrError / (fPercentToKeep*dat.iPoints);
}



//
// nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances()
//
// This is Sandy's idea, to test inter-point distances in 3D.
//	Goal: see if distances between points remain constant, to see if transform is rigid.
//
double 
nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];

	for( int i = 0; i < 7; i++ )
		vec2[i] = x[i];

	float fScale = 1;
	if( n == 7 )
	{
		fScale = exp(vec2[6]);
	}

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	convert_vector_to_matrix4x4( vec2, mat44 );

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	float *pfPoints3D = new float[3*dat.iPoints];

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

		mult_4x4_vector_homogenous( mat44_1_p, p_in, pfPoints3D+3*i );
	}

	// Now compute inter-point distances for a random subset of (20 points)
	FILE *outfile = fopen( "interpoint_distance.txt", "wt" );
	for( int i = 0; i < dat.iPoints; i += dat.iPoints/20 )
	{
		for( int j = i; j < dat.iPoints; j += dat.iPoints/20 )
		{
			float fErrorSqr = vec3D_distsqr_3d( pfPoints3D+3*i, pfPoints3D+3*j );
			fprintf( outfile, "%f\t", fErrorSqr );
		}
		fprintf( outfile, "\n" );
	}
	fclose( outfile );

	return 1.0;
}

//
// nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances_to_ground_truth()
//
// Here try to figure out what precisely is the error compared to ground truth (identity) calibration matrix.
//
double 
nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances_to_ground_truth(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data,
			char *pcOutFile
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];

	for( int i = 0; i < 7; i++ )
		vec2[i] = x[i];

	float fScale = 1;
	if( n == 7 )
	{
		fScale = exp(vec2[6]);
	}

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	convert_vector_to_matrix4x4( vec2, mat44 );

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	float *pfPoints3D1 = new float[3*dat.iPoints];
	float *pfPoints3D2 = new float[3*dat.iPoints];

	FILE *outfile = fopen( pcOutFile, "a+" );
	for( int i = 0; i < 16; i++ )
	{
		if( i % 4 == 0 )
		{
			fprintf( outfile, "\n" );
		}
		fprintf( outfile, "%f\t", mat44[i] );
	}
	fprintf( outfile, " ------ \n" );

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );


		// Point 2: multiply current calibration matrix
		float point_test_out[3];

		// Pass Point 2 through current calibration matrix, for autocalibration optimization
		mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
		mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );


		fprintf( outfile, "%f\t%f\t%f\n%f\t%f\t%f\n",
			point_test_in[0], point_test_in[1], point_test_in[2], 
			point_test_out[0], point_test_out[1], point_test_out[2]
			);
	}

	fclose( outfile );

	return 1.0;
}





//
// nlopt_func_Fitness_rigid_US_scaled_plane()
//
// This function was develop to determine the scaled rigid transform
// from 2D points in an ultrasound plane to 3D points in the world.
// June 10 2013.
// Here, argument *x is the transformation matrix for a particular
// ultrasound plane. In func_data (my_constraint_data), the
//
double 
nlopt_func_Fitness_rigid_US_scaled_plane(
		   unsigned n,
		   const double *x,
			double *gradient, /* NULL if not needed */
			void *func_data
			)
{
	double dSumSqrError = 0;

    my_constraint_data &dat = *(my_constraint_data *) func_data;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG];
	float mat[9];
	for( int i = 0; i < 7; i++ )
		vec2[i] = x[i];
	// Create rotation matrix
	convert_vector_to_rotation_matrix( vec2, mat );
	// Multiply scale - pretend it doesn't change
	//mult_3x3_scalar( mat, exp(vec2[6]), mat );
	float fScale = 1;//exp(vec2[6]);

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	memset( &mat44[0], 0, sizeof(mat44) );
	memcpy( &mat44[0], &mat[0], 3*sizeof(float) );
	memcpy( &mat44[4], &mat[3], 3*sizeof(float) );
	memcpy( &mat44[8], &mat[6], 3*sizeof(float) );
	// Copy translations to 4x4 homegenous matrix
	mat44[3] = vec2[3];
	mat44[7] = vec2[4];
	mat44[11] = vec2[5];
	mat44[15] = 1;

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	for( int i = 0; i < dat.iPoints; i++ )
	{
		// Note: each data record is 15 elements long (4x3=12 for transformation matrix, 3 for an xyz point).
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			// Determine scaled 3D point
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 ); // ***
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Copy in current estimated transformation matrix - should hopefully converge
		// to something like what is currently there ...
		// This is actually redundant from the step above *** ...
		memcpy( &mat44_1[0], mat44, sizeof(float)*16 );

		// Multiply current calibration matrix - embedded in data
		mult_4x4_matrix( mat44_1, dat.pfCalib, mat44_1_p );
		mult_4x4_matrix( mat44_2, dat.pfCalib, mat44_2_p );

		// Generate output 3D test points
		float point_test_in[3];
		float point_test_out[3];
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );
		mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		// Sum up 
		dSumSqrError += fErrorSqr;
		pfDistSqr[i] = fErrorSqr;
	}

	// Robust: keep top half of distance 
	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
	dSumSqrError = 0;
	//float fPercentToKeep = 3.0/4.0;
	float fPercentToKeep =0.75;
	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
	{
		dSumSqrError += sqrt(pfDistSqr[i]);
	}

	delete [] pfDistSqr;

	if( fPercentToKeep*dat.iPoints > 0 )
	{
		dSumSqrError /= fPercentToKeep*dat.iPoints;
	}
	return dSumSqrError;
}


int
main_test_rotation_minimization(
								)
{
		nlopt_stopping nlStopCond;


	double xtolabs[PARAMETER_DIM_ROT] = {0,0,0};

	memset( &nlStopCond, 0, sizeof(nlStopCond) );
	nlStopCond.xtol_abs = &xtolabs[0];
	nlStopCond.xtol_rel = -1;
	nlStopCond.ftol_abs = -1;
	nlStopCond.ftol_rel = -1;
	nlStopCond.maxeval = 300; // Go 300 iterations
	nlStopCond.n = PARAMETER_DIM_ROT;

	my_constraint_data dat;

	generate_minimization_test_data( dat, 1000 );

    int n = PARAMETER_DIM_ROT;
	double x_params[PARAMETER_DIM_ROT] = {1,1,1};
	double x_lb[PARAMETER_DIM_ROT] = {-10,-10,-10};
	double x_ub[PARAMETER_DIM_ROT] = {10,10,10};
	double x_step[PARAMETER_DIM_ROT] = {5,5,5};
	double minf;
	nldrmd_minimize(
		n, nlopt_func_Fitness, &dat,			    
							 &x_lb[0],
							 &x_ub[0],
							 &x_params[0],
							 &minf,
			     
							 &x_step[0], //const double *xstep, /* initial step sizes */
			     
							 &nlStopCond
							 );

	return 1;
}


//
// premultiply_data()
//
// Multiply data by a test matrix.
//
// 
//
int
random_multiply_data(
	my_constraint_data &dat,
	int iParamCount,
	int bPostMultiply = 1, // 0 for pre-multiply (alignment), 1 for post-multiply (calibration)
	int iSide = 0 // If 1, then only premultiple pfPoints1, if 2 then only pfPoints2, if 0 (default) then both
				 )
{
	float x_params[PARAMETER_DIM_RIG] = {1,2,3,0,0,0,2};

	float fNum;

	//x_params[0] = 30*(0.5f-(rand()/((float)RAND_MAX)));
	//x_params[1] = 30*(0.5f-(rand()/((float)RAND_MAX)));
	//x_params[2] = 30*(0.5f-(rand()/((float)RAND_MAX)));	
	x_params[3] = 200*2*(0.5f-(rand()/((float)RAND_MAX)));
	x_params[4] = 200*2*(0.5f-(rand()/((float)RAND_MAX)));
	x_params[5] = 200*2*(0.5f-(rand()/((float)RAND_MAX)));
	x_params[0] = 10*(0.5f-rand()/((float)RAND_MAX));
	x_params[1] = 10*(0.5f-rand()/((float)RAND_MAX));
	x_params[2] = 10*(0.5f-rand()/((float)RAND_MAX));
	if( iParamCount > 6 )
	{
		//x_params[6] = 2*(0.5f-rand()/((float)RAND_MAX));
		//x_params[6] = 2.3*2.0*(0.5f-rand()/((float)RAND_MAX));
		x_params[6] = 2.3*2.0*(0.5f-rand()/((float)RAND_MAX));
	}
	else
	{
		x_params[6] = 0;
	}

	printf( "Random pre-multiply: one-side %d:\n%f %f %f %f %f %f %f\n", 
		iSide,
		x_params[0],
		x_params[1],
		x_params[2],
		x_params[3],
		x_params[4],
		x_params[5],
		x_params[6]
	);

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float mat[9];

	// Create scaled rotation matrix
	convert_vector_to_rotation_matrix( x_params, mat );
	float fScale = exp(x_params[6]);
	mult_3x3_scalar( mat, fScale, mat );

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	memset( &mat44[0], 0, sizeof(mat44) );
	memcpy( &mat44[0], &mat[0], 3*sizeof(float) );
	memcpy( &mat44[4], &mat[3], 3*sizeof(float) );
	memcpy( &mat44[8], &mat[6], 3*sizeof(float) );
	// Copy translations to 4x4 homegenous matrix
	mat44[3] = x_params[3];
	mat44[7] = x_params[4];
	mat44[11] = x_params[5];
	mat44[15] = 1;

	// Two US tracker transform matrices - initialize constant entries
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		if( bPostMultiply ) // Post-multiply, random calibration maxtrix, typically applied to both sides to test calibration
		{
			// Multiply current calibration matrix
			mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
		}
		else // Pre-multiply, random alignment matrix, typically applied to one side to test alignment
		{
			// Multiply current calibration matrix
			mult_4x4_matrix( mat44, mat44_1, mat44_1_p );
			mult_4x4_matrix( mat44, mat44_2, mat44_2_p );
		}

		// Write out transformed 3 rows (length 3x4 = 12) of US tracker transform matrices
		if( iSide == 0 || iSide == 1 )
			memcpy( point_in,  &mat44_1_p[0], sizeof(float)*12 );
		if( iSide == 0 || iSide == 2 )
			memcpy( point_out, &mat44_2_p[0], sizeof(float)*12 );
	}

	return 1;
}


int
main_test_rigid_US_minimization(
	my_constraint_data &dat,
	float &fResultError,
	int bRandomStartPoint = 0 // For a random starting point
)
{
	nlopt_stopping nlStopCond;

	double xtolabs[PARAMETER_DIM_RIG] = {0,0,0,0,0,0,0};

	memset( &nlStopCond, 0, sizeof(nlStopCond) );
	nlStopCond.xtol_abs = &xtolabs[0]; // tolerances for individual  parameters, for assessing convergence
	nlStopCond.xtol_rel = -1; // tolerance
	nlStopCond.ftol_abs = -1;
	nlStopCond.ftol_rel = -1;
	nlStopCond.maxeval = 10000; // Go 300 iterations
	nlStopCond.n = PARAMETER_DIM_RIG;

	// Calibration matrix for one good AMIGO calibration
	float pfCalib[] = {0.104828101,	0.000981117,	0.004546451,	-36.313503,	0.004569771,	-0.002627669,	-0.104798743,	-1.8013641,	-0.000866027,	0.104893735,	-0.002667809,	-1.0628991,	0,	0,	0,	1 };
	// Calibration matrix: identity
	float pfCalibIdenity[] = 
	{1, 0, 0, 0, 
		0, 1, 0, 0, 
		0, 0, 1, 0, 
	0, 0, 0, 1 };

	// Set input data fields for optimization
	dat.fPercentToKeep = 0.2;//3.0/4.0;
	dat.iErrorCalib = 0;
	dat.bPerformAlign = 0;
	// Set to 
	memset( dat.pfAlign, 0, sizeof(dat.pfAlign) );
	memset( dat.pfCalib, 0, sizeof(dat.pfCalib) );

	//float vec[3] = {-.3,-.2,.1};
	//float mat[3*3] = {1,0,0,0,1,0,0,0,1};
	//convert_vector_to_rotation_matrix( vec, mat );
	//convert_rotation_matrix_to_vector( vec, mat );


    int n = PARAMETER_DIM_RIG;

	//
	// *** Note: we need to specify proper/effective bounds
	//   
	//
	//double x_params[PARAMETER_DIM_RIG] = {1,-4,3,0,0,0,1};
	double x_params[PARAMETER_DIM_RIG] = {0,0.0,0.0,0,0,0,0};
	//double x_lb[PARAMETER_DIM_RIG] = {-10,-10,-10,-500,-500,-500,-5};
	//double x_ub[PARAMETER_DIM_RIG] = { 10, 10, 10, 500, 500, 500, 5};
	double x_lb[PARAMETER_DIM_RIG] = {-10,-10,-10,-200,-200,-200,-5};
	double x_ub[PARAMETER_DIM_RIG] = { 10, 10, 10, 200, 200, 200, 5};
	double x_step[PARAMETER_DIM_RIG] = {2,2,2,100,100,100,1};

	//double x_params[PARAMETER_DIM_RIG] = {0.01,0.01,0.01,0,0,0,0};
	//double x_lb[PARAMETER_DIM_RIG] = {-7,-7,-7,-20000,-20000,-2000,-3};
	//double x_ub[PARAMETER_DIM_RIG] = { 7, 7, 7, 20000, 20000, 2000, 3};
	//double x_step[PARAMETER_DIM_RIG] = {1,1,1,1000,1000,1000,1};

	// Control wether last scale parameter is considered
	//n = 7;
	n = 6;
	
	// **********************************************
	// This is for testing random calibration matrices
	//		Comment this out when testing ground truth error
	//
	// random_multiply_data( dat, n, 1 ); // For testing calibration
	//
	// random_multiply_data( dat, n, 0, 1 ); // For testing alignment

	dat.bPerformAlign = 1;
	float fError;
	fError = nlopt_func_Fitness_rigid_US_scaled( n, &x_params[0],0,&dat);
	fError = nlopt_func_Fitness_rigid_US_scaled_test_alignment( n, &x_params[0],0,&dat);

	//x_params[4] = 20000;
	//fError = nlopt_func_Fitness_rigid_US_scaled( n, &x_params[0],0,&dat);

	if( bRandomStartPoint )
	{
		// Initialize algorithm starting point to somewhere within search range
		for( int i = 0; i < PARAMETER_DIM_RIG; i++ )
		{
			float fRange = x_ub[i]-x_lb[i];
			float fRand = rand()/(RAND_MAX+1.0f);
			x_params[i] = x_lb[i]+fRand*fRange;
		}
	}

	double minf;
	nlopt_result res;

	for( int i = 0; i < 100; i++ )
	{
		// 1) Perform alignment - should be independent of calibration

		dat.bPerformAlign = 0;
		dat.fPercentToKeep = 0.7;

		//n = 6;

		nlStopCond.nevals = 0;


		res = 
		nldrmd_minimize(
		//sbplx_minimize(
			n,
			nlopt_func_Fitness_rigid_US_scaled_test_alignment,
			//nlopt_func_Fitness_rigid_US_scaled, // ** US paper results
			//nlopt_func_Fitness_rigid_US,
			&dat,			    
								 &x_lb[0],
								 &x_ub[0],
								 &x_params[0],
								 &minf,
				     
								 &x_step[0], //const double *xstep, /* initial step sizes */
				     
								 &nlStopCond
								 ); 

		// See how it works ...
		fError = nlopt_func_Fitness_rigid_US_scaled( n, &x_params[0],0,&dat);
		fError = nlopt_func_Fitness_rigid_US_scaled_test_alignment( n, &x_params[0],0,&dat);
		fResultError = fError;

		//     save alignment matrix, perform calibration
		if( dat.bPerformAlign )
		{
			for( int j = 0; j < PARAMETER_DIM_RIG; j++ ) dat.pfAlign[j] = x_params[j];
		}
		else
		{
			for( int j = 0; j < PARAMETER_DIM_RIG; j++ ) dat.pfCalib[j] = x_params[j];
		}

		// 2) Perform calibration - should be independent of alignment, if alignment is correct...
		dat.bPerformAlign = 0;
		dat.fPercentToKeep = 0.2;

		res = 
		nldrmd_minimize(
		//sbplx_minimize(
			n,
			nlopt_func_Fitness_rigid_US_scaled_test_alignment,
			//nlopt_func_Fitness_rigid_US_scaled, // ** US paper results
			//nlopt_func_Fitness_rigid_US,
			&dat,			    
								 &x_lb[0],
								 &x_ub[0],
								 &x_params[0],
								 &minf,
				     
								 &x_step[0], //const double *xstep, /* initial step sizes */
				     
								 &nlStopCond
								 ); 

		//		save calibration matrix for alignment
		//     save alignment matrix, perform calibration
		if( dat.bPerformAlign )
		{
			for( int j = 0; j < PARAMETER_DIM_RIG; j++ ) dat.pfAlign[j] = x_params[j];
		}
		else
		{
			for( int j = 0; j < PARAMETER_DIM_RIG; j++ ) dat.pfCalib[j] = x_params[j];
		}

		// See how it works ...
		fError = nlopt_func_Fitness_rigid_US_scaled( n, &x_params[0],0,&dat);
		fError = nlopt_func_Fitness_rigid_US_scaled_test_alignment( n, &x_params[0],0,&dat);
		printf( "Iteration: %d\tError: %f\n", i, fError );
		printf( "Calib:\n%f %f %f %f %f %f %f\n", 
			dat.pfCalib[0],
			dat.pfCalib[1],
			dat.pfCalib[2],
			dat.pfCalib[3],
			dat.pfCalib[4],
			dat.pfCalib[5],
			dat.pfCalib[6]
		);
		printf( "Align:\n%f %f %f %f %f %f %f\n", 
			dat.pfAlign[0],
			dat.pfAlign[1],
			dat.pfAlign[2],
			dat.pfAlign[3],
			dat.pfAlign[4],
			dat.pfAlign[5],
			dat.pfAlign[6]
		);
		fResultError = fError;

	}



	// Test distances between transformed points - sanity check
	// nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances( n, &x_params[0],0,&dat);

	// Output points according to identified transform
	nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances_to_ground_truth( n, &x_params[0],0,&dat, "transformed_points_calibration.txt" );


	// **********************************************
	// Test ground truth calibration
	//
	// Works when random_postmultiply_data( dat, n ) above is disabled/commented out
	// This should be commented out when running randomized calibration,
	// because the ground truth transform is no longer identity.
	//
	set_vector_to_identity( &x_params[0] );
	fError = nlopt_func_Fitness_rigid_US_scaled( n, &x_params[0],0,&dat);
	//x_params[3] = 315;x_params[4] = 315;

	// Test distances between transformed points - sanity check
	// nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances( n, &x_params[0],0,&dat);

	// Output points according to identity transform
	nlopt_func_Fitness_rigid_US_scaled_compute_transformed_point_distances_to_ground_truth( n, &x_params[0],0,&dat, "transformed_points_identity.txt" );


	float fScale = exp( x_params[6] );

	return 1;
}

//
// main_test_plane_fitting()
//
// The idea here is to test fitting of individual ultrasound planes.
// June, 2013.
//
int
main_test_plane_fitting(
						)
{
	nlopt_stopping nlStopCond;

	double xtolabs[PARAMETER_DIM_RIG] = {0,0,0,0,0,0,0};

	memset( &nlStopCond, 0, sizeof(nlStopCond) );
	nlStopCond.xtol_abs = &xtolabs[0]; // tolerances for individual  parameters, for assessing convergence
	nlStopCond.xtol_rel = -1; // tolerance
	nlStopCond.ftol_abs = -1;
	nlStopCond.ftol_rel = -1;
	nlStopCond.maxeval = 10000; // Go 300 iterations
	nlStopCond.n = PARAMETER_DIM_RIG;

	my_constraint_data dat;
	// Calibration matrix for one good AMIGO calibration
	float pfCalib[] = {0.104828101,	0.000981117,	0.004546451,	-36.313503,	0.004569771,	-0.002627669,	-0.104798743,	-1.8013641,	-0.000866027,	0.104893735,	-0.002667809,	-1.0628991,	0,	0,	0,	1 };
	//dat.pfCalib = pfCalib;
	memset( dat.pfCalib, 0, sizeof(dat.pfCalib) );

	if( 1 )
	{
		// Read data from files
		int iRows, iCols;
		float ** p1, ** p2;
		char pf1[2], pf2[2];

		readArray( "C:\\writting\\cvs_stuff\\papers\\miccai-2013\\pics\\reconstruction\\non-linear.matches.X", &p1, iCols, iRows, pf1 );
		readArray( "C:\\writting\\cvs_stuff\\papers\\miccai-2013\\pics\\reconstruction\\non-linear.matches.Y", &p2, iCols, iRows, pf2 );

		//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\non-linear.matches.X", &p1, iCols, iRows, pf1 );
		//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\non-linear.matches.Y", &p2, iCols, iRows, pf2 );

		//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\mr-non-linear.matches.X", &p1, iCols, iRows, pf1 );
		//readArray( "C:\\Visual Studio Projects\\libs\\OpenCV2.2\\windows\\samples\\cpp\\mr-non-linear.matches.Y", &p2, iCols, iRows, pf2 );

		dat.iPoints = iRows;
		dat.pfPoints1 = p1[0];
		dat.pfPoints2 = p2[0];
	}
	else
	{
		generate_minimization_test_data_rigid_US( dat, 200 );
	}

    int n = PARAMETER_DIM_RIG;

	double x_params[PARAMETER_DIM_RIG] = {0.01,0.01,0.01,0,0,0,0};
	double x_lb[PARAMETER_DIM_RIG] = {-10,-10,-10,-200,-200,-200,-5};
	double x_ub[PARAMETER_DIM_RIG] = { 10, 10, 10, 200, 200, 200, 5};
	double x_step[PARAMETER_DIM_RIG] = {1,1,1,5,5,5,1};

	FILE *outfile = fopen( "correspondences-error.txt", "wt" );

	int iFrame = 0;

	for( int iIndex = 0; iIndex < dat.iPoints; iIndex++ )
	{
		// Start of current testing point
		int iStartIndex = iIndex;

		// Collect all correspondences with similar transform matrices
		while( iIndex < dat.iPoints
			&& memcmp( dat.pfPoints1+iStartIndex*15, dat.pfPoints1+iIndex*15, sizeof(float)*12 ) == 0
			)
		{
			iIndex++;
		}

		// Init optimization parameters close to pre-existing transform

		my_constraint_data dat_one_slice;
		dat_one_slice.iPoints = iIndex-iStartIndex;
		dat_one_slice.pfPoints1 = dat.pfPoints1 + 15*iStartIndex;
		dat_one_slice.pfPoints2 = dat.pfPoints2 + 15*iStartIndex;
		memcpy( dat_one_slice.pfCalib, dat.pfCalib, sizeof(dat.pfCalib) );

		// Determine plane

		n = 6; // Forget scale parameter for now

		double minf;
		nlopt_result res = 
			nldrmd_minimize(
			n,
			nlopt_func_Fitness_rigid_US_scaled_plane,
			&dat_one_slice,			    
			&x_lb[0],
			&x_ub[0],
			&x_params[0],
			&minf,

			&x_step[0], //const double *xstep, /* initial step sizes */

			&nlStopCond
			); 

		float pfPointIn[3] = {320,240,0};
		float pfPointOut1[3];
		float pfPointOut2[3];


		float vec[3];
		float mat44[16];
		vec[0] = x_params[0];
		vec[1] = x_params[1];
		vec[2] = x_params[2];
		convert_vector_to_rotation_matrix_4X4( vec, mat44 );
		mat44[3] = x_params[3];
		mat44[7] = x_params[4];
		mat44[11] = x_params[5];

		// Generate test point from estimated matrix
		mult_4x4_vector_homogenous( mat44, pfPointIn, pfPointOut1 );

		// Generate test point from ground truth matrix
		// Copy in ground truth
		memcpy( mat44, dat.pfPoints1+iStartIndex*15, sizeof(float)*12 );
		mult_4x4_vector_homogenous( mat44, pfPointIn, pfPointOut2 );

		float fErrorPlane = vec3D_dist_3d( pfPointOut1, pfPointOut2 );

		// Test error
		float fRobustError = nlopt_func_Fitness_rigid_US_scaled_plane( n, &x_params[0],0,&dat_one_slice);

		fprintf( outfile, "%d\t%d\t%f\t%f\n", iFrame, iIndex-iStartIndex, fRobustError, fErrorPlane );
		iFrame++;
	}

	fclose( outfile );

	return 1;
}


typedef struct _KEYPOINT
{
	float x;
	float y;
	float size;
	float angle;
	float response;
	int octave;
	int id;
} KEYPOINT;

typedef float DESCRIPTOR[64];

//
// read_in_features()
//
// Description:
//	Reads in arrays from two binary files:
//		1) *.surf.key : keypoints
//		2) *.surf.dsc : descriptors
//
// Input:
//	Text file name for ultrasound sequence/directory, e.g. "_02\pre\sweep_2a\2d\2a.txt"
//
// Outputs:
//  None for the moment, only internal data structures are filed up:
//		int				iImages; // Number of images
//		int			*	piKeysPerImage; // Array of keypoint counts for each image
//		KEYPOINT	**	ppkKeypoints;  // Array of keypoints for each image
//		float		**	ppfDescriptors; // Array of descriptors for each image
//
int
read_in_features(
		char *pcInputFile // "_02\\pre\\sweep_2a\\2d\\2a.txt"
		)
{
	// Data structures to read in 
	int				iImages; // Number of images
	int			*	piKeysPerImage; // Array of keypoint counts for each image
	KEYPOINT	**	ppkKeypoints;  // Array of keypoints for each image
	float		**	ppfDescriptors; // Array of descriptors for each image


	FILE *infile;
	int iReturn;
	int iKeypoints;


	char pcFileName[400];

	// Identify extension 
	//pcInputFile = "C:\\downloads\\data\\MNI\\bite-us\\group1\\_02\\pre\\sweep_2a\\2d\\2a.txt";

	sprintf( pcFileName, pcInputFile );

	char *pch = strstr( pcFileName, ".txt" );

	// Open keypoints file
	sprintf( pch, "%s", ".surf.key" );
	infile = fopen( pcFileName, "rb" );
	if( !infile )
	{
		printf( "Error: no key file found\n" );
		return -1;
	}
	fread( &iImages, sizeof(iImages), 1, infile ); // number of images

	ppkKeypoints = new KEYPOINT*[iImages];
	piKeysPerImage = new int[iImages];

	for( int i = 0; i < iImages; i++ )
	{
		fread( &iKeypoints, sizeof(iKeypoints), 1, infile ); // number of keypoints perimage
		ppkKeypoints[i] = new KEYPOINT[iKeypoints];
		piKeysPerImage[i] = iKeypoints;
		for( int j = 0; j < iKeypoints; j++ )
		{
			fread( &(ppkKeypoints[i][j]), sizeof(KEYPOINT), 1, infile );
		}
	}
	fclose( infile );

	// Open descriptor file
	sprintf( pch, "%s", ".surf.dsc" );
	infile = fopen( pcFileName, "rb" );
	fread( &iImages, sizeof(iImages), 1, infile ); // number of images, should be the same as previous
	ppfDescriptors = new float*[iImages];
	for( int i = 0; i < iImages; i++ )
	{
		int iRows;
		int iCols;
		int iType;
		fread( &iRows, sizeof(iRows), 1, infile ); // number of keypoints perimage
		fread( &iCols, sizeof(iCols), 1, infile ); // Should always be 64
		fread( &iType, sizeof(iType), 1, infile ); // Should always be 5
		if( iRows > 0 && iCols > 0 )
		{
			ppfDescriptors[i] = new float[iRows*iCols];
			iReturn = fread( ppfDescriptors[i], sizeof(float), iRows*iCols, infile );
			if( iReturn != iRows*iCols )
			{
				printf( "Something went wrong.\n" );
			}
		}
		else
		{
			printf( "Zero feature image: %d\n", i );
		}
	}
	fclose( infile );

	return 0;
}















template <class T, class DIV>
void
mult_3x3_matrix_template(
			T mat1[3][3],
			T mat2[3][3],
			T mat_out[3][3]
			)
{
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			mat_out[i][j]=0;
			for( int ii = 0; ii < 3; ii++ )
			{
				mat_out[i][j] += mat1[i][ii]*mat2[ii][j];
			}
		}
	}
}

template <class T, class DIV>
void
mult_3x3_template(
			T mat[3][3],
			T vec_in[3],
			T vec_out[3]
			)
{
	for( int i = 0; i < 3; i++ )
	{
		vec_out[i] = 0;
		for( int j = 0; j < 3; j++ )
		{
			vec_out[i] += mat[i][j]*vec_in[j];
		}
	}
}



template <class T, class DIV>
void
transpose_3x3(
			T mat_in[3][3],
			T mat_trans[3][3]
			)
{
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			mat_trans[i][j] = mat_in[j][i];
		}
	}
}

//
// determine_rotation_3point()
//
// Determines a 3D rotation matrix from 3 points to the cartesian frame.
// The first point pfP01 is at the origin (0,0,0)=(x,y,z)
// Vector pfP01->pfP02 defines the x-axis (1,0,0)=(x,y,z)
// The cross product of vectors pfP01->pfP02 and pfP01->pfP03 (pfV0_nm) defines the z-axis (0,0,1)=(x,y,z)
// The cross product of vectors pfP01->pfP02 and pfV0_nm definex the y-axis (0,1,0)=(x,y,z)
//
int
determine_rotation_3point(
					   float *pfP01, float *pfP02, float *pfP03, // points in image 1 (3d)
					   float *rot // rotation matrix  (3x3)
					   )
{
	float pfV0_12[3];
	float pfV0_13[3];
	float pfV0_nm[3];

	// Subtract point 1 to convert to vectors
	vec3D_diff_3d( pfP01, pfP02, pfV0_12 );
	vec3D_diff_3d( pfP01, pfP03, pfV0_13 );

	// Normalize vectors
	vec3D_norm_3d( pfV0_12 );
	vec3D_norm_3d( pfV0_13 );

	// Cross product between 2 vectors to get normal 
	vec3D_cross_3d( pfV0_12, pfV0_13, pfV0_nm );
	vec3D_norm_3d(  pfV0_nm );

	// Cross product between 1st vectors and normal to get orthogonal 3rd vector
	vec3D_cross_3d( pfV0_nm, pfV0_12, pfV0_13 );
	vec3D_norm_3d(  pfV0_13 );

	rot[0]=pfV0_12[0];
	rot[1]=pfV0_12[1];
	rot[2]=pfV0_12[2];

	rot[3]=pfV0_13[0];
	rot[4]=pfV0_13[1];
	rot[5]=pfV0_13[2];

	rot[6]=pfV0_nm[0];
	rot[7]=pfV0_nm[1];
	rot[8]=pfV0_nm[2];

	return 1;
}



//
// determine_similarity_transform_3point()
//
// Determines a 3 point similarity transform in 3D.
//
// Inputs:
//    Corresponding points
//		float *pfP01, float *pfP02, float *pfP03,
//		float *pfP11, float *pfP12, float *pfP13,
//
// Outputs
//		rot0: rotation matrix, about the first pair of corresponding points pfP01<->pfP11.
//      fScaleDiff: scale change.
//		(rot1: used internally)
//
int
determine_similarity_transform_3point(
					   float *pfP01, float *pfP02, float *pfP03,
					   float *pfP11, float *pfP12, float *pfP13,

					   // Output
					   // Rotation about point pfP01/pfP11
					   float *rot0,
					   float *rot1,

					   // Scale change (magnification) from image 1 to image 2
					   float &fScaleDiff
					   )
{
	// Points cannot be equal or colinear

	float fDist012 = vec3D_dist_3d( pfP01, pfP02 );
	float fDist013 = vec3D_dist_3d( pfP01, pfP03 );
	float fDist023 = vec3D_dist_3d( pfP02, pfP03 );
	float fDist112 = vec3D_dist_3d( pfP11, pfP12 );
	float fDist113 = vec3D_dist_3d( pfP11, pfP13 );
	float fDist123 = vec3D_dist_3d( pfP12, pfP13 );

	if(
		fDist012 == 0 ||
		fDist013 == 0 ||
		fDist023 == 0 ||
		fDist112 == 0 ||
		fDist113 == 0 ||
		fDist123 == 0 )
	{
		// Points have to be distinct
		return -1;
	}

	fScaleDiff = (fDist112+fDist113+fDist123)/(fDist012+fDist013+fDist023);

	// Determine rotation about point pfP01/pfP11

	int iReturn;

	iReturn = determine_rotation_3point( pfP01, pfP02, pfP03, rot0 );
	iReturn = determine_rotation_3point( pfP11, pfP12, pfP13, rot1 );

	// Create 
	float ori0[3][3];
	memcpy( &(ori0[0][0]), rot0, sizeof(float)*3*3 );
	float ori1[3][3];
	memcpy( &(ori1[0][0]), rot1, sizeof(float)*3*3 );

	float ori1trans[3][3];
	transpose_3x3<float,float>( ori1, ori1trans );

	float rot_one[3][3];
	mult_3x3_matrix_template<float,float>( ori1trans, ori0, rot_one );

	memcpy( rot0, &(rot_one[0][0]), sizeof(float)*3*3 );

	return iReturn;
}

//
// angular_difference()
//
// Measure angular difference.
//
float
consistency_angle(
					   float *pfP01, float *pfP02, float *pfP03,
					   float *pfP11, float *pfP12, float *pfP13
					   )
{
	float fDist012 = vec3D_dist_3d( pfP01, pfP02 );
	float fDist013 = vec3D_dist_3d( pfP01, pfP03 );

	float fDist023 = vec3D_dist_3d( pfP02, pfP03 );

	float fDist112 = vec3D_dist_3d( pfP11, pfP12 );
	float fDist113 = vec3D_dist_3d( pfP11, pfP13 );

	float fDist123 = vec3D_dist_3d( pfP12, pfP13 );
	
	if(
		fDist012 == 0 ||
		fDist013 == 0 ||
		fDist023 == 0 ||
		fDist112 == 0 ||
		fDist113 == 0 ||
		fDist123 == 0 )
	{
		// Points have to be distinct
		return -1;
	}

	float pfV012[3];
	float pfV013[3];
	vec3D_diff_3d( pfP01, pfP02, pfV012 );
	vec3D_diff_3d( pfP01, pfP03, pfV013 );
	vec3D_norm_3d(pfV012);
	vec3D_norm_3d(pfV013);

	float pfV112[3];
	float pfV113[3];
	vec3D_diff_3d( pfP11, pfP12, pfV112 );
	vec3D_diff_3d( pfP11, pfP13, pfV113 );
	vec3D_norm_3d(pfV112);
	vec3D_norm_3d(pfV113);

}

//
//  consistency_side_scale_difference()
//
// Perfect consistency is 0, side legth ratios are exact.
// Imperfect consistency higher or lower than 0.
//
float
consistency_side_scale_difference(
					   float *pfP01, float *pfP02, float *pfP03,
					   float *pfP11, float *pfP12, float *pfP13
					   )
{
	// Points cannot be equal or colinear

	float fDist012 = vec3D_dist_3d( pfP01, pfP02 );
	float fDist013 = vec3D_dist_3d( pfP01, pfP03 );
	float fDist023 = vec3D_dist_3d( pfP02, pfP03 );

	float fDist112 = vec3D_dist_3d( pfP11, pfP12 );
	float fDist113 = vec3D_dist_3d( pfP11, pfP13 );
	float fDist123 = vec3D_dist_3d( pfP12, pfP13 );

	if(
		fDist012 == 0 ||
		fDist013 == 0 ||
		fDist023 == 0 ||
		fDist112 == 0 ||
		fDist113 == 0 ||
		fDist123 == 0 )
	{
		// Points have to be distinct
		return -10000;
	} 

	float fR0 = fDist012/fDist013;
	float fR1 = fDist112/fDist113;

	return log(fR1/fR0);
}


void
vec3D_summ_3d(
			 float *pv1,
			 float *pv2,
			 float *pv12
			 )
{
	pv12[0] = pv2[0]+pv1[0];
	pv12[1] = pv2[1]+pv1[1];
	pv12[2] = pv2[2]+pv1[2];
}

int
similarity_transform_3point(
							float *pfP0,
							float *pfP1,

							float *pfCenter0,
							float *pfCenter1,
						float *rot01,
					   float fScaleDiff
					   )
{
	// Subtract origin from point
	float pfP0Diff[3];
	vec3D_diff_3d( pfCenter0, pfP0, pfP0Diff );

	// Rotate from image 0 to image 1
	float rot_one[3][3];
	memcpy( &(rot_one[0][0]), rot01, sizeof(float)*3*3 );
	mult_3x3_template<float,float>( rot_one, pfP0Diff, pfP1 );

	// Scale
	pfP1[0] *= fScaleDiff;
	pfP1[1] *= fScaleDiff;
	pfP1[2] *= fScaleDiff;

	// Add to origin
	vec3D_summ_3d( pfP1, pfCenter1, pfP1 );

	return 0;
}


int
determine_similarity_transform_ransac(
									  float *pf0, // Points in image 1
									  float *pf1, // Points in image 2
									  float *pfs0, // Scale for points in image 1
									  float *pfs1, // Scale for points in image 2
									  int iPoints, // number of points
						int iIterations,
						float fDistThreshold,
						int *piInlierFlags
						
					//	,									  
					//float fScaleDiffThreshold,
					//float fShiftThreshold,
					//float fCosineAngleThreshold
	)
{
	if( iPoints < 4 )
	{
		return 0;
	}

	int iMaxInliers = 0;
	float fMinErrorSum = 1000000;

	for( int i = 0; i < iIterations; i++ )
	{
		// Find a set of unique points
		int i1 = iPoints*(rand()/(RAND_MAX+1.0f));
		int i2 = iPoints*(rand()/(RAND_MAX+1.0f));
		int i3 = iPoints*(rand()/(RAND_MAX+1.0f));
		while( i1 == i2 || i1 == i3 || i2 == i3 )
		{
			i1 = iPoints*(rand()/(RAND_MAX+1.0f));
			i2 = iPoints*(rand()/(RAND_MAX+1.0f));
			i3 = iPoints*(rand()/(RAND_MAX+1.0f));
		}

		// Calculate transform
		float rot_here0[3][3];
		float rot_here1[3][3];
		float fScaleDiffHere;

		float fConsistencyScale = consistency_side_scale_difference(
					   pf0+i1*3, pf0+i2*3, pf0+i3*3,
					   pf1+i1*3, pf1+i2*3, pf1+i3*3
					   );

		//float fConsistencyAngle = consistency_side_scale_difference(
		//	   pf0+i1*3, pf0+i2*3, pf0+i3*3,
		//	   pf1+i1*3, pf1+i2*3, pf1+i3*3
		//	   );

		int iResult = determine_similarity_transform_3point(
					   pf0+i1*3, pf0+i2*3, pf0+i3*3,
					   pf1+i1*3, pf1+i2*3, pf1+i3*3,
					   (float*)(&(rot_here0[0][0])),
					   (float*)(&(rot_here1[0][0])),
					   fScaleDiffHere
					   );
		if( iResult < 0 )
		{
			// Didn't work
			continue;
		}

		float fErrorSum = 0;

		// Count the number of inliers
		int iInliers = 0;
		int iInputInliers = 0;
		float pfTest[3];
		for( int j = 0; j < iPoints; j++ )
		{
			similarity_transform_3point(
				pf0+j*3, pfTest,
				pf0+i1*3, pf1+i1*3,
				(float*)(&(rot_here0[0][0])), fScaleDiffHere );

			// Compute squared error
			float fErrorSqr = vec3D_distsqr_3d( pf1, pfTest );

			float fDist = sqrt( fErrorSqr );
			if( fDist < fDistThreshold )
			{
				fErrorSum += fDist;
				iInliers++;
				if( j == i1 || j == i2 || j == i3 )
				{
					// Count indices of input features
					//	- are the features used to generate the transform part of the solution?
					// Not necessarily in noisy error cases. The computed similarity transform 
					// does not necessarily operate well in noise conditions, but it should in valid conditions.
					iInputInliers++;
				}
			}
		}

		if( iInputInliers != 3 )
		{
			// Didn't work - similarity transform must
			// correctly map inliers used to compute it
			continue;
		}

		if( iMaxInliers < iInliers )
		{
			iMaxInliers = iInliers;
			fMinErrorSum = fErrorSum;

			// Set flags
			float pfTest[3];
			for( int j = 0; j < iPoints; j++ )
			{
				similarity_transform_3point(
					pf0+j*3, pfTest,
					pf0+i1*3, pf1+i1*3,
					(float*)(&(rot_here0[0][0])), fScaleDiffHere );

				// Compute squared error
				float fErrorSqr = vec3D_distsqr_3d( pf1, pfTest );

				float fDist = sqrt( fErrorSqr );
				if( fDist < fDistThreshold )
				{
					piInlierFlags[j]=1;
				}
				else
				{
					piInlierFlags[j] = 0;
				}
			}

			assert( piInlierFlags[i1]==1);
			assert( piInlierFlags[i2]==1);
			assert( piInlierFlags[i3]==1);
			piInlierFlags[i1]=2;
			piInlierFlags[i2]=3;
			piInlierFlags[i3]=4;
		}
		else if( iMaxInliers == iInliers && fErrorSum < fMinErrorSum && iInliers >= 3 )
		{
			iMaxInliers = iInliers;
			fMinErrorSum = fErrorSum;

				// Set flags
			float pfTest[3];
			for( int j = 0; j < iPoints; j++ )
			{
				similarity_transform_3point(
					pf0+j*3, pfTest,
					pf0+i1*3, pf1+i1*3,
					(float*)(&(rot_here0[0][0])), fScaleDiffHere );

				// Compute squared error
				float fErrorSqr = vec3D_distsqr_3d( pf1, pfTest );

				float fDist = sqrt( fErrorSqr );
				if( fDist < fDistThreshold )
				{
					piInlierFlags[j]=1;
				}
				else
				{
					piInlierFlags[j] = 0;
				}
			}
			
			assert( piInlierFlags[i1]==1);
			assert( piInlierFlags[i2]==1);
			assert( piInlierFlags[i3]==1);
			piInlierFlags[i1]=2;
			piInlierFlags[i2]=3;
			piInlierFlags[i3]=4;
		}
	}

	return iMaxInliers;
}




//
// test_navigation_error()
//
// Testing navigation error for MICCAI 2014
//
double 
test_navigation_error(
					  char *pcFileName,
					  my_constraint_data &dat
			)
{
	double dSumSqrError = 0;

	dat.iErrorCalib = 0;
	dat.fPercentToKeep = 1;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG] = {0,0,0,0,0,0,0};
	float mat[9];

	float fScale = 1;

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	convert_vector_to_matrix4x4( vec2, mat44 );

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	vector<int> vecFrameIndices;

	// 
	// Here we gather statistics about correspondences in data
	//
	// Identify correspondences part of the same solution
	int iStartCountIndex = 0;
	vecFrameIndices.push_back(iStartCountIndex);

	int *piInlierCountPerFeature = new int[dat.iPoints];
	int *piInlierFrequency = new int[dat.iPoints];
	float *pfDotProductZ = new float[dat.iPoints];
	memset( piInlierCountPerFeature, 0, sizeof(int)*dat.iPoints );
	memset( piInlierFrequency, 0, sizeof(int)*dat.iPoints );
	float fVecMag1 = dat.pfPoints1[ 2]*dat.pfPoints1[ 2] + dat.pfPoints1[ 6]*dat.pfPoints1[ 6] + dat.pfPoints1[10]*dat.pfPoints1[10];
	float fVecMag2 = dat.pfPoints2[ 2]*dat.pfPoints2[ 2] + dat.pfPoints2[ 6]*dat.pfPoints2[ 6] + dat.pfPoints2[10]*dat.pfPoints2[10];
	float fVecNorm = sqrt( fVecMag1*fVecMag2 );
	for( int i = 0; i < dat.iPoints; i++ )
	{
		// Compute dot product with z-vectors
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;
		float fDot =
			point_in[ 2]*point_out[ 2] +
			point_in[ 6]*point_out[ 6] +
			point_in[10]*point_out[10];
		pfDotProductZ[i] = fDot/fVecNorm;

		if( dat.pfPoints1[15*i] != dat.pfPoints1[15*iStartCountIndex] )
		{
			int iInlierCount = i-iStartCountIndex;
			piInlierFrequency[iInlierCount]++;
			for( int j = iStartCountIndex; j < i; j++ )
			{
				piInlierCountPerFeature[j] = iInlierCount;
			}
			iStartCountIndex = i;
			vecFrameIndices.push_back(iStartCountIndex);
		}
	}

	// Finish up
	int iInlierCount = dat.iPoints-iStartCountIndex;
	piInlierFrequency[dat.iPoints-iStartCountIndex]++;
	for( int j = iStartCountIndex; j < dat.iPoints; j++ )
	{
		piInlierCountPerFeature[j] = dat.iPoints-iStartCountIndex;
	}

	FILE *outXFormFit = fopen( "_transform_fit.txt", "a+" );

	//
	// Here we create vectors of points in 3D in the pre-operative model.
	//
	//
	vector< float > vecPoints3DImage1;
	vector< float > vecDotProduct;
	vector< float > vecPoints3DImage2;
	vector< float > vecPoints3DError;
	
	// Keep cooridinate errors
	vector< float > vecErrorX;
	vector< float > vecErrorY;
	vector< float > vecErrorZ;

	for( int i = 0; i < dat.iPoints; i++ )
	{
		if( piInlierCountPerFeature[i] <= 3 )
		{
			continue;
		}

		vecPoints3DImage1.clear();
		vecPoints3DImage2.clear();
		vecPoints3DError.clear();
		vecDotProduct.clear();
		vecErrorX.clear();
		vecErrorY.clear();
		vecErrorZ.clear();

		// Valid plane here, check error

		int iStartPointIndex = vecPoints3DImage2.size();
		int iInlierCount = piInlierCountPerFeature[i];



		for( int ik = i; ik < i + iInlierCount; ik++ )
		{
			float *point_in = dat.pfPoints1 + 15*ik;
			float *point_out = dat.pfPoints2 + 15*ik;

			// Copy in scaled points
			float p_in[3];
			float p_out[3];
			for( int j = 0; j < 3; j++ )
			{
				p_in[j] = fScale*point_in[12+j];
				p_out[j] = fScale*point_out[12+j];
			}

			// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
			memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
			memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

			// Point 1: multiply current calibration matrix
			float point_test_in[3];
			mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
			mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

			// Pass Point 2 through current calibration matrix
			float point_test_out[3];
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

			// Compute error: for debugging
			float fErrorSqr = sqrt( vec3D_distsqr_3d( point_test_in, point_test_out ) );

			vecErrorX.push_back( point_test_out[0]-point_test_in[0] );
			vecErrorY.push_back( point_test_out[1]-point_test_in[1] );
			vecErrorZ.push_back( point_test_out[2]-point_test_in[2] );

			vecDotProduct.push_back( pfDotProductZ[ik] );
			vecPoints3DError.push_back( fErrorSqr );
			for( int j = 0; j < 3; j++ )
			{
				vecPoints3DImage1.push_back( p_in[j] );				// In-plane pixel coordinates
				vecPoints3DImage2.push_back( point_test_out[j] );	// 3D coordinates
			}
			
		}

		int piInlierFlags[200];
		memset( piInlierFlags, 0, sizeof(piInlierFlags) );

		float pfPlane[4];
		// Estimate plane & error
		int iInliers;
		
		//iInliers = 
		//ransac_plane_estimation(
		//	iInlierCount,
		//	&vecPoints3DImage2[iStartPointIndex],
		//	100,
		//	1.5f,
		//	pfPlane );

		iInliers =
		determine_similarity_transform_ransac(			
			&vecPoints3DImage1[iStartPointIndex],
			&vecPoints3DImage2[iStartPointIndex],
			0, 0, iInlierCount, 100, 3.0f, piInlierFlags );

		if( iInliers >= 3 )
		{
			static int iFrameCount = 0;

			for( int j = 0; j < iInlierCount; j++ )
			{
				if( piInlierFlags[j] )
				{
					fprintf( outXFormFit, "%s\t%d\t%f\t%f\n", strrchr( pcFileName, '\\')+1, iFrameCount, vecPoints3DError[j], vecDotProduct[j] );
				}
			}			

			iInliers =
			determine_similarity_transform_ransac(			
				&vecPoints3DImage1[iStartPointIndex],
				&vecPoints3DImage2[iStartPointIndex],
				0, 0, iInlierCount, 100, 3.0f, piInlierFlags );
			iFrameCount++;
		}

		// Skip past current inliers
		i += iInlierCount-1;
	}

	fclose( outXFormFit );

	return 1;

	// Now loop through all frames and compute planes
	for( int i = 0; i < vecFrameIndices.size(); i++ )
	{
		float pfPlane[4];
		int iFrameStartIndex = vecFrameIndices[i];
		ransac_plane_estimation(
			piInlierCountPerFeature[iFrameStartIndex],
			&vecPoints3DImage2[3*iFrameStartIndex],
			100,
			1.5f,
			pfPlane );
	}

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

		// Point 2: multiply current calibration matrix
		float point_test_out[3];
		if( dat.iErrorCalib )
		{
			// Pass Point 1 through ground truth calibration, for evaluating calibratin error
			mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_in, point_test_out );

			// Alternative 1 - pass point 2 though ground truth calibration
			mult_4x4_matrix( mat44_2, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

			// Alternative 2 - for shifted image points, reverse shift,
			// pass through ground truth. 
			float point_out_shifted[3];
			point_out_shifted[0] = p_in[0] + 315;
			point_out_shifted[1] = p_in[1] + 315;
			point_out_shifted[2] = p_in[2];
			mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, point_out_shifted, point_test_out );
		}
		else
		{
			// Pass Point 2 through current calibration matrix, for autocalibration optimization
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );
		}

		// Rescale to pixel units??
		//vec3D_mult_scalar( point_test_in, 1.0f/fScale );
		//vec3D_mult_scalar( point_test_out, 1.0f/fScale );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		for( int k = 0; k < 3; k++ )
			pdAbsErrorXYZ[k] = fabs(point_test_in[k] - point_test_out[k]);

		// Sum up 
		dSumSqrError += fErrorSqr;
		pfDistSqr[i] = fErrorSqr;
	}

	FILE *outf = fopen( "inliers.txt", "a+" );
	//fprintf( outf, "InlierCountInImage\tErrorThisFeature\tDotProductZ\tInlierFrequency (not related to individual features)\n" );
	for( int i = 0; i < dat.iPoints; i++ )
	{
		if( piInlierCountPerFeature[i] > 3 )
		{
			fprintf( outf, "%s\t%d\t%f\t%f\t%d\n", strrchr( pcFileName, '\\'), piInlierCountPerFeature[i], sqrt(pfDistSqr[i]), pfDotProductZ[i], piInlierFrequency[i] );
		}
	}
	fclose( outf );

	delete [] pfDotProductZ;
	delete [] piInlierCountPerFeature;
	delete [] piInlierFrequency;

	// Robust: keep top half of distance 
	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
	dSumSqrError = 0;
	//float fPercentToKeep = 3.0/4.0;
	float fPercentToKeep = dat.fPercentToKeep;
	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
	{
		dSumSqrError += sqrt(pfDistSqr[i]);
	}



	delete [] pfDistSqr;

	return dSumSqrError / (fPercentToKeep*dat.iPoints);
}



//
// test_trackerless_reconstruction()
//
// This function was scavanged from test_navigation_error()
//
// test_trackerless_reconstruction
//
double 
test_trackerless_reconstruction(
					  char *pcFileName,
					  my_constraint_data &dat
			)
{
	double dSumSqrError = 0;

	dat.iErrorCalib = 0;
	dat.fPercentToKeep = 1;

	// Create rotation matrix from current vector (length 7)
	// Parameters: rotation vector (3), translation vector (3), isotropic scale (1)
	float vec2[PARAMETER_DIM_RIG] = {0,0,0,0,0,0,0};
	float mat[9];

	float fScale = 1;

	// Copy scaled rotation matrix to 4x4 homogenous matrix
	float mat44[16];
	convert_vector_to_matrix4x4( vec2, mat44 );

	// Two US tracker transform matrices
	float mat44_1[16];
	float mat44_2[16];
	// Set last row to: [0,0,0,1]
	memset( &mat44_1[0], 0, sizeof(mat44_1) );
	memset( &mat44_2[0], 0, sizeof(mat44_2) );
	mat44_1[15] = 1;
	mat44_2[15] = 1;

	// Two product transform matrices
	float mat44_1_p[16];
	float mat44_2_p[16];

	float *pfDistSqr = new float[dat.iPoints];

	double pdAbsErrorXYZ[3];
	pdAbsErrorXYZ[0] = pdAbsErrorXYZ[1] = pdAbsErrorXYZ[2] = 0;

	vector<int> vecFrameIndices;

	// 
	// Here we gather statistics about correspondences in data
	//
	// Identify correspondences part of the same solution
	int iStartCountIndex = 0;
	vecFrameIndices.push_back(iStartCountIndex);

	int *piInlierCountPerFeature = new int[dat.iPoints];
	int *piInlierFrequency = new int[dat.iPoints];
	float *pfDotProductZ = new float[dat.iPoints];
	memset( piInlierCountPerFeature, 0, sizeof(int)*dat.iPoints );
	memset( piInlierFrequency, 0, sizeof(int)*dat.iPoints );
	float fVecMag1 = dat.pfPoints1[ 2]*dat.pfPoints1[ 2] + dat.pfPoints1[ 6]*dat.pfPoints1[ 6] + dat.pfPoints1[10]*dat.pfPoints1[10];
	float fVecMag2 = dat.pfPoints2[ 2]*dat.pfPoints2[ 2] + dat.pfPoints2[ 6]*dat.pfPoints2[ 6] + dat.pfPoints2[10]*dat.pfPoints2[10];
	float fVecNorm = sqrt( fVecMag1*fVecMag2 );
	for( int i = 0; i < dat.iPoints; i++ )
	{
		// Compute dot product with z-vectors
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;
		float fDot =
			point_in[ 2]*point_out[ 2] +
			point_in[ 6]*point_out[ 6] +
			point_in[10]*point_out[10];
		pfDotProductZ[i] = fDot/fVecNorm;

		if( dat.pfPoints1[15*i] != dat.pfPoints1[15*iStartCountIndex] )
		{
			int iInlierCount = i-iStartCountIndex;
			piInlierFrequency[iInlierCount]++;
			for( int j = iStartCountIndex; j < i; j++ )
			{
				piInlierCountPerFeature[j] = iInlierCount;
			}
			iStartCountIndex = i;
			vecFrameIndices.push_back(iStartCountIndex);
		}
	}

	// Finish up
	int iInlierCount = dat.iPoints-iStartCountIndex;
	piInlierFrequency[dat.iPoints-iStartCountIndex]++;
	for( int j = iStartCountIndex; j < dat.iPoints; j++ )
	{
		piInlierCountPerFeature[j] = dat.iPoints-iStartCountIndex;
	}

	FILE *outXFormFit = fopen( "_transform_fit.txt", "a+" );

	//
	// Here we create vectors of points in 3D in the pre-operative model.
	//
	//
	vector< float > vecPoints3DImage1;
	vector< float > vecDotProduct;
	vector< float > vecPoints3DImage2;
	vector< float > vecPoints3DError;
	
	// Keep cooridinate errors
	vector< float > vecErrorX;
	vector< float > vecErrorY;
	vector< float > vecErrorZ;

	for( int i = 0; i < dat.iPoints; i++ )
	{
		if( piInlierCountPerFeature[i] <= 3 )
		{
			continue;
		}

		vecPoints3DImage1.clear();
		vecPoints3DImage2.clear();
		vecPoints3DError.clear();
		vecDotProduct.clear();
		vecErrorX.clear();
		vecErrorY.clear();
		vecErrorZ.clear();

		// Valid plane here, check error

		int iStartPointIndex = vecPoints3DImage2.size();
		int iInlierCount = piInlierCountPerFeature[i];



		for( int ik = i; ik < i + iInlierCount; ik++ )
		{
			float *point_in = dat.pfPoints1 + 15*ik;
			float *point_out = dat.pfPoints2 + 15*ik;

			// Copy in scaled points
			float p_in[3];
			float p_out[3];
			for( int j = 0; j < 3; j++ )
			{
				p_in[j] = fScale*point_in[12+j];
				p_out[j] = fScale*point_out[12+j];
			}

			// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
			memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
			memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

			// Point 1: multiply current calibration matrix
			float point_test_in[3];
			mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
			mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

			// Pass Point 2 through current calibration matrix
			float point_test_out[3];
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

			// Compute error: for debugging
			float fErrorSqr = sqrt( vec3D_distsqr_3d( point_test_in, point_test_out ) );

			vecErrorX.push_back( point_test_out[0]-point_test_in[0] );
			vecErrorY.push_back( point_test_out[1]-point_test_in[1] );
			vecErrorZ.push_back( point_test_out[2]-point_test_in[2] );

			vecDotProduct.push_back( pfDotProductZ[ik] );
			vecPoints3DError.push_back( fErrorSqr );
			for( int j = 0; j < 3; j++ )
			{
				vecPoints3DImage1.push_back( p_in[j] );				// In-plane pixel coordinates
				vecPoints3DImage2.push_back( point_test_out[j] );	// 3D coordinates
			}
			
		}

		int piInlierFlags[200];
		memset( piInlierFlags, 0, sizeof(piInlierFlags) );

		float pfPlane[4];
		// Estimate plane & error
		int iInliers;
		
		//iInliers = 
		//ransac_plane_estimation(
		//	iInlierCount,
		//	&vecPoints3DImage2[iStartPointIndex],
		//	100,
		//	1.5f,
		//	pfPlane );

		iInliers =
		determine_similarity_transform_ransac(			
			&vecPoints3DImage1[iStartPointIndex],
			&vecPoints3DImage2[iStartPointIndex],
			0, 0, iInlierCount, 100, 3.0f, piInlierFlags );

		if( iInliers >= 3 )
		{
			static int iFrameCount = 0;

			for( int j = 0; j < iInlierCount; j++ )
			{
				if( piInlierFlags[j] )
				{
					fprintf( outXFormFit, "%s\t%d\t%f\t%f\n", strrchr( pcFileName, '\\')+1, iFrameCount, vecPoints3DError[j], vecDotProduct[j] );
				}
			}			

			iInliers =
			determine_similarity_transform_ransac(			
				&vecPoints3DImage1[iStartPointIndex],
				&vecPoints3DImage2[iStartPointIndex],
				0, 0, iInlierCount, 100, 3.0f, piInlierFlags );
			iFrameCount++;
		}

		// Skip past current inliers
		i += iInlierCount-1;
	}

	fclose( outXFormFit );

	return 1;

	// Now loop through all frames and compute planes
	for( int i = 0; i < vecFrameIndices.size(); i++ )
	{
		float pfPlane[4];
		int iFrameStartIndex = vecFrameIndices[i];
		ransac_plane_estimation(
			piInlierCountPerFeature[iFrameStartIndex],
			&vecPoints3DImage2[3*iFrameStartIndex],
			100,
			1.5f,
			pfPlane );
	}

	for( int i = 0; i < dat.iPoints; i++ )
	{
		float *point_in = dat.pfPoints1 + 15*i;
		float *point_out = dat.pfPoints2 + 15*i;

		// Copy in scaled points
		float p_in[3];
		float p_out[3];
		for( int j = 0; j < 3; j++ )
		{
			p_in[j] = fScale*point_in[12+j];
			p_out[j] = fScale*point_out[12+j];
			if( j < 2 )
			{
				// Add 10 to each coordinate to see the effect
				//p_in[j] += 50;
				//p_out[j] += 50;
			}
		}

		// Copy in first 3 rows (length 3x4 = 12) of US tracker transform matrices
		memcpy( &mat44_1[0], point_in, sizeof(float)*12 );
		memcpy( &mat44_2[0], point_out, sizeof(float)*12 );

		// Point 1: multiply current calibration matrix
		float point_test_in[3];
		mult_4x4_matrix( mat44_1, mat44, mat44_1_p );
		mult_4x4_vector_homogenous( mat44_1_p, p_in, point_test_in );

		// Point 2: multiply current calibration matrix
		float point_test_out[3];
		if( dat.iErrorCalib )
		{
			// Pass Point 1 through ground truth calibration, for evaluating calibratin error
			mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_in, point_test_out );

			// Alternative 1 - pass point 2 though ground truth calibration
			mult_4x4_matrix( mat44_2, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );

			// Alternative 2 - for shifted image points, reverse shift,
			// pass through ground truth. 
			float point_out_shifted[3];
			point_out_shifted[0] = p_in[0] + 315;
			point_out_shifted[1] = p_in[1] + 315;
			point_out_shifted[2] = p_in[2];
			mult_4x4_matrix( mat44_1, &dat.pfCalib[0], mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, point_out_shifted, point_test_out );
		}
		else
		{
			// Pass Point 2 through current calibration matrix, for autocalibration optimization
			mult_4x4_matrix( mat44_2, mat44, mat44_2_p );
			mult_4x4_vector_homogenous( mat44_2_p, p_out, point_test_out );
		}

		// Rescale to pixel units??
		//vec3D_mult_scalar( point_test_in, 1.0f/fScale );
		//vec3D_mult_scalar( point_test_out, 1.0f/fScale );

		// Compute squared error
		float fErrorSqr = vec3D_distsqr_3d( point_test_in, point_test_out );

		for( int k = 0; k < 3; k++ )
			pdAbsErrorXYZ[k] = fabs(point_test_in[k] - point_test_out[k]);

		// Sum up 
		dSumSqrError += fErrorSqr;
		pfDistSqr[i] = fErrorSqr;
	}

	FILE *outf = fopen( "inliers.txt", "a+" );
	//fprintf( outf, "InlierCountInImage\tErrorThisFeature\tDotProductZ\tInlierFrequency (not related to individual features)\n" );
	for( int i = 0; i < dat.iPoints; i++ )
	{
		if( piInlierCountPerFeature[i] > 3 )
		{
			fprintf( outf, "%s\t%d\t%f\t%f\t%d\n", strrchr( pcFileName, '\\'), piInlierCountPerFeature[i], sqrt(pfDistSqr[i]), pfDotProductZ[i], piInlierFrequency[i] );
		}
	}
	fclose( outf );

	delete [] pfDotProductZ;
	delete [] piInlierCountPerFeature;
	delete [] piInlierFrequency;

	// Robust: keep top half of distance 
	qsort( pfDistSqr, dat.iPoints, sizeof(float), _compareFloat );
	dSumSqrError = 0;
	//float fPercentToKeep = 3.0/4.0;
	float fPercentToKeep = dat.fPercentToKeep;
	for( int i = 0; i < fPercentToKeep*dat.iPoints; i++ )
	{
		dSumSqrError += sqrt(pfDistSqr[i]);
	}

	delete [] pfDistSqr;

	return dSumSqrError / (fPercentToKeep*dat.iPoints);
}




//
// main()
//
// Main function for Autocalibration, TMI 2017
//
int
main(
	 int argc,
	 char **argv
	 )
{
	
	// Testing calibration
	// Test on rigid US minimization with simulated data

	my_constraint_data dat;

	if( read_point_data_US_multiple( argc-1, argv+1, dat ) < 0 )
	{
		printf( "Could not read data\n" );
		return -1;
	}

	float fError;
	FILE *outfile = fopen( "result.txt", "a+" );
	fprintf( outfile, "\n%s\t", argv[1] );
	for( int i = 0; i < 10/*10*/; i++ )
	{
		main_test_rigid_US_minimization( dat, fError, 1 );
		fprintf( outfile, "%f\t", fError );
	}
	fclose( outfile );

	return 1;
}


