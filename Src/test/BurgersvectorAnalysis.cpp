/*
 * BurgersvectorAnalysis.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */

#include <BurgersvectorAnalysis.h>


double find_BV_plane(double *bv, int *plane, int *bvl)
{
	//cout << "came here\n";
	double *vv;
	double result=2;
	double max_val=0.;
	int max_val_r=-1;
	*plane = -1000;
	*bvl = -2;
	vv   = (double *)malloc(3*sizeof(double));
	bool entered = false;
	for (int i=0;i< (VECTORS_ENUM);i++)
	{
		vv[0] = BVS[i][0];      vv[1] = BVS[i][1];      vv[2] = BVS[i][2];
		result = V3dot (vv, bv);

		// cout << fabs(result)<<"\n";
		if((fabs(max_val)< fabs(result))&&(fabs(result)>ANGLE_DEVIATION))
		{
			max_val = result;
			*bvl = i;
			if(!entered)
			{
				entered = true;
				//	cout << "first enter\t"<<i<<"\t"<<max_val<<"\n";
			}else
			{
				cout << "entering multiple times\t"<<i<<"\t"<<max_val<<"\n";
			}
		}



	}
	free(vv);
	return max_val;
	///if(max_val<0) *bvl = (*bvl)*-1;
	//  cout << BVS[*bvl][0]<<"**************\t"<<BVS[*bvl][1]<<"*****************\t"<<BVS[*bvl][2]<<"***********\t"<<max_val<<"*******************\n";

}


void determine_BV_plane(double rx,double ry,double rz, int *plane, double *angle)
{
	double *v,*bv;
	int pl=-2;
	int bvl=-100;
	v   = (double *)malloc(3*sizeof(double));
	bv   = (double *)malloc(3*sizeof(double));
	//  for (int id=0; id< 3; id++) v[id] = malloc(3*sizeof(double));
	v[0] = rx;
	v[1] = ry;
	v[2] = rz;
	V3norm(v);

	V3mulM3(v,H0_geo,bv);
	// cout << index<<"\n";
	*angle = find_BV_plane(bv,&pl,&bvl);
	*plane = bvl+1;
	// atom[index].pl = pl;
	// atom[index].tot = -1;

	//if(index==147545)
    //{

		double bv_max = max(fabs(bv[0]),fabs(bv[1]));
		bv_max = max(fabs(bv[2]),bv_max);
		// cout << bvl<<"\t"<<bv[0]<<"\t"<<bv[1]<<"\t"<<bv[2]<<"\t";
		// cout << atom[index].disrigistry[2]<<"\t"<<acos(fabs(atom[index].disrigistry[2]))*180/PI<<"\n";
		//}

		free(v);
		free(bv);

}

void determine_BV_plane(int index, int num_slip,double rx,double ry,double rz)
{
	double *v,*bv;
	int pl=-2;
	int bvl=-100;
	v   = (double *)malloc(3*sizeof(double));
	bv   = (double *)malloc(3*sizeof(double));
	//  for (int id=0; id< 3; id++) v[id] = malloc(3*sizeof(double));
	v[0] = rx;
	v[1] = ry;
	v[2] = rz;
	V3norm(v);

	V3mulM3(v,H0_geo,bv);
	// cout << index<<"\n";
	atom[index].disrigistry[2] = find_BV_plane(bv,&pl,&bvl);
	atom[index].BV = bvl+1;
	// atom[index].pl = pl;
	// atom[index].tot = -1;

	//if(index==147545)
    //{

		double bv_max = max(fabs(bv[0]),fabs(bv[1]));
		bv_max = max(fabs(bv[2]),bv_max);
		// cout << bvl<<"\t"<<bv[0]<<"\t"<<bv[1]<<"\t"<<bv[2]<<"\t";
		// cout << atom[index].disrigistry[2]<<"\t"<<acos(fabs(atom[index].disrigistry[2]))*180/PI<<"\n";
		//}

		free(v);
		free(bv);

}
