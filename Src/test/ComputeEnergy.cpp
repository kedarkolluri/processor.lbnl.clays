/*
 * ComputeEnergy.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: kedar
 */
#include <ComputeEnergy.h>
double compute_energy_difference( atomic_dat *atom_ref1, atomic_dat *atom_now1, int n_now, double lower_limit, bool replace, int *number_of_atoms)
{
	int rejected=0; int accepted=0;
	double energy_diff=0.0;
	for(int i=0;i<n_now;i++)
	{
		if(fabs(atom_now1[i].pe-atom_ref1[i].pe) > lower_limit)
		{
			energy_diff +=atom_now1[i].pe-atom_ref1[i].pe;
			accepted++;
		}else
		{
			rejected++;
		}

		//if(replace) atom_now1[i].pe = fabs(atom_now1[i].pe-atom_ref1[i].pe);
		if(replace) atom_now1[i].pe = (atom_now1[i].pe-atom_ref1[i].pe);

	}

	*number_of_atoms = accepted;
	return energy_diff;
}
double compute_energy( atomic_dat *atom_now, int n_now)
{
	double energy=0.0;
	for(int i=0;i<n_now;i++)
	{
		energy+=atom_now[i].pe;
	}
	return energy;
}
