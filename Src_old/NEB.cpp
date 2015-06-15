/*
 * NEB.cpp
 *
 *  Created on: Sep 14, 2009
 *      Author: kedar
 */

/*
Algorithm:

1. Compute the energies of all the states, and mark the states where you would like the method to
   be climbing image instead of springs (after few iterations of step 2)
2. take combination of 3 states with the middle state being the state for which properties are being
	calculated
	(a) Determine the tangent for these states - traditional (1998 paper) or hybrid (Mukul's suggestion) methods
	(b) compute force components as required and f_neb, if f_neb < tolerence, go to (e)
	(c) use velocity verlet to change the positions - store the velocities
	(d) compute the energies and forces for these states - go to (a)

*/

#include <NEB.h>


int NEB_SAVE_INFO = 3;
bool NEB_CLIMB = true;
int NEB_CLIMB_AFTER = 100;
int NEB_CLIMB_TOL = 1e-2;
int NEB_CLIMB_CONDITION = 1;
int CLIMB_TOL = 1e-3;
int NEB_MAX_ITER = 999999999;
double NEB_FORCE_TOL = 1e-5;
bool NEB_MODIFY_LINE_SEARCH = true;
bool NEB_SAVE_EACHSTEP = false;
double NEB_dt = 0.02;
double NEB_springK = 1;
double NEB_DELTOL = 1e-8;
bool NEB_BACKTRACK_M = false;
bool NEB_BACKTRACKFIRST = false;
double NEB_dt_FACTOR = 0.7;
double NEB_dt_FACTOR_DEFAULT = 0.7;
double NEB_BACKTRACK_ETOL = 0.051;
double NEB_BACKTRACK_dt_MIN = 6e-3;

int NEB_SLOT = 1;
//char NEB_SCRIPT="NEB.script"

bool *climb_states;
double **v;

double *energies;
double *neb_forces;
double *spring_forces;
double *true_forces;
double *reac_coord;
double *t_mag_all;
double *t_l_all;
double *t_r_all;
double max_move = 0.1;
double max_move_max = 0.2;
double dt_max = 0.08;
bool *isdynamic;
double *max_neb_forces;


int n_atoms=0;

int total_length=0;
double force_to_velocity = 1 / 1.0364269e-4;
double R_max = 0;

bool accept_change(atomic_dat *atom_d, int dont_first, int dont_second)
{

	bool return_val= false;
	int pp = 0;
	//if(dont_first==0) pp = dont_second+10000;
	save_lammps_specific("dat_backtrack",pp,atom_d, Hcry);
	char exec_[80]="";
	strcat(exec_,"./backtrack.script ");
	strcat(exec_,"dat_backtrack.");
	char str[80];
	sprintf(str,"%d",pp);
	strcat(exec_,str);
	if(UTILITIES_ZIP)
	{
		strcat(exec_,".gz");
	}


	execute_system_command(exec_);
	//cout << exec_<<"\n";

	int load_data = read_lammps("backtrack_dat.0",atom,true,true);
	double delE=0.0;
	for(int i =0;i<n;i++)
	{
		delE += atom[i].pe-atom_d[i].pe;
	}
	if((delE-NEB_BACKTRACK_ETOL)<=0) return_val = true;
	cout <<delE<<"\t"<<delE-NEB_BACKTRACK_ETOL<<"\t";
	free(atom); atom = NULL;

	return return_val;
}

void update_forces(int *lammps_filename_array,int chain_length, int iter)
{
	if(NEB_SAVE_EACHSTEP)
	{
		char exec_2[80] = "mkdir ";
		char str2[80];
		sprintf(str2,"%d",iter);
		strcat(exec_2,str2);
		execute_system_command(exec_2);
		char exec_3[80] = "cp dat_* ";
		char str3[80];
		sprintf(str3,"%d",iter);
		strcat(exec_3,str3);
		strcat(exec_3,"/");
		execute_system_command(exec_3);
	}

	if(NEB_SLOT==1)
	{

		for(int i =1;i<chain_length-1;i++)
		{
			char exec_[80]="";
			strcat(exec_,"./det_forces.lammps.script ");
			char str[80];
			sprintf(str,"%d",lammps_filename_array[i]);
			strcat(exec_,str);
			execute_system_command(exec_);
	//		cout << exec_<<"\n";

		}
	}else
	{

		if(((total_length-2)%NEB_SLOT)!=0)
		{
			cout << "not sufficient number of processors; not multiples;\n";
			exit(1);
		}

		for(int i = 1;i<chain_length-1;i=i+NEB_SLOT)
		{
			int one = lammps_filename_array[i];
			int last = lammps_filename_array[i+NEB_SLOT-1];
			char exec_[80]="";
			strcat(exec_,"./NEB.script ");
			char str[80];
			sprintf(str,"%d",one);
			strcat(exec_,str);
			strcat(exec_," ");
			sprintf(str,"%d",last);
			strcat(exec_,str);
			strcat(exec_," 1");
			cout << exec_<<"\n";

			execute_system_command(exec_);

		}
	}
}


void save_energies(int iter, int *filename_array, bool all)
{

	FILE *fptr;
	atomic_dat *atom_prev;
	double initial = energies[0];
	if(all)
	{

		double f_mag_1[total_length];

		for(int i=0;i<total_length;i++)
		{
			char filename[80]="dat.",str[80];
			sprintf(str,"%d",filename_array[i]);
			strcat(filename,str);
			int is_read = read_lammps(filename,atom,true,true);
			if(i == 0)
			{
				atom_prev = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				copy_atomstruct(atom, atom_prev,n);
			}

			reac_coord[i] =0.0;
			f_mag_1[i] = 0.0;
			for(int ii=0;ii<n;ii++)
			{
				double r[3],s[3];
				double delx = (atom[ii].sx-atom_prev[ii].sx);
				double dely = (atom[ii].sy-atom_prev[ii].sy);
				double delz = (atom[ii].sz-atom_prev[ii].sz);
				delx = delx-(int)(delx*2);
				dely = dely-(int)(dely*2);
				delz = delz-(int)(delz*2);
				s[0] = delx;s[1] =dely;s[2] = delz;
				V3mulM3(s,Hcry,r);
				reac_coord[i] +=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
				double ting = sqrt(atom[ii].fx*atom[ii].fx+atom[ii].fy*atom[ii].fy+atom[ii].fz*atom[ii].fz);
				if(ting > f_mag_1[i]) f_mag_1[i]= ting;

			}

			reac_coord[i] = sqrt(reac_coord[i]);

			if(i!=0)
			{
				reac_coord[i]+=reac_coord[i-1];
				copy_atomstruct(atom, atom_prev,n);
			}

		}
		free(atom_prev);
		char fn[80]="energy.",str[80];
		sprintf(str,"%d",iter);
					strcat(fn,str);
		fptr = fopen(fn,"w");
		fprintf(fptr, "\n");
		double R1_max = reac_coord[total_length-1]/R_max;
		for(int i =0;i< total_length;i++)
		{
			fprintf(fptr, "%d %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf\n",i,energies[i]-initial,neb_forces[i],true_forces[i],spring_forces[i],t_mag_all[i],t_l_all[i],t_r_all[i],iter, reac_coord[i]/reac_coord[total_length-1],reac_coord[i], reac_coord[total_length-1],sqrt(max_neb_forces[i]),f_mag_1[i]);
		}
		fclose(fptr);
	}
	char fn1[80]="neb_energy.info";
	fptr = fopen(fn1,"w");
	for(int i =0;i< total_length;i++)
	{
		fprintf(fptr, "%d %lf %lf %lf %lf %d %lf\n",i,energies[i]-initial,neb_forces[i],true_forces[i],spring_forces[i],iter,sqrt(max_neb_forces[i]));
	}
	fclose(fptr);


}

void determine_climbstates(int *filename_array, int chain_length, bool dummy_to_fill_energies)
{
	int max_climb_states = 0;

if(dummy_to_fill_energies)
{
	atomic_dat *atom_prev;

	char filename[80]="dat.",str_[80];
	sprintf(str_,"%d",filename_array[0]);
	strcat(filename,str_);
	int is_read = read_lammps(filename,atom,true,true);
	atom_prev = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));

	copy_atomstruct(atom, atom_prev,n);
	char filename1[80]="dat.",str1[80];
	sprintf(str1,"%d",filename_array[total_length-1]);
	strcat(filename1,str1);
	is_read = read_lammps(filename1,atom,true,true);
	if(is_read==0) cout << "read the second file\n";
	for(int ii=0;ii<n;ii++)
	{
		double r[3],s[3];
		double delx = (atom[ii].sx-atom_prev[ii].sx);
		double dely = (atom[ii].sy-atom_prev[ii].sy);
		double delz = (atom[ii].sz-atom_prev[ii].sz);
		delx = delx-(int)(delx*2);
		dely = dely-(int)(dely*2);
		delz = delz-(int)(delz*2);
		s[0] = delx;s[1] =dely;s[2] = delz;
		V3mulM3(s,Hcry,r);
		R_max +=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
		//cout << R_max<<"\n";

	}

	R_max = sqrt(R_max);
	cout << "R_max value is \t"<<R_max<<"\n";
	free(atom_prev);

	for(int i=0;i<chain_length;i++)
		{
			char filename[80]="dat.",str[80];
			sprintf(str,"%d",filename_array[i]);
			strcat(filename,str);

			int is_read = read_lammps(filename,atom,true,true);
			cout << "read lammps in climbstates\n";
			energies[i] = compute_energy(atom, n);
			cout << "finished computing energy\t"<<i<<"\t"<<energies[i]<<"\n";
			climb_states[i] = false;
		}
}else
{
		int max_energy_state = -1;
		double max_energy = -1e31;
		if(NEB_CLIMB_CONDITION ==1)
		{
			for(int i = 0;i<total_length;i++)
			{

					if(max_energy < energies[i])
					{
						max_energy = energies[i];
						max_energy_state = i;
						cout << max_energy_state<<"\n";
					}
			}
			if((max_energy_state>-1)&&(max_energy_state<total_length))
			{
				climb_states[max_energy_state] = true;
				cout << max_energy_state<< "\t set to max value\n";
				max_climb_states++;
			}
		}else if(NEB_CLIMB_CONDITION == 2)
		{
			for(int i=1;i<total_length-1;i++)
			{
				if((energies[i]> (CLIMB_TOL+ energies[i+1]))&&(energies[i]> (CLIMB_TOL+energies[i-1])))
				{
					climb_states[i] = true;
					max_climb_states++;
					climb_states[i-1] = climb_states[i+1] = false;
				}else if((energies[i]> (energies[i+1]))&&(energies[i]> (CLIMB_TOL+energies[i-1])))
				{
					climb_states[i] = true;
					max_climb_states++;
					climb_states[i-1]= false;
				//}
				//else if((energies[i]> (energies[i+1]+CLIMB_TOL))&&(energies[i]> (energies[i-1])))
				//{
				//	climb_states[i] = true;
				//	max_climb_states++;
				}else if((energies[i]== energies[i-1])&&(energies[i]> (energies[i+1]+CLIMB_TOL)))
				{
					climb_states[i] = true;
					max_climb_states++;
				}
			}
		}


 }

 cout << "\t"<< total_length<<"\t"<<n_atoms<<"\t"<< max_climb_states<<"\n";
}

int determine_tangent( atomic_dat * states[3],double en[3],double **t, double **t_l, double **t_r, double *t_mag, double *t_l_mag, double *t_r_mag, bool hybrid)
{
	int category = 0;
	double emax = max(fabs(en[0]-en[1]),fabs(en[1]-en[2]));
	double emin = min(fabs(en[0]-en[1]),fabs(en[1]-en[2]));
	double r[3],s[3];

	for(int i =0;i<n_atoms;i++)
	{
		double delx = ((states[2])[i].sx-(states[1])[i].sx);
		double dely = ((states[2])[i].sy-(states[1])[i].sy);
		double delz = ((states[2])[i].sz-(states[1])[i].sz);
		delx = delx-(int)(delx*2);
		dely = dely-(int)(dely*2);
		delz = delz-(int)(delz*2);
		s[0] = delx;s[1] =dely;s[2] = delz;
		V3mulM3(s,Hcry,r);
		t_r[i][0] = r[0];t_r[i][1]=r[1];t_r[i][2]=r[2];

		delx = ((states[1])[i].sx-(states[0])[i].sx);
		dely = ((states[1])[i].sy-(states[0])[i].sy);
		delz = ((states[1])[i].sz-(states[0])[i].sz);
		delx = delx-(int)(delx*2);
		dely = dely-(int)(dely*2);
		delz = delz-(int)(delz*2);

		s[0] = delx;s[1] =dely;s[2] = delz;
		V3mulM3(s,Hcry,r);
		t_l[i][0] = r[0];t_l[i][1]=r[1];t_l[i][2]=r[2];
		for(int j=0;j<3;j++)
		{
			*t_r_mag += t_r[i][j]*t_r[i][j];
			*t_l_mag += t_l[i][j]*t_l[i][j];
		}
	}
	*t_r_mag = sqrt(*t_r_mag);
	*t_l_mag = sqrt(*t_l_mag);


if(hybrid)
{
	if((en[0]<en[1])&&(en[1]<en[2]))
	{
		category = 1;

		for(int i =0;i<n_atoms;i++)
		{
			for(int j=0;j<3;j++)
			{
				t[i][j]=t_r[i][j];
			}
		}

		*t_mag = *t_r_mag;


	}else if((en[0]>en[1])&&(en[1]>en[2]))
	{
		category = 2;
		for(int i =0;i<n_atoms;i++)
		{
			for(int j=0;j<3;j++)
			{
				t[i][j]=t_l[i][j];
			}
		}
		*t_mag = *t_l_mag;
	}else
	{
		double e1=emax; double e2=emin;
		if(en[2]>en[0])
		{
			category = 3;
		}else
		{
			e1=emin; e2=emax;
			category = 4;
		}

		for(int i =0;i<n_atoms;i++)
		{

			for(int j=0;j<3;j++)
			{
				t[i][j]=e2*t_l[i][j]+e1*t_r[i][j];
				*t_mag += t[i][j]*t[i][j];
			}
		}

		*t_mag = sqrt(*t_mag);

	}
}else
{
	category = 0;

			for(int i =0;i<n_atoms;i++)
			{
				double delx = ((states[2])[i].sx-(states[0])[i].sx);
				double dely = ((states[2])[i].sy-(states[0])[i].sy);
				double delz = ((states[2])[i].sz-(states[0])[i].sz);
				delx = delx-(int)(delx*2);
				dely = dely-(int)(dely*2);
				delz = delz-(int)(delz*2);
				s[0] = delx;s[1] =dely;s[2] = delz;
				V3mulM3(s,Hcry,r);
				t[i][0] = r[0];t[i][1]=r[1];t[i][2]=r[2];
				for(int j=0;j<3;j++)
				{
					*t_mag += t[i][j]*t[i][j];
				}
			}

			*t_mag = sqrt(*t_mag);


}
	return category;
}

void initiate_others()
{
	energies = (double *) malloc(total_length*sizeof(double));
	reac_coord = (double *)malloc(total_length*sizeof(double));
	neb_forces = (double *)malloc(total_length*sizeof(double));
	max_neb_forces = (double *)malloc(total_length*sizeof(double));
	true_forces = (double *)malloc(total_length*sizeof(double));
	spring_forces = (double *)malloc(total_length*sizeof(double));
	climb_states = (bool *) malloc(total_length*sizeof(bool));

	t_mag_all = (double *)malloc(total_length*sizeof(double));
	t_r_all =(double *)malloc(total_length*sizeof(double));
	t_l_all = (double *)malloc(total_length*sizeof(double));

	for(int i=0;i<total_length;i++)
	{
		energies[i] = reac_coord[i] = neb_forces[i] = spring_forces[i] = true_forces[i] = 0.0;
		t_mag_all[i] = t_r_all[i] = t_l_all[i] = 0.0;
		max_neb_forces[i] = 0.0;
	}

}


int NEB3(int *filename_array,int chain_length, double spring_K, double dt)
{

	total_length = chain_length;
	initiate_others();
	FILE *fptr_convergence,*fptr_force_magnitudes;
	char convergence_filename[80] = "neb_convergence.info";
	fptr = fopen (convergence_filename,"w");
	fclose(fptr);

	char forcemag_filename[80] = "neb_force_magnitdues.info";
	fptr = fopen (forcemag_filename,"w");


	total_length = chain_length;

	for(int i =0;i<chain_length;i++)
	{
		energies[i]=0.0;
		reac_coord[i]=0.0;
		climb_states[i]=false;
	}
	// obtain number of atoms
	char filename[80]="dat.",str[80];
	sprintf(str,"%d",filename_array[0]);
	strcat(filename,str);
	int is_read = read_lammps(filename,atom,true,true);
	cout << is_read<<"\n";
	if(is_read==0)
		{
			n_atoms = n;

		}
	cout << n_atoms<<"\t"<<n<<"\n";
	cout << "before freeing atom\n";
	free(atom);atom = NULL;
	//end obtain number of atoms

	//initialize velocities array
	cout << "before setting velocites\n";
	v = (double **) malloc(n*chain_length*sizeof(double));
	for(int i=0;i<n*chain_length;i++)
	{
		v[i] = (double *) malloc(3*sizeof(double));
	}
	cout << "set velocities\n";
	for(int i =0; i<chain_length;i++)
	{
		for(int j=0;j<n;j++)
		{
			for(int k=0;k<3;k++)
			{
				v[i*n+j][k] = 0.0;
			}
		}
	}

	//end initialize velocities array
	// fill the energy array
	determine_climbstates(filename_array, chain_length, true);
	cout << energies[0]<<"\t"<<energies[1]<<"\t"<<energies[2]<<"\n";

	// End of all initializations

	// start take groups of 3 and determine energies, forces and steps
	double ftol = 1.0e-5;
	double fneb_val = 0.0;
	double fneb_prev = fneb_val;
	int iteration=0;
	int category=-1;
//	double spring_K = 0.001;
	do
	{
		cout << iteration<<"\t"<<fneb_prev<<"\t"<<fneb_val<<"\t"<< category<<" he he he he he \n";
		fptr_convergence = fopen(convergence_filename,"a");
		fprintf(fptr_convergence,"%d %lf %lf %12.9lf\n",iteration, fneb_prev,fneb_val, fabs(fneb_prev-fneb_val));
		fclose(fptr_convergence);
		fneb_prev = fneb_val;

		fneb_val = 0.0;
		if(iteration%3 == 0)
		{
			save_energies(iteration, filename_array,true);
		}else
		{
			save_energies(iteration, filename_array,false);
		}
		update_forces(filename_array, chain_length, iteration);

		for(int i =1;i<chain_length-1;i++)
		{
			double **t,**t_l,**t_r, t_mag,t_l_mag,t_r_mag;
			t = (double **)malloc(n_atoms*sizeof(double));
			t_l = (double **)malloc(n_atoms*sizeof(double));
			t_r = (double **)malloc(n_atoms*sizeof(double));

			for(int j=0;j<n_atoms;j++)
			{
				t[j] = (double *)malloc(n_atoms*sizeof(double));
				t_l[j] = (double *)malloc(n_atoms*sizeof(double));
				t_r[j] = (double *)malloc(n_atoms*sizeof(double));
				for(int jj=0;jj<3;jj++)
				{
					t[j][jj]=t_l[j][jj]=t_r[j][jj]=0.0;
				}
			}
			t_mag = t_l_mag = t_r_mag = 0.0;

			atomic_dat *neb_state[3];
			double en[3];

			en[0]=en[1]=en[2]=0.0;

// load sets of 3 states
			for(int ii=-1;ii<=1;ii++)
			{
				char filename[80]="dat.",str[80];
				sprintf(str,"%d",filename_array[i+ii]);
				strcat(filename,str);

				is_read = read_lammps(filename,atom,true,true);
				neb_state[ii+1] = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				en[ii+1] = compute_energy(atom, n);
				energies[i+ii] = en[ii+1];
				copy_atomstruct(atom, neb_state[ii+1],n);
			}

			category = determine_tangent(neb_state,en,t,t_l,t_r,&t_mag,&t_l_mag,&t_r_mag, true);
//u is the new f_neb in the atomic structure
			cout << category<< " here the category is \n";
			double ftrue_mag = 0; double fspring_mag = 0;double fneb_mag = 0;
			for(int j=0;j<n_atoms;j++)
			{
				//F_trueforce_component_perpendicular_to_tangent

				(neb_state[1])[j].ux = (neb_state[1])[j].fx*(1-t[j][0]*t[j][0]/t_mag/t_mag);
				(neb_state[1])[j].uy = (neb_state[1])[j].fy*(1-t[j][0]*t[j][1]/t_mag/t_mag);
				(neb_state[1])[j].uz = (neb_state[1])[j].fz*(1-t[j][0]*t[j][2]/t_mag/t_mag);
				ftrue_mag += ((neb_state[1])[j].ux)*((neb_state[1])[j].ux)+((neb_state[1])[j].uy)*((neb_state[1])[j].uy)+((neb_state[1])[j].uz)*((neb_state[1])[j].uz);

				// add F_spring
				(neb_state[1])[j].ux += spring_K*(t_r_mag-t_l_mag)*t[j][0]/t_mag;
				(neb_state[1])[j].uy += spring_K*(t_r_mag-t_l_mag)*t[j][1]/t_mag;
				(neb_state[1])[j].uz += spring_K*(t_r_mag-t_l_mag)*t[j][2]/t_mag;

				fneb_mag += ((neb_state[1])[j].ux)*((neb_state[1])[j].ux)+((neb_state[1])[j].uy)*((neb_state[1])[j].uy)+((neb_state[1])[j].uz)*((neb_state[1])[j].uz);

				(neb_state[1])[j].delr[0] = (neb_state[1])[j].ux*(neb_state[1])[j].ux+(neb_state[1])[j].uy*(neb_state[1])[j].uy +(neb_state[1])[j].uz*(neb_state[1])[j].uz;
			}
			fneb_mag = sqrt(fneb_mag);
			if(fneb_mag > fneb_val) fneb_val = fneb_mag;
			double force_to_velocity = 1 / 1.0364269e-4;
//			double dt=0.01;
			int zero_value=0;
			for(int j=0;j<n_atoms;j++)
			{
				//update velocities - first half step

				v[i*n+j][0] += neb_state[1][j].ux*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][1] += neb_state[1][j].uy*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][2] += neb_state[1][j].uz*force_to_velocity*dt*0.5/neb_state[1][j].ma;

				// velocity-verlet change in coordinates

				double r[3],s[3];
				s[0] = neb_state[1][j].sx;s[1] =neb_state[1][j].sy;s[2] = neb_state[1][j].sz;
				V3mulM3(s,Hcry,r);

				// change in coords due to velocity component in the direction of force and due to force itself
				double vel_component = 0;
				if(v[i*n+j][0]*neb_state[1][j].ux>0)
				{
					vel_component = v[i*n+j][0]*neb_state[1][j].ux/fneb_mag*neb_state[1][j].ux/fneb_mag;
					//vel_component = v[i*n+j][0]*neb_state[1][j].ux/fneb_mag*neb_state[1][j].ux*force_to_velocity/neb_state[1][j].ma;
				}
				r[0]+=vel_component*dt;
				r[0]+=dt*dt*0.5*neb_state[1][j].ux*force_to_velocity/neb_state[1][j].ma;

				if(v[i*n+j][1]*neb_state[1][j].uy>0)
				{
					vel_component = v[i*n+j][1]*neb_state[1][j].uy/fneb_mag*neb_state[1][j].uy/fneb_mag;
					//vel_component = v[i*n+j][0]*neb_state[1][j].uy/fneb_mag*neb_state[1][j].uy*force_to_velocity/neb_state[1][j].ma;

				}
				r[1]+=vel_component*dt;
				r[1]+=dt*dt*0.5*neb_state[1][j].uy*force_to_velocity/neb_state[1][j].ma;

				if(v[i*n+j][2]*neb_state[1][j].uz>0)
				{
					vel_component = v[i*n+j][2]*neb_state[1][j].uz/fneb_mag*neb_state[1][j].uz/fneb_mag;
					//vel_component = v[i*n+j][0]*neb_state[1][j].uz/fneb_mag*neb_state[1][j].uz*force_to_velocity/neb_state[1][j].ma;

				}
				r[2]+=vel_component*dt;
				r[1]+=dt*dt*0.5*neb_state[1][j].uz*force_to_velocity/neb_state[1][j].ma;


				V3mulM3(r,Hcry_inv,s);
				neb_state[1][j].sx = s[0];neb_state[1][j].sy = s[1];neb_state[1][j].sz = s[2];
				if(neb_state[1][j].sx>1.0) neb_state[1][j].sx = neb_state[1][j].sx-1;
				if(neb_state[1][j].sy>1.0) neb_state[1][j].sy = neb_state[1][j].sy-1;
				if(neb_state[1][j].sz>1.0) neb_state[1][j].sz = neb_state[1][j].sz-1;

				if(neb_state[1][j].sx<0.0) neb_state[1][j].sx = 1.0+neb_state[1][j].sx;
				if(neb_state[1][j].sy<0.0) neb_state[1][j].sy = 1.0+neb_state[1][j].sy;
				if(neb_state[1][j].sz<0.0) neb_state[1][j].sz = 1.0+neb_state[1][j].sz;

				//update velocities - second half step

				v[i*n+j][0] += neb_state[1][j].ux*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][1] += neb_state[1][j].uy*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][2] += neb_state[1][j].uz*force_to_velocity*dt*0.5/neb_state[1][j].ma;

						//v[i*n+j][0] = 0;v[i*n+j][1] = 0;v[i*n+j][2] = 0;


			}
			if(atom!=NULL){	free(atom); atom=NULL;}
			atom = neb_state[1];
			save_lammps(filename_array[i],Hcry);
//			update_forces(filename_array[i]);
			for(int j=0;j<n_atoms;j++)
			{
				free(t[j]);free(t_l[j]);free(t_r[j]);
			}

			free(t);free(t_l);free(t_r);
			free(neb_state[0]);free(neb_state[1]);free(neb_state[2]);
		}

		iteration++;
	}while (ftol< fneb_val);
}


int NEB2(int *filename_array,int chain_length, double spring_K, double dt)
{

	total_length = chain_length;
	initiate_others();
	FILE *fptr_convergence,*fptr_force_magnitudes;
	char convergence_filename[80] = "neb_convergence.info";
	fptr = fopen (convergence_filename,"w");
	fclose(fptr);

	char forcemag_filename[80] = "neb_force_magnitdues.info";
	fptr = fopen (forcemag_filename,"w");



	for(int i =0;i<chain_length;i++)
	{
		energies[i]=0.0;
		reac_coord[i]=0.0;
		climb_states[i]=false;
	}
	// obtain number of atoms
	char filename[80]="dat.",str[80];
	sprintf(str,"%d",filename_array[0]);
	strcat(filename,str);
	int is_read = read_lammps(filename,atom,true,true);
	cout << is_read<<"\n";
	if(is_read==0)
		{
			n_atoms = n;

		}
	cout << n_atoms<<"\t"<<n<<"\n";
	cout << "before freeing atom\n";
	free(atom);atom = NULL;
	//end obtain number of atoms

	//initialize velocities array
	cout << "before setting velocites\n";
	v = (double **) malloc(n*chain_length*sizeof(double));
	for(int i=0;i<n*chain_length;i++)
	{
		v[i] = (double *) malloc(3*sizeof(double));
	}
	cout << "set velocities\n";
	for(int i =0; i<chain_length;i++)
	{
		for(int j=0;j<n;j++)
		{
			for(int k=0;k<3;k++)
			{
				v[i*n+j][k] = 0.0;
			}
		}
	}

	//end initialize velocities array
	// fill the energy array
	determine_climbstates(filename_array, chain_length, true);
	cout << energies[0]<<"\t"<<energies[1]<<"\t"<<energies[2]<<"\n";

	// End of all initializations

	// start take groups of 3 and determine energies, forces and steps
	double ftol = 1.0e-5;
	double fneb_val = 0.0;
	double fneb_prev = fneb_val;
	int iteration=0;
	int category=-1;
//	double spring_K = 0.001;
	do
	{
		cout << iteration<<"\t"<<fneb_prev<<"\t"<<fneb_val<<"\t"<< category<<" he he he he he \n";
		fptr_convergence = fopen(convergence_filename,"a");
		fprintf(fptr_convergence,"%d %lf %lf %12.9lf\n",iteration, fneb_prev,fneb_val, fabs(fneb_prev-fneb_val));
		fclose(fptr_convergence);
		fneb_prev = fneb_val;

		fneb_val = 0.0;
		if(iteration%3 == 0)
		{
			save_energies(iteration, filename_array,true);
		}else
		{
			save_energies(iteration, filename_array, false);
		}
		update_forces(filename_array, chain_length,iteration);

		for(int i =1;i<chain_length-1;i++)
		{
			double **t,**t_l,**t_r, t_mag,t_l_mag,t_r_mag;
			t = (double **)malloc(n_atoms*sizeof(double));
			t_l = (double **)malloc(n_atoms*sizeof(double));
			t_r = (double **)malloc(n_atoms*sizeof(double));

			for(int j=0;j<n_atoms;j++)
			{
				t[j] = (double *)malloc(n_atoms*sizeof(double));
				t_l[j] = (double *)malloc(n_atoms*sizeof(double));
				t_r[j] = (double *)malloc(n_atoms*sizeof(double));
				for(int jj=0;jj<3;jj++)
				{
					t[j][jj]=t_l[j][jj]=t_r[j][jj]=0.0;
				}
			}
			t_mag = t_l_mag = t_r_mag = 0.0;

			atomic_dat *neb_state[3];
			double en[3];

			en[0]=en[1]=en[2]=0.0;

// load sets of 3 states
			for(int ii=-1;ii<=1;ii++)
			{
				char filename[80]="dat.",str[80];
				sprintf(str,"%d",filename_array[i+ii]);
				strcat(filename,str);

				is_read = read_lammps(filename,atom,true,true);
				neb_state[ii+1] = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				en[ii+1] = compute_energy(atom, n);
				energies[i+ii] = en[ii+1];
				copy_atomstruct(atom, neb_state[ii+1],n);
			}

			category = determine_tangent(neb_state,en,t,t_l,t_r,&t_mag,&t_l_mag,&t_r_mag, true);
//u is the new f_neb in the atomic structure
			//cout << category<< " here the category is \n";
			double ftrue_mag = 0; double fspring_mag = 0;double fneb_mag = 0;

			// Determine the dot product of force and tangent
			double fdott = 0.0;
			double vdotf = 0.0;
			double entire_force_magnitude = 0.0;

			for(int j=0;j<n_atoms;j++)
			{
				fdott+=(neb_state[1])[j].fx*t[j][0]/t_mag;
				fdott+=(neb_state[1])[j].fy*t[j][1]/t_mag;
				fdott+=(neb_state[1])[j].fz*t[j][2]/t_mag;
			}

			// determine force and velocity component in the direction of entire force
			// entire force on the image is stored in delr variable of the atom structure
			// magnitude of neb force is stored in delr[3] of the atom structure

			spring_forces[i] = 0.0;
			for(int j=0;j<n_atoms;j++)
			{
				//F_trueforce_component_perpendicular_to_tangent

				(neb_state[1])[j].ux = (neb_state[1])[j].fx-fdott*t[j][0]/t_mag;
				(neb_state[1])[j].uy = (neb_state[1])[j].fy-fdott*t[j][1]/t_mag;;
				(neb_state[1])[j].uz = (neb_state[1])[j].fz-fdott*t[j][2]/t_mag;;
				ftrue_mag += ((neb_state[1])[j].ux)*((neb_state[1])[j].ux)+((neb_state[1])[j].uy)*((neb_state[1])[j].uy)+((neb_state[1])[j].uz)*((neb_state[1])[j].uz);



				// add F_spring
				(neb_state[1])[j].ux += spring_K*(t_r_mag-t_l_mag)*t[j][0]/t_mag;
				(neb_state[1])[j].uy += spring_K*(t_r_mag-t_l_mag)*t[j][1]/t_mag;
				(neb_state[1])[j].uz += spring_K*(t_r_mag-t_l_mag)*t[j][2]/t_mag;

				fneb_mag += ((neb_state[1])[j].ux)*((neb_state[1])[j].ux)+((neb_state[1])[j].uy)*((neb_state[1])[j].uy)+((neb_state[1])[j].uz)*((neb_state[1])[j].uz);

				(neb_state[1])[j].delr[3] = (neb_state[1])[j].ux*(neb_state[1])[j].ux+(neb_state[1])[j].uy*(neb_state[1])[j].uy +(neb_state[1])[j].uz*(neb_state[1])[j].uz;

				//update velocities - first half step
				v[i*n+j][0] += neb_state[1][j].ux*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][1] += neb_state[1][j].uy*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][2] += neb_state[1][j].uz*force_to_velocity*dt*0.5/neb_state[1][j].ma;


/*
				//velocity component parallel to the entire force acting on the image
				neb_state[1][j].delr[0] = neb_state[1][j].fx+spring_K*(t_r[j][0]-t_l[j][0]);
				neb_state[1][j].delr[1] = neb_state[1][j].fy+spring_K*(t_r[j][1]-t_l[j][1]);
				neb_state[1][j].delr[2] = neb_state[1][j].fz+spring_K*(t_r[j][2]-t_l[j][2]);


				for(int jj=0;jj<3;jj++)
				{
					entire_force_magnitude += neb_state[1][j].delr[jj]*neb_state[1][j].delr[jj];
					vdotf += v[i*n+j][jj]*neb_state[1][j].delr[jj];
				}
*/

				//velocity component parallel to the NEB force acting on the image
				neb_state[1][j].delr[0] = neb_state[1][j].ux;
				neb_state[1][j].delr[1] = neb_state[1][j].uy;
				neb_state[1][j].delr[2] = neb_state[1][j].uz;
				for(int jj=0;jj<3;jj++)
				{
					entire_force_magnitude += neb_state[1][j].delr[jj]*neb_state[1][j].delr[jj];
					vdotf += v[i*n+j][jj]*neb_state[1][j].delr[jj];
				}

			}

			if(vdotf<0) vdotf = 0.0;

			entire_force_magnitude = sqrt(entire_force_magnitude);

			fneb_mag = sqrt(fneb_mag);

			neb_forces[i] = fneb_mag;

			if(fneb_mag > fneb_val) fneb_val = fneb_mag;

//			double dt=0.01;
			int zero_value=0;
			for(int j=0;j<n_atoms;j++)
			{




				// velocity-verlet change in coordinates

				double r[3],s[3];
				s[0] = neb_state[1][j].sx;s[1] =neb_state[1][j].sy;s[2] = neb_state[1][j].sz;
				V3mulM3(s,Hcry,r);
				if((j==3996)&&((i==9)||(i==10)))
				{
					cout << "***** before "<<i<<"\t"<<iteration<<" ******\n";
					cout << neb_state[1][j].ux<<neb_state[1][j].uy<<neb_state[1][j].uz<<"\n";
					cout << r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\n";
					cout << v[i*n+j][0]<<"\t"<<v[i*n+j][1]<<"\t"<<v[i*n+j][2]<<"\n";

				}
				// change in coords due to velocity component in the direction of force and due to force itself
				double vel_component = 0;
				double r_change[3];
				r_change[0] = r_change[1] = r_change[2] = 0;
				r_change[0] = dt*neb_state[1][j].ux;
				r_change[1] = dt*neb_state[1][j].uy;
				r_change[2] = dt*neb_state[1][j].uz;
				double r_change_mag = (r_change[0]*r_change[0]+r_change[1]*r_change[1]+r_change[2]*r_change[2]);
				if((r_change[0]*r_change[0]+r_change[1]*r_change[1]+r_change[2]*r_change[2]) > max_move*max_move)
				{
					for(int jj=0;jj<3;jj++)
					{
						r_change[jj] = r_change[jj]*max_move/r_change_mag;
					}
				}
				for(int jj=0;jj<3;jj++)
				{
					r[jj] += r_change[jj];
				}
/*
				// Velocity verlet based atom movement - Begin
				vel_component = v[i*n+j][0]*vdotf/entire_force_magnitude*neb_state[1][j].delr[0]/entire_force_magnitude;
			//	r[0]+=vel_component*dt;
				r[0]+=dt*dt*0.5*neb_state[1][j].ux*force_to_velocity/neb_state[1][j].ma;

				vel_component = v[i*n+j][1]*vdotf/entire_force_magnitude*neb_state[1][j].delr[1]/entire_force_magnitude;
			//	r[1]+=vel_component*dt;
				r[1]+=dt*dt*0.5*neb_state[1][j].uy*force_to_velocity/neb_state[1][j].ma;

				vel_component = v[i*n+j][2]*vdotf/entire_force_magnitude*neb_state[1][j].delr[2]/entire_force_magnitude;
			//	r[2]+=vel_component*dt;
				r[1]+=dt*dt*0.5*neb_state[1][j].uz*force_to_velocity/neb_state[1][j].ma;
				// Velocity verlet based atom movement - End
*/
				V3mulM3(r,Hcry_inv,s);
				neb_state[1][j].sx = s[0];neb_state[1][j].sy = s[1];neb_state[1][j].sz = s[2];
				if(neb_state[1][j].sx>1.0) neb_state[1][j].sx = neb_state[1][j].sx-1;
				if(neb_state[1][j].sy>1.0) neb_state[1][j].sy = neb_state[1][j].sy-1;
				if(neb_state[1][j].sz>1.0) neb_state[1][j].sz = neb_state[1][j].sz-1;

				if(neb_state[1][j].sx<0.0) neb_state[1][j].sx = 1.0+neb_state[1][j].sx;
				if(neb_state[1][j].sy<0.0) neb_state[1][j].sy = 1.0+neb_state[1][j].sy;
				if(neb_state[1][j].sz<0.0) neb_state[1][j].sz = 1.0+neb_state[1][j].sz;

				//update velocities - second half step

				v[i*n+j][0] += neb_state[1][j].ux*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][1] += neb_state[1][j].uy*force_to_velocity*dt*0.5/neb_state[1][j].ma;
				v[i*n+j][2] += neb_state[1][j].uz*force_to_velocity*dt*0.5/neb_state[1][j].ma;

				if((j==3996)&&((i==9)||(i==10)))
				{
					cout << "***** after "<<i<<"\t"<<iteration<<" ******\n";
					cout << r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\n";
					cout << v[i*n+j][0]<<"\t"<<v[i*n+j][1]<<"\t"<<v[i*n+j][2]<<"\n";

					cout << "***** end "<<i<<"\t"<<iteration<<" ******\n";
				}

						//v[i*n+j][0] = 0;v[i*n+j][1] = 0;v[i*n+j][2] = 0;


			}
			if(atom!=NULL){	free(atom); atom=NULL;}
			atom = neb_state[1];
			save_lammps(filename_array[i],Hcry);
//			update_forces(filename_array[i]);
			for(int j=0;j<n_atoms;j++)
			{
				free(t[j]);free(t_l[j]);free(t_r[j]);
			}

			free(t);free(t_l);free(t_r);
			free(neb_state[0]);free(neb_state[1]);free(neb_state[2]);
		}

		iteration++;
	}while (ftol< fneb_val);
}

int NEB(int *filename_array,int chain_length)
{

	if(NEB_dt_FACTOR >=1) NEB_dt_FACTOR = NEB_dt_FACTOR_DEFAULT ;

	total_length = chain_length;
	initiate_others();
	FILE *fptr_convergence,*fptr_force_magnitudes;
	char convergence_filename[80] = "neb_convergence.info";
	fptr = fopen (convergence_filename,"w");
	fclose(fptr);



	for(int i =0;i<chain_length;i++)
	{
		energies[i]=0.0;
		reac_coord[i]=0.0;
		climb_states[i]=false;
	}
	// obtain number of atoms
	char filename[80]="dat.",str[80];
	sprintf(str,"%d",filename_array[0]);
	strcat(filename,str);
	int is_read = read_lammps(filename,atom,true,true);
	cout << is_read<<"\n";
	if(is_read==0)
		{
			n_atoms = n;

		}
	cout << n_atoms<<"\t"<<n<<"\n";
	cout << "before freeing atom\n";
	free(atom);atom = NULL;
	//end obtain number of atoms

	// fill the energy, rc_distance array
	determine_climbstates(filename_array, chain_length, true);
	int *dynamic_list;
	dynamic_list = NULL;
	int single_len = 0;
	isdynamic = (bool *)malloc(n_atoms*sizeof(bool));
	for(int i=0;i<n_atoms;i++)
	{
		isdynamic[i] = true;
	}

	getdynamic(filename_array, chain_length, dynamic_list, &single_len, isdynamic);

	cout << energies[0]<<"\t"<<energies[1]<<"\t"<<energies[2]<<"\n";

	// End of all initializations

	// start take groups of 3 and determine energies, forces and steps
	double ftol = NEB_FORCE_TOL;
	double fneb_val = 0.0;
	double fneb_prev = fneb_val;
	int iteration=0;
	int category=-1;
//	double spring_K = 0.001;
	save_energies(iteration, filename_array, true); exit(1);
	double fneb_max=0.0;
	double fspring_max = 0.0;
	double ftrueforce_max = 0.0;
	do
	{

		if((NEB_CLIMB_AFTER>=0)&&(iteration==NEB_CLIMB_AFTER))
		{
			determine_climbstates(filename_array, chain_length, false);
		}

		cout << iteration<<"\t"<<fneb_prev<<"\t"<<fneb_val<<" \n";
		fptr_convergence = fopen(convergence_filename,"a");
		fprintf(fptr_convergence,"%d %lf %lf %12.9lf %lf %lf %lf %lf %lf \n",iteration, fneb_prev,fneb_val, fabs(fneb_prev-fneb_val),NEB_dt,max_move,sqrt(fneb_max),sqrt(fspring_max),sqrt(ftrueforce_max));
		fclose(fptr_convergence);
		fneb_max=0.0;
		fspring_max = 0.0;
		ftrueforce_max = 0.0;

		if(NEB_MODIFY_LINE_SEARCH)
		{

			if((fabs(fneb_prev-fneb_val) <3e-4) && (fneb_val > 5*ftol))
			{
				if (max_move < max_move_max)
				{
					max_move = min(max_move*1.5,max_move_max);
				}

				if(NEB_dt< dt_max)
				{
					NEB_dt= min(NEB_dt*2,dt_max);
				}
			}
		}

		fneb_prev = fneb_val;

		fneb_val = 0.0;
		if(iteration%NEB_SAVE_INFO == 0)
		{
			save_energies(iteration, filename_array, true);
		}else
		{
			save_energies(iteration, filename_array, false);
		}
		update_forces(filename_array, chain_length,iteration);
		for(int i =1;i<chain_length-1;i++)
		{
			double **t,**t_l,**t_r, t_mag,t_l_mag,t_r_mag;
			t = (double **)malloc(n_atoms*sizeof(double));
			t_l = (double **)malloc(n_atoms*sizeof(double));
			t_r = (double **)malloc(n_atoms*sizeof(double));

			for(int j=0;j<n_atoms;j++)
			{
				t[j] = (double *)malloc(n_atoms*sizeof(double));
				t_l[j] = (double *)malloc(n_atoms*sizeof(double));
				t_r[j] = (double *)malloc(n_atoms*sizeof(double));
				for(int jj=0;jj<3;jj++)
				{
					t[j][jj]=t_l[j][jj]=t_r[j][jj]=0.0;
				}
			}
			t_mag = t_l_mag = t_r_mag = 0.0;

			atomic_dat *neb_state[3];
			double en[3];

			en[0]=en[1]=en[2]=0.0;

// load sets of 3 states
			for(int ii=-1;ii<=1;ii++)
			{
				char filename[80]="dat.",str[80];
				sprintf(str,"%d",filename_array[i+ii]);
				strcat(filename,str);

				is_read = read_lammps(filename,atom,true,true);
				neb_state[ii+1] = (struct atomic_dat *) malloc((n+3)*sizeof(struct atomic_dat));
				en[ii+1] = compute_energy(atom, n);
				energies[i+ii] = en[ii+1];
				copy_atomstruct(atom, neb_state[ii+1],n);
			}

			category = determine_tangent(neb_state,en,t,t_l,t_r,&t_mag,&t_l_mag,&t_r_mag, true);
//u is the new f_neb in the atomic structure
			//cout << category<< " here the category is \n";

			double ftrue_mag = 0; double fspring_mag = 0;double fneb_mag = 0;
			t_mag_all[i]= t_mag;
			t_l_all[i] = t_l_mag;
			t_r_all[i] = t_r_mag;
			// Determine the dot product of force and tangent
			double fdott = 0.0;
			double vdotf = 0.0;
			double entire_force_magnitude = 0.0;

			for(int j=0;j<n_atoms;j++)
			{
				if(isdynamic[j])
				{
					fdott+=(neb_state[1])[j].fx*t[j][0]/t_mag;
					fdott+=(neb_state[1])[j].fy*t[j][1]/t_mag;
					fdott+=(neb_state[1])[j].fz*t[j][2]/t_mag;
				}
			}

			// determine force and velocity component in the direction of entire force
			// entire force on the image is stored in delr variable of the atom structure
			// magnitude of neb force is stored in delr[3] of the atom structure




			max_neb_forces[i] = 0.0;
			for(int j=0;j<n_atoms;j++)
			{
				//F_trueforce_component_perpendicular_to_tangent
				if(isdynamic[j])
				{
					(neb_state[1])[j].ux = (neb_state[1])[j].fx-fdott*t[j][0]/t_mag;
					(neb_state[1])[j].uy = (neb_state[1])[j].fy-fdott*t[j][1]/t_mag;
					(neb_state[1])[j].uz = (neb_state[1])[j].fz-fdott*t[j][2]/t_mag;
					double true_f = ((neb_state[1])[j].ux)*((neb_state[1])[j].ux)+((neb_state[1])[j].uy)*((neb_state[1])[j].uy)+((neb_state[1])[j].uz)*((neb_state[1])[j].uz);
					ftrue_mag += true_f;
					//fneb_max=0.0;
					//	double fspring_max = 0.0;
					if (ftrueforce_max < true_f) ftrueforce_max = true_f;

					// add F_spring
					double fspring = 0.0;
					if(!climb_states[i])
					{
						(neb_state[1])[j].ux += NEB_springK*(t_r_mag-t_l_mag)*t[j][0]/t_mag;
						(neb_state[1])[j].uy += NEB_springK*(t_r_mag-t_l_mag)*t[j][1]/t_mag;
						(neb_state[1])[j].uz += NEB_springK*(t_r_mag-t_l_mag)*t[j][2]/t_mag;
						spring_forces[i] = NEB_springK*(t_r_mag-t_l_mag);
						fspring = NEB_springK*(t_r_mag-t_l_mag)*t[j][0]/t_mag*NEB_springK*(t_r_mag-t_l_mag)*t[j][0]/t_mag;
						fspring += NEB_springK*(t_r_mag-t_l_mag)*t[j][1]/t_mag*NEB_springK*(t_r_mag-t_l_mag)*t[j][1]/t_mag;
						fspring +=NEB_springK*(t_r_mag-t_l_mag)*t[j][2]/t_mag*NEB_springK*(t_r_mag-t_l_mag)*t[j][2]/t_mag;
					}else
					{
						(neb_state[1])[j].ux += -1.0*fdott*t[j][0]/t_mag;
						(neb_state[1])[j].uy += -1.0*fdott*t[j][1]/t_mag;
						(neb_state[1])[j].uz += -1.0*fdott*t[j][2]/t_mag;
						spring_forces[i] = -1.0*fdott;
						fspring = 1.0*fdott*t[j][0]/t_mag*1.0*fdott*t[j][0]/t_mag;
						fspring += 1.0*fdott*t[j][1]/t_mag*1.0*fdott*t[j][1]/t_mag;
						fspring += 1.0*fdott*t[j][2]/t_mag*1.0*fdott*t[j][2]/t_mag;

					}
					if(fspring_max < fspring) fspring_max = fspring;

				}else
				{
					(neb_state[1])[j].ux = 0.0;
					(neb_state[1])[j].uy = 0.0;
					(neb_state[1])[j].uz = 0.0;
				}
				double local_neb_force = ((neb_state[1])[j].ux)*((neb_state[1])[j].ux)+((neb_state[1])[j].uy)*((neb_state[1])[j].uy)+((neb_state[1])[j].uz)*((neb_state[1])[j].uz);
				if (max_neb_forces[i] < local_neb_force) max_neb_forces[i] = local_neb_force;
				if(fneb_max< local_neb_force) fneb_max = local_neb_force;
				fneb_mag += local_neb_force;

				(neb_state[1])[j].delr[3] = local_neb_force;

			}

			fneb_mag = sqrt(fneb_mag);
			neb_forces[i] = fneb_mag;
			true_forces[i] = sqrt(ftrue_mag);

			if(fneb_mag > fneb_val) fneb_val = fneb_mag;

			int zero_value=0;
			double local_dt = NEB_dt;
			bool backtrack_exit = true;
			if (NEB_BACKTRACK_M) NEB_BACKTRACKFIRST = true;
			do
			{
				for(int j=0;j<n_atoms;j++)
				{
					if((NEB_BACKTRACK_M)&&(NEB_BACKTRACKFIRST))
					{
						neb_state[1][j].dummy[0] = neb_state[1][j].sx;
						neb_state[1][j].dummy[1] = neb_state[1][j].sy;
						neb_state[1][j].dummy[2] = neb_state[1][j].sz;

					}

					double r[3],s[3];

					s[0] = neb_state[1][j].sx;s[1] =neb_state[1][j].sy;s[2] = neb_state[1][j].sz;
					V3mulM3(s,Hcry,r);

					double r_change[3];
					r_change[0] = r_change[1] = r_change[2] = 0;
					if(isdynamic[j])
					{
						r_change[0] = local_dt*neb_state[1][j].ux;
						r_change[1] = local_dt*neb_state[1][j].uy;
						r_change[2] = local_dt*neb_state[1][j].uz;
						double r_change_mag = sqrt(r_change[0]*r_change[0]+r_change[1]*r_change[1]+r_change[2]*r_change[2]);
						if(r_change_mag > max_move)
						{
							for(int jj=0;jj<3;jj++)
							{
								r_change[jj] = r_change[jj]*max_move/r_change_mag;
							}
						}
						for(int jj=0;jj<3;jj++)
						{
							r[jj] += r_change[jj];
						}

						V3mulM3(r,Hcry_inv,s);
						neb_state[1][j].sx = s[0];neb_state[1][j].sy = s[1];neb_state[1][j].sz = s[2];
						if(neb_state[1][j].sx>=1.0) neb_state[1][j].sx = neb_state[1][j].sx-1;
						if(neb_state[1][j].sy>=1.0) neb_state[1][j].sy = neb_state[1][j].sy-1;
						if(neb_state[1][j].sz>=1.0) neb_state[1][j].sz = neb_state[1][j].sz-1;

						if(neb_state[1][j].sx<0.0) neb_state[1][j].sx = 1.0+neb_state[1][j].sx;
						if(neb_state[1][j].sy<0.0) neb_state[1][j].sy = 1.0+neb_state[1][j].sy;
						if(neb_state[1][j].sz<0.0) neb_state[1][j].sz = 1.0+neb_state[1][j].sz;
					}

				}
				if (NEB_BACKTRACK_M)
				{
					NEB_BACKTRACKFIRST = false;
					cout << iteration<<"\t"<< i<<"\t"<<local_dt<<"\t";
					backtrack_exit = accept_change(neb_state[1],filename_array[0], filename_array[total_length-1]);
					cout <<backtrack_exit<<"\t";
					if(!backtrack_exit)
					{

						for(int j=0;j<n_atoms;j++)
						{
							neb_state[1][j].sx = neb_state[1][j].dummy[0];
							neb_state[1][j].sy = neb_state[1][j].dummy[1];
							neb_state[1][j].sz = neb_state[1][j].dummy[2];

						}

						if(local_dt< NEB_BACKTRACK_dt_MIN)
						{
							backtrack_exit = true;
							cout << "exiting backtrack local_dt too small\n";
						}

						local_dt = local_dt*NEB_dt_FACTOR;

					}
					cout <<"\n";

				}
			}while(!backtrack_exit);
			if(atom!=NULL){	free(atom); atom=NULL;}
			atom = neb_state[1];
			save_lammps(filename_array[i],Hcry);
//			update_forces(filename_array[i]);
			for(int j=0;j<n_atoms;j++)
			{
				free(t[j]);free(t_l[j]);free(t_r[j]);
			}

			free(t);free(t_l);free(t_r);
			free(neb_state[0]);free(neb_state[1]);free(neb_state[2]);
			neb_state[0] = NULL; neb_state[1] = NULL; neb_state[2] = NULL;
			atom = NULL;
		}

		iteration++;
	}while ((ftol< fneb_val)&&(iteration<(NEB_MAX_ITER+1))&&(fabs(fneb_val-fneb_prev)>NEB_DELTOL));
}

