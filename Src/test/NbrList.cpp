/*
 * NbrList.cpp
 *
 *  Created on: Jul 19, 2009
 *      Author: kedar
 */



/* --------------------------------------------------------------------- */
#include <NbrList.h>

int cellno(int ix, int iy, int iz)
{
	int result;

	result = 1 + (ix - 1 + mx)%mx + ((iy - 1 + my)%my)*mx +
		((iz - 1 + mz)%mz)*mx*my;
	return(result);
}


void mapcells()
{
	int ix, iy, iz, imap;

	for (iz = 1; iz <= mz; iz++)
		for (iy = 1; iy <= my; iy++)
			for (ix = 1; ix <= mx; ix++)
			{

				imap = (cellno(ix,iy,iz) - 1)*13;

				map1[imap+ 1] = cellno(ix+1,iy  ,iz  );
				map1[imap+ 2] = cellno(ix+1,iy+1,iz  );
				map1[imap+ 3] = cellno(ix  ,iy+1,iz  );
				map1[imap+ 4] = cellno(ix-1,iy+1,iz  );
				map1[imap+ 5] = cellno(ix+1,iy  ,iz-1);
				map1[imap+ 6] = cellno(ix+1,iy+1,iz-1);
				map1[imap+ 7] = cellno(ix  ,iy+1,iz-1);
				map1[imap+ 8] = cellno(ix-1,iy+1,iz-1);
				map1[imap+ 9] = cellno(ix+1,iy  ,iz+1);
				map1[imap+10] = cellno(ix+1,iy+1,iz+1);
				map1[imap+11] = cellno(ix  ,iy+1,iz+1);
				map1[imap+12] = cellno(ix-1,iy+1,iz+1);
				map1[imap+13] = cellno(ix  ,iy  ,iz+1);

			}
}


void make_nbr_lst(double H1[3][3])
// This is a hybrid subroutine using cell-link-list algorithm
// to create the Verlet neighbor list. Uses extra memory for
// holding the Verlet neighbor list but eliminates the necessity
// to refresh the cell-link list at every call to force subroutine.
{
	TEST_BOOL = true;
	MAX_COORDNUM=0;
	int i,j,nlst,icell,jcell,nabor,jcell0;

	double sxi,syi,szi,sij[3],rij[3],rijsq;
	cout << ncell<<" VALUES \n";
	nlst = 0;
	// for (i = 1; i <= n; i++) nn[i]=0;
	// reset the head list of cell-link list
	for (i = 1; i <= ncell; i++) head[i] = -1;
	// reset the displacements, the neigbor list pointers. reconstruct
	//   head list, link-list until next neighbor list update
	for (i = 0; i < n; i++)
	{
		icell = 1 + ((int)((atom[i].sx)*(mx*1.0))) + ((int)((atom[i].sy)*(my*1.0)))*mx +((int)((atom[i].sz)*(mz*1.0)))*mx*my;
		// cout <<i<<"\t"<< atom[i].rx<<"\t"<<atom[i].ry<<"\t"<<atom[i].rz<<"\t"<<icell<<"\n";
		if(icell<=0)
		{
			cout <<" icell is less than or equal to zero - does not make sense ";
			cout << i<<"\t"<<atom[i].rx<<"\t"<<atom[i].ry<<"\t"<<atom[i].rz<<"\n";
		}
		if(icell>ncell)
			{
				cout <<"ICELL > N CELL WHAT IS WRONG\t"<<icell<<"\t\n";
				cout << i<<"\t"<<atom[i].rx<<"\t"<<atom[i].ry<<"\t"<<atom[i].rz<<"\n";
			}
		list[i] = head[icell];
		head[icell] = i;
		nbr_ptr[i] = 0;
		nbr_ptr1[i] = -1;


		//	if(icell==10000)
		//	{
		//		cout << icell<<"\t"<< i<<"\t"<<list[i]<<"\t"<<head[icell]<<"\n";
		//	}

	}

	if(TEST_BOOL) { cout << "first\n";cout << ncell<<"\n";}
	// loop over each cell for neighbor list update
	for (icell = 1; icell <= ncell; icell++) {
		i = head[icell];
		// loop over each atom in icell
		//    cout << icell<<"sd\n";
		while (i >= 0) {
			int i_start = nlst;
			nbr_ptr[i] = nlst + 1;
			sxi = atom[i].sx;
			syi = atom[i].sy;
			szi = atom[i].sz;
			j = list[i];
			//	  if(icell==10000)
			//	  {
			//	  cout << j<<"\n";
			//	  }

			// cout << i<< "\t"<<icell<<"\n";
			if((j>n )||(j<-1))
			{
				cout << "j GONE AWRY!!\t"<<j<<"\n";
			}
			// search neighbors of i in icell
			if(i!=j)
			{
				if(j==i) cout << "the same atom "<<i<<" is twice in the root cell\n";

				while (j >= 0) {
					sij[0] = sxi - atom[j].sx;
					sij[1] = syi - atom[j].sy;
					sij[2] = szi - atom[j].sz;
					sij[0] = sij[0]-(int)(sij[0]*2);
					sij[1] = sij[1]-(int)(sij[1]*2);
					sij[2] = sij[2]-(int)(sij[2]*2);

					V3mulM3(sij,H1,rij);
					rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];



					//	if(icell==10000) cout << rijsq<<"\t";
					if (rijsq <= rlstsq)
					{
						nlst++;
						nbr_lst[nlst] = j;
					}

					if(rijsq<=0.05)
					{
						cout << i<<"\t"<<atom[i].elem<<"\t"<<j<<"\t"<<atom[j].elem<<"\t"<<atom[i].sx-atom[j].sx;
						cout<<"\t"<<atom[i].sy-atom[j].sy<< "\t" << atom[i].sz-atom[j].sz<<"\t"<<rijsq<<" examples of too close atoms \n";
							//exit(1);
					}

					j = list[j];
					//	if(icell==10000) cout << i<<"\tnext j is\t"<<j<<"\n";
				}
			}
			// if(TEST_BOOL) cout << "second\n";
			jcell0 = 13 * (icell - 1);
			// loop over neighbor cells of icell to search neighbors of i
			for (nabor = 1; nabor <= 13; nabor++) {
				jcell = map1[jcell0 + nabor];
				j = head[jcell];
				// search neighbors of i in jcell
				if((j>n )||(j<-1))
				{
					cout << "j GONE AWRY!!\t"<<j<<"\n";
				}
				if(i!=j)
				{
					if(j==i) cout << "the same atom "<<i<<" is twice in the adjacent cells\n";

					while (j >= 0) {
						sij[0] = sxi - atom[j].sx;
						sij[1] = syi - atom[j].sy;
						sij[2] = szi - atom[j].sz;

						sij[0] = sij[0]-(int)(sij[0]*2);
						sij[1] = sij[1]-(int)(sij[1]*2);
						sij[2] = sij[2]-(int)(sij[2]*2);


						V3mulM3(sij,H1,rij);
						rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];

						if (rijsq <= rlstsq)
						{
							nlst++;
							nbr_lst[nlst] = j;
						}
						j = list[j];
					}
				}
			}
			nbr_ptr1[i] = nlst;
			if (MAX_COORDNUM < (nlst-i_start))
			{
				MAX_COORDNUM = nlst-i_start;
			}
			i = list[i];
		}
	}
	cout << nlst*1.0/n<<" avg nlst \n";
	cout << MAX_COORDNUM <<"\t is the maximum coordination number\n";
}
