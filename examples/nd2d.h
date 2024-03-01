#include <cmath>
#include <vector>

#define MESHMIN 3
#define mesh_f_(i,j)      mesh[nx*((j)-1)+(i)-1]

int nd2d(int nx, int ny, std::vector<int> &mesh, std::vector<int> &perm)
{
   int status = 0;
   int i, j, count;
   int nxtop, nxbot, nyleft, nyright;
   int midcol, midrow;

   if (nx <= MESHMIN && ny <= MESHMIN) {
      /* at a leaf, walk thru in column major */
      count = 0;
      for (j = 1; j <=ny; j++) 
         for (i=1; i<=nx; i++) {
            perm[count] = mesh_f_(i,j);
            count++;
         }  
   } 
   else if ( nx > ny) {
      /* cut mesh in the middle row */
      midrow = (int)ceil((double)nx/2.0);

      /* order the top half */
      nxtop = midrow-1;
      std::vector<int> meshtop(nxtop*ny);
      std::vector<int> ptop(nxtop*ny, 0);
      for (j = 1; j<=ny; j++) 
         for (i=1; i<=nxtop; i++) {
	    meshtop[nxtop*((j)-1)+(i)-1] = mesh_f_(i,j);
         }

      nd2d(nxtop, ny, meshtop, ptop);
      for (i = 0; i<nxtop*ny; i++)
         perm[i] = ptop[i];

      ptop.clear();
      meshtop.clear();

      /* order the bottom half */
      nxbot = nx - nxtop - 1;
      std::vector<int> meshbot(nxbot*ny); 
      for (j = 1; j<=ny; j++) 
	 for (i=1; i<=nxbot; i++) { 
	   meshbot[nxbot*((j)-1)+(i)-1] = mesh_f_(nxtop+1+i,j);
         }

      std::vector<int> pbot(nxbot*ny, 0);
      nd2d(nxbot, ny, meshbot, pbot);
      for (i = 0; i< nxbot*ny; i++)
         perm[nxtop*ny+i] = pbot[i];

      /* append the seperator */
      for (j = 0; j<ny; j++)
	 perm[(nxtop+nxbot)*ny+j]=mesh_f_(midrow,j+1); 

   }
   else {
      /* cut mesh in the middle column */
      midcol = (int)ceil((double)ny/2);
      nyleft = midcol-1;
      std::vector<int> meshleft(nx*nyleft); 
      for (j = 1; j<=nyleft; j++) {
         for (i=1; i<=nx; i++) 
	   meshleft[(j-1)*nx+i-1] = mesh_f_(i,j);
      }
      std::vector<int> pleft(nx*nyleft, 0);
      nd2d(nx, nyleft, meshleft, pleft);
      for (i = 0; i<nyleft*nx; i++)
         perm[i] = pleft[i];

      pleft.clear();
      meshleft.clear();

      nyright = ny - nyleft - 1;
      std::vector<int> meshright(nx*nyright); 
      for (j = 1; j<=nyright; j++) {
         for (i=1; i<=nx; i++) 
	   meshright[(j-1)*nx+i-1] = mesh_f_(i,nyleft+1+j);
      }
      std::vector<int> pright(nx*nyright, 0);
      nd2d(nx, nyright, meshright, pright);
      for (i = 0; i<nyright*nx; i++)
         perm[nyleft*nx+i] = pright[i];

      for (i = 0; i<nx; i++)
	perm[(nyleft+nyright)*nx+i] = mesh_f_(i+1,midcol);

   } 
   return status;
}
