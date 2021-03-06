
#include "allvars.h"
#include "grid.h"
#include "tools.h"
#include "profiles.h"

void compute_profiles()
{
   int          i,k,ic,jc,kc,l,ii,jj,kk,next,ibin,in,m,NumGrid;
   double       xc[3],xt[3],dx[3],vt[3],dist,Rho,GridSize[3];
   double       dR,DeltaDiff,DeltaCum,MinDist,MaxDist,GAP;
   double       CGal[NumProfileBins],Suma[NumProfileBins],CVel[NumProfileBins];
   double       VRad,rm,ri,rs,Vol,DeltaMax,Radius;
   vector <int> Indx; 
   char         OutFile[MAXCHAR];
   int          NumQuery;
   struct query Query;
   struct grid  *GridList;
   FILE         *fd;
   clock_t      t;

   fprintf(logfile,"\n COMPUTING VOID PROFILES \n");
   t = clock();

   NumGrid = (int)(BoxSize/ProxyGridSize);
   GridList = (struct grid *) malloc(NumGrid*NumGrid*NumGrid*sizeof(struct grid));
   build_grid_list(Tracer,NumTrac,GridList,NumGrid,GridSize,false);

   // Only for true voids

   for (i=0; i<NumVoid; i++) 
       if (Void[i].ToF) 
	  Indx.push_back(i);       

   dR = (log10(MaxProfileDist)-log10(MinProfileDist))/(double)NumProfileBins;
   
   // Selecciono grides

   MinDist = 0.0;
   MaxDist = 0.0;
   for (i=0; i<NumVoid; i++) {
       if (!Void[i].ToF) continue;
       if (Void[i].Rad > MaxDist) MaxDist = Void[i].Rad;
   }
   MaxDist *= MaxProfileDist;
   query_grid(&Query,GridSize,MinDist,MaxDist);
   NumQuery = Query.i.size();
   GAP = 0.5*sqrt(3.0)*max_grid_size(GridSize);
   
   fprintf(logfile," | MinDist - MaxDist = %5.3f - %5.3f [Mpc/h], %d grids \n",MinDist,MaxDist,NumQuery);
   fflush(logfile);

   #pragma omp parallel for default(none) schedule(dynamic)                   \
    shared(NumVoid,Void,Tracer,NumQuery,Query,dR,NumGrid,GridSize,GridList,   \
           Indx,MeanNumTrac,LBox,NumProfileBins,MinProfileDist,MaxProfileDist,\
	   WriteProfiles,PathProfiles,GAP)                                    \
    private(i,m,k,ii,jj,kk,l,Radius,ic,jc,kc,CVel,CGal,Suma,xc,xt,dx,vt,next, \
	    dist,VRad,ibin,DeltaMax,ri,rm,rs,Rho,Vol,DeltaDiff,DeltaCum,fd,in,\
	    OutFile)

   for (i=0; i<(int)Indx.size(); i++) {
       
       for (k=0; k<NumProfileBins; k++) {
           CVel[k] = 0.0;
           CGal[k] = 0.0;
           Suma[k] = 0.0;
       }
       
       Radius = Void[Indx[i]].Rad;
       for (k=0; k<3; k++) 
	   xc[k] = Void[Indx[i]].Pos[k];
       Void[Indx[i]].Dtype = 0.0;
      
       ic = (int)(xc[0]/GridSize[0]);
       jc = (int)(xc[1]/GridSize[1]);
       kc = (int)(xc[2]/GridSize[2]);

       for (in=0; in<NumQuery; in++) {

       	   ii = Query.i[in]; 
	   jj = Query.j[in]; 
	   kk = Query.k[in]; 

	   dist = (double)(ii*ii)*(GridSize[0]*GridSize[0])
	        + (double)(jj*jj)*(GridSize[1]*GridSize[1])
	        + (double)(kk*kk)*(GridSize[2]*GridSize[2]);
	   dist = sqrt(dist);

	   if (dist > MaxProfileDist*Radius+GAP) continue;

       	   ii = periodic_grid(ii + ic, NumGrid); 
	   jj = periodic_grid(jj + jc, NumGrid); 
	   kk = periodic_grid(kk + kc, NumGrid); 	       
	   
	   l = index_1d(ii,jj,kk,NumGrid);

           if (GridList[l].NumMem == 0) continue;

           for (m=0; m<GridList[l].NumMem; m++) {
		
	       next = GridList[l].Member[m];	

	       for (k=0; k<3; k++) {
	           xt[k] = (double)Tracer[next].Pos[k];	 
	           vt[k] = (double)(Tracer[next].Vel[k] - Void[Indx[i]].Vel[k]);
	           dx[k] = periodic_delta(xt[k] - xc[k],LBox[k])/Radius;
	       }

               dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	       if (dist > MinProfileDist && dist < MaxProfileDist) {
           
	          ibin = (int)((log10(dist)-log10(MinProfileDist))/dR);

	          VRad = vt[0]*dx[0] + vt[1]*dx[1] + vt[2]*dx[2];
	          VRad /= dist;

	          CVel[ibin] += VRad;
	          CGal[ibin] += 1.0; 

	       }
	   }
       }
       
       for (k=0; k<NumProfileBins; k++) {
	   if (CGal[k] < 3.0) {
	       CGal[k] = 0.0;
               CVel[k] = 0.0;
           } else {
               CVel[k] /= CGal[k];	   
	   }
       }

       if (WriteProfiles == 1) {
          sprintf(OutFile,"%s/void_%d.dat",PathProfiles,i);
          fd = safe_open(OutFile,"w");
       }

       DeltaMax = -1.0;
       for (k=0; k<NumProfileBins; k++) {

	   for (kk=0; kk<=k; kk++) 
	       Suma[k] += CGal[kk]; 

	   ri = (double)(k    )*dR + log10(MinProfileDist);
	   rm = (double)(k+0.5)*dR + log10(MinProfileDist);
	   rs = (double)(k+1.0)*dR + log10(MinProfileDist);

	   ri = pow(10.0,ri)*Radius;
	   rm = pow(10.0,rm)*Radius;
	   rs = pow(10.0,rs)*Radius;

	   Vol = (4.0/3.0)*PI*(pow(rs,3) - pow(ri,3));
	   Rho = CGal[k]/Vol;
	   DeltaDiff = Rho/MeanNumTrac - 1.0;

	   Vol = (4.0/3.0)*PI*pow(rs,3);
	   Rho = Suma[k]/Vol;
           DeltaCum = Rho/MeanNumTrac - 1.0;

    	   if (WriteProfiles == 1) 
	      fprintf(fd,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",
			  ri,rm,rs,CVel[k],DeltaDiff,DeltaCum,Radius);
	   
	   if (rs < 2.0*Radius || rs > 3.0*Radius) continue;

	   if (DeltaCum > DeltaMax) DeltaMax = DeltaCum;

       }
       
       if (WriteProfiles == 1)
          fclose(fd);

       Void[Indx[i]].Dtype = DeltaMax;

   } 

   Indx.clear();
   free_query_grid(&Query);
   free_grid_list(GridList,NumGrid);
   
   StepName.push_back("Computing profiles");
   StepTime.push_back(get_time(t,OMPcores));

}
