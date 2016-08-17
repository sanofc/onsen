/*
 *  advect.cpp
 *  smoke3D
 *
 */

#include "advect.h"
#include "utility.h"

static FLOAT u_ref( FLOAT ****u, int dir, int i, int j, int k, int N ) {
	if( dir == 0 )
		return u[0][MAX(0,MIN(N,i))][MAX(0,MIN(N-1,j))][MAX(0,MIN(N-1,k))];
	else if( dir == 1 )
		return u[1][MAX(0,MIN(N-1,i))][MAX(0,MIN(N,j))][MAX(0,MIN(N-1,k))];
	else
		return u[2][MAX(0,MIN(N-1,i))][MAX(0,MIN(N-1,j))][MAX(0,MIN(N,k))];
}

		   
void semiLagrangian( FLOAT ***d, FLOAT ***d0, int width, int height, int depth, FLOAT ****u, int N, FLOAT dt ) {
	OPENMP_FOR for( int n=0; n<width*height*depth; n++ ) {
		int i = (n%(width*height))%width;
		int j = (n%(width*height))/width;
		int k = n/(width*height);
		
		d[i][j][k] = interp( d0, width, height, depth, i-N*u[0][i][j][k]*dt, j-N*u[1][i][j][k]*dt, k-N*u[2][i][j][k]*dt );
	}
}

void advect::advect( FLOAT ****u, FLOAT ***c, FLOAT ***s, FLOAT ***t, int N, FLOAT dt ) {
	
	// Compute Fluid Velocity At Each Staggered Faces And Concentration Cell Centers
	static FLOAT ***ux[3] = { alloc3D(N+1,N,N), alloc3D(N+1,N,N), alloc3D(N+1,N,N) };
	static FLOAT ***uy[3] = { alloc3D(N,N+1,N), alloc3D(N,N+1,N), alloc3D(N,N+1,N) };
	static FLOAT ***uz[3] = { alloc3D(N,N,N+1), alloc3D(N,N,N+1), alloc3D(N,N,N+1) };
	static FLOAT ***uc[3] = { alloc3D(N,N,N), alloc3D(N,N,N), alloc3D(N,N,N) };
	static FLOAT ***out[6] = { alloc3D(N+1,N,N), alloc3D(N,N+1,N), alloc3D(N,N,N+1), alloc3D(N,N,N) , alloc3D(N,N,N) ,alloc3D(N,N,N)};
	
	FOR_EVERY_X_FLOW {
		ux[0][i][j][k] = u[0][i][j][k];
		ux[1][i][j][k] = (u_ref(u,1,i-1,j,k,N)+u_ref(u,1,i-1,j+1,k,N)+u_ref(u,1,i,j,k,N)+u_ref(u,1,i,j+1,k,N))/4.0;
		ux[2][i][j][k] = (u_ref(u,2,i-1,j,k,N)+u_ref(u,2,i-1,j,k+1,N)+u_ref(u,2,i,j,k,N)+u_ref(u,2,i,j,k+1,N))/4.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		uy[0][i][j][k] = (u_ref(u,0,i,j-1,k,N)+u_ref(u,0,i+1,j-1,k,N)+u_ref(u,0,i,j,k,N)+u_ref(u,0,i+1,j,k,N))/4.0;
		uy[1][i][j][k] = u[1][i][j][k];
		uy[2][i][j][k] = (u_ref(u,2,i,j-1,k,N)+u_ref(u,2,i,j-1,k+1,N)+u_ref(u,2,i,j,k,N)+u_ref(u,2,i,j,k+1,N))/4.0;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		uz[0][i][j][k] = (u_ref(u,0,i,j,k-1,N)+u_ref(u,0,i+1,j,k-1,N)+u_ref(u,0,i,j,k,N)+u_ref(u,0,i+1,j,k,N))/4.0;
		uz[1][i][j][k] = (u_ref(u,1,i,j,k-1,N)+u_ref(u,1,i,j+1,k-1,N)+u_ref(u,1,i,j,k,N)+u_ref(u,1,i,j+1,k,N))/4.0;
		uz[2][i][j][k] = u[2][i][j][k];
	} END_FOR
	
	FOR_EVERY_CELL {
		uc[0][i][j][k] = 0.5*u[0][i][j][k]+0.5*u[0][i+1][j][k];
		uc[1][i][j][k] = 0.5*u[1][i][j][k]+0.5*u[1][i][j+1][k];
		uc[2][i][j][k] = 0.5*u[2][i][j][k]+0.5*u[2][i][j][k+1];
	} END_FOR
	
	// BackTrace X Flow
	semiLagrangian( out[0], u[0], N+1, N, N, ux, N, dt );
	
	// BackTrace Y Flow
	semiLagrangian( out[1], u[1], N, N+1, N, uy, N, dt );
	
	// BackTrace Z Flow
	semiLagrangian( out[2], u[2], N, N, N+1, uz, N, dt );
	
	// BackTrace Vapor Concentration
	semiLagrangian( out[3], c, N, N, N, uc, N, dt );

	// BackTrace Steam Concentration
	semiLagrangian( out[4], s, N, N, N, uc, N, dt );

	// BackTrace Temperature
	semiLagrangian( out[5], t, N, N, N, uc, N, dt);
	
	// Copy Back To The Result
	copy3D(u[0],out[0],N);
	copy3D(u[1],out[1],N);
	copy3D(u[2],out[2],N);
	copy3D(c,out[3],N);
	copy3D(s,out[4],N);
	copy3D(t,out[5],N);
}
