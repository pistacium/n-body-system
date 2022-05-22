#include<stdio.h>
#include<math.h>

void accel(int N, float r[][3], float (*m), double G, float a[][3]){
  float distsx[N][N], distsy[N][N], distsz[N][N], ir3[N][N], multx[N][N], multy[N][N], multz[N][N];
  float x[N], y[N], z[N], ax[N], ay[N], az[N];
  
  for(int i=0; i<N; i++){
    //put Nx3 position matrix into N long arrays of x y and z positions for ease of calculations
    x[i] = r[i][0];
    y[i] = r[i][1];
    z[i] = r[i][2];
  }

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      //calulate r_j - r_i
      if(i!=j){
	distsx[i][j] = x[j] - x[i];
	distsy[i][j] = y[j] - y[i];
	distsz[i][j] = z[j] - z[i];

	//calculate inverse magnitude of (r_j - r_i) cubed and put into ir3
	//extra factor of 0.1^2 because otherwise particles that collide will cause code to malfunction
	ir3[i][j] = pow(distsx[i][j],2) + pow(distsy[i][j],2) + pow(distsz[i][j],2) + pow(0.1,2);
	if(ir3[i][j]>0){
	  //at this point all ir3 should be > 0, but just in case so that other values do not blow up, make sure of this
	  ir3[i][j] = pow(ir3[i][j],(-3.0/2));
	}
      }
      else{
	//place 0 into i=j, otherwise code will take values from memory allocated in those areas previously
	distsx[i][j] = 0;
	distsy[i][j] = 0;
	distsz[i][j] = 0;
	ir3[i][j] = 0;
      }
    }
  }
  
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      //putting r_j - r_i / |r_j - r_i|^3 into multx multy and multz
      if(i!=j){
	multx[i][j] = distsx[i][j] * ir3[i][j];
	multy[i][j] = distsy[i][j] * ir3[i][j];
	multz[i][j] = distsz[i][j] * ir3[i][j];
      }
      else{
	multx[i][j] = 0;
	multy[i][j] = 0;
	multz[i][j] = 0;
      }
    }
  }
  for(int i=0; i<N; i++){
    //initializing ax ay and az as 0
    ax[i]=0;
    ay[i]=0;
    az[i]=0;
  }
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      //matrix multiplication of multx multy and multz by m (N length array of masses of particles)
      ax[i] += multx[i][j] * m[j];
      ay[i] += multy[i][j] * m[j];
      az[i] += multz[i][j] * m[j];
    }
  }
    
  for(int i=0; i<N; i++){
    //finally, multiply all accelerations by G
    a[i][0] = G * ax[i];
    a[i][1] = G * ay[i];
    a[i][2] = G * az[i];
  }
  
  
}

int main(){
  int N;
  float dt=0.1, tf, pt;

  //N - number of bodies
  //dt - time interval of simulation
  //tf - total simulation time
  //pt - time interval to print data to console

  
  //get user inputs for parameters
  printf("\n");
  printf("Enter total simulation time: ");
  scanf("%f",&tf);
  printf("Enter interval of time to display information: ");
  scanf("%f", &pt);
  printf("Enter total number of particles: ");
  scanf("%d",&N);
  printf("\n");
  
  float m[N];
  float r[N][3], v[N][3], a[N][3];
  double G = 6.67e-11;

  //m - N length array of masses
  //r - Nx3 size matrix of positions in 3d space for N bodies
  //v - Nx3 size matrix of velocities in 3d space for N bodies
  //a - Nx3 size matrix of accelerations in 3d space for N bodies


  //get user inputs for each body
  for(int i=0; i<N; i++){
    printf("Enter mass [kg] of particle: ");
    scanf("%f", &m[i]);
    printf("Enter x position [m] of particle: ");
    scanf("%f", &r[i][0]);
    printf("Enter y position [m] of particle: ");
    scanf("%f", &r[i][1]);
    printf("Enter z position [m] of particle: ");
    scanf("%f", &r[i][2]);
    printf("\n");
  }
  
  for(int i=0; i<N; i++){
    //initialize velocities as 0
    //(otherwise program will grab whatever was allocated at the memory addresses previously)
    v[i][0] = 0;
    v[i][1] = 0;
    v[i][2] = 0;
  }

  //initialize accelerations
  accel(N,r,m,G,a);
  
  float t;
  for(t=0; t<tf; t+=dt){
    //simulation loop for each time step
    //if t is close enough to a multiple of tp or t==0, print time and positions
    if(fmod(t,pt)<0.09 || pt==0) printf("t = %4f \n", t);
    
    for(int i=0; i<N; i++){
      //leap-frog: half step increment to velocity and full step increment to positions
      
      if(fmod(t,pt)<0.09 || pt==0) printf("Body %3d : r = [ %4f , %4f , %4f ]\n", i, r[i][0], r[i][1], r[i][2]);

      v[i][0] = v[i][0] + (dt/2.0) * a[i][0];
      v[i][1] = v[i][1] + (dt/2.0) * a[i][1];
      v[i][2] = v[i][2] + (dt/2.0) * a[i][2];

      r[i][0] = r[i][0] + dt*v[i][0];
      r[i][1] = r[i][1] + dt*v[i][1];
      r[i][2] = r[i][2] + dt*v[i][2];
    }

    //calculate new accelerations
    accel(N,r,m,G,a);

    for(int i=0; i<N; i++){
      //half step increment to velocities
      v[i][0] = v[i][0] + (dt/2.0) * a[i][0];
      v[i][1] = v[i][1] + (dt/2.0) * a[i][1];
      v[i][2] = v[i][2] + (dt/2.0) * a[i][2];
    }
    if(fmod(t,pt)<0.09 || pt==0) printf("\n---------------\n");
    
  }

  //print final state
  printf("FINAL STATE: \n");
  printf("t = %4f \n", t);
  for(int i=0; i<N; i++){
    printf("Body %3d : r = [ %4f , %4f , %4f ]\n", i, r[i][0], r[i][1], r[i][2]);
  }
      

  
}
