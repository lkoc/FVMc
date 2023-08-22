__kernel void heat2d(
    __global const float *u0,
    __global float *u1,
    __global const float *Br, // Robin
    __global const float *alpha, // alpha
    __global const float *f, // q

// alpha= k(x,y) * delta_t/(Cp(x,y) * rho(x,y)) --> (x,y)) alpha es "a"
//is the thermal diffusivity, a material-specific quantity 
// f = q * delta_t/(rho*Cp) = q* delta_t *a/K

    float gamma, //,
//    float dt,
    float dx,
    float dy,
    float A_r)
 {
  //Get total number of cells
  int nx = get_global_size(0);
  int ny = get_global_size(1);
  int i = get_global_id(0);
  int j = get_global_id(1);
  int ii = get_local_id(0);
  int jj = get_local_id(1);

  //Calculate the four indices of our neighbouring cells
  // j es eje vertical de arriba hacia abajo
  // i es eje horizontal de izquierda a derecha
  int centro = j * nx + i;
  int sur    = (j + 1) * nx + i; //abajo
  int norte  = (j - 1) * nx + i; //arriba
  int este   = j * nx + (i + 1); //derecha
  int oeste   = j * nx + (i - 1); //izquierda

  //Internall cells
  if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
//       u1[centro] = u0[centro]
//       + gamma * (u0[oeste]  - 2 * u0[centro] + u0[este])
//       + gamma * (u0[sur]    - 2 * u0[centro] + u0[norte]);

u1[centro] = u0[centro] 
  + (alpha[este]/(2*dx) - alpha[oeste]/(2*dx))*(u0[este]/(2*dx) - u0[oeste]/(2*dx))
  + (alpha[sur]/(2*dy) - alpha[norte]/(2*dy))*(u0[sur]/(2*dy) - u0[norte]/(2*dy))
  + (u0[este]/(dx*dx) + u0[oeste]/(dx*dx) - 2*u0[centro]/(dx*dx))*alpha[centro]  
  + (u0[sur]/(dy*dy) + u0[norte]/(dy*dy) - 2*u0[centro]/(dy*dy))*alpha[centro] 
  + f[centro];
           }
  // Boundary conditions (Diritchlet)
  else {
       u1[centro] = u0[centro];
        };
//Condiciones de frontera
 // q=0 Neumman BC left
//barrier(CLK_GLOBAL_MEM_FENCE);
 if ( i ==0 && j > 0 && j < ny - 1 ){
         u1[centro] = u0[centro]
         + gamma * (u0[este]  - 2 * u0[centro] + u0[este])
         + gamma * (u0[sur] - 2 * u0[centro] + u0[norte]);
      //u1[centro] = u1[este];
           }
 // q=0 Neumman BC right
//  if ( i ==nx -1 && j > 0 && j < ny - 1  ){
//          u1[centro] = u0[centro]
//          + gamma * (u0[oeste]  - 2 * u0[centro] + u0[oeste])
//          + gamma * (u0[sur] - 2 * u0[centro] + u0[norte]);
//          //u1[centro] = u1[oeste];
//          }

 // q=0 Neumman BC top
  if (i > 0 && i < nx - 1 && j ==0){
          u1[centro] = u0[centro]
          + gamma * (u0[oeste]  - 2 * u0[centro] + u0[este])
          + gamma * (u0[sur] - 2 * u0[centro] + u0[sur]);
          //u1[centro] = u1[oeste];
          }

 // q= h*(t - t_inf) Robin BC Top
// if (i > 0 && i < nx - 1 && j ==0){
//      u1[centro] = u0[centro]
//    + gamma * ( u0[oeste]  - 2 * u0[centro] + u0[este])
//      + gamma * ( u0[sur] - 2*u0[centro] + (u0[sur] - A_r*u0[centro] + Br[centro]));
//        };

 // q=0 Neumman BC bottom
  if (i > 0 && i < nx - 1 && j ==ny-1){
          u1[centro] = u0[centro]
          + gamma * (u0[oeste]  - 2 * u0[centro] + u0[este])
          + gamma * (u0[norte] - 2 * u0[centro] + u0[norte]);
          //u1[centro] = u1[oeste];
          }
        
 // Esquina sup-izq
 if (i ==0 && j == 0){
         u1[centro]  = (0.5)*(u0[sur] + u0[este]);}
 //Esquina sup-der
 if (i == nx-1 && j == 0){
         u1[centro] =  (0.5)*(u0[sur] + u0[oeste]);}
 //Esquina inf-derecha
 if (i == nx-1 && j == ny-1){
         u1[centro] = (0.5)*(u0[norte] + u0[oeste]);}
 //Esquina inf-izq
 if (i ==0 && j == ny-1){
         u1[centro] = (0.5)*(u0[norte] + u0[este]);}
 }