//create a OpenCl RNG function to feed random numbers inside the kernel for optimized 2d random walk simulation 
// the RNG need to be able to feed the parallel threads with random numbers
// The RNG function only recieve a position i,j a the number of thread as seeds and return a random number
// the main kernel function only recieve of the external C code the position and return average a temperature  and their variance 

// first create the RNG function suitable for parallelization
// Use the Mersenne Twister RNG algorithm




