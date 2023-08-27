// Define constants
#define DOMAIN_SIZE 10.0f // Size of the square domain
#define NUM_DIRECTIONS 4   // North, South, East, West
#define LOCAL_WORKGROUP_SIZE 64 // Local workgroup size for local memory

// Initialize XOR-Shift constants
__constant uint4 xorshift_state = (uint4)(0u, 0u, 0u, 0u);

// Function prototype for xorshift128
uint4 xorshift128();


// Initialize the random number generator with skip ahead functionality
void init_xorshift_skip_ahead(uint4 seed, uint skip_count) {
    xorshift_state = seed;
    for (uint i = 0; i < skip_count; i++) {
        xorshift128();
        }
    }

    // Generate random numbers in batches using XOR-Shift
    uint4 xorshift128() {
        uint4 state = xorshift_state;
        uint t = state.x;
        t ^= t << 11;
        t ^= t >> 8;
        state.x = state.y;
        state.y = state.z;
        state.z = state.w;
        state.w = (state.w ^ (state.w >> 19)) ^ (t ^ (t >> 19));
        xorshift_state = state;
        return state;
    }

    // Main kernel function for random walk simulation
    __kernel void random_walk(__global float* temperatures,
                              uint numWalkers, uint numSteps, uint4 seed,
                              float h, float pN, float pS, float pE, float pW,
                              float T_N, float T_S, float T_E, float T_W,
                              float x_ini, float y_ini) { 
        size_t global_id = get_global_id(0);
        size_t local_id = get_local_id(0);
        size_t local_size = get_local_size(0);
        
        // Allocate local memory for a portion of the temperatures array
        __local float local_temperatures[LOCAL_WORKGROUP_SIZE];
        
        for (uint walker = 0; walker < numWalkers; walker++) {
            // Calculate the skip ahead value based on the walker index
            uint skip_count = walker * numSteps + global_id;
            
            // Initialize XOR-Shift state with skip ahead
            uint4 local_seed = seed + (uint4)(skip_count, skip_count, skip_count, skip_count);
            init_xorshift_skip_ahead(local_seed, skip_count);
            
            // Load a portion of the temperatures array into local memory
            uint local_index = local_id;
            while (local_index < numWalkers) {
                local_temperatures[local_id] = temperatures[local_index];
                local_index += local_size;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            
            // Initialize the starting position for each walker
            float2 walker_position = (float2)(x_ini, y_ini);
            
            for (uint step = 0; step < numSteps; step++) {
                uint4 random_nums = xorshift128();
                
                // Choose a random direction based on probabilities
                float prob = (float)random_nums.x / 4294967295.0f; // Convert to [0, 1]
                uint direction = 0;
                if (prob < pN) {
                    direction = 0; // North
                } else if (prob < pN + pS) {
                    direction = 1; // South
                } else if (prob < pN + pS + pE) {
                    direction = 2; // East
                } else {
                    direction = 3; // West
                }
                
                // Update walker's position based on chosen direction
                switch (direction) {
                    case 0: // North
                        walker_position.y += h;
                        break;
                    case 1: // South
                        walker_position.y -= h;
                        break;
                    case 2: // East
                        walker_position.x += h;
                        break;
                    case 3: // West
                        walker_position.x -= h;
                        break;
                }
                
                // Check if the walker is out of bounds and score temperature
                if (walker_position.x < 0.0f) {
                    temperatures[walker] = T_W;
                    break;
                }
                if (walker_position.x >= DOMAIN_SIZE) {
                temperatures[walker] = T_E;
                break;
            }
            if (walker_position.y < 0.0f) {
                temperatures[walker] = T_S;
                break;
            }
            if (walker_position.y >= DOMAIN_SIZE) {
                temperatures[walker] = T_N;
                break;
            }
        }
    }
}

// Global reduction kernel to calculate average temperature and variance
__kernel void calculate_average_and_variance(__global float* temperatures,
                                             uint numWalkers,
                                             __global float* averageTemperature,
                                             __global float* variance) {
    float sumTemperature = 0.0f;
    float sumSquaredTemperature = 0.0f;
    
    for (uint walker = 0; walker < numWalkers; walker++) {
        float temperature = temperatures[walker];
        sumTemperature += temperature;
        sumSquaredTemperature += temperature * temperature;
    }
    
    *averageTemperature = sumTemperature / numWalkers;
    
    float meanSquaredTemperature = sumSquaredTemperature / numWalkers;
    *variance = meanSquaredTemperature - (*averageTemperature) * (*averageTemperature);
}
