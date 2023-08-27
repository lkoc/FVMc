// #include <CL/cl.h>
#include "C:\\OpenCL_SDK\\OpenCL-SDK\\external\\OpenCL-Headers\\CL\cl.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Define the platform and device index

#define PLATFORM_INDEX 0
#define DEVICE_INDEX 0

// Define the number of walkers and steps
#define NUM_WALKERS 1000
#define NUM_STEPS 100

// Function to load OpenCL source code from a file
char* load_source_code(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        exit(1);
    }

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    char* source_code = (char*)malloc(file_size + 1);
    if (!source_code) {
        perror("Failed to allocate memory for source code");
        fclose(file);
        exit(1);
    }

    fread(source_code, 1, file_size, file);
    source_code[file_size] = '\0';

    fclose(file);
    return source_code;
}

int main() {
    // Initialize OpenCL and get platform/device
    cl_platform_id platform;
    clGetPlatformIDs(1, &platform, NULL);

    cl_device_id device;
    clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);

    // Get the maximum workgroup size
    size_t maxLocalWorkgroupSize;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &maxLocalWorkgroupSize, NULL);

    // Create a context and command queue
    cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, NULL);
    cl_command_queue command_queue = clCreateCommandQueue(context, device, 0, NULL);

    // Load OpenCL program source code from file
    const char* source_filename = "random_walk.cl";
    char* source_code = load_source_code(source_filename);

    // Create a program and build it
    cl_program program = clCreateProgramWithSource(context, 1, &source_code, NULL, NULL);
    clBuildProgram(program, 1, &device, NULL, NULL, NULL);

    // Create kernels
    cl_kernel random_walk_kernel = clCreateKernel(program, "random_walk", NULL);
    cl_kernel calc_avg_var_kernel = clCreateKernel(program, "calculate_average_and_variance", NULL);

    // Set kernel arguments
    // ... (set kernel arguments using clSetKernelArg)

    // Set LOCAL_WORKGROUP_SIZE to the maximum value
    uint LOCAL_WORKGROUP_SIZE = (uint) maxLocalWorkgroupSize;

    // Print device information
    char deviceName[128];
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
    printf("Device: %s\n", deviceName);
    printf("Max Local Workgroup Size: %zu\n", maxLocalWorkgroupSize);

    // Launch the random walk kernel
    size_t global_work_size = NUM_WALKERS;
    size_t local_work_size = LOCAL_WORKGROUP_SIZE;

    // Measure execution time
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    clEnqueueNDRangeKernel(command_queue, random_walk_kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);

    // Launch the average/variance calculation kernel
    clEnqueueNDRangeKernel(command_queue, calc_avg_var_kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);

    // Retrieve results
    float averageTemperature, variance;
    // ... (read results using clEnqueueReadBuffer)

    // Print results
    gettimeofday(&end_time, NULL);
    double execution_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) * 1e-6;
    printf("Average Temperature: %f\n", averageTemperature);
    printf("Variance: %f\n", variance);
    printf("Number of Walkers: %d\n", NUM_WALKERS);
    printf("Execution Time: %f seconds\n", execution_time);

    // Clean up resources
    clReleaseKernel(random_walk_kernel);
    clReleaseKernel(calc_avg_var_kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);

    return 0;
}
