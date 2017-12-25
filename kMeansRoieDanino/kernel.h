#ifndef __KERNEL_H
#define __KERNEL_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
	double *dev_productsTrans;
	double *dev_centersTrans;
	double *dev_distances;
	int* dev_ClusterIds;
}CudaMem;


cudaError_t getProductsClustersIdsWithCuda(double* productsTransactions, int numOfProducts, double* centersTransactions, int numOfCenters, int dimension, int* clusterIds,
	CudaMem* cudaMem);

void handleErrors(cudaError_t cudaStatus, const char* errorMessage, double* dev_arr, double* dev_histo, double* dev_localisto);

#endif