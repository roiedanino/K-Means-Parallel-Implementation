#include "kernel.h"


__global__ void distance(double* products, int numOfProducts, double* centers, int numOfCenters, int dimension, double* results, int productsPerBlock)
{
	int i, xid = threadIdx.x, yid = threadIdx.y, blockId = blockIdx.x;
	double result = 0, delta = 0;
	
	if (blockIdx.x == gridDim.x - 1 && numOfProducts % blockDim.x <= xid)
		return;
	
	for (i = 0; i < dimension; i++)
	{
		delta = (products[(blockId*productsPerBlock + xid)*dimension + i] - centers[yid*dimension + i]);
		
		result += delta*delta;
	}

	results[numOfProducts*yid + (blockId*productsPerBlock + xid)]  = result;
}

__global__ void findClosestCluster(double* distances, int numOfCenters, int numOfProducts, int productsPerBlock, int* centersIds)
{
	int i, xid = threadIdx.x, blockId = blockIdx.x;
	double minIndex = 0, minDistance, tempDistance;

	if (blockIdx.x == gridDim.x - 1 && numOfProducts % blockDim.x <= xid)
		return;
	
	minDistance = distances[productsPerBlock*blockId + xid];

	for (i = 1; i < numOfCenters; i++)
	{
		tempDistance = distances[productsPerBlock*blockId + xid + i*numOfProducts];
		
		if (minDistance > tempDistance)
		{
			minIndex = i;
			minDistance = tempDistance;
		}
	}
	centersIds[productsPerBlock*blockId + xid] = minIndex;
}


cudaError_t getProductsClustersIdsWithCuda(double* productsTransactions, int numOfProducts, double* centersTransactions, int numOfCenters, int dimension, int* clusterIds, 
	CudaMem* cudaMem)
{
	double *dev_productsTrans = cudaMem->dev_productsTrans;
	double *dev_centersTrans = cudaMem->dev_centersTrans;
	double *dev_distances = cudaMem->dev_distances;
	int* dev_ClusterIds = cudaMem->dev_ClusterIds;

	int productsPerBlock;
	int numBlocks;
	int extraBlock;
	cudaDeviceProp prop;
	cudaError_t cudaStatus;

	cudaStatus = cudaGetDeviceProperties(&prop, 0);
	handleErrors(cudaStatus, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?", dev_productsTrans, dev_centersTrans, dev_distances);

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_productsTrans, productsTransactions, numOfProducts*dimension * sizeof(double), cudaMemcpyHostToDevice);
	handleErrors(cudaStatus, "cudaMemcpy failed!", dev_productsTrans, dev_centersTrans, dev_distances);

	cudaStatus = cudaMemcpy(dev_centersTrans, centersTransactions, dimension * numOfCenters * sizeof(double), cudaMemcpyHostToDevice);
	handleErrors(cudaStatus, "cudaMemcpy failed!", dev_productsTrans, dev_centersTrans, dev_distances);



	productsPerBlock = prop.maxThreadsPerBlock / numOfCenters;
	dim3 dim(productsPerBlock, numOfCenters);
	numBlocks = numOfProducts / productsPerBlock;
	extraBlock = numOfProducts % productsPerBlock == 0 ? 0 : 1;

	// Launch a kernel on the GPU with one thread for each element.
	distance << <numBlocks + extraBlock, dim >> > (dev_productsTrans, numOfProducts, dev_centersTrans, numOfCenters, dimension, dev_distances, productsPerBlock);

	cudaStatus = cudaDeviceSynchronize();
	handleErrors(cudaStatus, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", dev_productsTrans, dev_centersTrans, dev_distances);


	productsPerBlock = prop.maxThreadsPerBlock;
	numBlocks = numOfProducts / productsPerBlock;
	extraBlock = numOfProducts % productsPerBlock == 0 ? 0 : 1;

	findClosestCluster << <numBlocks + extraBlock, productsPerBlock >> > (dev_distances, numOfCenters, numOfProducts, productsPerBlock, dev_ClusterIds);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	handleErrors(cudaStatus, "addKernel launch failed: %s\n", dev_productsTrans, dev_centersTrans, dev_distances);


	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	handleErrors(cudaStatus, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", dev_productsTrans, dev_centersTrans, dev_distances);



	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(clusterIds, dev_ClusterIds, numOfProducts * sizeof(int), cudaMemcpyDeviceToHost);
	handleErrors(cudaStatus, "cudaMemcpy failed!", dev_productsTrans, dev_centersTrans, dev_distances);

	return cudaStatus;
}

void handleErrors(cudaError_t cudaStatus, const char* errorMessage ,double* dev_arr, double* dev_histo, double* dev_localisto)
{
	 if (cudaStatus != cudaSuccess) 
	 {
		 cudaFree(dev_arr);
		 cudaFree(dev_histo);
		 cudaFree(dev_localisto);
	 }
}

