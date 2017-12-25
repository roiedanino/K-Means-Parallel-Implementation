#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <cstdlib>
#include "Cluster.h"
#include "Product.h"
#include "MPI_Functions.h"
#include "ProductDatabase.h"
#include "omp.h"


Cluster* dynamicKMeans(Product* originalProducts, ProductDatabase* myProducts, Info* info, int* numOfClusters);

double* gatherDiameters(Cluster* myClusters, ProductDatabase* productDatabase,int myNumOfClusters, int k, int numprocs, int myRank);

double calcQuality(Product* allCenters, int numOfCenters, double* diameters);

Product* kMeans(Product* originalProducts, ProductDatabase* productDatabase, Info* info, int k, int myNumOfClusters, int prevK, Cluster* myClusters, Product* firstCenters, int numprocs, int myRank, CudaMem* cudaMem);

void cudaMemAllocations(CudaMem* cmem, Info* info);

void cudaFreeAllocations(CudaMem* cmem);

void writeResultsToFile(Product* allCenters, int dimension, int numOfClusters, double qm);

int main(int argc,char *argv[])
{
    int numOfClusters = 2; 
    int numprocs, myid, myNumberOfProducts;

    Info info;
    Product *allProducts = NULL, *myProducts = NULL;

    ProductDatabase productDatabase;

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    MPI_Datatype InfoType, ProductType;

    createMPI_Info(&InfoType);


    if(myid ==  ROOT)
    {
   
        allProducts = readProductsFromFile(&info);
        printf("Num of products: %d\nDimension: %d\nMaxK: %d\nLimit Of Iterations: %d\nQuality:%lf\n",info.numOfProducts,
               info.dimension, info.maxK, info.limitOfIterations, info.qualityMeasure);
		fflush(stdout);
    }

	//Sharing Meta data such as QM, MaxK and dimension
    MPI_Bcast(&info, 1, InfoType, ROOT,MPI_COMM_WORLD);

    createMPI_Product(&ProductType, info.dimension);

    myNumberOfProducts = info.numOfProducts / numprocs;

    myProducts = (Product*)calloc(sizeof(Product),(size_t)myNumberOfProducts + (myid == ROOT ? info.numOfProducts % numprocs : 0));
    checkAllocation(myProducts);

	//Equally scattering products
    scatterProducts(allProducts, info.numOfProducts, info.dimension, myNumberOfProducts, ProductType, myProducts, myid, ROOT);

	//ProductDatabase - process's own set of products
    initProductDatabase(&productDatabase, myProducts, myNumberOfProducts, info.dimension, ProductType);

	//Root Process takes the remaining of the products
    if(myid == ROOT)
    {
        int remaining = info.numOfProducts % numprocs;
        int i;
        for(i = 0 ; i < remaining ; i++)
        {
            addProduct(&productDatabase, nullptr, &allProducts[info.numOfProducts - remaining + i]);
        }
    }

    Cluster* clusters = dynamicKMeans(allProducts, &productDatabase, &info, &numOfClusters);

	MPI_Finalize();

    freeAllClusters(clusters, numOfClusters);

	return 0;
}


Cluster* dynamicKMeans(Product* originalProducts, ProductDatabase* myProducts, Info* info, int* numOfClusters)
{
    double error = INT_MAX, t1,t2;
    int myRank, numprocs, stop = 0, k = 2, prevNumOfClusters = 0, physicalSize = 2;
	
	CudaMem cudaMemory;

    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 

    Cluster* myClusters = (Cluster*)malloc(physicalSize*sizeof(Cluster)); 
    checkAllocation(myClusters);

    Product* firstCenters = (Product*)malloc(sizeof(Product)*physicalSize); 
    checkAllocation(firstCenters);

	cudaMemAllocations(&cudaMemory, info);

    *numOfClusters = k / numprocs;

    if(k%numprocs > myRank)
        (*numOfClusters)++;

	if (myRank == ROOT)
	{
		t1 = MPI_Wtime();
	}

    do
    {
        Product* allCenters = kMeans(originalProducts, myProducts, info, k, *numOfClusters, prevNumOfClusters, myClusters, firstCenters, numprocs, myRank, &cudaMemory);

		//Each process calculate his own clusters diamters and sends it to the Root process
        double* diameters = gatherDiameters(myClusters, myProducts, *numOfClusters, k,numprocs, myRank);

        if(myRank == ROOT)
        {
            error = calcQuality(allCenters, k, diameters);
            stop = info->qualityMeasure >= error || k >= info->maxK;
            free(diameters);
        }

        MPI_Bcast(&stop,1,MPI_INT,ROOT,MPI_COMM_WORLD);


        if(!stop)
        {
            cleanClusters(myClusters, *numOfClusters);
            prevNumOfClusters = *numOfClusters;

            k++;

            *numOfClusters = k / numprocs;

			//Each process adds a cluster in his turn
            if(k % numprocs > myRank)
            {
                (*numOfClusters)++;
            }

			//dynamically increamenting arrays memory allocations
            if(*numOfClusters >= physicalSize)
            {
                physicalSize *= 2;

                myClusters = (Cluster *) realloc(myClusters, physicalSize * sizeof(Cluster));
                checkAllocation(myClusters);

                firstCenters = (Product *) realloc(firstCenters, physicalSize * sizeof(Product));
                checkAllocation(firstCenters);
            }

            free(allCenters);
        }

        else if(myRank == ROOT)
        {
			t2 = MPI_Wtime();
			printf("\nTime: %lf\n", (t2 - t1));
            writeResultsToFile(allCenters, info->dimension, k, error);
        }

    }while(!stop);

	printf("Num Of Clusters: %d   current quality: %lf   Desired Quality: %lf\n",k, error, info->qualityMeasure);
	fflush(stdout);

	cudaFreeAllocations(&cudaMemory);

    return myClusters;
}

Product* kMeans(Product* originalProducts, ProductDatabase* productDatabase, Info* info, int k, int myNumOfClusters, int prevK, Cluster* myClusters, Product* firstCenters, int numprocs, int myRank, CudaMem* cudaMem)
{
    int hasChanged = 0, tempHasChanged = 0, toClean = 1, dimension = info->dimension, iteration = 0;
    int allCentersBufferSize = ((k / numprocs) + (k % numprocs == 0? 0: 1))*numprocs;

    Product* allCenters = (Product*)calloc((size_t)allCentersBufferSize,sizeof(Product)); //TODO CHANGE HERE
    checkAllocation(allCenters);

	//scattering first k clusters centers to the other processes
    scatterVcenters(originalProducts, k, dimension, firstCenters, myNumOfClusters, productDatabase->ProductType, numprocs, myRank);

    initAllClusters(myClusters, myNumOfClusters, prevK, firstCenters, myRank);

    do
    {
		//Sharing clusters centers
        shareCenters(allCenters, allCentersBufferSize, myClusters, myNumOfClusters, k, myRank, productDatabase, numprocs);

		//Matching Points to clusters with CUDA
		hasChanged = setEveryProductToHisClusterWithCuda(productDatabase, allCenters, k, myRank, toClean, cudaMem);
		
        toClean = 0;

		//Send products to the processes that has the clusters they belong to
        sendProductsToTheirClusters(numprocs, productDatabase, myRank, myClusters, myNumOfClusters);
			
		tempHasChanged = hasChanged;
        MPI_Allreduce(&tempHasChanged, &hasChanged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		
		//Each process calculates his own clusters new centers
		if(hasChanged) 
			calculateCenters(myClusters,myNumOfClusters, productDatabase, allCenters);        
		
    }while(hasChanged && ++iteration < info->limitOfIterations);

    return allCenters;
}




double calcQuality(Product* allCenters, int numOfCenters, double* diameters)
{
    int i, j;
    double sum = 0;

#pragma omp parallel for  private(j) reduction(+:sum)
    for(i = 0 ; i < numOfCenters ; i++)
    {
        for(j = i+1 ; j < numOfCenters ; j++)
        {
            sum += (diameters[i]+diameters[j])/difference(&(allCenters[i]),&(allCenters[j]));
        }
    }

    return sum / (numOfCenters*(numOfCenters-1));
}

double* gatherDiameters(Cluster* myClusters, ProductDatabase* productDatabase,int myNumOfClusters, int k, int numprocs, int myRank)
{
    int i, *recvCounts = NULL, *displ = NULL;
    double* myDiameters = (double*)calloc((size_t)myNumOfClusters , sizeof(double));
    double* allDiameters = NULL;
    Product* myProductsMatrix = NULL;

    if(myRank == ROOT)//TODO CHANGE HERE
    {
        allDiameters = (double *) calloc((size_t) myNumOfClusters * numprocs, sizeof(double));
        checkAllocation(allDiameters);

        recvCounts = (int*)calloc((size_t)numprocs , sizeof(int));
        checkAllocation(recvCounts);

        displ = (int*)calloc((size_t)numprocs , sizeof(int));
        checkAllocation(displ);

        for (i = 0; i < numprocs; i++)
        {
            recvCounts[i] = k / numprocs;

            if (k % numprocs > i)
            {
                recvCounts[i] += 1;
            }

            displ[i] = i != 0 ? displ[i - 1] + recvCounts[i - 1] : 0;
        }
    }

    myProductsMatrix = fromDatabaseTo_ByCluster_Matrix(productDatabase, myNumOfClusters);

#pragma omp parallel for
    for(i = 0 ; i < myNumOfClusters ; i++)
    {
        calculateDiameter(&myClusters[i], &myProductsMatrix[i*productDatabase->numOfProducts]);

        myDiameters[i] = myClusters[i].diameter;
    }

    MPI_Gatherv(myDiameters, myNumOfClusters, MPI_DOUBLE, allDiameters, recvCounts, displ, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    free(myDiameters);
    free(recvCounts);
    free(displ);
    free(myProductsMatrix);

    return allDiameters;
}



void writeResultsToFile(Product* allCenters, int dimension, int numOfClusters, double qm)
{
    int i, j;
    FILE* resultFile = fopen("kmeans_results.txt", "w");

    fprintf(resultFile,"Number of clusters with the best measure: K = %d QM = %lf\n\nCenters Of Clusters:\n\n",numOfClusters,qm);

    for(i = 0; i < numOfClusters ; i++)
    {

        fprintf(resultFile, "C%d  \n", (i+1));

        for(j = 0 ; j < dimension ; j++)
        {
            fprintf(resultFile,"%4.3lf    ", allCenters[i].transactions[j]);
        }
        fprintf(resultFile,"\n\n");
    }
    fclose(resultFile);
}

// Pre - Allocations for cuda
void cudaMemAllocations(CudaMem* cmem, Info* info)
{
	cudaError_t cudaStatus;
	int dev_distancesSize = info->numOfProducts*info->dimension*info->maxK;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	handleErrors(cudaStatus, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?", cmem->dev_productsTrans, cmem->dev_centersTrans, cmem->dev_distances);

	cudaStatus = cudaMalloc((void**)&cmem->dev_productsTrans, info->numOfProducts*info->dimension* sizeof(double));
	handleErrors(cudaStatus, "cudaMalloc failed!", cmem->dev_productsTrans, cmem->dev_centersTrans, cmem->dev_distances);

	cudaStatus = cudaMalloc((void**)&cmem->dev_centersTrans, info->dimension *info->maxK* sizeof(double));
	handleErrors(cudaStatus, "cudaMalloc failed!", cmem->dev_productsTrans, cmem->dev_centersTrans, cmem->dev_distances);

	cudaStatus = cudaMalloc((void**)&cmem->dev_distances, dev_distancesSize * sizeof(double));
	handleErrors(cudaStatus, "cudaMalloc failed!", cmem->dev_productsTrans, cmem->dev_centersTrans, cmem->dev_distances);

	cudaStatus = cudaMalloc((void**)&cmem->dev_ClusterIds, info->numOfProducts * sizeof(int));
	handleErrors(cudaStatus, "cudaMalloc failed!", cmem->dev_productsTrans, cmem->dev_centersTrans, cmem->dev_distances);
}

void cudaFreeAllocations(CudaMem* cmem)
{
	cudaFree(cmem->dev_centersTrans);
	cudaFree(cmem->dev_ClusterIds);
	cudaFree(cmem->dev_distances);
	cudaFree(cmem->dev_productsTrans);
}