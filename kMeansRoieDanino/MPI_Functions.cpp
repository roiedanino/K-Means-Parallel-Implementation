//
// Created by Roie Danino on 15/09/2017.
//

#include "MPI_Functions.h"


 // #################################################### MPI FUNCTIONS #############################################################################


// Data Types Creations

void createMPI_Product(MPI_Datatype* Product_Type ,int dimension)
{
  
    MPI_Datatype type[] = {MPI_CHAR, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};

    int blocklen[] = { LEN, 1, 1, 1, 1, 1, 1};

    MPI_Aint disp[7];

    disp[0] = offsetof(Product, productName);
    disp[1] = offsetof(Product,transactions);
    disp[2] = offsetof(Product,dimension);
    disp[3] = offsetof(Product, clusterId);
    disp[4] = offsetof(Product, newClusterId);
    disp[5] = offsetof(Product, index);
    disp[6] = offsetof(Product, rank);

    MPI_Type_create_struct(7, blocklen, disp, type, Product_Type);
    MPI_Type_commit(Product_Type);

}

void createMPI_Info(MPI_Datatype* InfoType)
{
    MPI_Datatype type[] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};

    int blocklen[] = {1,1,1,1,1};

    MPI_Aint disp[5];

    disp[0] = offsetof(Info, numOfProducts);
    disp[1] = offsetof(Info, dimension);
    disp[2] = offsetof(Info, maxK);
    disp[3] = offsetof(Info, limitOfIterations);
    disp[4] = offsetof(Info, qualityMeasure);

    MPI_Type_create_struct(5, blocklen, disp, type, InfoType);
    MPI_Type_commit(InfoType);
}


// Scatter And Scatter-v for products

void scatterProducts(Product* products, int numOfProducts, int dimension, int numOfProductsEach, MPI_Datatype ProductType, Product* myProducts, int rank, int root)
{
    double* tempTransactions =  NULL;

    double* myProductsTransactions = (double*)calloc((size_t)numOfProductsEach*dimension, sizeof(double));
    checkAllocation(myProductsTransactions);

    if(rank == root)
    {
        tempTransactions = getProductsTransactions(products,numOfProducts);
    }

    MPI_Scatter(products,numOfProductsEach,ProductType,myProducts,numOfProductsEach,ProductType,root,MPI_COMM_WORLD);

    MPI_Scatter(tempTransactions, dimension*numOfProductsEach,MPI_DOUBLE,myProductsTransactions,dimension*numOfProductsEach,MPI_DOUBLE,root,MPI_COMM_WORLD);

    setProductsTransactions(myProducts,numOfProductsEach,myProductsTransactions,dimension);

    if(rank == root)
    {
        free(tempTransactions);
    }
}


void scatterVcenters(Product* centers, int k, int dimension, Product* myCenters, int myNumOfCenters, MPI_Datatype ProductType, int numprocs, int myRank)
{

    int i, *sendCounts = NULL, *displ = NULL;
    double *tempTransactions = NULL, *myProductsTransactions = (double *) calloc(
            (size_t) myNumOfCenters * dimension + 1, sizeof(double));
    checkAllocation(myProductsTransactions);

    if (myRank == ROOT)
    {
        tempTransactions = getProductsTransactions(centers, k);

        sendCounts = (int *) calloc((size_t) numprocs * 2, sizeof(int));
        checkAllocation(sendCounts);

        displ = (int *) calloc((size_t) numprocs * 2, sizeof(int));
        checkAllocation(displ);

        for (i = 0; i < numprocs; i++) 
		{
            sendCounts[i] = k / numprocs;

            if (k % numprocs > i) {
                sendCounts[i] += 1;
            }
            sendCounts[i + numprocs] = sendCounts[i] * centers[0].dimension;

            displ[i] = (i != 0 ? displ[i - 1] + sendCounts[i - 1] : 0);
            displ[i + numprocs] = (i != 0 ? displ[(i - 1) + numprocs] + sendCounts[(i - 1) + numprocs] : 0);
        }
    }


    MPI_Scatterv(centers, sendCounts, displ, ProductType, myCenters, myNumOfCenters, ProductType, ROOT, MPI_COMM_WORLD);
    MPI_Scatterv(tempTransactions, sendCounts + numprocs, displ + numprocs, MPI_DOUBLE, myProductsTransactions,
                 myNumOfCenters * dimension, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    if (myNumOfCenters > 0)
    {
        setProductsTransactions(myCenters, myNumOfCenters, myProductsTransactions, dimension);
    }

    if(myRank == ROOT)
    {
        free(tempTransactions);
        free(sendCounts);
        free(displ);
    }
}

void allGatherProducts(Product *products, int numOfProducts, int dimension, int numOfProductsEach,
                       MPI_Datatype ProductType, Product *myProducts, int rank)
{
    double* gatheredTransactions = (double *) malloc(sizeof(double) * dimension * numOfProducts);
    checkAllocation(gatheredTransactions);

    double* tempTransactions = getProductsTransactions(myProducts, numOfProductsEach);

    MPI_Gather(myProducts, numOfProductsEach, ProductType, products, numOfProductsEach, ProductType, ROOT, MPI_COMM_WORLD);

    MPI_Gather(tempTransactions,dimension*numOfProductsEach,MPI_DOUBLE, gatheredTransactions,dimension*numOfProductsEach,MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

    MPI_Bcast(products,numOfProducts,ProductType,ROOT,MPI_COMM_WORLD);

    MPI_Bcast(gatheredTransactions,numOfProducts*dimension, MPI_DOUBLE,ROOT,MPI_COMM_WORLD);

    setProductsTransactions(products, numOfProducts, gatheredTransactions, dimension);

    free(tempTransactions);
}


int setProductToCluster(Product* product, int clusterIndex, Product* centers)
{
	if (centers[clusterIndex].clusterId != product->clusterId || product->rank != centers[clusterIndex].rank)
	{
		product->newClusterId = centers[clusterIndex].clusterId;
		product->rank = centers[clusterIndex].rank;
		return 1;
	}
	return 0;
}

int setEveryProductToHisClusterWithCuda(ProductDatabase* productDatabase, Product* centers, int numOfCenters, int myRank, int toClean, CudaMem* cudaMem)
{
    int i, hasChanged = 0;
	double* allProductsTransactions = getProductsTransactions(productDatabase->products, productDatabase->numOfProducts);
	double* allCentersTransactions = getProductsTransactions(centers, numOfCenters);
	
	int* clustersIdsPerProduct = (int*)calloc(productDatabase->numOfProducts*productDatabase->dimension*numOfCenters, sizeof(double));
	checkAllocation(clustersIdsPerProduct);
	
	getProductsClustersIdsWithCuda(allProductsTransactions, productDatabase->numOfProducts, allCentersTransactions, numOfCenters, productDatabase->dimension, clustersIdsPerProduct, cudaMem);

#pragma omp parallel for reduction(+:hasChanged)
    for( i = 0 ; i < productDatabase->numOfProducts ; i++)
    {
        if(toClean)
        {
            cleanProduct(&productDatabase->products[i]);
        }

        hasChanged |= setProductToCluster(&productDatabase->products[i], clustersIdsPerProduct[i], centers);
    }

    return hasChanged;
}

void shareCenters(Product* centers, int centersBufferSize, Cluster* myClusters, int numOfClusters, int k,
                  int myRank, ProductDatabase* productDatabase, int numprocs)
{

    int i, mySendBufferSize = centersBufferSize / numprocs , size  = centersBufferSize;
	Product* tempCenters = (Product*)calloc(centersBufferSize, sizeof(Product));
	checkAllocation(tempCenters);
	
#pragma omp parallel for
    for( i = 0 ; i < numOfClusters; i++)
    {
		tempCenters[i].transactions = (double*)calloc(productDatabase->dimension, sizeof(double)); 
		memcpy(tempCenters[i].transactions, myClusters[i].centerTransactions,sizeof(double)*productDatabase->dimension);
        tempCenters[i].rank = myRank;
        tempCenters[i].clusterId = myClusters[i].id;
        tempCenters[i].newClusterId = -1;
        tempCenters[i].dimension = productDatabase->dimension;
    }

    if(k%numprocs != 0)
    {
		tempCenters[numOfClusters].transactions = (double*)calloc(productDatabase->dimension, sizeof(double)); 
		tempCenters[numOfClusters].dimension = productDatabase->dimension;
        tempCenters[numOfClusters].newClusterId = PADDING_CENTER;
    }

	
    allGatherProducts(centers, centersBufferSize, productDatabase->dimension, mySendBufferSize, productDatabase->ProductType, tempCenters, myRank);

	trimToLogicalSize(size, mySendBufferSize, centers);

	free(tempCenters);
}

void trimToLogicalSize(int size, int mySendBufferSize, Product* centers)
{
	int i;

	for(i = size - mySendBufferSize - 1; i >= 0; i--)
    {
        if(centers[i].newClusterId == PADDING_CENTER)
        {
            if (centers[size - 1].newClusterId == PADDING_CENTER)
            {
                size--;
            }

            centers[i] = centers[size - 1];
            size--;
        }
    }
}

void sendProductsToTheirClusters(int numprocs, ProductDatabase* productDatabase, int myRank, Cluster* myClusters, int numOfClusters)
{
    int i, j, dimension  = productDatabase->dimension, removeCounter = 0;

	//A Matrix : NumOfProcesses * Num Of Products To Send To each of them
    Product *productsToSend = (Product *)calloc(numprocs * productDatabase->numOfProducts, sizeof(Product)); 
    checkAllocation(productsToSend);

    int *countersForEachProcess = (int *) calloc((size_t) numprocs, sizeof(int));
    checkAllocation(countersForEachProcess);

    int* indicesToRemove = (int*)calloc((size_t)productDatabase->numOfProducts, sizeof(int));
    checkAllocation(indicesToRemove);

    Product** productsToReceive = (Product**)malloc((size_t) sizeof(Product*)*numprocs);
    checkAllocation(productsToReceive);

    int* amounts = (int*)calloc((size_t)numprocs, sizeof(int));
    checkAllocation(amounts);

	removeCounter = toSendOrToSave(productDatabase, productsToSend, countersForEachProcess, indicesToRemove, myClusters, myRank);

	paddWithDefaultProduct(productDatabase, productsToSend, countersForEachProcess, numprocs);

    for(i = 0 ; i < numprocs ; i++)
    {
        int BuffSize = productDatabase->numOfProducts;
        MPI_Bcast(&BuffSize,1,MPI_INT, i, MPI_COMM_WORLD);

        productsToReceive[i] = (Product*)calloc((size_t)BuffSize, sizeof(Product));

        MPI_Scatter(countersForEachProcess, 1, MPI_INT, &amounts[i], 1, MPI_INT, i, MPI_COMM_WORLD);

        scatterProducts(productsToSend, numprocs*BuffSize, productDatabase->dimension, BuffSize, productDatabase->ProductType, productsToReceive[i], myRank, i);                                  
    }

    for( i = removeCounter - 1 ; i >= 0 ; i--)
    {
        removeProduct(productDatabase, &myClusters[productDatabase->products[indicesToRemove[i]].clusterId], indicesToRemove[i]);
    }

    for(i = 0 ; i < numprocs ; i++)
    {
        for (j = 0; j < amounts[i]; j++)
        {
            addProduct(productDatabase, &myClusters[productsToReceive[i][j].newClusterId], &productsToReceive[i][j]);
        }
    }


    free(amounts);
    free(productsToReceive);
    free(indicesToRemove);
    free(countersForEachProcess);
    free(productsToSend);

}

void paddWithDefaultProduct(ProductDatabase* productDatabase, Product* productsToSend, int* countersForEachProcess, int numprocs)
{
	  int i, j;

#pragma omp parallel for private(i)
    for (j = 0; j < numprocs; j++) 
	{
        for (i = 0; i < productDatabase->numOfProducts; i++) 
		{
            if(i >= countersForEachProcess[j])
            {
                productsToSend[j * productDatabase->numOfProducts + i] = productDatabase->products[0];
            }
        }
    }
}

int toSendOrToSave(ProductDatabase* productDatabase,  Product* productsToSend, int* countersForEachProcess, int* indicesToRemove, Cluster* myClusters, int myRank)
{
	int i, toRank, removeCounter = 0;

	for(i = 0; i < productDatabase->numOfProducts; i++)
		{
            toRank = productDatabase->products[i].rank;

                if (toRank != myRank)
                {
                    productsToSend[toRank * productDatabase->numOfProducts + countersForEachProcess[toRank]++] = productDatabase->products[i];
                    indicesToRemove[removeCounter++] = i;
                }
                else if (productDatabase->products[i].newClusterId != NO_CLUSTER)
                {
                    myClusters[productDatabase->products[i].clusterId].numOfProducts--;
                    myClusters[productDatabase->products[i].newClusterId].numOfProducts++;
                    productDatabase->products[i].clusterId = productDatabase->products[i].newClusterId;
                    productDatabase->products[i].newClusterId = NO_CLUSTER;
                }
        }

	return removeCounter;
}

void calculateCenters(Cluster* myClusters, int myNumOfClusters, ProductDatabase* productDatabase, Product* allCenters)
{
	int i,j,t, dimension = productDatabase->dimension;

#pragma omp parallel for private(j, t)
        for(i = 0 ; i < myNumOfClusters ; i++)
        {
            if(myClusters[i].numOfProducts > 0)
            {
                memset(myClusters[i].centerTransactions, 0, dimension * sizeof(double));
            }

            for(j = 0 ; j < productDatabase->numOfProducts ; j++)
            {

                if(productDatabase->products[j].clusterId == myClusters[i].id)
                {
                    for(t = 0 ; t < dimension; t++)
                    {
                        myClusters[i].centerTransactions[t] += productDatabase->products[j].transactions[t] / myClusters[i].numOfProducts;
                    }
                }
            }
                memcpy(allCenters[i].transactions, myClusters[i].centerTransactions, sizeof(double)*dimension);
		}
}

//Arranging products into a by-cluster matrix 
Product* fromDatabaseTo_ByCluster_Matrix(ProductDatabase* productDatabase, int numOfClusters)
{
    int i,clusterId = 0, numOfProducts = productDatabase->numOfProducts;

    Product* productMat = (Product*)malloc(sizeof(Product) * numOfProducts * numOfClusters);
    checkAllocation(productMat);

    int *countersForEachCluster = (int *) calloc((size_t) numOfClusters, sizeof(int));
    checkAllocation(countersForEachCluster);

    for(i = 0 ; i < numOfProducts ; i++)
    {
        clusterId = productDatabase->products[i].clusterId;
        productMat[clusterId * numOfProducts + countersForEachCluster[clusterId]++] = productDatabase->products[i];
    }

    free(countersForEachCluster);

    return productMat;
}
