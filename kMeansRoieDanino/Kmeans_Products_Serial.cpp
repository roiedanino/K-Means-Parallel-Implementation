#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Cluster.h"

int setProductToTheClosestCluster(Product *p , Cluster** clusters, int numOfClusters);

Cluster** kMeans(int numOfClusters, Product** products , int numOfProducts);

Cluster** dynamicKMeans(Product** products, int numOfProducts, double maxError, int* numOfClusters);

double calcAvgError(Cluster** clusters, int numOfClusters, int numOfProducts);


int main(int argc,char *argv[])
{
	const double MAX_ERR = 1.9;
    int numOfProducts = 0, numOfClusters = 2;

    Product** allProducts = readProductsFromFile(&numOfProducts);

    Cluster** clusters = dynamicKMeans(allProducts,numOfProducts, MAX_ERR, &numOfClusters);

	printAllClusters((const Cluster**)clusters, numOfClusters);


    freeAllClusters(clusters, numOfClusters);

    freeAllProducts(allProducts,numOfProducts);

	return 0;
}

int setProductToTheClosestCluster(Product *p , Cluster** clusters, int numOfClusters)
{
	int i,minIndex=0;
	double minDist = MAXFLOAT;
    double distance = 0;

    for(i = 0 ; i < numOfClusters ; i++)
	{
		 distance = difference(p,&(clusters[i]->logicalCenter));

       if(distance == 0)
        {
            minDist = distance;
            minIndex = i;
            break;
        }

        if (distance < minDist)
        {
            minDist = distance;
            minIndex = i;
        }

	}

	if(clusters[minIndex] != p->cluster)
	{
		removeProduct(p->cluster,p->index);
		addProduct(clusters[minIndex],p);
        p->distance = minDist;
		return 1;
	}
	return 0;
}

Cluster** kMeans(int numOfClusters, Product** products , int numOfProducts)
{
	int i, hasChanged = 0, toClean = 1;
	Cluster** allClusters = (Cluster**)malloc(numOfClusters* sizeof(Cluster*));
	checkAllocation(allClusters);

    if(numOfClusters > numOfProducts)
        return NULL;

//#pragma omp parallel for
	for(i = 0 ; i < numOfClusters ; i++)
	{
		allClusters[i] = (Cluster*)malloc( sizeof(Cluster));
		checkAllocation(allClusters[i]);
		initCluster(allClusters[i],products[i]);
	}

	do
	{
		hasChanged = 0;

		for( i = 0 ; i < numOfProducts ; i++)
		{
            products[i]->index = toClean ? -1 : products[i]->index;
			hasChanged |= setProductToTheClosestCluster(products[i], allClusters, numOfClusters);
		}

        toClean = 0;

#pragma omp parallel for private(i)
		for(i = 0 ; i < numOfClusters ; i++)
		{
			setCenterProduct(allClusters[i]);
		}

	}while(hasChanged);

	return allClusters;
}

Cluster** dynamicKMeans(Product** products, int numOfProducts, double maxError, int* numOfClusters)
{
    *numOfClusters = 2;
    double error;

    Cluster** allClusters;

    do
    {

        allClusters = kMeans(*numOfClusters, products, numOfProducts);

        error = calcAvgError(allClusters, *numOfClusters, numOfProducts);

        if(maxError < error)
            freeAllClusters(allClusters,(*numOfClusters)++);

        printf("Num Of Clusters: %d   error: %lf   maxError: %lf\n",*numOfClusters,error,maxError);
        fflush(stdout);
    }while(maxError < error);

    return allClusters;
}

double calcAvgError(Cluster** clusters, int numOfClusters, int numOfProducts)
{
    int i, j;
    double sum = 0;

#pragma omp parallel for reduction(+ : sum) private(j)
    for(i = 0 ; i < numOfClusters ; i++)
    {
        for( j = 0 ; j < clusters[i]->numOfProducts ; j++)
        {
            sum += clusters[i]->products[j]->distance;
        }
    }

    return sum / numOfProducts;
}