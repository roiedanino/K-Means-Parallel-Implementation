//
// Created by Roie Danino on 15/09/2017.
//

#ifndef __MPI_FUNCTIONS_H
#define __MPI_FUNCTIONS_H

#include <mpi.h>
#include "Product.h"
#include "Cluster.h"
#include "ProductDatabase.h"
#include "kernel.h"
#include "TechnicalFunctions.h"

#include <stddef.h>
#include "omp.h"

#define ROOT 0

#define PADDING_CENTER -2

void createMPI_Info(MPI_Datatype* InfoType);

void createMPI_Product(MPI_Datatype* Product_Type, int dimension);

void scatterProducts(Product* products, int numOfProducts, int dimension, int numOfProductsEach, MPI_Datatype ProductType, Product* myProducts, int rank, int root);

void scatterVcenters(Product* centers, int k, int dimension, Product* myCenters, int myNumOfCenters, MPI_Datatype ProductType, int numprocs, int myRank);

void allGatherProducts(Product *products, int numOfProducts, int dimension, int numOfProductsEach,
                       MPI_Datatype ProductType, Product *myProducts, int rank);

int setProductToCluster(Product* product, int clusterIndex, Product* centers);

int setEveryProductToHisClusterWithCuda(ProductDatabase* productDatabase, Product* centers, int numOfCenters, int myRank, int toClean, CudaMem* cudaMem);

void shareCenters(Product* centers, int centersBufferSize, Cluster* myClusters, int numOfClusters, int k,
                  int myRank, ProductDatabase* productDatabase, int numprocs);

void sendProductsToTheirClusters(int numprocs, ProductDatabase* productDatabase, int myRank, Cluster* myClusters, int numOfClusters);

void calculateCenters(Cluster* myClusters, int myNumOfClusters, ProductDatabase* productDatabase, Product* allCenters);

Product* fromDatabaseTo_ByCluster_Matrix(ProductDatabase* productDatabase, int numOfClusters);

int toSendOrToSave(ProductDatabase* productDatabase,  Product* productsToSend, int* countersForEachProcess, int* indicesToRemove, Cluster* myClusters, int myRank);

//Technical Functions
void paddWithDefaultProduct(ProductDatabase* productDatabase, Product* productsToSend, int* countersForEachProcess, int numprocs);

void trimToLogicalSize(int size, int mySendBufferSize, Product* centers);

#endif //__MPI_FUNCTIONS_H
