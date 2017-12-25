//
// Created by Roie Danino on 13/09/2017.
//
#include "Product.h"

#ifndef __CLUSTER_H
#define __CLUSTER_H

typedef struct cluster_t
{
    int numOfProducts;

    double* centerTransactions;

    int id;

    double diameter;

}Cluster;


void initAllClusters(Cluster* allClusters, int numOfClusters,int prevNumOfClusters, Product* centers, int myRank);

void calculateDiameter(Cluster* cluster, Product* products);

void cleanClusters(Cluster* allClusters, int numOfClusters);

void freeAllClusters(Cluster* clusters, int numOfClusters);


#endif //__CLUSTER_H
