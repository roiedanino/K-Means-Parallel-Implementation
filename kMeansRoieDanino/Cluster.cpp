//
// Created by Roie Danino on 13/09/2017.
//
#include <stdio.h>
#include "Cluster.h"
#include "Product.h"

void initAllClusters(Cluster* allClusters, int numOfClusters,int prevNumOfClusters, Product* centers, int myRank)
{
    int i;

    for(i = 0 ; i < numOfClusters ; i++)
    {
        if(i >= prevNumOfClusters)
        {
            allClusters[i].id = i;
            allClusters[i].centerTransactions = (double*)malloc(sizeof(double)*centers[i].dimension);
        }

        allClusters[i].numOfProducts = 0;
        allClusters[i].diameter = 0;
        allClusters[i].id = i;
        memcpy(allClusters[i].centerTransactions, centers[i].transactions, sizeof(double)*centers[i].dimension);
    }
}

void calculateDiameter(Cluster* cluster, Product* products)
{
    int i, j;
    double diameter = 0, diff;


    for(i = 0 ; i < cluster->numOfProducts ; i++)
    {
        for(j = i + 1 ; j < cluster->numOfProducts ; j++)
        {
            if(products[i].clusterId == cluster->id && products[j].clusterId == cluster->id)
            {
                diff = difference(&products[i], &products[j],0);

                if (diff > diameter)
                {
                    diameter = diff;
                }
            }
        }
    }

    cluster->diameter = sqrt(diameter);
}

void cleanClusters(Cluster* allClusters, int numOfClusters)
{
    int i;

    for(i = 0 ; i < numOfClusters ; i++)
    {
        allClusters[i].numOfProducts = 0;
        allClusters[i].diameter = INT_MAX;
    }
}


void freeAllClusters(Cluster* clusters, int numOfClusters)
{
    int i;

    for(i = 0 ; i < numOfClusters ; i++)
    {
        free(clusters[i].centerTransactions);
	}

    free(clusters);
}
