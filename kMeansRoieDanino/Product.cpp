//
// Created by Roie Danino on 15/09/2017.
//
#define _CRT_SECURE_NO_WARNINGS
#include "Product.h"
#include "MPI_Functions.h"
#include "TechnicalFunctions.h"
#include "omp.h"

Product* readProductsFromFile(Info* info)
{
    int i,j;
    const char* FILE_NAME = "C:\\Users\\afeka\\Downloads\\Sales_Transactions_Dataset_Weekly.txt";
    FILE* productsFile = fopen(FILE_NAME,"r");

    if(!productsFile)
        failedToOpenFile(FILE_NAME);
    fscanf(productsFile,"%d %d %d %d %lf",&info->numOfProducts, &info->dimension, &info->maxK, &info->limitOfIterations, &info->qualityMeasure);

    Product* products = (Product*)malloc(sizeof(Product)*info->numOfProducts);
    checkAllocation(products);

    for(i = 0 ; i < info->numOfProducts ; i++)
    {
        products[i].clusterId = -1;
        products[i].index = -1;
        products[i].dimension = info->dimension;
        products[i].newClusterId = -1;

        fscanf(productsFile,"%s", products[i].productName);

        products[i].transactions = (double*)malloc(sizeof(double)*(info->dimension));
		checkAllocation(products[i].transactions);

        for(j = 0 ; j < info->dimension ; j++)
        {
            fscanf(productsFile,"%lf", &products[i].transactions[j]);
        }

    }

    fclose(productsFile);
    return products;
}

//Euclidian distance calculations
double difference(const Product* p1, const Product* p2, int toSqrt)
{
    int i;
    double sum = 0, partialSum = 0;

    for(i = 0 ; i < p1->dimension ; i++)
    {
        sum += (p1->transactions[i] - p2->transactions[i])*(p1->transactions[i] - p2->transactions[i]);
    }

    return toSqrt ? sqrt(sum) : sum;
}

double* getProductsTransactions(Product* products, int numOfProducts)
{
    int i, dimension = products[0].dimension;
    double* allTransactions = (double*)malloc(sizeof(double)*dimension*numOfProducts);

    for(i = 0 ; i < numOfProducts ; i++)
    {
        memcpy(&allTransactions[i*dimension],products[i].transactions, sizeof(double)*dimension);
    }

    return allTransactions;
}

void setProductsTransactions(Product* products, int numOfProducts, double* transactions, int dimension)
{
    int i;

    for(i = 0 ; i < numOfProducts; i++)
    {
        products[i].transactions = &transactions[i*dimension];
    }
}

void cleanProduct(Product* product)
{
    product->newClusterId = NO_CLUSTER;
    product->clusterId = NO_CLUSTER;
}

void freeAllProducts(Product* allProducts,int numOfProducts)
{
    int i;

    for(i = 0 ; i < numOfProducts; i++)
    {
        free(allProducts[i].transactions);
    }
}