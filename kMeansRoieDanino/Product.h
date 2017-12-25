//
// Created by Roie Danino on 15/09/2017.
//

#ifndef __PRODUCT_H
#define __PRODUCT_H
#include "Product.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define LEN 5

#define NO_CLUSTER -1


#define MAX_NUM_OF_PRODUCTS 5000

struct cluster_t;

typedef struct
{
    char productName[LEN];

    double* transactions;

    int dimension;

    int clusterId;

    int newClusterId;

    int index;

    int rank;
}Product;

typedef struct info_t
{
    int numOfProducts;

    int dimension;

    int maxK;

    int limitOfIterations;

    double qualityMeasure;

}Info;

Product* readProductsFromFile(Info* info);

double difference(const Product* p1, const Product* p2, int toSqrt = 1);


double* getProductsTransactions(Product* products, int numOfProducts);

void setProductsTransactions(Product* products, int numOfProducts, double* transactions, int dimension);

void cleanProduct(Product* product);

void freeAllProducts(Product* allProducts,int numOfProducts);



#endif //__PRODUCT_H
