//
// Created by Roie Danino on 26/10/2017.
//

#ifndef __PRODUCTDATABASE_H
#define __PRODUCTDATABASE_H

#include "Product.h"
#include "stdlib.h"
#include "Cluster.h"
#include "TechnicalFunctions.h"

typedef struct
{
    Product* products;

    int physicalSize;

    int numOfProducts;

    int dimension;

    MPI_Datatype ProductType;

}ProductDatabase;

void initProductDatabase(ProductDatabase* productDatabase, Product* products, int numOfProducts, int dimension, MPI_Datatype ProductType);

void addProduct(ProductDatabase* productDatabase, Cluster* cluster, Product* product);

void removeProduct(ProductDatabase* productDatabase, Cluster* cluster, int index);

#endif //__PRODUCTDATABASE_H
