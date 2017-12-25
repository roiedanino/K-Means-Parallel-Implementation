//
// Created by Roie Danino on 26/10/2017.
//

#include "ProductDatabase.h"
#include "Product.h"


void initProductDatabase(ProductDatabase* productDatabase, Product* products, int numOfProducts, int dimension, MPI_Datatype ProductType)
{
    int i;
    productDatabase->numOfProducts = numOfProducts;
    productDatabase->physicalSize = numOfProducts;
    productDatabase->products = products;
    productDatabase->dimension = dimension;
    productDatabase->ProductType = ProductType;

    for(i = 0 ; i < productDatabase->numOfProducts ; i++)
    {
        productDatabase->products[i].index = i;
    }
}

void addProduct(ProductDatabase* productDatabase, Cluster* cluster,Product* product)
{
    product->index = productDatabase->numOfProducts;

    if(productDatabase->numOfProducts >= productDatabase->physicalSize -1)
    {
        productDatabase->physicalSize *= 2;
        productDatabase->products = (Product*)realloc(productDatabase->products,sizeof(Product)*productDatabase->physicalSize);
        checkAllocation(productDatabase->products);
    }

    if(cluster)
    {
        cluster->numOfProducts++;
    }

    productDatabase->products[productDatabase->numOfProducts] = *product;
    productDatabase->products[productDatabase->numOfProducts].index = productDatabase->numOfProducts;

    if(cluster)
    {
        productDatabase->products[productDatabase->numOfProducts].clusterId = cluster->id;
        productDatabase->products[productDatabase->numOfProducts].newClusterId = NO_CLUSTER;
    }
    productDatabase->numOfProducts++;
}

void removeProduct(ProductDatabase* productDatabase, Cluster* cluster, int index)
{

    if(index >= 0 && index < productDatabase->numOfProducts)
    {
        if(productDatabase->products[index].clusterId == cluster->id)
            cluster->numOfProducts--;

        productDatabase->products[index].index = NO_CLUSTER;
        productDatabase->products[index].clusterId = NO_CLUSTER;
        productDatabase->products[index] = productDatabase->products[--(productDatabase->numOfProducts)];
        productDatabase->products[index].index = index;
    }else
    {
            printf("Invalid Index: %d\n", index);
    }
}
