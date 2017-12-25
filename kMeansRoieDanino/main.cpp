#include <iostream>
#include <string.h>

using namespace std;
int main() 
{
    const int numprocs = 7, maxK = 10, bufferSize = maxK / numprocs + 1;
    int ranksCounters[] = {0,0,0,0,0,0,0};
    int data[numprocs*bufferSize];
    int k = 2;
    char key[] = "BIRD";

    //int newNumOfCluster = numOfCluster / numprocs;
    //newNumOfCluster += (numOfCluster - myid  > 1 ? 1 : 0);


    for(k = 2 ; k <= maxK ; k++)
    {
        cout << "K = " << k << "   K / numprocs + 1= " << k/numprocs + 1 << endl;

        for(int i = 0 ; i < numprocs ; i++)
        {
            ranksCounters[i] = k / numprocs;

            if(k%numprocs > i)
                ranksCounters[i]++;

            cout << "Rank = " << i << " num of clusters: " << ranksCounters[i] << endl;
        }
        cout << endl;
    }

   /* for(int i = 0 ; i < numprocs ; i++)
    {
        ranksCounters[i] = k / numprocs ;

        if(!(k%numprocs < i + 1 || (i + 1 == numprocs && k%numprocs == 0)))
        {
            ranksCounters[i] += 1;

        }
        cout << "Rank = " << i << " num of clusters: " << ranksCounters[i] << endl;
    }*/

    for(int i = 0 ; i < numprocs ; i++ )
    {
        //if(k%numprocs <= i + 1)//if(k%numprocs <= i + 1 || (i + 1 == numprocs && k%numprocs == 0))
       // {
         //   cout << "Rank = " << i << " num of clusters: " << ranksCounters[i] << endl;
        //}
        for(int j = 0 ; j <  bufferSize; j++)
        {
            if(ranksCounters[i] > j)
            {
                data[i*bufferSize + j] = i;
            }
            else
            {
                data[i*bufferSize + j] = -1;
            }
            printf("DATA i:%d j:%d == %d\n",i,j, data[i*bufferSize + j]);
        }
    }

    int size = bufferSize*numprocs;
    for(int i = size - 1 - bufferSize; i >= 0 ; i -= 1)
    {
        if(data[i] == -1 )
        {
            if(data[size - 1] == -1)
                size--;
            data[i] = data[size - 1];
            size--;
        }
      //  *data[i + bufferSize - 1] = -2;
        //data[5*bufferSize + 1] = data[6*bufferSize];
        /*
         * [0,0,0,0,1,1,1,X,2,2,2,X,3,3,3,X]
         * [0,0,0,0,1,1,1,X,2,2,2,3,3,3,_,X]
         * [0,0,0,0,1,1,1,3,2,2,2,3,3,_,_,X]
         */
    }
    printf("After Fix:\n");
    for(int i = 0 ; i < bufferSize*numprocs ; i++ )
    {
        printf("DATA[%d]  == %d\n",i, data[i]);
    }
    printf("First K:\n");
    for(int i = 0 ; i < k ; i++ )
    {
        printf("DATA[%d]  == %d\n",i, data[i]);
    }

    return 0;
}