//
// Created by Roie Danino on 15/09/2017.
//

#include "TechnicalFunctions.h"





void checkAllocation(void* ptr)
{
    if(!ptr)
    {
        printf("OUT OF MEMORY!!!!");
		fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
        exit(1);
    }
}

void failedToOpenFile(const char* path)
{
    printf("ERROR! Couldn't open the file from: \"%s\", make sure that the path that defined in \"FILE_NAME\" variable is correct");
    MPI_Finalize();
    exit(1);
}

