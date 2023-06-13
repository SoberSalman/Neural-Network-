#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <semaphore.h>
#include <sys/wait.h>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

double** matMul(double** arr1, int arr1Rows, int arr1Cols, double** arr2, int arr2Rows, int arr2Cols) {

    double** result = new double*[arr1Rows];
    for (int i = 0; i < arr1Rows; ++i) {
        result[i] = new double[arr2Cols];
        for (int j = 0; j < arr2Cols; ++j) {
            result[i][j] = 0;
        }
    }

    for (int i = 0; i < arr1Rows; ++i) {
        for (int j = 0; j < arr2Cols; ++j) {
            for (int k = 0; k < arr1Cols; ++k) {
                result[i][j] += arr1[i][k] * arr2[k][j];
            }
        }
    }

    return result;
}


struct Node{

    int nodeNumber;
    int numWeights;
    double** weights;

    double* doAlgorithm(double** input){

        double** x = matMul(input, 1, numWeights, weights, numWeights, 1);

        double* answer = new double(x[0][0]);

        return answer;

    }

};

struct BarrierArgs 
{
    pthread_barrier_t* barrier;
};

struct FunkyFuncArgs {
    Node* node;
    pthread_mutex_t* mutex;
    double** inputData;
    double* outputArray;  // instead of returning result through the pointer , We could try writing the result directly inside an output array
    int outputIndex;
};

void* funkyFunc(void* args) 
{
    FunkyFuncArgs* funcArgs = (FunkyFuncArgs*)args;
    Node* n = funcArgs->node;
    pthread_mutex_t mutex = *(funcArgs->mutex);
    double** input = funcArgs->inputData;

    pthread_mutex_lock(&mutex);

    double* output = n->doAlgorithm(input); // i assume this is some algo doing some of your calculations and returning the result

    pthread_mutex_unlock(&mutex);
    
    funcArgs->outputArray[funcArgs->outputIndex] = *output; // writing the result in the output array,

    

    pthread_exit((void*)output);

}

struct Layer{

    double** layerData;
    int layerNumber;
    Node* Nodes;
    int numberNodes = 0;
    pthread_mutex_t nodeMutex;

    Layer()
    {

        pthread_mutex_init(&nodeMutex, NULL);

    }

    double** callNodes()
    {

        pthread_t tid[numberNodes];
        FunkyFuncArgs args[numberNodes];

        double** outputArr = new double*[1];
        outputArr[0] = new double[numberNodes];

        for (int i = 0; i < numberNodes; i++) 
        {
            args[i].node = &Nodes[i];
            args[i].mutex = &nodeMutex;
            args[i].inputData = layerData;
            args[i].outputArray = outputArr[0]; // Pass the output array to each thread so that instead of joining them threads and passing the array as a pointer to write the result, we just write them in the ThreadArgument
            args[i].outputIndex = i;                 // indexing is important so that we should obviously know where to write
            pthread_create(&tid[i], NULL, funkyFunc, &args[i]);
        }

        for(int i = 0; i < numberNodes; i++)
        {
            usleep(1000); // Passing null here :P , waiting for all threads to finsih 
        }

        return outputArr;

    }

    void backProp(){

        cout << "BackProp begun by layer " << layerNumber << endl;

        pthread_t tid[numberNodes];
        FunkyFuncArgs args[numberNodes];

        for (int i = 0; i < numberNodes; i++) {
            args[i].node = &Nodes[i];
            args[i].mutex = &nodeMutex;
            pthread_create(&tid[i], NULL, funkyFunc, &args[i]);
        }

        for(int i = 0; i < numberNodes; i++){

            pthread_join(tid[i], NULL);

        }

    }

    friend void* funkyFunc(void* argc);

};

int counter = 0;

struct NueralNetwork{

    Layer* Layers;
    int numLayers;

    NueralNetwork(int nL, int* nodesInLayer, vector<double> values, int inputDim, int* weights){

        counter = 0;

        numLayers = nL;

        Layers = new Layer[nL];

        for(int i = 0; i < nL; i++){

            Layers[i].layerNumber = i + 1;
            Layers[i].numberNodes = nodesInLayer[i];
            Layers[i].Nodes = new Node[nodesInLayer[i]];

            for(int j = 0; j < nodesInLayer[i]; j++){

                Layers[i].Nodes[j].nodeNumber = j + 1;
                Layers[i].Nodes[j].weights = new double*[weights[i]];
                Layers[i].Nodes[j].numWeights = weights[i];

                for(int k = 0; k < weights[i]; k++){
                    Layers[i].Nodes[j].weights[k] = new double[1];
                    Layers[i].Nodes[j].weights[k][0] = values[counter++];
                }

            }

        }

        for(int i = 0; i < nL; i++){

            string fifoString = "";
            fifoString += "DataFifo";
            fifoString += to_string(i);
            fifoString += to_string(i + 1);

            mkfifo(fifoString.c_str(), 0666);

        }

    }

    double* stuff(double** input){

        string fifoString = "DataFifo01";

        cout << "Ai stuff begun by nueral network" << endl;

        int value;

        pid_t initialFork = fork();

        if(initialFork){

            int fd = open(fifoString.c_str(), O_WRONLY);

            int rows = 1;
            int cols = 2;

            double flat_array[rows * cols];
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    flat_array[i * cols + j] = input[i][j];
                }
            }

            write(fd, flat_array, sizeof(flat_array));

            close(fd);

            exit(0);

        }

        sleep(1);

        for(int i = 0; i < numLayers; i++){

            pid_t layerProcess = fork();

            if(layerProcess == 0){

                string readData = "DataFifo";
                readData += to_string(i);
                readData += to_string(i + 1);

                string writeData = "DataFifo";
                writeData += to_string(i + 1);
                writeData += to_string(i + 2);

                int received_rows = 1;
                int received_cols = Layers[i].Nodes->numWeights;

                int fd2 = open(readData.c_str(), O_RDONLY);
                double received_flat_array[received_rows * received_cols];
                read(fd2, received_flat_array, sizeof(received_flat_array));
                close(fd2);

                double** received_array_2d = new double*[received_rows];
                
                for(int i = 0; i < received_rows; i++){

                    received_array_2d[i] = new double[received_cols];

                }
            
                for (int i = 0; i < received_rows; ++i) {
                    for (int j = 0; j < received_cols; ++j) {
                        received_array_2d[i][j] = received_flat_array[i * received_cols + j];
                    }
                }

                Layers[i].layerData = received_array_2d;

                double** output = Layers[i].callNodes();

                if(i != numLayers - 1){

                    int fd3 = open(writeData.c_str(), O_WRONLY);

                    int rows = 1;
                    int cols = Layers[i].numberNodes;

                    double flat_array[rows * cols];
                    for (int i = 0; i < rows; ++i) {
                        for (int j = 0; j < cols; ++j) {
                            flat_array[i * cols + j] = output[i][j];
                        }
                    }

                    write(fd3, flat_array, sizeof(flat_array));

                    close(fd3);

                }

                else{

                    double x = output[0][0];
                    double* fx = new double[2];
                    fx[0] = ((x*x) + x + 1)/2.0;
                    fx[1] = ((x*x) - x) /2.0;

                    cout << endl << output[0][0] << " Final answer" << endl;
                    cout << endl << fx[0] << " " << fx[1] << endl;

                    int fd6 = open(readData.c_str(), O_WRONLY);
                    write(fd6, fx, 16);
                    close(fd6);

                    exit(0);

                }

                int fd4 = open(writeData.c_str(), O_RDONLY);
                double* finalfx = new double[2];
                read(fd4, finalfx, 16);
                close(fd4);

                cout << "Layer " << i + 1 << ": " << finalfx[0] << " " << finalfx[1] << endl;

                int fd5 = open(readData.c_str(), O_WRONLY);
                write(fd5, finalfx, 16);

                if(i == 0){

                    int fd10 = open("DataFifo01", O_WRONLY);
                    write(fd10, finalfx, 16);
                    close(fd10);

                }

                exit(0);

            }

        }

        double* finalfinal = new double[2];

        int fd11 = open("DataFifo01", O_RDONLY);
        read(fd11, finalfinal, 16);

        cout << "Output recieved in main : " <<  finalfinal[0] << " " << finalfinal[1] << endl;

        cout << "Nueral Network Pass Completed" << endl;

        return finalfinal;

    }

};

double randomDouble(){

    return (float) rand()/RAND_MAX;
    
}

int main(){

    ifstream file("weights.txt");
    string line;
    std::vector<double> values;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            double value;
            while (ss >> value) {
                values.push_back(value);
                if (ss.peek() == ',') {
                    ss.ignore();
                }
            }
        }
        file.close();
    }

    srand((unsigned)time(NULL));

    int numLayers = 7;
    int* nodesPerLayer = new int[7];
    nodesPerLayer[0] = 8;
    nodesPerLayer[1] = 8;
    nodesPerLayer[2] = 8;
    nodesPerLayer[3] = 8;
    nodesPerLayer[4] = 8;
    nodesPerLayer[5] = 8;
    nodesPerLayer[6] = 1;

    int* weightsPerLayer = new int[7];
    weightsPerLayer[0] = 2;
    weightsPerLayer[1] = 8;
    weightsPerLayer[2] = 8;
    weightsPerLayer[3] = 8;
    weightsPerLayer[4] = 8;
    weightsPerLayer[5] = 8;
    weightsPerLayer[6] = 8;

    srand(time(0));

    NueralNetwork NN(7, nodesPerLayer, values, 2, weightsPerLayer);

    cout << values.size() << " " << counter << endl;
 
    double** input = new double*[1];
    input[0] = new double[values.size() - counter];

    int size = values.size() - counter;

    for(int i = 0; i < size; i++){

        input[0][i] = values[counter++];

    }

    //input[0][0] = 0.1;
    //input[0][1] = 0.2;

    double* output = NN.stuff(input);
    double ** newInput = new double*[1];
    newInput[0] = new double[2];

    newInput[0][0] = output[0];
    newInput[0][1] = output[1];

    cout << endl << endl;

    for(int i = 0; i < NN.numLayers; i++)
    {

        string fifoString = "";
        fifoString += "DataFifo";
        fifoString += to_string(i);
        fifoString += to_string(i + 1);

        unlink(fifoString.c_str());

    }

    NueralNetwork NN2(7, nodesPerLayer, values, 2, weightsPerLayer);

    NN2.stuff(newInput);

    for(int i = 0; i < NN2.numLayers; i++){

        string fifoString = "";
        fifoString += "DataFifo";
        fifoString += to_string(i);
        fifoString += to_string(i + 1);

        unlink(fifoString.c_str());

    }

}