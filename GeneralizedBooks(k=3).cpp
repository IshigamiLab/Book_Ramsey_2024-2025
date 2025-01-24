#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <string>
#include <csignal>
#include <cstdlib>

using namespace std;

void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int elem : row) {
            std::cout << elem << " ";
        }
        std::cout << endl;
    }
}

void printArray(const vector<int>& arr) {
    for (const auto& elem : arr){
      std::cout << elem << " ";
    }
    std::cout << endl;
}

int main() {
    
    int n = 38;
    if (n <= 1) {
        std::cout << "Size of the matrix should be greater than 1." << endl;
        return 1;
    }

    int a,b;
    a = 5;
    b = 5;

    int zeros1 = (n/2-1)/4;
    int ones1 = (n/2-1)/2 - zeros1;
    int zeros = (n/2-1)/2+1;
    int ones = n/2 - zeros;

    vector<vector<int>> matrix(n, vector<int>(n));

    auto start = std::chrono::high_resolution_clock::now();

    vector<int> D(n/2);

    vector<vector<int>> matrix1(n/2, vector<int>(n/2));
    vector<vector<int>> matrix2(n/2, vector<int>(n/2));
    vector<vector<int>> matrix3(n/2, vector<int>(n/2));

    unsigned long long totalArrays = pow(2, (n/2-1)/2-1);
    unsigned long long totalArrays2 = pow(2,n/2);

    unsigned long long lb = 0;
    unsigned long long ub = totalArrays2;

    std::cout << "idd runs from " << lb << " to " << ub << ", " << ub - lb << "times" << endl;

    for (unsigned long long id = 1; id < totalArrays; id++) 
    {

        int det_counter1=0;
        D[0] = 0;

        for (int j = 1; j <= (n/2-1)/2; ++j) {
            D[j] = (id >> (j - 1)) & 1;
            D[n/2-j] = D[j];
	          det_counter1+=D[j];
        }

	      if( det_counter1==((n/2-1)/4) ){

            for (int k = 0; k < n/2; ++k) {
                for (int l = 0; l < n/2; ++l) {
                    matrix1[k][l] = D[(l - k + n/2) % (n/2)];
                }
            }

            for (unsigned long long idd=lb;idd<ub;idd++){
	              int det_counter = 0;
	              for (int j = 0; j < n/2; j++) {
	                  D[j] = (idd & (1 << j)) ? 1 : 0;
		                det_counter+=D[j];
                }

	              if(det_counter==(n/4)){
	      
                for(int i = 0; i < n/2; i++) {
                    for(int j = 0; j < n/2; j++){
                        matrix2[i][j] = D[(j - i + n/2) % (n/2)];
                    }
                }

		            for (unsigned long long iddd=1;iddd<totalArrays;iddd++){
		                int det_counter2=0;
		                D[0] = 0;

		                for(int j=1;j<=(n/2-1)/2;j++){
		                    D[j] = (iddd>>(j-1)) & 1;
		                    D[n/2-j] = D[j];
		                    det_counter2+=D[j];
		                }

		                if(det_counter2==((n/2-1)/4+1)){
		                    for(int k=0;k<n/2;k++){
		                        for(int l=0;l<n/2;l++){
			                          matrix3[k][l] = D[(l-k+n/2)%(n/2)];
		                        }
		                    }

                        for (int i = 0; i < n/2; i++) {
                            for(int j = 0; j < n/2; j++){
                                matrix[i][j] = matrix1[i][j];
                                matrix[i+n/2][j] = matrix2[i][j];
                                matrix[i][j+n/2] = matrix2[j][i];
                                matrix[i+n/2][j+n/2] = matrix3[i][j];
                            }
                        }

                        int flag1 = 0;
                        int flag0 = 0;

                        int i,j,k;
		                    int idx;
                        int counter1 = 0;
                        int counter0 = 0;

		                    #pragma omp parallel for
                        for(i=0;i<n-2;i++){
                            for(j=i+1;j<n-1;j++){
                                for(k=j+1;k<n;k++){
			                              if( ((matrix[i][j]*matrix[j][k]*matrix[i][k])==1) ){
			                                  for(idx=0;idx<n;idx++){
			                                      if( (i!=idx) && (j!=idx) && (k!=idx) ){
				                                        if( matrix[i][idx]*matrix[j][idx]*matrix[k][idx] ){
				                                            counter1++;
				                                            if(counter1==a){
				                                                flag1 = 1;
				                                            }
				                                        } 
			                                      }
			                                  }
			                                  counter1=0;
			                                  counter0=0;
			                              } else if( ((matrix[i][j]+matrix[j][k]+matrix[i][k])==0) ){
                                        for(idx=0;idx<n;idx++){
			                                      if( (i!=idx) && (j!=idx) && (k!=idx) ){
				                                        if( (matrix[i][idx]+matrix[j][idx]+matrix[k][idx])==0 ){
				                                            counter0++;
				                                            if(counter0==b){
				                                                flag0 = 1;
				                                            }
				                                        }
			                                      }                                                            
			                                  }
                                        counter0=0;
                                        counter1=0;
                                    }
			                          }
		                        }
		                    }


                        if((flag0 == 0) && (flag1 == 0)){
                            std::cout << "Graph not found as below" << std::endl;
                            printMatrix(matrix);
                            std::cout << "_________________________________________________" << std::endl;
                            std::cout << "R(" << a << "," << b << ") bigger than " << n << std::endl;
		                        std::cout << "ID is " << id << ", and idd is " << idd << endl;
                            auto end = std::chrono::high_resolution_clock::now();
                            std::chrono::duration<double, std::milli> elapsed = end - start;
                            std::cout << "Time: " << elapsed.count() << " ms" << std::endl;
		                        std::string filename = "results_K3_GB";
		                        filename += std::to_string(a);
		                        filename += std::to_string(b);
		                        filename += "_N";
		                        filename += std::to_string(n);
		                        filename += ".txt";
		    
		                        std::ofstream outFile(filename);
		                        if(outFile.is_open()){
		                            for(int i=0; i<n; i++){
			                              for(int j=0; j<n; j++){
			                                  outFile << matrix[i][j] << " ";
			                              }
		                            }
		                            outFile.close();
		                        }
		                        string mail_command = "echo 'The process R(" + to_string(a) + "," + to_string(b) + ") is done!' | mail -s 'Done!' -A "+ filename + " " + "Your mail address";
		                        system(mail_command.c_str());
                            return 0;   
                        }
	                    }
	                  }
	                }
	          }
        }
	      std::cout << "ID " << id << " passed!" << endl;
  }

  std::cout << " Cannot find the matrix that prove that R(" << a << "," << b << ") is bigger than " << n << std::endl;    
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;
  std::cout << "Time: " << elapsed.count() << " ms" << std::endl;
  return 0;
}
