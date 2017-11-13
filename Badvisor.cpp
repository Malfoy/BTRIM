#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>



using namespace std;



void help(){
	cout<<"./badvisor reads.fa"<<endl;
}



void printInt(const int value){
	string numWithCommas = to_string(value);
	int insertPosition = numWithCommas.length() - 3;
	while (insertPosition > 0) {
		numWithCommas.insert(insertPosition, ",");
		insertPosition-=3;
	}
	cout<<numWithCommas<<endl;
}



int main(int argc, char *argv[]) {
	if(argc<3){
		help();
		exit(0);
	}
	string readsFile(argv[1]),line,lineF0,lineF1;
	uint64_t coverageAsked(stol(argv[2]));
	vector<vector<uint64_t>> histograms;
	vector<uint64_t> numberKmerDistinct,minimumList;
	double frac(1.3);
	if(argc>3){
		frac=stod(argv[3]);
	}
	for(uint k(21);k<201;k+=10){
		ifstream stream(readsFile+"_k"+to_string(k)+".hist");
		if(not stream.is_open()){
			//~ cout<<"no file "+readsFile+"_k"+to_string(k)+".hist"<<endl;
			break;
		}
		getline(stream,lineF1,'	');
		getline(stream,lineF1);
		uint64_t F1(stol(lineF1));
		if(F1==0){
			break;
		}
		getline(stream,lineF0,'	');
		lineF0="";
		getline(stream,lineF0);
		//~ uint64_t numberKmer(stoi(lineF1));
		numberKmerDistinct.push_back(stol(lineF0));
		vector<uint64_t> abundances;
		while(not stream.eof()){
			getline(stream,line,'	');
			getline(stream,line);
			if(line.size()>0){
				abundances.push_back(stol(line));
			}
		}
		histograms.push_back(abundances);
	}
	bool cont(true);
	for(uint i(0);i<histograms.size() and cont;++i){
		for(uint ii(0);ii<histograms[i].size() and cont;++ii){
			if(histograms[i][ii]<histograms[i][ii+1]*frac){
				if(ii>=coverageAsked){
					minimumList.push_back(ii);

				}else{
					cont=false;
				}
				break;
			}
		}
	}
	cout<<minimumList.size()*10+21;
    return 0;
}
