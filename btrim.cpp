//~ de Bruijn graph TRIMing tool
//~ Copyright (C) 2017  Limasset Antoine and Marchet Camille

//~ This program is free software: you can redistribute it and/or modify
//~ it under the terms of the GNU Affero General Public License as
//~ published by the Free Software Foundation, either version 3 of the
//~ License, or (at your option) any later version.

//~ This program is distributed in the hope that it will be useful,
//~ but WITHOUT ANY WARRANTY; without even the implied warranty of
//~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//~ GNU Affero General Public License for more details.

//~ You should have received a copy of the GNU Affero General Public License
//~ along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fstream>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>
#include <chrono>
#include <unistd.h>





using namespace std;
using namespace chrono;




string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



uint64_t str2num(const string& str){
	uint64_t res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}



vector<bool> str2bool(const string& str){
	vector<bool> result;
	for(uint i(0);i<str.size();i++){
		switch (str[i]){
			case 'A':result.push_back(false);result.push_back(false);break;
			case 'C':result.push_back(false);result.push_back(true);break;
			case 'G':result.push_back(true);result.push_back(false);break;
			default:result.push_back(true);result.push_back(true);break;
		}
	}
	return result;
}



vector<bool> int2bool(uint n){
	vector<bool> result;
	while(n>0){
		if(n%2==0){
			result.push_back(false);
		}else{
			result.push_back(true);
		}
		n>>=1;
	}
	if(result.size()%2==0){
		result.push_back(false);
	}
	return result;
}



uint bool2int(const vector<bool>& v){
	uint result(0);
	uint acc(1);
	for(uint i(0);i<v.size();++i){
		if(v[i]){
			result+=acc;
		}
		acc<<=1;
	}
	return result;
}



string bool2str(const vector<bool>& v){
	string result;
	for(uint i(0);i<v.size();i+=2){
		if(v[i]){
			if(v[i+1]){
				result.push_back('T');
			}else{
				result.push_back('G');
			}
		}else{
			if(v[i+1]){
				result.push_back('C');
			}else{
				result.push_back('A');
			}
		}
	}
	return result;
}



uint32_t xs(uint32_t y){
	y^=(y<<13); y^=(y>>17); return (y^=(y<<15));
}



char randNucle(char c){
	switch (rand()%4){
		case 0:
			if(c!='A'){
				return 'A';
			}
		case 1:
			if(c!='C'){
				return 'C';
			}
		case 2:
			if(c!='G'){
				return 'G';
			}
		case 3:
			if(c!='T'){
				return 'T';
			}
	}
	return randNucle(c);
}



string compactionNoRecur(const string& seq1,const string& seq2, uint k){
	//~ cout<<seq1<<" "<<seq2<<endl;
	uint s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string rc2(revComp(seq2)),end1(seq1.substr(s1-k,k)), beg2(seq2.substr(0,k));
	if(end1==beg2){return seq1+(seq2.substr(k));}
	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){return seq1+(rc2.substr(k));}
	//~ return compaction(revComp(seq1),seq2,k);
	cout<<"fail"<<endl;
	return "";
}


string compaction(const string& seq1,const string& seq2, uint k){
	//~ cout<<seq1<<" "<<seq2<<endl;
	uint s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string rc2(revComp(seq2)),end1(seq1.substr(s1-k,k)), beg2(seq2.substr(0,k));
	if(end1==beg2){return seq1+(seq2.substr(k));}
	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){return seq1+(rc2.substr(k));}
	return compactionNoRecur(revComp(seq1),seq2,k);
}



uint getPosition(const vector<vector<bool>>& unitigs,uint n){
	//~ cout<<"getPOsition"<<n<<endl;
	auto v(unitigs[n]);
	if(v.size()%2==1){
		uint newPos(bool2int(v));
		return getPosition(unitigs,newPos);
	}else{
		return n;
	}
	return 0;
}



void str2bin(const string& str, char* cstr){
	uint j(0),mult(1);
	for(uint i(0); i<str.size(); ++i){
		switch (str[i]){
			//~ case 'A':cstr[j]+=0;break;
			case 'C':cstr[j]+=1*mult;break;
			case 'G':cstr[j]+=2*mult;break;
			case 'T':cstr[j]+=3*mult;break;
		}
		mult<<=2;
		if(mult>64){
			++j;
			mult=1;
		}
	}
}



vector<pair<string,uint32_t>> uniqueOnly(vector<pair<string,uint32_t>>& v){
	bool unique(true);
	vector<pair<string,uint32_t>> result;
	for(uint i(0);i+1<v.size();++i){
		if(v[i].first==v[i+1].first){
			unique=false;
		}else{
			if(unique){
				result.push_back(v[i]);
			}
			unique=true;
		}
	}
	if(unique and v.size()>0){
		result.push_back(v[v.size()-1]);
	}
	return result;
}



double parseCoverage(const string& str){
	size_t pos(str.find("km:f:"));
	if(pos==string::npos){
		pos=(str.find("KM:f:"));
	}
	uint i(1);
	while(str[i+pos+5]!=' '){
		++i;
	}
	return stof(str.substr(pos+5,i));
}

void usage(){
	cout<<"Usage:"<<endl;
		cout<<"-u [Unitig file]\n"
		<<"-k [Kmer size]\n"
		<<"-t [Tipping length (none)]\n"
		<<"-T [Cleaning Step (1)]\n"
		<<"-c [Core used (1)]\n"
		<<"-h [Hash size, use 2^h files (8 for 256 files)]\n"
		<<"-f [Unitig min coverage (none)]\n"
		<<"-a [Edge filtering ratio (none)]\n"
		<<"-o [Output file (out_tipped)]\n"
		<<endl;
}


void cleaning(string outFile, string inputUnitig,int nbFiles,int tipingSize,int coreUsed,int unitigThreshold,int ratioCoverage,int kmerSize){
	auto start=system_clock::now();
	uint64_t tiping(0),compactions(0),unitigFiltered(0);
	ofstream out(outFile);

	ifstream inUnitigs(inputUnitig);
	vector<fstream> beginFiles(nbFiles),endFiles(nbFiles);
	for(uint i(0); i< nbFiles; ++i){
		beginFiles[i].open(".begin"+to_string(i),fstream::out|fstream::in|fstream::binary|fstream::trunc);
		endFiles[i].open(".end"+to_string(i),fstream::out|fstream::in|fstream::binary|fstream::trunc);
	}
	//THIS MEMORY USAGE CAN BE REMOVED
	vector<vector<bool>> unitigs;
	vector<uint> coverages;



	cout<<"\tParitioning"<<endl;
	//PARTITIONING
	string begin,end,unitig,useless,beginRc,endRc;
	uint64_t hashBegin, hashBeginRc, hashEnd, hashEndRc;
	uint32_t unitigIndice(0);
	while(not inUnitigs.eof()){
		getline(inUnitigs,useless);
		unitig="";
		getline(inUnitigs,unitig);
		if(unitig.size()<kmerSize){
			continue;
		}
		uint coverage=parseCoverage(useless);
		//~ if(unitigThreshold>1 and coverage<unitigThreshold){
			//~ unitigFiltered++;
			//~ continue;
		//~ }
		unitigs.push_back(str2bool(unitig));
		//~ cout<<coverage<<endl;
		coverages.push_back(coverage);


		begin=unitig.substr(0,kmerSize);
		beginRc=revComp(begin);
		hashBegin=str2num(begin);
		hashBeginRc=str2num(beginRc);
		if(hashBegin<=hashBeginRc){
			beginFiles[xs(hashBegin)%nbFiles]<<begin;
			beginFiles[xs(hashBegin)%nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}else{
			endFiles[xs(hashBeginRc)%nbFiles]<<beginRc;
			endFiles[xs(hashBeginRc)%nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}

		end=unitig.substr(unitig.size()-kmerSize,kmerSize);
		endRc=revComp(end);
		hashEnd=str2num(end);
		hashEndRc=str2num(endRc);

		if(hashEnd<=hashEndRc){
			endFiles[xs(hashEnd)%nbFiles]<<end;
			endFiles[xs(hashEnd)%nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}else{
			beginFiles[xs(hashEndRc)%nbFiles]<<endRc;
			beginFiles[xs(hashEndRc)%nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}
		++unitigIndice;
	}


	cout<<"\tTipping"<<endl;
	//TIPPING
	//FOREACH FILE
	#pragma omp parallel for num_threads(coreUsed)
	for(uint i=0; i< nbFiles; ++i){
		string content,wordSeq,wordInt,seqBegin,seqEnd;
		vector<pair<string,uint32_t>> beginVector,endVector;
		uint32_t numberRead;
		uint sizeFileBegin(beginFiles[i].tellg());
		beginFiles[i].seekg(0,ios::beg);
		content.resize(sizeFileBegin);
		beginFiles[i].read(&content[0], content.size());
		for(uint ii(0); ii<content.size(); ){
			wordSeq=content.substr(ii,kmerSize);
			wordInt=content.substr(ii+kmerSize,4);
			memcpy(&numberRead,wordInt.c_str(),4);
			beginVector.push_back({wordSeq,numberRead});
			ii+=kmerSize+4;
		}

		uint sizeFileEnd(endFiles[i].tellg());
		endFiles[i].seekg(0,ios::beg);
		content.resize(sizeFileEnd);
		endFiles[i].read(&content[0], content.size());
		for(uint ii(0); ii<content.size(); ){
			wordSeq=content.substr(ii,kmerSize);
			wordInt=content.substr(ii+kmerSize,4);
			memcpy(&numberRead,wordInt.c_str(),4);
			endVector.push_back({wordSeq,numberRead});
			ii+=kmerSize+4;
		}
		//~ cout<<"gosort"<<endl;
		sort(beginVector.begin(), beginVector.end());
		sort(endVector.begin(), endVector.end());
		uint indiceBegin(0), indiceEnd(0);
		vector<pair<uint,uint>> coverageComparison;
		while(indiceBegin<beginVector.size() and indiceEnd<endVector.size()){
			//~ cout<<"golopp"<<endl;
			seqBegin=beginVector[indiceBegin].first;
			seqEnd=endVector[indiceEnd].first;
			if(seqBegin==seqEnd){
				coverageComparison={};
				//if two begin and one is ten time larger remove the lower
				while(seqBegin==beginVector[indiceBegin].first){
					if(ratioCoverage>1){
						coverageComparison.push_back({coverages[beginVector[indiceBegin].second],beginVector[indiceBegin].second});
					}
					++indiceBegin;
					//~ cout<<"go1"<<endl;
					if(indiceBegin==beginVector.size()){
						break;
					}
					//~ cout<<coverages.size()<<" "<< beginVector[indiceBegin].second <<endl;

					//~ cout<<"go1.6"<<endl;
				}

				if(coverageComparison.size()>1){
					//~ cout<<"go2"<<endl;
					sort(coverageComparison.begin(),coverageComparison.end());
					//~ cout<<coverageComparison[0].first<<" "<<coverageComparison[1].first<<endl;
					if(ratioCoverage*coverageComparison[0].first<coverageComparison[1].first){
						//~ cout<<"go3"<<endl;
						for(uint iii(1);iii<coverageComparison.size();++iii){
							#pragma omp critical(dataupdate)
							{
								unitigs[coverageComparison[iii].second]={};
							}
						}
					}
				}

				coverageComparison={};
				while(seqEnd==endVector[indiceEnd].first){
					if(ratioCoverage>1){
						coverageComparison.push_back({coverages[endVector[indiceEnd].second],endVector[indiceEnd].second});
					}
					++indiceEnd;
					if(indiceEnd==endVector.size()){
						break;
					}

					//~ cout<<"go1.5"<<endl;
				}

				if(coverageComparison.size()>1){
					//~ cout<<"go2"<<endl;
					sort(coverageComparison.begin(),coverageComparison.end());
					//~ cout<<coverageComparison[0].first<<" "<<coverageComparison[1].first<<endl;
					if(ratioCoverage*coverageComparison[0].first<coverageComparison[1].first){
						//~ cout<<"go3"<<endl;
						for(uint iii(1);iii<coverageComparison.size();++iii){
							#pragma omp critical(dataupdate)
							{
								unitigs[coverageComparison[iii].second]={};
							}
						}
					}
				}
				continue;
			}
			//~ cout<<"golop2"<<endl;
			if(seqBegin<seqEnd){
				while(seqBegin==beginVector[indiceBegin].first){
					if(unitigs[beginVector[indiceBegin].second].size()<tipingSize and unitigs[beginVector[indiceBegin].second].size()>0){
						#pragma omp critical(dataupdate)
						{
							++tiping;
							unitigs[beginVector[indiceBegin].second]={};
						}
					}
					++indiceBegin;
					if(indiceBegin==beginVector.size()){
						break;
					}
				}
				continue;
			}
			//~ cout<<"golopp3"<<endl;
			if(seqBegin>seqEnd){
				while(seqEnd==endVector[indiceEnd].first){
					if(unitigs[endVector[indiceEnd].second].size()<tipingSize and unitigs[endVector[indiceEnd].second].size()>0){
						#pragma omp critical(dataupdate)
						{
							++tiping;
							unitigs[endVector[indiceEnd].second]={};
						}
					}
					++indiceEnd;
					if(indiceEnd==endVector.size()){
						break;
					}
				}
				continue;
			}
		}
	}



	cout<<"\tRecompaction"<<endl;
	//RECOMPACTION
	#pragma omp parallel for num_threads(1)
	for(uint i=0; i< nbFiles; ++i){
		string content,wordSeq,wordInt,seqBegin,seqEnd,str1,str2,compacted;
		vector<pair<string,uint32_t>> beginVector,endVector;
		uint numberRead;
		uint sizeFileBegin(beginFiles[i].tellg());
		beginFiles[i].seekg(0,ios::beg);
		content.resize(sizeFileBegin);
		beginFiles[i].read(&content[0], content.size());
		for(uint ii(0); ii<content.size(); ){
			wordSeq=content.substr(ii,kmerSize);
			wordInt=content.substr(ii+kmerSize,4);
			memcpy(&numberRead,wordInt.c_str(),4);
			if(not unitigs[numberRead].empty()){
				beginVector.push_back({wordSeq,numberRead});
			}else{
			}
			ii+=kmerSize+4;
		}
		uint sizeFileEnd(endFiles[i].tellg());
		endFiles[i].seekg(0,ios::beg);
		content.resize(sizeFileEnd);
		endFiles[i].read(&content[0], content.size());
		for(uint ii(0); ii<content.size(); ){
			wordSeq=content.substr(ii,kmerSize);
			wordInt=content.substr(ii+kmerSize,4);
			memcpy(&numberRead,wordInt.c_str(),4);
			if(not unitigs[numberRead].empty()){
				endVector.push_back({wordSeq,numberRead});
			}
			ii+=kmerSize+4;
		}
		sort(beginVector.begin(), beginVector.end());
		beginVector=uniqueOnly(beginVector);

		sort(endVector.begin(), endVector.end());
		endVector=uniqueOnly(endVector);
		uint indiceBegin(0), indiceEnd(0);
		while(indiceBegin<beginVector.size() and indiceEnd<endVector.size()){
			seqBegin=beginVector[indiceBegin].first;
			seqEnd=endVector[indiceEnd].first;
			if(seqBegin==seqEnd){
				#pragma omp critical(dataupdate)
				{
					uint position1(getPosition(unitigs,beginVector[indiceBegin].second));
					uint position2(getPosition(unitigs,endVector[indiceEnd].second));
					if(position1!=position2){
						str1=(bool2str(unitigs[position1]));
						str2=(bool2str(unitigs[position2]));
						compacted=(compaction(str1,str2,kmerSize));
						unitigs[position1]=str2bool(compacted);
						coverages[position1]=(coverages[position1]*str1.size()+coverages[position2]*str2.size())/compacted.size();
						unitigs[position2]=int2bool(position1);
						compactions++;
					}
				}
				++indiceBegin;
				++indiceEnd;
				continue;
			}
			if(seqBegin<seqEnd){
				++indiceBegin;
				continue;
			}
			if(seqBegin>seqEnd){
				++indiceEnd;
				continue;
			}
		}
		remove((".begin"+to_string(i)).c_str());
		remove((".end"+to_string(i)).c_str());
	}

	//OUTPUT
	for(uint i(0); i<unitigs.size(); ++i){
		if((not unitigs[i].empty()) and (unitigs[i].size()%2==0) and coverages[i]>=unitigThreshold){
			out<<">km:f:"<<coverages[i]<<"\n";
			out<<bool2str(unitigs[i])<<'\n';
		}else{
			unitigFiltered++;
		}
	}

	if(unitigThreshold>1){
		cout<<"\tUnitig filtered: "+intToString(unitigFiltered)<<endl;
	}
	cout<<"\tTips removed: "+intToString(tiping)<<endl;
	cout<<"\tUnitigs compacted: "+intToString(compactions)<<endl;
	auto endTime=system_clock::now();
    auto waitedFor=endTime-start;
    cout<<"\nCleaned in "<<duration_cast<seconds>(waitedFor).count()<<" seconds"<<endl;
}



//TODO multiple functions
//ONE FUNCTION TO RULE THEM ALL
int main(int argc, char ** argv){
	//INIT
	if(argc<3){
		usage();
		exit(0);
	}
	int unitigThreshold(1);
	string inputUnitig;
	//~ ifstream inUnitigs(inputUnitig);
	uint kmerSize;
	uint tipingSize(1);
	uint tipingStep(1);
	uint coreUsed(1);
	uint hashSize(8);//256 FILES
	uint nbFiles(1<<(hashSize-1));
	uint ratioCoverage(1);
	string outFile("out_tipped");
	char c;
	while ((c = getopt (argc, argv, "u:k:t:c:h:f:a:o:T:")) != -1){
		switch(c){
		case 'u':
			inputUnitig=optarg;
			break;
		case 'k':
			kmerSize=stoi(optarg);
			--kmerSize;
			break;
		case 't':
			tipingSize=2*stoi(optarg);
			break;
		case 'T':
			tipingStep=stoi(optarg);
			break;
		case 'c':
			coreUsed=stoi(optarg);
			break;
		case 'h':
			hashSize=stoi(optarg);
			break;
		case 'f':
			unitigThreshold=stoi(optarg);
			break;
		case 'o':
			outFile=optarg;
			break;
		case 'a':
			ratioCoverage=stoi(optarg);
			break;
		}
	}
	for(uint i(1);i<=tipingStep;++i){
		cout<<"Strep "<<i<<endl;
		if(i==tipingStep){
			cleaning( outFile,  inputUnitig, nbFiles, tipingSize, coreUsed, unitigThreshold, ratioCoverage, kmerSize);
		}else{
			cleaning( outFile+to_string(i),  inputUnitig, nbFiles, tipingSize, coreUsed, unitigThreshold, ratioCoverage, kmerSize);
			inputUnitig=outFile+to_string(i);
		}
	}

}
