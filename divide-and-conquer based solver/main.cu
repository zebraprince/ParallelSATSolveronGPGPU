#include "struct_def.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>

using namespace std;

__global__ void Expansion_Kernel(Formula, FormulaBatch, FormulaBatch, int);
__global__ void BCP_Kernel(Formula, FormulaBatch);
__global__ void UnitCls_Kernel(Formula, FormulaBatch);
__global__ void SATChecking_Kernel(Formula, FormulaBatch);

/********************Read CNF from files**************************/
char* readFile(FILE *fileIn){
	
	long lSize;
	char * data;
	size_t result;

	if(fileIn == NULL){fputs("File error", stderr); exit(0);}

	//obtain file size
	int m = fseek(fileIn, 0, SEEK_END);
	lSize = ftell(fileIn);
	rewind(fileIn);

	//allocate memory to contain the whole file
	data = (char *) malloc(sizeof(char) * lSize);
	if (data == NULL){fputs("Memory error", stderr); exit(1);}

	//copy the file into memory
	result = fread(data, 1, lSize, fileIn);

	//if (result != lSize){fputs("Reading error", stderr); exit(2);}

	return data;
}

/*********** skip white space **********/
static void skipWhitespace(char** in) {
    	while ( (**in >= 9 && **in <= 13) || **in == 32 )
        	(*in)++;
}


/**********  skip mark line **********/
static void skipLine(char** in) {
    	for (;;){
        	if (**in == 0) return;
        	if (**in == '\n') {
        	  	(*in)++; return;
        	}
        (*in)++;
    	}
}


static int parseInt(char** in) {
    	
	int val = 0;
    	int neg_flag = 0;

    	skipWhitespace(in);
    	if (**in == '-'){
		neg_flag = 1;
		(*in)++;
	}
    	else if (**in == '+') 
		(*in)++;
    	if (**in < '0' || **in > '9')
    	  	fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", **in), exit(1);
    	while (**in >= '0' && **in <= '9'){
        	val = val*10 + (**in - '0');
        	(*in)++;
	}
    	
	return neg_flag ? -val : val;
}


static int getPara(char** in){

	int val = 0;
	
	while (**in < '0' || **in > '9')
		(*in)++;
	while (**in >= '0' && **in <= '9'){
        	val = val*10 + (**in - '0');
        	(*in)++;
	}
	return val;
}


static void readClause(char **in, Formula &f, Clause *cls, vector<int> &parsedLit, int currentCls){
    	
	int tmp;
	int counter=0;

    	for(;;){
        	tmp = parseInt(in);
        	if (tmp == 0) {
			cls[currentCls].countUnassignment = counter;
			cls[currentCls].clsState = Unknown;
			cls[currentCls].clsSize = counter;
			if(currentCls < f.numofCls-1)
				f.numofLitBeforeCls[currentCls+1] = f.numofLitBeforeCls[currentCls] + counter;
			counter=0;
			break;
		}
		parsedLit.push_back(tmp);
		counter++;
    	}
}


static void parse_DIMACS(FILE * in, Formula &f){
	
	char * text = readFile(in);
	int clsCounter = 0;
	vector<int> parsedLit;
	int numofCls;
	
    	for(;;){
        	skipWhitespace(&text);
        	if (*text == 0)
            		break;
		else if (*text == 'p'){
			f.numofLit = getPara(&text);
			f.numofCls = getPara(&text);
			numofCls = f.numofCls;
			f.cls = (Clause *)malloc(sizeof(Clause)*numofCls);			
			f.expPath = (int *)malloc(sizeof(int)*batchSize);
			f.numofLitBeforeCls = (int *)malloc(sizeof(int) * f.numofCls);
			f.numofLitBeforeCls[0] = 0;
			f.formulaState = Unknown;
			skipLine(&text);
		}
        	else if (*text == 'c')
            		skipLine(&text);
        	else{
			if(clsCounter<f.numofCls){
            			readClause(&text, f, f.cls, parsedLit, clsCounter);
				clsCounter++;
			}
		}
    	}
	f.numofTotalLit = f.numofLitBeforeCls[numofCls-1] + f.cls[numofCls-1].clsSize;
	int allLit = f.numofTotalLit;
	f.literal = (int *)malloc(sizeof(int) * allLit);
	f.litState = (int *)malloc(sizeof(int) * allLit);
	for(int i=allLit-1; i>=0; i--){
		f.literal[i] = parsedLit[i];
		f.litState[i] = Unassigned;
		parsedLit.pop_back();
	}	
	//free(text);	
}



void satPrinter(Formula &answer, Formula &origin, float &GPU_time, clock_t &startTime, clock_t &endTime, float &CPU_time){
	
	endTime = clock();
	CPU_time = (double)( (endTime - startTime) / (double)CLOCKS_PER_SEC);
	printf("     GPU running time(ms): %f", GPU_time);
	printf("     CPU running time(s): %f", CPU_time);
}



void runKernel(Formula &origin_d, FormulaBatch &source, FormulaBatch &result, int expPara, int numofCls, int numofTotalLit, int *litState_r, int *formulaState_r, Clause* cls_r, Formula *tmp, float &time){

	size_t malSize;
	malSize = sizeof(Clause) * source.numofFormula * numofCls;
	Clause *cls_h = (Clause *)malloc(malSize);
	for(int i=source.numofFormula -1; i>=0; i--)
		for(int j=numofCls-1; j>=0; j--)
			cls_h[i*numofCls + j] = tmp[i].cls[j];
	cudaMalloc((void **)&source.cls_b, malSize);
	cudaMemcpy(source.cls_b, cls_h, malSize, cudaMemcpyHostToDevice);
	cudaMalloc((void **)&result.cls_b, expPara * malSize);
	free(cls_h);

	malSize = sizeof(int) * source.numofFormula * numofTotalLit;
	int *litState_h = (int *)malloc(malSize);
	for(int i=source.numofFormula -1; i>=0; i--)
		for(int j=numofTotalLit-1; j>=0; j--)
			litState_h[i*numofTotalLit + j] = tmp[i].litState[j];
	cudaMalloc((void **)&source.litState_b, malSize);
	cudaMemcpy(source.litState_b, litState_h, malSize, cudaMemcpyHostToDevice);
	cudaMalloc((void **)&result.litState_b, expPara * malSize);
	free(litState_h);

	malSize = sizeof(int) * source.numofFormula;
	int *formulaState_h = (int *)malloc(malSize);
	for(int i=source.numofFormula -1; i>=0; i--)
		formulaState_h[i] = tmp[i].formulaState;
	cudaMalloc((void **)&source.formulaState_b, malSize);
	cudaMemcpy(source.formulaState_b, formulaState_h, malSize, cudaMemcpyHostToDevice);
	cudaMalloc((void **)&result.formulaState_b, expPara * malSize);
	free(formulaState_h);

	malSize = sizeof(int) * source.numofFormula * batchSize;
	int *expPath_h = (int *)malloc(malSize);
	for(int i=source.numofFormula -1; i>=0; i--)
		for(int j=batchSize-1; j>=0; j--)
			expPath_h[i*batchSize + j] = tmp[i].expPath[j];
	cudaMalloc((void **)&source.expPath_b, malSize);
	cudaMemcpy(source.expPath_b, expPath_h, malSize, cudaMemcpyHostToDevice);
	cudaMalloc((void **)&result.expPath_b, expPara * malSize);
	free(expPath_h);
			
	malSize = sizeof(int) * result.numofFormula * batchSize;
	int *inferPath_h = (int *)malloc(malSize);
	for(int i=result.numofFormula * batchSize -1; i>=0; i--)
		inferPath_h[i] = 0;
	cudaMalloc((void **)&result.inferPath_b, malSize);			
	cudaMemcpy(result.inferPath_b, inferPath_h, malSize, cudaMemcpyHostToDevice);
	free(inferPath_h);

	malSize = sizeof(int) * result.numofFormula;

	int gridX = (int)((numofCls+Thread_per_Block-1)/Thread_per_Block);
	dim3 dimBlock(Thread_per_Block, 1);
	dim3 dimGridEXP(gridX, source.numofFormula);

float time_tmp; 
cudaEvent_t start;
cudaEvent_t stop;

cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start, 0);
	Expansion_Kernel<<<dimGridEXP, dimBlock>>>(origin_d, source, result, expPara);
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
time += time_tmp;

	dim3 dimGridBCP(gridX, result.numofFormula);
	int *notifyUnit_h = (int *)malloc(sizeof(int));
	notifyUnit_h[0] = 0;
	cudaMalloc((void**)&result.notifyUnit, sizeof(int));
	cudaMemcpy(result.notifyUnit, notifyUnit_h, sizeof(int), cudaMemcpyHostToDevice);
	

cudaEventRecord(start, 0);
	UnitCls_Kernel<<<dimGridBCP, dimBlock>>>(origin_d, result);
	cudaThreadSynchronize();
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
time += time_tmp;

	cudaMemcpy(notifyUnit_h, result.notifyUnit, sizeof(int), cudaMemcpyDeviceToHost);

	while(notifyUnit_h[0]){

cudaEventRecord(start, 0);
		BCP_Kernel<<<dimGridBCP, dimBlock>>>(origin_d, result);
		cudaThreadSynchronize();
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
time += time_tmp;

		notifyUnit_h[0] = 0;
		cudaMemcpy(result.notifyUnit, notifyUnit_h, sizeof(int), cudaMemcpyHostToDevice);

cudaEventRecord(start, 0);
		UnitCls_Kernel<<<dimGridBCP, dimBlock>>>(origin_d, result);
		cudaMemcpy(notifyUnit_h, result.notifyUnit, sizeof(int), cudaMemcpyDeviceToHost);
		cudaThreadSynchronize();
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
time += time_tmp;
	}
	free(notifyUnit_h);

	malSize = sizeof(int) * result.numofFormula;
	int *notifyUnknown_h = (int *)malloc(malSize);
	for(int i = 0; i<result.numofFormula; i++)
		notifyUnknown_h[i] = 0;
	cudaMalloc((void**)&result.notifyUnknown, malSize);
	cudaMemcpy(result.notifyUnknown, notifyUnknown_h, malSize, cudaMemcpyHostToDevice);

	int *notifyUnsat_h = (int *)malloc(malSize);
	for(int i = 0; i<result.numofFormula; i++)
		notifyUnsat_h[i] = 0;
	cudaMalloc((void**)&result.notifyUnsat, malSize);
	cudaMemcpy(result.notifyUnsat, notifyUnsat_h, malSize, cudaMemcpyHostToDevice);


cudaEventRecord(start, 0);
	SATChecking_Kernel<<<dimGridBCP, dimBlock>>>(origin_d, result);
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
time += time_tmp;
cudaEventDestroy(start);
cudaEventDestroy(stop);


	cudaMemcpy(formulaState_r, result.formulaState_b, malSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(notifyUnknown_h, result.notifyUnknown, malSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(notifyUnsat_h, result.notifyUnsat, malSize, cudaMemcpyDeviceToHost);

	for(int i = 0; i<result.numofFormula; i++){

		if(notifyUnsat_h[i] == 1)
			formulaState_r[i] = Unsat;		
		else if(notifyUnknown_h[i] == 1)
			formulaState_r[i] = Unknown;
		else
			formulaState_r[i] = Sat;
	}


	free(notifyUnknown_h);
	free(notifyUnsat_h);

	malSize = sizeof(Clause) * result.numofFormula * numofCls;
	cudaMemcpy(cls_r, result.cls_b, malSize, cudaMemcpyDeviceToHost);

	malSize = sizeof(int) * result.numofFormula * numofTotalLit;
	cudaMemcpy(litState_r, result.litState_b, malSize, cudaMemcpyDeviceToHost);

	cudaFree(source.cls_b);
	cudaFree(source.litState_b);
	cudaFree(source.formulaState_b);
	cudaFree(source.expPath_b);
	cudaFree(result.cls_b);
	cudaFree(result.litState_b);
	cudaFree(result.formulaState_b);
	cudaFree(result.expPath_b);
	cudaFree(result.inferPath_b);
	cudaFree(result.notifyUnknown);
	cudaFree(result.notifyUnit);
}

void logicalPro(FormulaBatch result, int * formulaState, Formula *tmp, float GPU_time, clock_t startTime, clock_t endTime, float CPU_time, Formula origin, int *litState, Clause *cls_r, vector<Formula> &vec2, vector<Formula> &vec3, vector<Formula> &vec4, vector<Formula> &vecLarge){

	for(int i=0, j=0; i<result.numofFormula; i++){
		switch(formulaState[i]){
			case Sat:
				printf("Sat");
				tmp[0].formulaState = formulaState[i];
				for(int k=0; k<origin.numofTotalLit; k++)
					tmp[0].litState[k] = litState[i*origin.numofTotalLit + k];
				satPrinter(tmp[0], origin, GPU_time, startTime, endTime, CPU_time);
				exit(0);
			case Unknown:
				int mark1=0;int mark2=0;
				for(int k=0; k<origin.numofCls; k++){
					tmp[j].cls[k] = cls_r[i*origin.numofCls+k];
						if(tmp[j].cls[k].clsState == Unknown && tmp[j].cls[k].countUnassignment>mark1 && tmp[j].cls[k].countUnassignment<=batchSize){
							mark1 = tmp[j].cls[k].countUnassignment;
							mark2 = k;
						}
					}
				for(int k=0; k<origin.numofTotalLit; k++)
					tmp[j].litState[k] = litState[i*origin.numofTotalLit+k];

				for(int k=0,x=0; k<origin.cls[mark2].clsSize; k++)
					if(tmp[j].litState[origin.numofLitBeforeCls[mark2]+ k] == Unassigned){
						tmp[j].expPath[x] = origin.literal[origin.numofLitBeforeCls[mark2]+ k];
						x++;
					}

				for(int k=mark1; k<batchSize; k++)
					tmp[j].expPath[k] = 0;

				switch(mark1){
					case 4: vec4.push_back(tmp[j]);break;
					case 3: vec3.push_back(tmp[j]);break;
					case 2: vec2.push_back(tmp[j]);break;
					default: vecLarge.push_back(tmp[j]);}

				j++;
		}
	}

}


int main(int argc, char** argv){

float GPU_time = 0;
clock_t startTime, endTime;
float CPU_time;

startTime = clock();

//stream of the input file
	FILE * fileIn;


	//temporary objects used while reading the file
	Formula f_origin;
	
	/**********************read file into a formula**********************/
	fileIn = fopen(argv[1], "rb");
	if(fileIn == NULL){
        	fprintf( stderr, "ERROR! Could not open the file: %s \n", argc == 1 ? "<srdin>" : argv[1]), exit(1);}
	else{
        	parse_DIMACS(fileIn, f_origin);
    	}
    	fclose(fileIn);

	//Work on the origin formula and initial origin_d which will be used whole life time of the appliation
	int chosenCls=0;//clause be chosen as the base of expansion
	int expPara=f_origin.cls[0].clsSize;
	int numofCls = f_origin.numofCls;
	int numofTotalLit = f_origin.numofTotalLit;
	size_t malSize;

	for(int i=0; i<numofCls; i++)
		if((f_origin.cls[i].clsSize > f_origin.cls[chosenCls].clsSize) && (f_origin.cls[i].clsSize <= batchSize)){
			chosenCls = i;
			expPara = f_origin.cls[i].clsSize;
		}
	for(int i=0; i<=expPara-1; i++)
		f_origin.expPath[i] = f_origin.literal[f_origin.numofLitBeforeCls[chosenCls]+i];

	Formula origin_cpy = f_origin;
	Formula origin_d;
	origin_d.formulaState = f_origin.formulaState;
	origin_d.numofCls = numofCls;
	origin_d.numofLit = f_origin.numofLit;
	origin_d.numofTotalLit = numofTotalLit;
	origin_d.cls = f_origin.cls;
	origin_d.expPath = f_origin.expPath;
	origin_d.litState = f_origin.litState;

	malSize = sizeof(int) * numofCls;
	cudaMalloc((void **)&origin_d.numofLitBeforeCls, malSize);
	cudaMemcpy(origin_d.numofLitBeforeCls, f_origin.numofLitBeforeCls, malSize, cudaMemcpyHostToDevice);
	malSize = sizeof(int) * numofTotalLit;
	cudaMalloc((void **)&origin_d.literal, malSize);
	cudaMemcpy(origin_d.literal, f_origin.literal, malSize, cudaMemcpyHostToDevice);
	//origin formula done

	vector<Formula> vec2, vec3, vec4, vecLarge;//vectors to store generated formulas with different expPara 2, 3, 4 or larger
	switch(expPara){
		case 2:
			vec2.push_back(origin_cpy);break;
		case 3:
			vec3.push_back(origin_cpy);break;
		case 4:
			vec4.push_back(origin_cpy);break;
		default:
			vecLarge.push_back(origin_cpy);break;
	}

	Formula *tmp = (Formula *)malloc(sizeof(Formula) * batchSize);
	for(int i=0; i<batchSize; i++){
		tmp[i].formulaState = Unknown;
		tmp[i].numofCls = numofCls;
		tmp[i].numofLit = f_origin.numofLit;
		tmp[i].numofTotalLit = numofTotalLit;
		tmp[i].numofLitBeforeCls = f_origin.numofLitBeforeCls;
		tmp[i].cls = (Clause *)malloc(sizeof(Clause) * numofCls);
		tmp[i].expPath = (int *)malloc(sizeof(int) * batchSize);
		tmp[i].litState = (int *)malloc(sizeof(int) * numofTotalLit);
	}

	FormulaBatch source;
	FormulaBatch result;

	while(vec3.size()!=0 || vec4.size()!=0 || vec2.size()!=0 || vecLarge.size()!=0){


		if(vecLarge.size()>0){
			tmp[0] = vecLarge[vecLarge.size()-1];
			for(expPara=0; expPara<batchSize; expPara++)
				if(tmp[0].expPath[expPara] == 0)
					break;
			vecLarge.pop_back();

			source.numofFormula = 1;
			result.numofFormula =  source.numofFormula*expPara;

			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);		

			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
		else if(vec2.size()>5){
			expPara = 2;		
			for(int i=0; i<6; i++){
				tmp[i] = vec2[vec2.size()-1];
				vec2.pop_back();}
			source.numofFormula = 6;
			result.numofFormula =  source.numofFormula*expPara;

			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);

			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
		else if(vec4.size()>2){
			expPara = 4;		
			for(int i=0; i<3; i++){
				tmp[i] = vec4[vec4.size()-1];
				vec4.pop_back();}
			source.numofFormula = 3;
			result.numofFormula =  source.numofFormula*expPara;

			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);

			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
		else if(vec3.size()>3){
			expPara = 3;		
			for(int i=0; i<4; i++){
				tmp[i] = vec3[vec3.size()-1];
				vec3.pop_back();}
			source.numofFormula = 4;
			result.numofFormula =  source.numofFormula*expPara;

			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);
	
			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
		else if(vec4.size()>0){
			expPara = 4;
			source.numofFormula = vec4.size();
			result.numofFormula =  source.numofFormula*expPara;	

			for(int i=vec4.size()-1; i>=0; i--){
				tmp[i] = vec4[vec4.size()-1];

				vec4.pop_back();
			}
	
			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);

			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
		else if(vec3.size()>0){
			expPara = 3;
			source.numofFormula = vec3.size();
			result.numofFormula =  source.numofFormula*expPara;	
			for(int i=vec3.size()-1; i>=0; i--){
				tmp[i] = vec3[vec3.size()-1];
				vec3.pop_back();}
			
			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);

			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
		else{
			expPara = 2;
			source.numofFormula = vec2.size();
			result.numofFormula =  source.numofFormula*expPara;	
			for(int i=vec2.size()-1; i>=0; i--){
				tmp[i] = vec2[vec2.size()-1];
				vec2.pop_back();}
			
			int *litState_r = (int *)malloc(sizeof(int) * result.numofFormula * numofTotalLit);
			int *formulaState_r = (int *)malloc(sizeof(int) * result.numofFormula);
			Clause *cls_r = (Clause *)malloc(sizeof(Clause) * result.numofFormula * numofCls);

			runKernel(origin_d, source, result, expPara, numofCls, numofTotalLit, litState_r, formulaState_r, cls_r, tmp, GPU_time);

			logicalPro(result, formulaState_r, tmp, GPU_time, startTime, endTime, CPU_time, f_origin, litState_r, cls_r, vec2, vec3, vec4, vecLarge);

			free(litState_r);
			free(formulaState_r);
			free(cls_r);
		}
	}
	printf("Unsat");
	printf("     GPU running time: %f", GPU_time);
	endTime = clock();
	CPU_time = (double)( (endTime - startTime) / (double)CLOCKS_PER_SEC);
	printf("     CPU running time: %f", CPU_time);
	exit(0);

};
