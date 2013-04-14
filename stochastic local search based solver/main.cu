#include "struct_def.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>

using namespace std;

__global__ void LocalSearch_K(int *, int *, int *, int, int, int);
__global__ void LocalSearchSum_K(int *, int *, int, int);

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


static void readClause(char **in, Formula &f, Clause *cls, vector<int> &parsedLit, int currentCls, int &maxSize){
    	
	int tmp;
	int counter=0;

    	for(;;){
        	tmp = parseInt(in);
        	if (tmp == 0) {
			cls[currentCls].clsSize = counter;
			maxSize = counter > maxSize ? counter : maxSize;
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
	
    	for(;;){
        	skipWhitespace(&text);
        	if (*text == 0)
            		break;
		else if (*text == 'p'){
			f.numofVar = getPara(&text);
			f.numofCls = getPara(&text);
			f.cls = (Clause *)malloc(sizeof(Clause)*f.numofCls);			
			f.numofLitBeforeCls = (int *)malloc(sizeof(int) * f.numofCls);
			f.numofLitBeforeCls[0] = 0;
			f.MAX_CLS_SIZE = 0;
			skipLine(&text);
		}
        	else if (*text == 'c')
            		skipLine(&text);
        	else{
			if(clsCounter<f.numofCls){
            			readClause(&text, f, f.cls, parsedLit, clsCounter, f.MAX_CLS_SIZE);
				clsCounter++;
			}
		}
    	}
	f.numofTotalLit = f.numofLitBeforeCls[f.numofCls-1] + f.cls[f.numofCls-1].clsSize;
	int allLit = f.numofTotalLit;
	f.literal = (int *)malloc(sizeof(int) * allLit);
	for(int i=allLit-1; i>=0; i--){
		f.literal[i] = parsedLit[i];
		parsedLit.pop_back();
	}	
	//free(text);	
}


//Following two methods are for additional guiding search, use them by uncommment them and recompile
/**
void prepareDPLL(int *UCls_record, int numofVar, int numofCls, int *ref_formula, int max_cls_size, int *ass_record, Formula &f, int * weight, GuideVariable * guide, int &chosenCls){

	int minNumofLocalMinCls = UCls_record[0];
	vector<int> LocalMinAss;
	LocalMinAss.push_back(0);
	for(int i = 1; i < RestartTimes; i++)
		if(minNumofLocalMinCls > UCls_record[i]){
			minNumofLocalMinCls = UCls_record[i];
			LocalMinAss.clear();
			LocalMinAss.push_back(i);
		}else if(minNumofLocalMinCls == UCls_record[i])
			LocalMinAss.push_back(i);
	int numofLocalMin = LocalMinAss.size();

	int size = numofLocalMin*minNumofLocalMinCls;
	int * UClsSet = (int *)malloc(sizeof(int) * size);
	for(int assLoop=0; assLoop<numofLocalMin; assLoop++){
		int base = LocalMinAss[assLoop] * numofVar;
		int UCls_counter = 0;
		for(int clsLoop=0; clsLoop<numofCls; clsLoop++){
			int clsState = 0;
			int base = clsLoop * max_cls_size;
			for(int litLoop=0; litLoop<max_cls_size; litLoop++){
				int tmp = ref_formula[base + litLoop];
				
				if(tmp<0){
					tmp = 0-tmp;
					clsState = clsState || (1-ass_record[base + tmp -1]);
				}else if(tmp>0)
					clsState = clsState || ass_record[base + tmp -1];
				else{
					if(clsState == UNSAT){
						UClsSet[assLoop*minNumofLocalMinCls + UCls_counter] = (clsLoop+1);
						UCls_counter++;
					}
					break;
				}
				if(litLoop == max_cls_size-1){
					if(clsState == UNSAT){
						UClsSet[assLoop*minNumofLocalMinCls + UCls_counter] = (clsLoop+1);
						UCls_counter++;
					}
					break;
				}
					
			}
			if(UCls_counter == minNumofLocalMinCls)
				break;
		}
	}

	int * hardClsCounter = (int *)malloc(sizeof(int)*size);
	for(int outer_loop = 0; outer_loop < size; outer_loop++){
		hardClsCounter[outer_loop] = 1;
		for(int inner_loop = outer_loop+1; inner_loop < size; inner_loop++)
			if(UClsSet[outer_loop] == UClsSet[inner_loop])
				hardClsCounter[outer_loop]++;
	}
	chosenCls = UClsSet[0];
	int chosenCls_counter = hardClsCounter[0];
	int chosenAss = 0;
	for(int idx=1; idx<size; idx++)
		if(chosenCls_counter < hardClsCounter[idx]){
			chosenCls = UClsSet[idx];
			chosenCls_counter = hardClsCounter[idx];
			chosenAss = LocalMinAss[idx/minNumofLocalMinCls];
		}
		else if(chosenCls_counter == hardClsCounter[idx] && (f.cls[chosenCls].clsSize < f.cls[UClsSet[idx]].clsSize) ){
			chosenCls = UClsSet[idx];
			chosenCls_counter = hardClsCounter[idx];
			chosenAss = LocalMinAss[idx/minNumofLocalMinCls];
		}

	//initial variable guide
	int base = chosenAss * numofVar;
	GuideVariable tmpVar;
	for(int i=0; i<numofVar; i++){
		tmpVar.var = i+1;
		tmpVar.val = ass_record[base + i];
		tmpVar.weight = weight[i];
		guide[i] = tmpVar;
	}
	//sorting by weight
	for(int outer_loop = numofVar-1; outer_loop > 0; outer_loop--)
		for(int inner_loop = 0; inner_loop < outer_loop; inner_loop++)
			if(guide[inner_loop].weight > guide[inner_loop+1].weight){
				GuideVariable tmp = guide[inner_loop+1];
				guide[inner_loop+1] = guide[inner_loop];
				guide[inner_loop] = tmp;
			}
}

int DPLLwithGuide(int numofVar, int numofCls, GuideVariable *guideInput, Formula &f, int initVar, int max_cls_size, int *ref_formula){
	
	int *varSet = (int *)malloc(sizeof(int) * numofVar);
	GuideVariable *guide = (GuideVariable *)malloc(sizeof(int) * numofVar * 30);
	int *clsState = (int *)malloc(sizeof(int) * numofCls);
	int *clsSize = (int *)malloc(sizeof(int) * numofCls);

	for(int i=0; i<numofVar; i++){
		varSet[i] = UNSET;
		guide[i] = guideInput[i];
	}
	for(int j=0; j<numofCls; j++){
		clsState[j] = UNKNOWN;
		clsSize[j] = f.cls[j].clsSize;
	}
	int guide_pointer = numofVar;

	if(initVar>0){
		guide[guide_pointer].var = initVar;
		guide[guide_pointer].val = 1;
	}else if(initVar<0){
		guide[guide_pointer].var = 0-initVar;
		guide[guide_pointer].val = 0;
	}
	
	int fStateSAT = 0;

	while(guide_pointer>=0){

		int variable;
		int value;

		while( varSet[(guide[guide_pointer].var-1)] == SET  && guide_pointer>=0)
			guide_pointer--;
		if(guide_pointer<0){
			return UNSAT;
		}

		variable = guide[guide_pointer].var;
		value = guide[guide_pointer].val;
		guide_pointer--;
		int offset = 0;
		for(int clsLoop=0; clsLoop<numofCls; clsLoop++){
			if( clsState[clsLoop] != SAT ){
				int base = clsLoop * max_cls_size;
				for(int litLoop=0; litLoop<max_cls_size; litLoop++){
					
					int tmp = ref_formula[base + litLoop];
				
					if(tmp<0){
						tmp = 0-tmp;
						if(varSet[tmp-1] == UNSET)
							if( tmp == variable && value == 0) {
								clsSize[clsLoop]--;
								clsState[clsLoop] = SAT;
								fStateSAT++;
								break;
							}else if( tmp == variable && value == 1) {
								clsSize[clsLoop]--;
								if(clsSize[clsLoop] == 1){
									offset++;
									for(int unit = 0; unit < max_cls_size; unit++)
										if(varSet[ref_formula[base+unit]] == UNSET)
											if(ref_formula[base+unit]<0){
												guide[guide_pointer+offset].var = 0-ref_formula[base+unit];
												guide[guide_pointer+offset].val = 0;
												break;
											}else{
												guide[guide_pointer+offset].var = ref_formula[base+unit];
												guide[guide_pointer+offset].val = 1;
												break;
											}
								}
							}
					}else if(tmp>0)
						if(varSet[tmp-1] == UNSET)
							if( tmp == variable && value == 1) {
								clsSize[clsLoop]--;
								clsState[clsLoop] = SAT;
								fStateSAT++;
								break;
							}else if( tmp == variable && value == 0) {
								clsSize[clsLoop]--;
								if(clsSize[clsLoop] == 1){
									offset++;
									for(int unit = 0; unit < max_cls_size; unit++)
										if(varSet[ref_formula[base+unit]] == UNSET)
											if(ref_formula[base+unit]<0){
												guide[guide_pointer+offset].var = 0-ref_formula[base+unit];
												guide[guide_pointer+offset].val = 0;
												break;
											}else{
												guide[guide_pointer+offset].var = ref_formula[base+unit];
												guide[guide_pointer+offset].val = 1;
												break;
											}
								}
							}
					else{
						if(clsSize[clsLoop] == 0){
							clsState[clsLoop] = UNSAT;
							return UNSAT;
						}
						break;
					}
					if( litLoop == max_cls_size-1 ){
						if(clsSize[clsLoop] == 0){
							clsState[clsLoop] = UNSAT;
							clsState[clsLoop] = UNSAT;
							return UNSAT;
						}
						break;
					}
				}
			}
			if(fStateSAT == numofCls){
				clsState[clsLoop] = UNSAT;
				return SAT;
			}
		}
		guide_pointer += offset;
		varSet[variable-1] = SET;	
	}
	return SAT;
}**/

void reportSAT(int *init_ass, int *ass_record, int *weight, int *d_ref_formula, int *d_ass, int *d_numofUNSATCls, int *d_sum_numofUNSATCls, int *sum_numofUNSATCls, int *UCls_record, float &GPU_time, clock_t &startTime, clock_t &endTime, float &CPU_time, cudaEvent_t &start, cudaEvent_t &stop){
	printf("SAT");
	endTime = clock();
	CPU_time = (double)( (endTime - startTime) / (double)CLOCKS_PER_SEC);
	printf("     CPU running time: %f", CPU_time);
printf("     GPU running time(ms): %f", GPU_time);
	printf("\n");
	cudaFree(d_ref_formula);
	cudaFree(d_ass);
	cudaFree(d_numofUNSATCls);
	cudaFree(d_sum_numofUNSATCls);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	exit(0);
}

void reportUNSAT(int *init_ass, int *ass_record, int *weight, int *d_ref_formula, int *d_ass, int *d_numofUNSATCls, int *d_sum_numofUNSATCls, int *sum_numofUNSATCls, int numofVar, int *UCls_record, float &GPU_time, clock_t &startTime, clock_t &endTime, float &CPU_time, cudaEvent_t &start, cudaEvent_t &stop){
	printf("UNSAT");
endTime = clock();
	CPU_time = (double)( (endTime - startTime) / (double)CLOCKS_PER_SEC);
	printf("     CPU running time: %f", CPU_time);
printf("     GPU running time(ms): %f", GPU_time);
	printf("\n");

	cudaFree(d_ref_formula);
	cudaFree(d_ass);
	cudaFree(d_numofUNSATCls);
	cudaFree(d_sum_numofUNSATCls);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	exit(0);
}

void randomWalk(int * init_ass, int * ass_record, int * UCls_record, int numofVar, int &var, int walk_counter, int &track_pointer, vector<int> UNSATCls_track, int *weight, int &trap_counter){
	
	trap_counter = 0;	

	int offset = walk_counter * numofVar;
	for(int var_counter = 0; var_counter < numofVar; var_counter++){
		ass_record[offset+var_counter] = init_ass[var_counter];
		if(weight[var_counter]<1)
			init_ass[var_counter]=1-init_ass[var_counter];
		weight[var_counter] = 0;
	}

	UCls_record[walk_counter] = UNSATCls_track.back();
	
	var = rand()%numofVar + 1;

	UNSATCls_track.clear();
	track_pointer = -1;
	weight[var-1]++;
}


int main(int argc, char** argv){

float GPU_time = 0;
clock_t startTime, endTime;
float CPU_time;

startTime = clock();

//Read in DIMAS file
	//stream of the input file
	FILE * fileIn;

	//temporary objects used while reading the file
	Formula f;
	
	/**********************read file into a formula**********************/
	fileIn = fopen(argv[1], "rb");
	if(fileIn == NULL){
        	fprintf( stderr, "ERROR! Could not open the file: %s \n", argc == 1 ? "<srdin>" : argv[1]), exit(1);}
	else{
        	parse_DIMACS(fileIn, f);
    	}
    	fclose(fileIn);

	//Prepare for kernels
	int numofVar = f.numofVar;
	int numofCls = f.numofCls;
	int max_cls_size = f.MAX_CLS_SIZE;
	int gridX = (int)((numofCls+thread_per_block-1)/thread_per_block);
	int gridY = (int)((numofVar+thread_per_block-1)/thread_per_block);	

	vector<int> UNSATCls_track;
	int track_pointer = -1;
	int inverseVar = 1;
	int inverse_tmp;

	int * init_ass;
	int * ass_record;
	int * UCls_record;
	int * weight;
	GuideVariable * guide;
	guide = (GuideVariable *)malloc(sizeof(GuideVariable)*numofVar);

	//Create the formula stored in Global memory
	int * ref_formula = (int *)malloc(sizeof(int) * gridX * thread_per_block * max_cls_size);
	for(int cls_index=0; cls_index<numofCls; cls_index++)
		for(int lit_index=0; lit_index<max_cls_size; lit_index++){
			int c_i = cls_index*max_cls_size;
			if(lit_index < f.cls[cls_index].clsSize){
				int l_i = f.numofLitBeforeCls[cls_index];
				ref_formula[c_i + lit_index] = f.literal[ l_i + lit_index];
			}else
				ref_formula[c_i + lit_index] = 0;
		}
	for(int cls_index = numofCls; cls_index < gridX * thread_per_block; cls_index++){
		int c_i = cls_index*max_cls_size;
		for(int lit_index = 0; lit_index<max_cls_size; lit_index++)
			ref_formula[c_i+lit_index] = 0;
	}

	//Create random initial assignment
	init_ass = (int *)malloc(sizeof(int) * numofVar);
	ass_record = (int *)malloc(sizeof(int) * numofVar * RestartTimes);
	UCls_record = (int *)malloc(sizeof(int) * RestartTimes);
	

	srand((unsigned)(time(NULL)));
	for(int var_index=0; var_index<numofVar; var_index++)
		init_ass[var_index] = rand()%2;

	//Create weight for variables
	weight = (int *)malloc(sizeof(int) * numofVar);
	for(int i=0; i<numofVar; i++)
		weight[i] = 0;

	//copy to the GPU side
	size_t malSize;
	int *d_ref_formula;
	int *d_ass;
	int *d_numofUNSATCls;
	int *d_sum_numofUNSATCls;
	int * sum_numofUNSATCls;
	int sharedMemSize = sizeof(int) * (thread_per_block * max_cls_size + numofVar);

	malSize = sizeof(int) * numofVar;
	sum_numofUNSATCls = (int *)malloc(malSize);
	cudaMalloc((void **)&d_sum_numofUNSATCls, malSize);

	malSize = sizeof(int) * numofCls * max_cls_size;
	cudaMalloc((void **)&d_ref_formula, malSize);
	cudaMemcpy(d_ref_formula, ref_formula, malSize, cudaMemcpyHostToDevice);

	malSize = sizeof(int) * numofVar;
	cudaMalloc((void **)&d_ass, malSize);
	cudaMemcpy(d_ass, init_ass, malSize, cudaMemcpyHostToDevice);
	
	dim3 grid(gridX, gridY);
	dim3 block(thread_per_block, 1);
	dim3 gridSum(gridY,1);

	malSize = sizeof(int) * gridY *gridX * thread_per_block;
	cudaMalloc((void **)&d_numofUNSATCls, malSize);

	int trap_counter = 0;

float time_tmp; 
cudaEvent_t start;
cudaEvent_t stop;

cudaEventCreate(&start);
cudaEventCreate(&stop);


	for(int walk_counter = 0; walk_counter < RestartTimes; ){
		for(; trap_counter < Trapped; trap_counter++){



cudaEventRecord(start, 0);
			LocalSearch_K<<<grid, block, sharedMemSize>>>(d_ref_formula, d_ass, d_numofUNSATCls, numofVar, numofCls, max_cls_size);
			cudaStreamSynchronize(0);
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
GPU_time += time_tmp;

			malSize = sizeof(int) * numofVar;

cudaEventRecord(start, 0);
			LocalSearchSum_K<<<gridSum, block>>>(d_numofUNSATCls, d_sum_numofUNSATCls, numofVar, gridX);
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_tmp, start, stop);
GPU_time += time_tmp;


			cudaMemcpy(sum_numofUNSATCls, d_sum_numofUNSATCls, malSize, cudaMemcpyDeviceToHost);

			//flip a variable to approach the local best
			int UNSAT_min = sum_numofUNSATCls[0];
			int differ_weight=0;
			inverse_tmp = 1;
			
			for(int var_decision = 0; var_decision<numofVar; var_decision++)
				if(UNSAT_min > sum_numofUNSATCls[var_decision]){
							UNSAT_min = sum_numofUNSATCls[var_decision];
							inverse_tmp = var_decision+1;
							differ_weight = 0;
				}else if (UNSAT_min == sum_numofUNSATCls[var_decision]){
							int flag_weight = weight[inverse_tmp-1] - weight[var_decision];
							if(flag_weight > 0){
								inverse_tmp = var_decision+1;
								differ_weight = flag_weight;
							}
							else if(flag_weight < 0)
								differ_weight = ((differ_weight+flag_weight)>0) ? (-flag_weight) : differ_weight;
							else
								inverse_tmp = (rand()%2==0) ? inverse_tmp : (var_decision+1);
				}		
			
			inverseVar = inverse_tmp;
			weight[inverseVar-1] += differ_weight;
			weight[inverseVar-1] ++;
			UNSATCls_track.push_back(UNSAT_min);
			track_pointer++;

			if(UNSAT_min == 0)
				reportSAT(init_ass, ass_record, weight, d_ref_formula, d_ass, d_numofUNSATCls, d_sum_numofUNSATCls, sum_numofUNSATCls, UCls_record, GPU_time, startTime, endTime, CPU_time, start, stop);
			else if( (track_pointer > 0) && (UNSAT_min > UNSATCls_track[track_pointer-1]) ){
				randomWalk(init_ass, ass_record, UCls_record, numofVar, inverseVar, walk_counter, track_pointer, UNSATCls_track, weight, trap_counter);
				
				walk_counter++;
			}
			else if( (track_pointer > 0) && (UNSAT_min < UNSATCls_track[track_pointer-1]) )
				
					trap_counter = 0;

			init_ass[inverseVar-1] = 1 - init_ass[inverseVar-1];
			cudaMemcpy(d_ass, init_ass, sizeof(int) * numofVar, cudaMemcpyHostToDevice);
		}
		//do random walk with a probability
		int rwPro = rand()%10;
		if(rwPro < P){
			randomWalk(init_ass, ass_record, UCls_record, numofVar, inverseVar, walk_counter, track_pointer, UNSATCls_track, weight, trap_counter);
			walk_counter++;
			init_ass[inverseVar-1] = 1 - init_ass[inverseVar-1];
			cudaMemcpy(d_ass, init_ass, sizeof(int) * numofVar, cudaMemcpyHostToDevice);
		}
	}
//for additional guiding search
/**
	int chosenCls;
	prepareDPLL(UCls_record, numofVar, numofCls, ref_formula, max_cls_size, ass_record, f, weight, guide, chosenCls);

	for(int DPLLLoop = 0; DPLLLoop < f.cls[chosenCls].clsSize; DPLLLoop++){
		int result = DPLLwithGuide(numofVar, numofCls, guide, f, ref_formula[(chosenCls-1)*max_cls_size+DPLLLoop], max_cls_size, ref_formula);
		if(result == SAT)
			reportSAT(init_ass, ass_record, weight, d_ref_formula, d_ass, d_numofUNSATCls, d_sum_numofUNSATCls, sum_numofUNSATCls,UCls_record, GPU_time, startTime, endTime, CPU_time, start, stop);
	}
**/

	reportUNSAT(init_ass, ass_record, weight, d_ref_formula, d_ass, d_numofUNSATCls, d_sum_numofUNSATCls, sum_numofUNSATCls, numofVar, UCls_record, GPU_time, startTime, endTime, CPU_time, start, stop);
}
	
