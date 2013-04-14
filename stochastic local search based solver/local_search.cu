#ifndef _LOCAL_SEARCH_H_
#define _LOCAL_SEARCH_H_

#include "struct_def.h"

__global__ void LocalSearch_K(int *d_ref_formula, int *d_ass, int *d_numofUNSATCls, int numofVar, int numofCls, int max_cls_size){

	int clsDiv = blockIdx.x;
	int varDiv = blockIdx.y;
	int tid = threadIdx.x;
	const int stepChangingIdx = varDiv * thread_per_block + tid + 1;
	const int rangementY = varDiv * thread_per_block  + tid;

	extern __shared__ int array[];

	int * cls_partial = array;
	int * s_ass = (int *)&cls_partial[max_cls_size * thread_per_block];
	
	//load from global memory to shared memory
	for(int lit_counter = 0; lit_counter < max_cls_size; lit_counter++){
		int idx_lit = tid + lit_counter*thread_per_block;
		cls_partial[idx_lit] = d_ref_formula[thread_per_block * clsDiv * max_cls_size + idx_lit];
	}

	for(int var_counter = 0; var_counter * thread_per_block < numofVar; var_counter++){
			int idx_var = tid + var_counter * thread_per_block;
			if(idx_var < numofVar)
				s_ass[idx_var] = d_ass[idx_var];
	}
	__syncthreads();

	if(rangementY < numofVar){
		
		int usefulPartial = thread_per_block;
		if(clsDiv == gridDim.x-1)
			usefulPartial = numofCls - clsDiv*thread_per_block;
		int numofUNSATCls = 0;

		//calculate number of UNSAT clauses for each 1-bitwise assignment
		for(int cls_loop = 0; cls_loop < usefulPartial; cls_loop++){
			int clsSAT = UNSAT;
			for(int lit_loop = 0; lit_loop < max_cls_size; lit_loop++){
				int Idx = cls_partial[lit_loop + cls_loop * max_cls_size];
				if(Idx == 0){}				
				else if(Idx<0){
					Idx = 0-Idx;
					if( stepChangingIdx != Idx)
						clsSAT = clsSAT || (1-s_ass[Idx-1]);
					else
						clsSAT = clsSAT || s_ass[Idx-1];
				}
				else{
					if( stepChangingIdx != Idx)
						clsSAT = clsSAT || s_ass[Idx-1];
					else
						clsSAT = clsSAT || (1-s_ass[Idx-1]);
				}
			}
			numofUNSATCls = numofUNSATCls + 1 - clsSAT;
		}

		d_numofUNSATCls[(gridDim.y * clsDiv + varDiv) * thread_per_block + tid] = numofUNSATCls;
	}
}

//sum up number of UNSAT cls, try to substitute with CPU realization next time if performance is bad
__global__ void LocalSearchSum_K(int *d_numofUNSATCls, int *d_sum_numofUNSATCls, int numofVar, int offset){
	
	int varDiv = blockIdx.x;
	int tid = threadIdx.x;

	const int Idx = tid + varDiv * thread_per_block;

	if(Idx < numofVar){
		int sum = 0;
		for(int offset_tmp = 0; offset_tmp < offset; offset_tmp++)
			sum = sum + d_numofUNSATCls[Idx + gridDim.x * thread_per_block * offset_tmp];
		d_sum_numofUNSATCls[Idx] = sum;
	}
}

#endif
