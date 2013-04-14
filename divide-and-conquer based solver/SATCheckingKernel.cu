#ifndef _SATChecking_KERNEL_H_
#define _SATChecking_KERNEL_H_

#include "struct_def.h"

__global__ void SATChecking_Kernel(Formula origin, FormulaBatch s){

	int formulaIdx = blockIdx.y;
	int blockCol = blockIdx.x;
	int clauseCol = threadIdx.x;

	
	int clsIdx_b = formulaIdx*origin.numofCls + blockCol*Thread_per_Block + clauseCol;

	if(blockCol*Thread_per_Block + clauseCol < origin.numofCls){

		if(s.cls_b[clsIdx_b].clsState == Unsat)
			s.notifyUnsat[formulaIdx] = 1;
		if(s.cls_b[clsIdx_b].clsState == Unknown)
			s.notifyUnknown[formulaIdx] = 1;
	}
}
#endif
