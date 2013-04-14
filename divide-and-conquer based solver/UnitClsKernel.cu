#ifndef _UnitCls_KERNEL_H_
#define _UnitCls_H_

#include "struct_def.h"

__global__ void UnitCls_Kernel(Formula origin, FormulaBatch s){

	int formulaIdx = blockIdx.y;
	int blockCol = blockIdx.x;
	int clauseCol = threadIdx.x;


	int clsIdx_b = formulaIdx*origin.numofCls + blockCol*Thread_per_Block + clauseCol;

	if((blockCol*Thread_per_Block + clauseCol < origin.numofCls) && (s.cls_b[clsIdx_b].clsState == Unknown) ){
		
		
		int litIdx_b = formulaIdx*origin.numofTotalLit + origin.numofLitBeforeCls[blockCol*Thread_per_Block + clauseCol];
		int litIdx = origin.numofLitBeforeCls[blockCol*Thread_per_Block + clauseCol];

		int remainLit = 0;

	
		if( (s.cls_b[clsIdx_b].countUnassignment == 1)&& (s.cls_b[clsIdx_b].clsState == Unknown)){
			
			for(int counterinCls=s.cls_b[clsIdx_b].clsSize-1; counterinCls>=0; counterinCls--)
				if(s.litState_b[litIdx_b+counterinCls] == Unassigned){
					remainLit = origin.literal[litIdx+counterinCls];
					break;
				}

			s.inferPath_b[formulaIdx*batchSize + clsIdx_b%batchSize] = remainLit;
			s.notifyUnit[0] = 1; 
		}

	}
}
#endif
