#ifndef _BCP_KERNEL_H_
#define _BCP_KERNEL_H_

#include "struct_def.h"

__global__ void BCP_Kernel(Formula origin, FormulaBatch s){

	int formulaIdx = blockIdx.y;
	int blockCol = blockIdx.x;
	int clauseCol = threadIdx.x;
	__shared__ int d_inferPath[batchSize];
	if(clauseCol<batchSize)
		d_inferPath[clauseCol] = s.inferPath_b[formulaIdx*batchSize + clauseCol];
	__syncthreads();

	int clsIdx_b = formulaIdx*origin.numofCls + blockCol*Thread_per_Block + clauseCol;

	if((blockCol*Thread_per_Block + clauseCol < origin.numofCls) && (s.cls_b[clsIdx_b].clsState == Unknown) ){
		
		int litIdx_b = formulaIdx*origin.numofTotalLit + origin.numofLitBeforeCls[blockCol*Thread_per_Block + clauseCol];
		int litIdx = origin.numofLitBeforeCls[blockCol*Thread_per_Block + clauseCol];

		for(int inferCounter=0; inferCounter<batchSize; inferCounter++){
			
			if(s.cls_b[clsIdx_b].clsState == Unknown){

				for(int counterinCls=s.cls_b[clsIdx_b].clsSize-1; counterinCls>=0; counterinCls--)

					if( (s.litState_b[litIdx_b+counterinCls]==Unassigned) && (origin.literal[litIdx+counterinCls]-d_inferPath[inferCounter]==0) ){
						s.litState_b[litIdx_b+counterinCls] = LitTrue;
						s.cls_b[clsIdx_b].countUnassignment--;
						s.cls_b[clsIdx_b].clsState = Sat;
						break;
						
					}else if( (s.litState_b[litIdx_b+counterinCls]==Unassigned) && (origin.literal[litIdx+counterinCls]+d_inferPath[inferCounter]==0)){
						s.litState_b[litIdx_b+counterinCls] = LitFalse;
						s.cls_b[clsIdx_b].countUnassignment--;
						if((s.cls_b[clsIdx_b].countUnassignment == 0) && (s.cls_b[clsIdx_b].clsState == Unknown) ){
							s.cls_b[clsIdx_b].clsState = Unsat;
							break;
						}
					}
			}else
				break;
		}
	}
}
#endif
