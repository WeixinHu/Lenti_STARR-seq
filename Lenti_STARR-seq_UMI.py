import os
import numpy as np
from tqdm import tqdm

def score_martix_fun(martix_type):
	score_martix={}
	lst=['A','G','C','T','N']
	if martix_type=='unitary matrix':
		score_m=[[1,0,0,0,-5],[0,1,0,0,-5],[0,0,1,0,-5],[0,0,0,1,-5],[-5,-5,-5,-5,-5]]
	elif martix_type=='BLAST matrix':
		score_m=[[5,-4,-4,-4,-5],[-4,5,-4,-4,-5],[-4,-4,5,-4,-5],[-4,-4,-4,5,-5],[-5,-5,-5,-5,-5]]
	elif martix_type=='transition-transversion matrix':
		score_m=[[1,-5,-5,-1,-5],[-5,1,-1,-5,-5],[-5,-1,1,-5,-5],[-1,-5,-5,1,-5],[-5,-5,-5,-5,-5]]	
	for i in range(0,len(lst)):
		for j in range(0,len(lst)):
			score_martix[lst[i]+lst[j]]=score_m[i][j]
	return score_martix
 
def Needleman_Wunsch_BLAST(score_martix,gap,long_seq,short_seq):
	score=np.zeros((len(short_seq)+1,len(long_seq)+1))
	for i in range(1,len(short_seq)+1):
		letter_shrot=short_seq[i-1]
		for j in range(1,len(long_seq)+1):
			letter_long=long_seq[j-1]
			diagonal_score_1 = score_martix[letter_shrot+letter_long] + score[i-1][j-1]
			diagonal_score_2 = gap + score[i-1][j]
			diagonal_score_3 = gap + score[i][j-1]
			score[i][j]=max([diagonal_score_1,diagonal_score_2,diagonal_score_3,0])
	current_pos=list(np.unravel_index(np.argmax(score),score.shape))
	alignment_pos=[]
	trace_direct=[]
	while True:
		alignment_pos.append(current_pos)
		trace_line=[score[current_pos[0]-1][current_pos[1]-1],score[current_pos[0]-1][current_pos[1]],score[current_pos[0]][current_pos[1]-1]]
		if max(trace_line)==score[current_pos[0]-1][current_pos[1]-1]:
			current_pos=[current_pos[0]-1,current_pos[1]-1]
			trace_direct.append(1)
		elif max(trace_line)==score[current_pos[0]-1][current_pos[1]]:
			current_pos=[current_pos[0]-1,current_pos[1]]
			trace_direct.append(2)
		elif max(trace_line)==score[current_pos[0]][current_pos[1]-1]:
			current_pos=[current_pos[0],current_pos[1]-1]
			trace_direct.append(3)		
		if score[current_pos[0]][current_pos[1]]==0:
			break		
	alignment_pos.reverse()
	trace_direct.reverse()
	alignment_1=[]
	alignment_2=[]	
	alignment_status=''
	match_num=0
	indel_num=0
	delet_num=0
	misma_num=0
	start_pos=alignment_pos[0][1]-1
	end_pos=alignment_pos[-1][1]
	for i in range(0,len(trace_direct)):
		pos=alignment_pos[i]
		if trace_direct[i]==1:
			alignment_1.append(long_seq[pos[1]-1])
			alignment_2.append(short_seq[pos[0]-1])
			if long_seq[pos[1]-1]==short_seq[pos[0]-1]:
				match_num+=1
			else:
				misma_num+=1
		elif trace_direct[i]==3:
			alignment_1.append(long_seq[pos[1]-1])
			alignment_2.append('-')
			delet_num+=1
		elif trace_direct[i]==2:
			alignment_1.append('-')
			alignment_2.append(short_seq[pos[0]-1])
			indel_num+=1
	alignment_1=''.join(alignment_1)
	alignment_2=''.join(alignment_2)
	return score,alignment_pos,trace_direct,alignment_1,alignment_2,start_pos,end_pos,match_num,delet_num,indel_num,misma_num



score_martix=score_martix_fun('BLAST matrix')
gap=-5

bar_seq_name={}
bar_seq_name['TAGGCTCT']='JK001'
bar_seq_name['GAAGACTG']='JK002'
bar_seq_name['CATTGCAC']='JK003'
bar_seq_name['CGGAAGAA']='JK004'
bar_seq_name['AATCCAGG']='JK005'
bar_seq_name['TGAGGAGA']='JK006'
bar_seq_name['GACTTGGA']='JK007'
bar_seq_name['TCTCACCA']='JK008'
bar_seq_list=['TAGGCTCT','GAAGACTG','CATTGCAC','CGGAAGAA','AATCCAGG','TGAGGAGA','GACTTGGA','TCTCACCA']	

gobal_loc=os.getcwd()
fq_file_list=os.popen('ls %s/*.fq'%gobal_loc).readlines()
for tp in fq_file_list:
	tp=tp.strip().split('\t')[0]
	fq_file=tp.split('/')[-1]
	tp_name=fq_file[:-3]
	os.system('mkdir %s/%s'%(gobal_loc,tp_name))
	f1=open(fq_file,'r')
	m1=f1.readlines()
	f1.close()	
	f0=open('%s/%s/%s_JK000.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f1=open('%s/%s/%s_JK001.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f2=open('%s/%s/%s_JK002.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f3=open('%s/%s/%s_JK003.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f4=open('%s/%s/%s_JK004.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f5=open('%s/%s/%s_JK005.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f6=open('%s/%s/%s_JK006.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f7=open('%s/%s/%s_JK007.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	f8=open('%s/%s/%s_JK008.UMI.fq'%(gobal_loc,tp_name,tp_name),'w')
	for i in range(0,len(m1)):
		if i%4 == 1:
			p1=m1[i].strip().split('\t')
			long_seq=p1[0]
			short_seq='CAGTCTAGTGAATTCG'
			try:
				all_start_pos=long_seq.index(short_seq)
				all_end_pos=all_start_pos+len(short_seq)
				long_seq=long_seq[:all_start_pos]
				error_num=0
				for seq in bar_seq_list:
					short_seq=seq
					score,alignment_pos,trace_direct,alignment_1,alignment_2,start_pos,end_pos,match_num,delet_num,indel_num,misma_num=Needleman_Wunsch_BLAST(score_martix,gap,long_seq,short_seq)
					mis_type=str(match_num)+'M'+str(delet_num)+'D'+str(indel_num)+'I'+str(misma_num)+'S'
					if mis_type[:2]=='8M':
						if seq==bar_seq_list[0]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f1.write(m1[i-1])
							f1.write(tp_seq[end_pos:all_start_pos-4])
							f1.write('\n')
							f1.write(m1[i+1])
							f1.write(qc[end_pos:all_start_pos-4])				
							f1.write('\n')
							break

						if seq==bar_seq_list[1]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f2.write(m1[i-1])
							f2.write(tp_seq[end_pos:all_start_pos-4])
							f2.write('\n')
							f2.write(m1[i+1])
							f2.write(qc[end_pos:all_start_pos-4])				
							f2.write('\n')
							break
							
						if seq==bar_seq_list[2]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f3.write(m1[i-1])
							f3.write(tp_seq[end_pos:all_start_pos-4])
							f3.write('\n')
							f3.write(m1[i+1])
							f3.write(qc[end_pos:all_start_pos-4])				
							f3.write('\n')													
							break

						if seq==bar_seq_list[3]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f4.write(m1[i-1])
							f4.write(tp_seq[end_pos:all_start_pos-4])
							f4.write('\n')
							f4.write(m1[i+1])
							f4.write(qc[end_pos:all_start_pos-4])				
							f4.write('\n')			
							break

						if seq==bar_seq_list[4]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f5.write(m1[i-1])
							f5.write(tp_seq[end_pos:all_start_pos-4])
							f5.write('\n')
							f5.write(m1[i+1])
							f5.write(qc[end_pos:all_start_pos-4])				
							f5.write('\n')
							break

						if seq==bar_seq_list[5]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f6.write(m1[i-1])
							f6.write(tp_seq[end_pos:all_start_pos-4])
							f6.write('\n')
							f6.write(m1[i+1])
							f6.write(qc[end_pos:all_start_pos-4])				
							f6.write('\n')
							break
							
						if seq==bar_seq_list[6]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f7.write(m1[i-1])
							f7.write(tp_seq[end_pos:all_start_pos-4])
							f7.write('\n')
							f7.write(m1[i+1])
							f7.write(qc[end_pos:all_start_pos-4])				
							f7.write('\n')													
							break
							
						if seq==bar_seq_list[7]:
							tp_seq=p1[0]
							qc=m1[i+2].strip().split('\t')[0]
							f8.write(m1[i-1])
							f8.write(tp_seq[end_pos:all_start_pos-4])
							f8.write('\n')
							f8.write(m1[i+1])
							f8.write(qc[end_pos:all_start_pos-4])				
							f8.write('\n')		
							break							
					else:
						error_num+=1
				if error_num==8:
					f0.write(m1[i-1])
					f0.write(m1[i])
					f0.write(m1[i+1])
					f0.write(m1[i+2])			
			except:
				no=1
	f1.close()
	f2.close()
	f3.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()
	f7.close()
	f8.close()
	f0.close()

