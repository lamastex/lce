
sliceBoxesFileName = 'AABCsimdata2DGR_Av_s_10_ns_25000_th_0.02_ni_10000_1.txt'
slDiv1 = regexp(sliceBoxesFileName, '_', 'split')
slDiv1
size(slDiv1,2)
slDiv2 = regexp(slDiv1{size(slDiv1,2)}, '\.txt', 'split')
index = slDiv2{1}
    