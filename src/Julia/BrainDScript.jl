#Pateient 4 Visit 1 default values
#default values
 #const ρ = [0.99; 1; 1]
 #const t= [60; 54; 100]
 #const tmax = [62; 54; 100]

 DataDir = "./Data"
 Patient = "Patient4"
 Visit = "Visit1"
 extra = "refpointpercent50_value26"
 gamma_coef_file_name=@sprintf("./%s/%s_%s_%s_Gamma_dist_func.jld2",Patient,Patient, Visit, extra)
 Basic_file=@sprintf("%s/%s/%s_%s_16beams_%s_notincludingdeadvoxels_20230719.mat",DataDir,Patient,Patient,Visit,extra)
 RegFileName=@sprintf("%s/%s/%s_%s_%s_Gamma_dist_func.jld2",DataDir,Patient,Patient,Visit,extra)
 D_file=@sprintf("%s/%s/%s_%s_D_formatted.jld2",DataDir,Patient,Patient,Visit)

 file_basic = matopen(Basic_file)
 γ = read(file_basic, "neighbors_Mat")
 ϕ = read(file_basic, "omf_Vec")
 close(file_basic)

D = []
firstIndices = []
dvrhs = []
OARSize = []
 if bLOAD_FROM_FILE
     D = FileIO.load(D_file, "D")
     firstIndices, dvrhs, OARSize = FileIO.load(D_file, "firstIndices", "dvrhs", "OARSize")
     if length(firstIndices)<length(t)
         global firstIndices=[firstIndices; size(D,1)]
    else
        global firstIndices= firstIndices
    end
 else
     file = matopen(Basic_file) #changed from 13 since it did not cover the PTV
     inD = read(file, "Dij")
     V = read(file, "V")
     close(file)
     println("initial D size: ",size(inD),"  OAR num (length(V)-1): ",length(V) - 1," D nnz: ",nnz(inD))
     #println(inD.nzval)
     #sumVec = sum(inD,2)
     nonzeroIdxs = unique(inD.rowval) #finall(sumVec->sumVec>0,sumVec)
	 #V might include additional organs such as dead cells
	 @assert(length(V)>=length(t)+1)
	 println("nonzero rows: ",length(nonzeroIdxs)," min: ",minimum(nonzeroIdxs)," max: ",maximum(nonzeroIdxs))
     nb = size(inD, 2)
     D = spzeros(0, nb)
     firstIndices = [] # vector of indices of first voxels for each of the stuctures
     dvrhs = zeros(length(t))
     OARSize = zeros(length(t))
     # skipping the last strucure that includes dead voxels?
     for k = 1:length(t)+1 #length(V)-1
         oarIndices = Array{Int}(vec(V[k]))
         idxs = intersect(nonzeroIdxs, oarIndices) #findall(in(nonzeroIdxs),oarIndices)
         println("struct: ",k," before removing rows: ", length(V[k]), " after removing rows: ", length(idxs), " ", typeof(nonzeroIdxs), " ", typeof(oarIndices))
         #println(minimum(oarIndices), " ", maximum(oarIndices), " ")
         appendD = inD[idxs, :]
         rowN, colN = size(appendD)
		 if k > 1 #size(D,1)>0
	         println(size(D))
	         global firstIndices = [firstIndices; size(D, 1) + 1]
	     end
         println("For struct: ", k, " appending submat of D of size: ", size(appendD))
         if rowN > 0
             global D = [D; appendD]
             if k > 1
	         OARSize[k-1] = length(oarIndices)
                 dvrhs[k-1] = floor((1 - ρ[k-1]) * length(V[k]))
                 println("DVRHS: ", dvrhs[k-1] , " for organ: ", k-1)
             end
         end
     end
     if bSAVE_FILES
         FileIO.save(D_file,"D",D,"firstIndices",firstIndices,"dvrhs",dvrhs,"OARSize",OARSize)
     end
     close(file)
#	 println(firstIndices)
 end
 println(OARSize)
 dvrhs = (-ρ.+1.0).*OARSize

 inD = []
 println(firstIndices)
 @assert(size(γ, 1) == firstIndices[1] - 1)#need to make sure these are the same

 Din=D
