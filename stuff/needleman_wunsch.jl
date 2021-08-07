#match = 1, mismatch = -1, gap = -1
function nw(seq1, seq2)
	matrix_ = fill(0,length(seq2)+1,length(seq1)+1)
	trace_ = Matrix{Union{Char, Nothing}}(nothing,length(seq2)+1,length(seq1)+1)

	#first row and column
	for i = 0:length(seq1)
		matrix_[1, i+1] = -i
		trace_[1, i+1] = 'L'
	end
	for i = 0:length(seq2)
		matrix_[i+1, 1] = -i
		trace_[i+1, 1] = 'U'
	end
	trace_[1,1] = '-' 
	
	#fill
	for i = 2:length(seq2)+1
		for j = 2:length(seq1)+1
			#x stores whether it was a match or mismatch/indel
			if seq1[j-1]==seq2[i-1] x = 1 else x = -1 end
			
			top = matrix_[i-1,j] - 1
			left = matrix_[i, j-1] - 1 
			di = matrix_[i-1,j-1] + x
			
			if di > top && di > left
				matrix_[i,j] = di
				trace_[i,j] = 'D'
			elseif top > left
				matrix_[i,j] = top
				trace_[i,j] = 'U'
			else 
				matrix_[i,j] = left
				trace_[i,j] = 'L'
			end
			
		end
	end

	#trace
	firstseq = ""
	secondseq = ""
	cell_ = [x for x in size(matrix_)]
	for i = 0:max(length(seq1), length(seq2))
		x = cell_[1]
		y = cell_[2]	
		if trace_[x,y] == 'U'
			secondseq *= seq2[x-1]
			firstseq *= "-"
			cell_ = [x-1,y]
		elseif trace_[x,y] == 'L'
			firstseq *= seq1[y-1]
			secondseq *= "-"
			cell_ = [x, y-1]
		elseif trace_[x,y] == 'D'
			secondseq *= seq2[x-1]
			firstseq *= seq1[y-1]
			cell_ = [x-1, y-1]
		end
		
	end
	
	println(reverse(firstseq))
	println(reverse(secondseq))
end
nw(ARGS[1], ARGS[2])
